package prepare

import (
	"bytes"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"strings"
	"testing"

	"github.com/ezyxip/lumo/codegen/internal/ctest"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

func testConfig(t *testing.T) contractv1.CompilerInput {
	t.Helper()
	data, err := os.Open("../../contracts/v1/testdata/compiler_input.json")
	if err != nil {
		t.Fatal(err)
	}
	defer data.Close()
	cfg, err := contractv1.ParseCompilerInput(data)
	if err != nil {
		t.Fatal(err)
	}
	return cfg
}

func TestGenerateAndRunC99(t *testing.T) {
	harness, err := os.ReadFile("testdata/check.c")
	if err != nil {
		t.Fatal(err)
	}
	for _, variant := range []string{"binary", "ternary", "multiple_feeds", "same_stage_feeds"} {
		t.Run(variant, func(t *testing.T) {
			cfg := testConfig(t)
			switch variant {
			case "ternary":
				third := cfg.Components[1]
				third.ID = "third"
				cfg.Components = append(cfg.Components, third)
			case "multiple_feeds":
				cfg.Column.Stages = 5
				cfg.Column.Feeds = append(cfg.Column.Feeds, contractv1.FeedLocation{ID: "second", Stage: 4})
			case "same_stage_feeds":
				cfg.Column.Feeds = append(cfg.Column.Feeds, contractv1.FeedLocation{ID: "second", Stage: 2})
			}
			files, err := Generate(cfg)
			if err != nil {
				t.Fatal(err)
			}
			again, err := Generate(cfg)
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(files.Header, again.Header) || !bytes.Equal(files.Source, again.Source) {
				t.Fatal("nondeterministic generation")
			}
			if files.Mesh.FeedIDs[0] != "main_feed" || files.Mesh.Thermo.ComponentIDs[0] != "benzene" {
				t.Fatal("input array order changed")
			}
			program := ctest.Build(t, sourceFiles(files, harness))
			for _, name := range []string{"reference", "conversions", "feed_states", "errors", "atomic_failure", "mesh"} {
				t.Run(name, func(t *testing.T) { program.Run(t, name) })
			}
		})
	}
}

func sourceFiles(f Files, harness []byte) map[string][]byte {
	return map[string][]byte{
		"lumo_prepare.h": f.Header, "lumo_prepare.c": f.Source,
		"lumo_mesh.h": f.Mesh.Header, "lumo_mesh.c": f.Mesh.Source,
		"lumo_thermo.h": f.Mesh.Thermo.Header, "lumo_thermo.c": f.Mesh.Thermo.Source,
		"lumo_jacobian.h": f.Mesh.Jacobian.Header, "lumo_jacobian.c": f.Mesh.Jacobian.Source,
		"check.c": harness,
	}
}

func TestGenerateRejectsInvalidInputs(t *testing.T) {
	for _, tc := range []struct {
		name string
		edit func(*contractv1.CompilerInput)
		want string
	}{
		{"mass", func(c *contractv1.CompilerInput) { c.Components[0].MolarMassGMol = 0 }, "molar_mass"},
		{"nonfinite mass", func(c *contractv1.CompilerInput) { c.Components[0].MolarMassGMol = math.NaN() }, "molar_mass"},
		{"builtin", func(c *contractv1.CompilerInput) {
			c.Components[0].Properties = contractv1.Properties{Source: "builtin", Key: "benzene", Revision: "1"}
		}, "resolved inline"},
		{"stage", func(c *contractv1.CompilerInput) { c.Column.Feeds[0].Stage = 1 }, "stage"},
		{"version", func(c *contractv1.CompilerInput) { c.SchemaVersion = 0 }, "schema_version"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			cfg := testConfig(t)
			tc.edit(&cfg)
			f, err := Generate(cfg)
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Fatalf("error=%v, want %q", err, tc.want)
			}
			if len(f.Header) != 0 || len(f.Source) != 0 {
				t.Fatal("invalid input produced code")
			}
		})
	}
}

// Full MESH reference is independently computed by the scalar Python model.
// Only this integration fixture, not code generation, fixes runtime values.
func TestIndependentMeshReference(t *testing.T) {
	data, err := os.ReadFile("../mesh/testdata/nonlinear/reference.json")
	if err != nil {
		t.Fatal(err)
	}
	var reference struct {
		State []float64 `json:"state"`
	}
	if err := json.Unmarshal(data, &reference); err != nil {
		t.Fatal(err)
	}
	files, err := Generate(testConfig(t))
	if err != nil {
		t.Fatal(err)
	}
	if len(reference.State) != files.Mesh.Variables {
		t.Fatal("wrong reference dimensions")
	}
	var source bytes.Buffer
	source.WriteString("#include \"lumo_prepare.h\"\n#include \"ctest.h\"\nstatic const double state[] = {\n")
	for _, x := range reference.State {
		fmt.Fprintf(&source, "%.17g,\n", x)
	}
	source.WriteString(`};
static int independent_reference(void) {
    const struct lumo_prepare_input input = {
        .feeds = {{.flow=360, .flow_unit=LUMO_FLOW_KMOL_H,
            .composition_basis=LUMO_COMPOSITION_MASS,
            .fractions={78.11/(78.11+92.14),92.14/(78.11+92.14)},
            .state=LUMO_FEED_SATURATED_LIQUID}},
        .distillate_flow=180, .distillate_unit=LUMO_FLOW_KMOL_H, .reflux_ratio=2.5
    };
    struct lumo_mesh_operating operating;
    double r[LUMO_MESH_SIZE];
    size_t i;
    CTEST_STATUS(lumo_prepare_operating(&input,&operating),LUMO_PREPARE_OK);
    CTEST_STATUS(lumo_mesh_residual(state,&operating,r),LUMO_MESH_OK);
    for(i=0;i<LUMO_MESH_SIZE;++i){
        double scale = i%LUMO_MESH_BLOCK==LUMO_MESH_H ? 1e7 : 100;
        CTEST_NEAR(r[i]/scale,0);
    }
    return 0;
}
int main(int argc,char **argv){
    const struct ctest_case cases[]={{"independent_reference",independent_reference}};
    return ctest_run(argc,argv,cases,1);
}
`)
	program := ctest.Build(t, sourceFiles(files, source.Bytes()))
	program.Run(t, "independent_reference")
}

func TestCloseKFeed(t *testing.T) {
	cfg := testConfig(t)
	for i, k := range []float64{1.000001, 0.999999} {
		cfg.Components[i].Properties.Antoine = &contractv1.Antoine{
			A: math.Log10(k*101325/133.3224) + 1000/(200+350-273.15),
			B: 1000, C: 200, PressureUnit: "mmHg", TemperatureUnit: "C",
		}
	}
	files, err := Generate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	harness := []byte(`#include "lumo_prepare.h"
#include "ctest.h"
static int close_k_feed(void){
    const struct lumo_prepare_input input={
        .feeds={{.flow=100,.flow_unit=LUMO_FLOW_MOL_S,
            .composition_basis=LUMO_COMPOSITION_MOLE,
            .fractions={0.50000025,0.49999975},
            .state=LUMO_FEED_TEMPERATURE,.state_value=350}},
        .distillate_flow=50,.distillate_unit=LUMO_FLOW_MOL_S,.reflux_ratio=2.5
    };
    struct lumo_mesh_operating operating;
    const double unit[]={1,1},zero[]={0,0};
    double e[2],hl,hv;
    long double a,b,beta,y0,y1,expected_h;
    CTEST_STATUS(lumo_thermo_eval(350,101325,unit,zero,e,&hl,&hv),LUMO_THERMO_OK);
    a=-(long double)e[0]-1;b=-(long double)e[1]-1;
    beta=-((long double)input.feeds[0].fractions[0]*a+(long double)input.feeds[0].fractions[1]*b)/(a*b);
    y0=input.feeds[0].fractions[0]*(1+a)/(1+beta*a);
    y1=input.feeds[0].fractions[1]*(1+b)/(1+beta*b);
    expected_h=beta*(30000*y0+33000*y1);
    CTEST_STATUS(lumo_prepare_operating(&input,&operating),LUMO_PREPARE_OK);
    CTEST_NEAR(operating.feeds[0].enthalpy_J_mol,(double)expected_h);
    CTEST_NEAR(operating.feeds[0].flow_mol_s,100);
    return 0;
}
int main(int argc,char **argv){
    const struct ctest_case cases[]={{"close_k_feed",close_k_feed}};
    return ctest_run(argc,argv,cases,1);
}
`)
	program := ctest.Build(t, sourceFiles(files, harness))
	program.Run(t, "close_k_feed")
}
