package initialguess

import (
	"bytes"
	"encoding/json"
	"fmt"
	"os"
	"strings"
	"testing"

	"github.com/ezyxip/lumo/codegen/internal/ctest"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

func config(t *testing.T) contractv1.CompilerInput {
	t.Helper()
	r, err := os.Open("../../contracts/v1/testdata/compiler_input.json")
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	cfg, err := contractv1.ParseCompilerInput(r)
	if err != nil {
		t.Fatal(err)
	}
	return cfg
}
func sources(f Files, h []byte) map[string][]byte {
	return map[string][]byte{"lumo_initialguess.h": f.Header, "lumo_initialguess.c": f.Source,
		"lumo_mesh.h": f.Mesh.Header, "lumo_mesh.c": f.Mesh.Source,
		"lumo_thermo.h": f.Mesh.Thermo.Header, "lumo_thermo.c": f.Mesh.Thermo.Source,
		"lumo_jacobian.h": f.Mesh.Jacobian.Header, "lumo_jacobian.c": f.Mesh.Jacobian.Source, "check.c": h}
}
func TestInitialGuessC99(t *testing.T) {
	h, err := os.ReadFile("testdata/check.c")
	if err != nil {
		t.Fatal(err)
	}
	for _, name := range []string{"binary", "five_stages", "ternary", "multiple_feeds", "same_stage"} {
		t.Run(name, func(t *testing.T) {
			cfg := config(t)
			switch name {
			case "five_stages":
				cfg.Column.Stages = 5
				cfg.Column.Feeds[0].Stage = 3
			case "ternary":
				c := cfg.Components[1]
				c.ID = "third"
				cfg.Components = append(cfg.Components, c)
			case "multiple_feeds":
				cfg.Column.Stages = 5
				cfg.Column.Feeds = append(cfg.Column.Feeds, contractv1.FeedLocation{ID: "second", Stage: 4})
			case "same_stage":
				cfg.Column.Feeds = append(cfg.Column.Feeds, contractv1.FeedLocation{ID: "second", Stage: 2})
			}
			f, err := Generate(cfg)
			if err != nil {
				t.Fatal(err)
			}
			again, err := Generate(cfg)
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(f.Header, again.Header) || !bytes.Equal(f.Source, again.Source) {
				t.Fatal("nondeterministic generation")
			}
			p := ctest.Build(t, sources(f, h))
			for _, test := range []string{"profiles", "zero_reflux", "pure", "errors", "enthalpy", "distillate_limits"} {
				t.Run(test, func(t *testing.T) { p.Run(t, test) })
			}
		})
	}
}
func TestRejectsInvalidConfiguration(t *testing.T) {
	cfg := config(t)
	cfg.Column.Feeds[0].Stage = 1
	if f, err := Generate(cfg); err == nil || len(f.Source) != 0 {
		t.Fatal("invalid configuration accepted")
	}
	cfg = config(t)
	cfg.Components[0].Properties = contractv1.Properties{Source: "builtin", Key: "x", Revision: "1"}
	if _, err := Generate(cfg); err == nil {
		t.Fatal("unresolved properties accepted")
	}
}
func TestNewtonFromInitialGuess(t *testing.T) {
	f, err := Generate(config(t))
	if err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile("../mesh/testdata/nonlinear/reference.json")
	if err != nil {
		t.Fatal(err)
	}
	var ref struct {
		State     []float64 `json:"state"`
		Operating struct {
			H float64 `json:"feed_enthalpy_J_mol"`
		} `json:"operating"`
	}
	if err := json.Unmarshal(data, &ref); err != nil {
		t.Fatal(err)
	}
	if len(ref.State) != f.Mesh.Variables {
		t.Fatal("invalid reference size")
	}
	var header bytes.Buffer
	header.WriteString("static const double reference_state[LUMO_MESH_SIZE]={\n")
	for _, v := range ref.State {
		fmt.Fprintf(&header, "%.17g,\n", v)
	}
	fmt.Fprintf(&header, "};\nstatic const struct lumo_mesh_operating reference_operating={.feeds={{100,{0.5,0.5},%.17g}},.distillate_mol_s=50,.reflux_ratio=2.5};\n", ref.Operating.H)
	h, err := os.ReadFile("../mesh/testdata/nonlinear/check.c")
	if err != nil {
		t.Fatal(err)
	}
	// Reuse the independent column's existing Newton probe and reference checks.
	if !strings.Contains(string(h), "LUMO_TEST_INITIALGUESS") {
		t.Fatal("Newton harness lacks initialguess integration")
	}
	files := sources(f, append([]byte("#define LUMO_TEST_INITIALGUESS\n"), h...))
	files["reference.h"] = header.Bytes()
	p := ctest.Build(t, files)
	p.Run(t, "newton_reference")
}

func TestNarrowTemperatureRange(t *testing.T) {
	cfg := config(t)
	for i := range cfg.Components {
		cfg.Components[i].Properties.TemperatureRangeK = &contractv1.TemperatureRange{Min: 365, Max: 365.5}
		cfg.Components[i].Properties.Enthalpy.ReferenceTemperatureK = 365
	}
	f, err := Generate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	h, err := os.ReadFile("testdata/check.c")
	if err != nil {
		t.Fatal(err)
	}
	p := ctest.Build(t, sources(f, h))
	p.Run(t, "profiles")
}
func TestBubblePointOutsideRange(t *testing.T) {
	cfg := config(t)
	for i := range cfg.Components {
		cfg.Components[i].Properties.TemperatureRangeK = &contractv1.TemperatureRange{Min: 350, Max: 355}
	}
	f, err := Generate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	h := []byte(`#include "lumo_initialguess.h"
#include "ctest.h"
static int range_error(void){
 const struct lumo_mesh_operating op={.feeds={{100,{0.5,0.5},2060}},.distillate_mol_s=50,.reflux_ratio=2.5};
 double state[LUMO_MESH_SIZE],before[LUMO_MESH_SIZE];
 memset(state,0x34,sizeof state);memcpy(before,state,sizeof state);
 CTEST_STATUS(lumo_initialguess_linear_T_x(&op,state),LUMO_THERMO_TEMPERATURE_RANGE);
 CTEST_STATUS(memcmp(state,before,sizeof state),0);return 0;
}
int main(int argc,char **argv){const struct ctest_case cases[]={{"range_error",range_error}};return ctest_run(argc,argv,cases,1);}
`)
	p := ctest.Build(t, sources(f, h))
	p.Run(t, "range_error")
}
