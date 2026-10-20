package thermo

import (
	"fmt"
	"math"
	"os"
	"testing"

	"github.com/ezyxip/lumo/codegen/internal/ctest"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

func TestFlashC99(t *testing.T) {
	harness, err := os.ReadFile("testdata/flash.c")
	if err != nil {
		t.Fatal(err)
	}
	for _, ks := range [][]float64{{2, 0.5}, {1.685, 0.742, 0.532}} {
		t.Run(fmt.Sprintf("%d components", len(ks)), func(t *testing.T) {
			cfg := contractv1.CompilerInput{SchemaVersion: 1, Thermo: contractv1.Thermo{PhaseEquilibrium: "antoine_raoult", Enthalpy: "ideal"}}
			for i, k := range ks {
				cfg.Components = append(cfg.Components, contractv1.Component{ID: fmt.Sprintf("c%d", i), Properties: contractv1.Properties{Source: "inline",
					Antoine:           &contractv1.Antoine{A: math.Log10(k*101325/133.3224) + 1000/(200+350-273.15), B: 1000, C: 200, PressureUnit: "mmHg", TemperatureUnit: "C"},
					Enthalpy:          &contractv1.ConstantCP{Type: "constant_cp", ReferenceTemperatureK: 350, LiquidCPJMolK: 10, VaporCPJMolK: 20, VaporizationEnthalpyAtReferenceJMol: 20000},
					TemperatureRangeK: &contractv1.TemperatureRange{Min: 330, Max: 400},
				}})
			}
			files, err := Generate(cfg)
			if err != nil {
				t.Fatal(err)
			}
			program := ctest.Build(t, map[string][]byte{"lumo_thermo.h": files.Header, "lumo_thermo.c": files.Source, "check.c": harness})
			for _, name := range []string{"reference", "state_roundtrip", "single_phase", "errors", "pure_component"} {
				t.Run(name, func(t *testing.T) { program.Run(t, name) })
			}
		})
	}
}

func TestFlashConstantKLimit(t *testing.T) {
	cfg := testConfig() // The second component has constant K=1 at 133322.4 Pa.
	files, err := Generate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	harness := []byte(`#include "lumo_thermo.h"
#include "ctest.h"
static int constant_k(void){
    const double z[]={0,1};
    struct lumo_flash_result result,before;
    memset(&result,0x23,sizeof result);memcpy(&before,&result,sizeof result);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_TEMPERATURE,350,133322.4,z,&result),LUMO_THERMO_FLASH_AMBIGUOUS);
    CTEST_STATUS(lumo_thermo_flash(LUMO_FEED_VAPOR_FRACTION,0.5,133322.4,z,&result),LUMO_THERMO_INVALID_INPUT);
    CTEST_STATUS(memcmp(&before,&result,sizeof result),0);
    return 0;
}
int main(int argc,char **argv){
    const struct ctest_case cases[]={{"constant_k",constant_k}};
    return ctest_run(argc,argv,cases,1);
}
`)
	program := ctest.Build(t, map[string][]byte{"lumo_thermo.h": files.Header, "lumo_thermo.c": files.Source, "check.c": harness})
	program.Run(t, "constant_k")
}

func TestFlashNumericalBoundaries(t *testing.T) {
	harness, err := os.ReadFile("testdata/flash_boundaries.c")
	if err != nil {
		t.Fatal(err)
	}
	for _, tc := range []struct {
		name string
		k, z [2]float64
	}{
		{"close_k", [2]float64{1.000001, 0.999999}, [2]float64{0.50000025, 0.49999975}},
		{"trace_vapor", [2]float64{1e8, 0.5}, [2]float64{1.1e-8, 1 - 1.1e-8}},
		{"trace_liquid", [2]float64{1e-8, 2}, [2]float64{1.1e-8, 1 - 1.1e-8}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			cfg := testConfig()
			for i, k := range tc.k {
				cfg.Components[i].Properties.Antoine = &contractv1.Antoine{
					A: math.Log10(k*101325/133.3224) + 1000/(200+350-273.15),
					B: 1000, C: 200, PressureUnit: "mmHg", TemperatureUnit: "C",
				}
			}
			files, err := Generate(cfg)
			if err != nil {
				t.Fatal(err)
			}
			header := []byte(fmt.Sprintf("static const double z[] = {%.17g, %.17g};\n#define TRACE_LIQUID %d\n", tc.z[0], tc.z[1], map[bool]int{false: 0, true: 1}[tc.name == "trace_liquid"]))
			program := ctest.Build(t, map[string][]byte{
				"lumo_thermo.h": files.Header, "lumo_thermo.c": files.Source,
				"check.c": harness, "case.h": header,
			})
			program.Run(t, "numerical_boundary")
		})
	}
}
