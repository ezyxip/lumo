package thermo

import (
	"bytes"
	"math"
	"os"
	"strings"
	"testing"

	"github.com/ezyxip/lumo/codegen/internal/ctest"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

// Only the water vapor-pressure correlation is physical reference data.
// NIST WebBook, water, Bridgeman and Aldrich 1964, valid at 344..373 K:
// log10(P/bar) = 5.08354 - 1663.125/(T_K - 45.622).
// https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4
// Conversion to the contract's mmHg/C form is explicit below.
// Heat capacities, reference enthalpies and the second component are synthetic
// values for checking algebra, not a validated water-mixture property package.
func testConfig() contractv1.CompilerInput {
	return contractv1.CompilerInput{
		SchemaVersion: 1,
		Thermo:        contractv1.Thermo{PhaseEquilibrium: "antoine_raoult", Enthalpy: "ideal"},
		Components: []contractv1.Component{
			{ID: "water", Properties: contractv1.Properties{
				Source: "inline",
				Antoine: &contractv1.Antoine{A: 5.08354 + math.Log10(100000/133.3224), B: 1663.125,
					C: 273.15 - 45.622, PressureUnit: "mmHg", TemperatureUnit: "C"},
				Enthalpy: &contractv1.ConstantCP{Type: "constant_cp", ReferenceTemperatureK: 350,
					LiquidCPJMolK: 100, VaporCPJMolK: 50, VaporizationEnthalpyAtReferenceJMol: 30000},
				TemperatureRangeK: &contractv1.TemperatureRange{Min: 344, Max: 373},
			}},
			{ID: "synthetic", Properties: contractv1.Properties{
				Source:  "inline",
				Antoine: &contractv1.Antoine{A: 3, B: 0, C: 1, PressureUnit: "mmHg", TemperatureUnit: "C"},
				Enthalpy: &contractv1.ConstantCP{Type: "constant_cp", ReferenceTemperatureK: 350,
					LiquidCPJMolK: 200, VaporCPJMolK: 80, VaporizationEnthalpyAtReferenceJMol: 40000},
				TemperatureRangeK: &contractv1.TemperatureRange{Min: 344, Max: 373},
			}},
		},
	}
}

func TestGenerateAndRunC99(t *testing.T) {
	files, err := Generate(testConfig())
	if err != nil {
		t.Fatal(err)
	}
	again, err := Generate(testConfig())
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(files.Header, again.Header) || !bytes.Equal(files.Source, again.Source) {
		t.Fatal("generation is not deterministic")
	}
	if strings.Join(files.ComponentIDs, ",") != "water,synthetic" {
		t.Fatalf("unexpected component order: %v", files.ComponentIDs)
	}
	harness, err := os.ReadFile("testdata/check.c")
	if err != nil {
		t.Fatal(err)
	}
	program := ctest.Build(t, map[string][]byte{
		"lumo_thermo.h": files.Header,
		"lumo_thermo.c": files.Source,
		"check.c":       harness,
	})
	for _, name := range []string{
		"water_reference", "enthalpy", "pressure", "reference_temperature",
		"newton_iterates", "temperature_range", "invalid_input",
		"nonfinite_composition", "numeric_overflow",
	} {
		t.Run(name, func(t *testing.T) {
			program.Run(t, name)
		})
	}
}

func TestGenerateRejectsInvalidProperties(t *testing.T) {
	for _, tc := range []struct {
		name string
		edit func(*contractv1.CompilerInput)
		want string
	}{
		{"schema", func(c *contractv1.CompilerInput) { c.SchemaVersion = 0 }, "schema_version"},
		{"phase model", func(c *contractv1.CompilerInput) { c.Thermo.PhaseEquilibrium = "nrtl" }, "unsupported model"},
		{"enthalpy model", func(c *contractv1.CompilerInput) { c.Thermo.Enthalpy = "other" }, "unsupported model"},
		{"components", func(c *contractv1.CompilerInput) { c.Components = nil }, "at least two"},
		{"duplicate", func(c *contractv1.CompilerInput) { c.Components[1].ID = "water" }, "duplicate ID"},
		{"id", func(c *contractv1.CompilerInput) { c.Components[0].ID = "bad id" }, "invalid"},
		{"builtin", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Source = "builtin" }, "resolved inline"},
		{"mixed source", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Key = "water" }, "cannot contain"},
		{"missing Antoine", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Antoine = nil }, "required"},
		{"missing enthalpy", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Enthalpy = nil }, "required"},
		{"missing range", func(c *contractv1.CompilerInput) { c.Components[0].Properties.TemperatureRangeK = nil }, "required"},
		{"NaN", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Antoine.A = math.NaN() }, "finite"},
		{"infinite range", func(c *contractv1.CompilerInput) { c.Components[0].Properties.TemperatureRangeK.Max = math.Inf(1) }, "finite"},
		{"pressure unit", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Antoine.PressureUnit = "bar" }, "mmHg"},
		{"temperature unit", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Antoine.TemperatureUnit = "K" }, "temperature_unit"},
		{"Cp model", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Enthalpy.Type = "polynomial" }, "unsupported model"},
		{"Cp", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Enthalpy.LiquidCPJMolK = 0 }, "positive"},
		{"range", func(c *contractv1.CompilerInput) { c.Components[0].Properties.TemperatureRangeK.Min = 373 }, "min < max"},
		{"reference", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Enthalpy.ReferenceTemperatureK = 300 }, "reference temperature"},
		{"pole", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Antoine.C = 273.15 - 350 }, "pole"},
		{"no common range", func(c *contractv1.CompilerInput) {
			c.Components[1].Properties.TemperatureRangeK = &contractv1.TemperatureRange{Min: 400, Max: 450}
			c.Components[1].Properties.Enthalpy.ReferenceTemperatureK = 420
		}, "common interval"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			cfg := testConfig()
			tc.edit(&cfg)
			files, err := Generate(cfg)
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Fatalf("Generate error = %v, want %q", err, tc.want)
			}
			if len(files.Header) != 0 || len(files.Source) != 0 {
				t.Fatal("invalid input produced code")
			}
		})
	}
}

func TestGenerateMultipleComponents(t *testing.T) {
	cfg := testConfig()
	third := cfg.Components[1]
	third.ID = "third"
	cfg.Components = append(cfg.Components, third)
	files, err := Generate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Contains(files.Header, []byte("#define LUMO_THERMO_COMPONENTS 3")) {
		t.Fatal("component count was not specialized")
	}
}
