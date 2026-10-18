package contractv1

import (
	"errors"
	"math"
	"os"
	"strings"
	"testing"
)

func TestValidateCompilerFixtures(t *testing.T) {
	if err := ValidateCompilerInput(fixtureCompiler(t)); err != nil {
		t.Fatal(err)
	}
	data, err := os.Open("../../examples/v1/benzene_toluene/compiler_input.json")
	if err != nil {
		t.Fatal(err)
	}
	defer data.Close()
	builtin, err := ParseCompilerInput(data)
	if err != nil {
		t.Fatal(err)
	}
	if err := ValidateCompilerInput(builtin); err != nil {
		t.Fatal(err)
	}
	c := fixtureCompiler(t)
	c.Column.Feeds = append(c.Column.Feeds, FeedLocation{ID: "second", Stage: 2})
	third := c.Components[1]
	third.ID = "third"
	c.Components = append(c.Components, third)
	c.Components[0].Properties.Antoine.A = -1
	c.Components[0].Properties.Antoine.B = 0
	c.Components[0].Properties.Antoine.C = -1
	if err := ValidateCompilerInput(c); err != nil {
		t.Fatal(err)
	}
}

func TestValidateCompilerErrors(t *testing.T) {
	for _, tc := range []struct {
		name string
		edit func(*CompilerInput)
		want string
	}{
		{"version", func(c *CompilerInput) { c.SchemaVersion = 2 }, "schema_version"},
		{"name", func(c *CompilerInput) { c.Name = "" }, "name"},
		{"few components", func(c *CompilerInput) { c.Components = c.Components[:1] }, "at least two"},
		{"invalid id", func(c *CompilerInput) { c.Components[0].ID = "1bad" }, "invalid ID"},
		{"duplicate id", func(c *CompilerInput) { c.Components[1].ID = c.Components[0].ID }, "duplicate ID"},
		{"component name", func(c *CompilerInput) { c.Components[0].Name = "" }, ".name"},
		{"mass", func(c *CompilerInput) { c.Components[0].MolarMassGMol = 0 }, "molar_mass"},
		{"mass nan", func(c *CompilerInput) { c.Components[0].MolarMassGMol = math.NaN() }, "molar_mass"},
		{"model", func(c *CompilerInput) { c.Thermo.PhaseEquilibrium = "nrtl" }, "thermo.phase_equilibrium"},
		{"enthalpy model", func(c *CompilerInput) { c.Thermo.Enthalpy = "x" }, "thermo.enthalpy"},
		{"source", func(c *CompilerInput) { c.Components[0].Properties.Source = "x" }, ".source"},
		{"mixed properties", func(c *CompilerInput) { c.Components[0].Properties.Key = "x" }, "cannot contain builtin"},
		{"builtin missing", func(c *CompilerInput) { c.Components[0].Properties = Properties{Source: "builtin"} }, ".key"},
		{"builtin mixed", func(c *CompilerInput) { c.Components[0].Properties.Source = "builtin" }, "cannot contain inline"},
		{"missing antoine", func(c *CompilerInput) { c.Components[0].Properties.Antoine = nil }, ".antoine"},
		{"missing enthalpy", func(c *CompilerInput) { c.Components[0].Properties.Enthalpy = nil }, ".enthalpy"},
		{"missing range", func(c *CompilerInput) { c.Components[0].Properties.TemperatureRangeK = nil }, "temperature_range"},
		{"nan coefficient", func(c *CompilerInput) { c.Components[0].Properties.Antoine.A = math.NaN() }, ".antoine.A"},
		{"infinite coefficient", func(c *CompilerInput) { c.Components[0].Properties.Antoine.B = math.Inf(1) }, ".antoine.B"},
		{"pressure units", func(c *CompilerInput) { c.Components[0].Properties.Antoine.PressureUnit = "bar" }, "pressure_unit"},
		{"temperature units", func(c *CompilerInput) { c.Components[0].Properties.Antoine.TemperatureUnit = "K" }, "temperature_unit"},
		{"pole", func(c *CompilerInput) { c.Components[0].Properties.Antoine.C = 273.15 - 360 }, "pole"},
		{"range reversed", func(c *CompilerInput) { c.Components[0].Properties.TemperatureRangeK.Min = 400 }, "min < max"},
		{"range zero", func(c *CompilerInput) { c.Components[0].Properties.TemperatureRangeK.Min = 0 }, ".min"},
		{"range infinity", func(c *CompilerInput) { c.Components[0].Properties.TemperatureRangeK.Max = math.Inf(1) }, ".max"},
		{"reference outside", func(c *CompilerInput) { c.Components[0].Properties.Enthalpy.ReferenceTemperatureK = 400 }, "inside temperature range"},
		{"caloric model", func(c *CompilerInput) { c.Components[0].Properties.Enthalpy.Type = "x" }, "enthalpy.type"},
		{"cp zero", func(c *CompilerInput) { c.Components[0].Properties.Enthalpy.LiquidCPJMolK = 0 }, "liquid_cp"},
		{"cp infinity", func(c *CompilerInput) { c.Components[0].Properties.Enthalpy.VaporCPJMolK = math.Inf(1) }, "vapor_cp"},
		{"latent zero", func(c *CompilerInput) { c.Components[0].Properties.Enthalpy.VaporizationEnthalpyAtReferenceJMol = 0 }, "vaporization_enthalpy"},
		{"latent upper", func(c *CompilerInput) {
			p := &c.Components[0].Properties
			p.Enthalpy.ReferenceTemperatureK = 350
			p.Enthalpy.VaporizationEnthalpyAtReferenceJMol = 1000
			p.Enthalpy.LiquidCPJMolK = 140
			p.Enthalpy.VaporCPJMolK = 90
			p.TemperatureRangeK.Min = 330
			p.TemperatureRangeK.Max = 400
		}, "throughout temperature_range"},
		{"latent lower", func(c *CompilerInput) {
			p := &c.Components[0].Properties
			p.Enthalpy.ReferenceTemperatureK = 380
			p.Enthalpy.VaporizationEnthalpyAtReferenceJMol = 1000
			p.Enthalpy.LiquidCPJMolK = 90
			p.Enthalpy.VaporCPJMolK = 140
		}, "throughout temperature_range"},
		{"disjoint", func(c *CompilerInput) {
			p := &c.Components[1].Properties
			p.TemperatureRangeK.Min = 390
			p.TemperatureRangeK.Max = 400
			p.Enthalpy.ReferenceTemperatureK = 390
		}, "no common interval"},
		{"stages", func(c *CompilerInput) { c.Column.Stages = 2 }, "column.stages"},
		{"condenser", func(c *CompilerInput) { c.Column.Condenser.Type = "partial" }, "condenser.type"},
		{"reboiler", func(c *CompilerInput) { c.Column.Reboiler.Type = "kettle" }, "reboiler.type"},
		{"pressure type", func(c *CompilerInput) { c.Column.Pressure.Type = "linear_drop" }, "pressure.type"},
		{"pressure negative", func(c *CompilerInput) { c.Column.Pressure.ValuePa = -1 }, "pressure.value_Pa"},
		{"no feeds", func(c *CompilerInput) { c.Column.Feeds = nil }, "at least one feed"},
		{"feed id", func(c *CompilerInput) { c.Column.Feeds[0].ID = "bad id" }, "invalid ID"},
		{"duplicate feeds", func(c *CompilerInput) { c.Column.Feeds = append(c.Column.Feeds, c.Column.Feeds[0]) }, "duplicate ID"},
		{"feed top", func(c *CompilerInput) { c.Column.Feeds[0].Stage = 1 }, "[2, N-1]"},
		{"feed bottom", func(c *CompilerInput) { c.Column.Feeds[0].Stage = c.Column.Stages }, "[2, N-1]"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			c := fixtureCompiler(t)
			tc.edit(&c)
			wantValidation(t, ValidateCompilerInput(c), tc.want)
		})
	}
}

func wantValidation(t *testing.T, err error, fragment string) {
	t.Helper()
	var validation *ValidationError
	if !errors.As(err, &validation) || !strings.Contains(err.Error(), fragment) {
		t.Fatalf("error = %v, want ValidationError containing %q", err, fragment)
	}
}

func TestCollectAllIssues(t *testing.T) {
	c := fixtureCompiler(t)
	c.Name = ""
	c.SchemaVersion = 0
	c.Components[0].MolarMassGMol = -1
	c.Column.Feeds[0].Stage = 0
	err := ValidateCompilerInput(c)
	var validation *ValidationError
	if !errors.As(err, &validation) || len(validation.Issues) != 4 {
		t.Fatalf("expected four issues, got %v", err)
	}
	want := "invalid compiler input: 4 issue(s):\n  - schema_version: must be 1\n  - name: required field is missing\n  - components[0].molar_mass_g_mol: must be finite and positive\n  - column.feeds[0].stage: must be in [2, N-1]"
	if err.Error() != want {
		t.Fatalf("unexpected error format %q", err.Error())
	}
}

func TestNonfiniteCompilerNumbers(t *testing.T) {
	for _, x := range []float64{math.NaN(), math.Inf(1), math.Inf(-1)} {
		for _, edit := range []func(*CompilerInput, float64){
			func(c *CompilerInput, x float64) { c.Components[0].MolarMassGMol = x },
			func(c *CompilerInput, x float64) { c.Column.Pressure.ValuePa = x },
			func(c *CompilerInput, x float64) { c.Components[0].Properties.Antoine.A = x },
			func(c *CompilerInput, x float64) { c.Components[0].Properties.Antoine.B = x },
			func(c *CompilerInput, x float64) { c.Components[0].Properties.Antoine.C = x },
			func(c *CompilerInput, x float64) { c.Components[0].Properties.Enthalpy.ReferenceTemperatureK = x },
			func(c *CompilerInput, x float64) { c.Components[0].Properties.Enthalpy.LiquidCPJMolK = x },
			func(c *CompilerInput, x float64) { c.Components[0].Properties.Enthalpy.VaporCPJMolK = x },
			func(c *CompilerInput, x float64) {
				c.Components[0].Properties.Enthalpy.VaporizationEnthalpyAtReferenceJMol = x
			},
			func(c *CompilerInput, x float64) { c.Components[0].Properties.TemperatureRangeK.Min = x },
			func(c *CompilerInput, x float64) { c.Components[0].Properties.TemperatureRangeK.Max = x },
		} {
			c := fixtureCompiler(t)
			edit(&c, x)
			if err := ValidateCompilerInput(c); err == nil {
				t.Fatal("accepted nonfinite number")
			}
		}
	}
}
