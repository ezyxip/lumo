package mesh

import (
	"bytes"
	"fmt"
	"math"
	"os"
	"strings"
	"testing"

	"github.com/ezyxip/lumo/codegen/internal/ctest"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

// Synthetic constant K=1 correlations isolate the MESH balance arithmetic.
// These fixtures have dependent E/S rows and must not be used to test Newton.
// Nonlinear full-rank coverage is provided separately in nonlinear_test.go.
func testConfig(stages int, feedStages ...int) contractv1.CompilerInput {
	cfg := contractv1.CompilerInput{
		SchemaVersion: 1, Name: "mesh_algebra",
		Thermo: contractv1.Thermo{PhaseEquilibrium: "antoine_raoult", Enthalpy: "ideal"},
		Column: contractv1.Column{Stages: stages,
			Condenser: contractv1.Condenser{Type: "total"}, Reboiler: contractv1.Reboiler{Type: "equilibrium"},
			Pressure: contractv1.Pressure{Type: "constant", ValuePa: 133322.4}},
	}
	for i := 0; i < 2; i++ {
		cfg.Components = append(cfg.Components, contractv1.Component{
			ID: fmt.Sprintf("component_%d", i), Name: "Synthetic", MolarMassGMol: 50,
			Properties: contractv1.Properties{Source: "inline",
				Antoine: &contractv1.Antoine{A: 3, B: 0, C: 1, PressureUnit: "mmHg", TemperatureUnit: "C"},
				Enthalpy: &contractv1.ConstantCP{Type: "constant_cp", ReferenceTemperatureK: 300,
					LiquidCPJMolK: float64(2 + 2*i), VaporCPJMolK: float64(3 + 2*i),
					VaporizationEnthalpyAtReferenceJMol: float64(1000 + 1000*i)},
				TemperatureRangeK: &contractv1.TemperatureRange{Min: 300, Max: 400}}})
	}
	for i, stage := range feedStages {
		cfg.Column.Feeds = append(cfg.Column.Feeds, contractv1.FeedLocation{ID: fmt.Sprintf("feed_%d", i), Stage: stage})
	}
	return cfg
}

func TestGenerateAndRunC99(t *testing.T) {
	ternary := testConfig(3, 2)
	third := ternary.Components[1]
	third.ID = "component_2"
	ternary.Components = append(ternary.Components, third)
	for _, tc := range []struct {
		name string
		cfg  contractv1.CompilerInput
	}{
		{"3x2", testConfig(3, 2)},
		{"5x2", testConfig(5, 3)},
		{"3x3", ternary},
		{"multiple_feeds", testConfig(5, 2, 4)},
		{"same_stage_feeds", testConfig(5, 3, 3)},
	} {
		t.Run(tc.name, func(t *testing.T) {
			files, err := Generate(tc.cfg)
			if err != nil {
				t.Fatal(err)
			}
			again, err := Generate(tc.cfg)
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(files.Header, again.Header) || !bytes.Equal(files.Source, again.Source) ||
				!bytes.Equal(files.Thermo.Header, again.Thermo.Header) || !bytes.Equal(files.Thermo.Source, again.Thermo.Source) {
				t.Fatal("generation is not deterministic")
			}
			if files.Variables != tc.cfg.Column.Stages*(2*len(tc.cfg.Components)+3)+1 {
				t.Fatalf("unexpected dimension %d", files.Variables)
			}
			for i, f := range tc.cfg.Column.Feeds {
				if files.FeedIDs[i] != f.ID {
					t.Fatal("feed order changed")
				}
			}
			harness, err := os.ReadFile("testdata/check.c")
			if err != nil {
				t.Fatal(err)
			}
			program := ctest.Build(t, map[string][]byte{
				"lumo_mesh.h": files.Header, "lumo_mesh.c": files.Source,
				"lumo_thermo.h": files.Thermo.Header, "lumo_thermo.c": files.Thermo.Source,
				"check.c": harness,
			})
			names := []string{"zero_residual", "summations", "heat_signs", "overall_balances", "closure", "errors"}
			if tc.cfg.Column.Stages == 3 && len(tc.cfg.Components) == 2 && len(tc.cfg.Column.Feeds) == 1 {
				names = append(names, "hand_calculated")
			}
			for _, name := range names {
				t.Run(name, func(t *testing.T) { program.Run(t, name) })
			}
		})
	}
}

func TestGenerateRejectsUnsupportedInputs(t *testing.T) {
	for _, tc := range []struct {
		name string
		edit func(*contractv1.CompilerInput)
		want string
	}{
		{"stages", func(c *contractv1.CompilerInput) { c.Column.Stages = 2 }, "stages"},
		{"condenser", func(c *contractv1.CompilerInput) { c.Column.Condenser.Type = "partial" }, "total condenser"},
		{"reboiler", func(c *contractv1.CompilerInput) { c.Column.Reboiler.Type = "other" }, "equilibrium reboiler"},
		{"pressure type", func(c *contractv1.CompilerInput) { c.Column.Pressure.Type = "linear_drop" }, "pressure"},
		{"pressure zero", func(c *contractv1.CompilerInput) { c.Column.Pressure.ValuePa = 0 }, "pressure"},
		{"pressure NaN", func(c *contractv1.CompilerInput) { c.Column.Pressure.ValuePa = math.NaN() }, "pressure"},
		{"missing feeds", func(c *contractv1.CompilerInput) { c.Column.Feeds = nil }, "feed"},
		{"feed id", func(c *contractv1.CompilerInput) { c.Column.Feeds[0].ID = "bad id" }, "id"},
		{"duplicate feed", func(c *contractv1.CompilerInput) { c.Column.Feeds = append(c.Column.Feeds, c.Column.Feeds[0]) }, "duplicate"},
		{"condenser feed", func(c *contractv1.CompilerInput) { c.Column.Feeds[0].Stage = 1 }, "stage"},
		{"reboiler feed", func(c *contractv1.CompilerInput) { c.Column.Feeds[0].Stage = 3 }, "stage"},
		{"thermo", func(c *contractv1.CompilerInput) { c.Thermo.PhaseEquilibrium = "unknown" }, "thermo"},
		{"unresolved properties", func(c *contractv1.CompilerInput) { c.Components[0].Properties.Source = "builtin" }, "inline"},
		{"huge system", func(c *contractv1.CompilerInput) { c.Column.Stages = 2147483647 }, "dimension"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			cfg := testConfig(3, 2)
			tc.edit(&cfg)
			files, err := Generate(cfg)
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Fatalf("error = %v, want %q", err, tc.want)
			}
			if len(files.Header) != 0 || len(files.Source) != 0 {
				t.Fatal("invalid input produced code")
			}
		})
	}
}
