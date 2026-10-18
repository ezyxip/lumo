package contractv1

import (
	"fmt"
	"math"
	"testing"
)

func editFeed(input *ModelInput, edit func(*Feed)) {
	feed := input.Feeds["main_feed"]
	edit(&feed)
	input.Feeds["main_feed"] = feed
}
func number(x float64) *float64 { return &x }

func TestValidateModelErrors(t *testing.T) {
	for _, tc := range []struct {
		name string
		edit func(*ModelInput)
		want string
	}{
		{"version", func(m *ModelInput) { m.SchemaVersion = 0 }, "schema_version"},
		{"missing feed", func(m *ModelInput) { delete(m.Feeds, "main_feed") }, "required feed is missing"},
		{"unknown feed", func(m *ModelInput) { m.Feeds["extra"] = m.Feeds["main_feed"] }, "unknown feed ID"},
		{"flow zero", func(m *ModelInput) { editFeed(m, func(f *Feed) { f.Flow.Value = 0 }) }, "flow.value"},
		{"flow unit", func(m *ModelInput) { editFeed(m, func(f *Feed) { f.Flow.Unit = "kg_s" }) }, "flow.unit"},
		{"basis", func(m *ModelInput) { editFeed(m, func(f *Feed) { f.Composition.Basis = "volume" }) }, "composition.basis"},
		{"component missing", func(m *ModelInput) { delete(m.Feeds["main_feed"].Composition.Fractions, "toluene") }, "missing component ID"},
		{"component extra", func(m *ModelInput) { m.Feeds["main_feed"].Composition.Fractions["other"] = 0 }, "unknown component ID"},
		{"fraction negative", func(m *ModelInput) { m.Feeds["main_feed"].Composition.Fractions["benzene"] = -0.1 }, "in [0, 1]"},
		{"fraction large", func(m *ModelInput) { m.Feeds["main_feed"].Composition.Fractions["benzene"] = 1.1 }, "in [0, 1]"},
		{"fraction sum", func(m *ModelInput) { m.Feeds["main_feed"].Composition.Fractions["benzene"] = 0.6 }, "sum to 1"},
		{"state type", func(m *ModelInput) { editFeed(m, func(f *Feed) { f.State.Type = "q" }) }, "state.type"},
		{"saturated extra", func(m *ModelInput) { editFeed(m, func(f *Feed) { f.State.TemperatureK = number(360) }) }, "cannot contain"},
		{"temperature absent", func(m *ModelInput) { editFeed(m, func(f *Feed) { f.State = FeedState{Type: "temperature"} }) }, "temperature_K"},
		{"temperature zero", func(m *ModelInput) {
			editFeed(m, func(f *Feed) { f.State = FeedState{Type: "temperature", TemperatureK: number(0)} })
		}, "finite and positive"},
		{"temperature outside", func(m *ModelInput) {
			editFeed(m, func(f *Feed) { f.State = FeedState{Type: "temperature", TemperatureK: number(349)} })
		}, "outside component"},
		{"temperature mixed", func(m *ModelInput) {
			editFeed(m, func(f *Feed) {
				f.State = FeedState{Type: "temperature", TemperatureK: number(360), VaporFraction: number(0)}
			})
		}, "cannot contain"},
		{"vapor absent", func(m *ModelInput) { editFeed(m, func(f *Feed) { f.State = FeedState{Type: "vapor_fraction"} }) }, "vapor_fraction"},
		{"vapor large", func(m *ModelInput) {
			editFeed(m, func(f *Feed) { f.State = FeedState{Type: "vapor_fraction", VaporFraction: number(1.1)} })
		}, "in [0, 1]"},
		{"vapor negative", func(m *ModelInput) {
			editFeed(m, func(f *Feed) { f.State = FeedState{Type: "vapor_fraction", VaporFraction: number(-0.1)} })
		}, "in [0, 1]"},
		{"vapor mixed", func(m *ModelInput) {
			editFeed(m, func(f *Feed) {
				f.State = FeedState{Type: "vapor_fraction", VaporFraction: number(0), TemperatureK: number(360)}
			})
		}, "cannot contain"},
		{"spec type", func(m *ModelInput) { m.Specifications.Type = "design" }, "specifications.type"},
		{"D zero", func(m *ModelInput) { m.Specifications.DistillateFlow.Value = 0 }, "distillate_flow.value"},
		{"D unit", func(m *ModelInput) { m.Specifications.DistillateFlow.Unit = "m3_h" }, "distillate_flow.unit"},
		{"D equals F", func(m *ModelInput) { m.Specifications.DistillateFlow.Value = 100 }, "less than total feed"},
		{"D exceeds F", func(m *ModelInput) { m.Specifications.DistillateFlow.Value = 101 }, "less than total feed"},
		{"R negative", func(m *ModelInput) { m.Specifications.RefluxRatio.Value = -1 }, "reflux_ratio.value"},
		{"R basis", func(m *ModelInput) { m.Specifications.RefluxRatio.Basis = "mass" }, "reflux_ratio.basis"},
		{"guess", func(m *ModelInput) { m.InitialGuess = &InitialGuess{Profile: "x"} }, "initial_guess.profile"},
		{"conversion overflow", func(m *ModelInput) {
			editFeed(m, func(f *Feed) {
				f.Flow = Flow{Value: math.MaxFloat64, Unit: "mol_s"}
				f.Composition.Fractions["benzene"] = 1
				f.Composition.Fractions["toluene"] = 0
			})
			m.Feeds["extra"] = m.Feeds["main_feed"]
		}, "total flows must be finite"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			c, m := fixtureCompiler(t), fixtureModel(t)
			tc.edit(&m)
			wantValidation(t, ValidateModelInput(c, m), tc.want)
		})
	}
}

func TestValidateModelVariants(t *testing.T) {
	for _, state := range []FeedState{
		{Type: "saturated_liquid"},
		{Type: "temperature", TemperatureK: number(350)},
		{Type: "temperature", TemperatureK: number(390)},
		{Type: "vapor_fraction", VaporFraction: number(0)},
		{Type: "vapor_fraction", VaporFraction: number(1)},
	} {
		c, m := fixtureCompiler(t), fixtureModel(t)
		editFeed(&m, func(f *Feed) { f.State = state })
		m.Specifications.RefluxRatio.Value = 0
		m.InitialGuess = &InitialGuess{Profile: "linear_T_x"}
		if err := ValidateModelInput(c, m); err != nil {
			t.Fatal(err)
		}
	}
	c, m := fixtureCompiler(t), fixtureModel(t)
	c.Column.Feeds = append(c.Column.Feeds, FeedLocation{ID: "second", Stage: 2})
	m.Feeds["second"] = m.Feeds["main_feed"]
	m.Specifications.DistillateFlow.Value = 150
	if err := ValidateModelInput(c, m); err != nil {
		t.Fatal(err)
	}
	delete(m.Feeds, "second")
	wantValidation(t, ValidateModelInput(c, m), "required feed is missing")
}

func TestFlowBasisConversions(t *testing.T) {
	for _, basis := range []string{"mole", "mass"} {
		mw := 150.0
		if basis == "mass" {
			mw = 1 / (0.5/100 + 0.5/200)
		}
		for _, unit := range []string{"mol_s", "kmol_h", "kg_h"} {
			t.Run(basis+"/"+unit, func(t *testing.T) {
				c, m := fixtureCompiler(t), fixtureModel(t)
				c.Components[0].MolarMassGMol = 100
				c.Components[1].MolarMassGMol = 200
				value := 2.0
				switch unit {
				case "kmol_h":
					value *= 3.6
				case "kg_h":
					value *= mw * 3.6
				}
				editFeed(&m, func(f *Feed) { f.Flow = Flow{Value: value, Unit: unit}; f.Composition.Basis = basis })
				for _, distillate := range []Flow{{Value: 1, Unit: "mol_s"}, {Value: 3.6, Unit: "kmol_h"}, {Value: mw * 3.6, Unit: "kg_h"}} {
					m.Specifications.DistillateFlow = distillate
					if err := ValidateModelInput(c, m); err != nil {
						t.Fatal(err)
					}
					m.Specifications.DistillateFlow.Value *= 2.01
					wantValidation(t, ValidateModelInput(c, m), "less than total feed")
				}
			})
		}
	}
}

func TestCompositionToleranceAndZeros(t *testing.T) {
	c, m := fixtureCompiler(t), fixtureModel(t)
	fractions := m.Feeds["main_feed"].Composition.Fractions
	fractions["benzene"] = 1
	fractions["toluene"] = 0
	if err := ValidateModelInput(c, m); err != nil {
		t.Fatal(err)
	}
	fractions["benzene"] = 0.5 + 5e-11
	fractions["toluene"] = 0.5
	if err := ValidateModelInput(c, m); err != nil {
		t.Fatal(err)
	}
	fractions["benzene"] = 0.5 + 2e-10
	wantValidation(t, ValidateModelInput(c, m), "sum to 1")
}

func TestModelErrorOrderAndAggregation(t *testing.T) {
	c, m := fixtureCompiler(t), fixtureModel(t)
	m.Feeds["z"] = Feed{}
	m.Feeds["a"] = Feed{}
	m.Specifications.Type = "x"
	first := ValidateModelInput(c, m)
	wantValidation(t, first, "specifications.type")
	for range 10 {
		if got := ValidateModelInput(c, m).Error(); got != first.Error() {
			t.Fatal("unstable error order")
		}
	}
	c.Name = ""
	wantValidation(t, ValidateModelInput(c, m), "name")
}

func TestNonfiniteModelNumbers(t *testing.T) {
	for _, x := range []float64{math.NaN(), math.Inf(1), math.Inf(-1)} {
		for i, edit := range []func(*ModelInput){
			func(m *ModelInput) { editFeed(m, func(f *Feed) { f.Flow.Value = x }) },
			func(m *ModelInput) { m.Feeds["main_feed"].Composition.Fractions["benzene"] = x },
			func(m *ModelInput) {
				editFeed(m, func(f *Feed) { f.State = FeedState{Type: "temperature", TemperatureK: &x} })
			},
			func(m *ModelInput) {
				editFeed(m, func(f *Feed) { f.State = FeedState{Type: "vapor_fraction", VaporFraction: &x} })
			},
			func(m *ModelInput) { m.Specifications.DistillateFlow.Value = x },
			func(m *ModelInput) { m.Specifications.RefluxRatio.Value = x },
		} {
			t.Run(fmt.Sprintf("%g/%d", x, i), func(t *testing.T) {
				c, m := fixtureCompiler(t), fixtureModel(t)
				edit(&m)
				if err := ValidateModelInput(c, m); err == nil {
					t.Fatal("accepted nonfinite number")
				}
			})
		}
	}
}
