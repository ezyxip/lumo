package contractv1

import (
	"fmt"
	"math"
)

// ValidateModelInput checks operating parameters against a valid description.
// It does not flash feeds or resolve builtin properties. Flow comparisons use
// total molar flow for molar D and total mass flow for mass D, so the unknown
// distillate composition is never replaced with a feed composition.
func ValidateModelInput(c CompilerInput, input ModelInput) error {
	if err := ValidateCompilerInput(c); err != nil {
		return err
	}
	v := &validator{}
	if input.SchemaVersion != SchemaVersion {
		v.add("schema_version", "must be 1")
	}
	expected := map[string]bool{}
	for _, feed := range c.Column.Feeds {
		expected[feed.ID] = true
		if _, ok := input.Feeds[feed.ID]; !ok {
			v.add(fmt.Sprintf("feeds[%q]", feed.ID), "required feed is missing")
		}
	}
	totalMol, totalMass := 0.0, 0.0
	comparable := true
	for _, id := range sortedKeys(input.Feeds) {
		path := fmt.Sprintf("feeds[%q]", id)
		if !expected[id] {
			v.add(path, "unknown feed ID")
		}
		feed := input.Feeds[id]
		start := len(v.issues)
		v.flow(path+".flow", feed.Flow)
		mw := v.composition(path+".composition", feed.Composition, c.Components)
		v.feedState(path+".state", feed.State, c.Components)
		if len(v.issues) != start {
			comparable = false
			continue
		}
		mol, mass := flowRates(feed.Flow, mw)
		if !finite(mol) || !finite(mass) || mol <= 0 || mass <= 0 {
			v.add(path+".flow", "conversion must produce finite positive molar and mass flows")
			comparable = false
		} else {
			totalMol += mol
			totalMass += mass
		}
	}
	if len(input.Feeds) == 0 {
		v.add("feeds", "at least one feed is required")
	}
	specs := input.Specifications
	v.equal("specifications.type", specs.Type, "rating")
	v.flow("specifications.distillate_flow", specs.DistillateFlow)
	v.equal("specifications.reflux_ratio.basis", specs.RefluxRatio.Basis, "mole")
	if !finite(specs.RefluxRatio.Value) || specs.RefluxRatio.Value < 0 {
		v.add("specifications.reflux_ratio.value", "must be finite and nonnegative")
	}
	if input.InitialGuess != nil {
		v.equal("initial_guess.profile", input.InitialGuess.Profile, "linear_T_x")
	}
	if !finite(totalMol) || !finite(totalMass) {
		v.add("feeds", "total flows must be finite")
		comparable = false
	}
	if comparable && len(input.Feeds) == len(expected) && finite(specs.DistillateFlow.Value) && specs.DistillateFlow.Value > 0 {
		d := specs.DistillateFlow
		value, total := d.Value, totalMol
		switch d.Unit {
		case "mol_s":
		case "kmol_h":
			value /= 3.6
		case "kg_h":
			value /= 3600
			total = totalMass
		default:
			return v.result("model input")
		}
		if !finite(value) || value <= 0 {
			v.add("specifications.distillate_flow", "converted flow must be finite and positive")
		} else if value >= total {
			v.add("specifications.distillate_flow", "must be less than total feed flow on the same basis")
		}
	}
	return v.result("model input")
}

func (v *validator) flow(path string, flow Flow) {
	v.positive(path+".value", flow.Value)
	switch flow.Unit {
	case "mol_s", "kmol_h", "kg_h":
	default:
		v.add(path+".unit", "must be mol_s, kmol_h or kg_h")
	}
}

// Returns the mixture molar mass in g/mol. The caller only uses it if valid.
func (v *validator) composition(path string, z Composition, components []Component) float64 {
	switch z.Basis {
	case "mole", "mass":
	default:
		v.add(path+".basis", "must be mole or mass")
	}
	expected := map[string]bool{}
	sum, molecularWeight := 0.0, 0.0
	for _, comp := range components {
		expected[comp.ID] = true
		fraction, ok := z.Fractions[comp.ID]
		if !ok {
			v.add(path+".fractions", fmt.Sprintf("missing component ID %q", comp.ID))
			continue
		}
		if !finite(fraction) || fraction < 0 || fraction > 1 {
			v.add(path+".fractions", fmt.Sprintf("fraction for %q must be finite and in [0, 1]", comp.ID))
			continue
		}
		sum += fraction
		if z.Basis == "mass" {
			molecularWeight += fraction / comp.MolarMassGMol
		} else {
			molecularWeight += fraction * comp.MolarMassGMol
		}
	}
	for _, id := range sortedKeys(z.Fractions) {
		if !expected[id] {
			v.add(path+".fractions", fmt.Sprintf("unknown component ID %q", id))
		}
	}
	if !finite(sum) || math.Abs(sum-1) > 1e-10 {
		v.add(path+".fractions", "fractions must sum to 1 within 1e-10")
	}
	if z.Basis == "mass" {
		molecularWeight = 1 / molecularWeight
	}
	if !finite(molecularWeight) || molecularWeight <= 0 {
		v.add(path, "mixture molar mass must be finite and positive")
	}
	return molecularWeight
}

func flowRates(flow Flow, mw float64) (mol, mass float64) {
	switch flow.Unit {
	case "mol_s":
		mol = flow.Value
		mass = mol * (mw / 1000)
	case "kmol_h":
		mol = flow.Value / 3.6
		mass = mol * (mw / 1000)
	case "kg_h":
		mass = flow.Value / 3600
		mol = mass * (1000 / mw)
	}
	return
}

func (v *validator) feedState(path string, state FeedState, components []Component) {
	switch state.Type {
	case "saturated_liquid":
		if state.TemperatureK != nil || state.VaporFraction != nil {
			v.add(path, "saturated_liquid cannot contain temperature_K or vapor_fraction")
		}
	case "temperature":
		if state.VaporFraction != nil {
			v.add(path, "temperature cannot contain vapor_fraction")
		}
		if state.TemperatureK == nil {
			v.add(path+".temperature_K", "required field is missing")
		} else if v.positive(path+".temperature_K", *state.TemperatureK) {
			for _, comp := range components {
				r := comp.Properties.TemperatureRangeK
				if comp.Properties.Source == "inline" && r != nil && (*state.TemperatureK < r.Min || *state.TemperatureK > r.Max) {
					v.add(path+".temperature_K", "outside component temperature range")
					break
				}
			}
		}
	case "vapor_fraction":
		if state.TemperatureK != nil {
			v.add(path, "vapor_fraction cannot contain temperature_K")
		}
		if state.VaporFraction == nil {
			v.add(path+".vapor_fraction", "required field is missing")
		} else if x := *state.VaporFraction; !finite(x) || x < 0 || x > 1 {
			v.add(path+".vapor_fraction", "must be finite and in [0, 1]")
		}
	default:
		v.add(path+".type", "must be saturated_liquid, temperature or vapor_fraction")
	}
}
