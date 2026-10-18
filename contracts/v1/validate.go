package contractv1

import (
	"fmt"
	"math"
	"regexp"
	"sort"
	"strings"
)

// ValidationError collects all detected semantic problems in a stable order.
type ValidationError struct {
	Input  string
	Issues []string
}

func (e *ValidationError) Error() string {
	return fmt.Sprintf("invalid %s: %d issue(s):\n  - %s", e.Input, len(e.Issues), strings.Join(e.Issues, "\n  - "))
}

type validator struct{ issues []string }

func (v *validator) add(path, message string) { v.issues = append(v.issues, path+": "+message) }
func (v *validator) result(input string) error {
	if len(v.issues) == 0 {
		return nil
	}
	return &ValidationError{Input: input, Issues: v.issues}
}
func finite(x float64) bool { return !math.IsNaN(x) && !math.IsInf(x, 0) }
func (v *validator) positive(path string, x float64) bool {
	if !finite(x) || x <= 0 {
		v.add(path, "must be finite and positive")
		return false
	}
	return true
}
func (v *validator) equal(path, got, want string) {
	if got != want {
		v.add(path, fmt.Sprintf("must be %q, got %q", want, got))
	}
}

var identifier = regexp.MustCompile(`^[A-Za-z][A-Za-z0-9_-]*$`)

func (v *validator) id(path, id string, seen map[string]bool) {
	if !identifier.MatchString(id) {
		v.add(path, "invalid ID")
	}
	if seen[id] {
		v.add(path, "duplicate ID")
	}
	seen[id] = true
}
func sortedKeys[V any](m map[string]V) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

// ValidateCompilerInput checks the v1 description, including inline physics.
// Builtin references are checked for shape only. A property database/resolver
// is not implemented; generators still require resolved inline properties.
func ValidateCompilerInput(c CompilerInput) error {
	v := &validator{}
	if c.SchemaVersion != SchemaVersion {
		v.add("schema_version", "must be 1")
	}
	if c.Name == "" {
		v.add("name", "required field is missing")
	}
	v.equal("thermo.phase_equilibrium", c.Thermo.PhaseEquilibrium, "antoine_raoult")
	v.equal("thermo.enthalpy", c.Thermo.Enthalpy, "ideal")
	if len(c.Components) < 2 {
		v.add("components", "need at least two components")
	}
	ids := map[string]bool{}
	commonMin, commonMax := 0.0, math.Inf(1)
	inline := 0
	for i, comp := range c.Components {
		path := fmt.Sprintf("components[%d]", i)
		v.id(path+".id", comp.ID, ids)
		if comp.Name == "" {
			v.add(path+".name", "required field is missing")
		}
		v.positive(path+".molar_mass_g_mol", comp.MolarMassGMol)
		v.properties(path+".properties", comp.Properties)
		if r := comp.Properties.TemperatureRangeK; comp.Properties.Source == "inline" && r != nil && finite(r.Min) && finite(r.Max) && r.Min > 0 && r.Min < r.Max {
			commonMin = math.Max(commonMin, r.Min)
			commonMax = math.Min(commonMax, r.Max)
			inline++
		}
	}
	if inline > 1 && commonMin >= commonMax {
		v.add("components", "temperature ranges have no common interval")
	}
	col := c.Column
	if col.Stages < 3 {
		v.add("column.stages", "must be >= 3")
	}
	v.equal("column.condenser.type", col.Condenser.Type, "total")
	v.equal("column.reboiler.type", col.Reboiler.Type, "equilibrium")
	v.equal("column.pressure.type", col.Pressure.Type, "constant")
	v.positive("column.pressure.value_Pa", col.Pressure.ValuePa)
	if len(col.Feeds) == 0 {
		v.add("column.feeds", "at least one feed is required")
	}
	ids = map[string]bool{}
	for i, feed := range col.Feeds {
		path := fmt.Sprintf("column.feeds[%d]", i)
		v.id(path+".id", feed.ID, ids)
		if feed.Stage < 2 || feed.Stage >= col.Stages {
			v.add(path+".stage", "must be in [2, N-1]")
		}
	}
	return v.result("compiler input")
}

func (v *validator) properties(path string, p Properties) {
	switch p.Source {
	case "builtin":
		if p.Key == "" {
			v.add(path+".key", "required field is missing")
		}
		if p.Revision == "" {
			v.add(path+".revision", "required field is missing")
		}
		if p.Antoine != nil || p.Enthalpy != nil || p.TemperatureRangeK != nil {
			v.add(path, "builtin properties cannot contain inline parameters")
		}
		return
	case "inline":
		if p.Key != "" || p.Revision != "" {
			v.add(path, "inline properties cannot contain builtin key or revision")
		}
	default:
		v.add(path+".source", "must be builtin or inline")
		return
	}
	if p.Antoine == nil {
		v.add(path+".antoine", "required field is missing")
	}
	if p.Enthalpy == nil {
		v.add(path+".enthalpy", "required field is missing")
	}
	if p.TemperatureRangeK == nil {
		v.add(path+".temperature_range_K", "required field is missing")
	}
	a, h, r := p.Antoine, p.Enthalpy, p.TemperatureRangeK
	validRange := false
	if r != nil {
		low := v.positive(path+".temperature_range_K.min", r.Min)
		high := v.positive(path+".temperature_range_K.max", r.Max)
		if low && high {
			if r.Min >= r.Max {
				v.add(path+".temperature_range_K", "require min < max")
			} else {
				validRange = true
			}
		}
	}
	if a != nil {
		for _, item := range []struct {
			name  string
			value float64
		}{{"A", a.A}, {"B", a.B}, {"C", a.C}} {
			if !finite(item.value) {
				v.add(path+".antoine."+item.name, "must be finite")
			}
		}
		v.equal(path+".antoine.pressure_unit", a.PressureUnit, "mmHg")
		v.equal(path+".antoine.temperature_unit", a.TemperatureUnit, "C")
		pole := 273.15 - a.C
		if validRange && finite(a.C) && pole >= r.Min && pole <= r.Max {
			v.add(path+".antoine", "pole inside temperature range")
		}
	}
	if h != nil {
		v.equal(path+".enthalpy.type", h.Type, "constant_cp")
		tref := v.positive(path+".enthalpy.reference_temperature_K", h.ReferenceTemperatureK)
		cl := v.positive(path+".enthalpy.liquid_cp_J_mol_K", h.LiquidCPJMolK)
		cv := v.positive(path+".enthalpy.vapor_cp_J_mol_K", h.VaporCPJMolK)
		latent := v.positive(path+".enthalpy.vaporization_enthalpy_at_reference_J_mol", h.VaporizationEnthalpyAtReferenceJMol)
		if validRange && tref {
			if h.ReferenceTemperatureK < r.Min || h.ReferenceTemperatureK > r.Max {
				v.add(path+".enthalpy.reference_temperature_K", "must be inside temperature range")
			}
			if cl && cv && latent {
				for _, t := range []float64{r.Min, r.Max} {
					delta := h.VaporizationEnthalpyAtReferenceJMol + (h.VaporCPJMolK-h.LiquidCPJMolK)*(t-h.ReferenceTemperatureK)
					if !finite(delta) || delta <= 0 {
						v.add(path+".enthalpy", "vaporization enthalpy must be finite and positive throughout temperature_range_K")
						break
					}
				}
			}
		}
	}
}
