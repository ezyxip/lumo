// Package thermo generates the C99 thermodynamic module used by MESH.
package thermo

import (
	"bytes"
	_ "embed"
	"fmt"
	"math"
	"regexp"
	"strconv"
	"text/template"

	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

//go:embed templates/thermo.h.tmpl
var headerTemplate string

//go:embed templates/antoine_raoult.c.tmpl
var sourceTemplate string

// Files contains a standalone module. Write Header as lumo_thermo.h and Source
// as lumo_thermo.c in the same directory. Link with the standard math library.
// Array indices in the C interface follow ComponentIDs, in input order.
type Files struct {
	Header       []byte
	Source       []byte
	ComponentIDs []string
}

// Generate selects a thermodynamic implementation and embeds its parameters.
// Currently only inline Antoine/Raoult + ideal constant-Cp enthalpies are supported.
// It validates thermodynamic inputs, not column topology or operating parameters.
// Builtin property references must be resolved to inline properties beforehand.
func Generate(cfg contractv1.CompilerInput) (Files, error) {
	if cfg.SchemaVersion != contractv1.SchemaVersion {
		return Files{}, fmt.Errorf("thermo: unsupported schema_version %d", cfg.SchemaVersion)
	}
	if cfg.Thermo.PhaseEquilibrium != "antoine_raoult" || cfg.Thermo.Enthalpy != "ideal" {
		return Files{}, fmt.Errorf("thermo: unsupported model %q / %q", cfg.Thermo.PhaseEquilibrium, cfg.Thermo.Enthalpy)
	}
	if err := validateComponents(cfg.Components); err != nil {
		return Files{}, err
	}
	data := struct{ Components []contractv1.Component }{cfg.Components}
	funcs := template.FuncMap{"literal": func(v float64) string {
		return strconv.FormatFloat(v, 'g', -1, 64)
	}}
	header, err := render("thermo.h", headerTemplate, funcs, data)
	if err != nil {
		return Files{}, err
	}
	source, err := render("thermo.c", sourceTemplate, funcs, data)
	if err != nil {
		return Files{}, err
	}
	ids := make([]string, len(cfg.Components))
	for i, c := range cfg.Components {
		ids[i] = c.ID
	}
	return Files{Header: header, Source: source, ComponentIDs: ids}, nil
}

func render(name, source string, funcs template.FuncMap, data any) ([]byte, error) {
	t, err := template.New(name).Funcs(funcs).Option("missingkey=error").Parse(source)
	if err != nil {
		return nil, fmt.Errorf("thermo template %s: %w", name, err)
	}
	var out bytes.Buffer
	if err := t.Execute(&out, data); err != nil {
		return nil, fmt.Errorf("thermo template %s: %w", name, err)
	}
	return out.Bytes(), nil
}

var componentID = regexp.MustCompile(`^[A-Za-z][A-Za-z0-9_-]*$`)

func finite(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }

func validateComponents(components []contractv1.Component) error {
	if len(components) < 2 {
		return fmt.Errorf("thermo: need at least two components")
	}
	seen := make(map[string]bool, len(components))
	commonMin, commonMax := 0.0, math.Inf(1)
	for i, c := range components {
		path := fmt.Sprintf("thermo: components[%d]", i)
		if !componentID.MatchString(c.ID) || seen[c.ID] {
			return fmt.Errorf("%s.id: invalid or duplicate ID %q", path, c.ID)
		}
		seen[c.ID] = true
		p := c.Properties
		if p.Source != "inline" {
			return fmt.Errorf("%s.properties.source: expected resolved inline properties, got %q", path, p.Source)
		}
		if p.Key != "" || p.Revision != "" {
			return fmt.Errorf("%s.properties: inline properties cannot contain builtin key or revision", path)
		}
		if p.Antoine == nil || p.Enthalpy == nil || p.TemperatureRangeK == nil {
			return fmt.Errorf("%s.properties: antoine, enthalpy and temperature_range_K are required", path)
		}
		a, h, r := p.Antoine, p.Enthalpy, p.TemperatureRangeK
		for _, v := range []float64{a.A, a.B, a.C, h.ReferenceTemperatureK, h.LiquidCPJMolK,
			h.VaporCPJMolK, h.VaporizationEnthalpyAtReferenceJMol, r.Min, r.Max} {
			if !finite(v) {
				return fmt.Errorf("%s.properties: parameters must be finite", path)
			}
		}
		if a.PressureUnit != "mmHg" || a.TemperatureUnit != "C" {
			return fmt.Errorf("%s.antoine: expected pressure_unit=mmHg and temperature_unit=C", path)
		}
		if h.Type != "constant_cp" {
			return fmt.Errorf("%s.enthalpy.type: unsupported model %q", path, h.Type)
		}
		if h.LiquidCPJMolK <= 0 || h.VaporCPJMolK <= 0 || h.VaporizationEnthalpyAtReferenceJMol <= 0 {
			return fmt.Errorf("%s.enthalpy: heat capacities and reference vaporization enthalpy must be positive", path)
		}
		if r.Min <= 0 || r.Min >= r.Max || h.ReferenceTemperatureK < r.Min || h.ReferenceTemperatureK > r.Max {
			return fmt.Errorf("%s.temperature_range_K: require 0 < min < max and reference temperature inside range", path)
		}
		// Do not require positive Antoine coefficients: their signs depend on
		// the fitted correlation. Only its pole must be outside the range.
		pole := 273.15 - a.C
		if pole >= r.Min && pole <= r.Max {
			return fmt.Errorf("%s.antoine: pole inside temperature range", path)
		}
		commonMin = math.Max(commonMin, r.Min)
		commonMax = math.Min(commonMax, r.Max)
	}
	if commonMin >= commonMax {
		return fmt.Errorf("thermo: component temperature ranges have no common interval")
	}
	return nil
}
