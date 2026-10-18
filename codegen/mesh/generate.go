// Package mesh generates residuals and numerical Jacobians of a steady-state ideal column.
package mesh

import (
	"bytes"
	_ "embed"
	"fmt"
	"math"
	"regexp"
	"strconv"
	"text/template"

	"github.com/ezyxip/lumo/codegen/jacobian"
	"github.com/ezyxip/lumo/codegen/thermo"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

//go:embed templates/mesh.h.tmpl
var headerTemplate string

//go:embed templates/mesh.c.tmpl
var sourceTemplate string

//go:embed templates/jacobian.c.tmpl
var jacobianTemplate string

// Files contains the complete C99 residual and Jacobian module and its thermodynamics.
// Write Header/Source as lumo_mesh.h/.c, Thermo as lumo_thermo.h/.c,
// and Jacobian as lumo_jacobian.h/.c.
// FeedIDs and Thermo.ComponentIDs define the order of C input arrays.
type Files struct {
	Header    []byte
	Source    []byte
	Thermo    thermo.Files
	Jacobian  jacobian.Files
	FeedIDs   []string
	Variables int
}

type feedLocation struct {
	Index int
	Stage int // zero-based C index
}

type templateData struct {
	Stages     int
	Components int
	Variables  int
	Pressure   string
	Feeds      []feedLocation
}

var feedID = regexp.MustCompile(`^[A-Za-z][A-Za-z0-9_-]*$`)

// Generate supports total condensers, equilibrium reboilers, constant pressure
// and a molar reflux-ratio closure. Runtime operating parameters are normalized
// C inputs, not raw model_input.json. Properties must be resolved inline data.
func Generate(cfg contractv1.CompilerInput) (Files, error) {
	col := cfg.Column
	if col.Stages < 3 {
		return Files{}, fmt.Errorf("mesh: column.stages must be >= 3")
	}
	if col.Condenser.Type != "total" || col.Reboiler.Type != "equilibrium" {
		return Files{}, fmt.Errorf("mesh: only total condenser and equilibrium reboiler are supported")
	}
	if col.Pressure.Type != "constant" || math.IsNaN(col.Pressure.ValuePa) || math.IsInf(col.Pressure.ValuePa, 0) || col.Pressure.ValuePa <= 0 {
		return Files{}, fmt.Errorf("mesh: require constant finite positive pressure in Pa")
	}
	if len(col.Feeds) == 0 {
		return Files{}, fmt.Errorf("mesh: at least one feed is required")
	}
	feeds := make([]feedLocation, len(col.Feeds))
	ids := make([]string, len(col.Feeds))
	seen := make(map[string]bool, len(col.Feeds))
	for i, f := range col.Feeds {
		if !feedID.MatchString(f.ID) || seen[f.ID] {
			return Files{}, fmt.Errorf("mesh: column.feeds[%d].id is invalid or duplicate", i)
		}
		seen[f.ID] = true
		if f.Stage < 2 || f.Stage >= col.Stages {
			return Files{}, fmt.Errorf("mesh: column.feeds[%d].stage must be in [2, N-1]", i)
		}
		feeds[i] = feedLocation{Index: i, Stage: f.Stage - 1}
		ids[i] = f.ID
	}
	// Keep dimensions representable in the C integer constant expressions.
	if len(cfg.Components) > (2147483647-3)/2 {
		return Files{}, fmt.Errorf("mesh: too many components")
	}
	block := 2*len(cfg.Components) + 3
	if col.Stages > (2147483647-1)/block {
		return Files{}, fmt.Errorf("mesh: system dimension exceeds C int range")
	}
	properties, err := thermo.Generate(cfg)
	if err != nil {
		return Files{}, fmt.Errorf("mesh: %w", err)
	}
	data := templateData{
		Stages: col.Stages, Components: len(cfg.Components), Variables: col.Stages*block + 1,
		Pressure: strconv.FormatFloat(col.Pressure.ValuePa, 'g', -1, 64), Feeds: feeds,
	}
	derivatives, err := jacobian.Generate(data.Variables)
	if err != nil {
		return Files{}, fmt.Errorf("mesh: %w", err)
	}
	header, err := render("lumo_mesh.h", headerTemplate, data)
	if err != nil {
		return Files{}, err
	}
	source, err := render("lumo_mesh.c", sourceTemplate+jacobianTemplate, data)
	if err != nil {
		return Files{}, err
	}
	return Files{Header: header, Source: source, Thermo: properties, Jacobian: derivatives, FeedIDs: ids, Variables: data.Variables}, nil
}

func render(name, source string, data templateData) ([]byte, error) {
	t, err := template.New(name).Option("missingkey=error").Parse(source)
	if err != nil {
		return nil, fmt.Errorf("mesh template %s: %w", name, err)
	}
	var out bytes.Buffer
	if err := t.Execute(&out, data); err != nil {
		return nil, fmt.Errorf("mesh template %s: %w", name, err)
	}
	return out.Bytes(), nil
}
