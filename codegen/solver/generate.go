// Package solver generates a damped Newton solver for MESH columns.
package solver

import (
	"bytes"
	_ "embed"
	"fmt"
	"text/template"

	"github.com/ezyxip/lumo/codegen/mesh"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

//go:embed templates/solver.h.tmpl
var headerTemplate string

//go:embed templates/solver.c.tmpl
var sourceTemplate string

// Files contains lumo_solver.h/.c and their generated MESH dependencies.
// Component and feed array orders are Mesh.Thermo.ComponentIDs and Mesh.FeedIDs.
type Files struct {
	Header []byte
	Source []byte
	Mesh   mesh.Files
}

// Generate specializes the dense Newton solver for the column dimensions.
func Generate(cfg contractv1.CompilerInput) (Files, error) {
	if err := contractv1.ValidateCompilerInput(cfg); err != nil {
		return Files{}, fmt.Errorf("solver: %w", err)
	}
	model, err := mesh.Generate(cfg)
	if err != nil {
		return Files{}, fmt.Errorf("solver: %w", err)
	}
	render := func(name, source string) ([]byte, error) {
		t, err := template.New(name).Option("missingkey=error").Parse(source)
		if err != nil {
			return nil, fmt.Errorf("solver template %s: %w", name, err)
		}
		var out bytes.Buffer
		if err := t.Execute(&out, cfg); err != nil {
			return nil, fmt.Errorf("solver template %s: %w", name, err)
		}
		return out.Bytes(), nil
	}
	header, err := render("lumo_solver.h", headerTemplate)
	if err != nil {
		return Files{}, err
	}
	source, err := render("lumo_solver.c", sourceTemplate)
	if err != nil {
		return Files{}, err
	}
	return Files{Header: header, Source: source, Mesh: model}, nil
}
