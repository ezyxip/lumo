// Package initialguess generates initial MESH states for ideal columns.
package initialguess

import (
	"bytes"
	_ "embed"
	"fmt"
	"text/template"

	"github.com/ezyxip/lumo/codegen/mesh"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

//go:embed templates/initialguess.h.tmpl
var headerTemplate string

//go:embed templates/initialguess.c.tmpl
var sourceTemplate string

// Files contains lumo_initialguess.h/.c and their generated MESH dependencies.
// Component and feed array orders are Mesh.Thermo.ComponentIDs and Mesh.FeedIDs.
type Files struct {
	Header []byte
	Source []byte
	Mesh   mesh.Files
}

// Generate builds the linear_T_x initialization strategy for the column.
func Generate(cfg contractv1.CompilerInput) (Files, error) {
	if err := contractv1.ValidateCompilerInput(cfg); err != nil {
		return Files{}, fmt.Errorf("initialguess: %w", err)
	}
	model, err := mesh.Generate(cfg)
	if err != nil {
		return Files{}, fmt.Errorf("initialguess: %w", err)
	}
	render := func(name, source string) ([]byte, error) {
		t, err := template.New(name).Option("missingkey=error").Parse(source)
		if err != nil {
			return nil, fmt.Errorf("initialguess template %s: %w", name, err)
		}
		var out bytes.Buffer
		if err := t.Execute(&out, cfg); err != nil {
			return nil, fmt.Errorf("initialguess template %s: %w", name, err)
		}
		return out.Bytes(), nil
	}
	header, err := render("lumo_initialguess.h", headerTemplate)
	if err != nil {
		return Files{}, err
	}
	source, err := render("lumo_initialguess.c", sourceTemplate)
	if err != nil {
		return Files{}, err
	}
	return Files{Header: header, Source: source, Mesh: model}, nil
}
