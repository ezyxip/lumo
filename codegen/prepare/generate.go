// Package prepare generates C99 conversion of runtime feeds into MESH inputs.
package prepare

import (
	"bytes"
	_ "embed"
	"fmt"
	"text/template"

	"github.com/ezyxip/lumo/codegen/mesh"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

//go:embed templates/prepare.h.tmpl
var headerTemplate string

//go:embed templates/prepare.c.tmpl
var sourceTemplate string

// Files contains lumo_prepare.h/.c and their generated MESH dependencies.
// Component and feed array orders are Mesh.Thermo.ComponentIDs and Mesh.FeedIDs.
type Files struct {
	Header []byte
	Source []byte
	Mesh   mesh.Files
}

// Generate compiles the description, not runtime feed values. Builtin properties
// must already be resolved. Mass distillate specifications require a different
// MESH closure and are reported as unsupported by the C preparation function.
func Generate(cfg contractv1.CompilerInput) (Files, error) {
	if err := contractv1.ValidateCompilerInput(cfg); err != nil {
		return Files{}, fmt.Errorf("prepare: %w", err)
	}
	model, err := mesh.Generate(cfg)
	if err != nil {
		return Files{}, fmt.Errorf("prepare: %w", err)
	}
	render := func(name, source string) ([]byte, error) {
		t, err := template.New(name).Option("missingkey=error").Parse(source)
		if err != nil {
			return nil, fmt.Errorf("prepare template %s: %w", name, err)
		}
		var out bytes.Buffer
		if err := t.Execute(&out, cfg); err != nil {
			return nil, fmt.Errorf("prepare template %s: %w", name, err)
		}
		return out.Bytes(), nil
	}
	header, err := render("lumo_prepare.h", headerTemplate)
	if err != nil {
		return Files{}, err
	}
	source, err := render("lumo_prepare.c", sourceTemplate)
	if err != nil {
		return Files{}, err
	}
	return Files{Header: header, Source: source, Mesh: model}, nil
}
