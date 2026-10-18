// Package jacobian generates model-independent C99 Jacobian evaluators.
package jacobian

import (
	"bytes"
	_ "embed"
	"fmt"
	"text/template"
)

//go:embed templates/jacobian.h.tmpl
var headerTemplate string

//go:embed templates/jacobian.c.tmpl
var sourceTemplate string

// Files contains lumo_jacobian.h and lumo_jacobian.c.
type Files struct {
	Header []byte
	Source []byte
}

// Generate specializes scratch arrays for a square system of the given size.
// The generated code accepts a residual callback and has no model dependencies.
func Generate(size int) (Files, error) {
	if size < 1 || size > 2147483647 {
		return Files{}, fmt.Errorf("jacobian: size must be in [1, 2147483647]")
	}
	render := func(name, source string) ([]byte, error) {
		t, err := template.New(name).Option("missingkey=error").Parse(source)
		if err != nil {
			return nil, fmt.Errorf("jacobian template %s: %w", name, err)
		}
		var out bytes.Buffer
		if err := t.Execute(&out, size); err != nil {
			return nil, fmt.Errorf("jacobian template %s: %w", name, err)
		}
		return out.Bytes(), nil
	}
	header, err := render("lumo_jacobian.h", headerTemplate)
	if err != nil {
		return Files{}, err
	}
	source, err := render("lumo_jacobian.c", sourceTemplate)
	if err != nil {
		return Files{}, err
	}
	return Files{Header: header, Source: source}, nil
}
