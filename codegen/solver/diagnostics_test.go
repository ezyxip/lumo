package solver

import (
	"bytes"
	"os"
	"strings"
	"testing"

	"github.com/ezyxip/lumo/codegen/internal/ctest"
)

func TestDiagnosticsModes(t *testing.T) {
	h, err := os.ReadFile("testdata/observation.c")
	if err != nil {
		t.Fatal(err)
	}
	for _, enabled := range []bool{true, false} {
		name := "disabled"
		if enabled {
			name = "enabled"
		}
		t.Run(name, func(t *testing.T) {
			cfg := config(t)
			f, err := GenerateWithOptions(cfg, GenerationOptions{Diagnostics: enabled})
			if err != nil {
				t.Fatal(err)
			}
			if !enabled && (bytes.Contains(f.Source, []byte("notify_observer")) || bytes.Contains(f.Source, []byte("options->observer")) || bytes.Contains(f.Source, []byte("worst_residual"))) {
				t.Fatal("diagnostic code remains in disabled source")
			}
			again, err := GenerateWithOptions(cfg, GenerationOptions{Diagnostics: enabled})
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(f.Source, again.Source) || !bytes.Equal(f.Header, again.Header) {
				t.Fatal("nondeterministic generation")
			}
			files := sourceFiles(t, cfg, h)
			files["lumo_solver.h"] = f.Header
			files["lumo_solver.c"] = f.Source
			p := ctest.Build(t, files)
			for _, name := range []string{"observation", "observation_failures", "backtracking_events"} {
				t.Run(name, func(t *testing.T) { p.Run(t, name) })
			}
		})
	}
}
func TestDiagnosticsDoNotChangeSolution(t *testing.T) {
	cfg := config(t)
	off, err := GenerateWithOptions(cfg, GenerationOptions{})
	if err != nil {
		t.Fatal(err)
	}
	h, err := os.ReadFile("testdata/diagnostics_equivalence.c")
	if err != nil {
		t.Fatal(err)
	}
	files := sourceFiles(t, cfg, h)
	renamed := strings.NewReplacer("lumo_solver_solve", "off_solve", "lumo_solver_default_options", "off_default_options", "lumo_solver_workspace_doubles", "off_workspace_doubles")
	files["solver_off.c"] = []byte(renamed.Replace(string(off.Source)))
	p := ctest.Build(t, files)
	p.Run(t, "equivalence")
}
