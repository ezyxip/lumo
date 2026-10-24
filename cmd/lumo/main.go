// Command lumo compiles a v1 column description into C99 sources.
package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"sort"

	"github.com/ezyxip/lumo/codegen/initialguess"
	"github.com/ezyxip/lumo/codegen/prepare"
	"github.com/ezyxip/lumo/codegen/solver"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

func main() {
	if err := run(os.Args[1:], os.Stderr); err != nil {
		fmt.Fprintln(os.Stderr, "lumo:", err)
		os.Exit(1)
	}
}

func run(args []string, stderr io.Writer) error {
	flags := flag.NewFlagSet("lumo", flag.ContinueOnError)
	flags.SetOutput(stderr)
	output := flags.String("o", "model", "output directory (must not exist)")
	diagnostics := flags.Bool("diagnostics", false, "generate solver diagnostics callbacks")
	flags.Usage = func() {
		fmt.Fprintln(stderr, "Usage: lumo [-diagnostics] [-o directory] compiler_input.json")
		flags.PrintDefaults()
	}
	if err := flags.Parse(args); err != nil {
		if errors.Is(err, flag.ErrHelp) {
			return nil
		}
		return err
	}
	if flags.NArg() != 1 {
		return fmt.Errorf("expected one compiler_input.json path; use -help for usage")
	}
	input, err := os.Open(flags.Arg(0))
	if err != nil {
		return fmt.Errorf("open input: %w", err)
	}
	cfg, err := contractv1.ParseCompilerInput(input)
	closeErr := input.Close()
	if err != nil {
		return err
	}
	if closeErr != nil {
		return closeErr
	}
	if err := contractv1.ValidateCompilerInput(cfg); err != nil {
		return err
	}
	files, err := generate(cfg, *diagnostics)
	if err != nil {
		return err
	}
	return writeFiles(*output, files)
}

func generate(cfg contractv1.CompilerInput, diagnostics bool) (map[string][]byte, error) {
	s, err := solver.GenerateWithOptions(cfg, solver.GenerationOptions{Diagnostics: diagnostics})
	if err != nil {
		return nil, err
	}
	p, err := prepare.Generate(cfg)
	if err != nil {
		return nil, err
	}
	g, err := initialguess.Generate(cfg)
	if err != nil {
		return nil, err
	}
	return map[string][]byte{
		"lumo_solver.h": s.Header, "lumo_solver.c": s.Source,
		"lumo_mesh.h": s.Mesh.Header, "lumo_mesh.c": s.Mesh.Source,
		"lumo_thermo.h": s.Mesh.Thermo.Header, "lumo_thermo.c": s.Mesh.Thermo.Source,
		"lumo_jacobian.h": s.Mesh.Jacobian.Header, "lumo_jacobian.c": s.Mesh.Jacobian.Source,
		"lumo_prepare.h": p.Header, "lumo_prepare.c": p.Source,
		"lumo_initialguess.h": g.Header, "lumo_initialguess.c": g.Source,
	}, nil
}

// Stage files beside the destination so generation/write errors leave no partial output.
func writeFiles(output string, files map[string][]byte) error {
	output = filepath.Clean(output)
	if _, err := os.Lstat(output); err == nil {
		return fmt.Errorf("output path already exists: %s", output)
	} else if !os.IsNotExist(err) {
		return err
	}
	staging, err := os.MkdirTemp(filepath.Dir(output), ".lumo-*")
	if err != nil {
		return fmt.Errorf("create output directory: %w", err)
	}
	defer os.RemoveAll(staging)
	names := make([]string, 0, len(files))
	for name := range files {
		names = append(names, name)
	}
	sort.Strings(names)
	for _, name := range names {
		if err := os.WriteFile(filepath.Join(staging, name), files[name], 0644); err != nil {
			return err
		}
	}
	// Recheck before publishing, including dangling symlinks.
	if _, err := os.Lstat(output); err == nil {
		return fmt.Errorf("output path already exists: %s", output)
	} else if !os.IsNotExist(err) {
		return err
	}
	if err := os.Rename(staging, output); err != nil {
		return fmt.Errorf("publish output: %w", err)
	}
	return nil
}
