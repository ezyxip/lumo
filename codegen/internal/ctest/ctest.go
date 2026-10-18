// Package ctest builds and runs C99 checks for code generators using cc.
package ctest

import (
	_ "embed"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

//go:embed ctest.h
var assertions []byte

// Program is a compiled C test executable and its debug artifacts.
type Program struct {
	Dir    string
	Binary string
}

// Build writes the supplied files and the shared ctest.h, then compiles every
// .c file using cc. File names must be local base names; ctest.h is reserved.
// Artifacts survive failures, or successful runs with LUMO_CODEGEN_KEEP=1.
func Build(t *testing.T, files map[string][]byte) *Program {
	t.Helper()
	cc, err := exec.LookPath("cc")
	if err != nil {
		t.Fatal("C99 integration tests require cc in PATH")
	}
	dir, err := os.MkdirTemp("", "lumo-codegen-*")
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() {
		if t.Failed() || os.Getenv("LUMO_CODEGEN_KEEP") == "1" {
			t.Logf("C debug artifacts kept in %s", dir)
			return
		}
		if err := os.RemoveAll(dir); err != nil {
			t.Errorf("remove C debug artifacts: %v", err)
		}
	})
	write(t, dir, "ctest.h", assertions, 0600)
	names := make([]string, 0, len(files))
	for name := range files {
		if !baseName(name) || name == "ctest.h" || name == "check" || name == "build.sh" || name == "build.log" {
			t.Fatalf("invalid or reserved C test file name %q", name)
		}
		names = append(names, name)
	}
	sort.Strings(names)
	args := []string{"-std=c99", "-g", "-O0", "-Wall", "-Wextra", "-Werror", "-pedantic"}
	for _, name := range names {
		write(t, dir, name, files[name], 0600)
		if strings.HasSuffix(name, ".c") {
			args = append(args, filepath.Join(dir, name))
		}
	}
	bin := filepath.Join(dir, "check")
	args = append(args, "-lm", "-o", bin)
	cmd := exec.Command(cc, args...)
	quoted := make([]string, len(cmd.Args))
	for i, arg := range cmd.Args {
		quoted[i] = "'" + strings.ReplaceAll(arg, "'", "'\"'\"'") + "'"
	}
	write(t, dir, "build.sh", []byte("#!/bin/sh\nset -eu\nexec "+strings.Join(quoted, " ")+"\n"), 0700)
	out, err := cmd.CombinedOutput()
	write(t, dir, "build.log", out, 0600)
	if err != nil {
		t.Fatalf("compile generated C99: %v\nartifacts: %s\n%s", err, dir, out)
	}
	return &Program{Dir: dir, Binary: bin}
}

// Run executes one named C test and reports its output on failure.
func (p *Program) Run(t *testing.T, name string) {
	t.Helper()
	out, err := p.run(t, name)
	if err != nil {
		t.Fatalf("C test %s: %v\nartifacts: %s\n%s", name, err, p.Dir, out)
	}
}

// RunFailure verifies an intentional assertion failure and its diagnostics.
func (p *Program) RunFailure(t *testing.T, name string, want ...string) {
	t.Helper()
	out, err := p.run(t, name)
	var exitErr *exec.ExitError
	if !errors.As(err, &exitErr) || exitErr.ExitCode() != 1 {
		t.Fatalf("C test %s: expected exit status 1, got %v\nartifacts: %s\n%s", name, err, p.Dir, out)
	}
	for _, fragment := range want {
		if !strings.Contains(string(out), fragment) {
			t.Errorf("C test %s: missing diagnostic %q\nartifacts: %s\n%s", name, fragment, p.Dir, out)
		}
	}
}

func (p *Program) run(t *testing.T, name string) ([]byte, error) {
	t.Helper()
	if !baseName(name) {
		t.Fatalf("invalid C test name %q", name)
	}
	out, err := exec.Command(p.Binary, name).CombinedOutput()
	write(t, p.Dir, name+".log", out, 0600)
	return out, err
}

func baseName(name string) bool {
	return name != "" && name != "." && name != ".." && filepath.Base(name) == name
}

func write(t *testing.T, dir, name string, data []byte, mode os.FileMode) {
	t.Helper()
	if err := os.WriteFile(filepath.Join(dir, name), data, mode); err != nil {
		t.Fatal(err)
	}
}
