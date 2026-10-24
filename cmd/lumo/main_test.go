package main

import (
	"bytes"
	"encoding/json"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

func inputFile(t *testing.T) string {
	t.Helper()
	data, err := os.ReadFile("../../codegen/mesh/testdata/nonlinear/case.json")
	if err != nil {
		t.Fatal(err)
	}
	var fixture struct {
		Input json.RawMessage `json:"compiler_input"`
	}
	if err := json.Unmarshal(data, &fixture); err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(t.TempDir(), "input.json")
	if err := os.WriteFile(path, fixture.Input, 0600); err != nil {
		t.Fatal(err)
	}
	return path
}

func TestGenerateCLI(t *testing.T) {
	for _, diagnostics := range []bool{false, true} {
		t.Run(map[bool]string{false: "off", true: "on"}[diagnostics], func(t *testing.T) {
			output := filepath.Join(t.TempDir(), "model")
			args := []string{"-o", output}
			if diagnostics {
				args = append(args, "-diagnostics")
			}
			args = append(args, inputFile(t))
			if err := run(args, &bytes.Buffer{}); err != nil {
				t.Fatal(err)
			}
			entries, err := os.ReadDir(output)
			if err != nil || len(entries) != 12 {
				t.Fatalf("files=%d err=%v", len(entries), err)
			}
			header, err := os.ReadFile(filepath.Join(output, "lumo_solver.h"))
			if err != nil {
				t.Fatal(err)
			}
			macro := "#define LUMO_SOLVER_HAS_DIAGNOSTICS 0"
			if diagnostics {
				macro = "#define LUMO_SOLVER_HAS_DIAGNOSTICS 1"
			}
			if !bytes.Contains(header, []byte(macro)) {
				t.Fatal("wrong diagnostics mode")
			}
			if err := run(args, &bytes.Buffer{}); err == nil {
				t.Fatal("overwrote existing output")
			}
			cc, err := exec.LookPath("cc")
			if err != nil {
				t.Skip("cc unavailable")
			}
			sources, err := filepath.Glob(filepath.Join(output, "*.c"))
			if err != nil {
				t.Fatal(err)
			}
			// Linking a shared library verifies all generated modules and dependencies together.
			ccArgs := append([]string{"-std=c99", "-Wall", "-Wextra", "-Werror", "-pedantic", "-shared", "-fPIC"}, sources...)
			ccArgs = append(ccArgs, "-lm", "-o", filepath.Join(t.TempDir(), "model.so"))
			if out, err := exec.Command(cc, ccArgs...).CombinedOutput(); err != nil {
				t.Fatalf("cc: %v\n%s", err, out)
			}
		})
	}
}

func TestErrors(t *testing.T) {
	input := inputFile(t)
	for _, args := range [][]string{nil, {"-unknown"}, {input, input}, {"missing.json"}} {
		if err := run(args, &bytes.Buffer{}); err == nil {
			t.Fatalf("accepted %v", args)
		}
	}
	for _, data := range []string{"{", `{"schema_version":1}`, strings.ReplaceAll(readInput(t, input), `"stages": 3`, `"stages": 1`)} {
		root := t.TempDir()
		bad := filepath.Join(root, "bad.json")
		if err := os.WriteFile(bad, []byte(data), 0600); err != nil {
			t.Fatal(err)
		}
		output := filepath.Join(root, "model")
		if err := run([]string{"-o", output, bad}, &bytes.Buffer{}); err == nil {
			t.Fatal("accepted bad input")
		}
		if _, err := os.Stat(output); !os.IsNotExist(err) {
			t.Fatalf("partial output: %v", err)
		}
	}
	var help bytes.Buffer
	if err := run([]string{"-help"}, &help); err != nil || !strings.Contains(help.String(), "Usage:") {
		t.Fatalf("help: %v %s", err, help.String())
	}
}

func readInput(t *testing.T, path string) string {
	t.Helper()
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	return string(data)
}
