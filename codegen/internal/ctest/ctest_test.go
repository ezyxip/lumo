package ctest

import (
	"os"
	"testing"
)

func TestDiagnostics(t *testing.T) {
	harness, err := os.ReadFile("testdata/check.c")
	if err != nil {
		t.Fatal(err)
	}
	program := Build(t, map[string][]byte{"check.c": harness})
	t.Run("success", func(t *testing.T) {
		program.Run(t, "success")
	})
	for _, tc := range []struct {
		name string
		want []string
	}{
		{"near_failure", []string{"near_failure", "check.c:", "actual=2", "expected=1", "error=1", "tolerance="}},
		{"status_failure", []string{"status_failure", "check.c:", "actual=3", "expected=0"}},
		{"nan_failure", []string{"nan_failure", "check.c:", "actual=nan", "expected=1"}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			program.RunFailure(t, tc.name, tc.want...)
		})
	}
}
