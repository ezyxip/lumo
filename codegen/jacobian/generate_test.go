package jacobian

import (
	"bytes"
	"fmt"
	"os"
	"testing"

	"github.com/ezyxip/lumo/codegen/internal/ctest"
)

func TestGenerateAndRunC99(t *testing.T) {
	harness, err := os.ReadFile("testdata/check.c")
	if err != nil {
		t.Fatal(err)
	}
	for _, size := range []int{1, 3} {
		files, err := Generate(size)
		if err != nil {
			t.Fatal(err)
		}
		again, err := Generate(size)
		if err != nil {
			t.Fatal(err)
		}
		if !bytes.Equal(files.Header, again.Header) || !bytes.Equal(files.Source, again.Source) {
			t.Fatal("generation is not deterministic")
		}
		program := ctest.Build(t, map[string][]byte{
			"lumo_jacobian.h": files.Header, "lumo_jacobian.c": files.Source,
			"check.c": harness,
		})
		for _, name := range []string{"polynomial", "domain", "errors"} {
			t.Run(fmt.Sprintf("%d/%s", size, name), func(t *testing.T) { program.Run(t, name) })
		}
	}
}

func TestInvalidSize(t *testing.T) {
	for _, size := range []int{-1, 0, 2147483648} {
		files, err := Generate(size)
		if err == nil || len(files.Header) != 0 || len(files.Source) != 0 {
			t.Fatalf("size %d produced code or no error", size)
		}
	}
}
