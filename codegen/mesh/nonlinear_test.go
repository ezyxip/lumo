package mesh

import (
	"bytes"
	"crypto/sha256"
	"encoding/json"
	"fmt"
	"os"
	"testing"

	"github.com/ezyxip/lumo/codegen/internal/ctest"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

func TestNonlinearReference(t *testing.T) {
	input, err := os.ReadFile("testdata/nonlinear/case.json")
	if err != nil {
		t.Fatal(err)
	}
	var testcase struct {
		CompilerInput json.RawMessage `json:"compiler_input"`
	}
	if err := json.Unmarshal(input, &testcase); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile("testdata/nonlinear/reference.json")
	if err != nil {
		t.Fatal(err)
	}
	var reference struct {
		InputSHA256 string    `json:"input_sha256"`
		State       []float64 `json:"state"`
		Operating   struct {
			Flow float64   `json:"feed_flow_mol_s"`
			Z    []float64 `json:"feed_mole_fractions"`
			H    float64   `json:"feed_enthalpy_J_mol"`
			D    float64   `json:"distillate_mol_s"`
			R    float64   `json:"reflux_ratio"`
		} `json:"operating"`
	}
	if err := json.Unmarshal(data, &reference); err != nil {
		t.Fatal(err)
	}
	if fmt.Sprintf("%x", sha256.Sum256(input)) != reference.InputSHA256 {
		t.Fatal("stale reference: run python codegen/mesh/testdata/nonlinear/reference.py")
	}
	cfg, err := contractv1.ParseCompilerInput(bytes.NewReader(testcase.CompilerInput))
	if err != nil {
		t.Fatal(err)
	}
	if err := contractv1.ValidateCompilerInput(cfg); err != nil {
		t.Fatal(err)
	}
	files, err := Generate(cfg)
	if err != nil {
		t.Fatal(err)
	}
	if len(reference.State) != files.Variables || files.Variables != 22 || len(reference.Operating.Z) != 2 {
		t.Fatal("nonlinear fixture must describe a binary three-stage column")
	}
	var header bytes.Buffer
	header.WriteString("/* Reference generated independently by reference.py. */\nstatic const double reference_state[LUMO_MESH_SIZE] = {\n")
	for _, v := range reference.State {
		fmt.Fprintf(&header, "%.17g,\n", v)
	}
	header.WriteString("};\nstatic const struct lumo_mesh_operating reference_operating = {\n")
	op := reference.Operating
	fmt.Fprintf(&header, ".feeds = {{.flow_mol_s = %.17g, .mole_fractions = {%.17g, %.17g}, .enthalpy_J_mol = %.17g}},\n.distillate_mol_s = %.17g, .reflux_ratio = %.17g\n};\n", op.Flow, op.Z[0], op.Z[1], op.H, op.D, op.R)
	harness, err := os.ReadFile("testdata/nonlinear/check.c")
	if err != nil {
		t.Fatal(err)
	}
	program := ctest.Build(t, map[string][]byte{
		"lumo_mesh.h": files.Header, "lumo_mesh.c": files.Source,
		"lumo_jacobian.h": files.Jacobian.Header, "lumo_jacobian.c": files.Jacobian.Source,
		"lumo_thermo.h": files.Thermo.Header, "lumo_thermo.c": files.Thermo.Source,
		"reference.h": header.Bytes(), "check.c": harness,
	})
	for _, name := range []string{"reference_residual", "temperature_dependence", "jacobian_values", "jacobian_limits", "jacobian_errors", "jacobian_rank", "newton_reference"} {
		t.Run(name, func(t *testing.T) { program.Run(t, name) })
	}
}
