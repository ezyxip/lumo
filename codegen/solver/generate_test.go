package solver

import (
	"bytes"
	"encoding/json"
	"fmt"
	"os"
	"regexp"
	"testing"

	"github.com/ezyxip/lumo/codegen/initialguess"
	"github.com/ezyxip/lumo/codegen/internal/ctest"
	"github.com/ezyxip/lumo/codegen/prepare"
	contractv1 "github.com/ezyxip/lumo/contracts/v1"
)

func config(t *testing.T) contractv1.CompilerInput {
	t.Helper()
	f, e := os.Open("../../contracts/v1/testdata/compiler_input.json")
	if e != nil {
		t.Fatal(e)
	}
	defer f.Close()
	c, e := contractv1.ParseCompilerInput(f)
	if e != nil {
		t.Fatal(e)
	}
	return c
}
func sourceFiles(t *testing.T, c contractv1.CompilerInput, h []byte) map[string][]byte {
	t.Helper()
	s, e := Generate(c)
	if e != nil {
		t.Fatal(e)
	}
	g, e := initialguess.Generate(c)
	if e != nil {
		t.Fatal(e)
	}
	p, e := prepare.Generate(c)
	if e != nil {
		t.Fatal(e)
	}
	return map[string][]byte{"lumo_solver.h": s.Header, "lumo_solver.c": s.Source,
		"lumo_mesh.h": s.Mesh.Header, "lumo_mesh.c": s.Mesh.Source, "lumo_jacobian.h": s.Mesh.Jacobian.Header, "lumo_jacobian.c": s.Mesh.Jacobian.Source,
		"lumo_thermo.h": s.Mesh.Thermo.Header, "lumo_thermo.c": s.Mesh.Thermo.Source,
		"lumo_initialguess.h": g.Header, "lumo_initialguess.c": g.Source, "lumo_prepare.h": p.Header, "lumo_prepare.c": p.Source, "check.c": h}
}
func TestSolverC99(t *testing.T) {
	c := config(t)
	f, e := Generate(c)
	if e != nil {
		t.Fatal(e)
	}
	again, e := Generate(c)
	if e != nil {
		t.Fatal(e)
	}
	if !bytes.Equal(f.Header, again.Header) || !bytes.Equal(f.Source, again.Source) {
		t.Fatal("nondeterministic generation")
	}
	h, e := os.ReadFile("testdata/check.c")
	if e != nil {
		t.Fatal(e)
	}
	data, e := os.ReadFile("../mesh/testdata/nonlinear/reference.json")
	if e != nil {
		t.Fatal(e)
	}
	var reference struct {
		State []float64 `json:"state"`
	}
	if e = json.Unmarshal(data, &reference); e != nil {
		t.Fatal(e)
	}
	if len(reference.State) != f.Mesh.Variables {
		t.Fatal("wrong reference dimensions")
	}
	var header bytes.Buffer
	header.WriteString("static const double reference_state[LUMO_MESH_SIZE]={\n")
	for _, x := range reference.State {
		fmt.Fprintf(&header, "%.17g,\n", x)
	}
	header.WriteString("};\n")
	files := sourceFiles(t, c, h)
	files["reference.h"] = header.Bytes()
	p := ctest.Build(t, files)
	for _, name := range []string{"reference", "already_converged", "limits", "callbacks", "errors", "scales", "operating_variants", "workspace_contract"} {
		t.Run(name, func(t *testing.T) { p.Run(t, name) })
	}
}
func TestSolverOtherDimensions(t *testing.T) {
	h, e := os.ReadFile("testdata/dimensions.c")
	if e != nil {
		t.Fatal(e)
	}
	for _, name := range []string{"five_stages", "ternary", "multiple_feeds", "same_stage"} {
		t.Run(name, func(t *testing.T) {
			c := config(t)
			switch name {
			case "five_stages":
				c.Column.Stages = 5
				c.Column.Feeds[0].Stage = 3
			case "ternary":
				third := c.Components[1]
				third.ID = "third"
				c.Components = append(c.Components, third)
			case "multiple_feeds":
				c.Column.Stages = 5
				c.Column.Feeds = append(c.Column.Feeds, contractv1.FeedLocation{ID: "second", Stage: 4})
			case "same_stage":
				c.Column.Feeds = append(c.Column.Feeds, contractv1.FeedLocation{ID: "second", Stage: 2})
			}
			p := ctest.Build(t, sourceFiles(t, c, h))
			p.Run(t, "dimensions")
		})
	}
}
func TestRejectInvalidConfiguration(t *testing.T) {
	c := config(t)
	c.Column.Stages = 1
	if f, e := Generate(c); e == nil || len(f.Source) != 0 {
		t.Fatal("invalid description accepted")
	}
}
func TestSingularColumn(t *testing.T) {
	c := config(t)
	c.Column.Pressure.ValuePa = 133322.4
	for i := range c.Components {
		c.Components[i].Properties.Antoine = &contractv1.Antoine{A: 3, B: 0, C: 1, PressureUnit: "mmHg", TemperatureUnit: "C"}
	}
	h := []byte(`#include "lumo_solver.h"
#include "ctest.h"
static double workspace[LUMO_SOLVER_WORKSPACE_DOUBLES];
static int singular_column(void){
 const struct lumo_mesh_operating op={.feeds={{100,{0.5,0.5},2000}},.distillate_mol_s=50,.reflux_ratio=2.5};
 double state[LUMO_MESH_SIZE]={0},before[LUMO_MESH_SIZE];struct lumo_solver_report report;size_t j;
 for(j=0;j<LUMO_MESH_STAGES;++j){size_t b=j*LUMO_MESH_BLOCK;state[b]=0.5;state[b+1]=0.5;state[b+2]=0.5;state[b+3]=0.5;state[b+LUMO_MESH_T]=360;state[b+LUMO_MESH_L]=j==2?50:125;state[b+LUMO_MESH_V]=175;}
 memcpy(before,state,sizeof state);
 CTEST_STATUS(lumo_solver_solve(&op,state,NULL,workspace,LUMO_SOLVER_WORKSPACE_DOUBLES,&report),LUMO_SOLVER_SINGULAR_JACOBIAN);
 CTEST_STATUS(report.iterations,0);CTEST_STATUS(memcmp(before,state,sizeof state),0);return 0;
}
int main(int argc,char **argv){const struct ctest_case cases[]={{"singular_column",singular_column}};return ctest_run(argc,argv,cases,1);}
`)
	p := ctest.Build(t, sourceFiles(t, c, h))
	p.Run(t, "singular_column")
}

func TestDenseLinearSolve(t *testing.T) {
	c := config(t)
	f, err := Generate(c)
	if err != nil {
		t.Fatal(err)
	}
	h := []byte(`#include "solver_under_test.inc"
#include "ctest.h"
static int linear_reference(void){
 double a[LUMO_MESH_SIZE*LUMO_MESH_SIZE]={0},b[LUMO_MESH_SIZE]={0},x[LUMO_MESH_SIZE],expected[LUMO_MESH_SIZE];size_t i,j;
 for(i=0;i<N;++i){a[i*N+i]=2;expected[i]=((int)(i%3)-1)*0.5;if(i+1<N)a[i*N+i+1]=0.1;if(i>0)a[i*N+i-1]=0.2;}
 for(j=0;j<N;++j){double v=a[j];a[j]=a[(N-1)*N+j];a[(N-1)*N+j]=v;}
 for(i=0;i<N;++i)for(j=0;j<N;++j)b[i]+=a[i*N+j]*expected[j];
 CTEST_STATUS(linear_solve(a,b,x),LUMO_SOLVER_CONVERGED);
 for(i=0;i<N;++i){CTEST_NEAR(x[i],expected[i]);}
 memset(a,0,sizeof a);for(i=0;i<N;++i)a[i*N+i]=1;a[(N-1)*N+N-1]=1e-20;
 CTEST_STATUS(linear_solve(a,b,x),LUMO_SOLVER_SINGULAR_JACOBIAN);return 0;
}
int main(int argc,char **argv){const struct ctest_case cases[]={{"linear_reference",linear_reference}};return ctest_run(argc,argv,cases,1);}
`)
	files := sourceFiles(t, c, h)
	delete(files, "lumo_solver.c")
	files["solver_under_test.inc"] = f.Source
	p := ctest.Build(t, files)
	p.Run(t, "linear_reference")
}

// A library build must remain free of heap-allocation calls in every generated
// module, including preparation and initialization, not just the solver body.
func TestGeneratedModulesHaveNoHeapCalls(t *testing.T) {
	files := sourceFiles(t, config(t), nil)
	allocator := regexp.MustCompile(`\b(malloc|calloc|realloc|free|aligned_alloc|posix_memalign|alloca|strdup|strndup|printf|fprintf|puts|fputs|fopen|fwrite|write)\s*\(`)
	for name, data := range files {
		if location := allocator.FindIndex(data); location != nil {
			t.Errorf("%s contains an allocation or output call: %s", name, data[location[0]:location[1]])
		}
	}
}
