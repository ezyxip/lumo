package topology

import (
	"os"
	"strings"
	"testing"
)

// mustLoadFixture читает fixture из testdata.
func mustLoadFixture(t *testing.T, path string) Config {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	c, err := Parse(f)
	if err != nil {
		t.Fatalf("Parse(%s): %v", path, err)
	}
	return c
}

// mustValidateFixture парсит и валидирует fixture, ожидая успех.
func mustValidateFixture(t *testing.T, path string) Config {
	c := mustLoadFixture(t, path)
	if err := Validate(c); err != nil {
		t.Fatalf("Validate(%s): unexpected error:\n%v", path, err)
	}
	return c
}

func TestValidateFixtures(t *testing.T) {
	mustValidateFixture(t, "testdata/01_gubkin_octane_decane.json")
	mustValidateFixture(t, "testdata/02_benzene_toluene.json")
	mustValidateFixture(t, "testdata/03_synthetic_wide_boiling.json")
}

// okConfig — валидный базовый конфиг, от которого наследуются ошибочные кейсы.
func okConfig() Config {
	f := 80.1
	return Config{
		Name: "test",
		Components: []Component{
			{Name: "benzene", M: 78.11, Antoine: Antoine{A: 6.90565, B: 1211.033, C: 220.79}, TBoilC: &f},
			{Name: "toluene", M: 92.14, Antoine: Antoine{A: 6.95464, B: 1344.8, C: 219.48}},
		},
		Thermo: Thermo{Model: ThermoModel},
		Column: Column{
			N: 10, Condenser: "total", Reboiler: "kettle",
			FeedStages: []int{5}, PressureMPa: 0.101325, PressureProfile: "const",
		},
	}
}

func TestValidateErrors(t *testing.T) {
	cases := []struct {
		name    string
		mutate  func(*Config)
		wantSub string // подстрока, которая должна быть в ошибке
	}{
		{"missing name", func(c *Config) { c.Name = "" }, "name: required field is missing"},

		{"too few components", func(c *Config) { c.Components = c.Components[:1] }, "components: need 2..3 components, got 1"},
		{"too many components", func(c *Config) {
			c.Components = append(c.Components,
				Component{Name: "x1", M: 1, Antoine: Antoine{A: 1, B: 1, C: 1}},
				Component{Name: "x2", M: 1, Antoine: Antoine{A: 1, B: 1, C: 1}})
		}, "components: need 2..3 components, got 4"},

		{"missing component name", func(c *Config) { c.Components[0].Name = "" }, "components[0].name: required field is missing"},
		{"duplicate component", func(c *Config) { c.Components[1].Name = "benzene" }, "components[1].name=\"benzene\": duplicate component name (already defined at components[0])"},
		{"zero molar mass", func(c *Config) { c.Components[1].M = 0 }, "components[1].M=0: must be > 0"},

		{"zero antoine A", func(c *Config) { c.Components[0].Antoine.A = 0 }, "components[0].antoine.A=0: must be > 0"},
		{"zero antoine B", func(c *Config) { c.Components[1].Antoine.B = 0 }, "components[1].antoine.B=0: must be > 0"},
		{"negative antoine C", func(c *Config) { c.Components[0].Antoine.C = -1 }, "components[0].antoine.C=-1: must be > 0"},

		{"wrong P units", func(c *Config) { c.Components[0].Antoine.PUnits = "kPa" }, "components[0].antoine.P_units=\"kPa\": unsupported unit"},
		{"wrong T units", func(c *Config) { c.Thermo.TUnits = "F" }, "thermo.T_units=\"F\": unsupported unit"},
		{"wrong equation", func(c *Config) { c.Components[1].Antoine.Equation = "ln(P)=A-B/(C+t)" }, "unsupported equation form"},

		{"unknown thermo model", func(c *Config) { c.Thermo.Model = "nrtl" }, "thermo.model=\"nrtl\": unsupported thermo model"},
		{"missing thermo model", func(c *Config) { c.Thermo.Model = "" }, "thermo.model: required field is missing"},

		{"missing N", func(c *Config) { c.Column.N = 0 }, "column.N=0: must be >= 3"},
		{"N too small", func(c *Config) { c.Column.N = 2 }, "column.N=2: must be >= 3"},
		{"missing condenser", func(c *Config) { c.Column.Condenser = "" }, "column.condenser: required field is missing"},
		{"bad condenser", func(c *Config) { c.Column.Condenser = "totalx" }, "column.condenser=\"totalx\""},
		{"missing reboiler", func(c *Config) { c.Column.Reboiler = "" }, "column.reboiler: required field is missing"},
		{"bad reboiler", func(c *Config) { c.Column.Reboiler = "internal" }, "column.reboiler=\"internal\""},
		{"zero pressure", func(c *Config) { c.Column.PressureMPa = 0 }, "column.pressure_MPa=0: must be > 0"},
		{"wrong pressure profile", func(c *Config) { c.Column.PressureProfile = "linear_rise" }, "column.pressure_profile=\"linear_rise\""},

		{"feed on condenser", func(c *Config) { c.Column.FeedStages = []int{1} }, "column.feed_stages[0]=1: stage out of range, must be in [2..N-1] (N=10)"},
		{"feed on reboiler", func(c *Config) { c.Column.FeedStages = []int{10} }, "column.feed_stages[0]=10: stage out of range, must be in [2..N-1] (N=10)"},
		{"feed beyond N", func(c *Config) { c.Column.FeedStages = []int{12} }, "column.feed_stages[0]=12: stage out of range, must be in [2..N-1] (N=10)"},
		{"duplicate feed", func(c *Config) { c.Column.FeedStages = []int{5, 5} }, "column.feed_stages[1]=5: duplicate feed stage"},
		{"side draw beyond N", func(c *Config) {
			c.Column.SideDraws = []SideDraw{{Stage: 11}}
		}, "column.side_draws[0].stage=11: stage out of range, must be in [1..N] (N=10)"},
		{"side draw zero", func(c *Config) {
			c.Column.SideDraws = []SideDraw{{Stage: 0}}
		}, "column.side_draws[0].stage=0: stage out of range"},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			c := okConfig()
			tc.mutate(&c)
			err := Validate(c)
			if err == nil {
				t.Fatalf("Validate() = nil, want error containing %q", tc.wantSub)
			}
			if !strings.Contains(err.Error(), tc.wantSub) {
				t.Errorf("Validate() error:\n%v\n\ndoes not contain %q", err, tc.wantSub)
			}
		})
	}
}

func TestValidateCollectsAllIssues(t *testing.T) {
	c := okConfig()
	c.Name = ""
	c.Components[1].M = -5
	c.Column.N = 2
	c.Column.FeedStages = []int{50}

	err := Validate(c)
	verr, ok := err.(*ValidationError)
	if !ok {
		t.Fatalf("Validate() = %T, want *ValidationError", err)
	}
	if got := len(verr.Issues); got != 4 {
		t.Fatalf("got %d issues, want 4:\n%v", got, err)
	}
}

func TestValidationErrorFormat(t *testing.T) {
	c := okConfig()
	c.Column.N = 2
	c.Column.FeedStages = nil // иначе feed-ошибка добавит второй issue
	err := Validate(c)
	verr := err.(*ValidationError)
	if len(verr.Issues) != 1 {
		t.Fatalf("got %d issues, want 1", len(verr.Issues))
	}
	want := "invalid compiler input: 1 issue(s):\n  - column.N=2: must be >= 3 (condenser + at least 1 stage + reboiler)"
	if err.Error() != want {
		t.Errorf("Error() =\n%q\nwant\n%q", err.Error(), want)
	}
}
