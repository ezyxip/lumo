package topology

import (
	"strings"
	"testing"
)

func TestParseFixture(t *testing.T) {
	c := mustLoadFixture(t, "testdata/02_benzene_toluene.json")

	if c.Name != "02_benzene_toluene" {
		t.Errorf("Name = %q, want %q", c.Name, "02_benzene_toluene")
	}
	wantDesc := "Канонический идеальный пример бензол/толуол, P=1 атм, для кросс-чека с Raoult. Константы из NIST/BCCampus."
	if c.Description != wantDesc {
		t.Errorf("Description = %q, want %q", c.Description, wantDesc)
	}
	wantSource := "pressbooks.bccampus.ca/chbe220, NIST Antoine"
	if c.Provenance.Source != wantSource {
		t.Errorf("Provenance.Source = %q, want %q", c.Provenance.Source, wantSource)
	}
	if c.Provenance.Pages != "Txy benzene-toluene 1 atm" {
		t.Errorf("Provenance.Pages = %q, want %q", c.Provenance.Pages, "Txy benzene-toluene 1 atm")
	}

	if len(c.Components) != 2 {
		t.Fatalf("len(Components) = %d, want 2", len(c.Components))
	}
	benzene := c.Components[0]
	if benzene.Name != "benzene" {
		t.Errorf("Components[0].Name = %q, want %q", benzene.Name, "benzene")
	}
	if benzene.M != 78.11 {
		t.Errorf("Components[0].M = %v, want 78.11", benzene.M)
	}
	wantAntoine := Antoine{A: 6.90565, B: 1211.033, C: 220.79, PUnits: "mmHg", TUnits: "C", Equation: "log10(P)=A-B/(C+t)"}
	if benzene.Antoine != wantAntoine {
		t.Errorf("Components[0].Antoine = %+v, want %+v", benzene.Antoine, wantAntoine)
	}
	if benzene.TBoilC == nil || *benzene.TBoilC != 80.1 {
		t.Errorf("Components[0].TBoilC = %v, want 80.1", benzene.TBoilC)
	}
	toluene := c.Components[1]
	if toluene.Name != "toluene" || toluene.M != 92.14 {
		t.Errorf("Components[1] = {%s, %v}, want {toluene, 92.14}", toluene.Name, toluene.M)
	}

	if c.Thermo.Model != "antoine_raoult" {
		t.Errorf("Thermo.Model = %q, want %q", c.Thermo.Model, "antoine_raoult")
	}

	col := c.Column
	if col.N != 10 {
		t.Errorf("Column.N = %d, want 10", col.N)
	}
	if col.Condenser != "total" || col.Reboiler != "kettle" {
		t.Errorf("Column types = %q/%q, want total/kettle", col.Condenser, col.Reboiler)
	}
	if len(col.FeedStages) != 1 || col.FeedStages[0] != 5 {
		t.Errorf("Column.FeedStages = %v, want [5]", col.FeedStages)
	}
	if len(col.SideDraws) != 0 {
		t.Errorf("Column.SideDraws = %v, want empty", col.SideDraws)
	}
	if col.PressureMPa != 0.101325 {
		t.Errorf("Column.PressureMPa = %v, want 0.101325", col.PressureMPa)
	}
	if col.PressureProfile != "const" {
		t.Errorf("Column.PressureProfile = %q, want %q", col.PressureProfile, "const")
	}
}

func TestParseBadJSON(t *testing.T) {
	cases := []struct {
		name string
		in   string
	}{
		{"truncated", `{"name": "x", "components": [`},
		{"not-an-object", `[1, 2, 3]`},
		{"garbage", `not json at all`},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			_, err := Parse(strings.NewReader(tc.in))
			if err == nil {
				t.Fatalf("Parse(%q) = nil error, want error", tc.in)
			}
		})
	}
}

func TestParseUnknownField(t *testing.T) {
	_, err := Parse(strings.NewReader(`{"name": "x", "oops": 1}`))
	if err == nil {
		t.Fatal("Parse with unknown field = nil error, want error")
	}
}
