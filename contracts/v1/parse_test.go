package contractv1

import (
	"encoding/json"
	"errors"
	"io"
	"os"
	"strings"
	"testing"
)

func fixtureCompiler(t *testing.T) CompilerInput {
	t.Helper()
	data, err := os.ReadFile("testdata/compiler_input.json")
	if err != nil {
		t.Fatal(err)
	}
	c, err := ParseCompilerInput(strings.NewReader(string(data)))
	if err != nil {
		t.Fatal(err)
	}
	return c
}
func fixtureModel(t *testing.T) ModelInput {
	t.Helper()
	data, err := os.ReadFile("../../examples/v1/benzene_toluene/model_input.json")
	if err != nil {
		t.Fatal(err)
	}
	c, err := ParseModelInput(strings.NewReader(string(data)))
	if err != nil {
		t.Fatal(err)
	}
	return c
}

func TestParseFixtures(t *testing.T) {
	inline := fixtureCompiler(t)
	if inline.Column.Stages != 3 || inline.Components[0].ID != "benzene" {
		t.Fatal("incorrect inline input")
	}
	data, err := os.ReadFile("../../examples/v1/benzene_toluene/compiler_input.json")
	if err != nil {
		t.Fatal(err)
	}
	builtin, err := ParseCompilerInput(strings.NewReader(string(data)))
	if err != nil {
		t.Fatal(err)
	}
	if builtin.Components[0].Properties.Source != "builtin" {
		t.Fatal("lost builtin reference")
	}
	input := fixtureModel(t)
	if input.Feeds["main_feed"].Flow.Value != 100 || input.Specifications.RefluxRatio.Value != 2.5 {
		t.Fatal("incorrect model input")
	}
	// Parsing and validation are deliberately separate operations.
	inline.Components[0].Properties.Antoine.B = 0
	inline.Components[0].Properties.Antoine.C = -1
	encoded, err := json.Marshal(inline)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := ParseCompilerInput(strings.NewReader(string(encoded))); err != nil {
		t.Fatal(err)
	}
}

func TestParseRejectsInvalidJSON(t *testing.T) {
	for _, input := range []string{
		"", `{"schema_version":`, `[1,2]`, `null`, `42`, `"text"`,
		`{}`, `{"schema_version":1,"schema_version":1}`, `{} {}`, `{} trailing`,
	} {
		t.Run(input, func(t *testing.T) {
			if _, err := ParseCompilerInput(strings.NewReader(input)); err == nil {
				t.Fatal("compiler parser accepted invalid input")
			}
			if _, err := ParseModelInput(strings.NewReader(input)); err == nil {
				t.Fatal("model parser accepted invalid input")
			}
		})
	}
}

func TestCompilerStructure(t *testing.T) {
	for _, tc := range []struct {
		name string
		edit func(map[string]any)
	}{
		{"unknown", func(m map[string]any) { m["oops"] = true }},
		{"wrong case", func(m map[string]any) { m["Name"] = m["name"]; delete(m, "name") }},
		{"null optional", func(m map[string]any) { m["provenance"] = nil }},
		{"null required", func(m map[string]any) { m["column"] = nil }},
		{"null number", func(m map[string]any) { component(m)["molar_mass_g_mol"] = nil }},
		{"missing coefficient", func(m map[string]any) { delete(properties(m)["antoine"].(map[string]any), "B") }},
		{"inline extra empty key", func(m map[string]any) { properties(m)["key"] = "" }},
		{"missing range", func(m map[string]any) { delete(properties(m), "temperature_range_K") }},
		{"null component", func(m map[string]any) { m["components"].([]any)[0] = nil }},
		{"wrong array", func(m map[string]any) { m["components"] = map[string]any{} }},
		{"wrong number", func(m map[string]any) { component(m)["molar_mass_g_mol"] = "78" }},
		{"fractional stages", func(m map[string]any) { m["column"].(map[string]any)["stages"] = 3.5 }},
		{"overflow", func(m map[string]any) { component(m)["molar_mass_g_mol"] = json.Number("1e999") }},
		{"builtin inline", func(m map[string]any) {
			properties(m)["source"] = "builtin"
			properties(m)["key"] = "x"
			properties(m)["revision"] = "1"
		}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			m := object(t, fixtureCompiler(t))
			tc.edit(m)
			data, err := json.Marshal(m)
			if err != nil {
				t.Fatal(err)
			}
			if _, err := ParseCompilerInput(strings.NewReader(string(data))); err == nil {
				t.Fatal("accepted invalid structure")
			}
		})
	}
	data, _ := json.Marshal(fixtureCompiler(t))
	duplicate := strings.Replace(string(data), `"stages":3`, `"stages":3,"stages":4`, 1)
	if _, err := ParseCompilerInput(strings.NewReader(duplicate)); err == nil {
		t.Fatal("accepted nested duplicate")
	}
}

func TestModelStructure(t *testing.T) {
	for _, tc := range []struct {
		name string
		edit func(map[string]any)
	}{
		{"missing zero reflux", func(m map[string]any) {
			delete(m["specifications"].(map[string]any)["reflux_ratio"].(map[string]any), "value")
		}},
		{"state extra zero", func(m map[string]any) { feedObject(m)["state"].(map[string]any)["vapor_fraction"] = 0.0 }},
		{"missing vapor fraction", func(m map[string]any) { feedObject(m)["state"] = map[string]any{"type": "vapor_fraction"} }},
		{"missing temperature", func(m map[string]any) { feedObject(m)["state"] = map[string]any{"type": "temperature"} }},
		{"null fraction", func(m map[string]any) {
			feedObject(m)["composition"].(map[string]any)["fractions"].(map[string]any)["benzene"] = nil
		}},
		{"unknown nested", func(m map[string]any) { feedObject(m)["flow"].(map[string]any)["pressure"] = 1.0 }},
		{"null feed", func(m map[string]any) { m["feeds"].(map[string]any)["main_feed"] = nil }},
	} {
		t.Run(tc.name, func(t *testing.T) {
			m := object(t, fixtureModel(t))
			tc.edit(m)
			data, err := json.Marshal(m)
			if err != nil {
				t.Fatal(err)
			}
			if _, err := ParseModelInput(strings.NewReader(string(data))); err == nil {
				t.Fatal("accepted invalid structure")
			}
		})
	}
	input := fixtureModel(t)
	feed := input.Feeds["main_feed"]
	zero := 0.0
	feed.State = FeedState{Type: "vapor_fraction", VaporFraction: &zero}
	input.Feeds["main_feed"] = feed
	input.Specifications.RefluxRatio.Value = 0
	data, _ := json.Marshal(input)
	parsed, err := ParseModelInput(strings.NewReader(string(data)))
	if err != nil || parsed.Feeds["main_feed"].State.VaporFraction == nil {
		t.Fatalf("zero parameter lost: %v", err)
	}
	duplicate := strings.Replace(string(data), `"benzene":0.5`, `"benzene":0.5,"benzene":0.5`, 1)
	if _, err := ParseModelInput(strings.NewReader(duplicate)); err == nil {
		t.Fatal("accepted duplicate composition key")
	}
}

func object(t *testing.T, value any) map[string]any {
	t.Helper()
	data, err := json.Marshal(value)
	if err != nil {
		t.Fatal(err)
	}
	var out map[string]any
	if err := json.Unmarshal(data, &out); err != nil {
		t.Fatal(err)
	}
	return out
}
func component(m map[string]any) map[string]any  { return m["components"].([]any)[0].(map[string]any) }
func properties(m map[string]any) map[string]any { return component(m)["properties"].(map[string]any) }
func feedObject(m map[string]any) map[string]any {
	return m["feeds"].(map[string]any)["main_feed"].(map[string]any)
}

type failedReader struct{}

func (failedReader) Read([]byte) (int, error) { return 0, io.ErrUnexpectedEOF }
func TestParseReadError(t *testing.T) {
	if _, err := ParseCompilerInput(failedReader{}); !errors.Is(err, io.ErrUnexpectedEOF) {
		t.Fatalf("lost read error: %v", err)
	}
	if _, err := ParseModelInput(failedReader{}); !errors.Is(err, io.ErrUnexpectedEOF) {
		t.Fatalf("lost read error: %v", err)
	}
}
