package contractv1

import (
	"bytes"
	"encoding/json"
	"os"
	"reflect"
	"testing"
)

// Ensure the published examples can be read and written through the Go wire
// types without silently dropping fields or introducing fields from other branches.
func TestExampleRoundTrip(t *testing.T) {
	for _, tc := range []struct {
		name  string
		value any
	}{
		{"compiler_input", &CompilerInput{}},
		{"model_input", &ModelInput{}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			data, err := os.ReadFile("../../examples/v1/benzene_toluene/" + tc.name + ".json")
			if err != nil {
				t.Fatal(err)
			}
			dec := json.NewDecoder(bytes.NewReader(data))
			dec.DisallowUnknownFields()
			if err := dec.Decode(tc.value); err != nil {
				t.Fatal(err)
			}
			encoded, err := json.Marshal(tc.value)
			if err != nil {
				t.Fatal(err)
			}
			var before, after any
			if err := json.Unmarshal(data, &before); err != nil {
				t.Fatal(err)
			}
			if err := json.Unmarshal(encoded, &after); err != nil {
				t.Fatal(err)
			}
			if !reflect.DeepEqual(before, after) {
				t.Fatalf("round trip changed the contract:\n%s", encoded)
			}
		})
	}
}

func TestZeroVaporFractionIsNotOmitted(t *testing.T) {
	zero := 0.0
	data, err := json.Marshal(FeedState{Type: "vapor_fraction", VaporFraction: &zero})
	if err != nil {
		t.Fatal(err)
	}
	if string(data) != `{"type":"vapor_fraction","vapor_fraction":0}` {
		t.Fatalf("zero vapor fraction lost: %s", data)
	}
}

func TestCondenserHasNoOutgoingVapor(t *testing.T) {
	data, err := json.Marshal(StageProfile{Stage: 1, TemperatureK: 350,
		Liquid: Stream{FlowMolS: 1, MoleFractions: map[string]float64{"a": 0.5, "b": 0.5}}})
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Contains(data, []byte(`"vapor":null`)) {
		t.Fatalf("missing explicit absent vapor: %s", data)
	}
}
