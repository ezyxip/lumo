package contractv1

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"reflect"
	"strings"
)

// ParseCompilerInput strictly decodes v1 JSON. Call ValidateCompilerInput next
// to check semantics. Parsing preserves valid zero-valued coefficients.
func ParseCompilerInput(r io.Reader) (CompilerInput, error) {
	return parse[CompilerInput](r, "compiler input")
}

// ParseModelInput strictly decodes operating parameters. ValidateModelInput
// checks their semantics and references against a compiler input.
func ParseModelInput(r io.Reader) (ModelInput, error) {
	return parse[ModelInput](r, "model input")
}

func parse[T any](r io.Reader, name string) (T, error) {
	var result T
	data, err := io.ReadAll(r)
	if err != nil {
		return result, fmt.Errorf("invalid %s: reading: %w", name, err)
	}
	dec := json.NewDecoder(bytes.NewReader(data))
	dec.UseNumber()
	value, err := readValue(dec)
	if err == nil {
		if _, tail := dec.Token(); tail != io.EOF {
			err = fmt.Errorf("unexpected data after JSON value")
		}
	}
	if err == nil {
		err = checkStructure(value, reflect.TypeFor[T](), "$")
	}
	if err == nil {
		err = json.Unmarshal(data, &result)
	}
	if err != nil {
		var zero T
		return zero, fmt.Errorf("invalid %s: %w", name, err)
	}
	return result, nil
}

// Read objects without silently accepting duplicate keys, including map keys.
func readValue(dec *json.Decoder) (any, error) {
	token, err := dec.Token()
	if err != nil {
		return nil, err
	}
	switch token {
	case json.Delim('{'):
		value := map[string]any{}
		for dec.More() {
			keyToken, err := dec.Token()
			if err != nil {
				return nil, err
			}
			key := keyToken.(string)
			if _, exists := value[key]; exists {
				return nil, fmt.Errorf("duplicate field %q", key)
			}
			child, err := readValue(dec)
			if err != nil {
				return nil, err
			}
			value[key] = child
		}
		_, err = dec.Token()
		return value, err
	case json.Delim('['):
		value := []any{}
		for dec.More() {
			child, err := readValue(dec)
			if err != nil {
				return nil, err
			}
			value = append(value, child)
		}
		_, err = dec.Token()
		return value, err
	default:
		return token, nil
	}
}

// Required fields follow the wire types' json tags. Unlike encoding/json,
// field names are case-sensitive and null cannot replace a value/object.
func checkStructure(value any, typ reflect.Type, path string) error {
	if value == nil {
		return fmt.Errorf("%s: null is not allowed", path)
	}
	if typ.Kind() == reflect.Pointer {
		typ = typ.Elem()
	}
	switch typ.Kind() {
	case reflect.Struct:
		object, ok := value.(map[string]any)
		if !ok {
			return fmt.Errorf("%s: expected object", path)
		}
		fields := map[string]reflect.StructField{}
		for i := 0; i < typ.NumField(); i++ {
			field := typ.Field(i)
			name, options, _ := strings.Cut(field.Tag.Get("json"), ",")
			fields[name] = field
			if _, present := object[name]; !present && options != "omitempty" {
				return fmt.Errorf("%s.%s: required field is missing", path, name)
			}
		}
		// Branch-specific fields are forbidden even if explicitly zero/empty.
		var allowed map[string]bool
		if typ == reflect.TypeFor[Properties]() {
			switch object["source"] {
			case "builtin":
				allowed = fieldSet("source", "key", "revision")
			case "inline":
				allowed = fieldSet("source", "antoine", "enthalpy", "temperature_range_K")
			}
		}
		if typ == reflect.TypeFor[FeedState]() {
			switch object["type"] {
			case "saturated_liquid":
				allowed = fieldSet("type")
			case "temperature":
				allowed = fieldSet("type", "temperature_K")
			case "vapor_fraction":
				allowed = fieldSet("type", "vapor_fraction")
			}
		}
		for _, name := range sortedKeys(allowed) {
			if _, present := object[name]; !present {
				return fmt.Errorf("%s.%s: required field is missing", path, name)
			}
		}
		for _, name := range sortedKeys(object) {
			field, ok := fields[name]
			if !ok || (allowed != nil && !allowed[name]) {
				return fmt.Errorf("%s.%s: unknown or incompatible field", path, name)
			}
			if err := checkStructure(object[name], field.Type, path+"."+name); err != nil {
				return err
			}
		}
	case reflect.Slice:
		array, ok := value.([]any)
		if !ok {
			return fmt.Errorf("%s: expected array", path)
		}
		for i, child := range array {
			if err := checkStructure(child, typ.Elem(), fmt.Sprintf("%s[%d]", path, i)); err != nil {
				return err
			}
		}
	case reflect.Map:
		object, ok := value.(map[string]any)
		if !ok {
			return fmt.Errorf("%s: expected object", path)
		}
		for _, key := range sortedKeys(object) {
			if err := checkStructure(object[key], typ.Elem(), fmt.Sprintf("%s[%q]", path, key)); err != nil {
				return err
			}
		}
	}
	return nil
}

func fieldSet(names ...string) map[string]bool {
	result := make(map[string]bool, len(names))
	for _, name := range names {
		result[name] = true
	}
	return result
}
