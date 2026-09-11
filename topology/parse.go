package topology

import (
	"bytes"
	"encoding/json"
	"fmt"
	"io"
)

// Parse читает compiler_input.json и возвращает Config.
//
// Ошибки делятся на два типа:
//   - "invalid compiler input: <деталь>" — битый JSON, неверная структура
//     (ожидали object, получили array), неизвестные поля;
//   - семантические ошибки (missing field, битые ссылки, диапазоны) здесь
//     НЕ проверяются — их ловит Validate (см. validate.go). Пропущенное
//     required-поле детектируется там как невалидное нулевое значение.
func Parse(r io.Reader) (Config, error) {
	data, err := io.ReadAll(r)
	if err != nil {
		return Config{}, fmt.Errorf("invalid compiler input: reading: %w", err)
	}

	dec := json.NewDecoder(bytes.NewReader(data))
	dec.DisallowUnknownFields()
	var cfg Config
	if err := dec.Decode(&cfg); err != nil {
		return Config{}, fmt.Errorf("invalid compiler input: %w", err)
	}
	if _, err := dec.Token(); err != io.EOF {
		return Config{}, fmt.Errorf("invalid compiler input: unexpected data after JSON object")
	}
	return cfg, nil
}
