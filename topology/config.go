// Package topology читает описание топологии колонны (compiler_input.json)
// во внутреннее представление компилятора и проверяет его корректность.
//
// Формат входа зафиксирован схемой examples/schema/compiler_input.schema.json.
package topology

// Config — внутреннее представление входа компилятора (Level 1).
// Фиксирует структуру MESH-задачи: компоненты, термодинамику, топологию колонны.
type Config struct {
	Name        string      `json:"name"`
	Description string      `json:"description,omitempty"`
	Provenance  Provenance  `json:"provenance"`
	Components  []Component `json:"components"`
	Thermo      Thermo      `json:"thermo"`
	Column      Column      `json:"column"`
}

// Provenance — источник данных примера (для трассировки эталонов).
type Provenance struct {
	Source string `json:"source,omitempty"`
	Pages  string `json:"pages,omitempty"`
}

// Component — компонент смеси с константами уравнения Антуана.
// Antoine: lg P = A - B / (C + t), P в мм рт. ст. (по умолчанию), t в °C.
type Component struct {
	Name    string   `json:"name"`
	M       float64  `json:"M"` // молярная масса, г/моль
	Antoine Antoine  `json:"antoine"`
	TBoilC  *float64 `json:"T_boil_C_at_P,omitempty"` // опционально, для контроля
}

// Antoine — константы уравнения Антуана для компонента.
type Antoine struct {
	A        float64 `json:"A"`
	B        float64 `json:"B"`
	C        float64 `json:"C"`
	PUnits   string  `json:"P_units,omitempty"`  // default "mmHg"
	TUnits   string  `json:"T_units,omitempty"`  // default "C"
	Equation string  `json:"equation,omitempty"` // default "log10(P)=A-B/(C+t)"
}

// Thermo — модель термодинамики. В V0 поддерживается только идеальная смесь.
type Thermo struct {
	Model  string `json:"model"` // enum: "antoine_raoult"
	PUnits string `json:"P_units,omitempty"`
	TUnits string `json:"T_units,omitempty"`
}

// Column — топология колонны.
type Column struct {
	N               int        `json:"N"`         // теоретические ступени, включая condenser и reboiler
	Condenser       string     `json:"condenser"` // enum: "total" | "partial"
	Reboiler        string     `json:"reboiler"`  // enum: "kettle" | "thermosiphon"
	FeedStages      []int      `json:"feed_stages,omitempty"`
	SideDraws       []SideDraw `json:"side_draws,omitempty"`
	PressureMPa     float64    `json:"pressure_MPa"`
	PressureProfile string     `json:"pressure_profile,omitempty"` // enum: "const" | "linear_drop", default "const"
}

// SideDraw — боковой отбор. Схема пока описывает его как object без полей;
// минимально фиксируем номер ступени. TODO(#1): уточнить схему side_draws.
type SideDraw struct {
	Stage int `json:"stage"`
}
