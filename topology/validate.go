package topology

import (
	"fmt"
	"strings"
)

// Допустимые значения, зашитые скоупом V0 (см. doc/adr/0001-thermo-models-scope.md).
// Уравнение Антуана и единицы зафиксированы: генерируемый код интерпретирует
// константы строго в этих единицах, другие значения не имеют смысла.
const (
	ThermoModel     = "antoine_raoult"
	PressureUnit    = "mmHg"
	TemperatureUnit = "C"
	AntoineEquation = "log10(P)=A-B/(C+t)"
)

// ValidationError — агрегированная ошибка валидации со всеми найденными проблемами.
type ValidationError struct {
	Issues []string
}

func (e *ValidationError) Error() string {
	var b strings.Builder
	fmt.Fprintf(&b, "invalid compiler input: %d issue(s):", len(e.Issues))
	for _, s := range e.Issues {
		b.WriteString("\n  - ")
		b.WriteString(s)
	}
	return b.String()
}

// Validate проверяет семантику Config после успешного Parse и собирает все
// проблемы за один проход (а не первую встреченную).
//
// Отсутствующие required-поля детектируются как невалидные нулевые значения:
// у всех required-полей схемы нулевое значение физически невалидно
// (N=0, M=0, pressure=0, пустой enum и т.д.).
//
// Нормализация единиц и другие «удобства» сознательно не входят сюда —
// см. doc/adr/0003-strict-units-validation.md.
func Validate(c Config) error {
	v := &validator{}
	v.name(c.Name)
	v.components(c.Components)
	v.thermo(c.Thermo)
	v.column(c.Column)
	if len(v.issues) == 0 {
		return nil
	}
	return &ValidationError{Issues: v.issues}
}

type validator struct {
	issues []string
}

func (v *validator) addf(format string, args ...any) {
	v.issues = append(v.issues, fmt.Sprintf(format, args...))
}

func (v *validator) name(name string) {
	if name == "" {
		v.addf("name: required field is missing")
	}
}

func (v *validator) components(comps []Component) {
	// examples/schema/compiler_input.schema.json: minItems 2, maxItems 3.
	if n := len(comps); n < 2 || n > 3 {
		v.addf("components: need 2..3 components, got %d", n)
	}
	seen := make(map[string]int, len(comps))
	for i, comp := range comps {
		path := fmt.Sprintf("components[%d]", i)
		switch comp.Name {
		case "":
			v.addf("%s.name: required field is missing", path)
		default:
			if prev, dup := seen[comp.Name]; dup {
				v.addf("%s.name=%q: duplicate component name (already defined at components[%d])", path, comp.Name, prev)
			} else {
				seen[comp.Name] = i
			}
		}
		if comp.M <= 0 {
			v.addf("%s.M=%v: must be > 0 (molar mass)", path, comp.M)
		}
		v.antoine(path+".antoine", comp.Antoine)
	}
}

func (v *validator) antoine(path string, a Antoine) {
	if a.A <= 0 {
		v.addf("%s.A=%v: must be > 0", path, a.A)
	}
	if a.B <= 0 {
		v.addf("%s.B=%v: must be > 0", path, a.B)
	}
	if a.C <= 0 {
		v.addf("%s.C=%v: must be > 0", path, a.C)
	}
	v.unit(path+".P_units", a.PUnits, PressureUnit)
	v.unit(path+".T_units", a.TUnits, TemperatureUnit)
	if a.Equation != "" && a.Equation != AntoineEquation {
		v.addf("%s=%q: unsupported equation form, must be %q or empty (default)", path, a.Equation, AntoineEquation)
	}
}

func (v *validator) thermo(t Thermo) {
	switch t.Model {
	case "":
		v.addf("thermo.model: required field is missing")
	case ThermoModel:
		// ok
	default:
		v.addf("thermo.model=%q: unsupported thermo model, only %q is supported in V0 (see doc/adr/0001-thermo-models-scope.md)", t.Model, ThermoModel)
	}
	v.unit("thermo.P_units", t.PUnits, PressureUnit)
	v.unit("thermo.T_units", t.TUnits, TemperatureUnit)
}

// unit — строгая проверка единиц: пусто (= default) или точное значение.
// См. doc/adr/0003-strict-units-validation.md.
func (v *validator) unit(field, got, want string) {
	if got == "" || got == want {
		return
	}
	v.addf("%s=%q: unsupported unit, must be %q or empty (default %q)", field, got, want, want)
}

func (v *validator) column(col Column) {
	switch {
	case col.N < 3:
		v.addf("column.N=%d: must be >= 3 (condenser + at least 1 stage + reboiler)", col.N)
	}

	switch col.Condenser {
	case "":
		v.addf("column.condenser: required field is missing")
	case "total", "partial":
		// ok
	default:
		v.addf("column.condenser=%q: must be \"total\" or \"partial\"", col.Condenser)
	}

	switch col.Reboiler {
	case "":
		v.addf("column.reboiler: required field is missing")
	case "kettle", "thermosiphon":
		// ok
	default:
		v.addf("column.reboiler=%q: must be \"kettle\" or \"thermosiphon\"", col.Reboiler)
	}

	if col.PressureMPa <= 0 {
		v.addf("column.pressure_MPa=%v: must be > 0", col.PressureMPa)
	}

	switch col.PressureProfile {
	case "", "const", "linear_drop":
		// ok, "" = default "const"
	default:
		v.addf("column.pressure_profile=%q: must be \"const\" or \"linear_drop\"", col.PressureProfile)
	}

	// Ссылки на ступени. Диапазон фида — [2..N-1]: ступени 1 и N в V0 заняты
	// конденсатором (total) и ребойлером (kettle), фид на них физически
	// бессмыслен. См. doc/adr/0002-feed-stage-range.md.
	seenFeed := make(map[int]bool, len(col.FeedStages))
	for i, fs := range col.FeedStages {
		if fs < 2 || fs > col.N-1 {
			v.addf("column.feed_stages[%d]=%d: stage out of range, must be in [2..N-1] (N=%d); stages 1 and N are condenser/reboiler (see doc/adr/0002-feed-stage-range.md)", i, fs, col.N)
		}
		if seenFeed[fs] {
			v.addf("column.feed_stages[%d]=%d: duplicate feed stage", i, fs)
		}
		seenFeed[fs] = true
	}

	for i, sd := range col.SideDraws {
		if sd.Stage < 1 || sd.Stage > col.N {
			v.addf("column.side_draws[%d].stage=%d: stage out of range, must be in [1..N] (N=%d)", i, sd.Stage, col.N)
		}
	}
}
