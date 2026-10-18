// Package contractv1 defines the versioned JSON contracts for the ideal-column
// model. These are wire types, not a solver or a property database.
//
// The structural rules live in examples/schema/v1. Cross-file and physical
// validation rules are specified in examples/v1/README.md. ParseCompilerInput
// and ParseModelInput check JSON structure; ValidateCompilerInput and
// ValidateModelInput check semantics. Decoding into types alone is not validation.
package contractv1

const SchemaVersion = 1

// CompilerInput fixes the structure and properties of the generated model.
type CompilerInput struct {
	SchemaVersion int         `json:"schema_version"`
	Name          string      `json:"name"`
	Description   string      `json:"description,omitempty"`
	Provenance    *Provenance `json:"provenance,omitempty"`
	Components    []Component `json:"components"`
	Thermo        Thermo      `json:"thermo"`
	Column        Column      `json:"column"`
}

type Provenance struct {
	Source string `json:"source,omitempty"`
	Pages  string `json:"pages,omitempty"`
}

type Component struct {
	ID            string     `json:"id"`
	Name          string     `json:"name"`
	MolarMassGMol float64    `json:"molar_mass_g_mol"`
	Properties    Properties `json:"properties"`
}

// Properties is a discriminated union. Source=builtin requires Key/Revision;
// source=inline requires Antoine/Enthalpy/TemperatureRangeK. Mixing branches
// is invalid. A builtin reference must be resolved before code generation.
type Properties struct {
	Source            string            `json:"source"`
	Key               string            `json:"key,omitempty"`
	Revision          string            `json:"revision,omitempty"`
	Antoine           *Antoine          `json:"antoine,omitempty"`
	Enthalpy          *ConstantCP       `json:"enthalpy,omitempty"`
	TemperatureRangeK *TemperatureRange `json:"temperature_range_K,omitempty"`
}

// Antoine uses log10(Psat/mmHg) = A - B/(C + t_C).
type Antoine struct {
	A               float64 `json:"A"`
	B               float64 `json:"B"`
	C               float64 `json:"C"`
	PressureUnit    string  `json:"pressure_unit"`
	TemperatureUnit string  `json:"temperature_unit"`
}

// ConstantCP defines ideal pure-component enthalpies, in J/mol:
// hL = CpL*(T-Tref), hV = deltaHvap_ref + CpV*(T-Tref).
// Tref is in K. The liquid reference enthalpy is zero.
type ConstantCP struct {
	Type                                string  `json:"type"`
	ReferenceTemperatureK               float64 `json:"reference_temperature_K"`
	LiquidCPJMolK                       float64 `json:"liquid_cp_J_mol_K"`
	VaporCPJMolK                        float64 `json:"vapor_cp_J_mol_K"`
	VaporizationEnthalpyAtReferenceJMol float64 `json:"vaporization_enthalpy_at_reference_J_mol"`
}

type TemperatureRange struct {
	Min float64 `json:"min"`
	Max float64 `json:"max"`
}

type Thermo struct {
	PhaseEquilibrium string `json:"phase_equilibrium"`
	Enthalpy         string `json:"enthalpy"`
}

// Column numbers stages from 1 at the top through Stages at the bottom.
// Both condenser and reboiler are included in Stages.
type Column struct {
	Stages    int            `json:"stages"`
	Condenser Condenser      `json:"condenser"`
	Reboiler  Reboiler       `json:"reboiler"`
	Pressure  Pressure       `json:"pressure"`
	Feeds     []FeedLocation `json:"feeds"`
}

type Condenser struct {
	Type string `json:"type"`
}

type Reboiler struct {
	Type string `json:"type"`
}

type Pressure struct {
	Type    string  `json:"type"`
	ValuePa float64 `json:"value_Pa"`
}

type FeedLocation struct {
	ID    string `json:"id"`
	Stage int    `json:"stage"`
}

// ModelInput contains operating parameters for one invocation. Feeds are keyed
// by the IDs declared in CompilerInput.Column.Feeds.
type ModelInput struct {
	SchemaVersion  int             `json:"schema_version"`
	Feeds          map[string]Feed `json:"feeds"`
	Specifications Specifications  `json:"specifications"`
	InitialGuess   *InitialGuess   `json:"initial_guess,omitempty"`
}

type Feed struct {
	Flow        Flow        `json:"flow"`
	Composition Composition `json:"composition"`
	State       FeedState   `json:"state"`
}

type Flow struct {
	Value float64 `json:"value"`
	Unit  string  `json:"unit"`
}

// Composition.Fractions is keyed by component ID, independent of array order.
type Composition struct {
	Basis     string             `json:"basis"`
	Fractions map[string]float64 `json:"fractions"`
}

// FeedState is specified at the feed-stage pressure. Type selects exactly one
// of saturated_liquid (no fields), temperature, or vapor_fraction. Pointers
// distinguish omitted parameters from a vapor fraction of zero.
type FeedState struct {
	Type          string   `json:"type"`
	TemperatureK  *float64 `json:"temperature_K,omitempty"`
	VaporFraction *float64 `json:"vapor_fraction,omitempty"`
}

type Specifications struct {
	Type           string      `json:"type"`
	DistillateFlow Flow        `json:"distillate_flow"`
	RefluxRatio    RefluxRatio `json:"reflux_ratio"`
}

type RefluxRatio struct {
	Value float64 `json:"value"`
	Basis string  `json:"basis"`
}

type InitialGuess struct {
	Profile string `json:"profile"`
}

// ModelOutput uses mol/s, K and W throughout, with mole fractions keyed by ID.
// Products and Profiles are required for a converged result, optional otherwise.
type ModelOutput struct {
	SchemaVersion int              `json:"schema_version"`
	Provenance    OutputProvenance `json:"provenance"`
	Products      *Products        `json:"products,omitempty"`
	Profiles      []StageProfile   `json:"profiles,omitempty"`
	Diagnostics   Diagnostics      `json:"diagnostics"`
}

type OutputProvenance struct {
	Solver        string `json:"solver"`
	CompilerInput string `json:"compiler_input"`
	ModelInput    string `json:"model_input"`
}

type Stream struct {
	FlowMolS      float64            `json:"flow_mol_s"`
	MoleFractions map[string]float64 `json:"mole_fractions"`
}

type Products struct {
	Distillate     Stream  `json:"distillate"`
	Bottoms        Stream  `json:"bottoms"`
	CondenserDutyW float64 `json:"condenser_duty_W"`
	ReboilerDutyW  float64 `json:"reboiler_duty_W"`
}

type StageProfile struct {
	Stage        int     `json:"stage"`
	TemperatureK float64 `json:"temperature_K"`
	Liquid       Stream  `json:"liquid"`
	Vapor        *Stream `json:"vapor"` // null at the total condenser
}

type Diagnostics struct {
	Converged        bool      `json:"converged"`
	Iterations       int       `json:"iterations"`
	ResidualNorm     float64   `json:"residual_norm"`
	MassBalanceError *float64  `json:"mass_balance_error,omitempty"`
	ResidualHistory  []float64 `json:"residual_history,omitempty"`
	Message          string    `json:"message,omitempty"`
}
