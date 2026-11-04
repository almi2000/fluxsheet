export type Value = number | string;

export interface ParamSchema {
  name: string;
  label: string;
  quantity: string | null;
  default: Value | null;
  help: string;
  min: number | null;
  max: number | null;
  options: string[] | null;
  requires?: string | null; // shown only when this other parameter is set
  tab?: string | null; // short label among the alternatives of a group
}

export interface ResultField {
  name: string;
  label: string;
  quantity: string | null;
}

export interface ProfileSeries {
  label: string;
  quantity: string | null;
  values: number[];
}

export interface Profile {
  label: string;
  x_label: string;
  x: number[];
  series: ProfileSeries[];
}

/** Scalar results by name, plus optional profiles (stage or length profiles). */
export type UnitResult = Record<string, unknown> & { type?: string; profiles?: Record<string, Profile> };

export interface PortSchema {
  name: string;
  label: string;
  multiple: boolean;
  separate?: boolean;
}

export interface UnitType {
  type: string;
  label: string;
  category: string;
  description: string;
  inlets: PortSchema[];
  outlets: PortSchema[];
  params: ParamSchema[];
  exclusive_groups: string[][];
  required_groups: string[][];
  result_fields: ResultField[];
  dynamic?: DynamicSchema; // has a dynamic model (hold-up)
}

export interface Method {
  key: string;
  label: string;
  description: string;
}

export interface Compound {
  cas: string;
  name: string;
  formula: string;
  mw: number;
  smiles: string;
  simulation_ready: boolean;
}

export type Params = Record<string, unknown>;

export interface UnitModel {
  id: string;
  type: string;
  params?: Params;
  property_method?: string; // overrides the flowsheet's property_method for this unit
  ui?: { x: number; y: number };
}

export interface StreamModel {
  id: string;
  from: string;
  to?: string;
}

/** Power law r = k0 exp(-Ea/RT) prod(c_i^order_i), per unit volume; SI units. */
export interface RateLaw {
  k0: number;
  Ea?: Value;
  orders?: Record<string, number>; // default: the stoichiometric coefficients
  T_ref?: Value; // if set, k0 is the rate constant at T_ref
}

export interface Kinetics extends RateLaw {
  basis?: string; // "concentration" | "partial_pressure"
  reverse?: RateLaw;
  equilibrium_limited?: boolean;
}

/** ln K = A + B/T + C ln T + D T */
export interface Keq {
  A?: number;
  B?: number;
  C?: number;
  D?: number;
}

export interface Reaction {
  id: string;
  stoichiometry: Record<string, number>; // negative for reactants
  kinetics?: Kinetics;
  keq?: Keq;
}

/** A heating or cooling medium (utilities.py); values may carry units, e.g. "25 degC". */
export interface Utility {
  id: string;
  type: "fluid" | "steam";
  T_in?: Value;
  T_out?: Value;
  cp?: Value;
  P?: Value;
}

/** A user-defined compound; missing constants are estimated (custom_compounds.py). */
export interface CustomCompound {
  name: string;
  formula: string;
  smiles?: string;
  joback_groups?: Record<string, number>;
  unifac_groups?: Record<string, number>;
  constants?: Record<string, number | string>;
  vapor_pressure?: string;
  antoine?: Record<string, number | string>;
}

export interface FlowsheetModel {
  format: string;
  name: string;
  description?: string;
  components: string[];
  property_method: string;
  henry_components?: string[]; // dissolved gases of activity methods; absent = automatic
  custom_compounds?: CustomCompound[];
  units: UnitModel[];
  streams: StreamModel[];
  reactions?: Reaction[];
  utilities?: Utility[];
  design_specs?: DesignSpec[];
  sensitivity?: SensitivityStudy[];
  calculators?: unknown[]; // calculator blocks (calculator.py); kept as they are
  dynamics?: Dynamics;
}

export interface StreamResult {
  T: number;
  P: number;
  vapor_fraction: number;
  phase: string;
  molar_flow: number;
  mass_flow: number;
  enthalpy: number;
  enthalpy_flow: number;
  component_flows: Record<string, number>;
  mole_fractions: Record<string, number>;
  property_method?: string; // the method of the unit producing the stream
}

export interface SimError {
  message: string;
  where: string | null;
  hint: string | null;
}

export interface SolveResult {
  ok: boolean;
  streams: Record<string, StreamResult>;
  units: Record<string, UnitResult>;
  property_method: string | null;
  property_method_reason: string | null;
  property_crossings?: unknown[]; // streams re-flashed between methods; see PropertyMethods.tsx
  // Where each binary parameter set came from ("ThermoML fit", "ChemSep", "UNIFAC estimate",
  // "none (ideal)"); `method` is set for pairs of a per-unit method other than the default.
  parameter_sources?: { pair: string; source: string; method?: string; kind?: "binary" | "henry" | "estimate" }[];
  tear_streams: string[];
  iterations: number;
  converged: boolean;
  warnings: string[];
  error: SimError | null;
  design_specs?: DesignSpecResult[];
}

export interface Example {
  id: string;
  name: string;
  description: string;
}

// ---- design specs and sensitivity studies (see src/fluxsheet/study_models.py) ----------

/** A computed value: a stream property, or a scalar result of a unit (e.g. its duty). */
export interface VariableRef {
  stream?: string;
  unit?: string;
  property: string;
  component?: string; // for mole_fraction and component_flow
}

/** A unit parameter changed between bounds; `component` picks one entry of e.g. Feed flows. */
export interface VariedParam {
  unit: string;
  param: string;
  component?: string;
  lower: Value;
  upper: Value;
}

export interface DesignSpec {
  id: string;
  target: VariableRef;
  value: Value;
  vary: VariedParam;
  tolerance?: Value;
  max_iterations?: number;
  active?: boolean;
}

export interface SensitivityVary extends VariedParam {
  points?: number;
}

export interface SensitivityStudy {
  id: string;
  vary: SensitivityVary[];
  outputs: VariableRef[];
}

export interface DesignSpecResult {
  id: string;
  converged: boolean;
  target: string;
  quantity: string | null;
  value: number;
  achieved: number;
  vary: string;
  vary_quantity: string | null;
  vary_value: number;
  iterations: number;
}

export interface SensitivityColumn {
  label: string;
  quantity: string | null;
}

export interface SensitivityRow {
  inputs: number[];
  outputs: (number | null)[];
  error: SimError | null;
}

export interface SensitivityResult {
  ok: boolean;
  id?: string;
  inputs?: SensitivityColumn[];
  outputs?: SensitivityColumn[];
  rows?: SensitivityRow[];
  error?: SimError;
}

// ---- dynamic simulation (see src/fluxsheet/dynamics/models.py) --------------------------

/** Settings and results of a unit's dynamic model, from /api/unit-types. */
export interface DynamicSchema {
  params: ParamSchema[];
  result_fields: ResultField[];
}

/** A unit parameter or dynamic setting that a controller or event changes. */
export interface ManipulatedParam {
  unit: string;
  param: string;
  component?: string;
}

export interface Controller {
  id: string;
  measure: VariableRef;
  manipulate: ManipulatedParam;
  setpoint?: Value; // default: the steady-state measurement
  action?: "reverse" | "direct";
  gain?: number;
  integral_time?: Value;
  derivative_time?: Value;
  measurement_lag?: Value;
  op_range: [Value, Value];
  pv_range?: [Value, Value];
  active?: boolean;
}

export interface DynamicEvent {
  time: Value;
  unit?: string;
  param?: string;
  component?: string;
  controller?: string; // changes this controller's setpoint
  value: Value;
  ramp?: Value;
}

export interface Dynamics {
  t_end: Value;
  points?: number;
  units: Record<string, Params>;
  controllers?: Controller[];
  events?: DynamicEvent[];
  rtol?: number;
}

export interface DynamicSeries {
  key: string;
  group: "Streams" | "Units" | "Controllers";
  owner: string;
  label: string;
  quantity: string | null;
  values: number[];
}

export interface DynamicResult {
  ok: boolean;
  time: number[];
  series: DynamicSeries[];
  warnings: string[];
  error?: SimError;
}
