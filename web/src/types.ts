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

export interface FlowsheetModel {
  format: string;
  name: string;
  description?: string;
  components: string[];
  property_method: string;
  units: UnitModel[];
  streams: StreamModel[];
  reactions?: Reaction[];
  design_specs?: DesignSpec[];
  sensitivity?: SensitivityStudy[];
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
