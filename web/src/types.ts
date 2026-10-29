export type Value = number | string;

export interface ParamSchema {
  name: string;
  label: string;
  quantity: string | null;
  default: Value | null;
  help: string;
  min: number | null;
  max: number | null;
}

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

export interface FlowsheetModel {
  format: string;
  name: string;
  description?: string;
  components: string[];
  property_method: string;
  units: UnitModel[];
  streams: StreamModel[];
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
  units: Record<string, Record<string, number | string>>;
  property_method: string | null;
  property_method_reason: string | null;
  tear_streams: string[];
  iterations: number;
  converged: boolean;
  warnings: string[];
  error: SimError | null;
}

export interface Example {
  id: string;
  name: string;
  description: string;
}
