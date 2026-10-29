// Design specs and sensitivity studies (flowsheet settings). Both name flowsheet values the
// same way: a stream property or a unit result to read (VariableRef), and a unit parameter to
// change between bounds (VariedParam). See src/fluxsheet/specs.py.
import { useEffect, useState } from "react";
import { UNITS, formatQuantity, joinValue, splitValue } from "./quantities";
import "./studies.css";
import type {
  DesignSpec,
  DesignSpecResult,
  SensitivityStudy,
  SensitivityVary,
  UnitType,
  Value,
  VariableRef,
  VariedParam,
} from "./types";

export type UnitInfo = { id: string; schema: UnitType };

/** Mirrors STREAM_PROPERTIES in specs.py. */
export const STREAM_PROPERTIES: { key: string; label: string; quantity: string | null; perComponent?: boolean }[] = [
  { key: "T", label: "Temperature", quantity: "temperature" },
  { key: "P", label: "Pressure", quantity: "pressure" },
  { key: "vapor_fraction", label: "Vapor fraction", quantity: "fraction" },
  { key: "molar_flow", label: "Molar flow", quantity: "molar_flow" },
  { key: "mass_flow", label: "Mass flow", quantity: "mass_flow" },
  { key: "enthalpy_flow", label: "Enthalpy flow", quantity: "power" },
  { key: "mole_fraction", label: "Mole fraction of…", quantity: "fraction", perComponent: true },
  { key: "component_flow", label: "Molar flow of…", quantity: "molar_flow", perComponent: true },
];

const DISCRETE = ["integer", "choice", "component", "string_list"];

/** Parameters a study can vary: numbers, or one entry of a per-component table. */
function continuousParams(schema: UnitType) {
  return schema.params.filter((p) => !DISCRETE.includes(p.quantity ?? ""));
}

function perComponent(schema: UnitType | undefined, param: string): boolean {
  const q = schema?.params.find((p) => p.name === param)?.quantity;
  return q === "composition" || q === "component_fractions";
}

export function variableQuantity(ref: VariableRef, units: UnitInfo[]): string | null {
  if (ref.stream !== undefined) return STREAM_PROPERTIES.find((p) => p.key === ref.property)?.quantity ?? null;
  const schema = units.find((u) => u.id === ref.unit)?.schema;
  return schema?.result_fields.find((f) => f.name === ref.property)?.quantity ?? null;
}

function paramQuantity(vary: VariedParam, units: UnitInfo[]): string | null {
  const q = units.find((u) => u.id === vary.unit)?.schema.params.find((p) => p.name === vary.param)?.quantity ?? null;
  return q === "composition" ? "molar_flow" : q === "component_fractions" ? "fraction" : q;
}

function variableLabel(ref: VariableRef, units: UnitInfo[]): string {
  if (ref.stream !== undefined) {
    const p = STREAM_PROPERTIES.find((s) => s.key === ref.property);
    const name = p?.perComponent ? `${ref.component ?? "?"} ${p.label.replace(" of…", "").toLowerCase()}` : p?.label.toLowerCase();
    return `${ref.stream} ${name ?? ref.property}`;
  }
  const field = units.find((u) => u.id === ref.unit)?.schema.result_fields.find((f) => f.name === ref.property);
  return `${ref.unit} ${field?.label ?? ref.property}`;
}

function paramLabel(vary: VariedParam, units: UnitInfo[]): string {
  const p = units.find((u) => u.id === vary.unit)?.schema.params.find((x) => x.name === vary.param);
  return `${vary.unit} ${p?.label ?? vary.param}${vary.component ? ` of ${vary.component}` : ""}`;
}

function nextId(prefix: string, taken: string[]) {
  let i = 1;
  while (taken.includes(`${prefix}${i}`)) i++;
  return `${prefix}${i}`;
}

// ---- inputs ------------------------------------------------------------------------------

/** A value with a unit selector, or a plain number for fractions and dimensionless values. */
function Qty({
  quantity,
  value,
  placeholder,
  label,
  onChange,
}: {
  quantity: string | null;
  value: Value | undefined;
  placeholder?: string;
  label: string;
  onChange: (v: Value) => void;
}) {
  const units = quantity ? UNITS[quantity] : undefined;
  const { num, unit } = splitValue(value ?? "", units ? quantity! : "");
  const [text, setText] = useState(num);
  useEffect(() => setText(num), [num]);
  if (!units)
    return (
      <input
        type="text"
        inputMode="decimal"
        aria-label={label}
        placeholder={placeholder}
        value={text}
        onChange={(e) => {
          setText(e.target.value);
          const t = e.target.value.trim();
          onChange(t !== "" && Number.isFinite(Number(t)) ? Number(t) : t);
        }}
      />
    );
  return (
    <div className="quantity">
      <input
        type="text"
        inputMode="decimal"
        aria-label={label}
        placeholder={placeholder}
        value={num}
        onChange={(e) => onChange(joinValue(e.target.value, unit) ?? "")}
      />
      <select aria-label={`${label} unit`} value={unit} onChange={(e) => onChange(joinValue(num, e.target.value) ?? "")}>
        {units.map((u) => (
          <option key={u.value} value={u.value}>
            {u.label}
          </option>
        ))}
        {!units.some((u) => u.value === unit) && <option value={unit}>{unit}</option>}
      </select>
    </div>
  );
}

function VariablePicker({
  value,
  units,
  streams,
  components,
  onChange,
}: {
  value: VariableRef;
  units: UnitInfo[];
  streams: string[];
  components: string[];
  onChange: (v: VariableRef) => void;
}) {
  const owner = value.stream !== undefined ? `stream:${value.stream}` : `unit:${value.unit}`;
  const unitsWithResults = units.filter((u) => u.schema.result_fields.length > 0);
  const schema = units.find((u) => u.id === value.unit)?.schema;
  const prop = STREAM_PROPERTIES.find((p) => p.key === value.property);
  return (
    <div className="study-pick">
      <select
        aria-label="Stream or unit"
        value={owner}
        onChange={(e) => {
          const [kind, id] = e.target.value.split(/:(.*)/s);
          if (kind === "stream") onChange({ stream: id, property: value.stream !== undefined ? value.property : "T", component: value.component });
          else {
            const first = units.find((u) => u.id === id)?.schema.result_fields[0]?.name ?? "duty";
            onChange({ unit: id, property: first });
          }
        }}
      >
        <optgroup label="Streams">
          {streams.map((s) => (
            <option key={s} value={`stream:${s}`}>
              {s}
            </option>
          ))}
        </optgroup>
        <optgroup label="Unit results">
          {unitsWithResults.map((u) => (
            <option key={u.id} value={`unit:${u.id}`}>
              {u.id}
            </option>
          ))}
        </optgroup>
        {!streams.includes(value.stream ?? "") && !units.some((u) => u.id === value.unit) && (
          <option value={owner}>{value.stream ?? value.unit} (missing)</option>
        )}
      </select>
      <select
        aria-label="Property"
        value={value.property}
        onChange={(e) => {
          const p = STREAM_PROPERTIES.find((x) => x.key === e.target.value);
          onChange({ ...value, property: e.target.value, component: p?.perComponent ? value.component ?? components[0] : undefined });
        }}
      >
        {value.stream !== undefined
          ? STREAM_PROPERTIES.map((p) => (
              <option key={p.key} value={p.key}>
                {p.label}
              </option>
            ))
          : (schema?.result_fields ?? []).map((f) => (
              <option key={f.name} value={f.name}>
                {f.label}
              </option>
            ))}
      </select>
      {value.stream !== undefined && prop?.perComponent && (
        <select aria-label="Component" value={value.component ?? ""} onChange={(e) => onChange({ ...value, component: e.target.value })}>
          {components.map((c) => (
            <option key={c} value={c}>
              {c}
            </option>
          ))}
        </select>
      )}
    </div>
  );
}

function ParamPicker<T extends VariedParam>({
  value,
  units,
  components,
  onChange,
}: {
  value: T;
  units: UnitInfo[];
  components: string[];
  onChange: (v: T) => void;
}) {
  const candidates = units.filter((u) => continuousParams(u.schema).length > 0);
  const schema = units.find((u) => u.id === value.unit)?.schema;
  const reset = (unit: string, param: string): T => ({
    ...value,
    unit,
    param,
    component: perComponent(units.find((u) => u.id === unit)?.schema, param) ? components[0] : undefined,
    lower: "",
    upper: "",
  });
  return (
    <div className="study-pick">
      <select
        aria-label="Unit"
        value={value.unit}
        onChange={(e) => {
          const u = candidates.find((c) => c.id === e.target.value);
          if (u) onChange(reset(u.id, continuousParams(u.schema)[0].name));
        }}
      >
        {candidates.map((u) => (
          <option key={u.id} value={u.id}>
            {u.id}
          </option>
        ))}
        {!schema && <option value={value.unit}>{value.unit} (missing)</option>}
      </select>
      <select aria-label="Parameter" value={value.param} onChange={(e) => onChange(reset(value.unit, e.target.value))}>
        {(schema ? continuousParams(schema) : []).map((p) => (
          <option key={p.name} value={p.name}>
            {p.label}
          </option>
        ))}
      </select>
      {perComponent(schema, value.param) && (
        <select aria-label="Component" value={value.component ?? ""} onChange={(e) => onChange({ ...value, component: e.target.value })}>
          {components.map((c) => (
            <option key={c} value={c}>
              {c}
            </option>
          ))}
        </select>
      )}
    </div>
  );
}

function Bounds<T extends VariedParam>({ value, units, onChange }: { value: T; units: UnitInfo[]; onChange: (v: T) => void }) {
  const q = paramQuantity(value, units);
  return (
    <div className="study-bounds">
      <label>
        <span className="field-label">From</span>
        <Qty quantity={q} label="Lower bound" value={value.lower} onChange={(lower) => onChange({ ...value, lower })} />
      </label>
      <label>
        <span className="field-label">To</span>
        <Qty quantity={q} label="Upper bound" value={value.upper} onChange={(upper) => onChange({ ...value, upper })} />
      </label>
    </div>
  );
}

function Card({ title, subtitle, children }: { title: string; subtitle: string; children: React.ReactNode }) {
  const [open, setOpen] = useState(false);
  return (
    <li className="study-card">
      <button className="study-summary" onClick={() => setOpen(!open)} aria-expanded={open}>
        <b>{title}</b>
        <span className="muted small">{subtitle}</span>
      </button>
      {open && <div className="study-body">{children}</div>}
    </li>
  );
}

function IdField({ value, taken, onChange }: { value: string; taken: string[]; onChange: (id: string) => void }) {
  const [draft, setDraft] = useState(value);
  useEffect(() => setDraft(value), [value]);
  const clash = draft !== value && taken.includes(draft.trim());
  return (
    <label className="field">
      <span className="field-label">Id</span>
      <input
        value={draft}
        onChange={(e) => setDraft(e.target.value)}
        onBlur={() => (draft.trim() && !clash ? onChange(draft.trim()) : setDraft(value))}
      />
      {clash && <span className="field-error">{draft} already exists</span>}
    </label>
  );
}

// ---- design specs ------------------------------------------------------------------------

type Shared = { units: UnitInfo[]; streams: string[]; components: string[] };

function defaultTarget({ units, streams }: Shared): VariableRef {
  if (streams.length) return { stream: streams[streams.length - 1], property: "T" };
  const u = units.find((x) => x.schema.result_fields.length);
  return { unit: u?.id ?? "", property: u?.schema.result_fields[0].name ?? "duty" };
}

function defaultVary({ units, components }: Shared): VariedParam {
  const u = units.find((x) => x.schema.type !== "Feed" && continuousParams(x.schema).length) ??
    units.find((x) => continuousParams(x.schema).length);
  const p = u ? continuousParams(u.schema)[0].name : "";
  return { unit: u?.id ?? "", param: p, component: perComponent(u?.schema, p) ? components[0] : undefined, lower: "", upper: "" };
}

export function DesignSpecsEditor({
  specs,
  results,
  onChange,
  ...shared
}: Shared & { specs: DesignSpec[]; results?: DesignSpecResult[]; onChange: (s: DesignSpec[]) => void }) {
  const { units } = shared;
  const set = (i: number, spec: DesignSpec) => onChange(specs.map((s, j) => (j === i ? spec : s)));
  return (
    <div className="studies">
      <p className="muted small">
        Reach a target value by adjusting one unit parameter between bounds. Every Run solves the flowsheet
        repeatedly until the target holds.
      </p>
      <ul className="study-list">
        {specs.map((spec, i) => {
          const q = variableQuantity(spec.target, units);
          const res = results?.find((r) => r.id === spec.id);
          return (
            <Card
              key={i}
              title={`${spec.id}${spec.active === false ? " (off)" : ""}`}
              subtitle={`${variableLabel(spec.target, units)} = ${spec.value === "" ? "?" : spec.value} by ${paramLabel(spec.vary, units)}`}
            >
              {res && (
                <p className="callout">
                  Met: {formatQuantity(res.achieved, res.quantity)} with {res.vary} ={" "}
                  {formatQuantity(res.vary_value, res.vary_quantity)} ({res.iterations} solves)
                </p>
              )}
              <IdField value={spec.id} taken={specs.map((s) => s.id)} onChange={(id) => set(i, { ...spec, id })} />
              <label className="check">
                <input
                  type="checkbox"
                  checked={spec.active !== false}
                  onChange={(e) => set(i, { ...spec, active: e.target.checked ? undefined : false })}
                />
                Active
              </label>
              <div className="field">
                <span className="field-label">Target</span>
                <VariablePicker
                  {...shared}
                  value={spec.target}
                  onChange={(target) =>
                    set(i, { ...spec, target, value: variableQuantity(target, units) === q ? spec.value : "" })
                  }
                />
              </div>
              <label className="field">
                <span className="field-label">Value</span>
                <Qty quantity={q} label="Target value" value={spec.value} onChange={(value) => set(i, { ...spec, value })} />
              </label>
              <div className="field">
                <span className="field-label">Vary</span>
                <ParamPicker {...shared} value={spec.vary} onChange={(vary) => set(i, { ...spec, vary })} />
              </div>
              <Bounds value={spec.vary} units={units} onChange={(vary) => set(i, { ...spec, vary })} />
              <button className="danger" onClick={() => onChange(specs.filter((_, j) => j !== i))}>
                Delete design spec
              </button>
            </Card>
          );
        })}
      </ul>
      <button
        disabled={units.length === 0}
        onClick={() =>
          onChange([...specs, { id: nextId("DS", specs.map((s) => s.id)), target: defaultTarget(shared), value: "", vary: defaultVary(shared) }])
        }
      >
        + Design spec
      </button>
    </div>
  );
}

// ---- sensitivity studies -----------------------------------------------------------------

export function SensitivityEditor({
  studies,
  onChange,
  ...shared
}: Shared & { studies: SensitivityStudy[]; onChange: (s: SensitivityStudy[]) => void }) {
  const { units } = shared;
  const set = (i: number, study: SensitivityStudy) => onChange(studies.map((s, j) => (j === i ? study : s)));
  return (
    <div className="studies">
      <p className="muted small">
        Solve the flowsheet over a range of one or two parameters and record outputs. Run a study from the
        Sensitivity tab of the results.
      </p>
      <ul className="study-list">
        {studies.map((study, i) => {
          const setVary = (k: number, v: SensitivityVary) => set(i, { ...study, vary: study.vary.map((x, j) => (j === k ? v : x)) });
          const setOutput = (k: number, o: VariableRef) =>
            set(i, { ...study, outputs: study.outputs.map((x, j) => (j === k ? o : x)) });
          return (
            <Card
              key={i}
              title={study.id}
              subtitle={`${study.vary.map((v) => paramLabel(v, units)).join(" × ")} → ${study.outputs.length} output${study.outputs.length === 1 ? "" : "s"}`}
            >
              <IdField value={study.id} taken={studies.map((s) => s.id)} onChange={(id) => set(i, { ...study, id })} />
              {study.vary.map((v, k) => (
                <div key={k} className="field">
                  <span className="field-label">
                    {k === 0 ? "Vary" : "and vary"}
                    {k > 0 && (
                      <button className="link remove" onClick={() => set(i, { ...study, vary: study.vary.slice(0, 1) })}>
                        remove
                      </button>
                    )}
                  </span>
                  <ParamPicker {...shared} value={v} onChange={(x) => setVary(k, x)} />
                  <Bounds value={v} units={units} onChange={(x) => setVary(k, x)} />
                  <label className="study-points">
                    <span className="field-label">Points</span>
                    <input
                      type="number"
                      min={2}
                      max={101}
                      value={v.points ?? 11}
                      onChange={(e) => setVary(k, { ...v, points: Math.max(2, Math.min(101, Math.round(Number(e.target.value) || 2))) })}
                    />
                  </label>
                </div>
              ))}
              {study.vary.length === 1 && (
                <button
                  className="link"
                  onClick={() => set(i, { ...study, vary: [...study.vary, { ...defaultVary(shared), points: 3 }] })}
                >
                  + Second parameter
                </button>
              )}
              <div className="field">
                <span className="field-label">Outputs</span>
                {study.outputs.map((o, k) => (
                  <div key={k} className="study-output">
                    <VariablePicker {...shared} value={o} onChange={(x) => setOutput(k, x)} />
                    <button
                      aria-label="Remove output"
                      disabled={study.outputs.length === 1}
                      onClick={() => set(i, { ...study, outputs: study.outputs.filter((_, j) => j !== k) })}
                    >
                      ×
                    </button>
                  </div>
                ))}
                <button className="link" onClick={() => set(i, { ...study, outputs: [...study.outputs, defaultTarget(shared)] })}>
                  + Output
                </button>
              </div>
              <button className="danger" onClick={() => onChange(studies.filter((_, j) => j !== i))}>
                Delete study
              </button>
            </Card>
          );
        })}
      </ul>
      <button
        disabled={units.length === 0}
        onClick={() =>
          onChange([
            ...studies,
            {
              id: nextId("SENS", studies.map((s) => s.id)),
              vary: [{ ...defaultVary(shared), points: 11 }],
              outputs: [defaultTarget(shared)],
            },
          ])
        }
      >
        + Sensitivity study
      </button>
    </div>
  );
}

/** Both editors as sections of the flowsheet settings. */
export function StudiesSettings({
  specs,
  studies,
  results,
  onChange,
  ...shared
}: Shared & {
  specs: DesignSpec[];
  studies: SensitivityStudy[];
  results?: DesignSpecResult[];
  onChange: (patch: { design_specs?: DesignSpec[]; sensitivity?: SensitivityStudy[] }) => void;
}) {
  return (
    <>
      <h3>Design specs</h3>
      <DesignSpecsEditor {...shared} specs={specs} results={results} onChange={(s) => onChange({ design_specs: s })} />
      <h3>Sensitivity</h3>
      <SensitivityEditor {...shared} studies={studies} onChange={(s) => onChange({ sensitivity: s })} />
    </>
  );
}
