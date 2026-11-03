// Costing settings and optimizations (flowsheet settings). An optimization names variables
// (flowsheet values, costing totals or parameter values), minimizes or maximizes an expression
// over them by varying unit parameters, subject to constraints. See src/fluxsheet/optimization.py
// and costing.py.
import { useEffect, useState } from "react";
import {
  Bounds,
  Card,
  continuousParams,
  IdField,
  nextId,
  ParamPicker,
  perComponent,
  Qty,
  type UnitInfo,
  VariablePicker,
} from "./StudiesEditor";
import "./studies.css";
import type { Constraint, Costing, Optimization, OptimizationVary, ParamRef, VariableRef } from "./types";

type Shared = { units: UnitInfo[]; streams: string[]; components: string[] };

const isParamRef = (r: VariableRef | ParamRef): r is ParamRef => "param" in r;

/** A number input that keeps what is typed until it is a number (or empty: undefined). */
function Num({
  value,
  label,
  placeholder,
  onChange,
}: {
  value: number | undefined;
  label: string;
  placeholder?: string;
  onChange: (v: number | undefined) => void;
}) {
  const [text, setText] = useState(value === undefined ? "" : String(value));
  useEffect(() => setText(value === undefined ? "" : String(value)), [value]);
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
        if (t === "") onChange(undefined);
        else if (Number.isFinite(Number(t))) onChange(Number(t));
      }}
    />
  );
}

function NumField(props: { label: string; value: number | undefined; placeholder: string; onChange: (v: number | undefined) => void }) {
  return (
    <label>
      <span className="field-label">{props.label}</span>
      <Num {...props} />
    </label>
  );
}

// ---- costing -------------------------------------------------------------------------------

export function CostingEditor({ costing, onChange }: { costing?: Costing; onChange: (c: Costing | undefined) => void }) {
  if (!costing)
    return (
      <div className="studies">
        <p className="muted small">
          Estimate purchased and installed equipment costs (Turton correlations), utility costs and the total annualized
          cost (TAC) after every Run. Objectives and studies can then use <code>tac</code>.
        </p>
        <button onClick={() => onChange({})}>+ Costing</button>
      </div>
    );
  const set = (patch: Partial<Costing>) => onChange({ ...costing, ...patch });
  return (
    <div className="studies">
      <p className="muted small">
        Study estimate (±30 %). Empty prices pick steam, cooling water or refrigerant by temperature; prices in $/GJ.
      </p>
      <div className="dyn-row">
        <NumField label="Cost index (CEPCI)" value={costing.cepci} placeholder="800" onChange={(cepci) => set({ cepci })} />
        <NumField label="Hours per year" value={costing.hours_per_year} placeholder="8000" onChange={(hours_per_year) => set({ hours_per_year })} />
        <NumField label="Interest rate" value={costing.interest_rate} placeholder="0.1" onChange={(interest_rate) => set({ interest_rate })} />
        <NumField label="Plant life [yr]" value={costing.plant_life} placeholder="10" onChange={(plant_life) => set({ plant_life })} />
        <NumField label="Electricity [$/kWh]" value={costing.electricity_price} placeholder="0.07" onChange={(electricity_price) => set({ electricity_price })} />
        <NumField label="Tray efficiency" value={costing.tray_efficiency} placeholder="0.7" onChange={(tray_efficiency) => set({ tray_efficiency })} />
        <NumField label="Heating [$/GJ]" value={costing.heating_price} placeholder="by temperature" onChange={(heating_price) => set({ heating_price })} />
        <NumField label="Cooling [$/GJ]" value={costing.cooling_price} placeholder="by temperature" onChange={(cooling_price) => set({ cooling_price })} />
        <label>
          <span className="field-label">Material</span>
          <select aria-label="Material" value={costing.material ?? "CS"} onChange={(e) => set({ material: e.target.value as "CS" | "SS" })}>
            <option value="CS">Carbon steel</option>
            <option value="SS">Stainless steel</option>
          </select>
        </label>
        <label>
          <span className="field-label">Exchanger U</span>
          <Qty quantity="heat_transfer_coefficient" label="Exchanger U" value={costing.U} placeholder="500" onChange={(U) => set({ U: U === "" ? undefined : U })} />
        </label>
      </div>
      <button className="danger" onClick={() => onChange(undefined)}>
        Remove costing
      </button>
    </div>
  );
}

// ---- optimizations ---------------------------------------------------------------------------

function ParamRefPicker({ value, units, components, onChange }: { value: ParamRef } & Omit<Shared, "streams"> & { onChange: (v: ParamRef) => void }) {
  const candidates = units.filter((u) => continuousParams(u.schema).length > 0);
  const schema = units.find((u) => u.id === value.unit)?.schema;
  const pick = (unit: string, param: string): ParamRef => ({
    unit,
    param,
    component: perComponent(units.find((u) => u.id === unit)?.schema, param) ? components[0] : undefined,
  });
  return (
    <div className="study-pick">
      <select
        aria-label="Unit"
        value={value.unit}
        onChange={(e) => {
          const u = candidates.find((c) => c.id === e.target.value);
          if (u) onChange(pick(u.id, continuousParams(u.schema)[0].name));
        }}
      >
        {candidates.map((u) => (
          <option key={u.id} value={u.id}>
            {u.id}
          </option>
        ))}
      </select>
      <select aria-label="Parameter" value={value.param} onChange={(e) => onChange(pick(value.unit, e.target.value))}>
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

/** One named variable: a flowsheet value (result) or the value of a parameter. */
function VariableRow({
  name,
  value,
  taken,
  costing,
  shared,
  onRename,
  onChange,
  onRemove,
}: {
  name: string;
  value: VariableRef | ParamRef;
  taken: string[];
  costing: boolean;
  shared: Shared;
  onRename: (n: string) => void;
  onChange: (v: VariableRef | ParamRef) => void;
  onRemove: () => void;
}) {
  const [draft, setDraft] = useState(name);
  useEffect(() => setDraft(name), [name]);
  const bad = draft !== name && (!/^[A-Za-z_]\w*$/.test(draft) || taken.includes(draft));
  const firstParam = (): ParamRef => {
    const u = shared.units.find((x) => continuousParams(x.schema).length);
    return { unit: u?.id ?? "", param: u ? continuousParams(u.schema)[0].name : "" };
  };
  return (
    <div className="opt-var">
      <input
        aria-label="Variable name"
        className={bad ? "invalid" : ""}
        value={draft}
        onChange={(e) => setDraft(e.target.value)}
        onBlur={() => (!bad && draft !== name ? onRename(draft) : setDraft(name))}
      />
      <select
        aria-label="Variable kind"
        value={isParamRef(value) ? "param" : "result"}
        onChange={(e) => onChange(e.target.value === "param" ? firstParam() : costing ? { property: "tac" } : { stream: shared.streams[0] ?? "", property: "T" })}
      >
        <option value="result">result</option>
        <option value="param">parameter</option>
      </select>
      <div className="opt-var-pick">
        {isParamRef(value) ? (
          <ParamRefPicker value={value} units={shared.units} components={shared.components} onChange={onChange} />
        ) : (
          <VariablePicker {...shared} costing={costing} value={value} onChange={onChange} />
        )}
      </div>
      <button aria-label={`Remove variable ${name}`} onClick={onRemove}>
        ×
      </button>
    </div>
  );
}

function freeName(taken: string[]) {
  for (const n of ["x", "y", "z", "a", "b", "c"]) if (!taken.includes(n)) return n;
  return nextId("v", taken);
}

function defaultVary({ units, components }: Shared): OptimizationVary {
  const u = units.find((x) => x.schema.type !== "Feed" && continuousParams(x.schema).length) ?? units.find((x) => continuousParams(x.schema).length);
  const p = u ? continuousParams(u.schema)[0].name : "";
  return { unit: u?.id ?? "", param: p, component: perComponent(u?.schema, p) ? components[0] : undefined, lower: "", upper: "" };
}

export function OptimizationsEditor({
  optimizations,
  costing,
  onChange,
  ...shared
}: Shared & { optimizations: Optimization[]; costing: boolean; onChange: (o: Optimization[]) => void }) {
  const { units } = shared;
  const set = (i: number, o: Optimization) => onChange(optimizations.map((x, j) => (j === i ? o : x)));
  return (
    <div className="studies">
      <p className="muted small">
        Minimize or maximize an expression of named variables (SI units, $ for costs) by varying unit parameters. Run it
        from the Optimization tab of the results.
      </p>
      <ul className="study-list">
        {optimizations.map((opt, i) => {
          const names = Object.keys(opt.variables);
          const setVars = (variables: Optimization["variables"]) => set(i, { ...opt, variables });
          const setVary = (k: number, v: OptimizationVary) => set(i, { ...opt, vary: opt.vary.map((x, j) => (j === k ? v : x)) });
          const cons = opt.constraints ?? [];
          const setCon = (k: number, c: Constraint) => set(i, { ...opt, constraints: cons.map((x, j) => (j === k ? c : x)) });
          return (
            <Card key={i} title={opt.id} subtitle={`${opt.sense ?? "minimize"} ${opt.objective || "?"} by ${opt.vary.length} parameter${opt.vary.length === 1 ? "" : "s"}`}>
              <IdField value={opt.id} taken={optimizations.map((o) => o.id)} onChange={(id) => set(i, { ...opt, id })} />
              <div className="field">
                <span className="field-label">Variables</span>
                {names.map((n) => (
                  <VariableRow
                    key={n}
                    name={n}
                    value={opt.variables[n]}
                    taken={names}
                    costing={costing}
                    shared={shared}
                    onRename={(m) => setVars(Object.fromEntries(Object.entries(opt.variables).map(([k, v]) => [k === n ? m : k, v])))}
                    onChange={(v) => setVars({ ...opt.variables, [n]: v })}
                    onRemove={() => setVars(Object.fromEntries(Object.entries(opt.variables).filter(([k]) => k !== n)))}
                  />
                ))}
                <button
                  className="link"
                  onClick={() => setVars({ ...opt.variables, [freeName(names)]: costing ? { property: "tac" } : { stream: shared.streams[0] ?? "", property: "T" } })}
                >
                  + Variable
                </button>
              </div>
              <div className="field">
                <span className="field-label">Objective</span>
                <div className="opt-objective">
                  <select aria-label="Sense" value={opt.sense ?? "minimize"} onChange={(e) => set(i, { ...opt, sense: e.target.value as Optimization["sense"] })}>
                    <option value="minimize">minimize</option>
                    <option value="maximize">maximize</option>
                  </select>
                  <input aria-label="Objective" value={opt.objective} placeholder="e.g. tac" onChange={(e) => set(i, { ...opt, objective: e.target.value })} />
                </div>
              </div>
              {opt.vary.map((v, k) => (
                <div key={k} className="field">
                  <span className="field-label">
                    {k === 0 ? "Vary" : "and vary"}
                    {opt.vary.length > 1 && (
                      <button className="link remove" onClick={() => set(i, { ...opt, vary: opt.vary.filter((_, j) => j !== k) })}>
                        remove
                      </button>
                    )}
                  </span>
                  <ParamPicker {...shared} value={v} onChange={(x) => setVary(k, x)} />
                  <Bounds value={v} units={units} onChange={(x) => setVary(k, x)} />
                </div>
              ))}
              <button className="link" onClick={() => set(i, { ...opt, vary: [...opt.vary, defaultVary(shared)] })}>
                + Parameter to vary
              </button>
              <div className="field">
                <span className="field-label">Constraints</span>
                {cons.map((c, k) => (
                  <div key={k} className="opt-constraint">
                    <input aria-label="Constraint expression" value={c.expression} onChange={(e) => setCon(k, { ...c, expression: e.target.value })} />
                    <select aria-label="Relation" value={c.relation} onChange={(e) => setCon(k, { ...c, relation: e.target.value as Constraint["relation"] })}>
                      <option value="<=">≤</option>
                      <option value=">=">≥</option>
                      <option value="=">=</option>
                    </select>
                    <Num label="Constraint value" value={c.value} onChange={(value) => setCon(k, { ...c, value: value ?? 0 })} />
                    <button aria-label="Remove constraint" onClick={() => set(i, { ...opt, constraints: cons.filter((_, j) => j !== k) })}>
                      ×
                    </button>
                  </div>
                ))}
                <button className="link" onClick={() => set(i, { ...opt, constraints: [...cons, { expression: names[0] ?? "", relation: "<=", value: 0 }] })}>
                  + Constraint
                </button>
              </div>
              <div className="dyn-row">
                <label>
                  <span className="field-label">Method</span>
                  <select aria-label="Method" value={opt.method ?? "auto"} onChange={(e) => set(i, { ...opt, method: e.target.value as Optimization["method"] })}>
                    <option value="auto">SLSQP, COBYLA fallback</option>
                    <option value="slsqp">SLSQP</option>
                    <option value="cobyla">COBYLA</option>
                  </select>
                </label>
                <NumField label="Max. iterations" value={opt.max_iterations} placeholder="50" onChange={(max_iterations) => set(i, { ...opt, max_iterations })} />
              </div>
              <button className="danger" onClick={() => onChange(optimizations.filter((_, j) => j !== i))}>
                Delete optimization
              </button>
            </Card>
          );
        })}
      </ul>
      <button
        disabled={units.length === 0}
        onClick={() =>
          onChange([
            ...optimizations,
            {
              id: nextId("OPT", optimizations.map((o) => o.id)),
              variables: costing ? { tac: { property: "tac" } } : {},
              objective: costing ? "tac" : "",
              sense: "minimize",
              vary: [defaultVary(shared)],
            },
          ])
        }
      >
        + Optimization
      </button>
    </div>
  );
}

/** Both editors as sections of the flowsheet settings. */
export function OptimizationSettings({
  costing,
  optimizations,
  onChange,
  ...shared
}: Shared & {
  costing?: Costing;
  optimizations: Optimization[];
  onChange: (patch: { costing?: Costing; optimizations?: Optimization[] }) => void;
}) {
  return (
    <>
      <h3>Costing</h3>
      <CostingEditor costing={costing} onChange={(c) => onChange({ costing: c })} />
      <h3>Optimization</h3>
      <OptimizationsEditor {...shared} optimizations={optimizations} costing={costing !== undefined} onChange={(o) => onChange({ optimizations: o })} />
    </>
  );
}
