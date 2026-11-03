// Dynamic simulation setup (flowsheet settings): run time, which units hold up material
// and how big they are, PID controllers and events (steps and ramps). Run it from the
// Dynamics tab of the results. See src/fluxsheet/dynamics/.
import {
  Card,
  IdField,
  Qty,
  VariablePicker,
  continuousParams,
  nextId,
  perComponent,
  variableLabel,
  type UnitInfo,
} from "./StudiesEditor";
import "./studies.css";
import type { Controller, DynamicEvent, Dynamics, ManipulatedParam, ParamSchema, Params, UnitType } from "./types";

type Shared = { units: UnitInfo[]; streams: string[]; components: string[] };

/** A unit's schema with its dynamic results added, so controllers can measure e.g. a level. */
function withDynamicResults(units: UnitInfo[], dynamic: Record<string, Params>): UnitInfo[] {
  return units.map((u) =>
    u.id in dynamic && u.schema.dynamic
      ? { ...u, schema: { ...u.schema, result_fields: [...u.schema.dynamic.result_fields, ...u.schema.result_fields] } }
      : u,
  );
}

/** Parameters a controller or event can move: the unit's numbers and its dynamic settings. */
function movable(schema: UnitType, dynamic: boolean): ParamSchema[] {
  const own = continuousParams(schema);
  const dyn = dynamic ? (schema.dynamic?.params ?? []).filter((p) => p.quantity !== "choice") : [];
  return [...dyn, ...own.filter((p) => !dyn.some((d) => d.name === p.name))];
}

function quantityOf(m: ManipulatedParam, units: UnitInfo[], dynamic: Record<string, Params>): string | null {
  const schema = units.find((u) => u.id === m.unit)?.schema;
  const q = schema ? movable(schema, m.unit in dynamic).find((p) => p.name === m.param)?.quantity ?? null : null;
  return q === "composition" ? "molar_flow" : q === "component_fractions" ? "fraction" : q;
}

function TargetPicker({
  value,
  units,
  components,
  dynamic,
  onChange,
}: {
  value: ManipulatedParam;
  units: UnitInfo[];
  components: string[];
  dynamic: Record<string, Params>;
  onChange: (v: ManipulatedParam) => void;
}) {
  const candidates = units.filter((u) => movable(u.schema, u.id in dynamic).length > 0);
  const schema = units.find((u) => u.id === value.unit)?.schema;
  const reset = (unit: string, param: string): ManipulatedParam => ({
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
          if (u) onChange(reset(u.id, movable(u.schema, u.id in dynamic)[0].name));
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
        {(schema ? movable(schema, value.unit in dynamic) : []).map((p) => (
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

function SettingInput({ param, value, onChange }: { param: ParamSchema; value: unknown; onChange: (v: unknown) => void }) {
  if (param.quantity === "choice")
    return (
      <select aria-label={param.label} value={String(value ?? param.default ?? "")} onChange={(e) => onChange(e.target.value)}>
        {(param.options ?? []).map((o) => (
          <option key={o} value={o}>
            {o}
          </option>
        ))}
      </select>
    );
  const q = param.quantity === "fraction" ? null : param.quantity;
  return (
    <Qty
      quantity={q}
      label={param.label}
      placeholder={param.default != null ? String(param.default) : "steady state"}
      value={value as string | number | undefined}
      onChange={onChange}
    />
  );
}

function HoldUps({
  units,
  dynamic,
  onChange,
}: {
  units: UnitInfo[];
  dynamic: Record<string, Params>;
  onChange: (d: Record<string, Params>) => void;
}) {
  const eligible = units.filter((u) => u.schema.dynamic);
  if (!eligible.length)
    return <p className="muted small">No unit in this flowsheet has a dynamic model yet (vessels, CSTRs and heaters do).</p>;
  return (
    <div className="studies">
      {eligible.map((u) => {
        const on = u.id in dynamic;
        const settings = dynamic[u.id] ?? {};
        return (
          <div key={u.id} className="dyn-unit">
            <label className="check">
              <input
                type="checkbox"
                checked={on}
                onChange={(e) => {
                  const next = { ...dynamic };
                  if (e.target.checked) next[u.id] = {};
                  else delete next[u.id];
                  onChange(next);
                }}
              />
              <b>{u.id}</b> <span className="muted small">{u.schema.label}</span>
            </label>
            {on && (
              <div className="dyn-settings">
                {u.schema.dynamic!.params.map((p) => (
                  <label key={p.name} className="field" title={p.help}>
                    <span className="field-label">{p.label}</span>
                    <SettingInput
                      param={p}
                      value={settings[p.name]}
                      onChange={(v) => {
                        const next = { ...settings };
                        if (v === "" || v === undefined) delete next[p.name];
                        else next[p.name] = v;
                        onChange({ ...dynamic, [u.id]: next });
                      }}
                    />
                  </label>
                ))}
              </div>
            )}
          </div>
        );
      })}
    </div>
  );
}

function ControllersEditor({
  controllers,
  dynamic,
  onChange,
  ...shared
}: Shared & { controllers: Controller[]; dynamic: Record<string, Params>; onChange: (c: Controller[]) => void }) {
  const units = withDynamicResults(shared.units, dynamic);
  const set = (i: number, c: Controller) => onChange(controllers.map((x, j) => (j === i ? c : x)));
  const firstTarget = (): ManipulatedParam => {
    const u = units.find((x) => x.id in dynamic) ?? units.find((x) => movable(x.schema, false).length);
    const p = u ? movable(u.schema, u.id in dynamic)[0].name : "";
    return { unit: u?.id ?? "", param: p };
  };
  const firstMeasure = () => {
    const u = units.find((x) => x.id in dynamic);
    return u ? { unit: u.id, property: u.schema.result_fields[0].name } : { stream: shared.streams[0] ?? "", property: "T" };
  };
  return (
    <div className="studies">
      <ul className="study-list">
        {controllers.map((c, i) => {
          const pvq = units.length ? variableQuantityOf(c, units) : null;
          const opq = quantityOf(c.manipulate, units, dynamic);
          return (
            <Card
              key={i}
              title={`${c.id}${c.active === false ? " (off)" : ""}`}
              subtitle={`${variableLabel(c.measure, units)} → ${c.manipulate.unit} ${c.manipulate.param}`}
            >
              <IdField value={c.id} taken={controllers.map((x) => x.id)} onChange={(id) => set(i, { ...c, id })} />
              <label className="check">
                <input type="checkbox" checked={c.active !== false} onChange={(e) => set(i, { ...c, active: e.target.checked ? undefined : false })} />
                Active
              </label>
              <div className="field">
                <span className="field-label">Measure</span>
                <VariablePicker {...shared} units={units} value={c.measure} onChange={(measure) => set(i, { ...c, measure, setpoint: undefined, pv_range: undefined })} />
              </div>
              <div className="field">
                <span className="field-label">Manipulate</span>
                <TargetPicker
                  value={c.manipulate}
                  units={units}
                  components={shared.components}
                  dynamic={dynamic}
                  onChange={(manipulate) => set(i, { ...c, manipulate, op_range: ["", ""] })}
                />
              </div>
              <div className="dyn-row">
                <label>
                  <span className="field-label">Output from</span>
                  <Qty quantity={opq} label="Output low" value={c.op_range[0]} onChange={(v) => set(i, { ...c, op_range: [v, c.op_range[1]] })} />
                </label>
                <label>
                  <span className="field-label">to</span>
                  <Qty quantity={opq} label="Output high" value={c.op_range[1]} onChange={(v) => set(i, { ...c, op_range: [c.op_range[0], v] })} />
                </label>
              </div>
              <label className="field">
                <span className="field-label">Setpoint</span>
                <Qty
                  quantity={pvq}
                  label="Setpoint"
                  placeholder="steady state"
                  value={c.setpoint}
                  onChange={(setpoint) => set(i, { ...c, setpoint: setpoint === "" ? undefined : setpoint })}
                />
              </label>
              <div className="dyn-row">
                <label>
                  <span className="field-label">Action</span>
                  <select value={c.action ?? "reverse"} onChange={(e) => set(i, { ...c, action: e.target.value as Controller["action"] })}>
                    <option value="reverse">reverse (heater on a temperature)</option>
                    <option value="direct">direct (outlet valve on a level)</option>
                  </select>
                </label>
                <label>
                  <span className="field-label">Gain [%/%]</span>
                  <Qty quantity={null} label="Gain" value={c.gain ?? 1} onChange={(v) => set(i, { ...c, gain: Number(v) || 0 })} />
                </label>
              </div>
              <div className="dyn-row">
                <label>
                  <span className="field-label">Integral time</span>
                  <Qty
                    quantity="time"
                    label="Integral time"
                    placeholder="none"
                    value={c.integral_time}
                    onChange={(v) => set(i, { ...c, integral_time: v === "" ? undefined : v })}
                  />
                </label>
                <label>
                  <span className="field-label">Derivative time</span>
                  <Qty
                    quantity="time"
                    label="Derivative time"
                    placeholder="none"
                    value={c.derivative_time}
                    onChange={(v) => set(i, { ...c, derivative_time: v === "" ? undefined : v })}
                  />
                </label>
              </div>
              <button className="danger" onClick={() => onChange(controllers.filter((_, j) => j !== i))}>
                Delete controller
              </button>
            </Card>
          );
        })}
      </ul>
      <button
        disabled={shared.units.length === 0}
        onClick={() =>
          onChange([
            ...controllers,
            { id: nextId("C", controllers.map((c) => c.id)), measure: firstMeasure(), manipulate: firstTarget(), gain: 1, integral_time: "5 min", op_range: ["", ""] },
          ])
        }
      >
        + Controller
      </button>
    </div>
  );
}

function variableQuantityOf(c: Controller, units: UnitInfo[]): string | null {
  if (c.measure.stream !== undefined) {
    const p = c.measure.property;
    return p === "T" ? "temperature" : p === "P" ? "pressure" : ["molar_flow", "component_flow"].includes(p) ? "molar_flow" :
      p === "mass_flow" ? "mass_flow" : p === "enthalpy_flow" ? "power" : null;
  }
  return units.find((u) => u.id === c.measure.unit)?.schema.result_fields.find((f) => f.name === c.measure.property)?.quantity ?? null;
}

function EventsEditor({
  events,
  controllers,
  dynamic,
  onChange,
  ...shared
}: Shared & { events: DynamicEvent[]; controllers: Controller[]; dynamic: Record<string, Params>; onChange: (e: DynamicEvent[]) => void }) {
  const units = withDynamicResults(shared.units, dynamic);
  const set = (i: number, e: DynamicEvent) => onChange(events.map((x, j) => (j === i ? e : x)));
  const feed = units.find((u) => u.schema.type === "Feed");
  return (
    <div className="studies">
      <ul className="study-list">
        {events.map((ev, i) => {
          const isSp = ev.controller !== undefined;
          const ctrl = controllers.find((c) => c.id === ev.controller);
          const q = isSp ? (ctrl ? variableQuantityOf(ctrl, units) : null)
            : quantityOf({ unit: ev.unit ?? "", param: ev.param ?? "", component: ev.component }, units, dynamic);
          return (
            <Card
              key={i}
              title={`at ${ev.time === "" ? "?" : ev.time}`}
              subtitle={`${isSp ? `${ev.controller} setpoint` : `${ev.unit} ${ev.param}${ev.component ? ` of ${ev.component}` : ""}`} → ${ev.value === "" ? "?" : ev.value}${ev.ramp ? ` over ${ev.ramp}` : ""}`}
            >
              <div className="dyn-row">
                <label>
                  <span className="field-label">Time</span>
                  <Qty quantity="time" label="Event time" value={ev.time} onChange={(time) => set(i, { ...ev, time })} />
                </label>
                <label>
                  <span className="field-label">Ramp over</span>
                  <Qty quantity="time" label="Ramp" placeholder="step" value={ev.ramp} onChange={(v) => set(i, { ...ev, ramp: v === "" ? undefined : v })} />
                </label>
              </div>
              <div className="field">
                <span className="field-label">Change</span>
                <select
                  aria-label="What changes"
                  value={isSp ? "controller" : "unit"}
                  onChange={(e) =>
                    set(
                      i,
                      e.target.value === "controller"
                        ? { time: ev.time, controller: controllers[0]?.id ?? "", value: "" }
                        : { time: ev.time, unit: feed?.id ?? units[0]?.id, param: "flows", component: shared.components[0], value: "" },
                    )
                  }
                >
                  <option value="unit">a unit parameter</option>
                  <option value="controller" disabled={!controllers.length}>
                    a controller setpoint
                  </option>
                </select>
                {isSp ? (
                  <select aria-label="Controller" value={ev.controller} onChange={(e) => set(i, { ...ev, controller: e.target.value })}>
                    {controllers.map((c) => (
                      <option key={c.id} value={c.id}>
                        {c.id}
                      </option>
                    ))}
                  </select>
                ) : (
                  <TargetPicker
                    value={{ unit: ev.unit ?? "", param: ev.param ?? "", component: ev.component }}
                    units={units}
                    components={shared.components}
                    dynamic={dynamic}
                    onChange={(m) => set(i, { ...ev, ...m, value: "" })}
                  />
                )}
              </div>
              <label className="field">
                <span className="field-label">New value</span>
                <Qty quantity={q} label="New value" value={ev.value} onChange={(value) => set(i, { ...ev, value })} />
              </label>
              <button className="danger" onClick={() => onChange(events.filter((_, j) => j !== i))}>
                Delete event
              </button>
            </Card>
          );
        })}
      </ul>
      <button
        disabled={!feed}
        onClick={() => onChange([...events, { time: "10 min", unit: feed!.id, param: "flows", component: shared.components[0], value: "" }])}
      >
        + Event
      </button>
    </div>
  );
}

export function DynamicsSettings({
  dynamics,
  onChange,
  ...shared
}: Shared & { dynamics: Dynamics | undefined; onChange: (d: Dynamics | undefined) => void }) {
  if (!dynamics)
    return (
      <div className="studies">
        <p className="muted small">
          Simulate how the plant moves over time: vessels and reactors hold up material, controllers act, feeds and
          setpoints change.
        </p>
        <button onClick={() => onChange({ t_end: "1 h", points: 121, units: {}, controllers: [], events: [] })}>
          + Dynamic simulation
        </button>
      </div>
    );
  const dynamic = dynamics.units;
  const controllers = dynamics.controllers ?? [];
  return (
    <div className="studies">
      <p className="muted small">Starts from the steady state; run it from the Dynamics tab of the results.</p>
      <div className="dyn-row">
        <label>
          <span className="field-label">Run for</span>
          <Qty quantity="time" label="End time" value={dynamics.t_end} onChange={(t_end) => onChange({ ...dynamics, t_end })} />
        </label>
        <label>
          <span className="field-label">Output points</span>
          <input
            type="number"
            min={2}
            max={2001}
            value={dynamics.points ?? 201}
            onChange={(e) => onChange({ ...dynamics, points: Math.max(2, Math.min(2001, Math.round(Number(e.target.value) || 2))) })}
          />
        </label>
      </div>
      <h4>Hold-ups</h4>
      <HoldUps units={shared.units} dynamic={dynamic} onChange={(units) => onChange({ ...dynamics, units })} />
      <h4>Controllers</h4>
      <ControllersEditor {...shared} controllers={controllers} dynamic={dynamic} onChange={(c) => onChange({ ...dynamics, controllers: c })} />
      <h4>Events</h4>
      <EventsEditor
        {...shared}
        events={dynamics.events ?? []}
        controllers={controllers}
        dynamic={dynamic}
        onChange={(events) => onChange({ ...dynamics, events })}
      />
      <button className="danger" onClick={() => onChange(undefined)}>
        Remove dynamic simulation
      </button>
    </div>
  );
}
