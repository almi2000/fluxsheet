import { useEffect, useRef, useState } from "react";
import { api } from "./api";
import { CustomCompoundForm } from "./CustomCompoundForm";
import { UnitIcon } from "./icons";
import { ProfileView } from "./ProfileView";
import { ReactionsEditor } from "./ReactionsEditor";
import { UNITS, fmt, formatQuantity, joinValue, splitValue } from "./quantities";
import type { Compound, CustomCompound, Method, ParamSchema, Params, Reaction, StreamResult, UnitResult, UnitType } from "./types";

// ---- generic inputs ----------------------------------------------------------------------

export function QuantityInput({
  quantity,
  value,
  placeholder,
  onChange,
}: {
  quantity: string | null;
  value: unknown;
  placeholder?: string;
  onChange: (v: unknown) => void;
}) {
  const units = quantity ? UNITS[quantity] : undefined;
  const { num, unit } = splitValue(value, quantity ?? "");
  return (
    <div className="quantity">
      <input
        type="text"
        inputMode="decimal"
        value={num}
        placeholder={placeholder}
        onChange={(e) => onChange(units ? joinValue(e.target.value, unit) : numOrNull(e.target.value))}
      />
      {units && (
        <select value={unit} onChange={(e) => onChange(joinValue(num, e.target.value))}>
          {units.map((u) => (
            <option key={u.value} value={u.value}>
              {u.label}
            </option>
          ))}
          {!units.some((u) => u.value === unit) && <option value={unit}>{unit}</option>}
        </select>
      )}
    </div>
  );
}

function numOrNull(text: string): number | null {
  return text.trim() === "" ? null : Number(text);
}

function Field({ param, children, group }: { param: ParamSchema; children: React.ReactNode; group?: boolean }) {
  // A group of several inputs (a table) is not one labelled control: a <label> would forward
  // clicks on its caption to the first button or input inside.
  const Tag = group ? "div" : "label";
  return (
    <Tag className="field">
      <span className="field-label" title={param.help}>
        {param.label}
      </span>
      {children}
      {param.help && <span className="field-help">{param.help}</span>}
    </Tag>
  );
}

function ComponentTable({
  components,
  value,
  quantity,
  onChange,
}: {
  components: string[];
  value: Record<string, unknown>;
  quantity: "molar_flow" | "fraction";
  onChange: (v: Record<string, unknown>) => void;
}) {
  if (components.length === 0)
    return <p className="muted small">Add components in the flowsheet settings first.</p>;
  return (
    <div className="component-table">
      {components.map((c) => (
        <div key={c} className="component-row">
          <span>{c}</span>
          <QuantityInput
            quantity={quantity === "molar_flow" ? "molar_flow" : null}
            value={value[c]}
            placeholder="0"
            onChange={(v) => {
              const next = { ...value };
              if (v === null || v === undefined) delete next[c];
              else next[c] = v;
              onChange(next);
            }}
          />
        </div>
      ))}
    </div>
  );
}

/** One stage number per feed stream; empty fields fall back to the common feed stage. */
function StreamStageTable({
  streams,
  value,
  fallback,
  onChange,
}: {
  streams: string[];
  value: Record<string, number>;
  fallback: unknown;
  onChange: (v: Record<string, number> | null) => void;
}) {
  if (streams.length < 2)
    return <p className="muted small">Connect more feed streams to give each its own stage.</p>;
  return (
    <div className="component-table">
      {streams.map((sid) => (
        <div key={sid} className="component-row">
          <span>{sid}</span>
          <input
            type="number"
            step={1}
            min={1}
            value={value[sid] === undefined ? "" : String(value[sid])}
            placeholder={fallback === undefined || fallback === null ? "" : String(fallback)}
            onChange={(e) => {
              // Keep only the streams still connected, so deleted streams leave nothing behind.
              const next = Object.fromEntries(Object.entries(value).filter(([k]) => streams.includes(k)));
              if (e.target.value === "") delete next[sid];
              else next[sid] = Math.round(Number(e.target.value));
              onChange(Object.keys(next).length ? next : null);
            }}
          />
        </div>
      ))}
    </div>
  );
}

/** A "records" parameter: one row per entry, one column per field; rows can be added and removed. */
function RecordTable({
  fields,
  value,
  onChange,
}: {
  fields: ParamSchema[];
  value: Record<string, unknown>[];
  onChange: (v: Record<string, unknown>[] | null) => void;
}) {
  const update = (i: number, name: string, v: unknown) => {
    const next = value.map((row) => ({ ...row }));
    if (v === null || v === undefined || v === "") delete next[i][name];
    else next[i][name] = v;
    onChange(next);
  };
  const cell = (f: ParamSchema, row: Record<string, unknown>, i: number) => {
    const current = row[f.name] ?? undefined;
    if (f.quantity === "choice")
      return (
        <select value={String(current ?? f.default ?? "")} onChange={(e) => update(i, f.name, e.target.value || null)}>
          {f.default === null && <option value="">–</option>}
          {(f.options ?? []).map((o) => (
            <option key={o} value={o}>
              {o}
            </option>
          ))}
        </select>
      );
    if (f.quantity === "integer")
      return (
        <input
          type="number"
          step={1}
          min={f.min ?? undefined}
          value={current === undefined ? "" : String(current)}
          onChange={(e) => update(i, f.name, e.target.value === "" ? null : Math.round(Number(e.target.value)))}
        />
      );
    return (
      <QuantityInput
        quantity={f.quantity === "fraction" ? null : f.quantity}
        value={current}
        placeholder={f.default === null ? "–" : String(f.default)}
        onChange={(v) => update(i, f.name, v)}
      />
    );
  };
  return (
    <div className="record-table">
      {value.map((row, i) => (
        <div key={i} className="record-row">
          {fields.map((f) => (
            <label key={f.name} className="record-cell" title={f.help}>
              <span>{f.label}</span>
              {cell(f, row, i)}
            </label>
          ))}
          <button
            className="record-remove"
            title="Remove this entry"
            onClick={(e) => {
              e.preventDefault();
              const next = value.filter((_, k) => k !== i);
              onChange(next.length ? next : null);
            }}
          >
            ×
          </button>
        </div>
      ))}
      <button
        className="record-add"
        onClick={(e) => {
          e.preventDefault();
          const blank = Object.fromEntries(fields.filter((f) => f.default !== null).map((f) => [f.name, f.default]));
          onChange([...value, blank]);
        }}
      >
        + Add
      </button>
    </div>
  );
}

export type CurveValue = {
  flow_unit?: string;
  head_unit?: string;
  efficiency_unit?: string;
  points: (number | null)[][];
};

const HEAD_UNITS = [...UNITS.length, ...UNITS.specific_energy];

/** Performance curve: rows of flow, head and (optionally) efficiency, with their units. */
export function CurveInput({ value, onChange }: { value: CurveValue; onChange: (v: CurveValue) => void }) {
  const points = value.points ?? [];
  const width = points.length ? points[0].length : value.efficiency_unit === undefined ? 2 : 3;
  const effUnit = width === 3 ? (value.efficiency_unit ?? "") : null;
  const update = (patch: Partial<CurveValue>) => onChange({ ...value, ...patch });
  const setCell = (row: number, col: number, text: string) =>
    update({ points: points.map((r, i) => (i === row ? r.map((v, j) => (j === col ? numOrNull(text) : v)) : r)) });
  const setEfficiency = (unit: string) => {
    if (unit === "none") {
      const { efficiency_unit: _, ...rest } = value;
      onChange({ ...rest, points: points.map((r) => r.slice(0, 2)) });
    } else update({ efficiency_unit: unit, points: points.map((r) => (r.length === 3 ? r : [...r, null])) });
  };
  const select = (current: string, options: { value: string; label: string }[], set: (v: string) => void) => (
    <select value={current} onChange={(e) => set(e.target.value)}>
      {options.map((u) => (
        <option key={u.value} value={u.value}>
          {u.label}
        </option>
      ))}
      {!options.some((u) => u.value === current) && <option value={current}>{current}</option>}
    </select>
  );
  return (
    <div className="curve-table">
      <div className="curve-row curve-head">
        <span>Flow</span>
        <span>Head</span>
        <span>Efficiency</span>
        <span />
      </div>
      <div className="curve-row">
        {select(value.flow_unit ?? "m3/h", UNITS.volumetric_flow, (u) => update({ flow_unit: u }))}
        {select(value.head_unit ?? "m", HEAD_UNITS, (u) => update({ head_unit: u }))}
        {select(
          effUnit === null ? "none" : effUnit,
          [
            { value: "none", label: "none" },
            { value: "", label: "fraction" },
            { value: "%", label: "%" },
          ],
          setEfficiency,
        )}
        <span />
      </div>
      {points.map((row, i) => (
        <div key={i} className="curve-row">
          {[0, 1, 2].map((j) =>
            j < row.length ? (
              <input
                key={j}
                type="number"
                step="any"
                aria-label={["flow", "head", "efficiency"][j] + ` of point ${i + 1}`}
                value={row[j] === null || row[j] === undefined ? "" : String(row[j])}
                onChange={(e) => setCell(i, j, e.target.value)}
              />
            ) : (
              <span key={j} />
            ),
          )}
          <button
            type="button"
            className="icon-button"
            title="Remove point"
            onClick={() => update({ points: points.filter((_, k) => k !== i) })}
          >
            ×
          </button>
        </div>
      ))}
      <button
        type="button"
        className="small"
        onClick={() => {
          const last = points[points.length - 1];
          update({ points: [...points, last ? last.map(() => null) : width === 3 ? [null, null, null] : [null, null]] });
        }}
      >
        + Add point
      </button>
    </div>
  );
}

// ---- unit parameters ---------------------------------------------------------------------

/** An empty curve table in the unit's usual head unit, to be filled point by point. */
function newCurve(p: ParamSchema): CurveValue {
  return { flow_unit: "m3/h", head_unit: p.options?.[0] ?? "m", points: [] };
}

function noneLabel(group: string[]): string {
  if (group.includes("dP") && group.includes("P")) return "No ΔP";
  if (group.includes("duty")) return "Adiabatic";
  return "Default";
}

export function UnitInspector({
  id,
  schema,
  params,
  components,
  utilities = [],
  inletStreams = {},
  result,
  onRename,
  onChange,
  onDelete,
  children,
}: {
  id: string;
  schema: UnitType;
  params: Params;
  components: string[];
  utilities?: string[]; // ids of the flowsheet's utilities
  inletStreams?: Record<string, string[]>; // connected stream ids per inlet port
  result?: UnitResult;
  onRename: (newId: string) => string | null;
  onChange: (p: Params) => void;
  onDelete: () => void;
  children?: React.ReactNode; // extra sections after the parameters, e.g. the property method
}) {
  const [draftId, setDraftId] = useState(id);
  const [idError, setIdError] = useState<string | null>(null);
  useEffect(() => {
    setDraftId(id);
    setIdError(null);
  }, [id]);

  const set = (name: string, v: unknown) => {
    const next = { ...params };
    if (v === null || v === undefined) delete next[name];
    else next[name] = v;
    onChange(next);
  };

  const grouped = new Set(schema.exclusive_groups.filter((g) => g.length > 1).flat());
  const byName = Object.fromEntries(schema.params.map((p) => [p.name, p]));

  const renderParam = (p: ParamSchema) => {
    const current = params[p.name] ?? (p.default !== null ? p.default : undefined);
    if (p.quantity === "integer")
      return (
        <Field key={p.name} param={p}>
          <input
            type="number"
            step={1}
            min={p.min ?? undefined}
            max={p.max ?? undefined}
            value={current === undefined ? "" : String(current)}
            placeholder={p.default === null ? "default" : String(p.default)}
            onChange={(e) => set(p.name, e.target.value === "" ? null : Math.round(Number(e.target.value)))}
          />
        </Field>
      );
    if (p.quantity === "choice")
      return (
        <Field key={p.name} param={p}>
          <select value={String(current ?? "")} onChange={(e) => set(p.name, e.target.value || null)}>
            {p.default === null && <option value="">–</option>}
            {(p.options ?? []).map((o) => (
              <option key={o} value={o}>
                {o}
              </option>
            ))}
          </select>
        </Field>
      );
    if (p.quantity === "component")
      return (
        <Field key={p.name} param={p}>
          <select value={String(current ?? "")} onChange={(e) => set(p.name, e.target.value || null)}>
            <option value="">Choose a component…</option>
            {components.map((c) => (
              <option key={c} value={c}>
                {c}
              </option>
            ))}
            {typeof current === "string" && current && !components.includes(current) && (
              <option value={current}>{current} (not in component list)</option>
            )}
          </select>
        </Field>
      );
    if (p.quantity === "utility")
      return (
        <Field key={p.name} param={p}>
          <select value={String(current ?? "")} onChange={(e) => set(p.name, e.target.value || null)}>
            <option value="">{utilities.length ? "None" : "None (add utilities in the flowsheet settings)"}</option>
            {utilities.map((u) => (
              <option key={u} value={u}>
                {u}
              </option>
            ))}
            {typeof current === "string" && current && !utilities.includes(current) && (
              <option value={current}>{current} (not defined)</option>
            )}
          </select>
        </Field>
      );
    if (p.quantity === "stream") {
      const streams = Object.values(inletStreams).flat();
      return (
        <Field key={p.name} param={p}>
          <select value={String(current ?? "")} onChange={(e) => set(p.name, e.target.value || null)}>
            <option value="">{streams.length ? "Choose a stream…" : "Connect inlet streams first"}</option>
            {streams.map((sid) => (
              <option key={sid} value={sid}>
                {sid}
              </option>
            ))}
            {typeof current === "string" && current && !streams.includes(current) && (
              <option value={current}>{current} (not connected)</option>
            )}
          </select>
        </Field>
      );
    }
    if (p.quantity === "curve")
      return (
        <div key={p.name} className="field" role="group" aria-label={p.label}>
          <span className="field-label" title={p.help}>
            {p.label}
          </span>
          <CurveInput value={(params[p.name] as CurveValue) ?? newCurve(p)} onChange={(v) => set(p.name, v)} />
          {p.help && <span className="field-help">{p.help}</span>}
        </div>
      );
    if (p.quantity === "records")
      return (
        <Field key={p.name} param={p} group>
          <RecordTable
            fields={p.fields ?? []}
            value={Array.isArray(params[p.name]) ? (params[p.name] as Record<string, unknown>[]) : []}
            onChange={(v) => set(p.name, v)}
          />
        </Field>
      );
    if (p.quantity === "string_list")
      return (
        <Field key={p.name} param={p}>
          <input
            type="text"
            placeholder="comma-separated, e.g. R1, R2"
            value={Array.isArray(current) ? current.join(", ") : ""}
            onChange={(e) => {
              const items = e.target.value.split(",").map((v) => v.trim()).filter(Boolean);
              set(p.name, items.length ? items : null);
            }}
          />
        </Field>
      );
    if (p.quantity === "stream_stages")
      return (
        <Field key={p.name} param={p}>
          <StreamStageTable
            streams={Object.values(inletStreams).flat()}
            value={(params[p.name] as Record<string, number>) ?? {}}
            fallback={params.feed_stage ?? byName.feed_stage?.default}
            onChange={(v) => set(p.name, v)}
          />
        </Field>
      );
    if (p.quantity === "composition" || p.quantity === "component_fractions")
      return (
        <Field key={p.name} param={p}>
          <ComponentTable
            components={components}
            value={(params[p.name] as Record<string, unknown>) ?? {}}
            quantity={p.quantity === "composition" ? "molar_flow" : "fraction"}
            onChange={(v) => set(p.name, v)}
          />
        </Field>
      );
    return (
      <Field key={p.name} param={p}>
        <QuantityInput
          quantity={p.quantity === "fraction" ? null : p.quantity}
          value={current}
          placeholder={p.default === null ? "default" : String(p.default)}
          onChange={(v) => set(p.name, v)}
        />
      </Field>
    );
  };

  return (
    <div className="inspector-body">
      <header className="inspector-head">
        <UnitIcon type={schema.type} size={38} />
        <div>
          <div className="eyebrow">{schema.label}</div>
          <input
            className="id-input"
            value={draftId}
            onChange={(e) => setDraftId(e.target.value)}
            onBlur={() => setIdError(draftId === id ? null : onRename(draftId))}
            onKeyDown={(e) => e.key === "Enter" && (e.target as HTMLInputElement).blur()}
          />
          {idError && <div className="field-error">{idError}</div>}
        </div>
      </header>
      <p className="muted small">{schema.description}</p>

      {schema.exclusive_groups
        .filter((g) => g.length > 1)
        .map((group) => {
          const required = schema.required_groups.some((r) => r.join() === group.join());
          const active = group.find((n) => params[n] !== undefined && params[n] !== null) ?? null;
          return (
            <div key={group.join()} className="spec-group">
              <div className="segmented">
                {!required && (
                  <button
                    className={active === null ? "on" : ""}
                    onClick={() => {
                      const next = { ...params };
                      group.forEach((n) => delete next[n]);
                      onChange(next);
                    }}
                  >
                    {noneLabel(group)}
                  </button>
                )}
                {group.map((n) => (
                  <button
                    key={n}
                    className={active === n ? "on" : ""}
                    onClick={() => {
                      if (active === n) return;
                      const next = { ...params };
                      group.forEach((g) => delete next[g]);
                      const p = byName[n];
                      next[n] =
                        p.quantity === "curve"
                          ? newCurve(p)
                          : (p.default ?? (p.quantity && UNITS[p.quantity] ? `0 ${UNITS[p.quantity][0].value}` : 0));
                      onChange(next);
                    }}
                  >
                    {byName[n].tab ?? byName[n].label}
                  </button>
                ))}
              </div>
              {active && renderParam(byName[active])}
            </div>
          );
        })}

      {schema.params
        .filter((p) => !grouped.has(p.name) && (!p.requires || (params[p.requires] ?? null) !== null))
        .map(renderParam)}
      {children}

      {result && schema.result_fields.some((f) => result[f.name] != null) && (
        <section className="unit-results">
          <h3>Results</h3>
          {schema.result_fields
            .filter((f) => result[f.name] != null)
            .map((f) => (
              <div key={f.name} className="kv">
                <span>{f.label}</span>
                <b>{formatQuantity(result[f.name], f.quantity)}</b>
              </div>
            ))}
        </section>
      )}
      {result?.profiles &&
        Object.entries(result.profiles).map(([key, profile]) => <ProfileView key={key} profile={profile} />)}

      <button className="danger" onClick={onDelete}>
        Delete unit
      </button>
    </div>
  );
}

// ---- stream results ----------------------------------------------------------------------

export function StreamInspector({ id, result }: { id: string; result?: StreamResult }) {
  return (
    <div className="inspector-body">
      <header className="inspector-head">
        <div>
          <div className="eyebrow">Stream</div>
          <h2>{id}</h2>
        </div>
      </header>
      {!result ? (
        <p className="muted">Run the simulation to see this stream's conditions.</p>
      ) : (
        <>
          <div className="kv"><span>Temperature</span><b>{fmt.T(result.T)} °C</b></div>
          <div className="kv"><span>Pressure</span><b>{fmt.P(result.P)} bar</b></div>
          <div className="kv"><span>Phase</span><b>{result.phase}</b></div>
          <div className="kv"><span>Vapor fraction</span><b>{fmt.frac(result.vapor_fraction)}</b></div>
          <div className="kv"><span>Molar flow</span><b>{fmt.molar(result.molar_flow)} kmol/h</b></div>
          <div className="kv"><span>Mass flow</span><b>{fmt.mass(result.mass_flow)} kg/h</b></div>
          <h3>Composition</h3>
          <table className="mini">
            <thead>
              <tr><th>Component</th><th>kmol/h</th><th>mol frac.</th></tr>
            </thead>
            <tbody>
              {Object.keys(result.component_flows).map((c) => (
                <tr key={c}>
                  <td>{c}</td>
                  <td>{fmt.molar(result.component_flows[c])}</td>
                  <td>{result.molar_flow > 0 ? fmt.frac(result.mole_fractions[c]) : "–"}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </>
      )}
    </div>
  );
}

// ---- flowsheet settings ------------------------------------------------------------------

function CompoundSearch({ onPick }: { onPick: (c: Compound) => void }) {
  const [q, setQ] = useState("");
  const [hits, setHits] = useState<Compound[]>([]);
  const [busy, setBusy] = useState(false);
  const seq = useRef(0);

  useEffect(() => {
    if (q.trim().length < 2) {
      setHits([]);
      return;
    }
    const n = ++seq.current;
    setBusy(true);
    const t = setTimeout(() => {
      api
        .searchCompounds(q)
        .then((r) => n === seq.current && setHits(r))
        .catch(() => n === seq.current && setHits([]))
        .finally(() => n === seq.current && setBusy(false));
    }, 180);
    return () => clearTimeout(t);
  }, [q]);

  return (
    <div className="compound-search">
      <input
        type="search"
        placeholder="Add compound: name, CAS or formula…"
        value={q}
        onChange={(e) => setQ(e.target.value)}
      />
      {(hits.length > 0 || busy) && (
        <ul className="dropdown">
          {busy && hits.length === 0 && <li className="muted">Searching…</li>}
          {hits.map((c) => (
            <li key={c.cas}>
              <button
                disabled={!c.simulation_ready}
                title={c.simulation_ready ? "" : "Missing critical data; cannot be simulated yet"}
                onClick={() => {
                  onPick(c);
                  setQ("");
                  setHits([]);
                }}
              >
                <span>{c.name}</span>
                <span className="muted small">
                  {c.formula} · {c.cas}
                </span>
              </button>
            </li>
          ))}
        </ul>
      )}
    </div>
  );
}

export function FlowsheetInspector({
  name,
  description,
  components,
  customCompounds = [],
  method,
  methods,
  reactions,
  onChange,
  children,
}: {
  name: string;
  description: string;
  components: string[];
  customCompounds?: CustomCompound[];
  method: string;
  methods: Method[];
  reactions: Reaction[];
  onChange: (patch: {
    name?: string;
    description?: string;
    components?: string[];
    custom_compounds?: CustomCompound[];
    property_method?: string;
    reactions?: Reaction[];
  }) => void;
  children?: React.ReactNode; // further settings sections, e.g. design specs
}) {
  const [suggestion, setSuggestion] = useState<{ method: string; reason: string } | null>(null);
  // The custom compound form: closed, a new compound ("new") or the one being edited.
  const [customForm, setCustomForm] = useState<CustomCompound | "new" | null>(null);
  const customOf = (c: string) => customCompounds.find((d) => d.name === c);
  useEffect(() => {
    if (components.length === 0) return setSuggestion(null);
    let live = true;
    const used = customCompounds.filter((d) => components.includes(d.name));
    api
      .suggestMethod(components, used.length ? used : undefined)
      .then((s) => live && setSuggestion(s))
      .catch(() => live && setSuggestion(null));
    return () => {
      live = false;
    };
  }, [components, customCompounds]);

  const label = (key: string) => methods.find((m) => m.key === key)?.label ?? key;

  return (
    <div className="inspector-body">
      <header className="inspector-head">
        <div>
          <div className="eyebrow">Flowsheet</div>
          <input className="id-input wide" value={name} onChange={(e) => onChange({ name: e.target.value })} />
        </div>
      </header>
      <textarea
        className="description"
        placeholder="What does this process do?"
        value={description}
        onChange={(e) => onChange({ description: e.target.value })}
      />

      <h3>Components</h3>
      <CompoundSearch
        onPick={(c) => !components.includes(c.name) && onChange({ components: [...components, c.name] })}
      />
      <ul className="chips">
        {components.map((c) => {
          const custom = customOf(c);
          return (
            <li key={c} className={custom ? "custom" : undefined}>
              {custom ? (
                <button className="chip-edit" title={`User-defined (${custom.formula}); click to edit`} onClick={() => setCustomForm(custom)}>
                  {c} <span className="chip-tag">custom</span>
                </button>
              ) : (
                c
              )}
              <button
                aria-label={`Remove ${c}`}
                onClick={() =>
                  onChange({
                    components: components.filter((x) => x !== c),
                    ...(custom ? { custom_compounds: customCompounds.filter((d) => d !== custom) } : {}),
                  })
                }
              >
                ×
              </button>
            </li>
          );
        })}
      </ul>
      {customForm ? (
        <CustomCompoundForm
          key={customForm === "new" ? "new" : customForm.name}
          components={components}
          editing={customForm === "new" ? null : customForm}
          onClose={() => setCustomForm(null)}
          onSave={(def, previous) => {
            onChange({
              custom_compounds: [...customCompounds.filter((d) => d.name !== previous), def],
              components: previous
                ? components.map((c) => (c === previous ? def.name : c))
                : [...components, def.name],
            });
            setCustomForm(null);
          }}
        />
      ) : (
        <button className="cc-add" onClick={() => setCustomForm("new")}>
          + Custom compound…
        </button>
      )}

      <h3>Property method</h3>
      <select value={method} onChange={(e) => onChange({ property_method: e.target.value })}>
        <option value="auto">Automatic{suggestion ? ` (${label(suggestion.method)})` : ""}</option>
        {methods.map((m) => (
          <option key={m.key} value={m.key}>
            {m.label}
          </option>
        ))}
      </select>
      {method === "auto" && suggestion && <p className="callout">{suggestion.reason}</p>}
      {method !== "auto" && (
        <p className="muted small">{methods.find((m) => m.key === method)?.description}</p>
      )}
      {method !== "auto" && suggestion && suggestion.method !== method && (
        <p className="callout warn">
          For these components, {label(suggestion.method)} is usually the better choice. {suggestion.reason}
        </p>
      )}

      <h3>Reactions</h3>
      <ReactionsEditor reactions={reactions} components={components} onChange={(r) => onChange({ reactions: r })} />
      {children}
    </div>
  );
}
