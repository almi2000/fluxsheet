import { useEffect, useRef, useState } from "react";
import { api } from "./api";
import { UnitIcon } from "./icons";
import { ProfileView } from "./ProfileView";
import { ReactionsEditor } from "./ReactionsEditor";
import { UNITS, fmt, formatQuantity, joinValue, splitValue } from "./quantities";
import type { Compound, Method, ParamSchema, Params, Reaction, StreamResult, UnitResult, UnitType } from "./types";

// ---- generic inputs ----------------------------------------------------------------------

function QuantityInput({
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

function Field({ param, children }: { param: ParamSchema; children: React.ReactNode }) {
  return (
    <label className="field">
      <span className="field-label" title={param.help}>
        {param.label}
      </span>
      {children}
      {param.help && <span className="field-help">{param.help}</span>}
    </label>
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

// ---- unit parameters ---------------------------------------------------------------------

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
  result,
  onRename,
  onChange,
  onDelete,
}: {
  id: string;
  schema: UnitType;
  params: Params;
  components: string[];
  result?: UnitResult;
  onRename: (newId: string) => string | null;
  onChange: (p: Params) => void;
  onDelete: () => void;
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
        <UnitIcon type={schema.type} size={34} />
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
                      next[n] = p.default ?? (p.quantity && UNITS[p.quantity] ? `0 ${UNITS[p.quantity][0].value}` : 0);
                      onChange(next);
                    }}
                  >
                    {byName[n].label}
                  </button>
                ))}
              </div>
              {active && renderParam(byName[active])}
            </div>
          );
        })}

      {schema.params.filter((p) => !grouped.has(p.name)).map(renderParam)}

      {result && schema.result_fields.some((f) => result[f.name] !== undefined) && (
        <section className="unit-results">
          <h3>Results</h3>
          {schema.result_fields
            .filter((f) => result[f.name] !== undefined)
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
  method,
  methods,
  reactions,
  onChange,
}: {
  name: string;
  description: string;
  components: string[];
  method: string;
  methods: Method[];
  reactions: Reaction[];
  onChange: (patch: {
    name?: string;
    description?: string;
    components?: string[];
    property_method?: string;
    reactions?: Reaction[];
  }) => void;
}) {
  const [suggestion, setSuggestion] = useState<{ method: string; reason: string } | null>(null);
  useEffect(() => {
    if (components.length === 0) return setSuggestion(null);
    let live = true;
    api
      .suggestMethod(components)
      .then((s) => live && setSuggestion(s))
      .catch(() => live && setSuggestion(null));
    return () => {
      live = false;
    };
  }, [components]);

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
        {components.map((c) => (
          <li key={c}>
            {c}
            <button aria-label={`Remove ${c}`} onClick={() => onChange({ components: components.filter((x) => x !== c) })}>
              ×
            </button>
          </li>
        ))}
      </ul>

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
    </div>
  );
}
