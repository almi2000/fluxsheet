// Flowsheet-level reactions: stoichiometry, optional power-law kinetics and an optional
// equilibrium-constant correlation. Reactors refer to these by id.
import { useEffect, useState } from "react";
import { UNITS, joinValue, splitValue } from "./quantities";
import type { Keq, Kinetics, RateLaw, Reaction, Value } from "./types";

// ---- small inputs ------------------------------------------------------------------------

/** Number input that keeps partial text such as "-" or "1e" while typing. */
function NumInput({
  value,
  placeholder,
  onChange,
  label,
}: {
  value: number | undefined | null;
  placeholder?: string;
  onChange: (v: number | null) => void;
  label?: string;
}) {
  const [text, setText] = useState(value === undefined || value === null ? "" : String(value));
  useEffect(() => {
    const parsed = text.trim() === "" ? null : Number(text);
    if ((value ?? null) !== parsed) setText(value === undefined || value === null ? "" : String(value));
  }, [value]); // eslint-disable-line react-hooks/exhaustive-deps
  return (
    <input
      type="text"
      inputMode="decimal"
      aria-label={label}
      value={text}
      placeholder={placeholder}
      onChange={(e) => {
        setText(e.target.value);
        const t = e.target.value.trim();
        if (t === "") onChange(null);
        else if (Number.isFinite(Number(t))) onChange(Number(t));
      }}
    />
  );
}

function QtyInput({
  quantity,
  value,
  placeholder,
  onChange,
}: {
  quantity: string;
  value: Value | undefined | null;
  placeholder?: string;
  onChange: (v: Value | null) => void;
}) {
  const { num, unit } = splitValue(value, quantity);
  return (
    <div className="quantity">
      <input
        type="text"
        inputMode="decimal"
        value={num}
        placeholder={placeholder}
        onChange={(e) => onChange(joinValue(e.target.value, unit))}
      />
      <select value={unit} onChange={(e) => onChange(joinValue(num, e.target.value))}>
        {UNITS[quantity].map((u) => (
          <option key={u.value} value={u.value}>
            {u.label}
          </option>
        ))}
        {!UNITS[quantity].some((u) => u.value === unit) && <option value={unit}>{unit}</option>}
      </select>
    </div>
  );
}

function Row({ label, help, children }: { label: string; help?: string; children: React.ReactNode }) {
  return (
    <label className="field">
      <span className="field-label">{label}</span>
      {children}
      {help && <span className="field-help">{help}</span>}
    </label>
  );
}

// ---- helpers -----------------------------------------------------------------------------

export function equation(r: Reaction): string {
  const side = (sign: number) =>
    Object.entries(r.stoichiometry)
      .filter(([, v]) => v * sign > 0)
      .map(([c, v]) => (Math.abs(v) === 1 ? c : `${Math.abs(v)} ${c}`))
      .join(" + ") || "…";
  const arrow = r.kinetics?.reverse || r.kinetics?.equilibrium_limited || !r.kinetics ? "⇌" : "→";
  return `${side(-1)} ${arrow} ${side(1)}`;
}

function nextReactionId(reactions: Reaction[]): string {
  let i = 1;
  while (reactions.some((r) => r.id === `R${i}`)) i++;
  return `R${i}`;
}

/** Orders per component, shown with the elementary default as placeholder. */
function OrdersTable({
  components,
  orders,
  defaults,
  onChange,
}: {
  components: string[];
  orders: Record<string, number> | null | undefined;
  defaults: Record<string, number>;
  onChange: (o: Record<string, number> | null) => void;
}) {
  const involved = components.filter((c) => orders?.[c] !== undefined || defaults[c] !== undefined);
  const rest = components.filter((c) => !involved.includes(c));
  const [extra, setExtra] = useState<string[]>([]);
  const shown = [...involved, ...extra.filter((c) => rest.includes(c))];
  const set = (c: string, v: number | null) => {
    // Editing any order switches from the defaults to an explicit list.
    const next: Record<string, number> = orders ? { ...orders } : { ...defaults };
    if (v === null) delete next[c];
    else next[c] = v;
    onChange(next);
  };
  return (
    <div className="component-table">
      {shown.map((c) => (
        <div key={c} className="component-row">
          <span>{c}</span>
          <NumInput
            label={`order of ${c}`}
            value={orders ? orders[c] : undefined}
            placeholder={orders ? "0" : String(defaults[c] ?? 0)}
            onChange={(v) => set(c, v)}
          />
        </div>
      ))}
      <div className="reaction-actions">
        {rest.filter((c) => !extra.includes(c)).length > 0 && (
          <select value="" onChange={(e) => e.target.value && setExtra([...extra, e.target.value])}>
            <option value="">+ order for…</option>
            {rest
              .filter((c) => !extra.includes(c))
              .map((c) => (
                <option key={c} value={c}>
                  {c}
                </option>
              ))}
          </select>
        )}
        {orders && (
          <button className="link" onClick={() => onChange(null)}>
            reset to elementary
          </button>
        )}
      </div>
    </div>
  );
}

function RateLawFields({
  law,
  components,
  defaults,
  basis,
  onChange,
}: {
  law: RateLaw;
  components: string[];
  defaults: Record<string, number>;
  basis: string;
  onChange: (l: RateLaw) => void;
}) {
  const unit = basis === "partial_pressure" ? "Pa" : "mol/m³";
  return (
    <>
      <Row label="Pre-exponential factor k0" help={`SI: r in mol/(m³·s) with ${basis === "partial_pressure" ? "partial pressures" : "concentrations"} in ${unit}.`}>
        <NumInput label="k0" value={law.k0} placeholder="e.g. 1e6" onChange={(v) => onChange({ ...law, k0: v ?? 0 })} />
      </Row>
      <Row label="Activation energy Ea">
        <QtyInput quantity="molar_energy" value={law.Ea} placeholder="0" onChange={(v) => onChange({ ...law, Ea: v ?? 0 })} />
      </Row>
      <Row label="Reference temperature" help="Optional. If set, k0 is the rate constant at this temperature.">
        <QtyInput
          quantity="temperature"
          value={law.T_ref ?? null}
          placeholder="none"
          onChange={(v) => {
            const next = { ...law };
            if (v === null) delete next.T_ref;
            else next.T_ref = v;
            onChange(next);
          }}
        />
      </Row>
      <Row label="Reaction orders" help="Blank uses the stoichiometric coefficients (elementary reaction).">
        <OrdersTable
          components={components}
          orders={law.orders}
          defaults={defaults}
          onChange={(o) => {
            const next = { ...law };
            if (o === null) delete next.orders;
            else next.orders = o;
            onChange(next);
          }}
        />
      </Row>
    </>
  );
}

// ---- one reaction ------------------------------------------------------------------------

function ReactionCard({
  reaction: r,
  components,
  takenIds,
  onChange,
  onDelete,
}: {
  reaction: Reaction;
  components: string[];
  takenIds: string[];
  onChange: (r: Reaction) => void;
  onDelete: () => void;
}) {
  const [open, setOpen] = useState(Object.keys(r.stoichiometry).length === 0);
  const [draftId, setDraftId] = useState(r.id);
  useEffect(() => setDraftId(r.id), [r.id]);
  const idError = draftId.trim() === "" ? "Give it an id" : takenIds.includes(draftId.trim()) ? "Id already used" : null;

  const reactants = Object.fromEntries(Object.entries(r.stoichiometry).filter(([, v]) => v < 0).map(([c, v]) => [c, -v]));
  const products = Object.fromEntries(Object.entries(r.stoichiometry).filter(([, v]) => v > 0));
  const unknown = Object.keys(r.stoichiometry).filter((c) => !components.includes(c));
  const kin = r.kinetics ?? null;
  const reverseMode = !kin ? "none" : kin.reverse ? "rate" : kin.equilibrium_limited ? "equilibrium" : "irreversible";

  const setKinetics = (k: Kinetics | null) => {
    const next = { ...r };
    if (k === null) delete next.kinetics;
    else next.kinetics = k;
    onChange(next);
  };
  const setKeq = (k: Keq | null) => {
    const next = { ...r };
    if (k === null) delete next.keq;
    else next.keq = k;
    onChange(next);
  };

  return (
    <li className="reaction-card">
      <button className="reaction-summary" onClick={() => setOpen(!open)} aria-expanded={open}>
        <b>{r.id}</b>
        <span className="small">{equation(r)}</span>
        <span className="muted small">{kin ? "kinetics" : ""}{kin && r.keq ? " · " : ""}{r.keq ? "K(T)" : ""}</span>
      </button>
      {open && (
        <div className="reaction-body">
          <Row label="Id">
            <input
              value={draftId}
              onChange={(e) => setDraftId(e.target.value)}
              onBlur={() => !idError && draftId.trim() !== r.id && onChange({ ...r, id: draftId.trim() })}
              onKeyDown={(e) => e.key === "Enter" && (e.target as HTMLInputElement).blur()}
            />
            {idError && draftId !== r.id && <span className="field-error">{idError}</span>}
          </Row>

          <Row label="Stoichiometry" help="Negative for reactants, positive for products. Elements must balance.">
            {components.length === 0 ? (
              <p className="muted small">Add components first.</p>
            ) : (
              <div className="component-table">
                {components.map((c) => (
                  <div key={c} className="component-row">
                    <span>{c}</span>
                    <NumInput
                      label={`coefficient of ${c}`}
                      value={r.stoichiometry[c]}
                      placeholder="–"
                      onChange={(v) => {
                        const s = { ...r.stoichiometry };
                        if (v === null || v === 0) delete s[c];
                        else s[c] = v;
                        onChange({ ...r, stoichiometry: s });
                      }}
                    />
                  </div>
                ))}
              </div>
            )}
          </Row>
          {unknown.length > 0 && (
            <p className="callout warn">Not in the component list: {unknown.join(", ")}</p>
          )}

          <Row label="Kinetics" help="Needed by CSTR and plug flow reactors.">
            <select
              value={kin ? "power" : "none"}
              onChange={(e) =>
                setKinetics(e.target.value === "none" ? null : { k0: 1, Ea: "50 kJ/mol", basis: "concentration" })
              }
            >
              <option value="none">None (stoichiometric or equilibrium reactors only)</option>
              <option value="power">Power law  r = k0·exp(−Ea/RT)·Π cᵢ^aᵢ</option>
            </select>
          </Row>
          {kin && (
            <div className="reaction-section">
              <Row label="Basis">
                <select value={kin.basis ?? "concentration"} onChange={(e) => setKinetics({ ...kin, basis: e.target.value })}>
                  <option value="concentration">Concentration (mol/m³)</option>
                  <option value="partial_pressure">Partial pressure (Pa), gas only</option>
                </select>
              </Row>
              <RateLawFields
                law={kin}
                components={components}
                defaults={reactants}
                basis={kin.basis ?? "concentration"}
                onChange={(l) => setKinetics({ ...kin, ...l, orders: l.orders, T_ref: l.T_ref })}
              />
              <Row label="Reverse reaction">
                <select
                  value={reverseMode}
                  onChange={(e) => {
                    const base: Kinetics = { ...kin };
                    delete base.reverse;
                    delete base.equilibrium_limited;
                    if (e.target.value === "equilibrium") setKinetics({ ...base, equilibrium_limited: true });
                    else if (e.target.value === "rate") setKinetics({ ...base, reverse: { k0: 1, Ea: "50 kJ/mol" } });
                    else setKinetics(base);
                  }}
                >
                  <option value="irreversible">Irreversible</option>
                  <option value="equilibrium">Equilibrium-limited  r·(1 − Q/K)</option>
                  <option value="rate">Own reverse rate law</option>
                </select>
              </Row>
              {kin.reverse && (
                <div className="reaction-section">
                  <RateLawFields
                    law={kin.reverse}
                    components={components}
                    defaults={products}
                    basis={kin.basis ?? "concentration"}
                    onChange={(l) => setKinetics({ ...kin, reverse: l })}
                  />
                </div>
              )}
            </div>
          )}

          <Row
            label="Equilibrium constant"
            help={
              r.keq
                ? "ln K = A + B/T + C·ln T + D·T, T in K; K on activities (gas: y·P/1 bar, liquid: x)."
                : "From Gibbs energies of formation (ideal gas). Liquid-phase equilibria need your own K(T)."
            }
          >
            <select value={r.keq ? "user" : "data"} onChange={(e) => setKeq(e.target.value === "user" ? {} : null)}>
              <option value="data">From formation data (ΔG°, ΔH°, Cp)</option>
              <option value="user">Own correlation ln K(T)</option>
            </select>
          </Row>
          {r.keq && (
            <div className="keq-grid">
              {(["A", "B", "C", "D"] as const).map((k) => (
                <label key={k}>
                  <span className="field-label">{k}</span>
                  <NumInput
                    label={`ln K coefficient ${k}`}
                    value={r.keq?.[k]}
                    placeholder="0"
                    onChange={(v) => {
                      const next: Keq = { ...r.keq };
                      if (v === null) delete next[k];
                      else next[k] = v;
                      setKeq(next);
                    }}
                  />
                </label>
              ))}
            </div>
          )}

          <button className="danger" onClick={onDelete}>
            Delete reaction
          </button>
        </div>
      )}
    </li>
  );
}

// ---- list --------------------------------------------------------------------------------

export function ReactionsEditor({
  reactions,
  components,
  onChange,
}: {
  reactions: Reaction[];
  components: string[];
  onChange: (r: Reaction[]) => void;
}) {
  return (
    <div className="reactions">
      {reactions.length === 0 && (
        <p className="muted small">No reactions. Define them here and select them by id in a reactor.</p>
      )}
      <ul className="reaction-list">
        {reactions.map((r, i) => (
          <ReactionCard
            key={i}
            reaction={r}
            components={components}
            takenIds={reactions.filter((_, j) => j !== i).map((x) => x.id)}
            onChange={(next) => onChange(reactions.map((x, j) => (j === i ? next : x)))}
            onDelete={() => onChange(reactions.filter((_, j) => j !== i))}
          />
        ))}
      </ul>
      <button onClick={() => onChange([...reactions, { id: nextReactionId(reactions), stoichiometry: {} }])}>
        + Add reaction
      </button>
    </div>
  );
}
