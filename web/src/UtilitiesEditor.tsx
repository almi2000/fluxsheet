// Flowsheet-level utilities: heating and cooling media that units refer to by id (see
// utilities.py). Presets fill in typical values; every value stays editable.
import { useEffect, useState } from "react";
import { QuantityInput } from "./Inspector";
import type { Utility } from "./types";

const PRESETS: { label: string; utility: Utility }[] = [
  { label: "Cooling water", utility: { id: "CW", type: "fluid", T_in: "25 degC", T_out: "35 degC", cp: "4.18 kJ/kgK" } },
  { label: "Chilled water", utility: { id: "CHW", type: "fluid", T_in: "5 degC", T_out: "15 degC", cp: "4.2 kJ/kgK" } },
  { label: "Low-pressure steam (3 bar)", utility: { id: "LPS", type: "steam", P: "3 bar" } },
  { label: "Medium-pressure steam (10 bar)", utility: { id: "MPS", type: "steam", P: "10 bar" } },
  { label: "High-pressure steam (40 bar)", utility: { id: "HPS", type: "steam", P: "40 bar" } },
  { label: "Thermal oil", utility: { id: "OIL", type: "fluid", T_in: "300 degC", T_out: "250 degC", cp: "2.3 kJ/kgK" } },
];

function freshId(stem: string, taken: string[]): string {
  if (!taken.includes(stem)) return stem;
  let i = 2;
  while (taken.includes(`${stem}${i}`)) i++;
  return `${stem}${i}`;
}

function summary(u: Utility): string {
  if (u.type === "steam") return `steam, saturated at ${u.P ?? "?"}`;
  return `${u.T_in ?? "?"} → ${u.T_out ?? "?"}`.replace(/degC/g, "°C");
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

function UtilityCard({
  utility: u,
  takenIds,
  onChange,
  onDelete,
}: {
  utility: Utility;
  takenIds: string[];
  onChange: (u: Utility) => void;
  onDelete: () => void;
}) {
  const [open, setOpen] = useState(false);
  const [draftId, setDraftId] = useState(u.id);
  useEffect(() => setDraftId(u.id), [u.id]);
  const idError = draftId.trim() === "" ? "Give it an id" : takenIds.includes(draftId.trim()) ? "Id already used" : null;
  const set = (key: keyof Utility, v: unknown) => {
    const next = { ...u } as Record<string, unknown>;
    if (v === null || v === undefined || v === "") delete next[key];
    else next[key] = v;
    onChange(next as unknown as Utility);
  };

  return (
    <li className="reaction-card">
      <button className="reaction-summary" onClick={() => setOpen(!open)} aria-expanded={open}>
        <b>{u.id}</b>
        <span className="small">{summary(u)}</span>
      </button>
      {open && (
        <div className="reaction-body">
          <Row label="Id">
            <input
              value={draftId}
              onChange={(e) => setDraftId(e.target.value)}
              onBlur={() => !idError && draftId.trim() !== u.id && onChange({ ...u, id: draftId.trim() })}
            />
            {idError && <span className="field-error">{idError}</span>}
          </Row>
          <Row label="Kind">
            <select
              value={u.type}
              onChange={(e) =>
                onChange(
                  e.target.value === "steam"
                    ? { id: u.id, type: "steam", P: "3 bar" }
                    : { id: u.id, type: "fluid", T_in: "25 degC", T_out: "35 degC", cp: "4.18 kJ/kgK" },
                )
              }
            >
              <option value="fluid">Fluid (warms up or cools down)</option>
              <option value="steam">Steam (condenses)</option>
            </select>
          </Row>
          {u.type === "steam" ? (
            <Row label="Pressure" help="Saturated steam; its temperature and heat of condensation follow (IAPWS-95).">
              <QuantityInput quantity="pressure" value={u.P} onChange={(v) => set("P", v)} />
            </Row>
          ) : (
            <>
              <Row label="Inlet temperature">
                <QuantityInput quantity="temperature" value={u.T_in} onChange={(v) => set("T_in", v)} />
              </Row>
              <Row label="Outlet temperature" help="Sets the flow a unit needs. A unit can fix its flow instead.">
                <QuantityInput quantity="temperature" value={u.T_out} onChange={(v) => set("T_out", v)} />
              </Row>
              <Row label="Heat capacity">
                <QuantityInput quantity="specific_heat" value={u.cp} placeholder="4.18" onChange={(v) => set("cp", v)} />
              </Row>
            </>
          )}
          <button className="danger" onClick={onDelete}>
            Delete utility
          </button>
        </div>
      )}
    </li>
  );
}

export function UtilitiesEditor({ utilities, onChange }: { utilities: Utility[]; onChange: (u: Utility[]) => void }) {
  return (
    <div className="reactions">
      {utilities.length === 0 && (
        <p className="muted small">No utilities. Add heating and cooling media here and choose them in a unit.</p>
      )}
      <ul className="reaction-list">
        {utilities.map((u, i) => (
          <UtilityCard
            key={i}
            utility={u}
            takenIds={utilities.filter((_, j) => j !== i).map((x) => x.id)}
            onChange={(next) => onChange(utilities.map((x, j) => (j === i ? next : x)))}
            onDelete={() => onChange(utilities.filter((_, j) => j !== i))}
          />
        ))}
      </ul>
      <select
        value=""
        onChange={(e) => {
          const preset = PRESETS[Number(e.target.value)];
          if (!preset) return;
          const id = freshId(preset.utility.id, utilities.map((u) => u.id));
          onChange([...utilities, { ...preset.utility, id }]);
        }}
      >
        <option value="">+ Add utility…</option>
        {PRESETS.map((p, i) => (
          <option key={p.label} value={i}>
            {p.label}
          </option>
        ))}
      </select>
    </div>
  );
}
