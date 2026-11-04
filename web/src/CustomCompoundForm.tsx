import { useState } from "react";
import { EMPTY_DRAFT, fromDraft, toDraft, type CompoundDraft } from "./customCompounds";
import type { CustomCompound } from "./types";
import "./custom-compounds.css";

/**
 * Define a compound the database lacks: name, formula and groups, optional known constants.
 * Missing properties are estimated on the server (Joback, Lee-Kesler) and flagged in results.
 */
export function CustomCompoundForm({
  components,
  editing,
  onSave,
  onClose,
}: {
  components: string[];
  editing: CustomCompound | null; // null: a new compound
  onSave: (c: CustomCompound, previousName: string | null) => void;
  onClose: () => void;
}) {
  const [draft, setDraft] = useState<CompoundDraft>(editing ? toDraft(editing) : EMPTY_DRAFT);
  const [error, setError] = useState<string | null>(null);
  const set = (patch: Partial<CompoundDraft>) => {
    setDraft((d) => ({ ...d, ...patch }));
    setError(null);
  };
  const taken = components.filter((c) => c !== editing?.name);

  function save() {
    const result = fromDraft(draft, taken);
    if (typeof result === "string") return setError(result);
    // Fields the form does not show (set in the file) survive an edit.
    onSave({ vapor_pressure: editing?.vapor_pressure, antoine: editing?.antoine, ...result }, editing?.name ?? null);
  }

  const field = (label: string, key: keyof CompoundDraft, placeholder: string, multiline = false) => (
    <label className="cc-field">
      <span>{label}</span>
      {multiline ? (
        <textarea rows={2} value={draft[key]} placeholder={placeholder} onChange={(e) => set({ [key]: e.target.value })} />
      ) : (
        <input
          value={draft[key]}
          placeholder={placeholder}
          // Feeds, reactions and specs refer to the name: rename by removing and adding again.
          disabled={key === "name" && editing !== null}
          title={key === "name" && editing ? "To rename, remove the compound and add it again" : undefined}
          onChange={(e) => set({ [key]: e.target.value })}
        />
      )}
    </label>
  );

  return (
    <div className="custom-compound-form" role="group" aria-label="Custom compound">
      {field("Name", "name", "my-ester")}
      {field("Formula", "formula", "C5H10O2")}
      {field("SMILES (optional)", "smiles", "CCCC(=O)OC")}
      {field("Joback groups", "joback", "-CH3: 2, -CH2-: 2, -COO- (ester): 1", true)}
      {field("UNIFAC groups (for UNIFAC, UNIQUAC)", "unifac", "CH3: 2, CH2: 1, CH2COO: 1", true)}
      {field("Known constants (optional)", "constants", "Tb: 102.7 degC, Tc: 554.5 K, Pc: 34.7 bar, omega: 0.38", true)}
      <p className="muted small">
        Missing constants are estimated: Tb, Tc, Pc, Hf, Gf and Cp by Joback groups, omega and the vapor pressure by
        Lee-Kesler. Results flag them as estimates.
      </p>
      {error && <p className="callout warn">{error}</p>}
      <div className="cc-actions">
        <button className="primary" onClick={save}>
          {editing ? "Update compound" : "Add compound"}
        </button>
        <button onClick={onClose}>Cancel</button>
      </div>
    </div>
  );
}
