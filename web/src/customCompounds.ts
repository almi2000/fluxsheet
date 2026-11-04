// Text <-> data for the custom compound form; see src/fluxsheet/custom_compounds.py.
import type { CustomCompound } from "./types";

/** Constants a definition may give (the backend estimates the rest). */
export const CONSTANT_KEYS = ["MW", "Tb", "Tm", "Tc", "Pc", "Vc", "omega", "Hf", "Gf", "dipole", "r", "q"];

/** "name: value" entries separated by commas, semicolons or new lines. */
function entries(text: string): [string, string][] {
  return text
    .split(/[;\n]|,(?![^()]*\))/)
    .map((s) => s.trim())
    .filter(Boolean)
    .map((s) => {
      const at = Math.max(s.lastIndexOf(":"), s.lastIndexOf("="));
      if (at <= 0) throw new Error(`"${s}" needs the form name: value`);
      return [s.slice(0, at).trim(), s.slice(at + 1).trim()];
    });
}

/** Groups such as "-CH3: 2, -OH (alcohol): 1" to {"-CH3": 2, ...}. */
export function parseGroups(text: string): Record<string, number> {
  const out: Record<string, number> = {};
  for (const [name, value] of entries(text)) {
    const n = Number(value);
    if (!Number.isInteger(n) || n < 0) throw new Error(`count of ${name} must be a whole number, got "${value}"`);
    if (n) out[name] = (out[name] ?? 0) + n;
  }
  return out;
}

/** Constants such as "Tb: 78.4 degC, Pc: 61 bar, omega: 0.64"; numbers stay numbers (SI). */
export function parseConstants(text: string): Record<string, number | string> {
  const out: Record<string, number | string> = {};
  for (const [key, value] of entries(text)) {
    if (!CONSTANT_KEYS.includes(key)) throw new Error(`unknown constant "${key}"; known: ${CONSTANT_KEYS.join(", ")}`);
    if (!value) continue;
    out[key] = /^[-+]?\d*\.?\d+(e[-+]?\d+)?$/i.test(value) ? Number(value) : value;
  }
  return out;
}

export function formatEntries(values: Record<string, number | string> | undefined): string {
  return Object.entries(values ?? {})
    .map(([k, v]) => `${k}: ${v}`)
    .join(", ");
}

export interface CompoundDraft {
  name: string;
  formula: string;
  smiles: string;
  joback: string;
  unifac: string;
  constants: string;
}

export const EMPTY_DRAFT: CompoundDraft = { name: "", formula: "", smiles: "", joback: "", unifac: "", constants: "" };

export function toDraft(c: CustomCompound): CompoundDraft {
  return {
    name: c.name,
    formula: c.formula,
    smiles: c.smiles ?? "",
    joback: formatEntries(c.joback_groups),
    unifac: formatEntries(c.unifac_groups),
    constants: formatEntries(c.constants),
  };
}

/** A definition from the form, or an error message for what is missing or malformed. */
export function fromDraft(d: CompoundDraft, taken: string[]): CustomCompound | string {
  const name = d.name.trim();
  const formula = d.formula.trim();
  if (!name) return "give the compound a name";
  if (taken.some((t) => t.toLowerCase() === name.toLowerCase())) return `"${name}" is already a component`;
  if (!formula) return "give the molecular formula, e.g. C2H6O";
  if (!/^[A-Z][A-Za-z0-9()]*$/.test(formula)) return `formula "${formula}" should look like C2H6O`;
  try {
    const joback = parseGroups(d.joback);
    const unifac = parseGroups(d.unifac);
    const constants = parseConstants(d.constants);
    const out: CustomCompound = { name, formula };
    if (d.smiles.trim()) out.smiles = d.smiles.trim();
    if (Object.keys(joback).length) out.joback_groups = joback;
    if (Object.keys(unifac).length) out.unifac_groups = unifac;
    if (Object.keys(constants).length) out.constants = constants;
    const known = ["Tb", "Tc", "Pc"].every((k) => k in constants);
    if (!out.joback_groups && !out.smiles && !known)
      return "give Joback groups (or a SMILES) to estimate the properties, or at least Tb, Tc and Pc";
    return out;
  } catch (e) {
    return (e as Error).message;
  }
}
