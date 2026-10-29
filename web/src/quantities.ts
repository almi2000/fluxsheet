// Display units per quantity. The first entry is the default for new inputs; the
// backend accepts every unit listed here (see src/fluxsheet/quantities.py).
export const UNITS: Record<string, { value: string; label: string }[]> = {
  temperature: [
    { value: "degC", label: "°C" },
    { value: "K", label: "K" },
    { value: "degF", label: "°F" },
  ],
  pressure: [
    { value: "bar", label: "bar" },
    { value: "kPa", label: "kPa" },
    { value: "MPa", label: "MPa" },
    { value: "atm", label: "atm" },
    { value: "psi", label: "psi" },
    { value: "Pa", label: "Pa" },
  ],
  pressure_difference: [
    { value: "bar", label: "bar" },
    { value: "kPa", label: "kPa" },
    { value: "mbar", label: "mbar" },
    { value: "Pa", label: "Pa" },
  ],
  molar_flow: [
    { value: "kmol/h", label: "kmol/h" },
    { value: "mol/s", label: "mol/s" },
    { value: "kmol/s", label: "kmol/s" },
  ],
  power: [
    { value: "kW", label: "kW" },
    { value: "W", label: "W" },
    { value: "MW", label: "MW" },
  ],
  volume: [
    { value: "m3", label: "m³" },
    { value: "L", label: "L" },
  ],
  length: [
    { value: "m", label: "m" },
    { value: "cm", label: "cm" },
    { value: "mm", label: "mm" },
  ],
  area: [
    { value: "m2", label: "m²" },
    { value: "cm2", label: "cm²" },
  ],
  time: [
    { value: "s", label: "s" },
    { value: "min", label: "min" },
    { value: "h", label: "h" },
  ],
  molar_energy: [
    { value: "kJ/mol", label: "kJ/mol" },
    { value: "J/mol", label: "J/mol" },
    { value: "kcal/mol", label: "kcal/mol" },
  ],
};

const SI_UNIT: Record<string, string> = {
  temperature: "K",
  pressure: "Pa",
  pressure_difference: "Pa",
  molar_flow: "mol/s",
  power: "W",
  volume: "m3",
  length: "m",
  area: "m2",
  time: "s",
  molar_energy: "J/mol",
};

/** Split a stored value ("80 degC" or 353.15) into number text and unit. */
export function splitValue(value: unknown, quantity: string): { num: string; unit: string } {
  if (value === null || value === undefined || value === "") {
    return { num: "", unit: UNITS[quantity]?.[0].value ?? "" };
  }
  if (typeof value === "number") return { num: String(value), unit: SI_UNIT[quantity] ?? "" };
  const m = String(value).match(/^\s*([-+]?[\d.eE+-]+)\s*(.*)$/);
  if (!m) return { num: String(value), unit: UNITS[quantity]?.[0].value ?? "" };
  const unit = m[2] || SI_UNIT[quantity] || "";
  const known = UNITS[quantity]?.find((u) => u.value.toLowerCase() === unit.toLowerCase());
  return { num: m[1], unit: known ? known.value : unit };
}

export function joinValue(num: string, unit: string): string | number | null {
  if (num.trim() === "") return null;
  return unit ? `${num.trim()} ${unit}` : Number(num);
}

export const fmt = {
  T: (k: number) => (k - 273.15).toFixed(1),
  P: (pa: number) => trim(pa / 1e5, 4),
  molar: (mols: number) => trim(mols * 3.6, 4), // kmol/h
  mass: (kgs: number) => trim(kgs * 3600, 4), // kg/h
  kW: (w: number) => trim(w / 1e3, 4),
  frac: (x: number) => x.toFixed(4),
};

// How results in SI are shown: [factor from SI, offset, unit label].
const DISPLAY: Record<string, [number, number, string]> = {
  temperature: [1, -273.15, "°C"],
  pressure: [1e-5, 0, "bar"],
  pressure_difference: [1e-5, 0, "bar"],
  molar_flow: [3.6, 0, "kmol/h"],
  mass_flow: [3600, 0, "kg/h"],
  power: [1e-3, 0, "kW"],
  volume: [1, 0, "m³"],
  length: [1, 0, "m"],
  area: [1, 0, "m²"],
  time: [1, 0, "s"],
  molar_energy: [1e-3, 0, "kJ/mol"],
};

/** Display unit label for a quantity ("" for dimensionless). */
export function displayUnit(quantity: string | null): string {
  return (quantity && DISPLAY[quantity]?.[2]) || "";
}

/** Convert an SI value to its display unit (number only). */
export function toDisplay(value: number, quantity: string | null): number {
  const d = quantity ? DISPLAY[quantity] : undefined;
  return d ? value * d[0] + d[1] : value;
}

/** Format an SI value with its display unit, e.g. 353.15 K -> "80.0 °C". */
export function formatQuantity(value: unknown, quantity: string | null): string {
  if (typeof value !== "number" || !Number.isFinite(value)) return String(value ?? "–");
  if (quantity === "temperature") return `${fmt.T(value)} °C`;
  const shown = toDisplay(value, quantity);
  const unit = displayUnit(quantity);
  return unit ? `${trim(shown, 4)} ${unit}` : quantity === "fraction" ? shown.toFixed(4) : trim(shown, 4);
}

function trim(x: number, digits: number): string {
  if (x === 0) return "0";
  const abs = Math.abs(x);
  if (abs >= 1e5 || abs < 1e-3) return x.toExponential(2);
  return Number(x.toPrecision(digits)).toLocaleString("en-US", { maximumFractionDigits: 3 });
}
