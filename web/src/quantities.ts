// Display units per quantity. The first entry is the default for new inputs; the
// backend accepts every unit listed here (see src/openflowsheet/quantities.py).
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
};

const SI_UNIT: Record<string, string> = {
  temperature: "K",
  pressure: "Pa",
  pressure_difference: "Pa",
  molar_flow: "mol/s",
  power: "W",
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

function trim(x: number, digits: number): string {
  if (x === 0) return "0";
  const abs = Math.abs(x);
  if (abs >= 1e5 || abs < 1e-3) return x.toExponential(2);
  return Number(x.toPrecision(digits)).toLocaleString("en-US", { maximumFractionDigits: 3 });
}
