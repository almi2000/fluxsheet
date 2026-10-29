// Property sections in the editor: a per-unit method override, the method tag on nodes, and
// the methods summary in the results. See src/fluxsheet/property_sections.py.
import type { Method, SolveResult, UnitResult } from "./types";

const SHORT: Record<string, string> = {
  "peng-robinson": "PR", srk: "SRK", unifac: "UNIFAC", nrtl: "NRTL", uniquac: "UNIQUAC",
  wilson: "Wilson", ideal: "Ideal",
};

export const shortMethod = (key: string) => SHORT[key] ?? key;

type Suggestion = { method: string; reason: string };

/** Small coloured tag under a node's id: which property method the unit calculates with. */
export function MethodTag({ method, override }: { method: string; override: boolean }) {
  return (
    <div
      className={`method-tag m-${method}${override ? " override" : ""}`}
      title={`${shortMethod(method)}${override ? " (set on this unit)" : " (flowsheet default)"}`}
    >
      {shortMethod(method)}
    </div>
  );
}

/** The "Property method" section of the unit inspector. */
export function UnitMethodField({
  value,
  flowsheetMethod,
  methods,
  result,
  isProduct,
  isSeparation,
  onChange,
}: {
  value: string | undefined;
  flowsheetMethod: string; // the resolved default (e.g. "auto" already replaced when solved)
  methods: Method[];
  result?: UnitResult;
  isProduct: boolean;
  isSeparation: boolean;
  onChange: (method: string | undefined) => void;
}) {
  const label = (key: string) => methods.find((m) => m.key === key)?.label ?? key;
  const suggestion = result?.property_method_suggestion as Suggestion | undefined;
  const effective = (result?.property_method as string | undefined) ?? value ?? flowsheetMethod;
  const defaultLabel = isProduct
    ? "Same as the unit feeding it"
    : `Flowsheet default (${flowsheetMethod === "auto" ? "automatic" : label(flowsheetMethod)})`;
  return (
    <section className="unit-method">
      <h3>Property method</h3>
      <select value={value ?? ""} onChange={(e) => onChange(e.target.value || undefined)}>
        <option value="">{defaultLabel}</option>
        {methods.map((m) => (
          <option key={m.key} value={m.key}>
            {m.label}
            {suggestion?.method === m.key ? " (suggested)" : ""}
          </option>
        ))}
      </select>
      {result?.property_method !== undefined && (
        <div className="kv">
          <span>Calculated with</span>
          <b>{label(effective)}</b>
        </div>
      )}
      {suggestion && suggestion.method === effective && (
        <p className="callout">
          {label(suggestion.method)} suits what this unit sees. {suggestion.reason}
        </p>
      )}
      {suggestion && suggestion.method !== effective && (
        <p className="callout warn">
          For this unit's feed, {label(suggestion.method)} is usually the better choice. {suggestion.reason}{" "}
          <button className="link" onClick={() => onChange(suggestion.method)}>
            Use {label(suggestion.method)}
          </button>
        </p>
      )}
      {!suggestion && isSeparation && !result && (
        <p className="muted small">Run the flowsheet to get a suggestion based on this unit's own feed.</p>
      )}
    </section>
  );
}

export type Crossing = {
  stream: string;
  from_unit: string;
  to_unit: string;
  from_method: string;
  to_method: string;
  vapor_fraction_before: number;
  vapor_fraction_after: number;
  enthalpy_change: number;
};

/** "Method PR" for one method; "Methods PR (8 units) · NRTL (DEGASSER, COLUMN)" for several. */
export function MethodsSummary({ result }: { result: SolveResult }) {
  const byMethod: Record<string, string[]> = {};
  for (const [id, u] of Object.entries(result.units)) {
    const m = u.property_method as string | undefined;
    if (m && u.type !== "Product") (byMethod[m] ??= []).push(id);
  }
  const entries = Object.entries(byMethod);
  if (entries.length <= 1) {
    return (
      <span>
        Method <b>{result.property_method}</b>
      </span>
    );
  }
  const crossings = (result.property_crossings ?? []) as Crossing[];
  return (
    <>
      <span>
        Methods{" "}
        {entries.map(([m, ids], i) => (
          <span key={m}>
            {i > 0 && " · "}
            <b className={`method-text m-${m}`}>{shortMethod(m)}</b>{" "}
            ({ids.length > 3 ? `${ids.length} units` : ids.join(", ")})
          </span>
        ))}
      </span>
      {crossings.length > 0 && (
        <span
          title={crossings
            .map(
              (c) =>
                `${c.stream}: ${shortMethod(c.from_method)} → ${shortMethod(c.to_method)}, vapor fraction ` +
                `${c.vapor_fraction_before.toFixed(3)} → ${c.vapor_fraction_after.toFixed(3)}, ` +
                `ΔH ${(c.enthalpy_change / 1000).toFixed(1)} kW`,
            )
            .join("\n")}
        >
          Re-flashed at method boundaries: <b>{crossings.map((c) => c.stream).join(", ")}</b>
        </span>
      )}
    </>
  );
}
