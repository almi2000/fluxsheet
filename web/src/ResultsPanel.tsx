import { fmt } from "./quantities";
import type { SolveResult } from "./types";

export function ResultsPanel({
  result,
  stale,
  onSelectStream,
}: {
  result: SolveResult | null;
  stale: boolean;
  onSelectStream: (id: string) => void;
}) {
  if (!result) {
    return (
      <div className="results empty muted">
        Build a flowsheet and press <b>Run</b>. Stream results appear here.
      </div>
    );
  }
  if (!result.ok) return null;
  const ids = Object.keys(result.streams);
  const components = ids.length ? Object.keys(result.streams[ids[0]].component_flows) : [];
  const rows: [string, (id: string) => string][] = [
    ["Temperature [°C]", (id) => fmt.T(result.streams[id].T)],
    ["Pressure [bar]", (id) => fmt.P(result.streams[id].P)],
    ["Vapor fraction", (id) => fmt.frac(result.streams[id].vapor_fraction)],
    ["Molar flow [kmol/h]", (id) => fmt.molar(result.streams[id].molar_flow)],
    ["Mass flow [kg/h]", (id) => fmt.mass(result.streams[id].mass_flow)],
    ...components.map(
      (c) => [`  ${c} [mol frac.]`, (id: string) =>
        result.streams[id].molar_flow > 0 ? fmt.frac(result.streams[id].mole_fractions[c]) : "–"] as [string, (id: string) => string],
    ),
  ];
  return (
    <div className={`results ${stale ? "stale" : ""}`}>
      <div className="results-meta">
        <span>
          Method <b>{result.property_method}</b>
        </span>
        {result.tear_streams.length > 0 && (
          <span>
            Recycle via <b>{result.tear_streams.join(", ")}</b> converged in <b>{result.iterations}</b> iterations
          </span>
        )}
        {stale && <span className="warn-text">Flowsheet changed since this run</span>}
      </div>
      <div className="table-scroll">
        <table>
          <thead>
            <tr>
              <th />
              {ids.map((id) => (
                <th key={id}>
                  <button className="link" onClick={() => onSelectStream(id)}>
                    {id}
                  </button>
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map(([label, get]) => (
              <tr key={label}>
                <th className={label.startsWith("  ") ? "sub" : ""}>{label.trim()}</th>
                {ids.map((id) => (
                  <td key={id}>{get(id)}</td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
