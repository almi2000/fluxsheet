import { useState } from "react";
import { Diagrams } from "./Diagrams";
import { MethodsSummary, shortMethod } from "./PropertyMethods";
import { fmt, formatQuantity } from "./quantities";
import { SensitivityView } from "./SensitivityView";
import type { FlowsheetModel, SensitivityStudy, SolveResult } from "./types";

// Streams that split into two liquids carry this (see Stream.to_dict); null otherwise.
type LiquidSplit = { light: { fraction: number }; heavy: { fraction: number } } | null | undefined;

export function ResultsPanel({
  result,
  stale,
  onSelectStream,
  studies = [],
  getModel,
}: {
  result: SolveResult | null;
  stale: boolean;
  onSelectStream: (id: string) => void;
  studies?: SensitivityStudy[];
  getModel?: () => FlowsheetModel;
}) {
  const [view, setView] = useState<"streams" | "diagram" | "sensitivity">("streams");
  const [diagramOpened, setDiagramOpened] = useState(false);
  const [sensitivityOpened, setSensitivityOpened] = useState(false);
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
  const split = (id: string) => (result.streams[id] as unknown as { liquid_split?: LiquidSplit }).liquid_split;
  const rows: [string, (id: string) => string][] = [
    ["Temperature [°C]", (id) => fmt.T(result.streams[id].T)],
    ["Pressure [bar]", (id) => fmt.P(result.streams[id].P)],
    ["Vapor fraction", (id) => fmt.frac(result.streams[id].vapor_fraction)],
    ...(ids.some((id) => split(id))
      ? [["Heavy-liquid fraction", (id: string) => {
          const s = split(id);
          return s ? fmt.frac(s.heavy.fraction) : "–";
        }] as [string, (id: string) => string]]
      : []),
    ...(new Set(ids.map((id) => result.streams[id].property_method)).size > 1
      ? [["Property method", (id: string) => shortMethod(result.streams[id].property_method ?? "")] as
          [string, (id: string) => string]]
      : []),
    ["Molar flow [kmol/h]", (id) => fmt.molar(result.streams[id].molar_flow)],
    ["Mass flow [kg/h]", (id) => fmt.mass(result.streams[id].mass_flow)],
    ...components.map(
      (c) => [`  ${c} [mol frac.]`, (id: string) =>
        result.streams[id].molar_flow > 0 ? fmt.frac(result.streams[id].mole_fractions[c]) : "–"] as [string, (id: string) => string],
    ),
  ];
  // The composition of the first feed stream is the natural default for a phase envelope.
  const feed = ids.find((id) => result.streams[id].molar_flow > 0);
  return (
    <div className={`results ${stale ? "stale" : ""}`}>
      <div className="results-meta">
        <div className="segmented results-tabs">
          <button className={view === "streams" ? "on" : ""} onClick={() => setView("streams")}>
            Streams
          </button>
          <button
            className={view === "diagram" ? "on" : ""}
            onClick={() => {
              setView("diagram");
              setDiagramOpened(true);
            }}
          >
            Phase diagram
          </button>
          {getModel && (
            <button
              className={view === "sensitivity" ? "on" : ""}
              onClick={() => {
                setView("sensitivity");
                setSensitivityOpened(true);
              }}
            >
              Sensitivity
            </button>
          )}
        </div>
        <MethodsSummary result={result} />
        {result.tear_streams.length > 0 && (
          <span>
            Recycle via <b>{result.tear_streams.join(", ")}</b> converged in <b>{result.iterations}</b> iterations
          </span>
        )}
        {(result.design_specs ?? []).map((d) => (
          <span key={d.id} className="design-spec-summary" title={`${d.iterations} flowsheet solves`}>
            Spec {d.id}: {d.target} = <b>{formatQuantity(d.achieved, d.quantity)}</b> with {d.vary} ={" "}
            <b>{formatQuantity(d.vary_value, d.vary_quantity)}</b>
          </span>
        ))}
        {stale && <span className="warn-text">Flowsheet changed since this run</span>}
      </div>
      {sensitivityOpened && getModel && (
        <div hidden={view !== "sensitivity"}>
          <SensitivityView studies={studies} getModel={getModel} />
        </div>
      )}
      {diagramOpened && (
        <div hidden={view !== "diagram"}>
          <Diagrams
            components={components}
            method={result.property_method}
            fractions={feed ? result.streams[feed].mole_fractions : null}
          />
        </div>
      )}
      <div className="table-scroll" hidden={view !== "streams"}>
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
