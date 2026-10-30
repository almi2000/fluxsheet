// Runs a sensitivity study (POST /api/sensitivity) and shows it as a line chart, in the style
// of ProfileView, plus a table that can be exported as CSV.
import { useState } from "react";
import { api } from "./api";
import { ticks } from "./ProfileView";
import { displayUnit, formatQuantity, toDisplay } from "./quantities";
import "./studies.css";
import type { FlowsheetModel, SensitivityColumn, SensitivityResult, SensitivityStudy } from "./types";

const W = 460;
const H = 190;
const PAD = { l: 52, r: 12, t: 10, b: 30 };

function heading(c: SensitivityColumn): string {
  const unit = displayUnit(c.quantity);
  return unit ? `${c.label} [${unit}]` : c.label;
}

function show(v: number | null, quantity: string | null): string {
  if (v === null) return "–";
  const d = toDisplay(v, quantity);
  return quantity === "fraction" ? d.toFixed(5) : Number(d.toPrecision(6)).toString();
}

function csvCell(text: string): string {
  return /[",\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}

function toCsv(r: SensitivityResult): string {
  const cols = [...r.inputs!, ...r.outputs!];
  const lines = [[...cols.map(heading), "error"].map(csvCell).join(",")];
  for (const row of r.rows!) {
    const values = [...row.inputs, ...row.outputs].map((v, i) => (v === null ? "" : String(toDisplay(v, cols[i].quantity))));
    lines.push([...values, row.error?.message ?? ""].map(csvCell).join(","));
  }
  return lines.join("\n") + "\n";
}

function Chart({ result, output }: { result: SensitivityResult; output: number }) {
  const [xc, gc] = result.inputs!;
  const yc = result.outputs![output];
  // One line per value of the second parameter (if any).
  const groups = new Map<number, { x: number; y: number }[]>();
  for (const row of result.rows!) {
    const key = gc ? row.inputs[1] : 0;
    if (!groups.has(key)) groups.set(key, []);
    const y = row.outputs[output];
    if (y !== null) groups.get(key)!.push({ x: toDisplay(row.inputs[0], xc.quantity), y: toDisplay(y, yc.quantity) });
  }
  const pts = [...groups.values()].flat();
  if (pts.length === 0) return <p className="muted small">No point of this study solved.</p>;
  const xs = pts.map((p) => p.x);
  const ys = pts.map((p) => p.y);
  const [x0, x1] = [Math.min(...xs), Math.max(...xs)];
  let [y0, y1] = [Math.min(...ys), Math.max(...ys)];
  if (y0 === y1) [y0, y1] = [y0 - 1, y1 + 1];
  const sx = (x: number) => PAD.l + ((x - x0) / (x1 - x0 || 1)) * (W - PAD.l - PAD.r);
  const sy = (y: number) => H - PAD.b - ((y - y0) / (y1 - y0)) * (H - PAD.t - PAD.b);
  return (
    <div className="sens-chart">
      <svg viewBox={`0 0 ${W} ${H}`} className="profile-chart" role="img" aria-label={`${yc.label} vs ${xc.label}`}>
        {ticks(y0, y1).map((t) => (
          <g key={`y${t}`}>
            <line x1={PAD.l} x2={W - PAD.r} y1={sy(t)} y2={sy(t)} className="grid" />
            <text x={PAD.l - 4} y={sy(t) + 3} textAnchor="end">
              {Number(t.toPrecision(4))}
            </text>
          </g>
        ))}
        {ticks(x0, x1).map((t) => (
          <text key={`x${t}`} x={sx(t)} y={H - PAD.b + 13} textAnchor="middle">
            {Number(t.toPrecision(4))}
          </text>
        ))}
        <text x={(PAD.l + W - PAD.r) / 2} y={H - 3} textAnchor="middle" className="axis-label">
          {heading(xc)}
        </text>
        {[...groups.entries()].map(([key, line], gi) => (
          <g key={key} className={`series s${gi % 6}`}>
            <path d={line.map((p, i) => `${i ? "L" : "M"}${sx(p.x).toFixed(1)},${sy(p.y).toFixed(1)}`).join("")} className="line" />
            {line.map((p, i) => (
              <circle key={i} cx={sx(p.x)} cy={sy(p.y)} r={2.4} className="dot">
                <title>{`${xc.label} ${Number(p.x.toPrecision(5))}: ${Number(p.y.toPrecision(5))} ${displayUnit(yc.quantity)}`}</title>
              </circle>
            ))}
          </g>
        ))}
      </svg>
      {gc && (
        <div className="sens-legend">
          {[...groups.keys()].map((key, gi) => (
            <span key={key} className={`series s${gi % 6}`}>
              <i /> {gc.label} = {formatQuantity(key, gc.quantity)}
            </span>
          ))}
        </div>
      )}
    </div>
  );
}

export function SensitivityView({
  studies,
  getModel,
}: {
  studies: SensitivityStudy[];
  getModel: () => FlowsheetModel;
}) {
  const [studyId, setStudyId] = useState<string | null>(null);
  const [running, setRunning] = useState(false);
  const [result, setResult] = useState<SensitivityResult | null>(null);
  const [output, setOutput] = useState(0);

  if (studies.length === 0)
    return (
      <p className="muted small">
        No sensitivity study yet. Click the empty canvas to open the flowsheet settings and add one under{" "}
        <b>Sensitivity</b>.
      </p>
    );
  const study = studies.find((s) => s.id === studyId) ?? studies[0];
  const points = study.vary.reduce((n, v) => n * (v.points ?? 11), 1);

  const run = async () => {
    setRunning(true);
    try {
      setResult(await api.sensitivity(getModel(), study.id));
      setOutput(0);
    } catch (err) {
      setResult({ ok: false, error: { message: (err as Error).message, where: null, hint: null } });
    } finally {
      setRunning(false);
    }
  };

  const exportCsv = () => {
    const blob = new Blob([toCsv(result!)], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${result!.id ?? "sensitivity"}.csv`;
    a.click();
    URL.revokeObjectURL(a.href);
  };

  const ok = result?.ok && result.rows;
  const failed = ok ? result.rows!.filter((r) => r.error).length : 0;
  return (
    <div className="sensitivity">
      <div className="sens-toolbar">
        {studies.length > 1 && (
          <select aria-label="Study" value={study.id} onChange={(e) => setStudyId(e.target.value)}>
            {studies.map((s) => (
              <option key={s.id} value={s.id}>
                {s.id}
              </option>
            ))}
          </select>
        )}
        <button className="primary" onClick={run} disabled={running}>
          {running ? `Solving ${points} points…` : `▶ Run ${study.id} (${points} points)`}
        </button>
        {ok && result.outputs!.length > 1 && (
          <select aria-label="Output" value={output} onChange={(e) => setOutput(Number(e.target.value))}>
            {result.outputs!.map((o, i) => (
              <option key={o.label} value={i}>
                {heading(o)}
              </option>
            ))}
          </select>
        )}
        {ok && <button onClick={exportCsv}>Export CSV</button>}
        {ok && result.id !== study.id && <span className="muted small">showing {result.id}</span>}
        {failed > 0 && <span className="warn-text">{failed} points failed</span>}
      </div>
      {result && !result.ok && result.error && (
        <div className="banner error">
          <b>{result.error.where ? `${result.error.where}: ` : ""}</b>
          {result.error.message}
          {result.error.hint && <div className="hint">💡 {result.error.hint}</div>}
        </div>
      )}
      {ok && (
        <div className="sens-body">
          <Chart result={result} output={Math.min(output, result.outputs!.length - 1)} />
          <div className="table-scroll sens-table">
            <table>
              <thead>
                <tr>
                  {[...result.inputs!, ...result.outputs!].map((c) => (
                    <th key={c.label}>{heading(c)}</th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {result.rows!.map((row, i) => (
                  <tr key={i} title={row.error ? `${row.error.message}${row.error.hint ? ` (${row.error.hint})` : ""}` : undefined}>
                    {row.inputs.map((v, k) => (
                      <td key={`i${k}`} className="sens-input">
                        {show(v, result.inputs![k].quantity)}
                      </td>
                    ))}
                    {row.error ? (
                      <td colSpan={row.outputs.length} className="warn-text">
                        {row.error.message}
                      </td>
                    ) : (
                      row.outputs.map((v, k) => <td key={`o${k}`}>{show(v, result.outputs![k].quantity)}</td>)
                    )}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
    </div>
  );
}
