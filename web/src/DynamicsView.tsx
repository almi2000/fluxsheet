// Runs the dynamic simulation (POST /api/dynamic) and plots the chosen time series, one
// chart per quantity, in the style of the sensitivity chart; CSV export of what is shown.
import { useMemo, useState } from "react";
import { api } from "./api";
import { ticks } from "./ProfileView";
import { displayUnit, toDisplay } from "./quantities";
import "./studies.css";
import type { DynamicResult, DynamicSeries, FlowsheetModel } from "./types";

const W = 460;
const H = 190;
const PAD = { l: 56, r: 12, t: 10, b: 30 };

/** Time axis in the unit that keeps the numbers readable. */
function timeAxis(tEnd: number): [number, string] {
  if (tEnd >= 3 * 3600) return [3600, "h"];
  if (tEnd >= 180) return [60, "min"];
  return [1, "s"];
}

function name(s: DynamicSeries): string {
  const unit = displayUnit(s.quantity);
  return `${s.owner} ${s.label}${unit ? ` [${unit}]` : ""}`;
}

function csvCell(text: string): string {
  return /[",\n]/.test(text) ? `"${text.replace(/"/g, '""')}"` : text;
}

function Chart({ time, series, scale, unit }: { time: number[]; series: DynamicSeries[]; scale: number; unit: string }) {
  const xs = time.map((t) => t / scale);
  const lines = series.map((s) => s.values.map((v) => toDisplay(v, s.quantity)));
  const all = lines.flat();
  const [x0, x1] = [xs[0], xs[xs.length - 1]];
  let [y0, y1] = [Math.min(...all), Math.max(...all)];
  if (y1 - y0 < 1e-9 * Math.max(Math.abs(y0), 1)) [y0, y1] = [y0 - 1e-3 * Math.max(Math.abs(y0), 1), y1 + 1e-3 * Math.max(Math.abs(y1), 1)];
  const sx = (x: number) => PAD.l + ((x - x0) / (x1 - x0 || 1)) * (W - PAD.l - PAD.r);
  const sy = (y: number) => H - PAD.b - ((y - y0) / (y1 - y0)) * (H - PAD.t - PAD.b);
  const yUnit = displayUnit(series[0].quantity);
  return (
    <div className="sens-chart">
      <svg viewBox={`0 0 ${W} ${H}`} className="profile-chart" role="img" aria-label={series.map(name).join(", ")}>
        {ticks(y0, y1).map((t) => (
          <g key={`y${t}`}>
            <line x1={PAD.l} x2={W - PAD.r} y1={sy(t)} y2={sy(t)} className="grid" />
            <text x={PAD.l - 4} y={sy(t) + 3} textAnchor="end">
              {Number(t.toPrecision(5))}
            </text>
          </g>
        ))}
        {ticks(x0, x1).map((t) => (
          <text key={`x${t}`} x={sx(t)} y={H - PAD.b + 13} textAnchor="middle">
            {Number(t.toPrecision(4))}
          </text>
        ))}
        <text x={(PAD.l + W - PAD.r) / 2} y={H - 3} textAnchor="middle" className="axis-label">
          Time [{unit}]
        </text>
        {yUnit && (
          <text x={4} y={PAD.t + 8} className="axis-label">
            {yUnit}
          </text>
        )}
        {lines.map((ys, gi) => (
          <g key={series[gi].key} className={`series s${gi % 6}`}>
            <path d={ys.map((y, i) => `${i ? "L" : "M"}${sx(xs[i]).toFixed(1)},${sy(y).toFixed(1)}`).join("")} className="line" />
          </g>
        ))}
      </svg>
      <div className="sens-legend">
        {series.map((s, gi) => (
          <span key={s.key} className={`series s${gi % 6}`}>
            <i /> {s.owner} {s.label}
          </span>
        ))}
      </div>
    </div>
  );
}

/** Series shown by default: controllers' measurement and output, levels and temperatures. */
function defaults(r: DynamicResult): string[] {
  const keys = r.series
    .filter((s) =>
      s.group === "Controllers" ? s.label !== "Output" : s.group === "Units" && ["level", "T"].includes(s.key.split(".").pop()!),
    )
    .map((s) => s.key);
  return keys.length ? keys.slice(0, 6) : r.series.slice(0, 2).map((s) => s.key);
}

export function DynamicsView({
  getModel,
  onRun,
}: {
  getModel: () => FlowsheetModel;
  /** A run finished, for the flowsheet `model` (the report can include it). */
  onRun?: (model: FlowsheetModel, result: DynamicResult) => void;
}) {
  const [running, setRunning] = useState(false);
  const [result, setResult] = useState<DynamicResult | null>(null);
  const [shown, setShown] = useState<string[]>([]);
  const [filter, setFilter] = useState("");
  const model = getModel();
  const dyn = model.dynamics;

  const run = async () => {
    setRunning(true);
    try {
      const sent = getModel();
      const r = await api.dynamic(sent);
      setResult(r);
      onRun?.(sent, r);
      if (r.series.length) setShown((prev) => (prev.some((k) => r.series.some((s) => s.key === k)) ? prev : defaults(r)));
    } catch (err) {
      setResult({ ok: false, time: [], series: [], warnings: [], error: { message: (err as Error).message, where: null, hint: null } });
    } finally {
      setRunning(false);
    }
  };

  const byKey = useMemo(() => new Map((result?.series ?? []).map((s) => [s.key, s])), [result]);
  const visible = shown.map((k) => byKey.get(k)).filter((s): s is DynamicSeries => !!s);
  // One chart per quantity, so a level and a temperature never share an axis.
  const charts = useMemo(() => {
    const groups = new Map<string, DynamicSeries[]>();
    for (const s of visible) {
      const q = s.quantity ?? "";
      if (!groups.has(q)) groups.set(q, []);
      groups.get(q)!.push(s);
    }
    return [...groups.values()];
  }, [visible]);

  if (!dyn)
    return (
      <p className="muted small">
        No dynamic simulation set up. Click the empty canvas to open the flowsheet settings and switch on hold-ups,
        controllers and events under <b>Dynamics</b>.
      </p>
    );

  const exportCsv = () => {
    const cols = visible.length ? visible : result!.series;
    const lines = [["time [s]", ...cols.map(name)].map(csvCell).join(",")];
    result!.time.forEach((t, i) => lines.push([t, ...cols.map((s) => toDisplay(s.values[i], s.quantity))].join(",")));
    const blob = new Blob([lines.join("\n") + "\n"], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${model.name.replace(/[^\w-]+/g, "_") || "dynamic"}.csv`;
    a.click();
    URL.revokeObjectURL(a.href);
  };

  const has = result && result.time.length > 1;
  const [scale, unit] = has ? timeAxis(result.time[result.time.length - 1]) : [1, "s"];
  const toggle = (key: string) => setShown((prev) => (prev.includes(key) ? prev.filter((k) => k !== key) : [...prev, key]));
  const matches = (result?.series ?? []).filter((s) => `${s.owner} ${s.label}`.toLowerCase().includes(filter.toLowerCase()));
  const groups = ["Controllers", "Units", "Streams"] as const;
  return (
    <div className="sensitivity">
      <div className="sens-toolbar">
        <button className="primary" onClick={run} disabled={running}>
          {running ? "Simulating…" : `▶ Run dynamics (${dyn.t_end})`}
        </button>
        {has && <button onClick={exportCsv}>Export CSV</button>}
        {result?.warnings.map((w) => (
          <span key={w} className="warn-text small">
            {w}
          </span>
        ))}
      </div>
      {result && !result.ok && result.error && (
        <div className="banner error">
          <b>{result.error.where ? `${result.error.where}: ` : ""}</b>
          {result.error.message}
          {result.error.hint && <div className="hint">💡 {result.error.hint}</div>}
          {has && <div className="hint">The charts show the run up to the failure.</div>}
        </div>
      )}
      {has && (
        <div className="sens-body">
          <div className="dyn-charts">
            {charts.length === 0 && <p className="muted small">Pick series on the right to plot them.</p>}
            {charts.map((series) => (
              <Chart key={series[0].key} time={result.time} series={series} scale={scale} unit={unit} />
            ))}
          </div>
          <div className="dyn-picker">
            <input type="search" placeholder="Filter series…" aria-label="Filter series" value={filter} onChange={(e) => setFilter(e.target.value)} />
            {groups.map((g) => {
              const items = matches.filter((s) => s.group === g);
              if (!items.length) return null;
              return (
                <fieldset key={g}>
                  <legend>{g}</legend>
                  {items.map((s) => (
                    <label key={s.key} className="check">
                      <input type="checkbox" checked={shown.includes(s.key)} onChange={() => toggle(s.key)} />
                      {name(s)}
                    </label>
                  ))}
                </fieldset>
              );
            })}
          </div>
        </div>
      )}
    </div>
  );
}
