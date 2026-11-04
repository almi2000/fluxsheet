// Runs an optimization (POST /api/optimize) and shows the objective per iteration, the optimum
// and the constraints; and the costing of a solved flowsheet as a table. See optimization.py.
import { useState } from "react";
import { api } from "./api";
import { ticks } from "./ProfileView";
import { formatQuantity } from "./quantities";
import "./studies.css";
import type { CostingResult, FlowsheetModel, Optimization, OptimizationResult } from "./types";

const W = 460;
const H = 190;
const PAD = { l: 64, r: 12, t: 10, b: 30 };

/** A number for display: money and plain values with up to 6 significant digits. */
export function money(v: number): string {
  return Math.abs(v) >= 1000 ? Math.round(v).toLocaleString("en-US") : Number(v.toPrecision(4)).toString();
}

function num(v: number): string {
  return Number(v.toPrecision(6)).toString();
}

function HistoryChart({ history }: { history: NonNullable<OptimizationResult["history"]> }) {
  const xs = history.map((h) => h.iteration);
  const ys = history.map((h) => h.objective);
  const [x0, x1] = [Math.min(...xs), Math.max(...xs)];
  let [y0, y1] = [Math.min(...ys), Math.max(...ys)];
  if (y0 === y1) [y0, y1] = [y0 - Math.abs(y0) * 0.01 - 1, y1 + Math.abs(y1) * 0.01 + 1];
  const sx = (x: number) => PAD.l + ((x - x0) / (x1 - x0 || 1)) * (W - PAD.l - PAD.r);
  const sy = (y: number) => H - PAD.b - ((y - y0) / (y1 - y0)) * (H - PAD.t - PAD.b);
  return (
    <div className="sens-chart">
      <svg viewBox={`0 0 ${W} ${H}`} className="profile-chart" role="img" aria-label="Objective per iteration">
        {ticks(y0, y1).map((t) => (
          <g key={`y${t}`}>
            <line x1={PAD.l} x2={W - PAD.r} y1={sy(t)} y2={sy(t)} className="grid" />
            <text x={PAD.l - 4} y={sy(t) + 3} textAnchor="end">
              {money(t)}
            </text>
          </g>
        ))}
        {ticks(x0, x1)
          .filter((t) => Number.isInteger(t))
          .map((t) => (
            <text key={`x${t}`} x={sx(t)} y={H - PAD.b + 13} textAnchor="middle">
              {t}
            </text>
          ))}
        <text x={(PAD.l + W - PAD.r) / 2} y={H - 3} textAnchor="middle" className="axis-label">
          Iteration
        </text>
        <g className="series s0">
          <path d={history.map((h, i) => `${i ? "L" : "M"}${sx(h.iteration).toFixed(1)},${sy(h.objective).toFixed(1)}`).join("")} className="line" />
          {history.map((h) => (
            <circle key={h.iteration} cx={sx(h.iteration)} cy={sy(h.objective)} r={2.6} className={`dot${h.max_violation > 1e-4 ? " infeasible" : ""}`}>
              <title>{`Iteration ${h.iteration}: ${num(h.objective)}${h.max_violation > 1e-4 ? ` (violates constraints by ${h.max_violation.toPrecision(2)})` : ""}`}</title>
            </circle>
          ))}
        </g>
      </svg>
    </div>
  );
}

export function OptimizationView({
  optimizations,
  getModel,
  onApply,
}: {
  optimizations: Optimization[];
  getModel: () => FlowsheetModel;
  onApply?: (values: NonNullable<OptimizationResult["apply"]>) => void;
}) {
  const [optId, setOptId] = useState<string | null>(null);
  const [running, setRunning] = useState(false);
  const [result, setResult] = useState<OptimizationResult | null>(null);
  const [applied, setApplied] = useState(false);

  if (optimizations.length === 0)
    return (
      <p className="muted small">
        No optimization yet. Click the empty canvas to open the flowsheet settings and add one under <b>Optimization</b>.
      </p>
    );
  const opt = optimizations.find((o) => o.id === optId) ?? optimizations[0];

  const run = async () => {
    setRunning(true);
    setApplied(false);
    try {
      setResult(await api.optimize(getModel(), opt.id));
    } catch (err) {
      setResult({ ok: false, error: { message: (err as Error).message, where: null, hint: null } });
    } finally {
      setRunning(false);
    }
  };

  const r = result;
  return (
    <div className="sensitivity">
      <div className="sens-toolbar">
        {optimizations.length > 1 && (
          <select aria-label="Optimization" value={opt.id} onChange={(e) => setOptId(e.target.value)}>
            {optimizations.map((o) => (
              <option key={o.id} value={o.id}>
                {o.id}
              </option>
            ))}
          </select>
        )}
        <button className="primary" onClick={run} disabled={running}>
          {running ? "Optimizing…" : `▶ Run ${opt.id}`}
        </button>
        {r?.ok && r.apply && onApply && (
          <button
            disabled={applied}
            onClick={() => {
              onApply(r.apply!);
              setApplied(true);
            }}
          >
            {applied ? "Optimum applied" : "Use optimum"}
          </button>
        )}
        {r?.history && (
          <span className="muted small">
            {r.method}: {r.iterations} iterations, {r.evaluations} flowsheet solves
            {r.failed_evaluations ? `, ${r.failed_evaluations} failed` : ""}
          </span>
        )}
        {r?.ok && !r.converged && <span className="warn-text">Not converged: {r.message}</span>}
      </div>
      {r && !r.ok && r.error && (
        <div className="banner error">
          <b>{r.error.where ? `${r.error.where}: ` : ""}</b>
          {r.error.message}
          {r.error.hint && <div className="hint">💡 {r.error.hint}</div>}
        </div>
      )}
      {r?.history && r.variables && (
        <div className="sens-body">
          <HistoryChart history={r.history} />
          <div className="table-scroll sens-table">
            <table>
              <thead>
                <tr>
                  <th>{r.ok ? (r.converged ? "Optimum" : "Best point") : "Least violation"}</th>
                  <th>Start</th>
                  <th>Result</th>
                  <th>Bounds</th>
                </tr>
              </thead>
              <tbody>
                <tr>
                  <th>
                    {r.sense} {r.objective?.expression}
                  </th>
                  <td>{r.objective?.start != null ? money(r.objective.start) : "–"}</td>
                  <td>
                    <b>{r.objective ? money(r.objective.value) : "–"}</b>
                  </td>
                  <td />
                </tr>
                {r.variables.map((v) => (
                  <tr key={v.label}>
                    <th>{v.label}</th>
                    <td>{formatQuantity(v.start, v.quantity)}</td>
                    <td>
                      <b>{formatQuantity(v.value, v.quantity)}</b>
                      {v.at_bound && <span className="muted small"> at bound</span>}
                    </td>
                    <td className="muted">
                      {formatQuantity(v.lower, v.quantity)} … {formatQuantity(v.upper, v.quantity)}
                    </td>
                  </tr>
                ))}
                {(r.constraints ?? []).map((c, k) => (
                  <tr key={`c${k}`}>
                    <th>
                      {c.expression} {c.relation === "<=" ? "≤" : c.relation === ">=" ? "≥" : "="} {num(c.value)}
                    </th>
                    <td />
                    <td>{num(c.achieved)}</td>
                    <td className={c.satisfied ? "muted" : "warn-text"}>{!c.satisfied ? "violated" : c.active ? "active" : "met"}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}
      {r?.ok && r.result?.costing && <CostingView costing={r.result.costing} title="Costing at the optimum" />}
    </div>
  );
}

export function CostingView({ costing, title }: { costing: CostingResult; title?: string }) {
  const units = Object.entries(costing.units);
  return (
    <div className="costing">
      <p className="small">
        {title && <b>{title}: </b>}
        Capital <b>{money(costing.capital_cost)} $</b> (bare module {money(costing.bare_module_cost)} $ + 18 %), operating{" "}
        <b>{money(costing.operating_cost)} $/yr</b>, TAC <b>{money(costing.tac)} $/yr</b>{" "}
        <span className="muted">
          (annualization {costing.annualization_factor.toPrecision(3)}/yr, CEPCI {costing.cepci}, {costing.hours_per_year} h/yr)
        </span>
      </p>
      <div className="table-scroll">
        <table>
          <thead>
            <tr>
              <th>Unit</th>
              <th>Equipment</th>
              <th>Size</th>
              <th>Purchased [$]</th>
              <th>Bare module [$]</th>
              <th>Utilities</th>
              <th>Operating [$/yr]</th>
            </tr>
          </thead>
          <tbody>
            {units.map(([id, u]) => (
              <tr key={id} title={u.notes.join("\n") || undefined}>
                <th>
                  {id}
                  {u.notes.length > 0 && <span className="muted"> *</span>}
                </th>
                <td>{u.items.map((i) => i.item).join(", ") || "–"}</td>
                <td>{u.items.map((i) => `${Number(i.size.toPrecision(3))} ${i.size_unit}${i.parallel > 1 ? ` ×${i.parallel}` : ""}`).join(", ") || "–"}</td>
                <td>{money(u.purchased_cost)}</td>
                <td>{money(u.bare_module_cost)}</td>
                <td>{u.utilities.map((x) => `${x.utility} ${formatQuantity(Math.abs(x.duty), "power")}`).join(", ") || "–"}</td>
                <td>{money(u.operating_cost)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      {costing.not_costed.length > 0 && <p className="muted small">Not costed: {costing.not_costed.join(", ")}</p>}
    </div>
  );
}
