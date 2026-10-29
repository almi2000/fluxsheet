import { useState } from "react";
import { displayUnit, formatQuantity, toDisplay } from "./quantities";
import type { Profile } from "./types";

const W = 280;
const H = 150;
const PAD = { l: 44, r: 10, t: 10, b: 26 };

export function ticks(lo: number, hi: number): number[] {
  if (lo === hi) return [lo];
  const step = Math.pow(10, Math.floor(Math.log10((hi - lo) / 3)));
  const nice = [1, 2, 5, 10].map((m) => m * step).find((s) => (hi - lo) / s <= 5) ?? step * 10;
  const out: number[] = [];
  for (let v = Math.ceil(lo / nice) * nice; v <= hi + 1e-9; v += nice) out.push(v);
  return out;
}

/** A profile (e.g. temperature per stage): one series at a time as a line chart plus table. */
export function ProfileView({ profile }: { profile: Profile }) {
  const [active, setActive] = useState(0);
  const [showTable, setShowTable] = useState(false);
  const series = profile.series[Math.min(active, profile.series.length - 1)];
  if (!series) return null;

  const ys = series.values.map((v) => toDisplay(v, series.quantity));
  const xs = profile.x;
  const [x0, x1] = [Math.min(...xs), Math.max(...xs)];
  let [y0, y1] = [Math.min(...ys), Math.max(...ys)];
  if (y0 === y1) [y0, y1] = [y0 - 1, y1 + 1];
  const sx = (x: number) => PAD.l + ((x - x0) / (x1 - x0 || 1)) * (W - PAD.l - PAD.r);
  const sy = (y: number) => H - PAD.b - ((y - y0) / (y1 - y0)) * (H - PAD.t - PAD.b);
  const path = xs.map((x, i) => `${i ? "L" : "M"}${sx(x).toFixed(1)},${sy(ys[i]).toFixed(1)}`).join("");
  const unit = displayUnit(series.quantity);

  return (
    <section className="profile">
      <h3>{profile.label}</h3>
      {profile.series.length > 1 && (
        <select value={active} onChange={(e) => setActive(Number(e.target.value))}>
          {profile.series.map((s, i) => (
            <option key={s.label} value={i}>
              {s.label}
            </option>
          ))}
        </select>
      )}
      <svg viewBox={`0 0 ${W} ${H}`} className="profile-chart" role="img" aria-label={`${series.label} vs ${profile.x_label}`}>
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
        <text x={(PAD.l + W - PAD.r) / 2} y={H - 2} textAnchor="middle" className="axis-label">
          {profile.x_label}
        </text>
        <path d={path} className="line" />
        {xs.map((x, i) => (
          <circle key={i} cx={sx(x)} cy={sy(ys[i])} r={2.2} className="dot">
            <title>{`${profile.x_label} ${x}: ${formatQuantity(series.values[i], series.quantity)}`}</title>
          </circle>
        ))}
      </svg>
      <div className="muted small">
        {series.label}
        {unit ? ` [${unit}]` : ""}
        {" · "}
        <button className="link" onClick={() => setShowTable((v) => !v)}>
          {showTable ? "hide table" : "show table"}
        </button>
      </div>
      {showTable && (
        <table className="mini">
          <thead>
            <tr>
              <th>{profile.x_label}</th>
              {profile.series.map((s) => (
                <th key={s.label}>{s.label}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {xs.map((x, i) => (
              <tr key={i}>
                <td>{x}</td>
                {profile.series.map((s) => (
                  <td key={s.label}>{formatQuantity(s.values[i], s.quantity)}</td>
                ))}
              </tr>
            ))}
          </tbody>
        </table>
      )}
    </section>
  );
}
