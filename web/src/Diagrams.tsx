import { useEffect, useState } from "react";
import { toDisplay } from "./quantities";
import "./diagrams.css";

// Response shapes of /api/diagrams/*; kept here so the shared types stay untouched.
interface DiagramError {
  message: string;
  where: string | null;
  hint: string | null;
}

interface LiquidGap {
  three_phase: number; // T (T-xy) or P (P-xy) of the three-phase line
  curve?: { T: number; x_light: number; x_heavy: number }[]; // T-xy: solubility curve
  x_light?: number; // P-xy: the two liquids at the diagram temperature
  x_heavy?: number;
}

interface BinaryDiagram {
  kind: "txy" | "pxy";
  components: string[];
  method: string;
  x: number[];
  y: number[];
  T: number | number[]; // T-xy: bubble temperatures; P-xy: the fixed temperature
  P: number | number[]; // P-xy: bubble pressures; T-xy: the fixed pressure
  azeotrope: { x: number; T?: number; P?: number; type: string } | null;
  liquid_liquid: LiquidGap | null;
  warnings: string[];
  parameter_sources: { pair: string; source: string }[];
}

interface Envelope {
  components: string[];
  method: string;
  bubble: { T: number; P: number }[];
  dew: { T: number; P: number }[];
  cricondenbar: { T: number; P: number } | null;
  cricondentherm: { T: number; P: number } | null;
  warnings: string[];
}

type Kind = "txy" | "pxy" | "envelope";
type Response<T> = ({ ok: true } & T) | { ok: false; error: DiagramError };

async function post<T>(path: string, body: unknown): Promise<Response<T>> {
  const res = await fetch(`/api/diagrams/${path}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!res.ok) throw new Error(`${res.status} ${res.statusText}`);
  return res.json();
}

const EOS = ["peng-robinson", "srk"];

/** Phase diagrams for the flowsheet's components: T-xy, P-xy and phase envelopes. */
export function Diagrams({
  components,
  method,
  fractions,
}: {
  components: string[];
  method: string | null;
  fractions: Record<string, number> | null;
}) {
  const [kind, setKind] = useState<Kind>(EOS.includes(method ?? "") && components.length > 2 ? "envelope" : "txy");
  const [pair, setPair] = useState<[string, string]>([components[0] ?? "", components[1] ?? ""]);
  const [chosen, setChosen] = useState<string>(method ?? "auto");
  const [pressure, setPressure] = useState("1 atm");
  const [temperature, setTemperature] = useState("80 degC");
  const [shares, setShares] = useState<Record<string, number>>(() => defaultShares(components, fractions));
  const [methods, setMethods] = useState<{ key: string; label: string }[]>([]);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<DiagramError | null>(null);
  const [binary, setBinary] = useState<BinaryDiagram | null>(null);
  const [envelope, setEnvelope] = useState<Envelope | null>(null);

  useEffect(() => {
    fetch("/api/methods")
      .then((r) => r.json())
      .then(setMethods)
      .catch(() => setMethods([]));
  }, []);

  // A new run may bring other components: keep the selection valid.
  useEffect(() => {
    setPair(([a, b]) =>
      components.includes(a) && components.includes(b) && a !== b ? [a, b] : [components[0] ?? "", components[1] ?? ""],
    );
    setShares((s) => (components.every((c) => c in s) ? s : defaultShares(components, fractions)));
  }, [components.join("|")]); // eslint-disable-line react-hooks/exhaustive-deps

  const envelopeMethod = EOS.includes(chosen) ? chosen : "peng-robinson";

  async function plot() {
    setBusy(true);
    setError(null);
    try {
      if (kind === "envelope") {
        const r = await post<Envelope>("envelope", {
          components,
          fractions: components.map((c) => shares[c] ?? 0),
          method: envelopeMethod,
        });
        if (r.ok) setEnvelope(r);
        else setError(r.error);
      } else {
        const r = await post<BinaryDiagram>("binary", {
          components: pair,
          method: chosen,
          ...(kind === "txy" ? { P: pressure } : { T: temperature }),
        });
        if (r.ok) setBinary(r);
        else setError(r.error);
      }
    } catch (e) {
      setError({ message: `could not reach the server (${(e as Error).message})`, where: null, hint: "is `fluxsheet serve` running?" });
    } finally {
      setBusy(false);
    }
  }

  // Plot once when the view opens, so there is something to look at right away.
  useEffect(() => {
    if (components.length >= 2) void plot();
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  const shown = kind === "envelope" ? envelope : binary && binary.kind === kind ? binary : null;

  return (
    <div className="diagrams">
      <form
        className="diagram-controls"
        onSubmit={(e) => {
          e.preventDefault();
          void plot();
        }}
      >
        <div className="segmented">
          {(
            [
              ["txy", "T-xy"],
              ["pxy", "P-xy"],
              ["envelope", "Envelope"],
            ] as [Kind, string][]
          ).map(([k, label]) => (
            <button type="button" key={k} className={kind === k ? "on" : ""} onClick={() => setKind(k)}>
              {label}
            </button>
          ))}
        </div>
        {kind !== "envelope" ? (
          <>
            <label className="diagram-field">
              <span>Components (x = first)</span>
              <div className="diagram-pair">
                {[0, 1].map((i) => (
                  <select
                    key={i}
                    value={pair[i]}
                    onChange={(e) => setPair((p) => (i === 0 ? [e.target.value, p[1]] : [p[0], e.target.value]))}
                  >
                    {components.map((c) => (
                      <option key={c} value={c}>
                        {c}
                      </option>
                    ))}
                  </select>
                ))}
              </div>
            </label>
            {kind === "txy" ? (
              <label className="diagram-field">
                <span>Pressure</span>
                <input value={pressure} onChange={(e) => setPressure(e.target.value)} placeholder="1 atm" />
              </label>
            ) : (
              <label className="diagram-field">
                <span>Temperature</span>
                <input value={temperature} onChange={(e) => setTemperature(e.target.value)} placeholder="80 degC" />
              </label>
            )}
          </>
        ) : (
          <div className="diagram-field">
            <span>Mole fractions</span>
            {components.map((c) => (
              <label key={c} className="diagram-share">
                <span>{c}</span>
                <input
                  type="number"
                  min={0}
                  step="any"
                  value={shares[c] ?? 0}
                  onChange={(e) => setShares((s) => ({ ...s, [c]: Number(e.target.value) }))}
                />
              </label>
            ))}
          </div>
        )}
        <label className="diagram-field">
          <span>Property method</span>
          <select value={kind === "envelope" ? envelopeMethod : chosen} onChange={(e) => setChosen(e.target.value)}>
            {kind !== "envelope" && <option value="auto">Auto</option>}
            {methods
              .filter((m) => kind !== "envelope" || EOS.includes(m.key))
              .map((m) => (
                <option key={m.key} value={m.key}>
                  {m.label}
                </option>
              ))}
          </select>
        </label>
        <button className="primary" type="submit" disabled={busy || components.length < 2}>
          {busy ? "Calculating…" : "Plot"}
        </button>
      </form>

      <div className="diagram-plot">
        {components.length < 2 && <div className="muted">Phase diagrams need at least two components.</div>}
        {error && (
          <div className="callout warn">
            {error.message}
            {error.hint && <div className="muted small">{error.hint}</div>}
          </div>
        )}
        {!error && shown && "x" in shown && <BinaryChart d={shown} />}
        {!error && shown && "bubble" in shown && <EnvelopeChart d={shown} />}
        {!error && !shown && busy && <div className="muted">Calculating…</div>}
      </div>
    </div>
  );
}

function defaultShares(components: string[], fractions: Record<string, number> | null): Record<string, number> {
  return Object.fromEntries(
    components.map((c) => [c, Number((fractions?.[c] ?? 1 / components.length).toFixed(4))]),
  );
}

// ---- charts ------------------------------------------------------------------------------

const W = 520;
const H = 210;
const PAD = { l: 46, r: 12, t: 10, b: 32 };

function ticks(lo: number, hi: number): number[] {
  if (lo === hi) return [lo];
  const step = Math.pow(10, Math.floor(Math.log10((hi - lo) / 3)));
  const nice = [1, 2, 5, 10].map((m) => m * step).find((s) => (hi - lo) / s <= 6) ?? step * 10;
  const out: number[] = [];
  for (let v = Math.ceil(lo / nice) * nice; v <= hi + 1e-9; v += nice) out.push(v);
  return out;
}

function padded(values: number[]): [number, number] {
  let [lo, hi] = [Math.min(...values), Math.max(...values)];
  if (lo === hi) [lo, hi] = [lo - 1, hi + 1];
  const m = (hi - lo) * 0.06;
  return [lo - m, hi + m];
}

interface Series {
  points: [number, number][];
  className: string;
  label: string;
}

interface Marker {
  at: [number, number];
  label: string;
}

function Chart({
  xLabel,
  yLabel,
  x,
  y,
  series,
  markers,
}: {
  xLabel: string;
  yLabel: string;
  x: [number, number];
  y: [number, number];
  series: Series[];
  markers: Marker[];
}) {
  const sx = (v: number) => PAD.l + ((v - x[0]) / (x[1] - x[0])) * (W - PAD.l - PAD.r);
  const sy = (v: number) => H - PAD.b - ((v - y[0]) / (y[1] - y[0])) * (H - PAD.t - PAD.b);
  const path = (pts: [number, number][]) =>
    pts.map(([a, b], i) => `${i ? "L" : "M"}${sx(a).toFixed(1)},${sy(b).toFixed(1)}`).join("");
  const legend = series.filter((s, i) => series.findIndex((t) => t.label === s.label) === i);
  return (
    <>
      <svg viewBox={`0 0 ${W} ${H}`} className="diagram-chart" role="img" aria-label={`${yLabel} vs ${xLabel}`}>
        {ticks(y[0], y[1]).map((t) => (
          <g key={`y${t}`}>
            <line x1={PAD.l} x2={W - PAD.r} y1={sy(t)} y2={sy(t)} className="grid" />
            <text x={PAD.l - 4} y={sy(t) + 3} textAnchor="end">
              {Number(t.toPrecision(4))}
            </text>
          </g>
        ))}
        {ticks(x[0], x[1]).map((t) => (
          <g key={`x${t}`}>
            <line x1={sx(t)} x2={sx(t)} y1={PAD.t} y2={H - PAD.b} className="grid" />
            <text x={sx(t)} y={H - PAD.b + 13} textAnchor="middle">
              {Number(t.toPrecision(4))}
            </text>
          </g>
        ))}
        <text x={(PAD.l + W - PAD.r) / 2} y={H - 3} textAnchor="middle" className="axis-label">
          {xLabel}
        </text>
        <text
          transform={`translate(11 ${(PAD.t + H - PAD.b) / 2}) rotate(-90)`}
          textAnchor="middle"
          className="axis-label"
        >
          {yLabel}
        </text>
        {series.map((s, i) => (
          <path key={i} d={path(s.points)} className={`curve ${s.className}`} />
        ))}
        {markers.map((m) => (
          <circle key={m.label} cx={sx(m.at[0])} cy={sy(m.at[1])} r={3.5} className="marker">
            <title>{m.label}</title>
          </circle>
        ))}
      </svg>
      <div className="diagram-legend small">
        {legend.map((s) => (
          <span key={s.label}>
            <i className={`swatch ${s.className}`} />
            {s.label}
          </span>
        ))}
        {markers.map((m) => (
          <span key={m.label}>
            <i className="swatch marker" />
            {m.label}
          </span>
        ))}
      </div>
    </>
  );
}

const degC = (k: number) => toDisplay(k, "temperature");
const bar = (pa: number) => toDisplay(pa, "pressure");

function BinaryChart({ d }: { d: BinaryDiagram }) {
  const txy = d.kind === "txy";
  const values = (txy ? d.T : d.P) as number[];
  const show = txy ? degC : bar;
  const v = values.map(show);
  const bubble: [number, number][] = d.x.map((x, i) => [x, v[i]]);
  const dew: [number, number][] = d.y.map((y, i): [number, number] => [y, v[i]]).sort((a, b) => a[0] - b[0]);
  const series: Series[] = [
    { points: bubble, className: "liquid", label: "Bubble line (liquid x)" },
    { points: dew, className: "vapor", label: "Dew line (vapor y)" },
  ];
  const all = [...v];
  const gap = d.liquid_liquid;
  if (gap) {
    const top = show(gap.three_phase);
    if (gap.curve && gap.curve.length) {
      const curve = gap.curve.map((p) => ({ ...p, t: degC(p.T) }));
      const light = curve.map((p): [number, number] => [p.x_light, p.t]);
      const heavy = curve.map((p): [number, number] => [p.x_heavy, p.t]);
      const lo = Math.min(curve[0].x_light, curve[0].x_heavy);
      const hi = Math.max(curve[0].x_light, curve[0].x_heavy);
      series.push(
        { points: light, className: "lle", label: "Two liquids" },
        { points: heavy, className: "lle", label: "Two liquids" },
        { points: [[lo, top], [hi, top]], className: "lle", label: "Two liquids" },
      );
      all.push(...curve.map((p) => p.t));
    } else if (gap.x_light !== undefined && gap.x_heavy !== undefined) {
      const ceiling = Math.max(...v) * 1.15;
      series.push(
        { points: [[gap.x_light, top], [gap.x_light, ceiling]], className: "lle", label: "Two liquids" },
        { points: [[gap.x_heavy, top], [gap.x_heavy, ceiling]], className: "lle", label: "Two liquids" },
        { points: [[gap.x_light, top], [gap.x_heavy, top]], className: "lle", label: "Two liquids" },
      );
      all.push(ceiling);
    }
  }
  const az = d.azeotrope;
  const markers: Marker[] = az
    ? [
        {
          at: [az.x, show((txy ? az.T : az.P) ?? 0)],
          label: `${az.type === "heterogeneous" ? "Heteroazeotrope" : "Azeotrope"}: x = ${az.x.toFixed(3)}, ${
            txy ? `${degC(az.T ?? 0).toFixed(2)} °C` : `${bar(az.P ?? 0).toFixed(4)} bar`
          }`,
        },
      ]
    : [];
  const fixed = txy ? `${bar(d.P as number).toPrecision(4)} bar` : `${degC(d.T as number).toFixed(1)} °C`;
  return (
    <section className="profile">
      <h3>
        {txy ? "T-xy" : "P-xy"} · {d.components[0]} / {d.components[1]} · {d.method} · {fixed}
      </h3>
      <Chart
        xLabel={`mole fraction ${d.components[0]}`}
        yLabel={txy ? "T [°C]" : "P [bar]"}
        x={[0, 1]}
        y={padded(all)}
        series={series}
        markers={markers}
      />
      <Provenance d={d} />
    </section>
  );
}

function Provenance({ d }: { d: BinaryDiagram }) {
  if (!d.parameter_sources.length && !d.warnings.length) return null;
  return (
    <div className="muted small">
      {d.parameter_sources.map((s) => (
        <span key={s.pair}>
          Binary parameters {s.pair}: <b>{s.source}</b>.{" "}
        </span>
      ))}
      {d.warnings.map((w) => (
        <div key={w} className="warn-text">
          ⚠ {w}
        </div>
      ))}
    </div>
  );
}

function EnvelopeChart({ d }: { d: Envelope }) {
  const pts = (c: { T: number; P: number }[]) => c.map((p): [number, number] => [degC(p.T), bar(p.P)]);
  const bubble = pts(d.bubble);
  const dew = pts(d.dew);
  const all = [...bubble, ...dew];
  if (!all.length) return <div className="muted">No saturation points converged for this mixture.</div>;
  const markers: Marker[] = [];
  if (d.cricondenbar)
    markers.push({
      at: [degC(d.cricondenbar.T), bar(d.cricondenbar.P)],
      label: `Cricondenbar ${bar(d.cricondenbar.P).toFixed(1)} bar`,
    });
  if (d.cricondentherm)
    markers.push({
      at: [degC(d.cricondentherm.T), bar(d.cricondentherm.P)],
      label: `Cricondentherm ${degC(d.cricondentherm.T).toFixed(1)} °C`,
    });
  return (
    <section className="profile">
      <h3>
        Phase envelope · {d.components.join(", ")} · {d.method}
      </h3>
      <Chart
        xLabel="T [°C]"
        yLabel="P [bar]"
        x={padded(all.map((p) => p[0]))}
        y={[0, padded(all.map((p) => p[1]))[1]]}
        series={[
          { points: bubble, className: "liquid", label: "Bubble line" },
          { points: dew, className: "vapor", label: "Dew line" },
        ]}
        markers={markers}
      />
      <div className="muted small">Traced towards the critical region until the flash stops converging.</div>
    </section>
  );
}
