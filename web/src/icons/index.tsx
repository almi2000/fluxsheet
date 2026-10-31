// Equipment symbols. Every file in this folder that exports `symbols` (unit type -> drawing
// with port positions) is picked up automatically, so each unit family keeps its own file.
import type { PortAt, UnitSymbol } from "./style";

export type { PortAt, UnitSymbol };

const modules = import.meta.glob<{ symbols?: Record<string, UnitSymbol> }>("./*.tsx", { eager: true });
const symbols: Record<string, UnitSymbol> = Object.assign({}, ...Object.values(modules).map((m) => m.symbols ?? {}));

const fallback: UnitSymbol = {
  w: 60,
  h: 50,
  body: <rect x="4" y="4" width="52" height="42" rx="5" fill="url(#fx-metal-h)" stroke="var(--sym-line)" strokeWidth={1.2} />,
  ports: {},
};

export function getSymbol(type: string): UnitSymbol {
  return symbols[type] ?? fallback;
}

/** Where port `name` sits; ports a symbol does not place are spread along its left or right side. */
export function portAt(sym: UnitSymbol, name: string, side: "in" | "out", index: number, count: number): PortAt {
  return sym.ports[name] ?? { x: side === "in" ? 0 : sym.w, y: ((index + 1) / (count + 1)) * sym.h };
}

/** The symbol drawn at its natural size (canvas) or fitted into a square (palette, inspector). */
export function UnitIcon({ type, size }: { type: string; size?: number }) {
  const sym = getSymbol(type);
  return (
    <svg
      width={size ?? sym.w}
      height={size ?? sym.h}
      viewBox={`0 0 ${sym.w} ${sym.h}`}
      overflow="visible"
      aria-hidden="true"
      className="unit-icon"
    >
      {sym.body}
    </svg>
  );
}

const stop = (offset: number, color: string) => <stop offset={offset} style={{ stopColor: `var(${color})` }} />;

/** Gradients, patterns and clip paths the symbols refer to. Mounted once per page. */
export function SymbolDefs() {
  return (
    <svg width="0" height="0" style={{ position: "absolute" }} aria-hidden="true">
      <defs>
        <linearGradient id="fx-metal">
          {stop(0, "--sym-m0")}
          {stop(0.35, "--sym-m2")}
          {stop(0.7, "--sym-m1")}
          {stop(1, "--sym-m0")}
        </linearGradient>
        <linearGradient id="fx-metal-h" href="#fx-metal" x2="0" y2="1" />
        <radialGradient id="fx-metal-r" cx="0.38" cy="0.35" r="0.7">
          {stop(0, "--sym-m2")}
          {stop(0.6, "--sym-m1")}
          {stop(1, "--sym-m0")}
        </radialGradient>
        <linearGradient id="fx-blue">
          {stop(0, "--sym-b0")}
          {stop(0.35, "--sym-b2")}
          {stop(0.7, "--sym-b1")}
          {stop(1, "--sym-b0")}
        </linearGradient>
        <linearGradient id="fx-blue-h" href="#fx-blue" x2="0" y2="1" />
        <linearGradient id="fx-liquid" x2="0" y2="1">
          {stop(0, "--sym-liq1")}
          {stop(1, "--sym-liq0")}
        </linearGradient>
        <linearGradient id="fx-oil" x2="0" y2="1">
          {stop(0, "--sym-oil1")}
          {stop(1, "--sym-oil0")}
        </linearGradient>
        <pattern id="fx-mesh" width="4" height="4" patternUnits="userSpaceOnUse">
          <path d="M0 4L4 0" stroke="var(--sym-line)" strokeWidth={0.7} />
        </pattern>
        <pattern id="fx-packing" width="7" height="7" patternUnits="userSpaceOnUse">
          <circle cx="3.5" cy="3.5" r="2" fill="none" stroke="var(--sym-line)" strokeWidth={0.7} opacity={0.55} />
        </pattern>
        <pattern id="fx-catalyst" width="6" height="6" patternUnits="userSpaceOnUse">
          <circle cx="3" cy="3" r="1.4" fill="var(--sym-line)" opacity={0.45} />
        </pattern>
        {Object.entries(symbols).map(([type, s]) => s.defs && <g key={type}>{s.defs}</g>)}
      </defs>
    </svg>
  );
}
