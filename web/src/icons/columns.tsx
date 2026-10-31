// Symbols for columns: shortcut column, distillation column with condenser, reflux drum and
// reboiler, and a packed absorber.
import { cold, hot, liquid, metal, metalH, metalR, outline, pipe, thin, type UnitSymbol } from "./style";

const range = (from: number, to: number, step: number) =>
  Array.from({ length: Math.floor((to - from) / step) + 1 }, (_, i) => from + i * step);

/** Sieve trays with alternating downcomers between x0 and x1. */
const trays = (x0: number, x1: number, ys: number[]) =>
  ys.map((y, i) => (i % 2 ? `M${x0 + 6} ${y}H${x1}` : `M${x0} ${y}H${x1 - 6}`)).join("");

const REFLUX_DRUM = "M68 36H88A5 7 0 0 1 88 50H68A5 7 0 0 1 68 36Z";

export const symbols: Record<string, UnitSymbol> = {
  ShortcutColumn: {
    w: 60,
    h: 146,
    body: (
      <>
        <path {...pipe} d="M1 72H14M26 10V5H59M38 126H59" />
        <path {...metal} d="M18 132h16v12H18z" />
        <path {...metal} d="M14 16A12 6 0 0 1 38 16V128A12 6 0 0 1 14 128Z" />
        <path {...outline} strokeDasharray="3 2.5" d={range(24, 120, 12).map((y) => `M14 ${y}H38`).join("")} />
      </>
    ),
    ports: { in: { x: 1, y: 72 }, distillate: { x: 59, y: 5 }, bottoms: { x: 59, y: 126 } },
  },
  DistillationColumn: {
    w: 110,
    h: 180,
    defs: (
      <clipPath id="fx-clip-reflux">
        <path d={REFLUX_DRUM} />
      </clipPath>
    ),
    body: (
      <>
        <path
          {...pipe}
          d="M1 90H26M40 19V8H74V14M74 30V36M68 50V58H54M86 36V30H109M92 46H109M40 157V164H66M81 160V140H54M96 167H109"
        />
        <path {...metal} d="M26 26A14 7 0 0 1 54 26V150A14 7 0 0 1 26 150Z" />
        <path {...thin} d={trays(26, 54, range(36, 144, 9))} />
        <circle {...metalR} cx="74" cy="22" r="8" />
        <path {...cold} d="M66 22h2l2-4 4 8 4-8 2 4h2" />
        <path {...metalH} d={REFLUX_DRUM} />
        <rect {...liquid} x="60" y="44" width="36" height="8" clipPath="url(#fx-clip-reflux)" />
        <path {...outline} d={REFLUX_DRUM} />
        <path {...metalH} d="M71 160H91A5 7 0 0 1 91 174H71A5 7 0 0 1 71 160Z" />
        <path {...hot} d="M70 168h3l2-3 3 6 3-6 3 6 3-6 2 3h3" />
      </>
    ),
    ports: {
      in: { x: 1, y: 90 },
      distillate: { x: 109, y: 46 },
      vapor_distillate: { x: 109, y: 30 },
      bottoms: { x: 109, y: 167 },
    },
  },
  Absorber: {
    w: 60,
    h: 156,
    body: (
      <>
        <path {...pipe} d="M1 24H16M1 130H16M30 11V5H59M30 143V150H59" />
        <path {...metal} d="M16 18A14 7 0 0 1 44 18V136A14 7 0 0 1 16 136Z" />
        <rect x="16" y="32" width="28" height="90" fill="url(#fx-packing)" />
        <path {...outline} d="M16 32H44M16 122H44" />
        <path {...thin} d="M18 26H42M22 26v3M30 26v3M38 26v3" />
      </>
    ),
    ports: {
      liquid_in: { x: 1, y: 24 },
      gas_in: { x: 1, y: 130 },
      gas_out: { x: 59, y: 5 },
      liquid_out: { x: 59, y: 150 },
    },
  },
};
