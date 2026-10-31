// Symbols for the reactors. The four balance-based reactors share a vessel with a sight
// window that shows the model; the CSTR is a stirred tank, the PFR a packed tube.
import { blue, liquid, metal, metalH, outline, pipe, type UnitSymbol } from "./style";

const mark = { fill: "none", stroke: "var(--sym-mark)", strokeWidth: 1.8, strokeLinecap: "round" as const, strokeLinejoin: "round" as const };

/** Vertical reactor on legs, with the model's mark in a round window. */
function vessel(sign: JSX.Element): UnitSymbol {
  return {
    w: 60,
    h: 92,
    body: (
      <>
        <path {...pipe} d="M1 24H12M48 66H59" />
        <path {...pipe} d="M17 76L15 90M43 76L45 90" />
        <path {...metal} d="M12 20A18 9 0 0 1 48 20V70A18 9 0 0 1 12 70Z" />
        <circle cx="30" cy="45" r="11" fill="var(--sym-window)" stroke="var(--sym-line)" strokeWidth={1.2} />
        {sign}
      </>
    ),
    ports: { in: { x: 1, y: 24 }, out: { x: 59, y: 66 } },
  };
}

const CSTR = "M12 26A20 9 0 0 1 52 26V72A20 9 0 0 1 12 72Z";

export const symbols: Record<string, UnitSymbol> = {
  StoichiometricReactor: vessel(<path {...mark} d="M23 45h14M33 41l4 4-4 4" />),
  YieldReactor: vessel(<path {...mark} d="M24 51v-5M30 51v-11M36 51v-8" />),
  EquilibriumReactor: vessel(<path {...mark} d="M23 42h14l-3-3M37 48H23l3 3" />),
  GibbsReactor: vessel(
    <>
      <path {...mark} d="M23 38Q30 62 37 38" />
      <circle cx="30" cy="49.5" r="1.8" fill="var(--sym-mark)" />
    </>,
  ),
  CSTR: {
    w: 64,
    h: 94,
    defs: (
      <clipPath id="fx-clip-cstr">
        <path d={CSTR} />
      </clipPath>
    ),
    body: (
      <>
        <path {...pipe} d="M1 30H12M52 74H63" />
        <path {...pipe} d="M17 78L15 92M47 78L49 92" />
        <path {...metal} d={CSTR} />
        <rect {...liquid} x="0" y="38" width="64" height="50" clipPath="url(#fx-clip-cstr)" />
        <path {...outline} d="M12 38H52" />
        <path {...outline} d={CSTR} />
        <path {...outline} strokeWidth={1.6} d="M32 11V64" />
        <path {...outline} strokeWidth={3} strokeLinecap="round" d="M22 64h20M22 50h20" />
        <rect {...metal} x="28" y="11" width="8" height="8" />
        <rect {...blue} x="25" y="2" width="14" height="10" rx="2" />
      </>
    ),
    ports: { in: { x: 1, y: 30 }, out: { x: 63, y: 74 } },
  },
  PFR: {
    w: 120,
    h: 44,
    body: (
      <>
        <path {...pipe} d="M1 22H8M112 22H119" />
        <rect {...metalH} x="10" y="10" width="100" height="24" rx="3" />
        <rect x="16" y="13" width="88" height="18" fill="url(#fx-catalyst)" />
        <rect {...metalH} x="7" y="6" width="5" height="32" rx="1" />
        <rect {...metalH} x="108" y="6" width="5" height="32" rx="1" />
      </>
    ),
    ports: { in: { x: 1, y: 22 }, out: { x: 119, y: 22 } },
  },
};
