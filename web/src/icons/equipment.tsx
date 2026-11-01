// Symbols for the mixing tank, thin-film evaporator, turbine, multistage compressor, pipe and
// extraction column.
import { blue, cold, liquid, metal, metalH, metalR, oil, outline, pipe, thin, type UnitSymbol } from "./style";

const TANK = "M12 16H58V68A23 8 0 0 1 12 68Z";
const EVAPORATOR = "M18 34H46V92H18Z";
const EXTRACTOR = "M18 18A14 7 0 0 1 46 18V132A14 7 0 0 1 18 132Z";

export const symbols: Record<string, UnitSymbol> = {
  MixingTank: {
    w: 70,
    h: 84,
    defs: (
      <clipPath id="fx-clip-tank">
        <path d={TANK} />
      </clipPath>
    ),
    body: (
      <>
        <path {...pipe} d="M1 26H12M50 16V6H69M58 66H69" />
        <path {...metal} d={TANK} />
        <rect {...liquid} x="0" y="32" width="70" height="52" clipPath="url(#fx-clip-tank)" />
        <path {...outline} d="M12 32H58" />
        <path {...outline} d={TANK} />
        <path {...outline} strokeWidth={1.6} d="M35 9V62" />
        <path {...outline} strokeWidth={2.4} strokeLinecap="round" d="M25 64l10-4 10 4M27 46l8-3 8 3" />
        <rect {...blue} x="28" y="1" width="14" height="9" rx="2" />
      </>
    ),
    ports: { in: { x: 1, y: 26 }, vapor: { x: 69, y: 6 }, liquid: { x: 69, y: 66 } },
    jacket: TANK,
  },
  ThinFilmEvaporator: {
    w: 64,
    h: 128,
    body: (
      <>
        <path {...pipe} d="M1 40H18M46 28H63M32 112V122H63" />
        <path {...metal} d="M18 20H46V100L37 112H27L18 100Z" />
        <rect x="14" y="34" width="36" height="58" rx="3" fill="none" stroke="var(--sym-hot)" strokeWidth={2} opacity={0.75} />
        <path {...outline} strokeWidth={1.4} d="M32 12V100" />
        <path {...thin} d="M22 44l8 6M42 44l-8 6M22 60l8 6M42 60l-8 6M22 76l8 6M42 76l-8 6" />
        <rect {...metalH} x="14" y="16" width="36" height="5" rx="1" />
        <rect {...blue} x="25" y="2" width="14" height="10" rx="2" />
      </>
    ),
    ports: { in: { x: 1, y: 40 }, vapor: { x: 63, y: 28 }, residue: { x: 63, y: 122 } },
    jacket: EVAPORATOR,
  },
  Turbine: {
    w: 60,
    h: 48,
    body: (
      <>
        <path {...pipe} d="M1 24H6M54 24H59" />
        <path {...blue} d="M6 15L54 4V44L6 33Z" />
        <path {...thin} stroke="var(--sym-shine)" d="M18 13v22M30 10v28M42 7v34" />
      </>
    ),
    ports: { in: { x: 1, y: 24 }, out: { x: 59, y: 24 } },
  },
  MultistageCompressor: {
    w: 112,
    h: 58,
    body: (
      <>
        <path {...pipe} d="M1 26H4M30 26H44M66 26H80M106 26H111M55 37V50H111" />
        <path {...blue} d="M4 6L30 14V38L4 46Z" />
        <path {...blue} d="M80 10L106 16V36L80 42Z" />
        <circle {...metalR} cx="55" cy="26" r="11" />
        <path {...cold} d="M44 26h2l2-4 4 8 4-8 4 8 2-4h2" />
      </>
    ),
    ports: { in: { x: 1, y: 26 }, out: { x: 111, y: 26 }, condensate: { x: 111, y: 50 } },
  },
  Pipe: {
    w: 112,
    h: 24,
    body: (
      <>
        <path {...pipe} d="M1 12H5M107 12H111" />
        <rect {...metalH} x="6" y="8" width="100" height="8" rx="2" />
        <rect {...metalH} x="4" y="4" width="4" height="16" rx="1" />
        <rect {...metalH} x="104" y="4" width="4" height="16" rx="1" />
      </>
    ),
    ports: { in: { x: 1, y: 12 }, out: { x: 111, y: 12 } },
  },
  ExtractionColumn: {
    w: 64,
    h: 150,
    defs: (
      <clipPath id="fx-clip-extractor">
        <path d={EXTRACTOR} />
      </clipPath>
    ),
    body: (
      <>
        <path {...pipe} d="M1 24H18M1 126H18M32 11V5H63M32 139V146H63" />
        <path {...metal} d={EXTRACTOR} />
        <g clipPath="url(#fx-clip-extractor)">
          <rect {...oil} x="0" y="0" width="64" height="34" opacity={0.85} />
          <rect {...liquid} x="0" y="116" width="64" height="40" opacity={0.85} />
          <rect x="18" y="38" width="28" height="74" fill="url(#fx-catalyst)" />
        </g>
        <path {...outline} d="M18 50H40M24 66H46M18 82H40M24 98H46" />
        <path {...outline} d={EXTRACTOR} />
      </>
    ),
    ports: {
      feed: { x: 1, y: 24 },
      solvent: { x: 1, y: 126 },
      extract: { x: 63, y: 5 },
      raffinate: { x: 63, y: 146 },
    },
  },
};
