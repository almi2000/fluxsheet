// Symbol for the decanter: a horizontal drum with a light and a heavy liquid layer, a weir
// and a boot for the heavy phase.
import { liquid, metalH, oil, outline, pipe, type UnitSymbol } from "./style";

const DRUM = "M14 8H86A10 19 0 0 1 86 46H14A10 19 0 0 1 14 8Z";

export const symbols: Record<string, UnitSymbol> = {
  Decanter: {
    w: 100,
    h: 54,
    defs: (
      <clipPath id="fx-clip-decanter">
        <path d={DRUM} />
      </clipPath>
    ),
    body: (
      <>
        <path {...pipe} d="M0 27H4M80 8V3H99M95 24H99M54 50H99" />
        <rect {...liquid} x="40" y="42" width="14" height="10" stroke="var(--sym-line)" strokeWidth={1.2} />
        <path {...metalH} d={DRUM} />
        <g clipPath="url(#fx-clip-decanter)">
          <rect {...oil} x="0" y="22" width="100" height="12" />
          <rect {...liquid} x="0" y="34" width="100" height="14" />
        </g>
        <path {...outline} d="M8 22H92M5 34H95M74 46V24" />
        <path {...outline} d={DRUM} />
      </>
    ),
    ports: { in: { x: 0, y: 27 }, vapor: { x: 99, y: 3 }, light: { x: 99, y: 24 }, heavy: { x: 99, y: 50 } },
  },
};
