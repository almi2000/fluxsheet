// Symbols for the logic blocks: multiplier, duplicator and selector. Drawn as plain boxes
// (they are not equipment), with a sign of what they do.
import { metalH, outline, pipe, thin, type UnitSymbol } from "./style";

const BOX = "M6 6H42V42H6Z";

export const symbols: Record<string, UnitSymbol> = {
  Mult: {
    w: 48,
    h: 48,
    body: (
      <>
        <path {...pipe} d="M1 24H6M42 24H47" />
        <path {...metalH} strokeDasharray="4 2" d={BOX} />
        <path {...outline} strokeWidth={2} strokeLinecap="round" d="M15 18l8 8M23 18l-8 8" />
        <path {...outline} strokeWidth={1.6} strokeLinecap="round" d="M28 31V17M28 25l6-6M30 23l5 8" />
      </>
    ),
    ports: { in: { x: 1, y: 24 }, out: { x: 47, y: 24 } },
  },
  Dupl: {
    w: 48,
    h: 56,
    body: (
      <>
        <path {...pipe} d="M1 28H6M42 10H47M42 22H47M42 34H47M42 46H47" />
        <path {...metalH} strokeDasharray="4 2" d="M6 6H42V50H6Z" />
        <path {...outline} strokeWidth={1.6} strokeLinecap="round" d="M12 28H22M22 28L36 10M22 28L36 22M22 28L36 34M22 28L36 46" />
        <circle {...thin} cx="22" cy="28" r="2" />
      </>
    ),
    ports: { in: { x: 1, y: 28 }, out1: { x: 47, y: 10 }, out2: { x: 47, y: 22 }, out3: { x: 47, y: 34 }, out4: { x: 47, y: 46 } },
  },
  Selector: {
    w: 48,
    h: 48,
    body: (
      <>
        <path {...pipe} d="M1 24H6M42 24H47" />
        <path {...metalH} strokeDasharray="4 2" d={BOX} />
        <path {...thin} d="M10 14H16M10 24H16M10 34H16" />
        <path {...outline} strokeWidth={1.8} strokeLinecap="round" d="M16 14L32 24H38" />
        <circle {...outline} cx="16" cy="14" r="1.6" />
        <circle {...outline} cx="16" cy="24" r="1.6" />
        <circle {...outline} cx="16" cy="34" r="1.6" />
      </>
    ),
    ports: { in: { x: 1, y: 24 }, out: { x: 47, y: 24 } },
  },
};
