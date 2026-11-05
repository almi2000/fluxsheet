// Symbol for the batch reactor: a stirred kettle with a dished bottom, charged from the top,
// emptied through the bottom valve, with a vent nozzle on the head.
import { blue, liquid, metal, outline, pipe, type UnitSymbol } from "./style";

const KETTLE = "M12 24A20 8 0 0 1 52 24V64Q52 84 32 84Q12 84 12 64Z";

export const symbols: Record<string, UnitSymbol> = {
  BatchReactor: {
    w: 64,
    h: 100,
    defs: (
      <clipPath id="fx-clip-batch">
        <path d={KETTLE} />
      </clipPath>
    ),
    body: (
      <>
        <path {...pipe} d="M1 30H12M32 84V92H63M44 17V8H63" />
        <path {...outline} d="M40 88h8v8h-8z" />
        <path {...metal} d={KETTLE} />
        <rect {...liquid} x="0" y="44" width="64" height="44" clipPath="url(#fx-clip-batch)" />
        <path {...outline} d="M12 44H52" />
        <path {...outline} d={KETTLE} />
        <path {...outline} strokeWidth={1.6} d="M32 11V70" />
        <path {...outline} strokeWidth={2.4} strokeLinecap="round" d="M23 70Q32 64 41 70" />
        <rect {...blue} x="25" y="2" width="14" height="10" rx="2" />
      </>
    ),
    ports: { in: { x: 1, y: 30 }, out: { x: 63, y: 92 }, vent: { x: 63, y: 8 } },
    jacket: KETTLE,
  },
};
