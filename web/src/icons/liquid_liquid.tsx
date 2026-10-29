// Symbols for liquid-liquid separators. Simplified P&ID style on a 40x40 grid.
import { S } from "./style";

export const shapes: Record<string, JSX.Element> = {
  Decanter: (
    <>
      <rect {...S} x="3" y="11" width="34" height="18" rx="9" />
      <path {...S} d="M8 21h24" strokeDasharray="3 3" />
      <path {...S} d="M27 15v14" strokeWidth={1.4} />
    </>
  ),
};
