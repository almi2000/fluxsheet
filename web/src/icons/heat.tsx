// Symbols for two-sided heat exchangers. Simplified P&ID style on a 40x40 grid.
import { S } from "./style";

export const shapes: Record<string, JSX.Element> = {
  // Shell and tube: the shell, with a U-tube bundle entering from the left.
  HeatExchanger: (
    <>
      <rect {...S} x="3" y="11" width="34" height="18" rx="9" />
      <path {...S} d="M3 16h24a4 4 0 0 1 0 8H3" strokeWidth={1.4} />
    </>
  ),
};
