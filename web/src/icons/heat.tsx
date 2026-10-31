// Symbol for the two-sided heat exchanger: shell and tube, hot stream on the shell side,
// cold stream through the tubes.
import { hot, cold, metalH, thin, type UnitSymbol } from "./style";

export const symbols: Record<string, UnitSymbol> = {
  HeatExchanger: {
    w: 110,
    h: 56,
    body: (
      <>
        <path {...hot} strokeWidth={2.2} strokeLinecap="butt" d="M34 14V6H1M82 42V50H109" />
        <path {...cold} strokeWidth={2.2} strokeLinecap="butt" d="M1 36H10M98 28H109" />
        <path {...metalH} d="M24 14H90A8 14 0 0 1 90 42H24Z" />
        <path {...thin} d="M28 21H88M28 28H92M28 35H88" />
        <rect {...metalH} x="10" y="12" width="13" height="32" rx="3" />
        <rect {...metalH} x="21" y="9" width="5" height="38" rx="1" />
      </>
    ),
    ports: {
      hot_in: { x: 1, y: 6 },
      cold_in: { x: 1, y: 36 },
      hot_out: { x: 109, y: 50 },
      cold_out: { x: 109, y: 28 },
    },
  },
};
