// Symbols for columns. Simplified P&ID style on a 40x40 grid.
import { S } from "./style";

export const shapes: Record<string, JSX.Element> = {
  ShortcutColumn: (
    <>
      <rect {...S} x="13" y="3" width="14" height="34" rx="7" />
      <path {...S} d="M13 14h14M13 20h14M13 26h14" strokeDasharray="3 2" strokeWidth={1.2} />
    </>
  ),
  DistillationColumn: (
    <>
      <rect {...S} x="10" y="3" width="13" height="34" rx="6.5" />
      <path {...S} d="M10 12h13M10 17h13M10 22h13M10 27h13" strokeWidth={1.2} />
      <path {...S} d="M17 3V1h14v4" strokeWidth={1.4} />
      <circle {...S} cx="31" cy="8" r="3" strokeWidth={1.4} />
      <path {...S} d="M17 37v2h14v-4" strokeWidth={1.4} />
      <circle {...S} cx="31" cy="32" r="3" strokeWidth={1.4} />
    </>
  ),
  Absorber: (
    <>
      <rect {...S} x="12" y="3" width="16" height="34" rx="8" />
      <path {...S} d="M12 10l16 4M12 16l16 4M12 22l16 4M12 28l16 4" strokeWidth={1.2} />
    </>
  ),
};
