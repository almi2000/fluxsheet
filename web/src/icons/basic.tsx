// Symbols for the basic units: stream tags, mixer and splitter, heater, flash drum,
// component separator, valve, pump and compressor.
import { blue, blueH, hot, liquid, metal, metalH, metalR, outline, pipe, thin, type UnitSymbol } from "./style";

const FLASH = "M9 14A16 8 0 0 1 41 14V74A16 8 0 0 1 9 74Z";

export const symbols: Record<string, UnitSymbol> = {
  Feed: {
    w: 46,
    h: 26,
    body: <path d="M2 4h30l12 9-12 9H2z" fill="var(--sym-feed)" stroke="var(--sym-feed-line)" strokeWidth={1.4} strokeLinejoin="round" />,
    ports: { out: { x: 44, y: 13 } },
  },
  Product: {
    w: 46,
    h: 26,
    body: <path d="M2 4h30l12 9-12 9H2l8-9z" fill="var(--sym-product)" stroke="var(--sym-product-line)" strokeWidth={1.4} strokeLinejoin="round" />,
    ports: { in: { x: 9, y: 13 } },
  },
  Mixer: {
    w: 44,
    h: 44,
    body: (
      <>
        <path {...metalH} d="M4 4L40 22L4 40z" />
        <path {...thin} d="M8 12L24 20" />
      </>
    ),
    ports: { in: { x: 4, y: 22 }, out: { x: 40, y: 22 } },
  },
  Splitter: {
    w: 44,
    h: 44,
    body: (
      <>
        <path {...metalH} d="M40 4L4 22L40 40z" />
        <path {...thin} d="M36 10L20 18" />
      </>
    ),
    ports: { in: { x: 4, y: 22 }, out1: { x: 40, y: 12 }, out2: { x: 40, y: 32 } },
  },
  Heater: {
    w: 48,
    h: 48,
    body: (
      <>
        <path {...pipe} d="M1 24h5M42 24h5" />
        <circle {...metalR} cx="24" cy="24" r="18" />
        <path {...hot} d="M6 24h6l4-8 6 16 6-16 6 16 4-8h4" />
      </>
    ),
    ports: { in: { x: 1, y: 24 }, out: { x: 47, y: 24 } },
  },
  Flash: {
    w: 50,
    h: 94,
    defs: (
      <clipPath id="fx-clip-flash">
        <path d={FLASH} />
      </clipPath>
    ),
    body: (
      <>
        <path {...pipe} d="M25 6V3H49M1 44H9M41 72H49" />
        <path {...metal} d="M15 80h20v12H15z" />
        <path {...metal} d={FLASH} />
        <rect {...liquid} x="0" y="56" width="50" height="30" clipPath="url(#fx-clip-flash)" />
        <path {...outline} d="M9 56H41" />
        <rect x="11" y="17" width="28" height="6" fill="url(#fx-mesh)" stroke="var(--sym-line)" strokeWidth={0.8} />
        <path {...outline} d={FLASH} />
      </>
    ),
    ports: { in: { x: 1, y: 44 }, vapor: { x: 49, y: 3 }, liquid: { x: 49, y: 72 } },
  },
  ComponentSplitter: {
    w: 48,
    h: 80,
    body: (
      <>
        <path {...pipe} d="M1 40H8M40 12H47M40 68H47" />
        <rect {...metal} x="8" y="4" width="32" height="72" rx="5" />
        <path {...outline} d="M8 40l6-4 7 8 6-8 7 8 6-4" />
        <path {...thin} d="M24 30V12M19 17l5-5 5 5M24 50v18M19 63l5 5 5-5" />
      </>
    ),
    ports: { in: { x: 1, y: 40 }, top: { x: 47, y: 12 }, bottom: { x: 47, y: 68 } },
  },
  Valve: {
    w: 44,
    h: 40,
    body: (
      <>
        <path {...outline} d="M22 28V13" />
        <path {...blue} d="M12 13a10 9 0 0 1 20 0z" />
        <path {...metalH} d="M4 20L22 28L4 36zM40 20L22 28L40 36z" />
      </>
    ),
    ports: { in: { x: 4, y: 28 }, out: { x: 40, y: 28 } },
  },
  Pump: {
    w: 52,
    h: 48,
    body: (
      <>
        <path {...metal} d="M12 46l4-8h18l4 8z" />
        <rect {...blueH} x="2" y="22" width="8" height="8" />
        <rect {...blueH} x="24" y="9" width="24" height="9" />
        <circle {...blue} cx="25" cy="26" r="16" />
        <circle {...metalR} cx="25" cy="26" r="5" />
        <path {...outline} d="M48 8v11" />
      </>
    ),
    ports: { in: { x: 2, y: 26 }, out: { x: 50, y: 13.5 } },
  },
  Compressor: {
    w: 58,
    h: 48,
    body: (
      <>
        <path {...pipe} d="M52 24h5" />
        <path {...blue} d="M4 4L52 15V33L4 44z" />
        <path {...thin} stroke="var(--sym-shine)" d="M16 9v30M28 12v24M40 14v20" />
      </>
    ),
    ports: { in: { x: 4, y: 24 }, out: { x: 57, y: 24 } },
  },
};
