// Simplified P&ID-style symbols, drawn with currentColor so they follow the theme.
const S = { fill: "none", stroke: "currentColor", strokeWidth: 2, strokeLinejoin: "round" as const };

const shapes: Record<string, JSX.Element> = {
  Feed: (
    <>
      <path {...S} d="M6 14h20l6 6-6 6H6z" />
    </>
  ),
  Product: (
    <>
      <path {...S} d="M8 14h20l6 6-6 6H8l6-6z" />
    </>
  ),
  Mixer: (
    <>
      <path {...S} d="M8 8l26 12-26 12z" />
    </>
  ),
  Splitter: (
    <>
      <path {...S} d="M32 8L6 20l26 12z" />
    </>
  ),
  Heater: (
    <>
      <circle {...S} cx="20" cy="20" r="14" />
      <path {...S} d="M6 20h6l3-6 5 12 5-12 3 6h6" />
    </>
  ),
  Flash: (
    <>
      <rect {...S} x="10" y="3" width="20" height="34" rx="10" />
      <path {...S} d="M10 24h20" strokeDasharray="3 3" />
    </>
  ),
  ComponentSplitter: (
    <>
      <rect {...S} x="12" y="2" width="16" height="36" rx="8" />
      <path {...S} d="M12 12h16M12 18h16M12 24h16M12 30h16" strokeWidth={1.2} />
    </>
  ),
  Valve: (
    <>
      <path {...S} d="M6 12l28 16V12L6 28z" />
    </>
  ),
  Pump: (
    <>
      <circle {...S} cx="20" cy="20" r="14" />
      <path {...S} d="M20 6h14M12 30l-4 6h24l-4-6" />
    </>
  ),
  Compressor: (
    <>
      <path {...S} d="M6 8l28 7v10L6 32z" />
    </>
  ),
};

export function UnitIcon({ type, size = 40 }: { type: string; size?: number }) {
  return (
    <svg width={size} height={size} viewBox="0 0 40 40" aria-hidden="true">
      {shapes[type] ?? <rect {...S} x="6" y="6" width="28" height="28" rx="4" />}
    </svg>
  );
}
