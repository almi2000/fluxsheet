// Shared paints for the equipment symbols. Colors come from the --sym-* tokens in
// styles.css, so the symbols follow the theme; gradients and patterns live in SymbolDefs.

/** Where a stream connects, in symbol pixels. Inlets sit on the left, outlets on the right. */
export type PortAt = { x: number; y: number };

export interface UnitSymbol {
  w: number;
  h: number;
  body: JSX.Element;
  ports: Record<string, PortAt>;
  /** Clip paths and the like this symbol needs, mounted once in SymbolDefs. */
  defs?: JSX.Element;
}

const line = { stroke: "var(--sym-line)", strokeWidth: 1.2, strokeLinejoin: "round" as const };

/** Vertical vessel shell (shaded across). */
export const metal = { ...line, fill: "url(#fx-metal)" };
/** Horizontal vessel shell (shaded top to bottom). */
export const metalH = { ...line, fill: "url(#fx-metal-h)" };
/** Round casings. */
export const metalR = { ...line, fill: "url(#fx-metal-r)" };
/** Rotating machinery. */
export const blue = { ...line, fill: "url(#fx-blue)" };
export const blueH = { ...line, fill: "url(#fx-blue-h)" };
export const liquid = { fill: "url(#fx-liquid)" };
export const oil = { fill: "url(#fx-oil)" };
export const outline = { ...line, fill: "none" };
export const thin = { ...outline, strokeWidth: 0.8, opacity: 0.6 };
/** Nozzles and connecting pipes. */
export const pipe = { fill: "none", stroke: "var(--sym-pipe)", strokeWidth: 2.2, strokeLinejoin: "round" as const };
export const hot = { fill: "none", stroke: "var(--sym-hot)", strokeWidth: 1.8, strokeLinejoin: "round" as const, strokeLinecap: "round" as const };
export const cold = { ...hot, stroke: "var(--sym-cold)" };
