/**
 * Result overlays on the canvas: what a stream label shows after a run, and which run
 * warnings belong to a unit. Pure functions, so the canvas only maps them over its items.
 */
import { formatQuantity } from "./quantities";
import type { SolveResult, StreamResult } from "./types";

/** What a stream label shows, from just the stream id to its full state. */
export const OVERLAYS = [
  { value: "id", label: "Stream IDs" },
  { value: "T", label: "Temperature" },
  { value: "TP", label: "T, P" },
  { value: "TPF", label: "T, P, flow" },
  { value: "all", label: "T, P, flow, composition" },
] as const;

export type Overlay = (typeof OVERLAYS)[number]["value"];

export const DEFAULT_OVERLAY: Overlay = "T";

const KEY = "fluxsheet:overlay";

export function isOverlay(v: unknown): v is Overlay {
  return OVERLAYS.some((o) => o.value === v);
}

/** The overlay this viewer chose last, or the default (storage may be missing or throw). */
export function loadOverlay(): Overlay {
  try {
    const v = window.localStorage.getItem(KEY);
    return isOverlay(v) ? v : DEFAULT_OVERLAY;
  } catch {
    return DEFAULT_OVERLAY;
  }
}

export function saveOverlay(v: Overlay): void {
  try {
    window.localStorage.setItem(KEY, v);
  } catch {
    // a convenience only
  }
}

/** The leading component by mole fraction, e.g. "92 % water" (nothing for an empty stream). */
function mainComponent(s: StreamResult): string | null {
  let best: [string, number] | null = null;
  for (const [name, x] of Object.entries(s.mole_fractions ?? {})) {
    if (best === null || x > best[1]) best = [name, x];
  }
  if (best === null || !(s.molar_flow > 0)) return null;
  const pct = best[1] * 100;
  return `${pct >= 99.95 ? "100" : pct.toFixed(pct >= 10 ? 0 : 1)} % ${best[0]}`;
}

/** The label of stream `id`: its id, then the chosen results of the last run. */
export function streamLabel(id: string, s: StreamResult | undefined, overlay: Overlay): string {
  if (!s || overlay === "id") return id;
  const parts = [id, formatQuantity(s.T, "temperature")];
  if (overlay !== "T") parts.push(formatQuantity(s.P, "pressure"));
  if (overlay === "TPF" || overlay === "all") parts.push(formatQuantity(s.molar_flow, "molar_flow"));
  if (overlay === "all") {
    const main = mainComponent(s);
    if (main) parts.push(main);
  }
  return parts.join(" · ");
}

/**
 * Warnings of the last run that concern unit `id`: unit notes (reported as "ID: note") and
 * balance checks scoped to the unit that did not close.
 */
export function unitWarnings(id: string, result: SolveResult | null | undefined): string[] {
  if (!result?.ok) return [];
  const prefix = `${id}: `;
  const out = result.warnings.filter((w) => w.startsWith(prefix)).map((w) => w.slice(prefix.length));
  for (const b of result.balance_checks ?? []) {
    if (b.scope !== id || b.status !== "warning") continue;
    const what = b.component ? `${b.quantity} balance of ${b.component}` : `${b.quantity} balance`;
    out.push(`${what} does not close${b.residual != null ? ` (residual ${b.residual.toPrecision(3)} ${b.unit})` : ""}`);
  }
  return out;
}
