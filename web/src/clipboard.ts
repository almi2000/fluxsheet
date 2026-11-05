import type { CustomCompound, FlowsheetModel, StreamModel, UnitModel } from "./types";

/**
 * Copy/paste of flowsheet fragments. A fragment is an ordinary Fluxsheet file holding the
 * copied units and the streams between them, so it can travel through the system clipboard
 * into another editor tab (or a text editor) and back.
 */

export type Point = { x: number; y: number };

const unitOf = (port: string) => port.split(".")[0];

/** The selected units and the streams that run between two of them. */
export function extractFragment(model: FlowsheetModel, ids: ReadonlySet<string>): FlowsheetModel {
  return {
    format: model.format,
    name: model.name,
    components: model.components,
    property_method: model.property_method,
    // The definitions of user-defined components travel along, for pasting into another flowsheet.
    ...(model.custom_compounds?.length && {
      custom_compounds: model.custom_compounds.filter((c) => model.components.includes(c.name)),
    }),
    units: model.units.filter((u) => ids.has(u.id)),
    streams: model.streams.filter((s) => s.to && ids.has(unitOf(s.from)) && ids.has(unitOf(s.to))),
  };
}

/** Reads a fragment from clipboard text; null if the text is not a Fluxsheet file. */
export function parseFragment(text: string): FlowsheetModel | null {
  try {
    const m = JSON.parse(text);
    if (!m || typeof m !== "object" || typeof m.format !== "string" || !m.format.startsWith("fluxsheet/")) return null;
    if (!Array.isArray(m.units) || m.units.some((u: UnitModel) => typeof u?.id !== "string" || typeof u?.type !== "string"))
      return null;
    const custom = Array.isArray(m.custom_compounds)
      ? m.custom_compounds.filter((c: CustomCompound) => typeof c?.name === "string" && typeof c?.formula === "string")
      : undefined;
    return {
      ...m,
      streams: Array.isArray(m.streams) ? m.streams : [],
      components: Array.isArray(m.components) ? m.components.filter((c: unknown) => typeof c === "string") : [],
      custom_compounds: custom,
    };
  } catch {
    return null;
  }
}

/** `MIX1` -> `MIX2`, `CHILLER` -> `CHILLER1`: same name stem, first number not taken. */
export function freshId(id: string, taken: ReadonlySet<string>, fallback = "U"): string {
  const stem = id.replace(/\d+$/, "") || fallback;
  let i = 1;
  while (taken.has(`${stem}${i}`)) i++;
  return `${stem}${i}`;
}

/** Top-left corner of the fragment's units; units without positions count as the origin. */
export function fragmentOrigin(units: UnitModel[]): Point {
  if (units.length === 0) return { x: 0, y: 0 };
  return {
    x: Math.min(...units.map((u) => u.ui?.x ?? 0)),
    y: Math.min(...units.map((u) => u.ui?.y ?? 0)),
  };
}

/**
 * Gives every unit and stream of a fragment an id not in `takenUnits` / `takenStreams`
 * (streams become S1, S2, ...), rewires the streams accordingly and moves the units so the
 * fragment's top-left corner lands on `at`. Streams to units outside the fragment are dropped.
 */
export function remapFragment(
  fragment: FlowsheetModel,
  takenUnits: ReadonlySet<string>,
  takenStreams: ReadonlySet<string>,
  at: Point,
): { units: UnitModel[]; streams: StreamModel[]; idMap: Map<string, string> } {
  const usedUnits = new Set(takenUnits);
  const usedStreams = new Set(takenStreams);
  const idMap = new Map<string, string>();
  const origin = fragmentOrigin(fragment.units);

  const units = fragment.units.map((u, i) => {
    const id = freshId(u.id, usedUnits);
    usedUnits.add(id);
    idMap.set(u.id, id);
    // Units without a position (hand-written JSON) get a simple grid.
    const pos = u.ui ?? { x: origin.x + (i % 5) * 200, y: origin.y + Math.floor(i / 5) * 180 };
    return { ...u, id, ui: { x: Math.round(pos.x - origin.x + at.x), y: Math.round(pos.y - origin.y + at.y) } };
  });

  const streams: StreamModel[] = [];
  const streamMap = new Map<string, string>();
  for (const s of fragment.streams) {
    if (!s.to) continue;
    const [fromUnit, fromPort] = s.from.split(".");
    const [toUnit, toPort] = s.to.split(".");
    const a = idMap.get(fromUnit);
    const b = idMap.get(toUnit);
    if (!a || !b) continue;
    const id = freshId("S", usedStreams);
    usedStreams.add(id);
    streamMap.set(s.id, id);
    streams.push({ ...s, id, from: `${a}.${fromPort}`, to: `${b}.${toPort}` });
  }
  // Parameters keyed by stream id (a column's feed stage per stream) follow the new ids.
  const oldStreams = new Set(fragment.streams.map((s) => s.id));
  const isStreamMap = (v: unknown): v is Record<string, unknown> =>
    !!v && typeof v === "object" && !Array.isArray(v) && Object.keys(v).length > 0 &&
    Object.keys(v).every((k) => oldStreams.has(k));
  for (const u of units) {
    for (const [name, v] of Object.entries(u.params ?? {})) {
      if (!isStreamMap(v)) continue;
      const moved = Object.entries(v).flatMap(([k, x]) => (streamMap.has(k) ? [[streamMap.get(k)!, x]] : []));
      u.params = { ...u.params, [name]: Object.fromEntries(moved) };
    }
  }
  return { units, streams, idMap };
}

/** JSON with sorted keys, so two definitions compare equal whatever their key order. */
export function canonical(x: unknown): string {
  if (Array.isArray(x)) return `[${x.map(canonical).join(",")}]`;
  if (x && typeof x === "object")
    return `{${Object.keys(x)
      .filter((k) => (x as Record<string, unknown>)[k] !== undefined)
      .sort()
      .map((k) => `${JSON.stringify(k)}:${canonical((x as Record<string, unknown>)[k])}`)
      .join(",")}}`;
  return JSON.stringify(x);
}

/**
 * The components and custom compound definitions of `target` after pasting `fragment`: the
 * fragment's missing components are added, with their definitions if they are user-defined.
 * A name the target already uses keeps the target's meaning; if the fragment means another
 * compound by it (another definition, or custom versus database), it is listed in `conflicts`.
 * Names are matched exactly, as the backend does.
 */
export function mergeComponents(
  target: { components: string[]; custom_compounds?: CustomCompound[] },
  fragment: FlowsheetModel,
  /** The components the pasted units refer to; only they are reported as conflicts. */
  used?: ReadonlySet<string>,
): { components: string[]; custom_compounds: CustomCompound[]; conflicts: string[] } {
  const components = [...target.components];
  const custom = [...(target.custom_compounds ?? [])];
  const conflicts: string[] = [];
  const theirs = new Map((fragment.custom_compounds ?? []).map((c) => [c.name, c]));
  for (const name of fragment.components ?? []) {
    const def = theirs.get(name);
    const mine = custom.find((c) => c.name === name);
    if (!components.includes(name)) {
      components.push(name);
      if (def && !mine) custom.push(def);
      else if (def && mine && canonical(def) !== canonical(mine)) conflicts.push(name);
    } else if (def ? !mine || canonical(def) !== canonical(mine) : !!mine) {
      conflicts.push(name);
    }
  }
  return { components, custom_compounds: custom, conflicts: used ? conflicts.filter((c) => used.has(c)) : conflicts };
}
