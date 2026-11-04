import type { FlowsheetModel } from "./types";

/**
 * The editor's memory in the browser (localStorage): an autosave of the open flowsheet,
 * restored after a reload, and a short list of recently opened or saved flowsheets.
 *
 * Storage can be missing, disabled, full or hold anything (an older version, another app's
 * data), so every helper here catches errors, validates what it reads and never throws.
 * The helpers take the storage as an argument, which keeps them testable.
 */

export type StorageLike = Pick<Storage, "getItem" | "setItem" | "removeItem">;

export const AUTOSAVE_KEY = "fluxsheet:autosave:v1";
export const RECENT_KEY = "fluxsheet:recent:v1";
/** Larger flowsheets are not autosaved; localStorage holds about 5 MB per site. */
export const MAX_AUTOSAVE_CHARS = 2_000_000;
/** Recent entries carry the whole model, so they stay small and few. */
export const MAX_RECENT_CHARS = 500_000;
export const MAX_RECENT = 6;

export interface Autosave {
  savedAt: number;
  /** The flowsheet had unsaved changes (not yet saved to a file). */
  dirty: boolean;
  model: FlowsheetModel;
}

export interface RecentEntry {
  name: string;
  at: number;
  model: FlowsheetModel;
}

/** window.localStorage, or null where it is unavailable (the getter itself may throw). */
export function safeStorage(): StorageLike | null {
  try {
    return typeof window !== "undefined" && window.localStorage ? window.localStorage : null;
  } catch {
    return null;
  }
}

function read(storage: StorageLike | null, key: string): unknown {
  if (!storage) return null;
  try {
    const text = storage.getItem(key);
    return text === null ? null : JSON.parse(text);
  } catch {
    remove(storage, key); // unreadable: forget it
    return null;
  }
}

function write(storage: StorageLike | null, key: string, text: string): boolean {
  if (!storage) return false;
  try {
    storage.setItem(key, text);
    return true;
  } catch {
    return false; // quota exceeded, storage disabled, ...
  }
}

function remove(storage: StorageLike | null, key: string): void {
  try {
    storage?.removeItem(key);
  } catch {
    /* nothing to do */
  }
}

/** Enough of a flowsheet to put on the canvas. */
export function isModel(x: unknown): x is FlowsheetModel {
  const m = x as FlowsheetModel | null;
  return (
    !!m && typeof m === "object" && typeof m.name === "string" && Array.isArray(m.components) &&
    Array.isArray(m.units) && Array.isArray(m.streams) && typeof m.property_method === "string"
  );
}

/** A new, empty flowsheet: nothing worth restoring. */
export function isBlank(m: FlowsheetModel): boolean {
  return (
    m.units.length === 0 && m.streams.length === 0 && m.components.length === 0 &&
    m.name === "Untitled flowsheet" && !m.description && !m.reactions?.length &&
    !m.utilities?.length && !m.design_specs?.length && !m.sensitivity?.length &&
    !m.calculators?.length && !m.dynamics
  );
}

// ---- autosave ------------------------------------------------------------------------------

export function readAutosave(storage: StorageLike | null): Autosave | null {
  const data = read(storage, AUTOSAVE_KEY) as (Autosave & { v?: number }) | null;
  if (data === null) return null;
  if (data?.v !== 1 || !isModel(data.model) || typeof data.savedAt !== "number") {
    remove(storage, AUTOSAVE_KEY);
    return null;
  }
  return { savedAt: data.savedAt, dirty: !!data.dirty, model: data.model };
}

/**
 * Store the open flowsheet. The autosave matters more than the recent list: if storage is
 * full, recent entries are dropped (oldest first) to make room. Returns whether it was stored.
 */
export function writeAutosave(storage: StorageLike | null, model: FlowsheetModel, dirty: boolean, now = Date.now()): boolean {
  if (!storage) return false;
  let text: string;
  try {
    text = JSON.stringify({ v: 1, savedAt: now, dirty, model });
  } catch {
    return false; // not serializable
  }
  if (text.length > MAX_AUTOSAVE_CHARS) {
    remove(storage, AUTOSAVE_KEY); // an older, smaller autosave would restore the wrong flowsheet
    return false;
  }
  if (write(storage, AUTOSAVE_KEY, text)) return true;
  const recent = readRecent(storage);
  while (recent.length > 0) {
    recent.pop();
    if (recent.length) write(storage, RECENT_KEY, JSON.stringify(recent));
    else remove(storage, RECENT_KEY);
    if (write(storage, AUTOSAVE_KEY, text)) return true;
  }
  remove(storage, AUTOSAVE_KEY);
  return false;
}

export function clearAutosave(storage: StorageLike | null): void {
  remove(storage, AUTOSAVE_KEY);
}

export const UNRESTORED_KEY = "fluxsheet:autosave:unrestored";

/**
 * Move an autosave this editor cannot restore (e.g. unknown unit types) out of the way, so the
 * next autosave does not overwrite it. Returns whether it was kept.
 */
export function keepUnrestored(storage: StorageLike | null): boolean {
  try {
    const text = storage?.getItem(AUTOSAVE_KEY);
    if (!storage || text == null) return false;
    storage.setItem(UNRESTORED_KEY, text);
    storage.removeItem(AUTOSAVE_KEY);
    return true;
  } catch {
    return false;
  }
}

// ---- recent flowsheets ---------------------------------------------------------------------

export function readRecent(storage: StorageLike | null): RecentEntry[] {
  const data = read(storage, RECENT_KEY);
  if (data === null) return [];
  if (!Array.isArray(data)) {
    remove(storage, RECENT_KEY);
    return [];
  }
  return data
    .filter((e): e is RecentEntry => !!e && typeof e.at === "number" && isModel(e.model))
    .map((e) => ({ name: typeof e.name === "string" ? e.name : e.model.name, at: e.at, model: e.model }))
    .slice(0, MAX_RECENT);
}

/**
 * Put a flowsheet at the top of the recent list. An entry with the same name and content
 * moves up instead of appearing twice. Returns the new list (unchanged if it could not be stored).
 */
export function addRecent(storage: StorageLike | null, model: FlowsheetModel, now = Date.now()): RecentEntry[] {
  const before = readRecent(storage);
  let body: string;
  try {
    body = JSON.stringify(model);
  } catch {
    return before;
  }
  if (!storage || body.length > MAX_RECENT_CHARS) return before;
  const same = (e: RecentEntry) => e.name === model.name && JSON.stringify(e.model) === body;
  const list = [{ name: model.name, at: now, model }, ...before.filter((e) => !same(e))].slice(0, MAX_RECENT);
  // When storage is full, keep fewer entries rather than none.
  for (let n = list.length; n > 0; n--) {
    if (write(storage, RECENT_KEY, JSON.stringify(list.slice(0, n)))) return list.slice(0, n);
  }
  return before;
}

export function clearRecent(storage: StorageLike | null): void {
  remove(storage, RECENT_KEY);
}
