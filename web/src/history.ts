import { useCallback, useReducer, useRef } from "react";

/**
 * Undo/redo as a stack of document snapshots.
 *
 * The present document lives in the caller's state; the history only keeps what came before
 * (`past`) and what was undone (`future`). Every snapshot carries a version number so the
 * editor can tell whether the present equals the last saved document, even after undo/redo.
 */

export interface Entry<T> {
  doc: T;
  version: number;
}

export interface History<T> {
  past: Entry<T>[];
  future: Entry<T>[];
  /** Version of the present document. */
  version: number;
  nextVersion: number;
  /** Key and time of the last recorded step, for coalescing rapid edits. */
  lastKey: string | null;
  lastAt: number;
}

export const HISTORY_LIMIT = 200;
/** Edits with the same key closer together than this merge into one undo step. */
export const COALESCE_MS = 1000;

export function emptyHistory<T>(): History<T> {
  return { past: [], future: [], version: 0, nextVersion: 1, lastKey: null, lastAt: 0 };
}

/**
 * Record that the document is about to change from `before`.
 * With a `key`, consecutive edits with the same key within COALESCE_MS form one step.
 */
export function record<T>(h: History<T>, before: T, key: string | null = null, now = Date.now()): History<T> {
  if (key !== null && key === h.lastKey && now - h.lastAt < COALESCE_MS) {
    return { ...h, future: [], lastAt: now };
  }
  const past = [...h.past, { doc: before, version: h.version }].slice(-HISTORY_LIMIT);
  return { past, future: [], version: h.nextVersion, nextVersion: h.nextVersion + 1, lastKey: key, lastAt: now };
}

/** Step back: returns the new history and the document to show, or null if there is nothing to undo. */
export function undo<T>(h: History<T>, current: T): [History<T>, T] | null {
  const prev = h.past[h.past.length - 1];
  if (!prev) return null;
  return [
    {
      ...h,
      past: h.past.slice(0, -1),
      future: [{ doc: current, version: h.version }, ...h.future],
      version: prev.version,
      lastKey: null,
    },
    prev.doc,
  ];
}

export function redo<T>(h: History<T>, current: T): [History<T>, T] | null {
  const next = h.future[0];
  if (!next) return null;
  return [
    {
      ...h,
      past: [...h.past, { doc: current, version: h.version }],
      future: h.future.slice(1),
      version: next.version,
      lastKey: null,
    },
    next.doc,
  ];
}

/** Ends any coalescing run, so the next edit starts a new step. */
export function seal<T>(h: History<T>): History<T> {
  return h.lastKey === null ? h : { ...h, lastKey: null };
}

/**
 * React binding: `current()` reads the latest document, `apply(doc)` puts one back on screen.
 * `version` changes whenever the present document does; compare it with a saved version to
 * know whether there are unsaved changes.
 */
export function useHistory<T>(current: () => T, apply: (doc: T) => void) {
  const ref = useRef<History<T>>(emptyHistory<T>());
  const [, rerender] = useReducer((n: number) => n + 1, 0);

  const set = useCallback((h: History<T>) => {
    ref.current = h;
    rerender();
  }, []);

  /** Call right before a model-changing action; returns the version of the coming document. */
  const commit = useCallback(
    (key: string | null = null) => {
      const h = record(ref.current, current(), key);
      set(h);
      return h.version;
    },
    [current, set],
  );
  /** Record a step whose "before" document was captured earlier (e.g. at drag start). */
  const commitFrom = useCallback((before: T) => set(record(ref.current, before)), [set]);

  const step = useCallback(
    (dir: "undo" | "redo") => {
      const r = (dir === "undo" ? undo : redo)(ref.current, current());
      if (!r) return false;
      set(r[0]);
      apply(r[1]);
      return true;
    },
    [current, apply, set],
  );

  const h = ref.current;
  return {
    commit,
    commitFrom,
    undo: useCallback(() => step("undo"), [step]),
    redo: useCallback(() => step("redo"), [step]),
    seal: useCallback(() => set(seal(ref.current)), [set]),
    /** Forgets all steps, keeping the present version. */
    reset: useCallback(() => set({ ...emptyHistory<T>(), version: ref.current.version, nextVersion: ref.current.nextVersion }), [set]),
    canUndo: h.past.length > 0,
    canRedo: h.future.length > 0,
    version: h.version,
  };
}
