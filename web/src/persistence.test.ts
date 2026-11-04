import { describe, expect, it } from "vitest";
import {
  AUTOSAVE_KEY, MAX_AUTOSAVE_CHARS, MAX_RECENT, RECENT_KEY, addRecent, clearAutosave, isBlank,
  UNRESTORED_KEY, keepUnrestored, readAutosave, readRecent, writeAutosave, type StorageLike,
} from "./persistence";
import type { FlowsheetModel } from "./types";

/** In-memory storage with an optional size limit (like a browser quota) and failure modes. */
function memory({ limit = Infinity, failGet = false, failSet = false } = {}) {
  const data = new Map<string, string>();
  const size = () => [...data.values()].reduce((n, v) => n + v.length, 0);
  const storage: StorageLike = {
    getItem: (k) => {
      if (failGet) throw new Error("SecurityError");
      return data.get(k) ?? null;
    },
    setItem: (k, v) => {
      if (failSet) throw new Error("SecurityError");
      if (size() - (data.get(k)?.length ?? 0) + v.length > limit) throw new Error("QuotaExceededError");
      data.set(k, v);
    },
    removeItem: (k) => {
      data.delete(k);
    },
  };
  return { storage, data };
}

function model(name = "Flash", units = 1): FlowsheetModel {
  return {
    format: "fluxsheet/1", name, components: ["water"], property_method: "auto",
    units: Array.from({ length: units }, (_, i) => ({ id: `U${i}`, type: "Feed", params: { note: "x".repeat(50) } })),
    streams: [],
  };
}

describe("autosave", () => {
  it("round-trips the flowsheet and its dirty flag", () => {
    const { storage } = memory();
    expect(writeAutosave(storage, model(), true, 1234)).toBe(true);
    expect(readAutosave(storage)).toEqual({ savedAt: 1234, dirty: true, model: model() });
    clearAutosave(storage);
    expect(readAutosave(storage)).toBeNull();
  });

  it("never throws when storage fails", () => {
    const { storage } = memory({ failGet: true, failSet: true });
    expect(writeAutosave(storage, model(), false)).toBe(false);
    expect(readAutosave(storage)).toBeNull();
    expect(addRecent(storage, model())).toEqual([]);
    expect(readRecent(storage)).toEqual([]);
    expect(writeAutosave(null, model(), false)).toBe(false);
    expect(readAutosave(null)).toBeNull();
  });

  it("drops corrupt or foreign data", () => {
    for (const bad of ["{not json", JSON.stringify({ v: 2, savedAt: 1, model: model() }),
      JSON.stringify({ v: 1, savedAt: 1, model: { name: "x" } }), "null", "42"]) {
      const { storage, data } = memory();
      data.set(AUTOSAVE_KEY, bad);
      expect(readAutosave(storage)).toBeNull();
      expect(data.has(AUTOSAVE_KEY)).toBe(bad === "null" ? true : false);
    }
  });

  it("skips flowsheets too large to store and removes the older autosave", () => {
    const { storage, data } = memory();
    writeAutosave(storage, model("small"), false);
    const huge = model("huge", Math.ceil(MAX_AUTOSAVE_CHARS / 80));
    expect(writeAutosave(storage, huge, true)).toBe(false);
    expect(data.has(AUTOSAVE_KEY)).toBe(false);
  });

  it("makes room by dropping recent entries when storage is full", () => {
    const { storage } = memory({ limit: 6000 });
    for (let i = 0; i < 5; i++) addRecent(storage, model(`R${i}`, 10), i);
    expect(readRecent(storage).length).toBeGreaterThan(1);
    expect(writeAutosave(storage, model("big", 40), true)).toBe(true);
    expect(readAutosave(storage)?.model.name).toBe("big");
    expect(readRecent(storage).length).toBeLessThan(5);
  });

  it("keeps an autosave it cannot restore aside", () => {
    const { storage, data } = memory();
    writeAutosave(storage, model("future"), true);
    expect(keepUnrestored(storage)).toBe(true);
    expect(data.has(AUTOSAVE_KEY)).toBe(false);
    expect(JSON.parse(data.get(UNRESTORED_KEY)!).model.name).toBe("future");
    expect(keepUnrestored(memory({ failGet: true }).storage)).toBe(false);
  });

  it("knows an empty flowsheet", () => {
    expect(isBlank({ ...model("Untitled flowsheet", 0), components: [] })).toBe(true);
    expect(isBlank({ ...model("Untitled flowsheet", 0), components: ["water"] })).toBe(false);
  });
});

describe("recent flowsheets", () => {
  it("keeps the newest first, without duplicates, at most MAX_RECENT", () => {
    const { storage } = memory();
    for (let i = 0; i < MAX_RECENT + 2; i++) addRecent(storage, model(`F${i}`), i);
    addRecent(storage, model("F7"), 100); // same name and content: moves to the top
    const list = readRecent(storage);
    expect(list.length).toBe(MAX_RECENT);
    expect(list.map((e) => e.name).slice(0, 2)).toEqual(["F7", "F6"]);
    expect(list[0].at).toBe(100);
  });

  it("keeps two different flowsheets of the same name", () => {
    const { storage } = memory();
    addRecent(storage, model("Untitled flowsheet", 1));
    addRecent(storage, model("Untitled flowsheet", 2));
    expect(readRecent(storage).length).toBe(2);
  });

  it("ignores invalid entries", () => {
    const { storage, data } = memory();
    data.set(RECENT_KEY, JSON.stringify([{ name: "ok", at: 1, model: model() }, { at: 2 }, null]));
    expect(readRecent(storage).map((e) => e.name)).toEqual(["ok"]);
    data.set(RECENT_KEY, "{}");
    expect(readRecent(storage)).toEqual([]);
  });

  it("stores fewer entries rather than none when storage is full", () => {
    const { storage } = memory({ limit: 3000 });
    for (let i = 0; i < 4; i++) addRecent(storage, model(`F${i}`, 5), i);
    const list = readRecent(storage);
    expect(list[0].name).toBe("F3");
    expect(list.length).toBeGreaterThanOrEqual(1);
  });
});
