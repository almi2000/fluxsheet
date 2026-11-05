import { describe, expect, it } from "vitest";
import {
  AUTOSAVE_KEY, MAX_AUTOSAVE_CHARS, MAX_RECENT, RECENT_KEY, addRecent, clearAutosave, isBlank,
  UNRESTORED_KEY, keepUnrestored, otherTabAutosave, readAutosave, readRecent, tabId, writeAutosave, type StorageLike,
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

describe("autosave with several tabs", () => {
  const known = (m: FlowsheetModel) => JSON.stringify(m);

  it("stores and reads the tab that wrote the autosave", () => {
    const { storage } = memory();
    writeAutosave(storage, model(), false, 1, "tab-a");
    expect(readAutosave(storage)?.tab).toBe("tab-a");
  });

  it("reports another tab's newer changes, not our own or ones we know", () => {
    const { storage } = memory();
    const a = model("A");
    writeAutosave(storage, a, true, 1, "tab-a");
    expect(otherTabAutosave(storage, "tab-a", known(a))).toBeNull(); // our own
    expect(otherTabAutosave(storage, "tab-b", known(a))).toBeNull(); // b restored the same flowsheet
    const conflict = otherTabAutosave(storage, "tab-b", known(model("B")));
    expect(conflict?.model.name).toBe("A");
    expect(conflict?.dirty).toBe(true);
    expect(otherTabAutosave(storage, "tab-b", null)?.tab).toBe("tab-a");
  });

  it("treats autosaves of older versions (no tab) and missing ones as no conflict", () => {
    const { storage, data } = memory();
    expect(otherTabAutosave(storage, "tab-b", null)).toBeNull();
    writeAutosave(storage, model("Old"), false, 1);
    expect(otherTabAutosave(storage, "tab-b", null)).toBeNull();
    data.set(AUTOSAVE_KEY, "{broken");
    expect(otherTabAutosave(storage, "tab-b", null)).toBeNull();
  });

  it("never throws when storage fails", () => {
    expect(otherTabAutosave(memory({ failGet: true }).storage, "t", null)).toBeNull();
    expect(otherTabAutosave(null, "t", null)).toBeNull();
    expect(typeof tabId()).toBe("string");
  });

  it("gives every page load its own tab id", () => {
    expect(tabId()).not.toBe(tabId());
  });

  it("lets a tab notice another tab's edit before overwriting it", () => {
    // Tab A autosaves; tab B opens the same flowsheet; A edits again; B's next autosave must stop.
    const { storage } = memory();
    const v1 = model("Shared");
    writeAutosave(storage, v1, false, 1, "a");
    const knownB = known(v1); // B restored v1
    expect(otherTabAutosave(storage, "b", knownB)).toBeNull();
    writeAutosave(storage, v1, false, 2, "b"); // B re-saves it after restoring
    expect(otherTabAutosave(storage, "a", known(v1))).toBeNull(); // A may overwrite: nothing new
    const v2 = { ...v1, description: "edited in A" };
    writeAutosave(storage, v2, true, 3, "a");
    expect(otherTabAutosave(storage, "b", knownB)?.model.description).toBe("edited in A");
  });
});
