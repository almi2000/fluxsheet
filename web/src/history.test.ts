import { describe, expect, it } from "vitest";
import { COALESCE_MS, HISTORY_LIMIT, emptyHistory, record, redo, seal, undo } from "./history";

describe("history", () => {
  it("undoes and redoes steps in order", () => {
    let h = emptyHistory<string>();
    h = record(h, "a"); // a -> b
    h = record(h, "b"); // b -> c
    const u1 = undo(h, "c")!;
    expect(u1[1]).toBe("b");
    const u2 = undo(u1[0], "b")!;
    expect(u2[1]).toBe("a");
    expect(undo(u2[0], "a")).toBeNull();
    const r1 = redo(u2[0], "a")!;
    expect(r1[1]).toBe("b");
    expect(redo(r1[0], "b")![1]).toBe("c");
  });

  it("drops the redo branch on a new edit", () => {
    let h = record(emptyHistory<string>(), "a");
    h = undo(h, "b")![0];
    expect(h.future).toHaveLength(1);
    h = record(h, "a");
    expect(h.future).toHaveLength(0);
  });

  it("coalesces edits with the same key within the window", () => {
    let h = emptyHistory<string>();
    h = record(h, "", "name", 0);
    h = record(h, "F", "name", 300);
    h = record(h, "Fl", "name", 600);
    expect(h.past.map((e) => e.doc)).toEqual([""]);
    // A pause, another key, or seal() start a new step.
    h = record(h, "Flu", "name", 600 + COALESCE_MS + 1);
    h = record(h, "Flux", "desc", 1700);
    h = record(seal(h), "Fluxs", "desc", 1800);
    expect(h.past.map((e) => e.doc)).toEqual(["", "Flu", "Flux", "Fluxs"]);
  });

  it("does not coalesce across an undo", () => {
    let h = record(emptyHistory<string>(), "a", "k", 0);
    h = undo(h, "ab")![0];
    h = record(h, "a", "k", 10);
    expect(h.past).toHaveLength(1);
    expect(h.future).toHaveLength(0);
  });

  it("tracks versions so undoing back to the saved document is clean", () => {
    let h = emptyHistory<string>();
    h = record(h, "a");
    const saved = h.version;
    h = record(h, "b");
    expect(h.version).not.toBe(saved);
    h = undo(h, "c")![0];
    expect(h.version).toBe(saved);
    h = redo(h, "b")![0];
    expect(h.version).not.toBe(saved);
  });

  it("keeps at most HISTORY_LIMIT steps", () => {
    let h = emptyHistory<number>();
    for (let i = 0; i < HISTORY_LIMIT + 20; i++) h = record(h, i);
    expect(h.past).toHaveLength(HISTORY_LIMIT);
    expect(h.past[0].doc).toBe(20);
  });
});
