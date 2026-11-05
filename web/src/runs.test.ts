import { describe, expect, it } from "vitest";
import { rowKeys } from "./Inspector";
import { currentRuns, emptyRuns, runKey } from "./runs";
import { idleSideOutlets, idleSideOutletText } from "./sideOutlets";
import type { FlowsheetModel, UnitType } from "./types";

const model: FlowsheetModel = {
  format: "fluxsheet/1",
  name: "Runs",
  components: ["water"],
  property_method: "auto",
  units: [{ id: "F1", type: "Feed", params: { T: 300, P: 1e5 }, ui: { x: 0, y: 0 } }],
  streams: [],
  sensitivity: [{ id: "A", vary: [{ unit: "F1", param: "T", lower: 300, upper: 350 }], outputs: [] }],
  optimizations: [{ id: "O", variables: {}, objective: "1", vary: [] }],
  dynamics: { t_end: 10, units: {} },
};

describe("run fingerprints", () => {
  const key = runKey(model, "optimization", "O");

  it("ignore positions, names, key order and the other runs' definitions", () => {
    const moved = { ...model, name: "x", units: [{ ...model.units[0], ui: { x: 50, y: 9 }, params: { P: 1e5, T: 300 } }] };
    expect(runKey(moved, "optimization", "O")).toBe(key);
    expect(runKey({ ...model, sensitivity: [] }, "optimization", "O")).toBe(key);
    expect(runKey({ ...model, costing: {} }, "dynamics")).toBe(runKey(model, "dynamics"));
  });

  it("change with the flowsheet or the run's own definition", () => {
    expect(runKey({ ...model, units: [{ ...model.units[0], params: { T: 310, P: 1e5 } }] }, "optimization", "O")).not.toBe(key);
    expect(runKey({ ...model, optimizations: [{ ...model.optimizations![0], objective: "2" }] }, "optimization", "O")).not.toBe(key);
    expect(runKey({ ...model, costing: {} }, "optimization", "O")).not.toBe(key);
  });

  it("select the runs that still match the flowsheet for the report", () => {
    const runs = emptyRuns();
    runs.sensitivity.A = { key: runKey(model, "sensitivity", "A"), result: { ok: true, id: "A", rows: [] } };
    runs.optimizations.O = { key: "old", result: { ok: true, id: "O" } };
    runs.dynamics = { key: runKey(model, "dynamics"), result: { ok: true, time: [0, 1 / 3], series: [], warnings: [] } };
    const r = currentRuns(runs, model);
    expect(r.sensitivity.map((s) => s.id)).toEqual(["A"]);
    expect(r.optimizations).toEqual([]);
    expect(r.dynamics?.time).toEqual([0, 0.333333]);
    expect(currentRuns(runs, model, false).dynamics).toBeNull();
    // A removed study is not reported.
    expect(currentRuns(runs, { ...model, sensitivity: [] }).sensitivity).toEqual([]);
  });

  it("leave the solve at the optimum out of the report", () => {
    const runs = emptyRuns();
    runs.optimizations.O = { key: runKey(model, "optimization", "O"), result: { ok: true, id: "O", result: {} as never } };
    expect(currentRuns(runs, model).optimizations[0]).toEqual({ ok: true, id: "O", result: undefined });
  });
});

describe("column side outlets", () => {
  const column = { params: [{ name: "side_draws" }] } as unknown as UnitType;
  const connected = [
    { outlet: "distillate", stream: "S2" },
    { outlet: "side1", stream: "S5" },
    { outlet: "side3", stream: "S7" },
  ];

  it("lists connected side outlets beyond the side draws", () => {
    expect(idleSideOutlets(column, { side_draws: [{}, {}] }, connected)).toEqual([{ outlet: "side3", stream: "S7" }]);
    expect(idleSideOutlets(column, {}, connected).map((x) => x.outlet)).toEqual(["side1", "side3"]);
    expect(idleSideOutlets(column, { side_draws: [{}, {}, {}] }, connected)).toEqual([]);
    expect(idleSideOutlets({ params: [] } as unknown as UnitType, {}, connected)).toEqual([]);
  });

  it("explains the renumbering", () => {
    expect(idleSideOutletText([{ outlet: "side3", stream: "S7" }], 2)).toContain(
      "side3 (S7) is connected, but the column has only 2 side draws",
    );
  });
});

describe("curve table row keys", () => {
  it("keep keys while the length matches and renew them otherwise", () => {
    let n = 0;
    const next = () => n++;
    const keys = rowKeys([], 3, next);
    expect(keys).toEqual([0, 1, 2]);
    expect(rowKeys(keys, 3, next)).toBe(keys);
    expect(rowKeys([0, 2], 2, next)).toEqual([0, 2]); // after removing the middle row
    expect(rowKeys(keys, 1, next)).toEqual([3]); // changed from outside (undo)
  });
});
