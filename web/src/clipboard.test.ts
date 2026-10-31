import { describe, expect, it } from "vitest";
import { extractFragment, freshId, parseFragment, remapFragment } from "./clipboard";
import type { FlowsheetModel } from "./types";

const model: FlowsheetModel = {
  format: "fluxsheet/1",
  name: "Test",
  components: ["water", "ethanol"],
  property_method: "auto",
  units: [
    { id: "FEED", type: "Feed", ui: { x: 0, y: 100 } },
    { id: "MIX1", type: "Mixer", ui: { x: 200, y: 100 } },
    { id: "HX1", type: "Heater", ui: { x: 400, y: 120 }, params: { T: "80 degC" } },
    { id: "PROD1", type: "Product", ui: { x: 600, y: 120 } },
  ],
  streams: [
    { id: "S1", from: "FEED.out", to: "MIX1.in" },
    { id: "S2", from: "MIX1.out", to: "HX1.in" },
    { id: "S3", from: "HX1.out", to: "PROD1.in" },
  ],
};

describe("clipboard", () => {
  it("keeps only the streams between selected units", () => {
    const f = extractFragment(model, new Set(["MIX1", "HX1"]));
    expect(f.units.map((u) => u.id)).toEqual(["MIX1", "HX1"]);
    expect(f.streams.map((s) => s.id)).toEqual(["S2"]);
    expect(f.components).toEqual(model.components);
  });

  it("finds the next free id with the same stem", () => {
    const taken = new Set(["MIX1", "MIX2", "CHILLER", "S1"]);
    expect(freshId("MIX1", taken)).toBe("MIX3");
    expect(freshId("CHILLER", taken)).toBe("CHILLER1");
    expect(freshId("S", taken)).toBe("S2");
    expect(freshId("42", taken)).toBe("U1");
  });

  it("remaps ids, rewires streams and moves the fragment", () => {
    const f = extractFragment(model, new Set(["MIX1", "HX1"]));
    const taken = new Set(model.units.map((u) => u.id));
    const { units, streams, idMap } = remapFragment(f, taken, new Set(["S1", "S2", "S3"]), { x: 1000, y: 50 });
    expect(units.map((u) => u.id)).toEqual(["MIX2", "HX2"]);
    expect(units.map((u) => u.ui)).toEqual([{ x: 1000, y: 50 }, { x: 1200, y: 70 }]);
    expect(units[1].params).toEqual({ T: "80 degC" });
    expect(streams).toEqual([{ id: "S4", from: "MIX2.out", to: "HX2.in" }]);
    expect(idMap.get("HX1")).toBe("HX2");
    // The source fragment is untouched.
    expect(f.units[0].id).toBe("MIX1");
  });

  it("moves a column's feed stages to the pasted stream ids", () => {
    const f: FlowsheetModel = {
      ...model,
      units: [
        { id: "F1", type: "Feed" },
        { id: "F2", type: "Feed" },
        { id: "COL", type: "DistillationColumn", params: { feed_stages: { S1: 4, S2: 9 } } },
      ],
      streams: [
        { id: "S1", from: "F1.out", to: "COL.in" },
        { id: "S2", from: "F2.out", to: "COL.in" },
      ],
    };
    const { units, streams } = remapFragment(f, new Set(["F1", "F2", "COL"]), new Set(["S1", "S2"]), { x: 0, y: 0 });
    expect(streams.map((s) => s.id)).toEqual(["S3", "S4"]);
    expect(units[2].params).toEqual({ feed_stages: { S3: 4, S4: 9 } });
  });

  it("gives each pasted unit a distinct id even with the same stem", () => {
    const f: FlowsheetModel = { ...model, units: [{ id: "A1", type: "Mixer" }, { id: "A2", type: "Mixer" }], streams: [] };
    const { units } = remapFragment(f, new Set(["A1", "A2"]), new Set(), { x: 0, y: 0 });
    expect(units.map((u) => u.id)).toEqual(["A3", "A4"]);
  });

  it("parses only Fluxsheet JSON", () => {
    expect(parseFragment(JSON.stringify(model))?.units).toHaveLength(4);
    expect(parseFragment("hello")).toBeNull();
    expect(parseFragment('{"format": "other", "units": []}')).toBeNull();
    expect(parseFragment('{"format": "fluxsheet/1", "units": [{"id": 3}]}')).toBeNull();
    expect(parseFragment('{"format": "fluxsheet/1", "units": []}')?.streams).toEqual([]);
  });
});
