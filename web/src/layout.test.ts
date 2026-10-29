import { describe, expect, it } from "vitest";
import { assignLayers, findRecycles, layeredLayout, type LayoutEdge, type LayoutNode } from "./layout";

const node = (id: string, role: LayoutNode["role"] = "unit", outlets = ["out"], inlets = ["in"]): LayoutNode => ({
  id, role, width: 80, height: 70, inlets: role === "feed" ? [] : inlets, outlets: role === "product" ? [] : outlets,
});
const edge = (source: string, target: string, sourceHandle = "out", targetHandle = "in"): LayoutEdge => ({
  source, target, sourceHandle, targetHandle,
});

// Gas recycle: FEED -> MIX -> CHILL -> SEP -vapor-> SPLIT -> COMP -> MIX, plus two products.
const nodes = [
  node("FUEL", "product"),
  node("COMP"),
  node("SPLIT", "unit", ["out1", "out2"]),
  node("SEP", "unit", ["vapor", "liquid"]),
  node("CHILL"),
  node("MIX"),
  node("NGL", "product"),
  node("FEED", "feed"),
];
const edges = [
  edge("FEED", "MIX"),
  edge("MIX", "CHILL"),
  edge("CHILL", "SEP"),
  edge("SEP", "SPLIT", "vapor"),
  edge("SEP", "NGL", "liquid"),
  edge("SPLIT", "COMP", "out1"),
  edge("COMP", "MIX"),
  edge("SPLIT", "FUEL", "out2"),
];

describe("layout", () => {
  it("identifies the recycle stream from the feed", () => {
    const back = [...findRecycles(nodes, edges)];
    expect(back).toEqual([edges[6]]); // COMP -> MIX
  });

  it("puts units left to right in flow order and products last", () => {
    const back = findRecycles(nodes, edges);
    const layer = assignLayers(nodes, edges.filter((e) => !back.has(e)));
    const order = ["FEED", "MIX", "CHILL", "SEP", "SPLIT", "COMP"];
    order.slice(1).forEach((id, i) => expect(layer.get(id)!).toBeGreaterThan(layer.get(order[i])!));
    expect(layer.get("NGL")).toBe(layer.get("FUEL"));
    expect(layer.get("NGL")).toBe(Math.max(...layer.values()));
  });

  it("moves a side feed next to the unit it feeds", () => {
    const ns = [node("F1", "feed"), node("A"), node("B", "unit", ["out"], ["in"]), node("F2", "feed"), node("P", "product")];
    const es = [edge("F1", "A"), edge("A", "B"), edge("F2", "B"), edge("B", "P")];
    const layer = assignLayers(ns, es);
    expect(layer.get("F2")).toBe(layer.get("B")! - 1);
  });

  it("produces non-overlapping, grid-snapped positions", () => {
    const pos = layeredLayout(nodes, edges);
    expect(pos.size).toBe(nodes.length);
    const boxes = [...pos.entries()];
    for (const [, p] of boxes) {
      expect(p.x % 10).toBe(0);
      expect(p.y % 10).toBe(0);
    }
    for (let i = 0; i < boxes.length; i++)
      for (let j = i + 1; j < boxes.length; j++) {
        const [a, b] = [boxes[i][1], boxes[j][1]];
        const overlap = Math.abs(a.x - b.x) < 80 && Math.abs(a.y - b.y) < 70;
        expect(overlap, `${boxes[i][0]} overlaps ${boxes[j][0]}`).toBe(false);
      }
    // Flash vapor goes up, liquid down.
    expect(pos.get("SPLIT")!.y).toBeLessThan(pos.get("NGL")!.y);
    expect(Math.min(...[...pos.values()].map((p) => p.y))).toBe(0);
  });

  it("handles empty graphs and pure cycles", () => {
    expect(layeredLayout([], []).size).toBe(0);
    const ring = [node("A"), node("B")];
    const pos = layeredLayout(ring, [edge("A", "B"), edge("B", "A")]);
    expect(pos.get("A")!.x).toBeLessThan(pos.get("B")!.x);
  });
});
