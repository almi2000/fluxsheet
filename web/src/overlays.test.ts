import { afterEach, describe, expect, it, vi } from "vitest";
import { DEFAULT_OVERLAY, loadOverlay, saveOverlay, streamLabel, unitWarnings } from "./overlays";
import type { SolveResult, StreamResult } from "./types";

const stream = (over: Partial<StreamResult> = {}): StreamResult => ({
  T: 353.15, P: 2e5, vapor_fraction: 0, phase: "liquid", molar_flow: 100 / 3.6, mass_flow: 1,
  enthalpy: 0, enthalpy_flow: 0, component_flows: {}, mole_fractions: { water: 0.92, ethanol: 0.08 },
  ...over,
});

const solved = (over: Partial<SolveResult> = {}): SolveResult => ({
  ok: true, streams: {}, units: {}, property_method: null, property_method_reason: null,
  tear_streams: [], iterations: 1, converged: true, warnings: [], ...over,
} as SolveResult);

describe("streamLabel", () => {
  it("shows only the id without results or when asked", () => {
    expect(streamLabel("S1", undefined, "all")).toBe("S1");
    expect(streamLabel("S1", stream(), "id")).toBe("S1");
  });

  it("adds temperature, pressure, flow and the main component step by step", () => {
    expect(streamLabel("S1", stream(), "T")).toBe("S1 · 80.0 °C");
    expect(streamLabel("S1", stream(), "TP")).toBe("S1 · 80.0 °C · 2 bar");
    expect(streamLabel("S1", stream(), "TPF")).toBe("S1 · 80.0 °C · 2 bar · 100 kmol/h");
    expect(streamLabel("S1", stream(), "all")).toBe("S1 · 80.0 °C · 2 bar · 100 kmol/h · 92 % water");
  });

  it("leaves out the composition of an empty stream and rounds a pure one to 100 %", () => {
    expect(streamLabel("S1", stream({ molar_flow: 0 }), "all")).toBe("S1 · 80.0 °C · 2 bar · 0 kmol/h");
    expect(streamLabel("S1", stream({ mole_fractions: { water: 0.99999 } }), "all")).toMatch(/100 % water$/);
  });
});

describe("unitWarnings", () => {
  it("picks the notes of this unit only, without the id prefix", () => {
    const r = solved({ warnings: ["PU: the flow is outside the pump curve", "PU2: other", "a global warning"] });
    expect(unitWarnings("PU", r)).toEqual(["the flow is outside the pump curve"]);
  });

  it("adds balance checks of the unit that did not close", () => {
    const r = solved({
      balance_checks: [
        { scope: "R1", quantity: "mass", component: null, residual: 0.0123, unit: "kg/s", tolerance: 1e-6, status: "warning" },
        { scope: "R1", quantity: "energy", component: null, residual: 0, unit: "W", tolerance: 1, status: "ok" },
        { scope: "R2", quantity: "mass", component: null, residual: 1, unit: "kg/s", tolerance: 1e-6, status: "warning" },
      ],
    });
    expect(unitWarnings("R1", r)).toEqual(["mass balance does not close (residual 0.0123 kg/s)"]);
  });

  it("has nothing to say about a failed run or no run", () => {
    expect(unitWarnings("PU", null)).toEqual([]);
    expect(unitWarnings("PU", solved({ ok: false, warnings: ["PU: x"] }))).toEqual([]);
  });
});

describe("overlay choice", () => {
  afterEach(() => vi.unstubAllGlobals());

  it("falls back to the default when storage throws or holds junk", () => {
    vi.stubGlobal("window", { localStorage: { getItem: () => { throw new Error("blocked"); }, setItem: () => { throw new Error("blocked"); } } });
    expect(loadOverlay()).toBe(DEFAULT_OVERLAY);
    expect(() => saveOverlay("TP")).not.toThrow();
    vi.stubGlobal("window", { localStorage: { getItem: () => "nonsense" } });
    expect(loadOverlay()).toBe(DEFAULT_OVERLAY);
  });

  it("remembers the last choice", () => {
    const store = new Map<string, string>();
    vi.stubGlobal("window", { localStorage: { getItem: (k: string) => store.get(k) ?? null, setItem: (k: string, v: string) => store.set(k, v) } });
    saveOverlay("TPF");
    expect(loadOverlay()).toBe("TPF");
  });
});
