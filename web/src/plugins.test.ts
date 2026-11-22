import { describe, expect, it } from "vitest";
import { methodLabel, pluginProblems } from "./plugins";
import type { PluginInfo } from "./types";

const plugin = (over: Partial<PluginInfo>): PluginInfo => ({
  group: "units",
  name: "example_units",
  value: "example.units",
  distribution: "example",
  version: "1.0",
  provides: ["DewPointCooler"],
  error: null,
  ...over,
});

describe("plug-ins in the editor", () => {
  it("marks plug-in property methods", () => {
    expect(methodLabel({ key: "nrtl", label: "NRTL", description: "" })).toBe("NRTL");
    expect(methodLabel({ key: "rs", label: "Regular solution", description: "", plugin: "x" })).toBe(
      "Regular solution (plug-in)",
    );
  });

  it("lists only the plug-ins that failed, by distribution or entry point", () => {
    expect(
      pluginProblems([
        plugin({}),
        plugin({ error: "ValueError: unit type already defined: Heater" }),
        plugin({ distribution: null, name: "local", value: "mod:X", error: "ImportError: no" }),
      ]),
    ).toEqual([
      "example (example_units): ValueError: unit type already defined: Heater",
      "mod:X (local): ImportError: no",
    ]);
  });
});
