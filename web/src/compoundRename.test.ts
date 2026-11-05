import { describe, expect, it } from "vitest";
import { componentsUsed, renameComponent } from "./compoundRename";
import type { FlowsheetModel, ParamSchema, UnitType } from "./types";

const param = (name: string, quantity: string | null, fields?: ParamSchema[]): ParamSchema => ({
  name, label: name, quantity, default: null, help: "", min: null, max: null, options: null, fields,
});
const type = (t: string, params: ParamSchema[], dynamic?: ParamSchema[]): UnitType => ({
  type: t, label: t, category: "", description: "", inlets: [], outlets: [], params,
  exclusive_groups: [], required_groups: [], result_fields: [],
  ...(dynamic && { dynamic: { params: dynamic, result_fields: [] } }),
});
const schemas: Record<string, UnitType> = {
  Feed: type("Feed", [param("flows", "composition"), param("T", "temperature")]),
  Column: type("Column", [
    param("light_key", "component"),
    param("heavy_key", "component"),
    param("side_draws", "records", [param("stage", "integer"), param("component", "component")]),
  ]),
  Reactor: type("Reactor", [param("keys", "string_list"), param("reactions", "string_list"), param("conversions", "string_list")]),
  Splitter: type("Splitter", [param("split", "component_fractions")], [param("x0", "component_fractions")]),
};

const model: FlowsheetModel = {
  format: "fluxsheet/1",
  name: "Rename",
  components: ["water", "my-ester", "ethanol"],
  henry_components: ["my-ester"],
  property_method: "NRTL",
  custom_compounds: [{ name: "my-ester", formula: "C5H10O2" }],
  units: [
    { id: "F1", type: "Feed", params: { flows: { water: 1, "my-ester": 2 }, T: "25 degC" } },
    { id: "C1", type: "Column", params: { light_key: "my-ester", heavy_key: "water", side_draws: [{ stage: 3, component: "my-ester" }] } },
    { id: "R1", type: "Reactor", params: { keys: ["my-ester"], reactions: ["my-ester"], conversions: ["0.5"] } },
    { id: "SP", type: "Splitter", params: { split: { "my-ester": 0.3 } } },
    { id: "X", type: "Unknown", params: { flows: { "my-ester": 1 } } },
  ],
  streams: [],
  reactions: [
    { id: "my-ester", stoichiometry: { ethanol: -1, "my-ester": 1 }, kinetics: { k0: 1, orders: { "my-ester": 1 }, reverse: { k0: 2, orders: { "my-ester": 2 } } } },
  ],
  design_specs: [{ id: "D", target: { stream: "S1", property: "mole_fraction", component: "my-ester" }, value: 0.9, vary: { unit: "F1", param: "flows", component: "my-ester", lower: 0, upper: 5 } }],
  sensitivity: [{ id: "ST", vary: [{ unit: "F1", param: "flows", component: "my-ester", lower: 0, upper: 1 }], outputs: [{ stream: "S1", property: "component_flow", component: "my-ester" }] }],
  calculators: [{ id: "CALC", variables: { x: { stream: "S1", property: "mole_fraction", component: "my-ester" } }, expression: "x * 2" }],
  optimizations: [{ id: "O", variables: { y: { stream: "S1", property: "component_flow", component: "my-ester" } }, objective: "y", vary: [{ unit: "F1", param: "flows", component: "my-ester", lower: 0, upper: 1 }] }],
  dynamics: {
    t_end: 10,
    units: { SP: { x0: { "my-ester": 0.5 } } },
    events: [{ time: 1, unit: "F1", param: "flows", component: "my-ester", value: 2 }],
  },
};

describe("renameComponent", () => {
  const r = renameComponent(model, "my-ester", "ester-B", schemas);

  it("renames the component list, Henry list and definition", () => {
    expect(r.components).toEqual(["water", "ester-B", "ethanol"]);
    expect(r.henry_components).toEqual(["ester-B"]);
    expect(r.custom_compounds).toEqual([{ name: "ester-B", formula: "C5H10O2" }]);
  });

  it("renames unit parameters by their kind and keeps entry order", () => {
    const p = (id: string) => r.units.find((u) => u.id === id)!.params!;
    expect(Object.keys(p("F1").flows as object)).toEqual(["water", "ester-B"]);
    expect(p("F1").T).toBe("25 degC");
    expect(p("C1")).toMatchObject({ light_key: "ester-B", heavy_key: "water", side_draws: [{ stage: 3, component: "ester-B" }] });
    expect(p("R1")).toEqual({ keys: ["ester-B"], reactions: ["my-ester"], conversions: ["0.5"] }); // reaction ids stay
    expect(p("SP").split).toEqual({ "ester-B": 0.3 });
    expect(p("X").flows).toEqual({ "my-ester": 1 }); // unknown unit type: left alone
  });

  it("renames reactions but not their ids", () => {
    const rx = r.reactions![0];
    expect(rx.id).toBe("my-ester");
    expect(rx.stoichiometry).toEqual({ ethanol: -1, "ester-B": 1 });
    expect(rx.kinetics?.orders).toEqual({ "ester-B": 1 });
    expect(rx.kinetics?.reverse?.orders).toEqual({ "ester-B": 2 });
  });

  it("renames references of specs, studies, calculators, optimizations and dynamics", () => {
    expect(r.design_specs![0].target.component).toBe("ester-B");
    expect(r.design_specs![0].vary.component).toBe("ester-B");
    expect(r.sensitivity![0].vary[0].component).toBe("ester-B");
    expect(r.sensitivity![0].outputs[0].component).toBe("ester-B");
    expect(JSON.stringify(r.calculators)).not.toContain("my-ester");
    expect((r.calculators![0] as { expression: string }).expression).toBe("x * 2");
    expect((r.optimizations![0].variables.y as { component: string }).component).toBe("ester-B");
    expect(r.optimizations![0].vary[0].component).toBe("ester-B");
    expect(r.dynamics!.events![0].component).toBe("ester-B");
    expect(r.dynamics!.units.SP.x0).toEqual({ "ester-B": 0.5 });
  });

  it("does not change the input model", () => {
    expect(model.units[0].params!.flows).toEqual({ water: 1, "my-ester": 2 });
  });

  it("finds the components the units use", () => {
    expect([...componentsUsed(model, schemas)].sort()).toEqual(["my-ester", "water"]);
  });
});
