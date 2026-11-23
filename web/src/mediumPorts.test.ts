import { describe, expect, it } from "vitest";
import { activePorts } from "./mediumPorts";
import type { UnitType } from "./types";

const port = (name: string) => ({ name, label: name, multiple: false, separate: false });
const pfr = {
  type: "PFR",
  params: [{ name: "medium", label: "Medium side", quantity: "choice", default: "utility", help: "", min: null, max: null, options: ["utility", "stream"] }],
  inlets: [port("in"), port("coolant_in")],
  outlets: [port("out"), port("coolant_out")],
} as unknown as UnitType;
const heater = { type: "Heater", params: [], inlets: [port("in")], outlets: [port("out")] } as unknown as UnitType;

describe("medium ports", () => {
  it("hides a tubular reactor's medium ports unless a stream flows there", () => {
    expect(activePorts(pfr, {}, pfr.inlets).map((p) => p.name)).toEqual(["in"]);
    expect(activePorts(pfr, { medium: "utility" }, pfr.outlets).map((p) => p.name)).toEqual(["out"]);
    expect(activePorts(pfr, { medium: "stream" }, pfr.inlets).map((p) => p.name)).toEqual(["in", "coolant_in"]);
  });

  it("leaves other units alone", () => {
    expect(activePorts(heater, {}, heater.inlets)).toEqual(heater.inlets);
    expect(activePorts(undefined, {}, heater.outlets)).toEqual(heater.outlets);
  });
});
