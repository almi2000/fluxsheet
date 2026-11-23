// A tubular reactor exchanges heat along its length with a utility or with a process stream.
// Only with a stream on its medium side (param `medium` = "stream") does it have the ports
// that stream flows through; otherwise the solver knows no such ports, so the editor hides
// them.
import type { Params, UnitType } from "./types";

export const MEDIUM_PORTS = new Set(["coolant_in", "coolant_out"]);

/** The ports of `ports` the unit has with these parameters. */
export function activePorts<T extends { name: string }>(schema: UnitType | undefined, params: Params, ports: T[]): T[] {
  const medium = schema?.params.find((p) => p.name === "medium");
  if (!medium) return ports;
  const value = params.medium ?? medium.default;
  return value === "stream" ? ports : ports.filter((p) => !MEDIUM_PORTS.has(p.name));
}
