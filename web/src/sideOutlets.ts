// Column side draws leave through the outlets side1, side2, ... in the order of the side draw
// table, so deleting a draw in the middle renumbers the later ones. The editor warns about
// connected side outlets that no draw feeds (the column leaves them empty).
import type { Params, UnitType } from "./types";

export interface IdleSideOutlet {
  outlet: string;
  stream: string;
}

/** Connected `sideN` outlets of a unit beyond its number of side draws. */
export function idleSideOutlets(
  schema: UnitType | undefined,
  params: Params,
  connected: { outlet: string; stream: string }[],
): IdleSideOutlet[] {
  if (!schema?.params.some((p) => p.name === "side_draws")) return [];
  const draws = Array.isArray(params.side_draws) ? params.side_draws.length : 0;
  return connected.filter(({ outlet }) => {
    const m = /^side(\d+)$/.exec(outlet);
    return !!m && Number(m[1]) > draws;
  });
}

export function idleSideOutletText(idle: IdleSideOutlet[], draws: number): string {
  const list = idle.map((x) => `${x.outlet} (${x.stream})`).join(", ");
  const has = draws === 0 ? "no side draws" : draws === 1 ? "only 1 side draw" : `only ${draws} side draws`;
  return (
    `${idle.length === 1 ? "Outlet" : "Outlets"} ${list} ${idle.length === 1 ? "is" : "are"} connected, but the column has ${has}, ` +
    "so no product leaves there. Side draws use the outlets in table order: deleting one renumbers the later ones. " +
    "Reconnect the stream or add the draw."
  );
}
