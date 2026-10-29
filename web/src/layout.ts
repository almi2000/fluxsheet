/**
 * Left-to-right layered (Sugiyama-style) layout for flowsheets.
 *
 * 1. Recycle streams are found by a depth-first search from the feeds and set aside, which
 *    leaves the flow direction as a DAG.
 * 2. Units get a column by longest path; feeds move right next to the unit they feed, and
 *    products all go into the last column.
 * 3. A few barycenter sweeps order the units within each column to reduce crossings, taking
 *    port positions into account (vapor outlets above liquid outlets, ...).
 * 4. Rows: each unit is pulled toward the ports it connects to, alternating forward and
 *    backward passes, while keeping a minimum gap and the column order.
 */

export type Point = { x: number; y: number };

export interface LayoutNode {
  id: string;
  width: number;
  height: number;
  role: "feed" | "product" | "unit";
  inlets: string[];
  outlets: string[];
}

export interface LayoutEdge {
  source: string;
  sourceHandle?: string | null;
  target: string;
  targetHandle?: string | null;
}

export interface LayoutOptions {
  gapX?: number;
  gapY?: number;
  grid?: number;
}

/** Edges that close a cycle when walking downstream from the feeds (recycle streams). */
export function findRecycles(nodes: LayoutNode[], edges: LayoutEdge[]): Set<LayoutEdge> {
  const byId = new Map(nodes.map((n) => [n.id, n]));
  const out = new Map<string, LayoutEdge[]>(nodes.map((n) => [n.id, []]));
  const indeg = new Map<string, number>(nodes.map((n) => [n.id, 0]));
  for (const e of edges) {
    if (!byId.has(e.source) || !byId.has(e.target)) continue;
    out.get(e.source)!.push(e);
    indeg.set(e.target, indeg.get(e.target)! + 1);
  }
  // Follow outlets in port order so the main product path is explored first.
  for (const [id, list] of out) {
    const ports = byId.get(id)!.outlets;
    list.sort((a, b) => ports.indexOf(a.sourceHandle ?? "") - ports.indexOf(b.sourceHandle ?? ""));
  }

  const roots = [
    ...nodes.filter((n) => n.role === "feed"),
    ...nodes.filter((n) => n.role !== "feed" && indeg.get(n.id) === 0),
    ...nodes,
  ];
  const state = new Map<string, 0 | 1 | 2>(); // 0 new, 1 on the stack, 2 done
  const back = new Set<LayoutEdge>();
  for (const root of roots) {
    if (state.get(root.id)) continue;
    const stack: { id: string; i: number }[] = [{ id: root.id, i: 0 }];
    state.set(root.id, 1);
    while (stack.length) {
      const top = stack[stack.length - 1];
      const e = out.get(top.id)![top.i++];
      if (!e) {
        state.set(top.id, 2);
        stack.pop();
        continue;
      }
      const s = state.get(e.target) ?? 0;
      if (s === 1) back.add(e);
      else if (s === 0) {
        state.set(e.target, 1);
        stack.push({ id: e.target, i: 0 });
      }
    }
  }
  return back;
}

/** Column index per node, flow left to right, products in the last column. */
export function assignLayers(nodes: LayoutNode[], forward: LayoutEdge[]): Map<string, number> {
  const preds = new Map<string, string[]>(nodes.map((n) => [n.id, []]));
  const succs = new Map<string, string[]>(nodes.map((n) => [n.id, []]));
  for (const e of forward) {
    preds.get(e.target)?.push(e.source);
    succs.get(e.source)?.push(e.target);
  }
  // Kahn's algorithm gives a topological order of the (acyclic) forward graph.
  const indeg = new Map(nodes.map((n) => [n.id, preds.get(n.id)!.length]));
  const queue = nodes.filter((n) => indeg.get(n.id) === 0).map((n) => n.id);
  const order: string[] = [];
  while (queue.length) {
    const id = queue.shift()!;
    order.push(id);
    for (const t of succs.get(id)!) {
      indeg.set(t, indeg.get(t)! - 1);
      if (indeg.get(t) === 0) queue.push(t);
    }
  }
  const layer = new Map<string, number>();
  for (const id of order) layer.set(id, Math.max(0, ...preds.get(id)!.map((p) => layer.get(p)! + 1)));

  const role = new Map(nodes.map((n) => [n.id, n.role]));
  // Sources (feeds) sit just left of the first unit they feed instead of in column 0.
  for (const id of [...order].reverse()) {
    if (preds.get(id)!.length > 0 || role.get(id) === "product") continue;
    const next = succs.get(id)!.map((t) => layer.get(t)!);
    if (next.length) layer.set(id, Math.max(0, Math.min(...next) - 1));
  }
  const last = Math.max(-1, ...nodes.filter((n) => n.role !== "product").map((n) => layer.get(n.id)!));
  for (const n of nodes) if (n.role === "product") layer.set(n.id, last + 1);
  return layer;
}

/** Position of an edge's port along a node's left/right side, -0.5 (top) .. 0.5 (bottom). */
function portOffset(ports: string[], handle: string | null | undefined): number {
  const i = Math.max(0, ports.indexOf(handle ?? ""));
  return (i + 1) / (ports.length + 1) - 0.5;
}

export function layeredLayout(
  nodes: LayoutNode[],
  edges: LayoutEdge[],
  { gapX = 110, gapY = 50, grid = 10 }: LayoutOptions = {},
): Map<string, Point> {
  const byId = new Map(nodes.map((n) => [n.id, n]));
  const valid = edges.filter((e) => byId.has(e.source) && byId.has(e.target) && e.source !== e.target);
  const recycles = findRecycles(nodes, valid);
  const forward = valid.filter((e) => !recycles.has(e));
  const layer = assignLayers(nodes, forward);

  const nLayers = Math.max(0, ...layer.values()) + 1;
  const columns: string[][] = Array.from({ length: nLayers }, () => []);
  for (const n of nodes) columns[layer.get(n.id)!].push(n.id);

  // Neighbors with the port offsets on both ends, for ordering and row placement.
  type Link = { other: string; mine: number; theirs: number };
  const ins = new Map<string, Link[]>(nodes.map((n) => [n.id, []]));
  const outs = new Map<string, Link[]>(nodes.map((n) => [n.id, []]));
  for (const e of forward) {
    const s = byId.get(e.source)!;
    const t = byId.get(e.target)!;
    const so = portOffset(s.outlets, e.sourceHandle);
    const to = portOffset(t.inlets, e.targetHandle);
    outs.get(s.id)!.push({ other: t.id, mine: so, theirs: to });
    ins.get(t.id)!.push({ other: s.id, mine: to, theirs: so });
  }

  // ---- order within columns (barycenter sweeps) ----
  const rank = new Map<string, number>();
  const setRanks = () => columns.forEach((col) => col.forEach((id, i) => rank.set(id, i)));
  setRanks();
  const sortColumn = (col: string[], links: Map<string, Link[]>) => {
    const key = new Map(
      col.map((id) => {
        const ls = links.get(id)!;
        if (!ls.length) return [id, rank.get(id)!] as const;
        return [id, ls.reduce((sum, l) => sum + rank.get(l.other)! + l.theirs, 0) / ls.length] as const;
      }),
    );
    col.sort((a, b) => key.get(a)! - key.get(b)! || rank.get(a)! - rank.get(b)!);
    col.forEach((id, i) => rank.set(id, i));
  };
  for (let sweep = 0; sweep < 4; sweep++) {
    for (let l = 1; l < nLayers; l++) sortColumn(columns[l], ins);
    for (let l = nLayers - 2; l >= 0; l--) sortColumn(columns[l], outs);
  }

  // ---- rows ----
  const cy = new Map<string, number>(); // vertical center of each node
  for (const col of columns) {
    const total = col.reduce((s, id) => s + byId.get(id)!.height, 0) + gapY * (col.length - 1);
    let y = -total / 2;
    for (const id of col) {
      const h = byId.get(id)!.height;
      cy.set(id, y + h / 2);
      y += h + gapY;
    }
  }
  const portY = (id: string, offset: number) => cy.get(id)! + offset * byId.get(id)!.height;
  /** Where a node's center would line its ports up with its neighbors' ports. */
  const wantY = (id: string, links: Map<string, Link[]>) => {
    const ls = links.get(id)!;
    if (!ls.length) return cy.get(id)!;
    const h = byId.get(id)!.height;
    return ls.reduce((s, l) => s + portY(l.other, l.theirs) - l.mine * h, 0) / ls.length;
  };
  /** Moves each node of a column toward its neighbors' ports, keeping order and gaps. */
  const placeColumn = (col: string[], links: Map<string, Link[]>, reorder = false) => {
    if (!col.length) return;
    // Once rows exist, real heights order a column better than ranks in columns of other sizes.
    if (reorder) {
      const key = new Map(col.map((id) => [id, wantY(id, links)]));
      col.sort((a, b) => key.get(a)! - key.get(b)!);
    }
    const want = col.map((id) => wantY(id, links));
    // Pack top to bottom, then shift the column so it deviates least from what it wants.
    const got: number[] = [];
    col.forEach((id, i) => {
      const h = byId.get(id)!.height;
      const min = i === 0 ? -Infinity : got[i - 1] + byId.get(col[i - 1])!.height / 2 + gapY + h / 2;
      got.push(Math.max(want[i], min));
    });
    const shift = want.reduce((s, w, i) => s + w - got[i], 0) / col.length;
    col.forEach((id, i) => cy.set(id, got[i] + shift));
  };
  for (let pass = 0; pass < 6; pass++) {
    for (let l = 1; l < nLayers; l++) placeColumn(columns[l], ins, pass < 2);
    for (let l = nLayers - 2; l >= 0; l--) placeColumn(columns[l], outs);
  }
  for (let l = 1; l < nLayers; l++) placeColumn(columns[l], ins);

  // ---- columns ----
  const snap = (v: number) => Math.round(v / grid) * grid;
  const result = new Map<string, Point>();
  let x = 0;
  for (const col of columns) {
    const width = Math.max(0, ...col.map((id) => byId.get(id)!.width));
    for (const id of col) {
      const n = byId.get(id)!;
      result.set(id, { x: x + (width - n.width) / 2, y: cy.get(id)! - n.height / 2 });
    }
    x += width + gapX;
  }
  const minY = Math.min(...[...result.values()].map((p) => p.y));
  for (const [id, p] of result) result.set(id, { x: snap(p.x), y: snap(p.y - minY) });
  return result;
}
