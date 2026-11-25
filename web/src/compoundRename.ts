// Renaming a component everywhere the flowsheet refers to it by name.
import type { FlowsheetModel, ParamSchema, Params, Reaction, UnitType } from "./types";

/** {old: 1, b: 2} -> {new: 1, b: 2}, keeping the order of the entries. */
function renameKey<T>(map: Record<string, T>, from: string, to: string): Record<string, T> {
  if (!map || typeof map !== "object" || Array.isArray(map) || !(from in map)) return map;
  return Object.fromEntries(Object.entries(map).map(([k, v]) => [k === from ? to : k, v]));
}

const renameIn = (list: string[] | undefined, from: string, to: string) => list?.map((c) => (c === from ? to : c));

/** "string_list" parameters (units/base.py) that hold something other than component names. */
const NOT_COMPONENTS = new Set(["reactions", "conversions"]);

/** Parameter kinds (units/base.py) that hold component names. */
function renameParam(p: ParamSchema, value: unknown, from: string, to: string): unknown {
  switch (p.quantity) {
    case "component":
      return value === from ? to : value;
    case "composition":
    case "component_fractions":
      return renameKey(value as Record<string, unknown>, from, to);
    case "string_list": // component lists (key components, inerts), but not reaction ids or conversions
      return !NOT_COMPONENTS.has(p.name) && Array.isArray(value) ? renameIn(value as string[], from, to) : value;
    case "records":
      return Array.isArray(value) && p.fields
        ? value.map((r) => (r && typeof r === "object" ? renameParams(p.fields!, r as Params, from, to) : r))
        : value;
    default:
      return value;
  }
}

function renameParams(schema: ParamSchema[], params: Params, from: string, to: string): Params {
  const out: Params = { ...params };
  for (const p of schema) if (p.name in out) out[p.name] = renameParam(p, out[p.name], from, to);
  return out;
}

/** Any `component` field (variable references, varied parameters, events, ...) in `x`. */
function renameRefs<T>(x: T, from: string, to: string): T {
  if (Array.isArray(x)) return x.map((v) => renameRefs(v, from, to)) as T;
  if (!x || typeof x !== "object") return x;
  return Object.fromEntries(
    Object.entries(x).map(([k, v]) => [k, k === "component" && v === from ? to : renameRefs(v, from, to)]),
  ) as T;
}

function renameReaction(r: Reaction, from: string, to: string): Reaction {
  const k = r.kinetics;
  return {
    ...r,
    stoichiometry: renameKey(r.stoichiometry, from, to),
    ...(k && {
      kinetics: {
        ...k,
        ...(k.orders && { orders: renameKey(k.orders, from, to) }),
        ...(k.reverse && { reverse: { ...k.reverse, ...(k.reverse.orders && { orders: renameKey(k.reverse.orders, from, to) }) } }),
      },
    }),
  };
}

/**
 * The flowsheet with component `from` called `to`: the component list, Henry components,
 * custom compound definitions, unit parameters (by their kind in `schemas`), reactions, and
 * the references of design specs, studies, calculators, dynamics and optimizations.
 * Expressions name their variables, not components, so they stay as they are.
 */
export function renameComponent(model: FlowsheetModel, from: string, to: string, schemas: Record<string, UnitType>): FlowsheetModel {
  if (from === to) return model;
  const dyn = model.dynamics;
  return {
    ...model,
    components: renameIn(model.components, from, to)!,
    henry_components: renameIn(model.henry_components, from, to),
    custom_compounds: model.custom_compounds?.map((c) => (c.name === from ? { ...c, name: to } : c)),
    property_data: model.property_data && {
      ...model.property_data,
      binary_parameters: model.property_data.binary_parameters?.map((p) => ({
        ...p, components: [p.components[0] === from ? to : p.components[0], p.components[1] === from ? to : p.components[1]],
      })),
      unifac_groups: model.property_data.unifac_groups?.map((g) => g.component === from ? { ...g, component: to } : g),
    },
    units: model.units.map((u) => {
      const schema = schemas[u.type]?.params;
      return schema && u.params ? { ...u, params: renameParams(schema, u.params, from, to) } : u;
    }),
    reactions: model.reactions?.map((r) => renameReaction(r, from, to)),
    design_specs: renameRefs(model.design_specs, from, to),
    sensitivity: renameRefs(model.sensitivity, from, to),
    calculators: renameRefs(model.calculators, from, to),
    optimizations: renameRefs(model.optimizations, from, to),
    dynamics: dyn && {
      ...renameRefs(dyn, from, to),
      // A unit's dynamic settings (e.g. an initial composition), by the kinds of its dynamic model.
      units: Object.fromEntries(
        Object.entries(dyn.units ?? {}).map(([id, params]) => {
          const type = model.units.find((u) => u.id === id)?.type;
          const schema = type ? schemas[type]?.dynamic?.params : undefined;
          return [id, schema ? renameParams(schema, params, from, to) : params];
        }),
      ),
    },
  };
}

/**
 * The components the parameters of `units` refer to, found the way `renameComponent` finds
 * them: a component renamed to a fresh name shows up as that name.
 */
export function componentsUsed(model: FlowsheetModel, schemas: Record<string, UnitType>): Set<string> {
  const used = new Set<string>();
  const marker = "\u0000used\u0000";
  for (const c of model.components) {
    const units = renameComponent({ ...model, reactions: undefined, dynamics: undefined }, c, marker, schemas).units;
    if (JSON.stringify(units) !== JSON.stringify(model.units)) used.add(c);
  }
  return used;
}
