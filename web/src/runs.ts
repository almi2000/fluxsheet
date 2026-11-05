// Which flowsheet a sensitivity study, optimization or dynamic run was computed for, so the
// editor can tell when its result no longer matches the flowsheet on screen.
import { useEffect, useState } from "react";
import { canonical } from "./clipboard";
import type { DynamicResult, FlowsheetModel, OptimizationResult, SensitivityResult } from "./types";

export type RunKind = "sensitivity" | "optimization" | "dynamics";

/**
 * A fingerprint of everything in `model` that can change the answer of one run: the model
 * without unit positions, name and description, and without the other runs' definitions (for
 * an optimization: the studies, the dynamics and the other optimizations; for a dynamic run
 * also the costing).
 */
export function runKey(model: FlowsheetModel, kind: RunKind, id?: string): string {
  const { name: _n, description: _d, sensitivity, optimizations, dynamics, costing, ...rest } = model;
  const units = rest.units.map(({ ui: _ui, ...u }) => u);
  // Key order does not matter: a parameter cleared and set again is the same flowsheet.
  return canonical({
    ...rest,
    costing: kind === "dynamics" ? undefined : costing,
    units,
    sensitivity: kind === "sensitivity" ? (sensitivity ?? []).filter((s) => s.id === id) : undefined,
    optimizations: kind === "optimization" ? (optimizations ?? []).filter((o) => o.id === id) : undefined,
    dynamics: kind === "dynamics" ? dynamics : undefined,
  });
}

/**
 * `compute()` once the inputs stopped changing for `ms` (a drag changes the model every
 * frame); null while none was computed yet. Callers check again before acting on it.
 */
export function useSettled<T>(compute: () => T, deps: unknown[], ms = 200): T | null {
  const [value, setValue] = useState<T | null>(null);
  useEffect(() => {
    const t = window.setTimeout(() => setValue(compute()), ms);
    return () => window.clearTimeout(t);
  }, deps);
  return value;
}

/** Which runs of `runs` still match `model`, by kind and id (no copies of the results). */
export function currentRunNames(runs: Runs, model: FlowsheetModel): string[] {
  const ok = (r: Run<unknown> | null | undefined, kind: RunKind, id?: string) => !!r && r.key === runKey(model, kind, id);
  return [
    ...(model.sensitivity ?? []).filter((s) => ok(runs.sensitivity[s.id], "sensitivity", s.id)).map((s) => `study ${s.id}`),
    ...(model.optimizations ?? []).filter((o) => ok(runs.optimizations[o.id], "optimization", o.id)).map((o) => `optimization ${o.id}`),
    ...(ok(runs.dynamics, "dynamics") ? ["dynamic run (Excel)"] : []),
  ];
}

export interface Run<T> {
  key: string;
  result: T;
}

/** The last result of every run, as the results panel keeps them for the report. */
export interface Runs {
  sensitivity: Record<string, Run<SensitivityResult>>;
  optimizations: Record<string, Run<OptimizationResult>>;
  dynamics: Run<DynamicResult> | null;
}

export const emptyRuns = (): Runs => ({ sensitivity: {}, optimizations: {}, dynamics: null });

/** The results of `runs` still valid for `model`, as /api/report takes them. */
export function currentRuns(runs: Runs, model: FlowsheetModel, withDynamics = true) {
  const valid = <T,>(r: Run<T> | null | undefined, kind: RunKind, id?: string) =>
    r && r.key === runKey(model, kind, id) ? r.result : null;
  const sensitivity = (model.sensitivity ?? []).flatMap((s) => {
    const r = valid(runs.sensitivity[s.id], "sensitivity", s.id);
    return r ? [r] : [];
  });
  const optimizations = (model.optimizations ?? []).flatMap((o) => {
    const r = valid(runs.optimizations[o.id], "optimization", o.id);
    // The solve at the optimum is the bulk of the result and not part of the report.
    return r ? [{ ...r, result: undefined }] : [];
  });
  const dyn = withDynamics ? valid(runs.dynamics, "dynamics") : null;
  const dynamics = dyn && {
    ...dyn,
    series: dyn.series.map((s) => ({ ...s, values: s.values.map(round) })),
    time: dyn.time.map(round),
  };
  return { sensitivity, optimizations, dynamics };
}

/** 6 significant digits are plenty for a report and keep the request small. */
const round = (v: number) => (Number.isFinite(v) ? Number(v.toPrecision(6)) : v);
