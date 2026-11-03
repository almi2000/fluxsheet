"""Optimization: minimize or maximize an expression by varying unit parameters.

Like Aspen's Optimization block, the flowsheet is a black box: each trial point is a complete
solve (with design specs, calculators and costing), and the objective and constraints are
expressions (calculator.py) over flowsheet values, e.g. ``tac`` from the costing or
``0.3 * product - tac``. Run on request; see ``OptimizationModel`` for the file format.

Method:
- The decision variables are scaled to [0, 1] between their bounds, the objective by its
  magnitude at the start, each constraint by the larger of its value and its start, so
  every function is of order one.
- SLSQP (``scipy.optimize.minimize``) with forward differences of ``step`` (a fraction of the
  range; backward at the upper bound or where the forward point fails). Each trial is
  solved once: SLSQP asks for the objective and the constraints at the same point, and the
  results are cached. Inside the trials the recycle tolerance is tightened to 1e-8 and that
  of design specs to 1e-7 relative, so the differences see the change of the variables,
  not the convergence noise.
- When a trial flowsheet fails during the search, or SLSQP gives up for another reason
  than its iteration limit, the search goes on with COBYLA (derivative free; failed trials
  count as worse than any point seen) from the best feasible point so far (``method:
  "auto"``).
- The answer is the best feasible point evaluated: constraints met to 1e-4 relative.
  If none is feasible the run fails with the least violated point.

Each iteration's objective goes into the history. On the demo server the time budget
(limits.py) stops the run; the best point found so far is reported then.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from . import limits
from .calculator import CONSTANTS, FUNCTIONS, Expression, _read_param
from .errors import SimulationError
from .flowsheet import FlowsheetModel
from .solver import SolveResult, solve
from .specs import Handle, Variable, _Spec
from .study_models import OptimizationModel, ParamRef

FEASIBILITY = 1e-4  # largest scaled constraint violation of a feasible point
ACTIVE = 1e-3  # a constraint within this (scaled) of its limit is reported as active
RECYCLE_TOLERANCE = 1e-8
SPEC_TOLERANCE = 1e-7  # relative to the target value


@dataclass
class _Point:
    z: tuple[float, ...]
    result: SolveResult
    objective: float | None = None  # raw, as the expression gives it
    constraints: list[float] | None = None  # raw left-hand sides
    error: dict | None = None

    @property
    def ok(self) -> bool:
        return self.error is None


class _Abort(Exception):
    """Leaves SLSQP when a trial fails; the search goes on with COBYLA."""

    def __init__(self, point: _Point):
        super().__init__((point.error or {}).get("message"))
        self.point = point


class _Problem:
    def __init__(self, model: FlowsheetModel, opt: OptimizationModel):
        self.opt = opt
        where = opt.id
        work = model.model_copy(deep=True, update={"optimizations": None})
        work.solver.tolerance = min(work.solver.tolerance, RECYCLE_TOLERANCE)
        for spec in work.design_specs or ():
            if spec.active and spec.tolerance is None:
                value = _Spec(spec, work).value
                spec.tolerance = SPEC_TOLERANCE * max(abs(value), 1e-3)
        self.model = work
        self.handles = [Handle.resolve(v, work, where) for v in opt.vary]
        self._check_conflicts()
        self.results: dict[str, Variable] = {}
        self.params: dict[str, ParamRef] = {}
        for name, ref in opt.variables.items():
            if not name.isidentifier() or name in FUNCTIONS or name in CONSTANTS:
                raise SimulationError(f"{name!r} cannot be a variable name", where=where,
                                      hint="use letters, digits and _, not a function name")
            if isinstance(ref, ParamRef):
                _read_param(work, ref, where)  # checks the reference
                self.params[name] = ref
            else:
                self.results[name] = Variable.resolve(ref, work, where)
        names = set(opt.variables)
        self.objective = Expression(opt.objective, names, where)
        self.constraints = [Expression(c.expression, names, f"{where}, constraint {k + 1}")
                            for k, c in enumerate(opt.constraints)]
        self.sign = -1.0 if opt.sense == "maximize" else 1.0
        self.cache: dict[tuple[float, ...], _Point] = {}
        self.f_scale = 1.0
        self.g_scale = [1.0] * len(self.constraints)
        self.worst = 0.0  # largest scaled objective seen, for COBYLA's failed trials

    def _check_conflicts(self) -> None:
        taken = {}
        for spec in self.model.design_specs or ():
            if spec.active:
                taken[(spec.vary.unit, spec.vary.param, spec.vary.component)] = \
                    f"design spec {spec.id}"
        for calc in self.model.calculators or ():
            if calc.active:
                for a in calc.assign:
                    taken[(a.unit, a.param, a.component)] = f"calculator {calc.id}"
        seen = set()
        for h in self.handles:
            key = (h.vary.unit, h.vary.param, h.vary.component)
            if key in taken:
                raise SimulationError(f"{h.label} is already set by {taken[key]}",
                                      where=self.opt.id,
                                      hint="vary a parameter nothing else adjusts")
            if key in seen:
                raise SimulationError(f"{h.label} is varied twice", where=self.opt.id)
            seen.add(key)

    # ---- points --------------------------------------------------------------------------

    def x(self, z) -> list[float]:
        return [h.lower + min(max(zi, 0.0), 1.0) * (h.upper - h.lower)
                for h, zi in zip(self.handles, z)]

    def start(self) -> list[float]:
        z = []
        for h, v in zip(self.handles, self.opt.vary):
            x = h.current
            if v.start is not None:
                x = Variable(None, h.label, h.quantity).to_si(v.start, "start", self.opt.id)
            x = 0.5 * (h.lower + h.upper) if x is None else min(max(x, h.lower), h.upper)
            z.append((x - h.lower) / (h.upper - h.lower))
        return z

    def evaluate(self, z) -> _Point:
        key = tuple(float(min(max(v, 0.0), 1.0)) for v in z)
        if key in self.cache:
            return self.cache[key]
        limits.check_deadline(self.opt.id)
        trial = self.model
        for h, x in zip(self.handles, self.x(key)):
            trial = h.apply(trial, x)
        result = solve(trial)
        point = _Point(key, result)
        if not result.ok:
            point.error = result.error or {"message": "the flowsheet failed"}
            limits.check_deadline(self.opt.id)  # a timeout is no failed trial
        else:
            try:
                values = {n: v.read(result) for n, v in self.results.items()}
                for name, ref in self.params.items():
                    value = _read_param(trial, ref, self.opt.id)
                    if value is None:
                        raise SimulationError(f"variable {name}: {ref.param} of {ref.unit} is "
                                              "not set", where=self.opt.id)
                    values[name] = value
                point.objective = self.objective.evaluate(values)
                point.constraints = [c.evaluate(values) for c in self.constraints]
            except SimulationError as exc:
                point.error = exc.located(self.opt.id).to_dict()
        if point.ok:
            self.worst = max(self.worst, self.f(point))
        self.cache[key] = point
        return point

    def scale(self, start: _Point) -> None:
        self.f_scale = abs(start.objective) if start.objective else 1.0
        self.g_scale = [max(abs(c.value), abs(lhs)) or 1.0
                        for c, lhs in zip(self.opt.constraints, start.constraints)]

    # ---- scaled functions ----------------------------------------------------------------

    def f(self, p: _Point) -> float:
        return self.sign * p.objective / self.f_scale

    def g(self, p: _Point) -> list[float]:
        """Scaled constraints, each >= 0 (inequalities) or = 0 (equalities) when met."""
        out = []
        for c, lhs, s in zip(self.opt.constraints, p.constraints, self.g_scale):
            out.append((c.value - lhs) / s if c.relation == "<=" else (lhs - c.value) / s)
        return out

    def violation(self, p: _Point) -> float:
        worst = 0.0
        for c, gi in zip(self.opt.constraints, self.g(p)):
            worst = max(worst, abs(gi) if c.relation == "=" else -gi)
        return worst

    def gradients(self, z) -> tuple[np.ndarray, np.ndarray]:
        """Forward differences of the scaled objective and constraints at ``z``."""
        base = self.evaluate(z)
        if not base.ok:
            raise _Abort(base)
        f0, g0 = self.f(base), np.array(self.g(base))
        h = self.opt.step
        df = np.zeros(len(z))
        dg = np.zeros((len(g0), len(z)))
        for i in range(len(z)):
            steps = (h, -h) if z[i] + h <= 1.0 else (-h, h)
            for s in steps:
                zz = list(z)
                zz[i] = z[i] + s
                if not 0.0 <= zz[i] <= 1.0:
                    continue
                p = self.evaluate(zz)
                if p.ok:
                    break
            else:
                raise _Abort(p)
            df[i] = (self.f(p) - f0) / s
            dg[:, i] = (np.array(self.g(p)) - g0) / s
        return df, dg


# ---- running ---------------------------------------------------------------------------------

def run_optimization(model: FlowsheetModel, opt_id: str | None = None) -> dict:
    """Run one optimization of ``model`` (the first if ``opt_id`` is None); see the module
    docstring. Returns ``{"ok": True, ...}`` with the optimum, or ``{"ok": False, "error"}``
    (with the history so far when the search ran)."""
    try:
        opt = _pick(model, opt_id)
        problem = _Problem(model, opt)
        return _search(problem)
    except SimulationError as exc:
        return {"ok": False, "error": exc.to_dict()}


def _search(pb: _Problem) -> dict:
    opt = pb.opt
    z0 = pb.start()
    first = pb.evaluate(z0)
    if not first.ok:  # try the middle of the bounds before giving up
        mid = [0.5] * len(z0)
        if pb.evaluate(mid).ok:
            z0, first = mid, pb.evaluate(mid)
        else:
            raise SimulationError(
                f"the flowsheet fails at the start point: {_describe(first.error)}",
                where=opt.id, hint=(first.error or {}).get("hint")
                or "choose a start value or bounds where the flowsheet solves")
    pb.scale(first)
    history = [_entry(pb, 0, first)]
    n = len(z0)
    ineq = [k for k, c in enumerate(opt.constraints) if c.relation != "="]
    eq = [k for k, c in enumerate(opt.constraints) if c.relation == "="]
    used, message, converged, retry, stopped = [], "", False, False, None

    def record(zk) -> None:
        p = pb.evaluate(zk)
        if p.ok and pb.x(p.z) != history[-1]["x"]:  # a last call may repeat the point
            history.append(_entry(pb, len(history), p))

    try:
        if opt.method in ("auto", "slsqp"):
            used.append("SLSQP")

            def fun(z):
                p = pb.evaluate(z)
                if not p.ok:
                    raise _Abort(p)
                return pb.f(p)

            def con(idx):
                def g(z):
                    p = pb.evaluate(z)
                    if not p.ok:
                        raise _Abort(p)
                    return np.array([pb.g(p)[k] for k in idx])

                def jac(z):
                    return pb.gradients(z)[1][idx, :]
                return g, jac

            constraints = []
            for kind, idx in (("ineq", ineq), ("eq", eq)):
                if idx:
                    g, jac = con(idx)
                    constraints.append({"type": kind, "fun": g, "jac": jac})
            try:
                res = minimize(fun, np.array(z0), jac=lambda z: pb.gradients(z)[0],
                               method="SLSQP", bounds=[(0.0, 1.0)] * n,
                               constraints=constraints, callback=record,
                               options={"maxiter": opt.max_iterations, "ftol": opt.tolerance})
                message, converged = res.message, res.status == 0
                retry = res.status not in (0, 9)
            except _Abort as exc:
                message = f"a trial failed ({_describe(exc.point.error)})"
                retry = True
            if retry and opt.method == "auto":
                message = f"SLSQP: {message}; continued with COBYLA"
            elif retry:
                message = f"SLSQP: {message}"
        if opt.method == "cobyla" or (opt.method == "auto" and not converged and retry):
            used.append("COBYLA")
            best = _best(pb)
            start = best.z if best is not None else z0

            def fun_c(z):
                p = pb.evaluate(z)
                return pb.f(p) if p.ok else pb.worst + 1.0

            def con_c(z):
                p = pb.evaluate(z)
                if not p.ok:
                    return -np.ones(len(ineq) + 2 * len(eq))
                g = pb.g(p)
                return np.array([g[k] for k in ineq] + [g[k] for k in eq]
                                + [-g[k] for k in eq])

            cons = [{"type": "ineq", "fun": con_c}] if opt.constraints else []
            res = minimize(fun_c, np.array(start), method="COBYLA",
                           bounds=[(0.0, 1.0)] * n, constraints=cons, callback=record,
                           options={"maxiter": max(opt.max_iterations * (n + 1), 20),
                                    "rhobeg": 0.1, "tol": opt.tolerance})
            converged = bool(res.success)
            message = f"{message}: {res.message}" if message else res.message
    except SimulationError as exc:  # the demo server's time limit
        stopped = exc.to_dict()
        message = exc.message

    best = _best(pb)
    report = {
        "id": opt.id, "sense": opt.sense, "method": " → ".join(used),
        "converged": converged and stopped is None, "message": str(message),
        "iterations": len(history) - 1, "evaluations": len(pb.cache),
        "failed_evaluations": sum(1 for p in pb.cache.values() if not p.ok),
        "history": history,
    }
    if best is None:
        least = min((p for p in pb.cache.values() if p.ok), key=pb.violation)
        return {"ok": False, **report, "error": {
            "message": "no feasible point found: the constraints are violated everywhere the "
            f"search went (least violation {pb.violation(least):.3g}, relative)",
            "where": opt.id,
            "hint": "check that the constraints can hold together within the bounds, "
            "widen the bounds, or start from a feasible point"} if stopped is None else stopped,
            **_describe_point(pb, least, z0)}
    return {"ok": True, **report, **_describe_point(pb, best, z0)}


def _best(pb: _Problem) -> _Point | None:
    feasible = [p for p in pb.cache.values() if p.ok and pb.violation(p) <= FEASIBILITY]
    return min(feasible, key=pb.f) if feasible else None


def _entry(pb: _Problem, k: int, p: _Point) -> dict:
    return {"iteration": k, "objective": p.objective, "x": pb.x(p.z),
            "max_violation": pb.violation(p)}


def _describe_point(pb: _Problem, p: _Point, z0) -> dict:
    x, x0 = pb.x(p.z), pb.x(z0)
    start = pb.cache.get(tuple(float(v) for v in z0))
    span = [h.upper - h.lower for h in pb.handles]
    g = pb.g(p)
    return {
        "objective": {"expression": pb.opt.objective, "value": p.objective,
                      "start": start.objective if start is not None else None},
        "variables": [
            {"label": h.label, "quantity": h.quantity, "unit": h.vary.unit,
             "param": h.vary.param, "component": h.vary.component, "lower": h.lower,
             "upper": h.upper, "start": s, "value": v,
             "at_bound": min(v - h.lower, h.upper - v) <= 1e-6 * d}
            for h, s, v, d in zip(pb.handles, x0, x, span)
        ],
        "constraints": [
            {"expression": c.expression, "relation": c.relation, "value": c.value,
             "achieved": lhs, "satisfied": (abs(gi) if c.relation == "=" else -gi)
             <= FEASIBILITY, "active": abs(gi) <= ACTIVE}
            for c, lhs, gi in zip(pb.opt.constraints, p.constraints, g)
        ],
        "apply": [{"unit": h.vary.unit, "param": h.vary.param, "component": h.vary.component,
                   "value": v} for h, v in zip(pb.handles, x)],
        "result": p.result.to_dict(),
    }


def _describe(error: dict | None) -> str:
    error = error or {}
    where = f"{error['where']}: " if error.get("where") else ""
    return f"{where}{error.get('message')}"


def _pick(model: FlowsheetModel, opt_id: str | None) -> OptimizationModel:
    opts = model.optimizations or []
    if not opts:
        raise SimulationError("the flowsheet has no optimization",
                              hint="add one under 'optimizations' (in the editor: flowsheet "
                              "settings, Optimization)")
    if opt_id is None:
        return opts[0]
    for o in opts:
        if o.id == opt_id:
            return o
    raise SimulationError(f"there is no optimization {opt_id!r}",
                          hint=f"choose one of: {', '.join(o.id for o in opts)}")
