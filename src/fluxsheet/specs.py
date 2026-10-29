"""Design specifications: reach a target value by adjusting one unit parameter.

A design spec says "make this stream property (or unit result) equal to a value by varying
that parameter between two bounds", e.g. 99.5 % benzene in the distillate by varying the
reflux ratio. It is solved by an outer loop around the sequential-modular solver: each trial
value of the parameter is a complete flowsheet solve.

The loop starts with secant steps from the parameter's current value. As soon as two trials
lie on both sides of the target the root is bracketed and Brent's method finishes it, so
the search cannot leave the bracket again. Several specs are solved one after another and
repeated until all hold together.

Stream and unit variables are read from the solve results, so the same references serve
sensitivity studies (sensitivity.py).
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from scipy.optimize import brentq

from .errors import SimulationError
from .flowsheet import FlowsheetModel
from .quantities import UnitError, to_si
from .solver import Solver, SolveResult
from .study_models import DesignSpecModel, VariableRef, VariedParam
from .units import UNIT_TYPES

# Stream properties a spec can target: name -> (label, quantity, needs a component).
STREAM_PROPERTIES: dict[str, tuple[str, str | None, bool]] = {
    "T": ("temperature", "temperature", False),
    "P": ("pressure", "pressure", False),
    "vapor_fraction": ("vapor fraction", "fraction", False),
    "molar_flow": ("molar flow", "molar_flow", False),
    "mass_flow": ("mass flow", "mass_flow", False),
    "enthalpy_flow": ("enthalpy flow", "power", False),
    "mole_fraction": ("mole fraction", "fraction", True),
    "component_flow": ("molar flow", "molar_flow", True),
}
# Keys of Stream.to_dict for the per-component properties.
_STREAM_KEYS = {"mole_fraction": "mole_fractions", "component_flow": "component_flows"}

MAX_PASSES = 10  # rounds over all specs when there are several


# ---- reading results ---------------------------------------------------------------------

@dataclass
class Variable:
    """A validated reference to a computed value, with its label and quantity."""

    ref: VariableRef
    label: str
    quantity: str | None

    @classmethod
    def resolve(cls, ref: VariableRef, model: FlowsheetModel, where: str) -> Variable:
        if ref.stream is not None:
            streams = [s.id for s in model.streams]
            # Unconnected outlets appear in the results as "UNIT.port".
            if ref.stream not in streams and ref.stream.split(".")[0] not in _unit_ids(model):
                raise SimulationError(f"there is no stream called {ref.stream!r}", where=where,
                                      hint=f"choose one of: {', '.join(streams)}")
            if ref.property not in STREAM_PROPERTIES:
                raise SimulationError(
                    f"streams have no property {ref.property!r}", where=where,
                    hint=f"use one of: {', '.join(STREAM_PROPERTIES)}",
                )
            label, quantity, per_component = STREAM_PROPERTIES[ref.property]
            if per_component:
                _check_component(ref.component, model, where)
                label = f"{ref.component} {label}"
            return cls(ref, f"{ref.stream} {label}", quantity)

        cls_ = _unit_class(ref.unit, model, where)
        fields = {f.name: f for f in cls_.result_fields}
        if ref.property not in fields:
            raise SimulationError(
                f"{ref.unit} ({cls_.label}) reports no result {ref.property!r}", where=where,
                hint=f"use one of: {', '.join(fields) or 'none'}",
            )
        field = fields[ref.property]
        return cls(ref, f"{ref.unit} {_lower(field.label)}", field.quantity)

    def read(self, result: SolveResult) -> float:
        ref = self.ref
        if ref.stream is not None:
            stream = result.streams.get(ref.stream)
            if stream is None:
                raise SimulationError(f"stream {ref.stream} is not in the results",
                                      hint="connect the outlet, or pick another stream")
            value = stream[_STREAM_KEYS.get(ref.property, ref.property)]
            if ref.component is not None:
                value = value[ref.component]
        else:
            value = result.units.get(ref.unit, {}).get(ref.property)
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            raise SimulationError(f"{self.label} has no value in this solution",
                                  hint="check the unit's specification")
        return float(value)

    def to_si(self, value: float | str, what: str, where: str) -> float:
        try:
            return to_si(value, self.quantity) if self.quantity else float(value)
        except (UnitError, ValueError) as exc:
            raise SimulationError(f"{what}: {exc}", where=where) from exc


def _unit_ids(model: FlowsheetModel) -> list[str]:
    return [u.id for u in model.units]


def _unit_class(unit_id: str | None, model: FlowsheetModel, where: str):
    um = next((u for u in model.units if u.id == unit_id), None)
    if um is None:
        raise SimulationError(f"there is no unit called {unit_id!r}", where=where,
                              hint=f"choose one of: {', '.join(_unit_ids(model))}")
    cls = UNIT_TYPES.get(um.type)
    if cls is None:
        raise SimulationError(f"unknown unit type {um.type!r}", where=unit_id)
    return cls


def _check_component(name: str | None, model: FlowsheetModel, where: str) -> None:
    if name not in model.components:
        raise SimulationError(
            f"{'no component given' if name is None else repr(name) + ' is not a component'}",
            where=where, hint=f"set 'component' to one of: {', '.join(model.components)}",
        )


# ---- changing parameters -----------------------------------------------------------------

@dataclass
class Handle:
    """A validated unit parameter that can be set to a number (SI) in a copy of the model."""

    vary: VariedParam
    label: str
    quantity: str | None
    lower: float
    upper: float
    current: float | None  # value in the model, SI, if it is set

    @classmethod
    def resolve(cls, vary: VariedParam, model: FlowsheetModel, where: str) -> Handle:
        unit_cls = _unit_class(vary.unit, model, where)
        params = {p.name: p for p in unit_cls.params}
        param = params.get(vary.param)
        if param is None:
            raise SimulationError(
                f"{vary.unit} ({unit_cls.label}) has no parameter {vary.param!r}", where=where,
                hint=f"use one of: {', '.join(params)}",
            )
        q = param.quantity
        label = f"{vary.unit} {_lower(param.label)}"
        if q in ("composition", "component_fractions"):
            _check_component(vary.component, model, where)
            q = "molar_flow" if q == "composition" else "fraction"
            label = f"{label} of {vary.component}"
        elif q in ("integer", "choice", "component", "string_list"):
            raise SimulationError(
                f"{param.label} of {vary.unit} is not a continuous number", where=where,
                hint="vary a parameter such as a temperature, flow, duty or ratio; whole "
                "numbers like stage counts cannot be adjusted smoothly",
            )

        def si(value, what):
            try:
                return to_si(value, q) if q else float(value)
            except (UnitError, ValueError) as exc:
                raise SimulationError(f"{what} of {label}: {exc}", where=where) from exc

        lower, upper = si(vary.lower, "lower bound"), si(vary.upper, "upper bound")
        if not lower < upper:
            raise SimulationError(f"the lower bound of {label} is not below the upper bound",
                                  where=where, hint="swap or widen the bounds")
        for bound, name in ((lower, "lower"), (upper, "upper")):
            if (param.min is not None and bound < param.min) or (
                param.max is not None and bound > param.max
            ):
                raise SimulationError(
                    f"the {name} bound {bound:g} is outside the allowed range of {label}",
                    where=where,
                    hint=f"keep the bounds within {param.min} … {param.max}"
                    if param.max is not None else f"keep the bounds at or above {param.min}",
                )
        raw = next(u for u in model.units if u.id == vary.unit).params.get(vary.param)
        if vary.component is not None:
            raw = (raw or {}).get(vary.component)
        current = si(raw, "value") if raw is not None else None
        return cls(vary, label, q, lower, upper, current)

    def apply(self, model: FlowsheetModel, value: float) -> FlowsheetModel:
        """A copy of ``model`` with the parameter set to ``value`` (SI).

        Other members of the parameter's specification group are removed, so varying, say,
        a heater's duty replaces its outlet temperature specification.
        """
        out = model.model_copy(deep=True)
        unit = next(u for u in out.units if u.id == self.vary.unit)
        unit_cls = UNIT_TYPES[unit.type]
        for group in unit_cls.exclusive_groups:
            if self.vary.param in group:
                for other in group:
                    if other != self.vary.param:
                        unit.params.pop(other, None)
        if self.vary.component is not None:
            unit.params[self.vary.param] = {**(unit.params.get(self.vary.param) or {}),
                                            self.vary.component: value}
        else:
            unit.params[self.vary.param] = value
        return out


# ---- solving -----------------------------------------------------------------------------

class _Done(Exception):
    """Raised inside Brent's method once the target is met, to stop early."""


class _Spec:
    def __init__(self, spec: DesignSpecModel, model: FlowsheetModel):
        self.spec = spec
        self.target = Variable.resolve(spec.target, model, spec.id)
        self.handle = Handle.resolve(spec.vary, model, spec.id)
        self.value = self.target.to_si(spec.value, "target value", spec.id)
        if spec.tolerance is not None:
            # A temperature tolerance is a difference: "0.1 degC" means 0.1 K.
            q = self.target.quantity
            diff = Variable(spec.target, "", "temperature_difference" if q == "temperature" else q)
            self.tol = abs(diff.to_si(spec.tolerance, "tolerance", spec.id))
        else:
            self.tol = 1e-5 * max(abs(self.value), 1e-3)
        self.x: float | None = None
        self.achieved: float | None = None
        self.iterations = 0

    def error(self, result: SolveResult) -> float:
        return self.target.read(result) - self.value

    def ok(self, result: SolveResult) -> bool:
        return abs(self.error(result)) <= self.tol

    def solve(self, model: FlowsheetModel) -> tuple[FlowsheetModel, SolveResult]:
        """Adjust the parameter in ``model`` until the target holds; the other specs' parameters
        stay as they are."""
        h, spec = self.handle, self.spec
        lo, hi = h.lower, h.upper
        seen: dict[float, tuple[float, SolveResult]] = {}
        best: list = []  # [x, result] of the last trial that met the target

        def f(x: float) -> float:
            x = min(max(x, lo), hi)
            if x not in seen:
                if len(seen) >= spec.max_iterations:
                    raise self._not_converged(seen)
                trial = h.apply(model, x)
                result = Solver(trial).solve()
                if not result.ok:
                    raise _TrialFailed(x, result.error or {})
                try:
                    seen[x] = (self.error(result), result)
                except SimulationError as exc:
                    raise exc.located(spec.id)
                self.iterations += 1
            err, result = seen[x]
            if abs(err) <= self.tol:
                best[:] = [x, result]
                raise _Done
            return err

        try:
            self._search(f, seen, h.current)
        except _Done:
            x, result = best
            self.x, self.achieved = x, self.target.read(result)
            return h.apply(model, x), result
        except _TrialFailed as exc:
            raise self._trial_failed(exc) from None
        raise AssertionError("unreachable")  # _search ends by _Done or an error

    def _search(self, f, seen, start: float | None) -> None:
        lo, hi = self.handle.lower, self.handle.upper
        span = hi - lo
        x0 = start if start is not None and lo <= start <= hi else 0.5 * (lo + hi)
        x0, f0 = self._first(f, x0, lo, hi)
        x1 = x0 + 0.05 * span if x0 + 0.05 * span <= hi else x0 - 0.05 * span
        f1 = f(x1)
        while True:
            bracket = self._bracket(seen)
            if bracket:
                a, b = bracket
                brentq(f, a, b, xtol=1e-12 * span, maxiter=self.spec.max_iterations, disp=False)
                raise self._not_converged(seen)  # Brent stopped on x without meeting tol
            if f1 == f0:
                raise SimulationError(
                    f"{self.target.label} does not change with {self.handle.label}",
                    where=self.spec.id,
                    hint="vary a parameter that influences the target, or check that the "
                    "unit is connected upstream of it",
                )
            x2 = min(max(x1 - f1 * (x1 - x0) / (f1 - f0), lo), hi)
            if x2 in seen:
                # The secant returns to a known point, typically a bound that does not bracket
                # the target; the bounds decide whether the target is reachable at all.
                unseen = [b for b in (lo, hi) if b not in seen]
                if not unseen:
                    raise self._not_bracketed(seen)
                x2 = unseen[0]
            x0, f0, x1, f1 = x1, f1, x2, f(x2)

    def _first(self, f, x0, lo, hi) -> tuple[float, float]:
        """The first trial; if the flowsheet fails there, try the middle of the bounds."""
        try:
            return x0, f(x0)
        except _TrialFailed:
            mid = 0.5 * (lo + hi)
            if x0 == mid:
                raise
            return mid, f(mid)

    @staticmethod
    def _bracket(seen) -> tuple[float, float] | None:
        xs = sorted(seen)
        pairs = [(a, b) for a, b in itertools.pairwise(xs) if seen[a][0] * seen[b][0] < 0]
        return min(pairs, key=lambda p: p[1] - p[0]) if pairs else None

    def _fmt(self, value: float, quantity: str | None) -> str:
        if quantity == "temperature":
            return f"{value - 273.15:.4g} °C"
        scale, unit = {"pressure": (1e-5, " bar"), "molar_flow": (3.6, " kmol/h"),
                       "power": (1e-3, " kW"), "mass_flow": (3600.0, " kg/h")}.get(
            quantity or "", (1.0, ""))
        return f"{value * scale:.6g}{unit}"

    def _not_bracketed(self, seen) -> SimulationError:
        h, t = self.handle, self.target
        lo_val, hi_val = seen[h.lower][0] + self.value, seen[h.upper][0] + self.value
        return SimulationError(
            f"{t.label} cannot reach {self._fmt(self.value, t.quantity)}: with {h.label} "
            f"between {self._fmt(h.lower, h.quantity)} and {self._fmt(h.upper, h.quantity)} it "
            f"goes from {self._fmt(lo_val, t.quantity)} to {self._fmt(hi_val, t.quantity)}",
            where=self.spec.id,
            hint="the bounds do not bracket the target: widen them, or choose a target "
            "within that range",
        )

    def _not_converged(self, seen) -> SimulationError:
        x = min(seen, key=lambda k: abs(seen[k][0])) if seen else None
        detail = "" if x is None else (
            f"; closest: {self._fmt(seen[x][0] + self.value, self.target.quantity)} at "
            f"{self.handle.label} = {self._fmt(x, self.handle.quantity)}")
        return SimulationError(
            f"design spec did not converge in {len(seen)} flowsheet solves{detail}",
            where=self.spec.id,
            hint="raise max_iterations, loosen the tolerance, or narrow the bounds around the "
            "solution",
        )

    def _trial_failed(self, exc: _TrialFailed) -> SimulationError:
        err = exc.error
        place = f"{err.get('where')}: " if err.get("where") else ""
        return SimulationError(
            f"the flowsheet failed with {self.handle.label} = "
            f"{self._fmt(exc.x, self.handle.quantity)} ({place}{err.get('message')})",
            where=self.spec.id,
            hint=err.get("hint") or "narrow the bounds to values where the flowsheet solves",
        )

    def report(self) -> dict:
        h, t = self.handle, self.target
        return {
            "id": self.spec.id, "converged": True,
            "target": t.label, "quantity": t.quantity, "value": self.value,
            "achieved": self.achieved,
            "vary": h.label, "vary_quantity": h.quantity, "vary_value": self.x,
            "iterations": self.iterations,
        }


class _TrialFailed(Exception):
    def __init__(self, x: float, error: dict):
        super().__init__(error.get("message"))
        self.x, self.error = x, error


def solve_with_specs(model: FlowsheetModel) -> SolveResult:
    """Solve ``model`` with its active design specs met; see the module docstring."""
    work = model.model_copy(deep=True, update={"design_specs": None})
    try:
        specs = [_Spec(s, work) for s in model.design_specs or () if s.active]
        _check_distinct(specs)
        result = None
        for _ in range(MAX_PASSES):
            for spec in specs:
                work, result = spec.solve(work)
            if all(spec.ok(result) for spec in specs):
                break
        else:
            raise SimulationError(
                f"the {len(specs)} design specs did not hold together after {MAX_PASSES} rounds",
                where=", ".join(s.spec.id for s in specs),
                hint="each spec undoes another; make sure every spec varies a parameter that "
                "mainly affects its own target",
            )
    except SimulationError as exc:
        failed = SolveResult(ok=False, error=exc.to_dict())
        return failed
    for spec in specs:
        spec.achieved = spec.target.read(result)
    result.design_specs = [spec.report() for spec in specs]
    return result


def _check_distinct(specs: list[_Spec]) -> None:
    seen: dict[tuple, str] = {}
    for s in specs:
        key = (s.spec.vary.unit, s.spec.vary.param, s.spec.vary.component)
        if key in seen:
            raise SimulationError(
                f"{s.handle.label} is already varied by design spec {seen[key]}",
                where=s.spec.id, hint="let each design spec vary a different parameter",
            )
        seen[key] = s.spec.id


def _lower(label: str) -> str:
    """'Reflux ratio' -> 'reflux ratio' inside a sentence, but keep 'UA' and 'LMTD'."""
    return label[0].lower() + label[1:] if label[1:2].islower() else label
