"""Dynamic simulation of a flowsheet.

1. Solve the steady state (the ordinary sequential-modular solve).
2. Units listed under ``dynamics.units`` get their dynamic model (units.py) and start from
   their steady state; controllers start bumpless at the steady measurement and output.
3. Integrate with a stiff BDF method. Every right-hand side evaluation sets the event and
   controller values, takes the outlets of the hold-ups from their states, recalculates the
   algebraic units in flow order from there (cached, so finite-difference Jacobians only redo
   what a perturbed state reaches) and collects the derivatives. A recycle loop through a
   hold-up needs no iteration; a loop of algebraic units only is iterated as in the steady
   solve. The Jacobian's sparsity follows from which states can reach which units.
4. Events split the time axis, so the integrator restarts cleanly at each step or ramp.

Results are time series of the streams (T, flow, mole fractions), the dynamic units (level,
temperature, ...) and the controllers (measurement, setpoint, output), all in SI units.
"""

from __future__ import annotations

import csv
import io
import itertools
import math
from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

from ..calculator import _param_quantity
from ..errors import SimulationError
from ..flowsheet import FlowsheetModel
from ..limits import check_deadline
from ..quantities import UnitError, to_si
from ..solver import Solver
from ..specs import STREAM_PROPERTIES, Variable
from ..stream import Stream
from ..units.base import mix_streams
from .controllers import Controller
from .models import ControllerModel, EventModel, ManipulatedParam
from .units import DYNAMIC_MODELS, DynamicModel

CACHE_SIZE = 64


def simulate(model: FlowsheetModel | dict) -> dict:
    """Run ``model.dynamics``. Returns ``{"ok", "time", "series", "warnings"}`` or, on an
    error, ``{"ok": False, "error", ...}`` with the series up to the failure, if any."""
    if isinstance(model, dict):
        model = FlowsheetModel.model_validate(model)
    return Simulator(model).run()


@dataclass
class _Event:
    t: float
    ramp: float
    value: float
    target: tuple  # ("controller", id) or ("unit", unit, param, component)


class Simulator(Solver):
    def run(self) -> dict:
        out: dict = {"ok": False, "time": [], "series": [], "warnings": []}
        dyn = self.model.dynamics
        if dyn is None:
            out["error"] = SimulationError(
                "the flowsheet has no dynamic simulation set up",
                hint="add a 'dynamics' section with t_end and the units that hold up material",
            ).to_dict()
            return out
        if any(s.active for s in self.model.design_specs or ()):
            self.warnings.append("design specs are not applied in a dynamic simulation; it "
                                 "starts from the parameters as they are")
        steady = self.solve()
        if not steady.ok:
            out["error"], out["warnings"] = steady.error, steady.warnings
            return out
        out["steady"] = {"iterations": steady.iterations}
        try:
            self._setup()
            self._integrate(out)
            out["ok"] = True
        except SimulationError as exc:
            out["error"] = exc.to_dict()
        out["warnings"] = self.warnings
        return out

    # ---- setup -------------------------------------------------------------------------

    def _setup(self) -> None:
        dyn = self.model.dynamics
        self.t_end = _si(dyn.t_end, "time", "dynamics", "end time")
        if self.t_end <= 0.0:
            raise SimulationError("the end time must be positive", where="dynamics")
        self.by_source = {(c.src_unit, c.src_port): c.id for c in self.connections}
        self.snapshot: dict[str, Stream] = dict(self.streams)
        self.cache: dict[str, OrderedDict] = {uid: OrderedDict() for uid in self.units}

        self.dynamic: dict[str, DynamicModel] = {}
        blocks: list[tuple[str, int]] = []  # (owner, size)
        x0: list[np.ndarray] = []
        scales: list[np.ndarray] = []
        for uid, settings in dyn.units.items():
            unit = self.units.get(uid)
            if unit is None:
                raise SimulationError(f"there is no unit called {uid!r}", where="dynamics")
            cls = DYNAMIC_MODELS.get(unit.type)
            if cls is None:
                raise SimulationError(
                    f"{unit.label} has no dynamic model yet", where=uid,
                    hint=f"dynamic models exist for: {', '.join(DYNAMIC_MODELS)}",
                )
            model = cls(unit, self.sections.package(uid), settings)
            if not model.source:
                model.steady = self._steady_of(uid)
            x = model.initialize(self._inlets(uid), self._outlets_of(uid))
            self.dynamic[uid] = model
            blocks.append((uid, len(x)))
            x0.append(x)
            scales.append(model.scale(x))
        for (uid, model), x in zip(self.dynamic.items(), x0):  # their results at the start
            model.outlets(x, None if model.source else self._inlets(uid))

        self.controllers: dict[str, Controller] = {}
        for cm in dyn.controllers:
            if cm.active:
                ctrl = self._controller(cm)
                self.controllers[cm.id] = ctrl
                x = ctrl.initial(self._read(ctrl.measure))
                blocks.append((cm.id, len(x)))
                x0.append(x)
                scales.append(ctrl.scale())

        self.slices: dict[str, slice] = {}
        start = 0
        for owner, size in blocks:
            self.slices[owner] = slice(start, start + size)
            start += size
        if start == 0:
            raise SimulationError(
                "nothing in the flowsheet changes over time", where="dynamics",
                hint="list the units with a hold-up under dynamics.units, e.g. "
                f"{{\"TANK\": {{}}}}; dynamic models exist for {', '.join(DYNAMIC_MODELS)}",
            )
        self.x0 = np.concatenate(x0)
        self.xscale = np.concatenate(scales)

        self.events = [self._event(ev) for ev in dyn.events]
        self.base_values = {ev.target: self._current(ev.target) for ev in self.events}
        self._order_algebraic()
        self.sparsity = self._sparsity()

    def _steady_of(self, uid: str):
        def steady(inlets):
            return self._calc(uid, inlets)
        return steady

    def _inlets(self, uid: str) -> dict[str, list[Stream]]:
        unit = self.units[uid]
        inlets: dict[str, list[Stream]] = {p.name: [] for p in unit.inlets}
        unit.inlet_ids = {p.name: [] for p in unit.inlets}
        for c in self.connections:
            if c.dst_unit == uid:
                inlets[c.dst_port].append(self.sections.inlet(c, self.snapshot[c.id]))
                unit.inlet_ids[c.dst_port].append(c.id)
        return inlets

    def _outlets_of(self, uid: str) -> dict[str, Stream]:
        unit = self.units[uid]
        return {p.name: self.snapshot[self.by_source.get((uid, p.name), f"{uid}.{p.name}")]
                for p in unit.outlets
                if self.by_source.get((uid, p.name), f"{uid}.{p.name}") in self.snapshot}

    def _controller(self, cm: ControllerModel) -> Controller:
        where = cm.id
        self.measure_quantity = None
        ref = cm.measure
        if ref.stream is None and ref.unit is None:
            raise SimulationError(f"a controller cannot measure the costing total {ref.property}",
                                  where=where, hint="measure a stream or unit value")
        if ref.unit is not None and ref.unit in self.model.dynamics.units:
            dyn_cls = DYNAMIC_MODELS.get(self.units[ref.unit].type)
            fields = {f.name: f for f in dyn_cls.result_fields} if dyn_cls else {}
            if ref.property in fields:
                quantity = fields[ref.property].quantity
            else:
                quantity = Variable.resolve(ref, self.model, where).quantity
        else:
            quantity = Variable.resolve(ref, self.model, where).quantity
        measured = self._read(ref)
        sp = measured if cm.setpoint is None else _si(cm.setpoint, quantity, where, "setpoint")
        if cm.pv_range is not None:
            lo, hi = (_si(v, quantity, where, "measurement range") for v in cm.pv_range)
        elif quantity == "temperature":
            lo, hi = sp - 50.0, sp + 50.0
        else:
            half = 0.5 * abs(sp) if sp else 1.0
            lo, hi = sp - half, sp + half
        if not lo < hi:
            raise SimulationError("the measurement range is empty", where=where)
        mq = self._manipulated_quantity(cm.manipulate, where)
        op_lo, op_hi = (_si(v, mq, where, "output range") for v in cm.op_range)
        if not op_lo < op_hi:
            raise SimulationError("the output range is empty", where=where,
                                  hint="give op_range as [low, high]")
        current = self._current(_target(cm.manipulate))
        bias = (current - op_lo) / (op_hi - op_lo) * 100.0
        if not 0.0 <= bias <= 100.0:
            self.warnings.append(f"{where}: the steady-state output lies outside the output "
                                 "range; the controller starts at its limit")
        ti = _si(cm.integral_time, "time", where, "integral time") if cm.integral_time else None
        td = (_si(cm.derivative_time, "time", where, "derivative time")
              if cm.derivative_time else None)
        lag = max(_si(cm.measurement_lag, "time", where, "measurement lag"), 1e-3)
        ctrl = Controller(cm, lo, hi, op_lo, op_hi, sp, bias, ti, td, lag)
        ctrl.measure = ref
        ctrl.quantity = quantity
        ctrl.op_quantity = mq
        return ctrl

    def _manipulated_quantity(self, m: ManipulatedParam, where: str) -> str | None:
        unit = self.units.get(m.unit)
        if unit is None:
            raise SimulationError(f"there is no unit called {m.unit!r}", where=where)
        dyn_cls = DYNAMIC_MODELS.get(unit.type) if m.unit in self.model.dynamics.units else None
        if dyn_cls is not None:
            for p in dyn_cls.params:
                if p.name == m.param:
                    return p.quantity
        return _param_quantity(self.model, m.unit, m.param, m.component, where)

    def _event(self, ev: EventModel) -> _Event:
        where = "dynamics.events"
        t = _si(ev.time, "time", where, "event time")
        ramp = _si(ev.ramp, "time", where, "ramp") if ev.ramp is not None else 0.0
        if ev.controller is not None:
            ctrl = self.controllers.get(ev.controller)
            if ctrl is None:
                raise SimulationError(f"there is no active controller {ev.controller!r}",
                                      where=where)
            value = _si(ev.value, ctrl.quantity, where, "setpoint")
            return _Event(t, ramp, value, ("controller", ev.controller))
        m = ManipulatedParam(unit=ev.unit, param=ev.param, component=ev.component)
        q = self._manipulated_quantity(m, where)
        return _Event(t, ramp, _si(ev.value, q, where, "value"), _target(m))

    # ---- values that controllers and events change -------------------------------------

    def _current(self, target: tuple) -> float:
        if target[0] == "controller":
            return self.controllers[target[1]].setpoint
        _, uid, param, comp = target
        model = self.dynamic.get(uid)
        if model is not None and param in model.s:
            return model.s[param]
        unit = self.units[uid]
        value = unit.p.get(param)
        if comp is not None:
            value = (value or {}).get(comp, 0.0)
        if value is None:  # e.g. a heater's duty when it is specified by temperature
            value = unit.results.get(param)
        if not isinstance(value, (int, float)):
            raise SimulationError(f"{param} of {uid} has no value to start from", where=uid)
        return float(value)

    def _set(self, target: tuple, value: float) -> None:
        if target[0] == "controller":
            self.controllers[target[1]].setpoint = value
            return
        _, uid, param, comp = target
        model = self.dynamic.get(uid)
        if model is not None and model.set(param, value):
            return
        unit = self.units[uid]
        for group in unit.exclusive_groups:
            if param in group:
                for other in group:
                    if other != param:
                        unit.p[other] = None
        if comp is not None:
            unit.p[param] = {**(unit.p.get(param) or {}), comp: value}
        else:
            unit.p[param] = value

    def _event_value(self, target: tuple, t: float) -> float | None:
        evs = sorted((e for e in self.events if e.target == target), key=lambda e: e.t)
        if not evs or t < evs[0].t:
            return None
        v = self.base_values[target]
        for ev in evs:
            if t < ev.t:
                break
            base = v
            frac = 1.0 if ev.ramp <= 0.0 else min((t - ev.t) / ev.ramp, 1.0)
            v = base + (ev.value - base) * frac
        return v

    # ---- evaluating the flowsheet at one instant ---------------------------------------

    def _order_algebraic(self) -> None:
        """Flow order of the units that are recalculated at every instant, and the streams
        torn in loops without a hold-up."""
        sources = {uid for uid, m in self.dynamic.items() if m.source}
        nodes = [u for u in self.units if u not in sources]
        out_edges = {u: [c for c in self.connections if c.src_unit == u and c.dst_unit
                         and c.dst_unit not in sources] for u in nodes}
        tears: set[str] = set()
        state: dict[str, int] = {}

        def visit(u: str) -> None:
            state[u] = 1
            for c in out_edges[u]:
                if state.get(c.dst_unit) == 1:
                    tears.add(c.id)
                elif c.dst_unit not in state:
                    visit(c.dst_unit)
            state[u] = 2

        for u in nodes:
            if u not in state:
                visit(u)
        indegree = {u: 0 for u in nodes}
        for u in nodes:
            for c in out_edges[u]:
                if c.id not in tears:
                    indegree[c.dst_unit] += 1
        ready = [u for u in nodes if indegree[u] == 0]
        order = []
        while ready:
            u = ready.pop(0)
            order.append(u)
            for c in out_edges[u]:
                if c.id not in tears:
                    indegree[c.dst_unit] -= 1
                    if indegree[c.dst_unit] == 0:
                        ready.append(c.dst_unit)
        self.sources, self.order, self.tears = sorted(sources), order, tears
        if tears:
            self.warnings.append(
                f"the recycle through {', '.join(sorted(tears))} has no hold-up and is "
                "iterated at every time step; give a unit in it a dynamic model to speed "
                "this up")

    def _calc(self, uid: str, inlets: dict[str, list[Stream]] | None = None
              ) -> dict[str, Stream]:
        """Steady-state outlets of an algebraic unit, cached on its inlets and parameters."""
        unit = self.units[uid]
        if inlets is None:
            inlets = self._inlets(uid)
        key = (tuple((port, tuple((tuple(s.flows), s.T, s.P) for s in streams))
                     for port, streams in inlets.items()),
               repr(sorted(unit.p.items())))
        cache = self.cache[uid]
        hit = cache.get(key)
        if hit is not None:
            cache.move_to_end(key)
            unit.results.clear()
            unit.results.update(hit[1])
            return hit[0]
        check_deadline(uid)
        unit_pkg = self.sections.package(uid)
        from .. import utilities

        unit.notes = []  # not reported during the integration
        try:
            for port in unit.inlets:
                if len(inlets[port.name]) > 1 and not port.separate:
                    inlets = {**inlets, port.name: [mix_streams(unit_pkg, inlets[port.name])]}
            outputs = utilities.calculate(unit, unit_pkg, inlets)
        except SimulationError as exc:
            raise exc.located(uid)
        except Exception as exc:
            raise SimulationError(f"calculation failed: {exc}", where=uid) from exc
        cache[key] = (outputs, dict(unit.results))
        if len(cache) > CACHE_SIZE:
            cache.popitem(last=False)
        return outputs

    def _store(self, uid: str, outputs: dict[str, Stream]) -> None:
        for port, stream in outputs.items():
            self.snapshot[self.by_source.get((uid, port), f"{uid}.{port}")] = stream

    def _evaluate(self, t: float, x: np.ndarray) -> None:
        """Bring every stream and unit result to time ``t`` and state ``x``."""
        for target in {ev.target for ev in self.events}:
            v = self._event_value(target, t)
            self._set(target, v if v is not None else self.base_values[target])
        for cid, ctrl in self.controllers.items():
            value, _ = ctrl.output(x[self.slices[cid]])
            self._set(_target(ctrl.model.manipulate), value)
        for uid in self.sources:
            self._store(uid, self.dynamic[uid].outlets(x[self.slices[uid]]))
        if not self.tears:
            self._sweep_once(x)
            return
        tol = self.model.solver.tolerance
        for _ in range(self.model.solver.max_iterations):
            before = {sid: Solver._vector(self.snapshot[sid]) for sid in self.tears}
            self._sweep_once(x)
            if all(Solver._close(before[s], Solver._vector(self.snapshot[s]), tol)
                   for s in self.tears):
                return
        raise SimulationError(f"the recycle did not converge at t = {t:g} s",
                              where=", ".join(sorted(self.tears)))

    def _sweep_once(self, x: np.ndarray) -> None:
        for uid in self.order:
            model = self.dynamic.get(uid)
            if model is not None:
                self._store(uid, model.outlets(x[self.slices[uid]], self._inlets(uid)))
            else:
                self._store(uid, self._calc(uid))

    def _rhs(self, t: float, z: np.ndarray) -> np.ndarray:
        x = z * self.xscale
        self.t_now = t
        check_deadline("dynamics")
        self._evaluate(t, x)
        dx = np.zeros_like(x)
        for uid, model in self.dynamic.items():
            dx[self.slices[uid]] = model.derivatives(x[self.slices[uid]], self._inlets(uid))
        for cid, ctrl in self.controllers.items():
            dx[self.slices[cid]] = ctrl.derivatives(x[self.slices[cid]],
                                                   self._read(ctrl.measure))
        return dx / self.xscale

    # ---- reading values ----------------------------------------------------------------

    def _read(self, ref) -> float:
        if ref.stream is not None:
            s = self.snapshot.get(ref.stream)
            if s is None:
                raise SimulationError(f"stream {ref.stream} is not in the results")
            return _stream_value(s, ref.property, ref.component, self.pkg)
        model = self.dynamic.get(ref.unit)
        if model is not None and ref.property in model.results:
            return float(model.results[ref.property])
        value = self.units[ref.unit].results.get(ref.property)
        if not isinstance(value, (int, float)) or not math.isfinite(value):
            raise SimulationError(f"{ref.unit} reports no {ref.property} now")
        return float(value)

    # ---- integration -------------------------------------------------------------------

    def _sparsity(self):
        """Which state blocks can change which derivatives."""
        owners = list(self.slices)
        manipulators: dict[str, set[str]] = {}
        for cid, ctrl in self.controllers.items():
            manipulators.setdefault(ctrl.model.manipulate.unit, set()).add(cid)
        inlet_srcs: dict[str, list[str]] = {u: [] for u in self.units}
        for c in self.connections:
            if c.dst_unit is not None:
                inlet_srcs[c.dst_unit].append(c.src_unit)
        memo: dict[str, set[str]] = {}

        def affects(uid: str, seen: frozenset = frozenset()) -> set[str]:
            """State blocks that a unit's outlets and results depend on."""
            if uid in memo:
                return memo[uid]
            if uid in seen:
                return set()
            blocks = set(manipulators.get(uid, ()))
            if uid in self.dynamic:
                blocks.add(uid)
            if uid not in self.dynamic or not self.dynamic[uid].source:
                for src in inlet_srcs[uid]:
                    blocks |= affects(src, seen | {uid})
            if not seen:
                memo[uid] = blocks
            return blocks

        n = len(self.x0)
        S = lil_matrix((n, n), dtype=int)
        for owner in owners:
            deps = {owner}
            if owner in self.dynamic:
                deps |= affects(owner)
                for src in inlet_srcs[owner]:
                    deps |= affects(src)
            else:
                ref = self.controllers[owner].measure
                unit = ref.unit if ref.unit is not None else next(
                    (c.src_unit for c in self.connections if c.id == ref.stream),
                    ref.stream.split(".")[0])
                deps |= affects(unit)
                for src in inlet_srcs[unit]:
                    deps |= affects(src)
            rows = self.slices[owner]
            for d in deps:
                cols = self.slices[d]
                S[rows, cols] = 1
        return S.tocsr()

    def _integrate(self, out: dict) -> None:
        dyn = self.model.dynamics
        grid = np.linspace(0.0, self.t_end, dyn.points)
        breaks = {0.0, self.t_end}
        for ev in self.events:
            for t in (ev.t, ev.t + ev.ramp):
                if 0.0 < t < self.t_end:
                    breaks.add(t)
        breaks = sorted(breaks)
        z = self.x0 / self.xscale
        times: list[float] = []
        states: list[np.ndarray] = []
        try:
            for a, b in itertools.pairwise(breaks):
                pts = [t for t in grid if a <= t <= b and (not times or t > times[-1])]
                sol = solve_ivp(self._rhs, (a, b), z, method="BDF",
                                t_eval=np.unique(np.r_[pts, b]), rtol=dyn.rtol,
                                atol=dyn.rtol * 1e-2, jac_sparsity=self.sparsity,
                                first_step=min(1.0, (b - a) / 100))
                wanted = set(pts)
                for t, y in zip(sol.t, sol.y.T):
                    if t in wanted:
                        times.append(float(t))
                        states.append(y)
                if sol.status < 0:
                    raise SimulationError(
                        f"the integration stopped at t = {sol.t[-1] if sol.t.size else a:g} s: "
                        f"{sol.message}",
                        hint="a value may run out of its range (a vessel empties or "
                        "overfills); check the controllers and events")
                z = sol.y[:, -1]
        except SimulationError as exc:
            self._record(out, times, states)
            raise exc.located(exc.where or f"t = {getattr(self, 't_now', 0.0):g} s")
        self._record(out, times, states)

    def _record(self, out: dict, times: list[float], states: list[np.ndarray]) -> None:
        series: dict[str, dict] = {}

        def add(key, group, owner, label, quantity, value):
            entry = series.setdefault(key, {"key": key, "group": group, "owner": owner,
                                            "label": label, "quantity": quantity,
                                            "values": []})
            entry["values"].append(value)

        for t, z in zip(times, states):
            x = z * self.xscale
            self._evaluate(t, x)
            for c in self.connections:
                s = self.snapshot[c.id]
                add(f"stream.{c.id}.T", "Streams", c.id, "Temperature", "temperature", s.T)
                add(f"stream.{c.id}.molar_flow", "Streams", c.id, "Molar flow", "molar_flow",
                    s.total_flow)
                for name, zi in zip(self.pkg.names, s.zs):
                    add(f"stream.{c.id}.x.{name}", "Streams", c.id, f"Mole fraction {name}",
                        "fraction", zi)
            for uid, model in self.dynamic.items():
                for f in model.result_fields:
                    if f.name in model.results:
                        add(f"unit.{uid}.{f.name}", "Units", uid, f.label, f.quantity,
                            float(model.results[f.name]))
            for cid, ctrl in self.controllers.items():
                op, _ = ctrl.output(x[self.slices[cid]])
                add(f"controller.{cid}.pv", "Controllers", cid, "Measurement", ctrl.quantity,
                    float(x[self.slices[cid]][1]))
                add(f"controller.{cid}.sp", "Controllers", cid, "Setpoint", ctrl.quantity,
                    ctrl.setpoint)
                add(f"controller.{cid}.op", "Controllers", cid, "Output", ctrl.op_quantity, op)
        out["time"] = list(times)
        out["series"] = list(series.values())


def to_csv(out: dict, keys: list[str] | None = None) -> str:
    """The time series as CSV: time in s, then one column per series (SI units)."""
    series = [s for s in out["series"] if keys is None or s["key"] in keys]
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(["time [s]", *(_column(s) for s in series)])
    for i, t in enumerate(out["time"]):
        writer.writerow([f"{t:.10g}", *(f"{s['values'][i]:.10g}" for s in series)])
    return buf.getvalue()


def _column(s: dict) -> str:
    unit = SI_UNITS.get(s["quantity"] or "", "")
    return f"{s['owner']} {s['label']}" + (f" [{unit}]" if unit else "")


SI_UNITS = {"temperature": "K", "pressure": "Pa", "molar_flow": "mol/s", "power": "W",
            "length": "m", "volume": "m3", "amount": "mol", "time": "s", "mass_flow": "kg/s"}


def _target(m: ManipulatedParam) -> tuple:
    return ("unit", m.unit, m.param, m.component)


def _si(value, quantity: str | None, where: str, what: str) -> float:
    try:
        return to_si(value, quantity) if quantity else float(value)
    except (UnitError, ValueError) as exc:
        raise SimulationError(f"{what}: {exc}", where=where) from exc


def _stream_value(s: Stream, prop: str, component: str | None, pkg) -> float:
    if prop not in STREAM_PROPERTIES:
        raise SimulationError(f"streams have no property {prop!r}")
    if prop == "T":
        return s.T
    if prop == "P":
        return s.P
    if prop == "vapor_fraction":
        return s.state.VF
    if prop == "molar_flow":
        return s.total_flow
    if prop == "mass_flow":
        return sum(f * mw for f, mw in zip(s.flows, pkg.mws)) / 1000.0
    if prop == "enthalpy_flow":
        return s.enthalpy_flow
    i = pkg.names.index(component)
    return s.zs[i] if prop == "mole_fraction" else s.flows[i]
