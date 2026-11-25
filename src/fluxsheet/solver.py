"""Sequential-modular flowsheet solver.

Units are calculated one after another in flow order. Recycle loops are found
automatically: the solver picks tear streams, guesses them, and iterates until
the guess matches what comes back around the loop (Wegstein or Broyden acceleration).
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from math import isfinite
from typing import Any

import numpy as np

from . import utilities
from .errors import SimulationError
from .flowsheet import CirculationModel, FlowsheetModel, RecycleGuessModel
from .limits import check_deadline
from .property_sections import PropertySections
from .quantities import UnitError, to_si
from .reactions import parse_reactions
from .stream import Stream
from .thermo_package import METHODS, PropertyPackage, get_package, package_options, suggest_method
from .units import UNIT_TYPES, UnitOperation
from .units.base import mix_streams

# These unit sides exchange heat across a wall, so the ports share calculation
# dependencies while their material paths stay separate.
_SEPARATE_SIDES = {
    "HeatExchanger": {"hot_out": "hot_in", "cold_out": "cold_in"},
    "PFR": {"out": "in", "coolant_out": "coolant_in"},
}


@dataclass
class _Connection:
    id: str
    src_unit: str
    src_port: str
    dst_unit: str | None
    dst_port: str | None
    circulation: CirculationModel | None
    recycle_guess: RecycleGuessModel | None


def _mass_flow(flows, mws) -> float:
    return sum(flow * mw for flow, mw in zip(flows, mws)) / 1000.0


def _ancestors(predecessors: dict[str, set[str]], targets: set[str]) -> set[str]:
    seen, todo = set(targets), list(targets)
    while todo:
        for upstream in predecessors[todo.pop()] - seen:
            seen.add(upstream)
            todo.append(upstream)
    return seen


@dataclass
class _Circulation:
    mass_flow: float
    mole_fractions: list[float] | None

    def flows_for(self, flows: list[float], mws) -> list[float]:
        basis = self.mole_fractions if self.mole_fractions is not None else flows
        scale = self.mass_flow / _mass_flow(basis, mws)
        return [flow * scale for flow in basis]


@dataclass
class SolveResult:
    ok: bool
    streams: dict[str, dict] = field(default_factory=dict)
    units: dict[str, dict] = field(default_factory=dict)
    property_method: str | None = None
    property_method_reason: str | None = None
    # Where each binary parameter set came from: ThermoML fit, ChemSep or UNIFAC estimate.
    parameter_sources: list[dict] = field(default_factory=list)
    # Streams re-flashed where they pass from one property method to another.
    property_crossings: list[dict] = field(default_factory=list)
    balance_checks: list[dict] = field(default_factory=list)
    tear_streams: list[str] = field(default_factory=list)
    iterations: int = 0
    converged: bool = False
    warnings: list[str] = field(default_factory=list)
    error: dict | None = None
    design_specs: list[dict] = field(default_factory=list)  # filled by specs.py
    calculators: list[dict] = field(default_factory=list)  # filled by calculator.py
    costing: dict | None = None  # filled by costing.py
    specification_checks: list[dict] = field(default_factory=list)
    specification_status: str = "not_run"
    balance_status: str = "not_run"

    def to_dict(self) -> dict[str, Any]:
        return self.__dict__.copy()


def _check_status(checks: list[dict]) -> str:
    """Summarize a completed check pass, prioritizing a known violation."""
    statuses = {check["status"] for check in checks}
    if "warning" in statuses:
        return "warning"
    if "not_checkable" in statuses:
        return "not_checkable"
    return "ok" if checks else "not_applicable"


class Solver:
    def __init__(self, model: FlowsheetModel):
        self.model = model
        self.units: dict[str, UnitOperation] = {}
        self.connections: list[_Connection] = []
        self.warnings: list[str] = []
        self.balance_io: dict[str, tuple[list[Stream], list[Stream]]] = {}

    # ---- setup ---------------------------------------------------------------------------

    def _build(self) -> None:
        for um in self.model.units:
            if um.id in self.units:
                raise SimulationError("two units share this id", where=um.id,
                                      hint="give every unit a unique id")
            cls = UNIT_TYPES.get(um.type)
            if cls is None:
                raise SimulationError(f"unknown unit type {um.type!r}", where=um.id,
                                      hint=f"available types: {', '.join(UNIT_TYPES)}")
            self.units[um.id] = cls(um.id, dict(um.params))

        seen_ids: set[str] = set()
        used: dict[tuple[str, str], str] = {}
        for sm in self.model.streams:
            if sm.id in seen_ids:
                raise SimulationError("two streams share this id", where=sm.id,
                                      hint="give every stream a unique id")
            seen_ids.add(sm.id)
            src_unit, src_port = sm.source.split(".")
            dst_unit, dst_port = sm.target.split(".") if sm.target else (None, None)
            self._check_endpoint(sm.id, src_unit, src_port, "outlets")
            if dst_unit is not None:
                self._check_endpoint(sm.id, dst_unit, dst_port, "inlets")
            if (src_unit, src_port) in used:
                raise SimulationError(
                    f"outlet {sm.source} is already connected to stream "
                    f"{used[(src_unit, src_port)]}", where=sm.id,
                    hint="an outlet feeds exactly one stream; add a Splitter to divide it",
                )
            used[(src_unit, src_port)] = sm.id
            self.connections.append(_Connection(sm.id, src_unit, src_port, dst_unit, dst_port,
                                                 sm.circulation, sm.recycle_guess))

        for unit in self.units.values():
            for port in unit.inlets:
                n = sum(1 for c in self.connections
                        if c.dst_unit == unit.id and c.dst_port == port.name)
                if n == 0:
                    raise SimulationError(f"inlet {port.name!r} is not connected", where=unit.id,
                                          hint="connect a stream to it, or delete the unit")
                if n > 1 and not port.multiple:
                    raise SimulationError(
                        f"inlet {port.name!r} has {n} streams but accepts only one",
                        where=unit.id, hint="put a Mixer in front of it",
                    )
            for port in unit.outlets:
                if (unit.id, port.name) not in used and port.name not in _idle_outlets(unit):
                    self.warnings.append(f"{unit.id}: outlet {port.name!r} is not connected")

    def _check_endpoint(self, sid: str, unit_id: str, port: str, side: str) -> None:
        unit = self.units.get(unit_id)
        if unit is None:
            raise SimulationError(f"connects to unknown unit {unit_id!r}", where=sid)
        names = [p.name for p in getattr(unit, side)]
        if port not in names:
            kind = "outlet" if side == "outlets" else "inlet"
            raise SimulationError(
                f"{unit_id} has no {kind} called {port!r}", where=sid,
                hint=f"{unit.label} {side}: {', '.join(names) or 'none'}",
            )

    def _package(self) -> tuple[PropertyPackage, str | None]:
        method, reason = self.model.property_method, None
        options = package_options(self.model)
        if method == "auto":
            suggestion = suggest_method(self.model.components, options["custom"])
            method, reason = suggestion.method, suggestion.reason
        elif method not in METHODS:
            raise SimulationError(f"unknown property method {method!r}",
                                  hint=f"use 'auto' or one of: {', '.join(METHODS)}")
        return get_package(tuple(self.model.components), method, **options), reason

    # ---- ordering ------------------------------------------------------------------------

    def _order(self) -> tuple[list[str], set[str]]:
        """Calculation order, plus the stream ids torn to break recycle loops."""
        out_edges: dict[str, list[_Connection]] = {u: [] for u in self.units}
        for c in self.connections:
            if c.dst_unit is not None:
                out_edges[c.src_unit].append(c)

        preferred = {c.id for c in self.connections
                     if c.circulation is not None or c.recycle_guess is not None}
        for c in self.connections:
            if c.id in preferred and c.dst_unit not in self._upstream_of({c.src_unit}):
                raise SimulationError("stream is not part of a recycle loop", where=c.id,
                                      hint="use circulation and recycle guesses on a stream "
                                      "with a path back from its destination to its source")

        tears = set(preferred)
        state: dict[str, int] = {}  # 1 = on stack, 2 = done

        def visit(u: str) -> None:
            state[u] = 1
            for c in out_edges[u]:
                if c.id in preferred:
                    continue
                if state.get(c.dst_unit) == 1:
                    tears.add(c.id)  # back edge closes a loop
                elif c.dst_unit not in state:
                    visit(c.dst_unit)
            state[u] = 2

        roots = [u for u, unit in self.units.items() if unit.type == "Feed"]
        for u in roots + list(self.units):
            if u not in state:
                visit(u)

        if not roots and not tears:
            raise SimulationError("the flowsheet has no feed",
                                  hint="add a Feed to bring material into the process")
        indegree = {u: 0 for u in self.units}
        for c in self.connections:
            if c.dst_unit is not None and c.id not in tears:
                indegree[c.dst_unit] += 1
        ready = [u for u in self.units if indegree[u] == 0]
        order: list[str] = []
        while ready:
            u = ready.pop(0)
            order.append(u)
            for c in out_edges[u]:
                if c.id in tears:
                    continue
                indegree[c.dst_unit] -= 1
                if indegree[c.dst_unit] == 0:
                    ready.append(c.dst_unit)
        return order, tears

    def _recycle_material_ancestors(self, tears: set[str]) -> dict[str, set[str]]:
        if not tears:
            return {}
        incoming: dict[str, list[_Connection]] = {uid: [] for uid in self.units}
        for c in self.connections:
            if c.dst_unit is not None:
                incoming[c.dst_unit].append(c)
        predecessors = {}
        sources = set()
        specified = set()
        for c in self.connections:
            unit = self.units[c.src_unit]
            side = _SEPARATE_SIDES.get(unit.type, {}).get(c.src_port)
            predecessors[c.id] = {s.id for s in incoming[c.src_unit]
                                  if side is None or s.dst_port == side}
            if unit.type == "Feed":
                sources.add(c.id)
            if c.circulation is not None:
                specified.add(c.id)

        anchors = sources | specified
        ancestors = {sid: _ancestors(predecessors, predecessors[sid]) for sid in sorted(tears)}
        for sid, upstream in ancestors.items():
            if sid in specified and sid not in upstream:
                raise SimulationError(
                    "circulation specification requires a material recycle", where=sid,
                    hint="use recycle_guess on a thermal return with no material circulation",
                )
            if sid in upstream and anchors.isdisjoint(upstream):
                raise SimulationError(
                    "closed recycle has no circulation mass-flow specification", where=sid,
                    hint="set circulation.mass_flow on an internal recycle stream to define "
                    "the circulating mass flow. recycle_guess supplies an initial estimate",
                )
        return ancestors

    # ---- solving -------------------------------------------------------------------------

    def solve(self) -> SolveResult:
        result = SolveResult(ok=False)
        try:
            self._build()
            order, tears = self._order()
            result.tear_streams = sorted(tears)
            material_upstream = self._recycle_material_ancestors(tears)
            pkg, reason = self._package()
            result.property_method = pkg.method
            result.property_method_reason = reason
            result.parameter_sources = list(pkg.parameter_sources)
            self.warnings.extend(pkg.warnings)  # e.g. binary parameters estimated, not measured
            self.sections = PropertySections(self.model, pkg, self.units, self.connections)
            self.warnings.extend(self.sections.warnings())
            for method, other in self.sections.packages.items():  # per-unit methods
                if other is not pkg:
                    result.parameter_sources.extend({**s, "method": method}
                                                    for s in other.parameter_sources)
            reactions = parse_reactions(self.model.reactions, pkg)
            for unit in self.units.values():
                if hasattr(unit, "reactions"):  # reactors run the flowsheet's reactions
                    unit.reactions = reactions
            media = utilities.parse_utilities(self.model.utilities)
            for unit in self.units.values():
                if "utility" in unit.p:  # units that can be heated or cooled by a utility
                    utilities.attach(unit, media)
            streams = self._iterate(pkg, order, tears, material_upstream, result)
            for uid, unit in self.units.items():  # from the converged calculation only
                self.warnings.extend(f"{uid}: {note}" for note in unit.notes)
            self.pkg, self.streams = pkg, streams  # the solution, for dynamics/
        except SimulationError as exc:
            result.error = exc.to_dict()
            for uid, unit in self.units.items():  # often the reason, e.g. a surging compressor
                self.warnings.extend(f"{uid}: {note} (in the last, failed calculation)"
                                     for note in unit.notes)
            result.warnings = self.warnings
            return result

        result.ok = True
        ordered = [c.id for c in self.connections] + [k for k in streams if k not in
                                                       {c.id for c in self.connections}]
        result.streams = {sid: streams[sid].to_dict(pkg) for sid in ordered}
        result.units = {
            uid: {"type": u.type, **u.results, **({"profiles": u.profiles} if u.profiles else {})}
            for uid, u in self.units.items()
        }
        self._check_balances(pkg, streams)
        states = {k: s.state for k, s in streams.items()}
        for p in self.sections.packages.values():  # the default and per-unit methods
            self.warnings.extend(p.henry_validity(states))
        result.warnings = self.warnings
        self.sections.annotate(result, streams, self.connections)
        from .balances import check_balances, check_specifications

        result.balance_checks = check_balances(self, pkg, streams, result.property_crossings)
        result.balance_status = _check_status(result.balance_checks)
        result.specification_checks = check_specifications(self, pkg, streams)
        result.specification_status = _check_status(result.specification_checks)
        for check in result.balance_checks:
            if check["status"] == "warning":
                name = f" ({check['component']})" if check["component"] else ""
                self.warnings.append(
                    f"{check['scope']}: {check['quantity']}{name} balance residual "
                    f"{check['residual']:.6g} {check['unit']} exceeds "
                    f"{check['tolerance']:.6g} {check['unit']}")
        for check in result.specification_checks:
            if check["status"] == "warning":
                unit_of_measure = check["unit"]
                self.warnings.append(
                    f"{check['scope']}: {check['side']} duty specification "
                    f"target {check['target']:.6g} {unit_of_measure}, "
                    f"actual {check['actual']:.6g} {unit_of_measure}, "
                    f"residual {check['residual']:.6g} {unit_of_measure} "
                    f"(tolerance {check['tolerance']:.6g} {unit_of_measure})")
        if self.model.costing is not None:
            from .costing import estimate  # sizes and prices the solved equipment

            try:
                estimate(self.model, result, self.sections.package)
            except SimulationError as exc:
                result.ok, result.error = False, exc.to_dict()
        return result

    def _run_unit(self, pkg, unit: UnitOperation, streams: dict[str, Stream]) -> dict[str, Stream]:
        check_deadline(unit.id)  # demo server time budget; no-op otherwise
        unit_pkg = self.sections.package(unit.id)
        inlets: dict[str, list[Stream]] = {p.name: [] for p in unit.inlets}
        unit.inlet_ids = {p.name: [] for p in unit.inlets}
        for c in self.connections:
            if c.dst_unit == unit.id:  # re-flashed if it comes from another property method
                inlets[c.dst_port].append(self.sections.inlet(c, streams[c.id]))
                unit.inlet_ids[c.dst_port].append(c.id)
        unit.notes = []
        try:
            for port in unit.inlets:  # several feeds on one inlet enter as one mixed stream
                if len(inlets[port.name]) > 1 and not port.separate:
                    inlets[port.name] = [mix_streams(unit_pkg, inlets[port.name])]
            outputs = utilities.calculate(unit, unit_pkg, inlets)
        except SimulationError as exc:
            raise exc.located(unit.id)
        except Exception as exc:  # keep unexpected failures readable too
            raise SimulationError(f"calculation failed: {exc}", where=unit.id) from exc
        self.balance_io[unit.id] = ([s for group in inlets.values() for s in group],
                                    list(outputs.values()))
        self._check_unit(unit, inlets)
        return outputs

    def _check_unit(self, unit: UnitOperation, inlets: dict[str, list[Stream]]) -> None:
        if unit.type == "Pump" and inlets["in"][0].state.VF > 1e-6:
            self.warnings.append(f"{unit.id}: pump inlet contains vapor; it may cavitate")
        if unit.type == "Compressor" and inlets["in"][0].state.VF < 1.0 - 1e-6:
            self.warnings.append(f"{unit.id}: compressor inlet contains liquid")

    def _sweep(self, pkg, order, streams, tears=frozenset()) -> None:
        by_source = {(c.src_unit, c.src_port): c.id for c in self.connections}
        # Every consuming unit sees the original guess, even if the tear's producer ran first.
        returned = {}
        for uid in order:
            unit = self.units[uid]
            for port, stream in self._run_unit(pkg, unit, streams).items():
                sid = by_source.get((uid, port))
                if sid is None and getattr(unit, "drop_unconnected", False):
                    continue  # e.g. an unused copy of a duplicator
                if sid is None:
                    if port in _idle_outlets(unit):
                        if stream.total_flow <= 0.0:
                            continue
                        note = f"{uid}: outlet {port!r} carries flow but is not connected"
                        if note not in self.warnings:  # once, not per recycle iteration
                            self.warnings.append(note)
                    sid = f"{uid}.{port}"
                if sid in tears:
                    returned[sid] = stream
                else:
                    streams[sid] = stream
        streams.update(returned)

    @staticmethod
    def _recycle_quantity(value, quantity: str, sid: str, field: str, *, zero=False) -> float:
        try:
            number = to_si(value, quantity)
        except UnitError as exc:
            raise SimulationError(f"invalid {field}: {exc}", where=sid) from exc
        if not isfinite(number) or number < 0 or (number == 0 and not zero):
            bound = "nonnegative" if zero else "positive"
            raise SimulationError(f"{field} must be finite and {bound}", where=sid)
        return number

    def _recycle_defaults(self, source: str, fields: dict[str, float]) -> dict[str, float]:
        """Find the nearest specified T and P upstream of a tear without a feed state."""
        pending = deque([source])
        seen = {source}
        values = {}
        while pending and len(values) < len(fields):
            uid = pending.popleft()
            for key in fields.keys() - values.keys():
                value = self.units[uid].p.get(key)
                if isinstance(value, (int, float)) and value > 0:
                    values[key] = value
            for c in self.connections:
                if c.dst_unit == uid and c.src_unit not in seen:
                    seen.add(c.src_unit)
                    pending.append(c.src_unit)
        return {**fields, **values}

    def _initialize_tears(self, pkg, tears, src_pkg, material_upstream):
        streams: dict[str, Stream] = {}
        circulations: dict[str, _Circulation] = {}
        feed_units = {uid for uid, unit in self.units.items() if unit.type == "Feed"}
        feeds = [unit.p for uid, unit in self.units.items() if uid in feed_units]
        feed_state = {"T": (feeds[0]["T"] or 298.15) if feeds else 298.15,
                      "P": max((p["P"] for p in feeds), default=1e5)}
        composition_sources = {
            c.id for c in self.connections if c.src_unit in feed_units
            or (c.circulation is not None and c.circulation.mole_fractions is not None)
        }
        for c in self.connections:
            sid = c.id
            if sid not in tears:
                continue
            spec, guess = c.circulation, c.recycle_guess
            flows = [0.0] * pkg.n
            state = feed_state.copy()
            if not feeds or spec is not None:
                state = self._recycle_defaults(c.src_unit, state)
            if spec is not None:
                mass = self._recycle_quantity(spec.mass_flow, "mass_flow", sid,
                                              "circulation.mass_flow")
                fractions = spec.mole_fractions
                composition = None
                if fractions is not None:
                    values = self._named_recycle_values(pkg, sid, "circulation.mole_fractions",
                                                        fractions, "fraction")
                    if abs(sum(values) - 1.0) > 1e-8:
                        raise SimulationError("circulation mole fractions must sum to one",
                                              where=sid, hint="enter mole fractions from 0 to 1")
                    composition = values
                    flows = values
                elif pkg.n > 1 and composition_sources.isdisjoint(material_upstream[sid]):
                    raise SimulationError(
                        "closed recycle has no circulating composition specification",
                        where=sid, hint="set circulation.mole_fractions for the mixture",
                    )
                else:
                    flows = [1.0] * pkg.n
                circulation = circulations[sid] = _Circulation(mass, composition)
                flows = circulation.flows_for(flows, pkg.mws)
            if guess is not None:
                for field, quantity in (("T", "temperature"), ("P", "pressure")):
                    value = getattr(guess, field)
                    if value is not None:
                        state[field] = self._recycle_quantity(value, quantity, sid,
                                                              f"recycle_guess.{field}")
                if guess.flows is not None:
                    flows = self._named_recycle_values(pkg, sid, "recycle_guess.flows",
                                                       guess.flows, "molar_flow")
                    if spec is not None and sum(flows) == 0:
                        raise SimulationError("a specified circulation needs a positive start flow",
                                              where=sid, hint="give a positive recycle_guess.flows "
                                              "or use the default from the circulation specification")
            if not isfinite(sum(flows)):
                raise SimulationError("recycle start flow is too large", where=sid,
                                      hint="check the circulation and start flow units")
            try:
                streams[sid] = Stream.from_flash(src_pkg[sid], flows, **state)
            except SimulationError as exc:
                raise SimulationError(f"recycle start flash failed: {exc.message}", where=sid,
                                      hint="adjust recycle_guess.T, recycle_guess.P or "
                                      "recycle_guess.flows") from exc
        return streams, circulations

    def _named_recycle_values(self, pkg, sid, field, values, quantity) -> list[float]:
        unknown = values.keys() - set(pkg.names)
        if unknown:
            raise SimulationError(f"{field} names unknown components: {', '.join(sorted(unknown))}",
                                  where=sid, hint="use the flowsheet's component names")
        return [self._recycle_quantity(values.get(name, 0.0), quantity, sid, field, zero=True)
                for name in pkg.names]

    def _project_circulation(self, pkg, spec, stream, previous, sid):
        basis = stream.flows if stream.total_flow > 0 else previous.flows
        flows = spec.flows_for(basis, pkg.mws)
        total = sum(flows)
        if all(abs(f / total - z) < 1e-14 for f, z in zip(flows, stream.zs)):
            return Stream(flows, stream.state)
        try:
            return Stream.from_flash(pkg, flows, P=stream.P, H=stream.state.H)
        except SimulationError as exc:
            raise SimulationError(f"circulation iteration flash failed: {exc.message}", where=sid,
                                  hint="check the specified composition and recycle guess") from exc

    def _check_circulation_state(self, sid: str, upstream: set[str]) -> None:
        units = [self.units[uid] for uid in {c.src_unit for c in self.connections if c.id in upstream}]
        # These units only shift or pass on inlet pressure and either set outlet temperature
        # or shift molar enthalpy at fixed circulation. A uniform shift can survive a sweep
        # unchanged when no unit supplies the corresponding absolute condition.
        simple_units = {"Heater", "Valve", "Mixer", "Splitter"}
        if any(unit.type not in simple_units or unit.p.get("area") is not None for unit in units):
            return
        if len({self.sections.methods[unit.id] for unit in units}) > 1:
            return  # a change in property method can change the enthalpy map
        if all(unit.p.get("P") is None for unit in units):
            raise SimulationError(
                "closed recycle has no absolute pressure specification", where=sid,
                hint="set an absolute outlet pressure on a unit in the loop; "
                "recycle_guess.P only starts the iteration",
            )
        if all(unit.type != "Heater" or unit.p["duty"] is not None for unit in units):
            raise SimulationError(
                "closed recycle has no temperature (or enthalpy) specification", where=sid,
                hint="specify a unit's outlet temperature or vapor fraction; "
                "recycle_guess.T only starts the iteration",
            )

    def _iterate(self, pkg, order, tears, material_upstream, result: SolveResult) -> dict[str, Stream]:
        streams: dict[str, Stream] = {}
        if not tears:
            self._sweep(pkg, order, streams)
            result.iterations, result.converged = 1, True
            return streams

        src_pkg = {c.id: self.sections.package(c.src_unit) for c in self.connections if c.id in tears}
        streams, circulations = self._initialize_tears(pkg, tears, src_pkg, material_upstream)

        # Only units that feed back into a tear stream take part in the iterations. Units
        # downstream of the loop (e.g. a column on the crude product) run once at the end, so
        # they never see the unconverged intermediate states.
        loop = self._upstream_of({c.src_unit for c in self.connections if c.id in tears})
        loop_order = [u for u in order if u in loop]

        opts = self.model.solver
        torn = sorted(tears)
        broyden = _Broyden(pkg.n + 2) if opts.acceleration == "broyden" else None
        prev_x: dict[str, list[float]] = {}
        prev_g: dict[str, list[float]] = {}
        for iteration in range(1, opts.max_iterations + 1):
            previous = {sid: streams[sid] for sid in circulations}
            x = {sid: self._vector(streams[sid]) for sid in tears}
            self._sweep(pkg, loop_order, streams, tears)
            returned = {sid: streams[sid] for sid in circulations}
            for sid, spec in circulations.items():
                streams[sid] = self._project_circulation(src_pkg[sid], spec, returned[sid],
                                                         previous[sid], sid)
            g = {sid: self._vector(streams[sid]) for sid in tears}

            if all(self._close(x[s], g[s], opts.tolerance) for s in tears):
                # Verify the real return before replacing projected iterates with results.
                for sid, spec in circulations.items():
                    actual = _mass_flow(returned[sid].flows, pkg.mws)
                    if abs(actual - spec.mass_flow) > opts.tolerance * spec.mass_flow:
                        result.iterations = iteration
                        raise SimulationError(
                            f"circulation mass flow does not close (target {spec.mass_flow:.6g} "
                            f"kg/s, returned {actual:.6g} kg/s)", where=sid,
                            hint="check that the loop returns the specified flow and that any "
                            "makeup and losses balance",
                        )
                    if not self._close(g[sid], self._vector(returned[sid]), opts.tolerance):
                        result.iterations = iteration
                        raise SimulationError("circulating composition does not close", where=sid,
                                              hint="check the composition specification and "
                                              "the component balances around the loop")
                result.iterations = iteration
                for sid in circulations:
                    self._check_circulation_state(sid, material_upstream[sid])
                streams.update(returned)
                result.converged = True
                self._sweep(pkg, [u for u in order if u not in loop], streams)
                return streams

            new: dict[str, list[float]] = {}
            if opts.acceleration == "wegstein" and iteration >= 3 and prev_x:
                new = {sid: self._wegstein(prev_x[sid], prev_g[sid], x[sid], g[sid])
                       for sid in tears}
            elif broyden is not None:
                flat = broyden.step([v for s in torn for v in x[s]],
                                    [v for s in torn for v in g[s]])
                size = pkg.n + 2
                new = {s: flat[k * size:(k + 1) * size] for k, s in enumerate(torn)}
            for sid, vec in new.items():
                # Keep the physical sweep result when an accelerated guess is invalid.
                if not all(np.isfinite(v) for v in vec) or min(vec[:-2]) < 0.0 or vec[-1] <= 0:
                    continue
                if sum(vec[:-2]) <= 1e-14:
                    continue
                spec = circulations.get(sid)
                if spec is not None:
                    vec = [*spec.flows_for(vec[:-2], pkg.mws), *vec[-2:]]
                try:
                    streams[sid] = Stream.from_flash(src_pkg[sid], vec[:-2], H=vec[-2], P=vec[-1])
                except (SimulationError, ValueError, OverflowError):
                    pass
            prev_x, prev_g = x, g

        result.iterations = opts.max_iterations
        raise SimulationError(
            f"recycle did not converge in {opts.max_iterations} iterations",
            where=", ".join(sorted(tears)),
            hint=("check the circulation specifications and starting values, or raise "
                  "solver.max_iterations" if circulations else "check that every component has "
                  "a way out of the loop (a purge or product), or raise solver.max_iterations"),
        )

    def _upstream_of(self, targets: set[str]) -> set[str]:
        """``targets`` plus every unit with a path of streams leading to one of them."""
        feeds_into: dict[str, set[str]] = {u: set() for u in self.units}
        for c in self.connections:
            if c.dst_unit is not None:
                feeds_into[c.dst_unit].add(c.src_unit)
        return _ancestors(feeds_into, targets)

    @staticmethod
    def _vector(s: Stream) -> list[float]:
        return [*s.flows, s.state.H, s.P]

    @staticmethod
    def _close(a: list[float], b: list[float], tol: float) -> bool:
        scale = max(sum(abs(v) for v in a[:-2]), 1e-10)
        flows_ok = all(abs(x - y) <= tol * scale + 1e-12 for x, y in zip(a[:-2], b[:-2]))
        if not flows_ok:
            return False
        if max(sum(a[:-2]), sum(b[:-2])) <= 1e-14:
            return True
        enthalpy_ok = abs(a[-2] - b[-2]) <= max(tol * max(abs(a[-2]), abs(b[-2]), 1e3), 1e-4)
        pressure_ok = abs(a[-1] - b[-1]) <= max(tol * max(abs(a[-1]), abs(b[-1])), 1e-3)
        return enthalpy_ok and pressure_ok

    @staticmethod
    def _wegstein(x0, g0, x1, g1) -> list[float]:
        out = []
        for a0, b0, a1, b1 in zip(x0, g0, x1, g1):
            dx = a1 - a0
            if abs(dx) < 1e-14:
                out.append(b1)
                continue
            s = (b1 - b0) / dx
            q = s / (s - 1.0) if abs(s - 1.0) > 1e-12 else 0.0
            q = min(max(q, -5.0), 0.0)
            out.append(q * a1 + (1.0 - q) * b1)
        return out

    def _check_balances(self, pkg: PropertyPackage, streams: dict[str, Stream]) -> None:
        def mass(s: Stream) -> float:
            return sum(f * mw for f, mw in zip(s.flows, pkg.mws))

        # Units that do not conserve mass by design (logic blocks such as a multiplier or a
        # duplicator) are boundaries: what enters them leaves the balance, what they put out
        # enters it, as at a Product and a Feed.
        def source(uid: str) -> bool:
            unit = self.units.get(uid)
            return unit is not None and (unit.type == "Feed"
                                         or not getattr(unit, "conserves_mass", True))

        def sink(uid: str) -> bool:
            unit = self.units[uid]
            return unit.type == "Product" or not getattr(unit, "conserves_mass", True)

        by_id = {c.id: c for c in self.connections}
        feeds_in = leaving = 0.0
        for sid, s in streams.items():
            conn = by_id.get(sid)
            if source(conn.src_unit if conn else sid.rsplit(".", 1)[0]):
                feeds_in += mass(s)
            if conn is None or conn.dst_unit is None or sink(conn.dst_unit):
                leaving += mass(s)
        if feeds_in > 0 and abs(feeds_in - leaving) > 1e-5 * feeds_in:
            self.warnings.append(
                f"mass balance does not close: in {feeds_in / 1000:.6g} kg/s, "
                f"out {leaving / 1000:.6g} kg/s"
            )


class _Broyden:
    """Broyden's method on the tear residual f(x) = g(x) - x, with the inverse Jacobian
    updated by Sherman-Morrison ("good" Broyden). It starts as direct substitution and resets
    to it whenever a step makes the residual grow a lot, so it never does worse for long.

    Each tear stream contributes component flows, molar enthalpy and pressure.
    Scaling keeps the variables comparable.
    """

    def __init__(self, block: int) -> None:
        self.block = block  # numbers per tear stream: its component flows, then T
        self.H: np.ndarray | None = None  # approximate inverse Jacobian of f, scaled
        self.prev: tuple[np.ndarray, np.ndarray] | None = None
        self.scale: np.ndarray | None = None

    def step(self, x: list[float], g: list[float]) -> list[float]:
        xa, ga = np.asarray(x, float), np.asarray(g, float)
        if self.scale is None:
            parts = []
            for k in range(0, len(ga), self.block):
                flows = ga[k:k + self.block - 2]
                parts += [max(abs(flows).sum(), 1e-10)] * len(flows)
                parts += [max(abs(ga[k + self.block - 2]), 1e3),
                          max(abs(ga[k + self.block - 1]), 1.0)]
            self.scale = np.asarray(parts)
        xs, fs = xa / self.scale, (ga - xa) / self.scale
        n = len(xs)
        if self.H is None:
            self.H = -np.eye(n)
        elif self.prev is not None:
            dx, df = xs - self.prev[0], fs - self.prev[1]
            if np.linalg.norm(fs) > 2.0 * np.linalg.norm(self.prev[1]):
                self.H = -np.eye(n)  # the model went wrong; start again from substitution
            else:
                Hdf = self.H @ df
                denom = dx @ Hdf
                if abs(denom) > 1e-14:
                    self.H += np.outer(dx - Hdf, dx @ self.H) / denom
        self.prev = (xs, fs)
        dx = -self.H @ fs
        # Limit a step to a doubling or halving of any scaled variable.
        limit = np.max(np.abs(dx) / np.maximum(np.abs(xs), 1.0))
        if limit > 1.0:
            dx /= limit
        return list((xs + dx) * self.scale)


def _idle_outlets(unit: UnitOperation) -> tuple[str, ...]:
    """Outlets a unit leaves empty in its configuration (e.g. a total condenser's vapor).
    Units may define ``idle_outlets()``; unconnected idle outlets are neither listed nor
    warned about."""
    idle = getattr(unit, "idle_outlets", None)
    return idle() if idle else ()


def solve(model: FlowsheetModel | dict) -> SolveResult:
    if isinstance(model, dict):
        model = FlowsheetModel.model_validate(model)
    if any(spec.active for spec in model.design_specs or ()):
        from .specs import solve_with_specs  # an outer loop around this solver

        return solve_with_specs(model)
    return solve_flowsheet(model)


def solve_flowsheet(model: FlowsheetModel) -> SolveResult:
    """One flowsheet solve with its calculator blocks, but without the design specs."""
    if any(calc.active for calc in model.calculators or ()):
        from .calculator import solve_with_calculators

        return solve_with_calculators(model)
    return Solver(model).solve()
