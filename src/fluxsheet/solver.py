"""Sequential-modular flowsheet solver.

Units are calculated one after another in flow order. Recycle loops are found
automatically: the solver picks tear streams, guesses them, and iterates until
the guess matches what comes back around the loop (Wegstein acceleration).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .errors import SimulationError
from .flowsheet import FlowsheetModel
from .reactions import parse_reactions
from .stream import Stream
from .thermo_package import METHODS, PropertyPackage, get_package, suggest_method
from .units import UNIT_TYPES, UnitOperation


@dataclass
class _Connection:
    id: str
    src_unit: str
    src_port: str
    dst_unit: str | None
    dst_port: str | None


@dataclass
class SolveResult:
    ok: bool
    streams: dict[str, dict] = field(default_factory=dict)
    units: dict[str, dict] = field(default_factory=dict)
    property_method: str | None = None
    property_method_reason: str | None = None
    tear_streams: list[str] = field(default_factory=list)
    iterations: int = 0
    converged: bool = False
    warnings: list[str] = field(default_factory=list)
    error: dict | None = None

    def to_dict(self) -> dict[str, Any]:
        return self.__dict__.copy()


class Solver:
    def __init__(self, model: FlowsheetModel):
        self.model = model
        self.units: dict[str, UnitOperation] = {}
        self.connections: list[_Connection] = []
        self.warnings: list[str] = []

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
            self.connections.append(_Connection(sm.id, src_unit, src_port, dst_unit, dst_port))

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

        if not any(u.type == "Feed" for u in self.units.values()):
            raise SimulationError("the flowsheet has no feed",
                                  hint="add a Feed to bring material into the process")

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
        if method == "auto":
            suggestion = suggest_method(self.model.components)
            method, reason = suggestion.method, suggestion.reason
        elif method not in METHODS:
            raise SimulationError(f"unknown property method {method!r}",
                                  hint=f"use 'auto' or one of: {', '.join(METHODS)}")
        return get_package(tuple(self.model.components), method), reason

    # ---- ordering ------------------------------------------------------------------------

    def _order(self) -> tuple[list[str], set[str]]:
        """Calculation order, plus the stream ids torn to break recycle loops."""
        out_edges: dict[str, list[_Connection]] = {u: [] for u in self.units}
        for c in self.connections:
            if c.dst_unit is not None:
                out_edges[c.src_unit].append(c)

        tears: set[str] = set()
        state: dict[str, int] = {}  # 1 = on stack, 2 = done

        def visit(u: str) -> None:
            state[u] = 1
            for c in out_edges[u]:
                if state.get(c.dst_unit) == 1:
                    tears.add(c.id)  # back edge closes a loop
                elif c.dst_unit not in state:
                    visit(c.dst_unit)
            state[u] = 2

        roots = [u for u, unit in self.units.items() if unit.type == "Feed"]
        for u in roots + list(self.units):
            if u not in state:
                visit(u)

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

    # ---- solving -------------------------------------------------------------------------

    def solve(self) -> SolveResult:
        result = SolveResult(ok=False)
        try:
            self._build()
            pkg, reason = self._package()
            result.property_method = pkg.method
            result.property_method_reason = reason
            self.warnings.extend(pkg.warnings)  # e.g. binary parameters estimated, not measured
            reactions = parse_reactions(self.model.reactions, pkg)
            for unit in self.units.values():
                if hasattr(unit, "reactions"):  # reactors run the flowsheet's reactions
                    unit.reactions = reactions
            order, tears = self._order()
            result.tear_streams = sorted(tears)
            streams = self._iterate(pkg, order, tears, result)
        except SimulationError as exc:
            result.error = exc.to_dict()
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
        result.warnings = self.warnings
        return result

    def _run_unit(self, pkg, unit: UnitOperation, streams: dict[str, Stream]) -> dict[str, Stream]:
        inlets: dict[str, list[Stream]] = {p.name: [] for p in unit.inlets}
        for c in self.connections:
            if c.dst_unit == unit.id:
                inlets[c.dst_port].append(streams[c.id])
        try:
            outputs = unit.calculate(pkg, inlets)
        except SimulationError as exc:
            raise exc.located(unit.id)
        except Exception as exc:  # keep unexpected failures readable too
            raise SimulationError(f"calculation failed: {exc}", where=unit.id) from exc
        self._check_unit(unit, inlets)
        return outputs

    def _check_unit(self, unit: UnitOperation, inlets: dict[str, list[Stream]]) -> None:
        if unit.type == "Pump" and inlets["in"][0].state.VF > 1e-6:
            self.warnings.append(f"{unit.id}: pump inlet contains vapor; it may cavitate")
        if unit.type == "Compressor" and inlets["in"][0].state.VF < 1.0 - 1e-6:
            self.warnings.append(f"{unit.id}: compressor inlet contains liquid")

    def _sweep(self, pkg, order, streams) -> None:
        by_source = {(c.src_unit, c.src_port): c.id for c in self.connections}
        for uid in order:
            unit = self.units[uid]
            for port, stream in self._run_unit(pkg, unit, streams).items():
                sid = by_source.get((uid, port))
                if sid is None:
                    if port in _idle_outlets(unit):
                        continue
                    sid = f"{uid}.{port}"
                streams[sid] = stream

    def _iterate(self, pkg, order, tears, result: SolveResult) -> dict[str, Stream]:
        streams: dict[str, Stream] = {}
        if not tears:
            self._sweep(pkg, order, streams)
            result.iterations, result.converged = 1, True
            return streams

        feed_T = next(u.p["T"] for u in self.units.values() if u.type == "Feed") or 298.15
        feed_P = max(u.p["P"] for u in self.units.values() if u.type == "Feed")
        for sid in tears:
            streams[sid] = Stream.from_flash(pkg, [0.0] * pkg.n, T=feed_T, P=feed_P)

        # Only units that feed back into a tear stream take part in the iterations. Units
        # downstream of the loop (e.g. a column on the crude product) run once at the end, so
        # they never see the unconverged intermediate states.
        loop = self._upstream_of({c.src_unit for c in self.connections if c.id in tears})
        loop_order = [u for u in order if u in loop]

        opts = self.model.solver
        prev_x: dict[str, list[float]] = {}
        prev_g: dict[str, list[float]] = {}
        for iteration in range(1, opts.max_iterations + 1):
            x = {sid: self._vector(streams[sid]) for sid in tears}
            self._sweep(pkg, loop_order, streams)
            g = {sid: self._vector(streams[sid]) for sid in tears}

            if all(self._close(x[s], g[s], opts.tolerance) for s in tears):
                result.iterations, result.converged = iteration, True
                self._sweep(pkg, [u for u in order if u not in loop], streams)
                return streams

            if opts.acceleration == "wegstein" and iteration >= 3 and prev_x:
                for sid in tears:
                    new = self._wegstein(prev_x[sid], prev_g[sid], x[sid], g[sid])
                    calc = streams[sid]
                    streams[sid] = Stream.from_flash(pkg, new[:-1], T=new[-1], P=calc.P)
            prev_x, prev_g = x, g

        result.iterations = opts.max_iterations
        raise SimulationError(
            f"recycle did not converge in {opts.max_iterations} iterations",
            where=", ".join(sorted(tears)),
            hint="check that every component has a way out of the loop (a purge or product), "
            "or raise solver.max_iterations",
        )

    def _upstream_of(self, targets: set[str]) -> set[str]:
        """``targets`` plus every unit with a path of streams leading to one of them."""
        feeds_into: dict[str, set[str]] = {u: set() for u in self.units}
        for c in self.connections:
            if c.dst_unit is not None:
                feeds_into[c.dst_unit].add(c.src_unit)
        seen, todo = set(targets), list(targets)
        while todo:
            for src in feeds_into[todo.pop()]:
                if src not in seen:
                    seen.add(src)
                    todo.append(src)
        return seen

    @staticmethod
    def _vector(s: Stream) -> list[float]:
        return [*s.flows, s.T]

    @staticmethod
    def _close(a: list[float], b: list[float], tol: float) -> bool:
        scale = max(sum(abs(v) for v in a[:-1]), 1e-10)
        flows_ok = all(abs(x - y) <= tol * scale + 1e-12 for x, y in zip(a[:-1], b[:-1]))
        return flows_ok and abs(a[-1] - b[-1]) <= max(tol * a[-1], 1e-4)

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

        by_id = {c.id: c for c in self.connections}
        feeds_in = sum(mass(streams[c.id]) for c in self.connections
                       if self.units[c.src_unit].type == "Feed")
        leaving = 0.0
        for sid, s in streams.items():
            conn = by_id.get(sid)
            if conn is None or conn.dst_unit is None or self.units[conn.dst_unit].type == "Product":
                leaving += mass(s)
        if feeds_in > 0 and abs(feeds_in - leaving) > 1e-5 * feeds_in:
            self.warnings.append(
                f"mass balance does not close: in {feeds_in / 1000:.6g} kg/s, "
                f"out {leaving / 1000:.6g} kg/s"
            )


def _idle_outlets(unit: UnitOperation) -> tuple[str, ...]:
    """Outlets a unit leaves empty in its configuration (e.g. a total condenser's vapor).
    Units may define ``idle_outlets()``; unconnected idle outlets are neither listed nor
    warned about."""
    idle = getattr(unit, "idle_outlets", None)
    return idle() if idle else ()


def solve(model: FlowsheetModel | dict) -> SolveResult:
    if isinstance(model, dict):
        model = FlowsheetModel.model_validate(model)
    return Solver(model).solve()
