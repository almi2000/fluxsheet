"""Property sections: different property methods in different parts of a flowsheet.

The flowsheet's ``property_method`` is the default. A unit can override it with its own
``property_method``, e.g. NRTL for a methanol/water column in a plant whose high-pressure
synthesis loop uses Peng-Robinson. Each unit then calculates with its own package, so its
flashes and energy balance are consistent within one method.

Where a stream crosses from one method to another, its state is re-established on the
consuming side: component flows, temperature and pressure are kept from the producing
unit, and the stream is flashed again at that T and P with the consumer's method. Mass is
conserved exactly; the enthalpy may jump by the difference between the two methods'
departure functions (both are relative to the ideal gas at 298.15 K, so the jump is small
unless the phase changes). The jump is reported with each crossing, and a warning is added
when the crossing changes the phase noticeably. Keeping T and P (rather than H and P) is what
commercial simulators do at section boundaries: the stream table shows the same temperature
on both sides, and the energy difference is visible instead of hidden in a shifted T.

Product units without an override take the method of the unit that feeds them, so a product
stream is never re-flashed only to be reported.
"""

from __future__ import annotations

import difflib
import re
from dataclasses import dataclass

from thermo.interaction_parameters import IPDB

from .errors import SimulationError
from .flowsheet import FlowsheetModel
from .stream import Stream
from .thermo_package import METHODS, PropertyPackage, Suggestion, get_package, suggest_method
from .units import UnitOperation

# A crossing that changes the vapor fraction by more than this gets a warning.
VF_WARNING = 0.02
# Unit categories whose result depends mostly on phase equilibrium; these get a suggestion.
SEPARATION_CATEGORIES = ("Separators", "Columns")
# Mole fraction below which a component is a trace and does not drive the method choice.
TRACE = 1e-3
# Above this pressure, dissolved supercritical gases make ideal-vapor activity models doubtful.
HIGH_PRESSURE = 10e5


def check_method(method: str, where: str, default: str) -> str:
    """Validate a unit's ``property_method`` override and return it."""
    if method in METHODS:
        return method
    norm = method.strip().lower().replace(" ", "-").replace("_", "-")
    close = norm if norm in METHODS else next(
        iter(difflib.get_close_matches(norm, list(METHODS), n=1, cutoff=0.6)), None)
    hint = f"did you mean {close!r}? " if close else ""
    if method == "auto":
        hint = "'auto' is only for the flowsheet default; "
    raise SimulationError(
        f"unknown property method {method!r}",
        where=where,
        hint=f"{hint}use one of: {', '.join(METHODS)}, or remove property_method to use "
        f"the flowsheet default ({default})",
    )


@dataclass
class _Crossing:
    stream: str
    src_unit: str
    dst_unit: str
    produced: Stream
    consumed: Stream


class PropertySections:
    """Which package each unit uses, and the re-flash of streams between them."""

    def __init__(self, model: FlowsheetModel, default: PropertyPackage,
                 units: dict[str, UnitOperation], connections: list):
        self.default = default
        self.units = units
        self.methods: dict[str, str] = {}
        self.sources: dict[str, str] = {}  # "unit", "flowsheet" or "inlet"
        components = tuple(model.components)
        for um in model.units:
            if um.property_method is not None:
                self.methods[um.id] = check_method(um.property_method, um.id, default.method)
                self.sources[um.id] = "unit"
            else:
                self.methods[um.id] = default.method
                self.sources[um.id] = "flowsheet"
        for c in connections:  # products report in the method of the unit feeding them
            if (c.dst_unit is not None and units[c.dst_unit].type == "Product"
                    and self.sources[c.dst_unit] == "flowsheet"):
                self.methods[c.dst_unit] = self.methods[c.src_unit]
                self.sources[c.dst_unit] = "inlet"
        self.packages: dict[str, PropertyPackage] = {default.method: default}
        for method in sorted(set(self.methods.values()) - {default.method}):
            self.packages[method] = get_package(components, method)
        self._crossings: dict[str, _Crossing] = {}

    @property
    def mixed(self) -> bool:
        return len(self.packages) > 1

    def warnings(self) -> list[str]:
        """Warnings of every package in use beyond the default (e.g. estimated parameters),
        with the units that use it; one line per kind of warning rather than per pair."""
        out = []
        for method, pkg in self.packages.items():
            if method == self.default.method:
                continue
            label = METHODS[method].label
            users = [u for u, m in self.methods.items()
                     if m == method and self.sources[u] != "inlet"]
            where = f"{label} in {', '.join(users)}"
            for w in _grouped(pkg.warnings):
                out.append(f"{where}: {w.removeprefix(label + ': ')}")
        return out

    def package(self, unit_id: str) -> PropertyPackage:
        return self.packages[self.methods[unit_id]]

    def inlet(self, conn, stream: Stream) -> Stream:
        """``stream`` as the unit at the end of ``conn`` sees it: re-flashed at the same T
        and P when it comes from a unit with a different method."""
        if conn.dst_unit is None or self.methods[conn.src_unit] == self.methods[conn.dst_unit]:
            return stream
        known = self._crossings.get(conn.id)
        if known is not None and known.produced is stream:
            return known.consumed
        pkg = self.package(conn.dst_unit)
        try:
            consumed = Stream.from_flash(pkg, stream.flows, T=stream.T, P=stream.P)
        except SimulationError as exc:
            raise SimulationError(
                f"re-flash with {pkg.method} failed where the stream enters {conn.dst_unit}: "
                f"{exc.message}", where=conn.id,
                hint=f"{conn.src_unit} uses {self.methods[conn.src_unit]} and {conn.dst_unit} "
                f"uses {pkg.method}; give both the same method, or move the boundary to a "
                "stream that is clearly one phase",
            ) from exc
        self._crossings[conn.id] = _Crossing(conn.id, conn.src_unit, conn.dst_unit,
                                             stream, consumed)
        return consumed

    def annotate(self, result, streams: dict[str, Stream], connections: list) -> None:
        """Add methods, suggestions and crossings to a solved ``SolveResult``."""
        src = {c.id: c.src_unit for c in connections}
        for sid, d in result.streams.items():  # a stream is reported in its producer's method
            d["property_method"] = self.methods.get(src.get(sid, sid.split(".")[0]),
                                                    self.default.method)
        for uid, d in result.units.items():
            inlets = [streams[c.id] for c in connections if c.dst_unit == uid and c.id in streams]
            d.update(self.unit_info(uid, inlets))
        result.property_crossings = self.crossings()
        result.warnings.extend(self.crossing_warnings())

    def crossings(self) -> list[dict]:
        """Every boundary stream with its state on both sides, in the converged solution."""
        out = []
        for x in self._crossings.values():
            a, b = x.produced, x.consumed
            out.append({
                "stream": x.stream,
                "from_unit": x.src_unit,
                "to_unit": x.dst_unit,
                "from_method": self.methods[x.src_unit],
                "to_method": self.methods[x.dst_unit],
                "T": a.T,
                "P": a.P,
                "vapor_fraction_before": a.state.VF,
                "vapor_fraction_after": b.state.VF,
                "phase_before": a.state.phase,
                "phase_after": b.state.phase,
                "enthalpy_change": (b.state.H - a.state.H) * a.total_flow,  # W
            })
        return out

    def crossing_warnings(self) -> list[str]:
        out = []
        for c in self.crossings():
            split_changed = (self._crossings[c["stream"]].produced.state.liquid_split
                             != self._crossings[c["stream"]].consumed.state.liquid_split)
            if abs(c["vapor_fraction_after"] - c["vapor_fraction_before"]) > VF_WARNING \
                    or split_changed:
                out.append(
                    f"{c['stream']}: switching from {c['from_method']} ({c['from_unit']}) to "
                    f"{c['to_method']} ({c['to_unit']}) changes the stream from "
                    f"{_phase(c, 'before')} to {_phase(c, 'after')} at the same T and P; "
                    "check which method describes this stream better"
                )
        return out

    def unit_info(self, unit_id: str, inlets: list[Stream]) -> dict:
        """What the solve result reports per unit: method, where it came from, suggestion."""
        info: dict = {"property_method": self.methods[unit_id],
                      "property_method_source": self.sources[unit_id]}
        unit = self.units[unit_id]
        if unit.category in SEPARATION_CATEGORIES and inlets:
            try:
                s = suggest_for_unit(self.package(unit_id), inlets, unit.p.get("P"))
            except Exception:  # noqa: BLE001 - a suggestion never fails the solve
                s = None
            if s is not None:
                info["property_method_suggestion"] = {"method": s.method, "reason": s.reason}
        return info


_PAIR = re.compile(r" for (.+?/.+?)(?=[,;] )")


def _grouped(warnings: list[str]) -> list[str]:
    """Merge warnings that differ only in the binary pair they name."""
    groups: dict[str, list[str | None]] = {}
    for w in warnings:
        m = _PAIR.search(w)
        key = w if m is None else w[:m.start(1)] + "{pair}" + w[m.end(1):]
        groups.setdefault(key, []).append(m.group(1) if m else None)
    out = []
    for key, pairs in groups.items():
        if len(pairs) == 1:
            out.append(key.replace("{pair}", pairs[0] or "{pair}"))
        else:
            out.append(key.replace("{pair}", f"{len(pairs)} pairs ({', '.join(pairs)})")
                       .replace("this pair is treated as an ideal solution",
                                "these pairs are treated as ideal solutions")
                       .replace("this pair is", "these pairs are")
                       .replace("for this pair", "for these pairs")
                       .replace("data for it", "data for them"))
    return out


def _phase(c: dict, side: str) -> str:
    vf = c[f"vapor_fraction_{side}"]
    return f"{c[f'phase_{side}']} (vapor fraction {vf:.3f})"


def suggest_for_unit(pkg: PropertyPackage, inlets: list[Stream],
                     P_unit: float | None = None) -> Suggestion | None:
    """Suggest a method for one unit from what it actually sees: the components present in
    its feed (traces ignored) and its pressure."""
    flows = [sum(s.flows[i] for s in inlets) for i in range(pkg.n)]
    total = sum(flows)
    if total <= 0.0:
        return None
    zs = [f / total for f in flows]
    present = [pkg.names[i] for i in range(pkg.n) if zs[i] >= TRACE]
    P = P_unit or max(s.P for s in inlets)
    T = max(s.T for s in inlets)
    base = suggest_method(present)

    supercritical = [pkg.names[i] for i in range(pkg.n)
                     if zs[i] >= 0.01 and (pkg.consts.Tcs[i] or 0.0) < T]
    if P > HIGH_PRESSURE and supercritical and base.method != "peng-robinson":
        return Suggestion(
            "peng-robinson",
            f"At {P / 1e5:.0f} bar the feed carries dissolved gases above their critical "
            f"temperature ({', '.join(supercritical)}). Activity-coefficient models assume an "
            "ideal vapor and cannot describe supercritical components, so an equation of state "
            "(Peng-Robinson) is the safer choice here, even with polar liquids present.",
        )
    if base.method == "unifac":
        # The liquid's non-ideality comes from the condensable compounds; gases above their
        # critical temperature are nearly ideal in the vapor at moderate pressure.
        idx = [i for i in range(pkg.n) if zs[i] >= 0.01 and pkg.names[i] not in supercritical]
        cas = [pkg.compounds[i].cas for i in idx]
        pairs = [(a, b) for k, a in enumerate(cas) for b in cas[k + 1:]]
        if pairs and all(IPDB.has_ip_specific("ChemSep NRTL", [a, b], "bij") for a, b in pairs):
            names = [pkg.names[i] for i in idx]
            main = ", ".join(names[:-1]) + f" and {names[-1]}"
            gases = (f" The dissolved gases ({', '.join(supercritical)}) matter little for the "
                     "liquid at this pressure; their interactions are approximated." if supercritical
                     else " Traces are estimated.")
            return Suggestion(
                "nrtl",
                f"The liquid is mainly {main}, a non-ideal mixture at {P / 1e5:.1f} bar. ChemSep "
                "has NRTL parameters regressed to measured data for every pair of these, which "
                "is more accurate than the predictive UNIFAC model." + gases,
            )
    return base
