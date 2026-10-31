"""Base class and parameter schema for unit operations.

Each unit declares its ports and parameters as data. The same declaration drives
input validation, unit conversion, and the forms in the web editor.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, ClassVar

from ..errors import SimulationError
from ..quantities import UnitError, to_si
from ..stream import Stream
from ..thermo_package import PropertyPackage

# Parameter kinds besides the physical quantities in quantities.py:
#   None                  plain number
#   "integer"             whole number (stage counts, stage positions)
#   "choice"              one of `options`
#   "component"           name of one flowsheet component
#   "string_list"         list of names, e.g. reaction ids
#   "composition"         {component: molar flow}
#   "component_fractions" {component: fraction}
#   "stream_stages"       {inlet stream id: stage number}, e.g. column feed stages
SPECIAL_KINDS = ("integer", "choice", "component", "string_list", "composition",
                 "component_fractions", "stream_stages")


@dataclass(frozen=True)
class Param:
    name: str
    label: str
    quantity: str | None = None
    default: Any = None
    help: str = ""
    min: float | None = None
    max: float | None = None
    options: tuple[str, ...] | None = None  # for quantity="choice"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["options"] = list(self.options) if self.options else None
        return d


@dataclass(frozen=True)
class ResultField:
    """A scalar a unit reports after solving, e.g. its duty. Drives the editor display."""

    name: str
    label: str
    quantity: str | None = None  # "power", "temperature", ... or None for dimensionless


@dataclass(frozen=True)
class Port:
    name: str
    label: str
    # Accepts any number of streams. Unless `separate`, the solver mixes them adiabatically
    # into one before the unit sees them; a separate port gets every stream on its own.
    multiple: bool = False
    separate: bool = False


@dataclass
class UnitOperation:
    id: str
    raw: dict[str, Any] = field(default_factory=dict)

    type: ClassVar[str]
    label: ClassVar[str]
    category: ClassVar[str]
    description: ClassVar[str]
    inlets: ClassVar[tuple[Port, ...]] = ()
    outlets: ClassVar[tuple[Port, ...]] = ()
    params: ClassVar[tuple[Param, ...]] = ()
    # Each group lists alternative specifications; at most one per group may be set.
    # Groups marked as required in `required_groups` need exactly one.
    exclusive_groups: ClassVar[tuple[tuple[str, ...], ...]] = ()
    required_groups: ClassVar[tuple[tuple[str, ...], ...]] = ()
    result_fields: ClassVar[tuple[ResultField, ...]] = ()

    def __post_init__(self) -> None:
        self.p: dict[str, Any] = {}
        self.results: dict[str, Any] = {}
        # Tabular results along stages or reactor length; see add_profile().
        self.profiles: dict[str, dict] = {}
        # Ids of the streams on each inlet port, in the order calculate() receives them.
        self.inlet_ids: dict[str, list[str]] = {}
        known = {param.name for param in self.params}
        for name in self.raw:
            if name not in known:
                raise SimulationError(
                    f"unknown parameter {name!r}",
                    where=self.id,
                    hint=f"{self.label} accepts: {', '.join(sorted(known)) or 'no parameters'}",
                )
        for param in self.params:
            value = self.raw.get(param.name, param.default)
            if value is None:
                self.p[param.name] = None
                continue
            try:
                self.p[param.name] = self._convert(param, value)
            except UnitError as exc:
                raise SimulationError(f"{param.label}: {exc}", where=self.id) from exc
            self._check_bounds(param)
        self._check_groups()

    def _convert(self, param: Param, value: Any) -> Any:
        q = param.quantity
        if q == "integer":
            if isinstance(value, bool) or not float(value).is_integer():
                raise UnitError(f"expected a whole number, got {value!r}")
            return int(value)
        if q == "choice":
            if value not in (param.options or ()):
                raise UnitError(f"{value!r} is not an option; use one of: "
                                f"{', '.join(param.options or ())}")
            return value
        if q == "component":
            if not isinstance(value, str) or not value:
                raise UnitError(f"expected a component name, got {value!r}")
            return value
        if q == "string_list":
            if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
                raise UnitError(f"expected a list of names, got {value!r}")
            return list(value)
        if q == "stream_stages":
            if not isinstance(value, dict):
                raise UnitError(f"expected a mapping of stream to stage, got {value!r}")
            stages = {}
            for sid, stage in value.items():
                if isinstance(stage, bool) or not float(stage).is_integer():
                    raise UnitError(f"stage of {sid} must be a whole number, got {stage!r}")
                stages[str(sid)] = int(stage)
            return stages
        if q in ("composition", "component_fractions"):
            if not isinstance(value, dict):
                raise UnitError(f"expected a mapping of component to value, got {value!r}")
            quantity = "molar_flow" if param.quantity == "composition" else "fraction"
            return {str(k): to_si(v, quantity) for k, v in value.items()}
        if param.quantity is None:
            return float(value)
        return to_si(value, param.quantity)

    def _check_bounds(self, param: Param) -> None:
        value = self.p[param.name]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return
        if (param.min is not None and value < param.min) or (
            param.max is not None and value > param.max
        ):
            bounds = f"between {param.min} and {param.max}" if param.max is not None else (
                f"at least {param.min}"
            )
            raise SimulationError(f"{param.label} must be {bounds}, got {value:g}", where=self.id)

    def _check_groups(self) -> None:
        for group in self.exclusive_groups:
            given = [name for name in group if self.p.get(name) is not None]
            labels = " or ".join(self._label(n) for n in group)
            if len(given) > 1:
                raise SimulationError(
                    f"over-specified: {' and '.join(self._label(n) for n in given)} are both set",
                    where=self.id,
                    hint=f"set only one of {labels}",
                )
            if not given and group in self.required_groups:
                raise SimulationError(
                    f"under-specified: {labels} is missing", where=self.id,
                    hint=f"set one of {labels}",
                )

    @classmethod
    def _label(cls, name: str) -> str:
        return next(p.label for p in cls.params if p.name == name)

    def component_index(self, pkg: PropertyPackage, name: str, label: str = "component") -> int:
        """Index of a component parameter in the package, with a helpful error if absent."""
        if name not in pkg.names:
            raise SimulationError(
                f"{label} {name!r} is not in the flowsheet's component list", where=self.id,
                hint=f"choose one of: {', '.join(pkg.names)}",
            )
        return pkg.names.index(name)

    def add_profile(
        self,
        name: str,
        label: str,
        x_label: str,
        x: list[float],
        series: dict[str, tuple[str | None, list[float]]],
    ) -> None:
        """Record a profile, e.g. temperature per stage. ``series``: label -> (quantity, values)."""
        self.profiles[name] = {
            "label": label,
            "x_label": x_label,
            "x": list(x),
            "series": [
                {"label": k, "quantity": q, "values": list(v)} for k, (q, v) in series.items()
            ],
        }

    def outlet_pressure(self, p_in: float) -> float:
        """Outlet pressure from an optional absolute ``P`` or pressure drop ``dP`` parameter."""
        if self.p.get("P") is not None:
            return self.p["P"]
        return p_in - (self.p.get("dP") or 0.0)

    def calculate(
        self, pkg: PropertyPackage, inlets: dict[str, list[Stream]]
    ) -> dict[str, Stream]:
        raise NotImplementedError

    @classmethod
    def schema(cls) -> dict:
        return {
            "type": cls.type,
            "label": cls.label,
            "category": cls.category,
            "description": cls.description,
            "inlets": [asdict(p) for p in cls.inlets],
            "outlets": [asdict(p) for p in cls.outlets],
            "params": [p.to_dict() for p in cls.params],
            "exclusive_groups": [list(g) for g in cls.exclusive_groups],
            "required_groups": [list(g) for g in cls.required_groups],
            "result_fields": [asdict(f) for f in cls.result_fields],
        }


def empty_stream(pkg: PropertyPackage, T: float, P: float) -> Stream:
    """A zero-flow stream; carries T and P so downstream units have something to work with."""
    return Stream.from_flash(pkg, [0.0] * pkg.n, T=T, P=P)


def single_inlet(inlets: dict[str, list[Stream]], port: str = "in") -> Stream:
    return inlets[port][0]


def mix_streams(pkg: PropertyPackage, streams: list[Stream], P: float | None = None) -> Stream:
    """Adiabatic mixing; the pressure defaults to the lowest inlet pressure."""
    P = P if P is not None else min(s.P for s in streams)
    flows = [sum(s.flows[i] for s in streams) for i in range(pkg.n)]
    total = sum(flows)
    if total <= 0.0:
        return empty_stream(pkg, streams[0].T, P)
    return Stream.from_flash(pkg, flows, P=P, H=sum(s.enthalpy_flow for s in streams) / total)


# Building blocks shared by many units. Most units take any number of feeds on their inlet,
# mixed before they enter, as a Mixer in front would.
IN = Port("in", "Inlet", multiple=True)
OUT = Port("out", "Outlet")
P_OUT = Param("P", "Outlet pressure", "pressure", help="Absolute outlet pressure.", min=1.0)
DP = Param("dP", "Pressure drop", "pressure_difference",
           help="Pressure loss across the unit. Defaults to zero.", min=0.0)
DUTY = ResultField("duty", "Duty", "power")
POWER = ResultField("power", "Power", "power")
