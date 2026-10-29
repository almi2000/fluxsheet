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


@dataclass(frozen=True)
class Param:
    name: str
    label: str
    quantity: str | None = None  # see quantities.py; "composition" = {component: molar flow}
    default: Any = None
    help: str = ""
    min: float | None = None
    max: float | None = None

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class Port:
    name: str
    label: str
    multiple: bool = False  # e.g. a mixer accepts any number of inlet streams


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

    def __post_init__(self) -> None:
        self.p: dict[str, Any] = {}
        self.results: dict[str, Any] = {}
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
        if param.quantity in ("composition", "component_fractions"):
            if not isinstance(value, dict):
                raise UnitError(f"expected a mapping of component to value, got {value!r}")
            quantity = "molar_flow" if param.quantity == "composition" else "fraction"
            return {str(k): to_si(v, quantity) for k, v in value.items()}
        if param.quantity is None:
            return float(value)
        return to_si(value, param.quantity)

    def _check_bounds(self, param: Param) -> None:
        value = self.p[param.name]
        if not isinstance(value, float):
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
        }


def empty_stream(pkg: PropertyPackage, T: float, P: float) -> Stream:
    """A zero-flow stream; carries T and P so downstream units have something to work with."""
    return Stream.from_flash(pkg, [0.0] * pkg.n, T=T, P=P)
