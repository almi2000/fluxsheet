"""Plug-ins: unit operations, property methods and compounds from other installed packages.

A package registers them in its ``pyproject.toml``::

    [project.entry-points."fluxsheet.units"]
    my_units = "my_package.units"              # a module of UnitOperation subclasses

    [project.entry-points."fluxsheet.property_methods"]
    regular = "my_package.thermo:REGULAR"      # an ActivityMethod or EquationOfState

    [project.entry-points."fluxsheet.components"]
    oils = "my_package.compounds:COMPOUNDS"    # custom-compound definitions

An entry point may name a module, an object, a list of objects, or a function without
arguments that returns either. Plug-ins are loaded once, when Fluxsheet's unit and property
registries are first imported. A plug-in that fails to load, or that reuses a name Fluxsheet
or another plug-in already has, is skipped and reported (``fluxsheet plugins``,
``GET /api/plugins``); the rest of Fluxsheet keeps working. Set the environment variable
``FLUXSHEET_DISABLE_PLUGINS=1`` to load none.

Plug-ins are ordinary installed Python packages and run with Fluxsheet's permissions; only
install ones you trust. See ``examples/plugin/`` for a complete one.
"""

from __future__ import annotations

import inspect
import os
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from importlib import metadata

GROUPS = {
    "units": "fluxsheet.units",
    "property_methods": "fluxsheet.property_methods",
    "components": "fluxsheet.components",
}
DISABLE_ENV = "FLUXSHEET_DISABLE_PLUGINS"


@dataclass(frozen=True)
class ActivityMethod:
    """A property method from an activity-coefficient model: Fluxsheet supplies the ideal
    gas, vapor pressures, liquid volumes, heats of vaporization and the flash; the plug-in
    supplies the excess Gibbs energy.

    ``gibbs_excess(consts, components)`` returns a ``thermo`` ``GibbsExcess`` model for the
    components with indices ``components`` (all of them, in order). ``consts`` is the
    ``thermo.ChemicalConstantsPackage`` of the flowsheet's compounds. ``liquid_liquid``:
    whether flashes should look for two liquid phases.
    """

    key: str
    label: str
    description: str
    gibbs_excess: Callable
    liquid_liquid: bool = False


@dataclass(frozen=True)
class EquationOfState:
    """A property method from a cubic equation of state for both phases.

    ``eos_mix`` is a ``thermo`` mixture class (a ``GCEOSMIX`` subclass such as ``PRMIX``),
    built with ``Tcs``, ``Pcs``, ``omegas`` and ``kijs`` plus ``extra(consts)``, a function
    returning further keyword arguments. ``kijs(consts)`` returns the binary interaction
    matrix; zeros if omitted.
    """

    key: str
    label: str
    description: str
    eos_mix: type
    kijs: Callable | None = None
    extra: Callable | None = None


@dataclass
class Plugin:
    """One loaded (or failed) entry point."""

    group: str
    name: str
    value: str
    distribution: str | None
    version: str | None
    provides: list[str] = field(default_factory=list)
    error: str | None = None

    def to_dict(self) -> dict:
        return {"group": self.group, "name": self.name, "value": self.value,
                "distribution": self.distribution, "version": self.version,
                "provides": list(self.provides), "error": self.error}


PLUGINS: list[Plugin] = []


def disabled() -> bool:
    return os.environ.get(DISABLE_ENV, "").strip().lower() not in ("", "0", "false", "no")


def entry_points(group: str) -> list[metadata.EntryPoint]:
    """The installed entry points of one of the ``GROUPS``, sorted for a stable order."""
    if disabled():
        return []
    return sorted(metadata.entry_points(group=GROUPS[group]),
                  key=lambda ep: (ep.name, ep.value))


def _record(group: str, ep: metadata.EntryPoint) -> Plugin:
    dist = getattr(ep, "dist", None)
    return Plugin(group, ep.name, ep.value, dist.name if dist else None,
                  dist.version if dist else None)


def objects(value) -> list:
    """What an entry point provides: a callable without arguments is called first, and a
    list or tuple gives several objects."""
    if inspect.isroutine(value):
        value = value()
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def load(group: str, register: Callable[[Plugin, object], Iterable[str]],
         eps: Iterable[metadata.EntryPoint] | None = None) -> list[Plugin]:
    """Load each entry point of ``group`` and pass what it gives to ``register``, which
    returns the names it added or raises ``ValueError`` to reject it. Failures are recorded,
    never raised."""
    loaded = []
    for ep in entry_points(group) if eps is None else eps:
        plugin = _record(group, ep)
        try:
            plugin.provides = list(register(plugin, ep.load()))
        except Exception as exc:  # noqa: BLE001 - a broken plug-in must not stop Fluxsheet
            plugin.error = f"{type(exc).__name__}: {exc}"
        PLUGINS.append(plugin)
        loaded.append(plugin)
    return loaded


def report() -> list[dict]:
    return [p.to_dict() for p in PLUGINS]
