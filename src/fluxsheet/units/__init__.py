"""Unit operations, discovered automatically.

Every module in this package is imported, and every concrete ``UnitOperation`` subclass
that defines ``type`` is registered. To add a family of units, add a module here; nothing
else needs to change. Other packages add units through the ``fluxsheet.units`` entry point
(see ``fluxsheet.plugins``); they come after Fluxsheet's own and cannot replace them.
"""

import importlib
import inspect
import pkgutil
import sys

from .. import plugins
from .base import Param, Port, ResultField, UnitOperation

for _mod in pkgutil.iter_modules(__path__):
    if _mod.name != "base":
        importlib.import_module(f"{__name__}.{_mod.name}")


def _concrete(cls: type[UnitOperation]):
    for sub in cls.__subclasses__():
        if "type" in vars(sub):
            yield sub
        yield from _concrete(sub)


# Palette order in the editor; categories not listed here come last, alphabetically.
CATEGORY_ORDER = ["Streams", "Mixers & Splitters", "Heat Exchangers", "Separators", "Columns",
                  "Reactors", "Pressure Changers", "Manipulators"]


def _sort_key(cls: type[UnitOperation]):
    cat = cls.category
    rank = CATEGORY_ORDER.index(cat) if cat in CATEGORY_ORDER else len(CATEGORY_ORDER)
    return rank, cat, cls.__module__, inspect.getsourcelines(cls)[1]


UNIT_TYPES: dict[str, type[UnitOperation]] = {}
for _cls in sorted(_concrete(UnitOperation), key=_sort_key):
    if _cls.type in UNIT_TYPES and UNIT_TYPES[_cls.type] is not _cls:
        raise RuntimeError(f"two unit operations share the type {_cls.type!r}")
    UNIT_TYPES[_cls.type] = _cls

__all__ = ["UNIT_TYPES", "Param", "Port", "ResultField", "UnitOperation"]


def _register_plugin(plugin, value) -> list[str]:
    """Register the unit operations an entry point brings: a module (its classes and those
    of the modules it imports from its own package), a class, or a list of classes."""
    if inspect.ismodule(value):
        root = value.__name__.split(".")[0]
        # Classes still defined in a loaded module of the package (not a replaced one's).
        new = [c for c in _concrete(UnitOperation) if c.__module__.split(".")[0] == root
               and getattr(sys.modules.get(c.__module__), c.__name__, None) is c]
    else:
        new = [c for v in plugins.objects(value) for c in (v, *_concrete(v))
               if isinstance(c, type) and "type" in vars(c)]
        if not new or not all(issubclass(c, UnitOperation) for c in new):
            raise ValueError("expected a module, UnitOperation subclasses, or a list of them")
    new = [c for c in dict.fromkeys(new) if UNIT_TYPES.get(c.type) is not c]
    if not new:
        raise ValueError("defines no unit operation")
    taken = [c.type for c in new if c.type in UNIT_TYPES]
    if taken or len({c.type for c in new}) < len(new):
        raise ValueError(f"unit type already defined: {', '.join(sorted(set(taken))) or 'twice'}")
    for c in new:
        for attr in ("label", "category", "description"):
            if not isinstance(getattr(c, attr, None), str):
                raise TypeError(f"{c.__name__} has no {attr}")
    for c in new:
        c.plugin = plugin.distribution or plugin.name
        UNIT_TYPES[c.type] = c
    return [c.type for c in new]


plugins.load("units", _register_plugin)
