"""Unit operations, discovered automatically.

Every module in this package is imported, and every concrete ``UnitOperation`` subclass
that defines ``type`` is registered. To add a family of units, add a module here; nothing
else needs to change.
"""

import importlib
import inspect
import pkgutil

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
                  "Reactors", "Pressure Changers"]


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
