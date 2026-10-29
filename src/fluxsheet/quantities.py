"""Engineering units in, SI out.

Flowsheet files may give any physical value either as a plain number (SI) or as a
string with a unit, e.g. ``"80 degC"``, ``"1.5 bar"`` or ``"10 kmol/h"``. Internally
everything is SI: K, Pa, mol/s, W, J/mol, m3, m, s.
"""

from __future__ import annotations

import re

# quantity -> unit -> (factor, offset): SI = value * factor + offset
_UNITS: dict[str, dict[str, tuple[float, float]]] = {
    "temperature": {
        "k": (1.0, 0.0),
        "c": (1.0, 273.15),
        "degc": (1.0, 273.15),
        "°c": (1.0, 273.15),
        "f": (5.0 / 9.0, 273.15 - 32.0 * 5.0 / 9.0),
        "degf": (5.0 / 9.0, 273.15 - 32.0 * 5.0 / 9.0),
        "°f": (5.0 / 9.0, 273.15 - 32.0 * 5.0 / 9.0),
    },
    "pressure": {
        "pa": (1.0, 0.0),
        "kpa": (1e3, 0.0),
        "mpa": (1e6, 0.0),
        "bar": (1e5, 0.0),
        "bara": (1e5, 0.0),
        "mbar": (1e2, 0.0),
        "atm": (101325.0, 0.0),
        "psi": (6894.757293168, 0.0),
        "psia": (6894.757293168, 0.0),
    },
    "pressure_difference": {
        "pa": (1.0, 0.0),
        "kpa": (1e3, 0.0),
        "mpa": (1e6, 0.0),
        "bar": (1e5, 0.0),
        "mbar": (1e2, 0.0),
        "psi": (6894.757293168, 0.0),
    },
    "molar_flow": {
        "mol/s": (1.0, 0.0),
        "kmol/s": (1e3, 0.0),
        "mol/h": (1.0 / 3600.0, 0.0),
        "kmol/h": (1e3 / 3600.0, 0.0),
        "mol/min": (1.0 / 60.0, 0.0),
    },
    "power": {
        "w": (1.0, 0.0),
        "kw": (1e3, 0.0),
        "mw": (1e6, 0.0),
    },
    "volume": {
        "m3": (1.0, 0.0),
        "l": (1e-3, 0.0),
        "cm3": (1e-6, 0.0),
    },
    "length": {
        "m": (1.0, 0.0),
        "cm": (1e-2, 0.0),
        "mm": (1e-3, 0.0),
    },
    "area": {
        "m2": (1.0, 0.0),
        "cm2": (1e-4, 0.0),
    },
    "time": {
        "s": (1.0, 0.0),
        "min": (60.0, 0.0),
        "h": (3600.0, 0.0),
    },
    "molar_energy": {
        "j/mol": (1.0, 0.0),
        "kj/mol": (1e3, 0.0),
        "kcal/mol": (4184.0, 0.0),
    },
    "fraction": {
        "": (1.0, 0.0),
        "%": (0.01, 0.0),
    },
}

_PATTERN = re.compile(r"^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)\s*(.*?)\s*$")


class UnitError(ValueError):
    """Raised when a value cannot be converted to SI."""


def known_units(quantity: str) -> list[str]:
    return list(_UNITS[quantity])


def to_si(value: float | str, quantity: str) -> float:
    """Convert ``value`` (number in SI, or ``"<number> <unit>"``) to SI."""
    if isinstance(value, bool):
        raise UnitError(f"expected a {quantity.replace('_', ' ')}, got {value!r}")
    if isinstance(value, (int, float)):
        return float(value)
    if not isinstance(value, str):
        raise UnitError(f"expected a {quantity.replace('_', ' ')}, got {value!r}")
    match = _PATTERN.match(value)
    if not match:
        raise UnitError(f"could not read {value!r} as a {quantity.replace('_', ' ')}")
    number, unit = float(match.group(1)), match.group(2).lower().replace(" ", "")
    table = _UNITS[quantity]
    if unit == "" and "" not in table:
        return number
    if unit not in table:
        raise UnitError(
            f"unknown unit {match.group(2)!r} for {quantity.replace('_', ' ')}; "
            f"use one of: {', '.join(u for u in table if u)}"
        )
    factor, offset = table[unit]
    return number * factor + offset
