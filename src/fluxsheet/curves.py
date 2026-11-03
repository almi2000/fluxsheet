"""Performance curves of pumps and compressors: head and efficiency against flow.

A curve is given as a table of points at the machine's rated speed::

    {"flow_unit": "m3/h", "head_unit": "m", "efficiency_unit": "",
     "points": [[0, 52], [20, 50, 0.55], [40, 45, 0.72], [60, 36, 0.70]]}

Rows are ``[flow, head]`` or ``[flow, head, efficiency]`` (every row alike); flows are actual
volumetric flows at the inlet. The head may be a length (m, ft: a liquid column) or a specific
energy (J/kg, kJ/kg); it is stored as specific energy either way. A zero efficiency is allowed
at zero flow only, where it is ignored.

Between the points the curve is a monotone piecewise cubic (PCHIP: smooth, no overshoot);
two points give a straight line. Outside the points it continues as a straight line with the
end slope, and the head never drops below zero.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.interpolate import PchipInterpolator

from .quantities import UnitError, known_units, to_si

G = 9.80665  # m/s2

DEFAULT_UNITS = {"flow_unit": "m3/h", "head_unit": "m", "efficiency_unit": ""}


@dataclass(frozen=True)
class Curve:
    flows: tuple[float, ...]  # m3/s, strictly increasing
    heads: tuple[float, ...]  # J/kg
    eff_flows: tuple[float, ...] | None  # m3/s; None: the curve has no efficiency column
    effs: tuple[float, ...] | None  # fractions

    def inside(self, q: float) -> bool:
        return self.flows[0] <= q <= self.flows[-1]

    def head(self, q: float) -> float:
        """Head in J/kg at actual volumetric flow ``q`` (m3/s), at the rated speed."""
        return max(_evaluate(self.flows, self.heads, q), 0.0)

    def efficiency(self, q: float) -> float | None:
        if self.effs is None:
            return None
        return min(max(_evaluate(self.eff_flows, self.effs, q), 0.01), 1.0)


@lru_cache(maxsize=256)
def _interpolator(xs: tuple[float, ...], ys: tuple[float, ...]):
    if len(xs) == 2:
        slope = (ys[1] - ys[0]) / (xs[1] - xs[0])
        return None, slope, slope
    f = PchipInterpolator(np.array(xs), np.array(ys), extrapolate=False)
    d = f.derivative()
    return f, float(d(xs[0])), float(d(xs[-1]))


def _evaluate(xs: tuple[float, ...], ys: tuple[float, ...], x: float) -> float:
    f, slope_lo, slope_hi = _interpolator(xs, ys)
    if x < xs[0]:
        return ys[0] + slope_lo * (x - xs[0])
    if x > xs[-1]:
        return ys[-1] + slope_hi * (x - xs[-1])
    if f is None:
        return ys[0] + slope_lo * (x - xs[0])
    return float(f(x))


def parse_curve(value) -> Curve:
    """Read and check a curve table; raises UnitError with what is wrong."""
    if not isinstance(value, dict) or "points" not in value:
        raise UnitError("expected a table with 'points' (rows of flow, head[, efficiency])")
    unknown = set(value) - {"points", *DEFAULT_UNITS}
    if unknown:
        raise UnitError(f"unknown curve entries: {', '.join(sorted(unknown))}")
    units = {k: value.get(k) or v for k, v in DEFAULT_UNITS.items()}
    for key in DEFAULT_UNITS:
        if not isinstance(units[key], str):
            raise UnitError(f"{key} must be a unit name, got {units[key]!r}")

    def convert(number, quantity: str, unit: str, what: str) -> float:
        if isinstance(number, bool) or not isinstance(number, (int, float)):
            raise UnitError(f"{what} {number!r} is not a number")
        try:
            return to_si(f"{float(number)!r} {unit}", quantity)
        except UnitError:
            raise UnitError(f"unknown {what} unit {unit!r}; use one of: "
                            f"{', '.join(u for u in known_units(quantity) if u)}") from None

    head_unit = units["head_unit"]
    if head_unit.lower().replace(" ", "") in known_units("length"):
        def head_si(h):
            return convert(h, "length", head_unit, "head") * G
    elif head_unit.lower().replace(" ", "") in known_units("specific_energy"):
        def head_si(h):
            return convert(h, "specific_energy", head_unit, "head")
    else:
        raise UnitError(f"unknown head unit {head_unit!r}; use a length "
                        f"({', '.join(known_units('length'))}) or an energy per mass "
                        f"({', '.join(known_units('specific_energy'))})")

    rows = value["points"]
    if not isinstance(rows, list) or len(rows) < 2:
        n = len(rows) if isinstance(rows, list) else 0
        raise UnitError(f"a curve needs at least 2 points, got {n}")
    width = len(rows[0]) if isinstance(rows[0], (list, tuple)) else 0
    if width not in (2, 3):
        raise UnitError("each point is [flow, head] or [flow, head, efficiency]")
    flows, heads, effs = [], [], []
    for k, row in enumerate(rows, start=1):
        if not isinstance(row, (list, tuple)) or len(row) != width:
            raise UnitError(f"point {k} has {len(row) if isinstance(row, (list, tuple)) else 0}"
                            f" values; every point needs {width}")
        q = convert(row[0], "volumetric_flow", units["flow_unit"], "flow")
        h = head_si(row[1])
        if q < 0.0:
            raise UnitError(f"point {k}: the flow is negative")
        if h < 0.0:
            raise UnitError(f"point {k}: the head is negative")
        if flows and q <= flows[-1]:
            raise UnitError(f"the flows must increase from point to point; point {k} does not "
                            f"({row[0]} after {rows[k - 2][0]})")
        flows.append(q)
        heads.append(h)
        if width == 3:
            eta = convert(row[2], "fraction", units["efficiency_unit"], "efficiency")
            if eta == 0.0 and q == 0.0:
                continue  # shut-off: no efficiency to speak of
            if not 0.0 < eta <= 1.0:
                raise UnitError(f"point {k}: efficiency {row[2]} is not between 0 and 1"
                                + ("" if units["efficiency_unit"] else
                                   "; for percent set efficiency_unit to '%'"))
            effs.append((q, eta))
    if width == 3 and len(effs) < 2:
        raise UnitError("the efficiency column needs at least 2 points with an efficiency")
    return Curve(
        flows=tuple(flows), heads=tuple(heads),
        eff_flows=tuple(q for q, _ in effs) if effs else None,
        effs=tuple(e for _, e in effs) if effs else None,
    )
