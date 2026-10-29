"""Phase diagrams: T-xy and P-xy of a binary, and phase envelopes of mixtures.

Everything is computed with the same property package as the flowsheet, so a diagram
shows exactly what the simulation will do with these compounds and this method.
"""

from __future__ import annotations

import math
from itertools import pairwise

from .errors import SimulationError
from .thermo_package import METHODS, get_package

EOS_METHODS = ("peng-robinson", "srk")


def binary_diagram(components: list[str], method: str, *, P: float | None = None,
                   T: float | None = None, points: int = 41) -> dict:
    """T-xy at fixed ``P`` or P-xy at fixed ``T`` for a binary; x and y of the first component.

    Returns the bubble curve (liquid x vs T or P), the dew curve (vapor y vs T or P), an
    azeotrope if the curves touch inside the diagram, and, for methods that predict a
    liquid-liquid split, the compositions of the two liquids.
    """
    if (P is None) == (T is None):
        raise SimulationError("give either a pressure (T-xy) or a temperature (P-xy)")
    pkg = _package(components, method, n=2)
    kind, fixed = ("txy", {"P": P}) if P is not None else ("pxy", {"T": T})
    key = "T" if kind == "txy" else "P"
    xs = [i / (points - 1) for i in range(points)]
    bubble, guess, gap = [], None, []
    for x in xs:
        xc = min(max(x, 1e-9), 1.0 - 1e-9)  # some models divide by mole fractions
        sat = pkg.saturation_point([xc, 1.0 - xc], "bubble", guess=guess, **fixed)
        value = getattr(sat, key)
        bubble.append((x, sat.y[0], value))
        # Inside a liquid-liquid gap every bubble point is the same three-phase point, so the
        # previous answer is the best guess; elsewhere the homogeneous estimate is better.
        het = pkg.supports_lle and 0.0 < x < 1.0 and _splits(pkg, x, sat, kind)
        guess = value if het else None
        if het:
            gap.append((x, value))

    return {
        "kind": kind,
        "components": list(components),
        "method": pkg.method,
        "P": P,
        "T": T,
        "x": [b[0] for b in bubble],
        "y": [b[1] for b in bubble],
        key: [b[2] for b in bubble],
        "azeotrope": _azeotrope(pkg, bubble, key, fixed),
        "liquid_liquid": _liquid_gap(pkg, gap, kind, fixed) if gap else None,
        "warnings": list(pkg.warnings),
        "parameter_sources": list(pkg.parameter_sources),
    }


def _package(components: list[str], method: str, n: int | None = None):
    if n is not None and len(components) != n:
        raise SimulationError(
            f"a binary diagram needs exactly two components, got {len(components)}",
            hint="pick two of the flowsheet's components",
        )
    if len(set(components)) != len(components):
        raise SimulationError("the same component is listed twice", hint="pick two different ones")
    if method not in METHODS:
        raise SimulationError(f"unknown property method {method!r}",
                              hint=f"choose one of: {', '.join(METHODS)}")
    return get_package(tuple(components), method)


def _splits(pkg, x: float, sat, kind: str) -> bool:
    """Whether the liquid at this bubble point has split, judged just below the boiling point."""
    spec = {"T": sat.T - 0.05, "P": sat.P} if kind == "txy" else {"T": sat.T, "P": sat.P * 1.002}
    return pkg.flash([x, 1.0 - x], **spec).liquid_split


def _azeotrope(pkg, bubble, key: str, fixed: dict) -> dict | None:
    """Where y - x changes sign between grid points; refined by the secant method."""
    d = [(x, y - x) for x, y, _ in bubble[1:-1]]
    for (x0, d0), (x1, d1) in pairwise(d):
        if d0 == 0.0 or d0 * d1 >= 0.0:
            continue
        for _ in range(4):
            x = x0 - d0 * (x1 - x0) / (d1 - d0)
            sat = pkg.saturation_point([x, 1.0 - x], "bubble", **fixed)
            dx = sat.y[0] - x
            x0, d0, x1, d1 = x1, d1, x, dx
            if abs(dx) < 1e-6:
                break
        split = pkg.supports_lle and pkg.flash(
            [x, 1.0 - x], **({"T": sat.T - 0.05, "P": sat.P} if key == "T"
                             else {"T": sat.T, "P": sat.P * 1.002})).liquid_split
        return {"x": x, key: getattr(sat, key),
                "type": "heterogeneous" if split else "homogeneous"}
    return None


def _liquid_gap(pkg, gap: list[tuple[float, float]], kind: str, fixed: dict) -> dict | None:
    """The two coexisting liquids (mole fraction of the first component) below the boiling line.

    ``gap`` holds the grid points whose boiling liquid splits, with their three-phase T or P.
    For T-xy the result is a solubility curve from the three-phase temperature downwards; for
    P-xy, the two compositions at the diagram's temperature (nearly pressure-independent).
    """
    z = 0.5 * (gap[0][0] + gap[-1][0])  # a feed between the two liquid compositions
    three_phase = sorted(v for _, v in gap)[len(gap) // 2]
    if kind == "pxy":
        st = pkg.flash([z, 1.0 - z], T=fixed["T"], P=three_phase * 1.01)
        if not st.liquid_split:
            return None
        return {"three_phase": three_phase, "x_light": st.x_light[0], "x_heavy": st.x_heavy[0]}
    curve = []
    for i in range(10):
        T = three_phase - 0.05 - 6.0 * i
        if T < 250.0:
            break
        st = pkg.flash([z, 1.0 - z], T=T, P=fixed["P"])
        if st.liquid_split and st.VF == 0.0:
            curve.append({"T": T, "x_light": st.x_light[0], "x_heavy": st.x_heavy[0]})
    return {"three_phase": three_phase, "curve": curve} if curve else None


def phase_envelope(components: list[str], fractions: list[float], method: str,
                   points: int = 24) -> dict:
    """Bubble and dew curves (T vs P) of a mixture with an equation of state.

    The curves are traced by pressure from 0.5 bar towards the critical region, where the
    flash stops converging; the last points are refined by halving the pressure step.
    """
    if method not in EOS_METHODS:
        raise SimulationError(
            f"a phase envelope needs an equation of state, not {method!r}",
            hint="use Peng-Robinson or SRK; activity models do not describe the critical region",
        )
    if len(fractions) != len(components) or any(f < 0 for f in fractions) or sum(fractions) <= 0:
        raise SimulationError("give one non-negative mole fraction per component",
                              hint="for example 0.6, 0.3, 0.1")
    pkg = _package(components, method)
    total = sum(fractions)
    zs = [f / total for f in fractions]
    curves = {kind: _trace(pkg, zs, kind, points) for kind in ("bubble", "dew")}
    allpts = curves["bubble"] + curves["dew"]
    return {
        "components": list(components),
        "fractions": zs,
        "method": pkg.method,
        "bubble": curves["bubble"],
        "dew": curves["dew"],
        "cricondenbar": max(allpts, key=lambda p: p["P"]) if allpts else None,
        "cricondentherm": max(allpts, key=lambda p: p["T"]) if allpts else None,
        "warnings": list(pkg.warnings),
    }


def _trace(pkg, zs, kind: str, points: int) -> list[dict]:
    ln_lo, ln_hi = math.log(5e4), math.log(3e7)
    step = (ln_hi - ln_lo) / points
    out: list[dict] = []
    ln_p, halvings = ln_lo, 0
    while ln_p <= ln_hi and halvings <= 3:
        try:
            sat = pkg.saturation_point(zs, kind, P=math.exp(ln_p))
            ok = sat.T > 0 and (not out or abs(sat.T - out[-1]["T"]) < 0.25 * out[-1]["T"])
        except SimulationError:
            ok = False
        if ok:
            out.append({"T": sat.T, "P": sat.P})
            ln_p += step
        elif out:  # approach the critical region in smaller steps
            ln_p -= step
            step /= 2.0
            ln_p += step
            halvings += 1
        else:
            ln_p += step
    return out
