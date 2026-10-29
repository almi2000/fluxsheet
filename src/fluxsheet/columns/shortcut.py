"""Shortcut column design: Fenske, Underwood, Gilliland (Molokanov) and Kirkbride.

Stage counts follow the textbook convention: equilibrium stages including the reboiler (and
a partial condenser, which is an equilibrium stage too), excluding a total condenser.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from ..errors import SimulationError
from ..thermo_package import PropertyPackage
from .saturation import bubble_point, dew_point


def fenske_min_stages(d_lk: float, b_lk: float, d_hk: float, b_hk: float, alpha: float) -> float:
    """Minimum equilibrium stages at total reflux for a light/heavy key split."""
    return math.log((d_lk / b_lk) * (b_hk / d_hk)) / math.log(alpha)


def underwood_theta(alpha, z, q: float, lo: float, hi: float) -> float:
    """Root of sum(alpha_i z_i / (alpha_i - theta)) = 1 - q between two adjacent volatilities."""
    alpha, z = np.asarray(alpha), np.asarray(z)

    def f(theta: float) -> float:
        return float(np.sum(alpha * z / (alpha - theta))) - (1.0 - q)

    span = hi - lo
    return brentq(f, lo + 1e-10 * span, hi - 1e-10 * span, xtol=1e-14, rtol=1e-14)


def gilliland_molokanov(n_min: float, reflux: float, reflux_min: float) -> float:
    """Equilibrium stages at reflux ``reflux`` from the Molokanov form of Gilliland's chart."""
    X = (reflux - reflux_min) / (reflux + 1.0)
    Y = 1.0 - math.exp((1.0 + 54.4 * X) / (11.0 + 117.2 * X) * (X - 1.0) / math.sqrt(X))
    return (n_min + Y) / (1.0 - Y)


def kirkbride_ratio(z_lk: float, z_hk: float, xb_lk: float, xd_hk: float, B: float,
                    D: float) -> float:
    """Ratio of stages above the feed to stages below it (Kirkbride)."""
    return ((z_hk / z_lk) * (xb_lk / xd_hk) ** 2 * (B / D)) ** 0.206


@dataclass
class ShortcutResult:
    alpha: np.ndarray  # relative volatility to the heavy key, geometric mean of top and bottom
    min_stages: float
    min_reflux: float
    reflux_ratio: float
    stages: float  # equilibrium stages at reflux_ratio (Gilliland), not rounded
    rectifying_stages: float  # equilibrium stages above the feed (Kirkbride)
    q: float  # feed thermal condition: 1 = saturated liquid, 0 = saturated vapor
    distillate: np.ndarray  # component flows, mol/s
    bottoms: np.ndarray


def _split(feed, lk, hk, rec_lk, rec_hk, alpha, n_min):
    """Distillate flows: keys from their recoveries, the rest from the Fenske distribution."""
    r_hk = rec_hk / (1.0 - rec_hk)  # d/b of the heavy key
    log_ratio = np.log(r_hk) + n_min * np.log(alpha)
    frac = 1.0 / (1.0 + np.exp(-np.clip(log_ratio, -700.0, 700.0)))
    frac[lk], frac[hk] = rec_lk, rec_hk
    d = feed * frac
    return d, feed - d


def shortcut_design(
    pkg: PropertyPackage, feed, H_feed: float, P: float, lk: int, hk: int, rec_lk: float,
    rec_hk: float, *, reflux_ratio: float | None = None, reflux_factor: float | None = None,
) -> ShortcutResult:
    """Fenske-Underwood-Gilliland-Kirkbride design for one feed at column pressure ``P``.

    ``H_feed`` is the feed's molar enthalpy; ``rec_lk`` and ``rec_hk`` are the fractions of
    the light and heavy key that go to the distillate.
    """
    names = pkg.names
    feed = np.asarray(feed, dtype=float)
    F = feed.sum()
    z = feed / F
    if lk == hk:
        raise SimulationError("light and heavy key are the same component",
                              hint="choose two different components as keys")
    for i in (lk, hk):
        if feed[i] <= 0.0:
            raise SimulationError(f"key component {names[i]!r} is not in the feed",
                                  hint="choose keys that are present in the feed")
    if not 0.0 < rec_hk < rec_lk < 1.0:
        raise SimulationError(
            f"infeasible key recoveries: {names[lk]} {rec_lk:.4g} and {names[hk]} {rec_hk:.4g} "
            "to the distillate",
            hint="the light key's recovery to the distillate must be above the heavy key's, "
            "and both strictly between 0 and 1 (e.g. 0.99 and 0.01)",
        )

    _, K = bubble_point(pkg, z, P)
    alpha = K / K[hk]
    if alpha[lk] <= 1.0:
        raise SimulationError(
            f"{names[lk]} is not more volatile than {names[hk]} at column conditions "
            f"(relative volatility {alpha[lk]:.3g})",
            hint=f"swap the keys: make {names[hk]} the light key and {names[lk]} the heavy key",
        )

    # Non-key distribution and volatilities depend on each other; two passes are plenty.
    n_min = fenske_min_stages(rec_lk, 1 - rec_lk, rec_hk, 1 - rec_hk, alpha[lk])
    for _ in range(3):
        d, b = _split(feed, lk, hk, rec_lk, rec_hk, alpha, n_min)
        _, K_top = dew_point(pkg, d / d.sum(), P)
        _, K_bot = bubble_point(pkg, b / b.sum(), P)
        alpha = np.sqrt((K_top / K_top[hk]) * (K_bot / K_bot[hk]))
        if alpha[lk] <= 1.0:
            raise SimulationError(
                f"{names[lk]} is not more volatile than {names[hk]} between top and bottom "
                f"(mean relative volatility {alpha[lk]:.3g})",
                hint="check the key order, or whether the mixture forms an azeotrope",
            )
        n_min = fenske_min_stages(d[lk], b[lk], d[hk], b[hk], alpha[lk])
    D, B = d.sum(), b.sum()

    # Feed thermal condition from the enthalpies of saturated liquid and vapor of the feed.
    h_liq = pkg.flash(list(z), P=P, VF=0.0).H
    h_vap = pkg.flash(list(z), P=P, VF=1.0).H
    q = (h_vap - H_feed) / (h_vap - h_liq)

    # Underwood: theta lies between the heavy key and the next more volatile component.
    present = feed > 0.0
    above = alpha[present & (alpha > 1.0 + 1e-9)]
    theta = underwood_theta(alpha[present], z[present], q, 1.0, float(above.min()))
    x_d = d / D
    r_min = float(np.sum(alpha[present] * x_d[present] / (alpha[present] - theta))) - 1.0
    if r_min <= 0.0:
        raise SimulationError(
            f"Underwood gives a non-positive minimum reflux ({r_min:.3g}) for this split",
            hint="the separation is very easy or the feed is far from saturated; "
            "give an absolute reflux ratio, or bring the feed closer to its boiling point",
        )

    if reflux_ratio is None:
        if reflux_factor is None or reflux_factor <= 1.0:
            raise SimulationError(
                f"reflux must be more than the minimum; got {reflux_factor} times minimum",
                hint="use a reflux factor of 1.1 to 1.5 times minimum",
            )
        reflux_ratio = reflux_factor * r_min
    if reflux_ratio <= r_min:
        raise SimulationError(
            f"reflux ratio {reflux_ratio:.4g} is at or below the minimum {r_min:.4g}; "
            "the separation would need infinitely many stages",
            hint=f"use a reflux ratio of at least {1.1 * r_min:.3g} (1.1 times minimum), or "
            "relax the key recoveries",
        )
    n = gilliland_molokanov(n_min, reflux_ratio, r_min)
    ratio = kirkbride_ratio(z[lk], z[hk], b[lk] / B, d[hk] / D, B, D)
    n_rect = n * ratio / (1.0 + ratio)
    return ShortcutResult(alpha, n_min, r_min, reflux_ratio, n, n_rect, q, d, b)
