"""Binary activity-coefficient models and bubble points, vectorised over data points.

The phase-equilibrium model is the one Fluxsheet's activity-coefficient property packages
use: modified Raoult's law, ``y_i P = x_i gamma_i(T, x) Psat_i(T)``, with an ideal vapor
and the vapor-pressure correlations from ``thermo``. Fitting with the same model keeps the
regressed parameters consistent with the simulations that use them.

Parameters are always ``p = (a12, a21, b12, b21)`` for the component order of the pair:

* NRTL: ``tau_ij = a_ij + b_ij / T``, ``G_ij = exp(-alpha tau_ij)``
* UNIQUAC: ``tau_ij = exp(a_ij + b_ij / T)``

which is the convention of ``thermo.NRTL``/``thermo.UNIQUAC`` (``tau_as``, ``tau_bs``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

MODELS = ("nrtl", "uniquac")


def nrtl_ln_gammas(p, alpha: float, T, x1):
    a12, a21, b12, b21 = p
    x1 = np.asarray(x1, dtype=float)
    x2 = 1.0 - x1
    t12, t21 = a12 + b12 / T, a21 + b21 / T
    G12, G21 = np.exp(-alpha * t12), np.exp(-alpha * t21)
    d1, d2 = x1 + x2 * G21, x2 + x1 * G12
    ln1 = x2**2 * (t21 * (G21 / d1) ** 2 + t12 * G12 / d2**2)
    ln2 = x1**2 * (t12 * (G12 / d2) ** 2 + t21 * G21 / d1**2)
    return ln1, ln2


def uniquac_ln_gammas(p, rq, T, x1):
    a12, a21, b12, b21 = p
    (r1, r2), (q1, q2) = rq
    x1 = np.clip(np.asarray(x1, dtype=float), 1e-12, 1.0 - 1e-12)
    x2 = 1.0 - x1
    t12, t21 = np.exp(a12 + b12 / T), np.exp(a21 + b21 / T)
    phi1, phi2 = r1 * x1 / (r1 * x1 + r2 * x2), r2 * x2 / (r1 * x1 + r2 * x2)
    th1, th2 = q1 * x1 / (q1 * x1 + q2 * x2), q2 * x2 / (q1 * x1 + q2 * x2)
    l1, l2 = 5.0 * (r1 - q1) - (r1 - 1.0), 5.0 * (r2 - q2) - (r2 - 1.0)
    s1, s2 = th1 + th2 * t21, th2 + th1 * t12
    ln1 = (np.log(phi1 / x1) + 5.0 * q1 * np.log(th1 / phi1) + phi2 * (l1 - r1 / r2 * l2)
           - q1 * np.log(s1) + th2 * q1 * (t21 / s1 - t12 / s2))
    ln2 = (np.log(phi2 / x2) + 5.0 * q2 * np.log(th2 / phi2) + phi1 * (l2 - r2 / r1 * l1)
           - q2 * np.log(s2) + th1 * q2 * (t12 / s2 - t21 / s1))
    return ln1, ln2


def redlich_kister_ln_gammas(A, T, x1):
    """G^E/RT = x1 x2 sum_k A_k (x1 - x2)^k, temperature independent (for consistency tests)."""
    x1 = np.asarray(x1, dtype=float)
    x2, d = 1.0 - x1, 2.0 * x1 - 1.0
    S = sum(a * d**k for k, a in enumerate(A))
    dS = sum(k * a * d ** (k - 1) for k, a in enumerate(A) if k)
    g = x1 * x2 * S
    dg = (x2 - x1) * S + 2.0 * x1 * x2 * dS  # d(g)/d(x1) with x2 = 1 - x1
    return g + x2 * dg, g - x1 * dg


@dataclass
class Pair:
    """Pure-component data of a binary: vapor pressures, UNIQUAC sizes, names, CAS numbers."""

    names: tuple[str, str]
    cas: tuple[str, str]
    vapor_pressures: tuple  # thermo VaporPressure objects
    rq: tuple[tuple[float, float], tuple[float, float]] | None  # ((r1, r2), (q1, q2))
    Tbs: tuple[float | None, float | None]

    @classmethod
    def from_cas(cls, cas_a: str, cas_b: str) -> Pair:
        from ..thermo_package import _constants

        consts, props = _constants((cas_a, cas_b))
        rs, qs = consts.UNIFAC_Rs, consts.UNIFAC_Qs
        rq = None if None in (*rs, *qs) else (tuple(rs), tuple(qs))
        return cls(tuple(consts.names), (cas_a, cas_b), tuple(props.VaporPressures), rq,
                   tuple(consts.Tbs))

    def __post_init__(self):
        self._splines: dict[int, tuple[float, float, CubicSpline]] = {}

    def ln_psat(self, i: int, T):
        """ln Psat [Pa] of component ``i`` at ``T`` [K] and its T derivative, vectorised.

        The correlation is tabulated once on a 0.5 K grid and splined, which keeps the
        thousands of bubble-point evaluations of a fit fast without changing the values.
        """
        T = np.asarray(T, dtype=float)
        lo, hi, spline = self._splines.get(i, (np.inf, -np.inf, None))
        t_min, t_max = float(np.min(T)), float(np.max(T))
        if spline is None or t_min < lo or t_max > hi:
            lo, hi = min(lo, t_min - 60.0), max(hi, t_max + 60.0)
            lo = max(lo, 30.0)
            grid = np.arange(lo, hi + 0.5, 0.5)
            vp = self.vapor_pressures[i]
            values = np.log([max(vp(t), 1e-30) for t in grid])
            spline = CubicSpline(grid, values)
            self._splines[i] = (lo, hi, spline)
        return spline(T), spline(T, 1)


class Model:
    """An activity-coefficient model for one pair, evaluated with parameters ``p``."""

    def __init__(self, kind: str, pair: Pair, alpha: float | None = None):
        if kind not in (*MODELS, "redlich-kister"):
            raise ValueError(f"unknown model {kind!r}")
        if kind == "uniquac" and pair.rq is None:
            raise ValueError(f"UNIQUAC needs r and q (UNIFAC groups) for {' and '.join(pair.names)}")
        self.kind, self.pair, self.alpha = kind, pair, 0.3 if alpha is None else alpha

    def ln_gammas(self, p, T, x1):
        if self.kind == "nrtl":
            return nrtl_ln_gammas(p, self.alpha, T, x1)
        if self.kind == "uniquac":
            return uniquac_ln_gammas(p, self.pair.rq, T, x1)
        return redlich_kister_ln_gammas(p, T, x1)

    def bubble_P(self, p, T, x1):
        """Bubble pressure [Pa] and vapor mole fraction y1 at given T and x1."""
        T, x1 = np.asarray(T, dtype=float), np.asarray(x1, dtype=float)
        ln1, ln2 = self.ln_gammas(p, T, x1)
        p1 = x1 * np.exp(ln1 + self.pair.ln_psat(0, T)[0])
        p2 = (1.0 - x1) * np.exp(ln2 + self.pair.ln_psat(1, T)[0])
        P = p1 + p2
        return P, p1 / P

    def bubble_T(self, p, P, x1, T0=None, iterations: int = 50):
        """Bubble temperature [K] and y1 at given P [Pa] and x1, by damped Newton steps.

        ``T0`` (e.g. the measured temperatures) is the start; otherwise a mole-fraction
        average of the pure boiling points at P. The gamma T-dependence is left out of the
        derivative, which costs a few iterations but never the answer.
        """
        P, x1 = np.asarray(P, dtype=float), np.asarray(x1, dtype=float)
        if T0 is None:
            tb = [t if t is not None else 350.0 for t in self.pair.Tbs]
            T = x1 * tb[0] + (1.0 - x1) * tb[1]
        else:
            T = np.array(T0, dtype=float)
        T = np.broadcast_to(T, np.broadcast(P, x1).shape).copy()
        lnP = np.log(P)
        for _ in range(iterations):
            ln1, ln2 = self.ln_gammas(p, T, x1)
            (s1, d1), (s2, d2) = self.pair.ln_psat(0, T), self.pair.ln_psat(1, T)
            p1, p2 = x1 * np.exp(ln1 + s1), (1.0 - x1) * np.exp(ln2 + s2)
            total = p1 + p2
            f = np.log(total) - lnP
            slope = (p1 * d1 + p2 * d2) / total
            step = np.clip(-f / slope, -25.0, 25.0)
            T = np.clip(T + step, 30.0, 2000.0)
            if np.all(np.abs(step) < 1e-7):
                break
        ln1, ln2 = self.ln_gammas(p, T, x1)
        p1 = x1 * np.exp(ln1 + self.pair.ln_psat(0, T)[0])
        p2 = (1.0 - x1) * np.exp(ln2 + self.pair.ln_psat(1, T)[0])
        return T, p1 / (p1 + p2)

    def in_miscibility_gap(self, p, T: float, x1, n: int = 401):
        """Whether each liquid composition lies inside the model's liquid-liquid gap at T.

        A point is in the gap when the Gibbs energy of mixing there lies above its lower
        convex hull (the common tangent of the two coexisting liquids).
        """
        xs = np.linspace(1e-6, 1.0 - 1e-6, n)
        ln1, ln2 = self.ln_gammas(p, T, xs)
        g = xs * (np.log(xs) + ln1) + (1.0 - xs) * (np.log(1.0 - xs) + ln2)
        hull = _lower_hull(xs, g)
        above = g - np.interp(xs, xs[hull], g[hull]) > 1e-9
        return np.interp(np.asarray(x1, dtype=float), xs, above.astype(float)) > 0.5


def _lower_hull(x, y) -> list[int]:
    hull: list[int] = []
    for i in range(len(x)):
        while len(hull) >= 2:
            j, k = hull[-2], hull[-1]
            if (y[k] - y[j]) * (x[i] - x[j]) >= (y[i] - y[j]) * (x[k] - x[j]):
                hull.pop()
            else:
                break
        hull.append(i)
    return hull
