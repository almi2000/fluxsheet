"""Bubble and dew points with the K-values that belong to them."""

from __future__ import annotations

import math

import numpy as np

from ..thermo_package import PropertyPackage


def _normalized(v: np.ndarray) -> np.ndarray:
    return v / v.sum()


def _bubble_secant(pkg: PropertyPackage, x: np.ndarray, P: float, T: float):
    """Bubble point from a nearby temperature: secant on ln(sum K x) against 1/T, which is
    nearly linear. Returns None if it does not converge or lands on the trivial K = 1."""
    y = x
    inv_prev = f_prev = None
    for _ in range(25):
        K = np.asarray(pkg.k_values(T, P, list(x), list(y)))
        s = float(np.sum(K * x))
        y = K * x / s
        f = math.log(s)
        if abs(f) < 1e-10:
            return None if np.max(np.abs(np.log(K))) < 1e-3 else (T, K)
        inv = 1.0 / T
        # Clausius-Clapeyron slope for ~30 kJ/mol until there are two points for a secant.
        slope = (f - f_prev) / (inv - inv_prev) if f_prev is not None and inv != inv_prev \
            else -3600.0
        if slope >= 0.0:
            return None
        inv_prev, f_prev = inv, f
        T = 1.0 / (inv - f / slope)
        if not 0.5 * (1.0 / inv) < T < 2.0 * (1.0 / inv):
            return None
    return None


def bubble_point(pkg: PropertyPackage, x, P: float,
                 T_guess: float | None = None) -> tuple[float, np.ndarray]:
    """Bubble temperature of liquid ``x`` at ``P``, and the K-values there. A ``T_guess``
    nearby makes it about ten times faster than a full flash."""
    x = np.asarray(x, dtype=float)
    if T_guess is not None:
        try:
            found = _bubble_secant(pkg, x, P, T_guess)
        except Exception:  # noqa: BLE001 - property failures fall back to the flash below
            found = None
        if found is not None:
            return found
    T = pkg.flash(list(x), P=P, VF=0.0).T
    K = np.asarray(pkg.k_values(T, P, list(x), list(x)))
    # The vapor fugacities depend on y; one more pass with the equilibrium vapor suffices.
    K = np.asarray(pkg.k_values(T, P, list(x), list(_normalized(K * x))))
    return T, K


def dew_point(pkg: PropertyPackage, y, P: float) -> tuple[float, np.ndarray]:
    """Dew temperature of vapor ``y`` at ``P``, and the K-values there."""
    y = np.asarray(y, dtype=float)
    T = pkg.flash(list(y), P=P, VF=1.0).T
    K = np.asarray(pkg.k_values(T, P, list(y), list(y)))
    K = np.asarray(pkg.k_values(T, P, list(_normalized(y / K)), list(y)))
    return T, K
