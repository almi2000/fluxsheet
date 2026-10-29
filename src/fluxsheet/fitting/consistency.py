"""Thermodynamic consistency tests for binary T-P-x-y data.

Both tests check the measurements against the Gibbs-Duhem equation, which ties the two
activity coefficients of a binary together.

* **Point test** (Van Ness, Byer and Gibbs 1973; criterion of Fredenslund et al. 1977):
  a four-term Redlich-Kister expansion of G^E/RT is fitted to the bubble temperatures
  (isobaric data) or pressures (isothermal data) alone, by Barker's method, and the vapor
  compositions it predicts are compared with the measured ones. A set passes when the mean
  absolute deviation in y is at most 0.01.
* **Area test** (Redlich and Kister 1948; Herington 1951 for isobaric data): with
  experimental activity coefficients ``gamma_i = y_i P / (x_i Psat_i)``, the integral of
  ``ln(gamma1/gamma2)`` over x1 from 0 to 1 must vanish. With ``D = 100 |A+ - A-| / (A+ + A-)``,
  isothermal data pass when D <= 10, isobaric data when D - J <= 10 with
  ``J = 150 (Tmax - Tmin) / Tmin``. It needs data spanning most of the composition range.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

from .models import Model, Pair
from .thermoml import VLEDataset

POINT_TEST_MAX_AAD_Y = 0.01
AREA_TEST_MAX = 10.0


@dataclass(frozen=True)
class TestResult:
    applicable: bool
    passed: bool | None  # None when not applicable
    value: float | None  # point test: AAD y; area test: D (isothermal) or D - J (isobaric)
    note: str = ""

    def to_dict(self) -> dict:
        return {"passed": self.passed, "value": None if self.value is None else round(self.value, 4),
                **({"note": self.note} if self.note else {})}


def _interior(ds: VLEDataset):
    pts = [p for p in ds.points if 0.0 < p.x1 < 1.0 and p.y1 is not None and 0.0 < p.y1 < 1.0]
    return (np.array([p.T for p in pts]), np.array([p.P for p in pts]),
            np.array([p.x1 for p in pts]), np.array([p.y1 for p in pts]))


def point_test(pair: Pair, ds: VLEDataset, max_aad_y: float = POINT_TEST_MAX_AAD_Y) -> TestResult:
    if not ds.has_vapor:
        return TestResult(False, None, None, "no vapor compositions")
    T, P, x, y = _interior(ds)
    if len(x) < 5:
        return TestResult(False, None, None, "fewer than 5 mixture points")
    model = Model("redlich-kister", pair)
    n_terms = 4 if len(x) >= 8 else 3
    isothermal = ds.kind == "isothermal"

    def calc(A):
        if isothermal:
            return model.bubble_P(A, T, x)
        return model.bubble_T(A, P, x, T0=T)

    def residuals(A):
        value, _ = calc(A)
        return (value / P - 1.0) * 100.0 if isothermal else value - T

    try:
        fit = least_squares(residuals, np.zeros(n_terms), method="lm")
        _, y_calc = calc(fit.x)
    except (ValueError, FloatingPointError, OverflowError):
        return TestResult(False, None, None, "Redlich-Kister fit failed")
    aad = float(np.mean(np.abs(y_calc - y)))
    if not np.isfinite(aad):
        return TestResult(False, None, None, "Redlich-Kister fit failed")
    return TestResult(True, bool(aad <= max_aad_y), aad)


def area_test(pair: Pair, ds: VLEDataset) -> TestResult:
    if not ds.has_vapor:
        return TestResult(False, None, None, "no vapor compositions")
    T, P, x, y = _interior(ds)
    if len(x) < 5 or x.min() > 0.1 or x.max() < 0.9:
        return TestResult(False, None, None, "data do not span x1 = 0.1 to 0.9")
    order = np.argsort(x)
    T, P, x, y = T[order], P[order], x[order], y[order]
    ln_ratio = (np.log(y * P / (x * np.exp(pair.ln_psat(0, T)[0])))
                - np.log((1 - y) * P / ((1 - x) * np.exp(pair.ln_psat(1, T)[0]))))
    # Trapezoids from x = 0 to 1, holding the end values flat to the pure-component limits.
    xs = np.concatenate([[0.0], x, [1.0]])
    fs = np.concatenate([[ln_ratio[0]], ln_ratio, [ln_ratio[-1]]])
    area_pos = area_neg = 0.0
    for i in range(len(xs) - 1):
        x0, x1, f0, f1 = xs[i], xs[i + 1], fs[i], fs[i + 1]
        if f0 * f1 >= 0:
            a = 0.5 * (f0 + f1) * (x1 - x0)
            area_pos, area_neg = (area_pos + a, area_neg) if a >= 0 else (area_pos, area_neg - a)
        else:  # split the trapezoid where the curve crosses zero
            xc = x0 + (x1 - x0) * f0 / (f0 - f1)
            for a in (0.5 * f0 * (xc - x0), 0.5 * f1 * (x1 - xc)):
                area_pos, area_neg = (area_pos + a, area_neg) if a >= 0 else (area_pos, area_neg - a)
    if area_pos + area_neg <= 0:
        return TestResult(False, None, None, "ideal solution")
    D = 100.0 * abs(area_pos - area_neg) / (area_pos + area_neg)
    if ds.kind == "isothermal":
        return TestResult(True, bool(D <= AREA_TEST_MAX), float(D))
    J = 150.0 * (T.max() - T.min()) / T.min()
    return TestResult(True, bool(D - J <= AREA_TEST_MAX), float(D - J), f"D = {D:.1f}, J = {J:.1f}")
