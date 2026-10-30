"""Least-squares regression of NRTL and UNIQUAC binary parameters to VLE data.

Objective
    For isobaric points the model's bubble temperature at the measured P and x is compared
    with the measured T; for isothermal points the bubble pressure at T and x with the
    measured P. Both also compare the vapor composition. Residuals are scaled by typical
    measurement uncertainties (``sigma_T`` 0.2 K, ``sigma_P`` 0.5 %, ``sigma_y`` 0.005) and
    minimised with ``scipy.optimize.least_squares``.

Parameters
    ``b12`` and ``b21`` (K) are always fitted. ``a12`` and ``a21`` are fitted as well only
    when the data span at least ``temperature_span`` (30 K); over a narrower range a and b
    cannot be told apart and ``a`` stays 0, as in the ChemSep tables. NRTL's alpha is fixed:
    0.3, or 0.2 when UNIFAC predicts a liquid-liquid split (Renon and Prausnitz).

Data selection, in this order
    1. Only complete T-P-x-y sets at pressures up to ``max_pressure`` (300 kPa, where the
       ideal-vapor model of the property packages holds) with at least 3 mixture points.
       Pure-component points (x = 0 or 1) carry no binary information and are left out.
    2. Each set must pass the point consistency test (mean |dy| <= 0.01); the area test is
       reported but does not decide, because it is unreliable for near-ideal mixtures.
    3. After a first fit, points are rejected when a scaled residual exceeds three times
       its robust spread (1.4826 x median absolute residual) and three measurement
       uncertainties, or when the liquid lies inside the model's liquid-liquid gap (the
       measured x is then an overall composition, not a liquid). The fit is then repeated.
"""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .consistency import TestResult, area_test, point_test
from .models import MODELS, Model, Pair
from .thermoml import VLEDataset

PARAM_NAMES = ("a12", "a21", "b12", "b21")
_BOUNDS = {"nrtl": (20.0, 6000.0), "uniquac": (10.0, 4000.0)}  # |a|, |b| limits


@dataclass
class FitOptions:
    model: str = "nrtl"
    alpha: float | None = None
    sigma_T: float = 0.2  # K
    sigma_P: float = 0.005  # relative
    sigma_y: float = 0.005
    temperature_span: float = 30.0  # K of data needed to fit a_ij as well as b_ij
    outlier_factor: float = 3.0
    max_pressure: float = 300e3  # Pa
    max_point_test_aad: float = 0.01
    require_consistency: bool = True


@dataclass
class DatasetReport:
    dataset: VLEDataset
    citation: str
    point_test: TestResult
    area_test: TestResult
    used: bool
    reason: str = ""
    n_rejected: int = 0  # outlier points of this set

    def to_dict(self) -> dict:
        ds = self.dataset
        return {"doi": ds.doi, "dataset": ds.number, "conditions": ds.label,
                "points": len(ds.points), "used": self.used, **({"reason": self.reason}
                                                                if self.reason else {}),
                "rejected_points": self.n_rejected, "point_test": self.point_test.to_dict(),
                "area_test": self.area_test.to_dict()}


@dataclass
class Points:
    T: np.ndarray
    P: np.ndarray
    x: np.ndarray
    y: np.ndarray
    iso: np.ndarray  # True: isothermal (P is computed), False: T is computed
    owner: np.ndarray  # index of the data set report

    def subset(self, mask) -> Points:
        return Points(*(getattr(self, f)[mask] for f in ("T", "P", "x", "y", "iso", "owner")))

    def __len__(self) -> int:
        return len(self.x)


@dataclass
class FitResult:
    model: str
    pair: Pair
    params: tuple[float, float, float, float]
    alpha: float | None
    free: tuple[str, ...]
    datasets: list[DatasetReport]
    points: Points  # the points the final fit used
    n_rejected: int
    stats: dict
    comparison: dict = field(default_factory=dict)  # other parameter sets on the same points

    def entry(self) -> dict:
        """The record stored in the parameter file."""
        a12, a21, b12, b21 = self.params
        used = [r for r in self.datasets if r.used]
        dois: dict[str, dict] = {}
        for r in used:
            d = dois.setdefault(r.dataset.doi, {"doi": r.dataset.doi, "citation": r.citation,
                                                "datasets": []})
            d["datasets"].append(r.dataset.label)
        return {
            "components": [{"name": n, "cas": c} for n, c in zip(self.pair.names, self.pair.cas)],
            "a12": round(a12, 6), "a21": round(a21, 6), "b12": round(b12, 4),
            "b21": round(b21, 4), **({"alpha": self.alpha} if self.model == "nrtl" else {}),
            "fitted": list(self.free),
            "T_range_K": [round(float(self.points.T.min()), 2),
                          round(float(self.points.T.max()), 2)],
            "P_range_kPa": [round(float(self.points.P.min()) / 1e3, 3),
                            round(float(self.points.P.max()) / 1e3, 3)],
            "n_points": len(self.points), "n_rejected_points": self.n_rejected,
            "n_datasets": len(used), "fit": self.stats, "comparison": self.comparison,
            "sources": list(dois.values()),
            "rejected_datasets": [r.to_dict() for r in self.datasets if not r.used],
            "fitted_on": dt.date.today().isoformat(),
        }


def _points(reports: list[DatasetReport]) -> Points:
    cols: dict[str, list] = {k: [] for k in ("T", "P", "x", "y", "iso", "owner")}
    for i, r in enumerate(reports):
        if not r.used:
            continue
        for p in r.dataset.points:
            if 0.0 < p.x1 < 1.0 and p.y1 is not None:
                for k, v in zip(cols, (p.T, p.P, p.x1, p.y1, r.dataset.kind == "isothermal", i)):
                    cols[k].append(v)
    return Points(*(np.array(v, dtype=bool if k == "iso" else None) for k, v in cols.items()))


def calculate(model: Model, p, pts: Points):
    """Model T or P and y at the measured conditions: (T_calc, P_calc, y_calc)."""
    T_calc, P_calc, y_calc = pts.T.copy(), pts.P.copy(), np.empty_like(pts.y)
    iso = pts.iso
    if iso.any():
        P_calc[iso], y_calc[iso] = model.bubble_P(p, pts.T[iso], pts.x[iso])
    if (~iso).any():
        T_calc[~iso], y_calc[~iso] = model.bubble_T(p, pts.P[~iso], pts.x[~iso], T0=pts.T[~iso])
    return T_calc, P_calc, y_calc


def _scaled(model: Model, p, pts: Points, opt: FitOptions):
    """Scaled residuals: (temperature-or-pressure, vapor composition), one of each per point."""
    T_calc, P_calc, y_calc = calculate(model, p, pts)
    r_tp = np.where(pts.iso, (P_calc / pts.P - 1.0) / opt.sigma_P, (T_calc - pts.T) / opt.sigma_T)
    r_y = (y_calc - pts.y) / opt.sigma_y
    return np.nan_to_num(r_tp, nan=1e3, posinf=1e3, neginf=-1e3), np.nan_to_num(
        r_y, nan=1e3, posinf=1e3, neginf=-1e3)


def statistics(model: Model, p, pts: Points) -> dict:
    """Average absolute and root-mean-square deviations in T (K), P (%) and y."""
    T_calc, P_calc, y_calc = calculate(model, p, pts)
    out: dict = {}
    iso = pts.iso
    if (~iso).any():
        dT = T_calc[~iso] - pts.T[~iso]
        out |= {"AAD_T_K": float(np.mean(np.abs(dT))), "RMS_T_K": float(np.sqrt(np.mean(dT**2)))}
    if iso.any():
        dP = 100.0 * (P_calc[iso] / pts.P[iso] - 1.0)
        out |= {"AAD_P_pct": float(np.mean(np.abs(dP))),
                "RMS_P_pct": float(np.sqrt(np.mean(dP**2)))}
    dy = y_calc - pts.y
    out |= {"AAD_y": float(np.mean(np.abs(dy))), "RMS_y": float(np.sqrt(np.mean(dy**2)))}
    return {k: round(v, 4 if k.endswith("_y") else 3) if np.isfinite(v) else None
            for k, v in out.items()}


def screen(pair: Pair, datasets: list[tuple[VLEDataset, str]],
           opt: FitOptions) -> list[DatasetReport]:
    """Apply the data-selection rules 1 and 2 to (data set, citation) pairs."""
    reports = []
    for ds, citation in datasets:
        pt, at = point_test(pair, ds, opt.max_point_test_aad), area_test(pair, ds)
        mixture = [p for p in ds.points if 0.0 < p.x1 < 1.0]
        reason = ""
        if not ds.has_vapor:
            reason = "no vapor compositions"
        elif max(p.P for p in ds.points) > opt.max_pressure:
            reason = f"pressure above {opt.max_pressure / 1e3:.0f} kPa"
        elif len(mixture) < 3:
            reason = "fewer than 3 mixture points"
        elif opt.require_consistency and not pt.passed:
            reason = ("failed the point consistency test" if pt.applicable
                      else f"point consistency test not possible ({pt.note})")
        reports.append(DatasetReport(ds, citation, pt, at, not reason, reason))
    return reports


def _start_points(model: Model, opt: FitOptions, fit_a: bool, T_mean: float,
                  extra: list) -> list[np.ndarray]:
    from ..thermo_package import _estimate_binary

    starts = [np.zeros(4)]
    est = _estimate_binary(model.kind, *model.pair.cas)
    for guess in ([est[:4]] if est else []) + list(extra):
        a12, a21, b12, b21 = guess
        if not fit_a:  # same tau at the mean temperature with a = 0
            a12, a21, b12, b21 = 0.0, 0.0, b12 + a12 * T_mean, b21 + a21 * T_mean
        starts.append(np.array([a12, a21, b12, b21], dtype=float))
    return starts


def _fit(model: Model, pts: Points, opt: FitOptions, fit_a: bool, starts) -> np.ndarray:
    free = [0, 1, 2, 3] if fit_a else [2, 3]
    a_max, b_max = _BOUNDS[model.kind]
    lim = np.array([a_max, a_max, b_max, b_max])[free]

    def residuals(q):
        p = np.zeros(4)
        p[free] = q
        return np.concatenate(_scaled(model, p, pts, opt))

    best = None
    for s in starts:
        q0 = np.clip(s[free], -0.99 * lim, 0.99 * lim)
        try:
            res = least_squares(residuals, q0, bounds=(-lim, lim),
                                x_scale=np.where(np.array(free) < 2, 1.0, 300.0))
        except (ValueError, FloatingPointError):
            continue
        if best is None or res.cost < best.cost:
            best = res
    if best is None:
        raise RuntimeError("the regression did not converge from any starting point")
    p = np.zeros(4)
    p[free] = best.x
    return p


def fit_binary(pair: Pair, datasets: list[tuple[VLEDataset, str]],
               opt: FitOptions | None = None, extra_starts=()) -> FitResult:
    """Screen the data sets, fit, reject outliers, refit. ``datasets`` are already oriented
    to ``pair`` (first component = component 1)."""
    opt = opt or FitOptions()
    if opt.model not in MODELS:
        raise ValueError(f"model must be one of {MODELS}")
    alpha = None
    if opt.model == "nrtl":
        from ..thermo_package import _binary_splits

        alpha = opt.alpha if opt.alpha is not None else (0.2 if _binary_splits(*pair.cas)
                                                          else 0.3)
    model = Model(opt.model, pair, alpha)
    reports = screen(pair, datasets, opt)
    pts = _points(reports)
    if len(pts) < 4:
        raise ValueError(f"only {len(pts)} usable points for {'/'.join(pair.names)}; "
                         "need at least 4")
    fit_a = float(pts.T.max() - pts.T.min()) >= opt.temperature_span
    starts = _start_points(model, opt, fit_a, float(pts.T.mean()), extra_starts)
    p = _fit(model, pts, opt, fit_a, starts)

    keep, n_rejected = _inliers(model, p, pts, opt), 0
    if not keep.all():
        kept = pts.subset(keep)
        if len(kept) >= 4:
            for i in pts.owner[~keep]:
                reports[i].n_rejected += 1
            p = _fit(model, kept, opt, fit_a, [p, *starts])
            pts, n_rejected = kept, int((~keep).sum())
    free = PARAM_NAMES if fit_a else PARAM_NAMES[2:]
    return FitResult(opt.model, pair, tuple(float(v) for v in p), alpha, free, reports, pts,
                     n_rejected, statistics(model, p, pts))


def _inliers(model: Model, p, pts: Points, opt: FitOptions) -> np.ndarray:
    r_tp, r_y = _scaled(model, p, pts, opt)
    keep = np.ones(len(pts), dtype=bool)
    for r in (r_tp, r_y):
        spread = 1.4826 * np.median(np.abs(r))
        keep &= np.abs(r) <= opt.outlier_factor * max(spread, 1.0)
    for T in np.unique(np.round(pts.T, 1)):
        at = np.abs(pts.T - T) < 0.05
        keep[at] &= ~model.in_miscibility_gap(p, float(T), pts.x[at])
    return keep


def compare(result: FitResult, label: str, params, alpha: float | None = None) -> None:
    """Add another parameter set's deviations on the same points to ``result.comparison``."""
    model = Model(result.model, result.pair, alpha if alpha is not None else result.alpha)
    result.comparison[label] = statistics(model, params, result.points)
