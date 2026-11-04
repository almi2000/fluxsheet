"""Henry's law for light gases dissolved in an activity-coefficient liquid.

Above (or near) its critical temperature a gas such as nitrogen has no vapor pressure; the
extrapolated correlation that Raoult's law would use is meaningless. Its liquid fugacity is
instead ``f_i = x_i H_i`` with a Henry constant ``H_i`` measured in the solvent. Fluxsheet
hands the constants to ``thermo``'s ``GibbsExcessLiquid`` (``henry_components``), which puts
``H_i`` in place of the vapor pressure, mixes it over the solvents present
(``ln H_i = sum_j x_j ln H_ij / sum_j x_j``, Gmehling) and takes the heat of solution from
``d ln H_i / dT``. The activity model is wrapped by :class:`HenryDilutedModel`, which keeps
the solute's activity coefficient at 1 (the unsymmetric convention at infinite dilution).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cache

import numpy as np
from chemicals.iapws import iapws95_Psat
from thermo import GibbsExcessLiquid
from thermo.activity import GibbsExcess, IdealSolution
from thermo.interaction_parameters import IPDB
from thermo.phases import gibbs_excess as _ge

# Components whose critical temperature lies below this are Henry components by default in
# activity-coefficient methods: H2, He, N2, O2, Ar, CO, CH4, ethylene, CO2 (304 K) and ethane
# (305 K). They are supercritical, or nearly so, at any temperature where such a liquid model
# applies. H2S (373 K) and NH3 (405 K) condense and keep their vapor pressure unless listed.
HENRY_TC = 310.0

WATER = "7732-18-5"

# Tables shipped with thermo, best first (after IAPWS G7-04 for gases in water). Keys are
# "solute-CAS solvent-CAS"; ln H [Pa] = A + B/T + C ln T + D T + E/T^2 + F T^2 on the
# mole-fraction basis (H = f / x).
_TABLES = (
    ("Sander T dep", "Sander (2015) compilation, temperature-dependent"),
    ("ChemSep Henry", "ChemSep databank"),
    ("Sander Const", "Sander (2015) compilation, 25 °C value only (no heat of solution)"),
    ("PRTranslated Henry Best", "predicted with Peng-Robinson (thermo), not measured"),
)
_KEYS = ("A", "B", "C", "D", "E", "F")

# IAPWS G7-04 (2004), "Guideline on the Henry's constant and vapor-liquid distribution
# constant for gases in H2O and D2O at high temperatures" (Fernandez-Prini et al. 2003):
# ln(H / p_w*) = A / Tr + B tau^0.355 / Tr + C Tr^-0.41 exp(tau), Tr = T / 647.096 K,
# tau = 1 - Tr, p_w* the vapor pressure of water. Valid from the freezing point to near
# water's critical point. Gas CAS number -> (A, B, C).
IAPWS_G704 = {
    "7440-59-7": (-3.52839, 7.12983, 4.47770),  # helium
    "7440-01-9": (-3.18301, 5.31448, 5.43774),  # neon
    "7440-37-1": (-8.40954, 4.29587, 10.52779),  # argon
    "7439-90-9": (-8.97358, 3.61508, 11.29963),  # krypton
    "7440-63-3": (-14.21635, 4.00041, 15.60999),  # xenon
    "1333-74-0": (-4.73284, 6.08954, 6.06066),  # hydrogen
    "7727-37-9": (-9.67578, 4.72162, 11.70585),  # nitrogen
    "7782-44-7": (-9.44833, 4.43822, 11.42005),  # oxygen
    "630-08-0": (-10.52862, 5.13259, 12.01421),  # carbon monoxide
    "124-38-9": (-8.55445, 4.01195, 9.52345),  # carbon dioxide
    "7783-06-4": (-4.51499, 5.23538, 4.42126),  # hydrogen sulfide
    "74-82-8": (-10.44708, 4.66491, 12.12986),  # methane
    "74-84-0": (-19.67563, 4.51222, 20.62567),  # ethane
    "2551-62-4": (-16.56118, 2.15289, 20.35440),  # sulfur hexafluoride
}
_IAPWS_SOURCE = "IAPWS G7-04 (2004), temperature-dependent to 250 °C"


def iapws_henry(cas: str, T: float) -> float:
    """Henry constant of a gas in water by IAPWS G7-04, Pa."""
    A, B, C = IAPWS_G704[cas]
    Tr = T / 647.096
    tau = 1.0 - Tr
    return iapws95_Psat(T) * math.exp(A / Tr + B * tau**0.355 / Tr
                                      + C * Tr**-0.41 * math.exp(tau))


@cache
def _iapws_fit(cas: str) -> tuple[float, ...]:
    """IAPWS G7-04 as ln H = A + B/T + C ln T + D T, the form thermo evaluates; fitted from
    0 to 250 °C (deviation below 1 % there)."""
    Ts = np.linspace(273.16, 523.15, 60)
    X = np.column_stack([np.ones_like(Ts), 1.0 / Ts, np.log(Ts), Ts])
    y = np.log([iapws_henry(cas, T) for T in Ts])
    coeffs = np.linalg.lstsq(X, y, rcond=None)[0]
    return (*(float(c) for c in coeffs), 0.0, 0.0)


@dataclass(frozen=True)
class HenryConstant:
    coeffs: tuple[float, float, float, float, float, float]
    source: str

    def value(self, T: float) -> float:
        a, b, c, d, e, f = self.coeffs
        return math.exp(a + b / T + c * math.log(T) + d * T + e / T**2 + f * T**2)


def henry_constant(solute_cas: str, solvent_cas: str) -> HenryConstant | None:
    """The best available Henry constant of a solute in a solvent, or None."""
    if solvent_cas == WATER and solute_cas in IAPWS_G704:
        return HenryConstant(_iapws_fit(solute_cas), _IAPWS_SOURCE)
    key = f"{solute_cas} {solvent_cas}"
    found = []
    for table, source in _TABLES:
        # The Peng-Robinson predictions are far off in water (N2: 1e12 Pa against 8.6e9).
        if table == "PRTranslated Henry Best" and solvent_cas == WATER:
            continue
        row = IPDB.tables[table].get(key)
        if row is not None:
            found.append(HenryConstant(tuple(float(row[k]) for k in _KEYS), source))
    # Sander's single van't Hoff term misses the maximum of H near 100 °C (O2 in water is
    # 40 % too high at 80 °C); ChemSep's form with ln T has it. ChemSep has errors of a
    # factor of 5 to 10 for some gases (N2, H2), so it is taken only where it agrees with
    # Sander at 25 °C.
    sander_vant_hoff = len(found) > 1 and found[0].coeffs[2:] == (0.0,) * 4
    if (sander_vant_hoff and found[1].source == "ChemSep databank"
            and abs(math.log(found[1].value(298.15) / found[0].value(298.15))) < math.log(1.1)):
        return found[1]
    return found[0] if found else None


def raoult_estimate(vapor_pressure) -> HenryConstant:
    """A stand-in Henry constant from the extrapolated vapor pressure (ideal solution).

    ``ln H = A + B / T`` through the correlation at 0 and 100 °C. Only for solvents that
    have no measured constant, so that they still count in the mixing rule.
    """
    T1, T2 = 273.15, 373.15
    l1, l2 = math.log(vapor_pressure(T1)), math.log(vapor_pressure(T2))
    B = (l1 - l2) / (1.0 / T1 - 1.0 / T2)
    return HenryConstant((l1 - B / T1, B, 0.0, 0.0, 0.0, 0.0),
                         "estimated from the extrapolated vapor pressure (ideal solution)")


class HenryLiquid(GibbsExcessLiquid):
    """GibbsExcessLiquid whose Henry constants survive a liquid without solvent.

    thermo weights the solvents' constants by their mole fractions and divides by their sum,
    which fails for a liquid of dissolved gases only (a trial phase in a stability test, a
    pure-gas stream). There the solvents are weighted equally instead.
    """

    def _henry_weights(self) -> list[float]:
        zs, h = self.zs, self.henry_components
        if sum(z for z, hi in zip(zs, h) if not hi) > 1e-12:
            return zs
        return [0.0 if hi else 1.0 for hi in h]

    def Henry_constants(self):
        Hs = [0.0] * self.N
        _ge.Henry_constants(self.lnHenry_matrix(), self._henry_weights(), self.henry_components,
                            True, Hs)
        return Hs

    def dHenry_constants_dT(self):
        dHs = [0.0] * self.N
        _ge.dHenry_constants_dT(self.lnHenry_matrix(), self.dlnHenry_matrix_dT(),
                                self._henry_weights(), self.henry_components, True, dHs)
        return dHs

    def d2Henry_constants_dT2(self):
        d2Hs = [0.0] * self.N
        _ge.d2Henry_constants_dT2(self.lnHenry_matrix(), self.dlnHenry_matrix_dT(),
                                  self.d2lnHenry_matrix_dT2(), self._henry_weights(),
                                  self.henry_components, True, d2Hs)
        return d2Hs


class HenryDilutedModel(GibbsExcess):
    """An activity model of the solvents, extended by solutes that behave as ideal dilute.

    ``G^E(x) = S g(x_s / S)``, with ``S`` the solvents' total mole fraction and ``g`` the
    model of the solute-free solvent mixture. ``n G^E`` then depends on the solvent amounts
    only, so ``ln gamma`` of a solute is exactly 0 (``gamma* = gamma / gamma_inf = 1``: its
    Henry constant already holds the solute-solvent interaction), and the solvents keep the
    activity coefficients of their solute-free mixture. This is what an unsymmetric model
    with no solute-solvent parameters gives in a single solvent; salting-out and the solute's
    effect on solvent activities, both second order in a dilute liquid, are neglected.
    """

    def __init__(self, solvent_model, solvents: list[int], N: int, T: float = 298.15,
                 xs: list[float] | None = None):
        self.inner = solvent_model
        self.solvents = list(solvents)
        self.N = N
        self.vectorized = False
        self.model_id = solvent_model.model_id
        xs = xs if xs is not None else [1.0 / N] * N
        self.T, self.xs = T, list(xs)
        self.S = sum(xs[i] for i in self.solvents)
        self.ys = ([xs[i] / self.S for i in self.solvents] if self.S > 0.0
                   else [1.0 / len(self.solvents)] * len(self.solvents))
        self.inner = solvent_model.to_T_xs(T, self.ys)

    def __repr__(self) -> str:
        return f"HenryDilutedModel({self.inner!r}, solvents={self.solvents}, N={self.N})"

    def model_hash(self) -> int:
        return hash((type(self).__name__, self.N, tuple(self.solvents), self.inner.model_hash()))

    def to_T_xs(self, T: float, xs: list[float]) -> HenryDilutedModel:
        return HenryDilutedModel(self.inner, self.solvents, self.N, T, xs)

    def _full(self, values) -> list[float]:
        out = [0.0] * self.N
        for i, v in zip(self.solvents, values):
            out[i] = float(v)
        return out

    def GE(self) -> float:
        return self.S * self.inner.GE()

    def dGE_dT(self) -> float:
        return self.S * self.inner.dGE_dT()

    def d2GE_dT2(self) -> float:
        return self.S * self.inner.d2GE_dT2()

    def d3GE_dT3(self) -> float:
        return self.S * self.inner.d3GE_dT3()

    def _dx(self, g: float, g_k) -> list[float]:
        # d/dx_k of S g(x_s/S) = g + g_k - sum_m y_m g_m for a solvent k, 0 for a solute.
        mean = sum(y * gk for y, gk in zip(self.ys, g_k))
        return self._full([g + gk - mean for gk in g_k])

    def dGE_dxs(self) -> list[float]:
        return self._dx(self.inner.GE(), self.inner.dGE_dxs())

    def d2GE_dTdxs(self) -> list[float]:
        return self._dx(self.inner.dGE_dT(), self.inner.d2GE_dTdxs())

    def d2GE_dxixjs(self) -> list[list[float]]:
        # (g_kl - sum_n g_kn y_n - sum_m y_m g_ml + sum_mn y_m y_n g_mn) / S for solvents k, l.
        h, ys, n = self.inner.d2GE_dxixjs(), self.ys, len(self.ys)
        row = [sum(h[k][m] * ys[m] for m in range(n)) for k in range(n)]
        col = [sum(ys[m] * h[m][k] for m in range(n)) for k in range(n)]
        both = sum(ys[k] * row[k] for k in range(n))
        out = [[0.0] * self.N for _ in range(self.N)]
        if self.S <= 1e-12:
            return out
        for a, k in enumerate(self.solvents):
            for b, m in enumerate(self.solvents):
                out[k][m] = (h[a][b] - row[a] - col[b] + both) / self.S
        return out


def ideal_solvent_model(n: int, T: float = 298.15):
    """The solvent model when the method has none (a single solvent)."""
    return IdealSolution(T=T, xs=[1.0 / n] * n)
