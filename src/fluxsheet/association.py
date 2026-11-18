"""Vapor-phase dimerization of carboxylic acids (chemical theory) for the activity methods.

Carboxylic acids form hydrogen-bonded dimers in the vapor, 2 A = A2. Near the boiling point
most of the acid in the saturated vapor is dimerized (acetic acid at 1 atm: 75 % of the
molecules), so an ideal gas of monomers gets the acid's heat of vaporization 40 to 65 % too
high, its vapor density too low and, through the formation enthalpies, liquid-phase heats of
reaction about 10 kJ per mole of acid off.

Model (Prausnitz, Lichtenthaler and de Azevedo, *Molecular Thermodynamics of Fluid-Phase
Equilibria*, ch. 7): the vapor is an ideal gas of true species, the monomers of every
component plus the homodimers of the acids, in equilibrium

    p_D = K(T) p_A^2,   log10 K = A + B / T  (K in 1/mmHg),

with the published constants in ``DIMERIZATION``. Its dimerization enthalpy is
dH = -R ln(10) B, about -60 kJ per mole of dimer. Per mole of the apparent (stoichiometric)
composition z, with n true moles and m_i, d_i the monomer and dimer amounts of acid i:

    z_i = m_i + 2 d_i,   n = 1 - sum d_i,   phi_i = y_mi / z_i (acids), 1 / n (others)
    H = H_ig + sum d_i dH_i,   G = G_ig + R T sum z_i ln phi_i,   V = n R T / P.

The fugacity of an acid is that of its monomer, f_i = y_mi P. The pure liquid is therefore
in equilibrium with the monomer pressure of the saturated vapor, p_m with p_m + K p_m^2 =
Psat, not with Psat itself: :class:`MonomerVaporPressure` gives the liquid that reference
fugacity, and :class:`AssociatedLatentHeat` the enthalpy between the liquid and the ideal
monomer gas, Hvap - d_sat dH (the measured heat of vaporization plus the bond enthalpy of
the dimers that form). The saturation pressure, the measured heat of vaporization and the
formation enthalpies are all met at once.

Not described: cross dimers between different acids and acid-water complexes; a mixture
of several acids is treated as each dimerizing with itself only.
"""

from __future__ import annotations

import math

from scipy.optimize import brentq
from thermo import FlashPureVLS, IdealGas

R = 8.314462618  # J/(mol K)
LN10 = math.log(10.0)
MMHG = 133.322368  # Pa

# CAS -> (A, B) of log10 K = A + B/T, K in 1/mmHg, from vapor density measurements.
DIMERIZATION = {
    # Marek and Standart (1954), Collect. Czech. Chem. Commun. 19, 1074, the values of
    # Ritter and Simons; dH = -60.6 kJ/mol, about 75 % of the molecules dimerized in the
    # saturated vapor at the normal boiling point.
    "64-19-7": (-10.4205, 3166.0),  # acetic acid
    # Coolidge (1928), J. Am. Chem. Soc. 50, 2166; dH = -59.0 kJ/mol.
    "64-18-6": (-10.743, 3083.0),  # formic acid
    # Vapor-density fit 50-150 degC (dissociation constant inverted); dH = -63.5 kJ/mol.
    "79-09-4": (-10.834, 3316.0),  # propionic acid
}


def ln_k(constants: tuple[float, float], T: float) -> float:
    """ln K of the dimerization, K in 1/Pa."""
    A, B = constants
    return LN10 * (A - math.log10(MMHG) + B / T)


def dimerization_enthalpy(constants: tuple[float, float]) -> float:
    """J per mole of dimer formed (negative)."""
    return -R * LN10 * constants[1]


def _ln_k_slope(constants: tuple[float, float], T: float) -> float:
    return -LN10 * constants[1] / (T * T)


def saturation_monomer(constants, T: float, Psat: float) -> tuple[float, float]:
    """Monomer pressure (Pa) of the pure saturated vapor and its dimers per mole of acid."""
    K = math.exp(ln_k(constants, T))
    p_m = 2.0 * Psat / (1.0 + math.sqrt(1.0 + 4.0 * K * Psat))
    u = K * p_m
    return p_m, u / (1.0 + 2.0 * u)


class MonomerVaporPressure:
    """Reference fugacity of the pure liquid acid: the monomer pressure of its saturated
    vapor. Stands in for the vapor pressure correlation in the liquid phase."""

    def __init__(self, vapor_pressure, constants: tuple[float, float]):
        self.base, self.constants = vapor_pressure, constants

    def _values(self, T: float) -> tuple[float, float, float]:
        """p_m and its first and second temperature derivatives."""
        base = self.base
        Ps = base.T_dependent_property(T)
        dPs = base.T_dependent_property_derivative(T)
        d2Ps = base.T_dependent_property_derivative(T, order=2)
        K = math.exp(ln_k(self.constants, T))
        a = _ln_k_slope(self.constants, T)  # d ln K / dT
        dK = K * a
        d2K = K * (a * a + 2.0 * LN10 * self.constants[1] / T**3)
        p = 2.0 * Ps / (1.0 + math.sqrt(1.0 + 4.0 * K * Ps))
        # Psat = p + K p^2, differentiated once and twice.
        denominator = 1.0 + 2.0 * K * p
        dp = (dPs - dK * p * p) / denominator
        d2p = (d2Ps - d2K * p * p - 4.0 * dK * p * dp - 2.0 * K * dp * dp) / denominator
        return p, dp, d2p

    def __call__(self, T: float) -> float:
        return self._values(T)[0]

    def T_dependent_property(self, T: float) -> float:
        return self._values(T)[0]

    def T_dependent_property_derivative(self, T: float, order: int = 1) -> float:
        return self._values(T)[order]

    def T_dependent_property_transform(self, T: float, transform: str) -> float:
        p, dp, d2p = self._values(T)
        if transform == "log_transform":
            return math.log(p)
        if transform in ("logarithmic_derivative", "derivative_over_value"):
            return dp / p
        if transform == "second_logarithmic_derivative":
            return d2p / p - (dp / p) ** 2
        if transform == "second_derivative_over_value":
            return d2p / p
        raise ValueError(f"unknown transform {transform!r}")

    def solve_property(self, value: float) -> float:
        return brentq(lambda T: self(T) - value, 0.3 * self.base.Tc, 1.2 * self.base.Tc)


class AssociatedLatentHeat:
    """Enthalpy from the liquid acid to its ideal monomer gas, J/mol: the measured heat of
    vaporization plus the bond enthalpy of the dimers in the saturated vapor."""

    def __init__(self, enthalpy_vaporization, vapor_pressure, constants: tuple[float, float]):
        self.base, self.vapor_pressure, self.constants = (enthalpy_vaporization,
                                                          vapor_pressure, constants)
        self.dH = dimerization_enthalpy(constants)
        self._monomer = MonomerVaporPressure(vapor_pressure, constants)

    def _dimers(self, T: float) -> tuple[float, float]:
        """Dimers per mole of acid in the saturated vapor, and their T derivative."""
        p, dp, _ = self._monomer._values(T)
        K = math.exp(ln_k(self.constants, T))
        u = K * p
        du = K * _ln_k_slope(self.constants, T) * p + K * dp
        return u / (1.0 + 2.0 * u), du / (1.0 + 2.0 * u) ** 2

    def __call__(self, T: float) -> float:
        measured = self.base(T) or 0.0
        return measured - self._dimers(T)[0] * self.dH

    def T_dependent_property_derivative(self, T: float) -> float:
        slope = self.base.T_dependent_property_derivative(T) if self.base(T) else 0.0
        return (slope or 0.0) - self._dimers(T)[1] * self.dH


class AssociatingGas(IdealGas):
    """Ideal gas of monomers and acid dimers in chemical equilibrium, on the apparent
    (stoichiometric) composition. ``acids`` maps component index to dimerization constants.
    """

    def __init__(self, acids: dict[int, tuple[float, float]], **kwargs):
        super().__init__(**kwargs)
        self.acids = tuple(sorted(acids.items()))

    def to_TP_zs(self, T, P, zs):
        new = super().to_TP_zs(T, P, zs)
        new.acids = self.acids
        return new

    def to(self, zs, T=None, P=None, V=None):
        if V is not None:  # n R T / P = V has no closed form in T or P
            if T is not None and P is None:
                P = brentq(lambda P: self.to_TP_zs(T, P, zs).V() - V, 1e-3 * self.R * T / V,
                           self.R * T / V)
            elif P is not None and T is None:
                T = brentq(lambda T: self.to_TP_zs(T, P, zs).V() - V, P * V / self.R,
                           2.0 * P * V / self.R)
            elif T is None or P is None:
                raise ValueError("Two of T, P, or V are needed")
        return self.to_TP_zs(T, P, zs)

    # ---- chemical equilibrium ------------------------------------------------------------

    def _equilibrium_at(self, T: float, P: float):
        """(true moles n per apparent mole, {acid: (monomer fraction y_m, dimers d)})."""
        zs = self.zs
        terms = [(i, math.exp(ln_k(c, T)) * P, zs[i]) for i, c in self.acids if zs[i] > 0.0]
        if not terms:
            return 1.0, {}

        def monomer(kp, z, n):
            q = z / n
            return 2.0 * q / (1.0 + math.sqrt(1.0 + 8.0 * kp * q))

        def excess(n):
            return n * (1.0 + sum(kp * monomer(kp, z, n) ** 2 for _, kp, z in terms)) - 1.0

        lowest = 1.0 - 0.5 * sum(z for _, _, z in terms)
        n = 1.0 if excess(1.0) <= 0.0 else brentq(excess, lowest, 1.0, xtol=1e-16,
                                                   rtol=1e-15, maxiter=200)
        species = {}
        for i, kp, z in terms:
            y = monomer(kp, z, n)
            species[i] = (y, n * kp * y * y)
        return n, species

    def _equilibrium(self):
        try:
            return self._eq
        except AttributeError:
            self._eq = self._equilibrium_at(self.T, self.P)
            return self._eq

    def _lnphis_from(self, eq) -> list[float]:
        n, species = eq
        out = [-math.log(n)] * self.N
        for i, (y, _) in species.items():
            out[i] = math.log(y / self.zs[i])
        return out

    def _H_dep_from(self, eq) -> float:
        _, species = eq
        dH = {i: dimerization_enthalpy(c) for i, c in self.acids}
        return sum(d * dH[i] for i, (_, d) in species.items())

    def _difference(self, f, variable: str):
        """Central difference of f(eq) in T or P."""
        T, P = self.T, self.P
        if variable == "T":
            h = 1e-5 * T
            return [(a - b) / (2 * h) for a, b in zip(
                _listed(f(self._equilibrium_at(T + h, P))),
                _listed(f(self._equilibrium_at(T - h, P))))]
        h = 1e-5 * P
        return [(a - b) / (2 * h) for a, b in zip(
            _listed(f(self._equilibrium_at(T, P + h))),
            _listed(f(self._equilibrium_at(T, P - h))))]

    # ---- fugacities ----------------------------------------------------------------------

    def lnphis(self):
        try:
            return self._lnphis
        except AttributeError:
            self._lnphis = self._lnphis_from(self._equilibrium())
            return self._lnphis

    lnphis_G_min = lnphis
    lnphis_lowest_Gibbs = lnphis

    def phis(self):
        return [math.exp(v) for v in self.lnphis()]

    def fugacities(self):
        P = self.P
        return [P * z * phi for z, phi in zip(self.zs, self.phis())]

    def dlnphis_dT(self):
        return self._difference(self._lnphis_from, "T")

    def dlnphis_dP(self):
        return self._difference(self._lnphis_from, "P")

    def dphis_dT(self):
        return [phi * d for phi, d in zip(self.phis(), self.dlnphis_dT())]

    def dphis_dP(self):
        return [phi * d for phi, d in zip(self.phis(), self.dlnphis_dP())]

    # ---- volume --------------------------------------------------------------------------

    def Z(self):
        return self._equilibrium()[0]

    def V(self):
        return self._equilibrium()[0] * self.R * self.T / self.P

    def dV_dT(self):
        dn = self._difference(lambda eq: eq[0], "T")[0]
        n = self._equilibrium()[0]
        return self.R * (n + self.T * dn) / self.P

    def dV_dP(self):
        dn = self._difference(lambda eq: eq[0], "P")[0]
        n = self._equilibrium()[0]
        return self.R * self.T * (dn / self.P - n / self.P**2)

    dV_dT_P = dV_dT
    dV_dP_T = dV_dP

    # ---- energies ------------------------------------------------------------------------

    def H_dep(self):
        return self._H_dep_from(self._equilibrium())

    def G_dep(self):
        return self.R * self.T * sum(z * ln for z, ln in zip(self.zs, self.lnphis()))

    def S_dep(self):
        return (self.H_dep() - self.G_dep()) / self.T

    def U_dep(self):
        return self.H_dep() - self.R * self.T * (self._equilibrium()[0] - 1.0)

    def A_dep(self):
        return self.U_dep() - self.T * self.S_dep()

    def H(self):
        try:
            return self._H
        except AttributeError:
            self._H = IdealGas.H(self) + self.H_dep()
            return self._H

    def S(self):
        try:
            return self._S
        except AttributeError:
            self._S = IdealGas.S(self) + self.S_dep()
            return self._S

    def Cp(self):
        try:
            return self._Cp
        except AttributeError:
            self._Cp = IdealGas.Cp(self) + self._difference(self._H_dep_from, "T")[0]
            return self._Cp

    dH_dT = Cp
    dH_dT_P = Cp

    def dS_dT(self):
        return self.Cp() / self.T

    dS_dT_P = dS_dT

    def dH_dP(self):
        return self._difference(self._H_dep_from, "P")[0]

    dH_dP_T = dH_dP

    def dS_dP(self):
        return -self.dV_dT()

    dS_dP_T = dS_dP


def _listed(value) -> list[float]:
    return list(value) if isinstance(value, (list, tuple)) else [value]


class FlashPureAssociating(FlashPureVLS):
    """FlashPureVLS for one dimerizing acid. Its shortcuts for an ideal gas over an activity
    liquid read the saturation pressure from the liquid, whose reference fugacity is the
    monomer pressure here; they take the measured vapor pressure instead."""

    def flash_TPV(self, T, P, V, zs=None, solution=None, hot_start=None, solution_target=None):
        if (self.VL_IG_activity and self.VL_IG_hack and V is None and solution is None
                and T is not None and P is not None):
            if P > self.correlations.VaporPressures[0](T):
                return None, [self.liquid.to(zs=zs, T=T, P=P)], [], [1.0], None
            return self.gas.to(zs=zs, T=T, P=P), [], [], [1.0], None
        return super().flash_TPV(T, P, V, zs=zs, solution=solution, hot_start=hot_start,
                                 solution_target=solution_target)

    def flash_TVF(self, T, VF=None, zs=None, hot_start=None):
        if self.VL_IG_activity:
            zs = [1.0]
            Psat = self.correlations.VaporPressures[0](T)
            return (Psat, self.liquid.to_TP_zs(T, Psat, zs), self.gas.to_TP_zs(T, Psat, zs),
                    0, 0.0)
        return super().flash_TVF(T, VF=VF, zs=zs, hot_start=hot_start)


def vapor_from_monomers(acids: dict[int, tuple[float, float]], T: float,
                        pressures: list[float]) -> tuple[float, list[float]]:
    """Total pressure and apparent vapor composition from the monomer partial pressures."""
    units, total = list(pressures), sum(pressures)
    for i, c in acids.items():
        dimer = math.exp(ln_k(c, T)) * pressures[i] ** 2
        units[i] += 2.0 * dimer
        total += dimer
    apparent = sum(units)
    return total, [u / apparent for u in units]
