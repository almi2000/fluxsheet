"""The liquid phase of the activity-coefficient methods, with measured heats of vaporization.

``thermo``'s ``GibbsExcessLiquid`` can take the liquid enthalpy from the vapor pressure curve
(``H_L = H_ig - R T^2 dlnPsat/dT``, Clausius-Clapeyron). That assumes the saturated vapor is an
ideal gas, so the heat of vaporization comes out too large by the factor 1 / (Z_V - Z_L):
about 5 % at 1 bar, 15 % at 5 bar and 25 to 40 % at 10 bar, and by 40 to 65 % for carboxylic
acids, whose vapor dimerizes. Condenser, reboiler and vaporizer duties inherit the error.

Here the liquid enthalpy is ``H_L = H_ig - Hvap(T)`` with the heat of vaporization
correlations of the compounds (``caloric_basis="Hvap"``), as process simulators do for
activity-coefficient methods. The phase equilibrium is untouched: fugacities still come from
the vapor pressures and activity coefficients.

``thermo`` evaluates H and S on this basis but not their temperature derivatives, and it
knows no heat of solution for Henry components. Both are added:

- a Henry component keeps ``R T^2 dlnH/dT`` (its heat of solution, as before),
- a compound without a correlation falls back to Clausius-Clapeyron, and so does one given
  :class:`ClausiusClapeyron` in place of its correlation (carboxylic acids, see
  ``PropertyPackage._warn_association``),
- above its critical temperature a compound's correlation gives zero.
"""

from __future__ import annotations

from thermo import GibbsExcessLiquid

R = 8.314462618  # J/(mol K)


class ClausiusClapeyron:
    """Stands in for a heat of vaporization correlation that is not to be used."""

    def __call__(self, T: float) -> None:
        return None

    def T_dependent_property_derivative(self, T: float) -> None:
        return None


class ActivityLiquid(GibbsExcessLiquid):
    """``GibbsExcessLiquid`` on the ``"Hvap"`` caloric basis, with consistent derivatives."""

    def _latent(self) -> tuple[list[float], list[float]]:
        """Heat of vaporization (or of solution) per component and its T derivative."""
        try:
            return self._Hvaps, self._dHvaps_dT
        except AttributeError:
            pass
        T = self.T
        henry = self.henry_components if self.has_henry_components else None
        values, slopes = [], []
        ratio = second = None
        for i, correlation in enumerate(self.EnthalpyVaporizations):
            value = None if henry is not None and henry[i] else correlation(T)
            if value is None:
                if ratio is None:
                    ratio = self.dPsats_dT_over_Psats()
                    second = self.d2Psats_dT2_over_Psats()
                a = ratio[i]
                values.append(R * T * T * a)
                slopes.append(R * T * (2.0 * a + T * (second[i] - a * a)))
            else:
                slope = correlation.T_dependent_property_derivative(T)
                values.append(value)
                slopes.append(slope if slope is not None else 0.0)
        self._Hvaps, self._dHvaps_dT = values, slopes
        return values, slopes

    def Hvaps(self):
        return self._latent()[0]

    def dHvaps_dT(self):
        return self._latent()[1]

    def Cp(self):
        try:
            return self._Cp
        except AttributeError:
            pass
        slopes = self._latent()[1]
        Cp = sum(z * (cp - d) for z, cp, d in zip(self.zs, self.Cpigs_pure(), slopes))
        if not self.composition_independent:
            Cp += self.GibbsExcessModel.CpE()
        self._Cp = Cp
        return Cp

    dH_dT = dH_dT_P = Cp  # thermo's pure-component flashes look up the _P names

    def dS_dT(self):
        try:
            return self._dS_dT
        except AttributeError:
            pass
        T = self.T
        values, slopes = self._latent()
        ratio = self.dPsats_dT_over_Psats()
        dS_dT = sum(z * (cp / T - R * a - d / T + h / (T * T))
                    for z, cp, a, d, h in zip(self.zs, self.Cpigs_pure(), ratio, slopes, values))
        if not self.composition_independent:
            dS_dT += self.GibbsExcessModel.dSE_dT()
        self._dS_dT = dS_dT
        return dS_dT

    dS_dT_P = dS_dT
