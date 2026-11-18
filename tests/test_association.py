"""Vapor-phase dimerization of carboxylic acids in the activity-coefficient methods.

Reference values: the compounds' heat of vaporization correlations, their formation
enthalpies as gas and liquid, the closed form of the dimerization equilibrium, and two
published facts about acid/water mixtures at 1 atm (acetic acid/water boils without an
azeotrope; formic acid/water has a maximum-boiling azeotrope at 107.3 °C and 77.5 wt %
formic acid, x_formic = 0.574).
"""

import math
from itertools import pairwise

import numpy as np
import pytest

import fluxsheet.thermo_package as tp
from fluxsheet.association import (
    DIMERIZATION,
    AssociatingGas,
    dimerization_enthalpy,
    ln_k,
    saturation_monomer,
)
from fluxsheet.fitting.models import Model, Pair
from fluxsheet.thermo_package import PropertyPackage, _estimate_binary, get_package

ACIDS = {"acetic acid": "64-19-7", "formic acid": "64-18-6", "propionic acid": "79-09-4"}
R = 8.314462618


@pytest.mark.parametrize("method", ["nrtl", "uniquac", "unifac"])
@pytest.mark.parametrize("name", list(ACIDS))
@pytest.mark.parametrize("P", [1.01325e5, 5e5, 10e5])
def test_acid_heat_of_vaporization_and_boiling_point_are_the_measured_ones(method, name, P):
    # An ideal monomer vapor gave 40 to 65 % too much.
    pkg = get_package((name,), method)
    liquid, vapor = pkg.flash([1.0], P=P, VF=0.0), pkg.flash([1.0], P=P, VF=1.0)
    assert liquid.T == pytest.approx(pkg.props.VaporPressures[0].solve_property(P), abs=1e-6)
    measured = pkg.props.EnthalpyVaporizations[0](liquid.T)
    assert vapor.H - liquid.H == pytest.approx(measured, rel=1e-6)


def test_saturated_acetic_acid_vapor_follows_the_dimerization_equilibrium():
    pkg = get_package(("acetic acid",), "nrtl")
    vapor = pkg.flash([1.0], P=101325.0, VF=1.0)
    constants = DIMERIZATION["64-19-7"]
    K = math.exp(ln_k(constants, vapor.T))
    p_m = 2.0 * 101325.0 / (1.0 + math.sqrt(1.0 + 4.0 * K * 101325.0))
    p_d = 101325.0 - p_m
    assert p_d == pytest.approx(K * p_m**2, rel=1e-12)
    n = (p_m + p_d) / (p_m + 2.0 * p_d)  # true molecules per acid unit
    assert vapor.V == pytest.approx(n * R * vapor.T / 101325.0, rel=1e-9)
    assert 2.0 * p_d / (p_m + 2.0 * p_d) == pytest.approx(0.745, abs=0.005)  # units in dimers
    density = pkg.mws[0] / 1000.0 / vapor.V
    assert density == pytest.approx(1.594 * 101325.0 * 0.0600520 / (R * vapor.T), rel=0.005)
    assert dimerization_enthalpy(constants) == pytest.approx(-60.6e3, rel=0.002)


def test_esterification_heat_of_reaction_follows_the_formation_data():
    # The liquid acid now lies below the monomer gas by the heat of vaporization plus the
    # dimers' bond enthalpy: 50.7 kJ/mol at 25 °C against 51.4 from the formation data.
    names = ("acetic acid", "ethanol", "ethyl acetate", "water")
    pkg = get_package(names, "nrtl")
    liquid = pkg._liquid_template().to(T=298.15, P=1e5, zs=[0.25] * 4)
    latent, c = liquid.Hvaps(), pkg.consts
    assert latent[0] == pytest.approx(c.Hfgs[0] - c.Hf_STPs[0], abs=1.5e3)
    nu = (-1, -1, 1, 1)
    reaction = sum(v * (c.Hfgs[i] - latent[i]) for i, v in enumerate(nu))
    assert reaction == pytest.approx(sum(v * c.Hf_STPs[i] for i, v in enumerate(nu)), abs=2e3)


@pytest.mark.parametrize("method", ["nrtl", "uniquac", "unifac"])
def test_acetic_acid_water_boils_without_an_azeotrope(method):
    # With an ideal monomer vapor every method predicted a minimum-boiling azeotrope below
    # 100 °C near pure water (372.3 K at x_water = 0.8).
    pkg = get_package(("water", "acetic acid"), method)
    xs = np.linspace(0.02, 0.98, 25)
    T = [pkg.saturation_point([x, 1.0 - x], "bubble", P=101325.0).T for x in xs]
    assert all(b < a for a, b in pairwise(T))
    assert 373.15 < T[-1] and T[0] < 391.1


@pytest.mark.parametrize("method", ["nrtl", "unifac"])
def test_formic_acid_water_maximum_boiling_azeotrope(method):
    # Measured: 107.3 °C (380.45 K) at x_formic = 0.574. Ideal monomer vapor: 377.5 K.
    pkg = get_package(("water", "formic acid"), method)
    xs = np.linspace(0.25, 0.65, 41)
    T = [pkg.saturation_point([x, 1.0 - x], "bubble", P=101325.0).T for x in xs]
    k = int(np.argmax(T))
    assert T[k] == pytest.approx(380.45, abs=2.5)
    assert 1.0 - xs[k] == pytest.approx(0.574, abs=0.06)
    assert 0 < k < len(xs) - 1  # a true maximum


@pytest.mark.parametrize("components,zs", [
    (("acetic acid",), [1.0]),
    (("water", "acetic acid"), [0.6, 0.4]),
    (("acetic acid", "formic acid", "water"), [0.3, 0.2, 0.5]),
])
@pytest.mark.parametrize("T,P", [(390.0, 1e5), (430.0, 3e5), (360.0, 2e4)])
def test_associating_vapor_is_thermodynamically_consistent(components, zs, T, P):
    gas = get_package(components, "nrtl")._flasher.gas
    assert isinstance(gas, AssociatingGas)

    def at(t=T, p=P):
        return gas.to(T=t, P=p, zs=zs)

    phase, h, dp = at(), 1e-3, 1e-3 * P
    assert phase.Cp() == pytest.approx((at(T + h).H() - at(T - h).H()) / (2 * h), rel=1e-6)
    assert phase.dS_dT() == pytest.approx((at(T + h).S() - at(T - h).S()) / (2 * h), rel=1e-6)
    # Gibbs-Helmholtz: the dimers' enthalpy is the temperature dependence of G_dep / RT.
    g = [at(t).G_dep() / (R * t) for t in (T - h, T + h)]
    assert phase.H_dep() == pytest.approx(-R * T * T * (g[1] - g[0]) / (2 * h), rel=1e-6)
    assert phase.G_dep() == pytest.approx(phase.H_dep() - T * phase.S_dep(), abs=1e-8)
    # Maxwell and the volume: dS/dP = -dV/dT, V = (true moles) R T / P.
    s = [at(p=P + s_).S() for s_ in (-dp, dp)]
    assert (s[1] - s[0]) / (2 * dp) == pytest.approx(-phase.dV_dT(), rel=1e-5)
    assert phase.V() == pytest.approx(phase.Z() * R * T / P, rel=1e-9)
    assert phase.dlnphis_dT() == pytest.approx(
        [(a - b) / (2 * h) for a, b in zip(at(T + h).lnphis(), at(T - h).lnphis())], rel=1e-5)


def test_acid_fugacity_is_its_monomer_pressure():
    gas = get_package(("water", "acetic acid"), "nrtl")._flasher.gas
    phase = gas.to(T=400.0, P=1e5, zs=[0.5, 0.5])
    n, species = phase._equilibrium()
    y_m, d = species[1]
    assert phase.fugacities()[1] == pytest.approx(y_m * 1e5, rel=1e-12)
    assert phase.fugacities()[0] == pytest.approx(0.5 / n * 1e5, rel=1e-12)
    K = math.exp(ln_k(DIMERIZATION["64-19-7"], 400.0))
    assert d / n == pytest.approx(K * 1e5 * y_m**2, rel=1e-12)  # y_D = K P y_m^2
    assert 0.5 == pytest.approx(n * y_m + 2.0 * d, rel=1e-12)


@pytest.mark.parametrize("T", [360.0, 391.0, 420.0])
def test_pure_acid_phase_and_round_trips(T):
    pkg = get_package(("acetic acid",), "nrtl")
    Psat = pkg.props.VaporPressures[0](T)
    p_m, _ = saturation_monomer(DIMERIZATION["64-19-7"], T, Psat)
    # The liquid's reference fugacity is the monomer pressure p_m, but the phase changes at
    # the vapor pressure: between the two the acid is vapor.
    assert p_m < 0.6 * Psat
    assert pkg.flash([1.0], T=T, P=0.5 * (p_m + Psat)).phase == "vapor"
    assert pkg.flash([1.0], T=T, P=0.99 * Psat).phase == "vapor"
    assert pkg.flash([1.0], T=T, P=1.01 * Psat).phase == "liquid"
    for P in (0.5 * Psat, 2.0 * Psat):
        state = pkg.flash([1.0], T=T, P=P)
        assert pkg.flash([1.0], P=P, H=state.H).T == pytest.approx(T, abs=1e-6)
        assert pkg.flash([1.0], P=P, S=state.S).T == pytest.approx(T, abs=1e-6)


def test_mixture_flashes_round_trip_in_the_two_phase_region():
    pkg = get_package(("water", "acetic acid"), "nrtl")
    state = pkg.flash([0.5, 0.5], T=380.0, P=101325.0)
    assert 0.0 < state.VF < 1.0
    for spec in ("H", "S"):
        back = pkg.flash([0.5, 0.5], P=101325.0, **{spec: getattr(state, spec)})
        assert back.T == pytest.approx(380.0, abs=1e-4)
        assert back.VF == pytest.approx(state.VF, abs=1e-5)


def test_binary_fit_model_uses_the_same_vapor():
    pair = Pair.from_cas("7732-18-5", "64-19-7")
    a12, a21, b12, b21, alpha = _estimate_binary("nrtl", "7732-18-5", "64-19-7")
    model = Model("nrtl", pair, alpha=alpha)
    xs = np.array([0.05, 0.5, 0.95])
    T, y = model.bubble_T((a12, a21, b12, b21), 101325.0, xs)
    pkg = get_package(("water", "acetic acid"), "nrtl")
    for x, t, y1 in zip(xs, T, y):
        point = pkg.saturation_point([x, 1.0 - x], "bubble", P=101325.0)
        assert t == pytest.approx(point.T, abs=0.01)
        assert y1 == pytest.approx(point.y[0], abs=1e-4)
    P, _ = model.bubble_P((a12, a21, b12, b21), T, xs)
    assert P == pytest.approx([101325.0] * 3, rel=1e-6)


def test_mixtures_without_acids_are_untouched():
    pkg = PropertyPackage(["water", "ethanol"], "nrtl")
    assert type(pkg._flasher.gas).__name__ == "IdealGas"
    assert pkg._liquid_template().VaporPressures == list(pkg.props.VaporPressures)
    assert not any("dimers" in w for w in pkg.warnings)


def test_acid_without_dimerization_constants_is_flagged(monkeypatch):
    monkeypatch.setattr(tp, "DIMERIZATION", {})
    pkg = PropertyPackage(["acetic acid", "water"], "nrtl")
    assert type(pkg._flasher.gas).__name__ == "IdealGas"
    flagged = [w for w in pkg.warnings if "dimers" in w]
    assert len(flagged) == 1 and "no dimerization constants" in flagged[0]


def test_described_acids_warn_only_with_an_equation_of_state():
    assert not any("dimers" in w for w in PropertyPackage(["acetic acid", "water"],
                                                           "uniquac").warnings)
    eos = [w for w in PropertyPackage(["acetic acid", "water"], "peng-robinson").warnings
           if "dimers" in w]
    assert len(eos) == 1 and "NRTL, UNIQUAC or UNIFAC" in eos[0]
