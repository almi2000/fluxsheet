"""Hand calculations and textbook results behind the process-engineering audit: latent heats,
liquid densities, the multistage compressor's discharge pressure and compressible pipe flow.

Reference values are either the compound correlations the property data come with (heats of
vaporization, liquid volumes) or closed-form results (Fanno flow of a perfect gas).
"""

import math

import numpy as np
import pytest
from scipy.optimize import brentq

from fluxsheet import solve
from fluxsheet.contents_flash import ContentsFlash
from fluxsheet.dynamics import simulate
from fluxsheet.thermo_package import PropertyPackage, get_package

R = 8.314462618


def one_unit(components, method, feed, unit, outlets=("out",)):
    model = {
        "components": list(components), "property_method": method,
        "units": [{"id": "F", "type": "Feed", "params": feed}, {"id": "U", **unit}],
        "streams": [{"id": "IN", "from": "F.out", "to": "U.in"},
                    *({"id": port.upper(), "from": f"U.{port}"} for port in outlets)],
    }
    return solve(model)


# --- heat of vaporization of the activity-coefficient methods -------------------------------

@pytest.mark.parametrize("method", ["nrtl", "ideal"])
@pytest.mark.parametrize("name", ["water", "benzene", "n-hexane", "ethanol"])
@pytest.mark.parametrize("P", [1.01325e5, 5e5, 10e5])
def test_latent_heat_is_the_measured_one_at_any_pressure(method, name, P):
    # Clausius-Clapeyron with an ideal vapor gave +4 % at 1 atm and +20 to 40 % at 10 bar.
    pkg = get_package((name,), method)
    liquid, vapor = pkg.flash([1.0], P=P, VF=0.0), pkg.flash([1.0], P=P, VF=1.0)
    measured = pkg.props.EnthalpyVaporizations[0](liquid.T)
    assert vapor.H - liquid.H == pytest.approx(measured, rel=1e-9)


@pytest.mark.parametrize("name,cp", [("water", 75.3), ("benzene", 135.7), ("n-hexane", 195.5),
                                      ("toluene", 157.1)])
def test_liquid_heat_capacity_at_25_degC(name, cp):
    pkg = get_package((name,), "nrtl")
    h = [pkg.phase_props("liquid", T, 101325.0, [1.0]).H for T in (297.15, 299.15)]
    assert (h[1] - h[0]) / 2.0 == pytest.approx(cp, rel=0.03)  # J/mol/K, literature


@pytest.mark.parametrize("components,method,zs", [
    (("water", "ethanol"), "nrtl", [0.6, 0.4]),
    (("benzene", "toluene"), "ideal", [0.5, 0.5]),
    (("water", "carbon dioxide", "nitrogen"), "nrtl", [0.98, 0.015, 0.005]),
    (("acetic acid", "water"), "uniquac", [0.5, 0.5]),
])
def test_liquid_enthalpy_and_entropy_derivatives_are_consistent(components, method, zs):
    liquid = get_package(components, method)._liquid_template()
    T, P, h = 330.0, 2e5, 1e-3

    def at(t):
        return liquid.to(T=t, P=P, zs=zs)

    phase = at(T)
    assert phase.Cp() == pytest.approx((at(T + h).H() - at(T - h).H()) / (2 * h), rel=1e-7)
    assert phase.dS_dT() == pytest.approx((at(T + h).S() - at(T - h).S()) / (2 * h), rel=1e-7)


def test_dissolved_gas_keeps_its_heat_of_solution():
    # CO2 in water at 25 °C: about -19.4 kJ/mol (from the temperature dependence of the
    # Henry constant); a heat of vaporization correlation would give zero above 304 K.
    liquid = get_package(("water", "carbon dioxide"), "nrtl")._liquid_template()
    cold = liquid.to(T=298.15, P=1e5, zs=[0.999, 0.001])
    hot = liquid.to(T=320.0, P=1e5, zs=[0.999, 0.001])
    assert -cold.Hvaps()[1] == pytest.approx(-19.4e3, rel=0.05)
    assert hot.Hvaps()[1] > 5e3  # above the critical temperature of CO2 as well


def test_isentropic_flash_into_the_two_phase_region_with_an_activity_method():
    pkg = get_package(("water", "ethanol"), "nrtl")
    state = pkg.flash([0.6, 0.4], T=355.0, P=1.01325e5)
    assert 0.0 < state.VF < 1.0
    back = pkg.flash([0.6, 0.4], P=1.01325e5, S=state.S)
    assert back.T == pytest.approx(state.T, abs=1e-4)
    assert back.VF == pytest.approx(state.VF, abs=1e-5)


def test_carboxylic_acids_are_flagged_with_an_equation_of_state():
    # The activity methods describe the dimers (test_association.py); a cubic equation
    # of state cannot.
    pkg = PropertyPackage(["acetic acid", "water"], "peng-robinson")
    flagged = [w for w in pkg.warnings if "acetic acid" in w and "dimers" in w]
    assert len(flagged) == 1 and "without association" in flagged[0]
    for method in ("peng-robinson", "nrtl"):
        assert not any("dimers" in w for w in PropertyPackage(["ethanol", "water"],
                                                              method).warnings)


# --- liquid volumes of the equations of state -----------------------------------------------

@pytest.mark.parametrize("method", ["peng-robinson", "srk"])
@pytest.mark.parametrize("name,rho", [("water", 997.0), ("methanol", 787.0), ("acetone", 784.0),
                                       ("benzene", 874.0), ("n-hexane", 655.0)])
def test_equation_of_state_liquid_density_at_25_degC(method, name, rho):
    # Untranslated: Peng-Robinson -15 % for water, SRK -9 to -24 % for all of them.
    pkg = get_package((name,), method)
    volume = pkg.phase_props("liquid", 298.15, 101325.0, [1.0]).V
    assert pkg.mws[0] / 1000.0 / volume == pytest.approx(rho, rel=0.005)


@pytest.mark.parametrize("method", ["peng-robinson", "srk"])
def test_volume_translation_leaves_the_phase_equilibrium_alone(method, monkeypatch):
    components, zs, P = ["methane", "propane", "n-hexane"], [0.3, 0.4, 0.3], 20e5
    translated = PropertyPackage(components, method)
    monkeypatch.setattr(PropertyPackage, "_volume_translations",
                        lambda self, eos: [0.0] * self.n)
    plain = PropertyPackage(components, method)
    for vf in (0.0, 1.0):
        assert translated.flash(zs, P=P, VF=vf).T == pytest.approx(
            plain.flash(zs, P=P, VF=vf).T, abs=1e-6)
    a, b = translated.flash(zs, P=P, T=300.0), plain.flash(zs, P=P, T=300.0)
    assert a.VF == pytest.approx(b.VF, abs=1e-8)
    assert a.y == pytest.approx(b.y, abs=1e-8)
    assert a.H == pytest.approx(b.H, abs=1e-9)  # the equation itself is untouched
    liquid = translated.flash(zs, P=P, T=190.0).V
    assert liquid != pytest.approx(plain.flash(zs, P=P, T=190.0).V, rel=1e-3)


def test_pump_power_of_water_no_longer_depends_on_the_method():
    feed = {"T": 298.15, "P": "1 atm", "flows": {"water": 1000.0 / 18.01528}}
    power = [one_unit(["water"], method, feed, {"type": "Pump", "params": {
        "dP": "10 bar", "efficiency": 0.75}}).units["U"]["power"]
             for method in ("nrtl", "peng-robinson", "srk")]
    assert power[0] == pytest.approx(1e6 / 997.0 / 0.75, rel=0.002)  # V dP / efficiency
    assert power[1:] == pytest.approx([power[0]] * 2, rel=0.005)


# --- multistage compressor ------------------------------------------------------------------

@pytest.mark.parametrize("aftercooler", ["no", "yes"])
@pytest.mark.parametrize("stages", [2, 4])
def test_multistage_compressor_reaches_its_outlet_pressure_with_cooler_drops(aftercooler,
                                                                            stages):
    r = one_unit(["nitrogen"], "peng-robinson",
                 {"T": 300, "P": "1 bar", "flows": {"nitrogen": 10}},
                 {"type": "MultistageCompressor", "params": {
                     "P": "9 bar", "stages": stages, "T_cooler": 310, "dP_cooler": "0.2 bar",
                     "aftercooler": aftercooler}}, outlets=("out", "condensate"))
    assert r.ok, r.error
    assert r.streams["OUT"]["P"] == pytest.approx(9e5, rel=1e-10)
    ratio, P = r.units["U"]["pressure_ratio"], 1e5
    for k in range(1, stages + 1):  # every stage has the same pressure ratio
        P = P * ratio - (0.2e5 if k < stages or aftercooler == "yes" else 0.0)
    assert P == pytest.approx(9e5, rel=1e-10)
    assert ratio > 9.0 ** (1.0 / stages)


# --- compressible pipe flow -----------------------------------------------------------------

def gas_pipe(flow, length=30.0, **params):
    return one_unit(["nitrogen"], "peng-robinson",
                    {"T": 300, "P": "5 bar", "flows": {"nitrogen": flow}},
                    {"type": "Pipe", "params": {"length": length, "diameter": 0.05,
                                                "segments": 200, **params}})


def fanno(M, gamma=1.4):
    """4 f L* / D of adiabatic pipe flow of a perfect gas (Darcy f: f L* / D)."""
    return (1 - M * M) / (gamma * M * M) + (gamma + 1) / (2 * gamma) * math.log(
        (gamma + 1) * M * M / (2 + (gamma - 1) * M * M))


def test_gas_pipe_follows_fanno_flow():
    r = gas_pipe(20)
    assert r.ok, r.error
    u = r.units["U"]
    gamma = 1.4
    M1 = u["velocity_in"] / math.sqrt(gamma * R * 300.0 / 0.0280134)
    remaining = fanno(M1) - u["friction_factor"] * 30.0 / 0.05
    M2 = brentq(lambda M: fanno(M) - remaining, M1, 1.0 - 1e-12)
    P2 = 5e5 * (M1 / M2) * math.sqrt((1 + 0.2 * M1 * M1) / (1 + 0.2 * M2 * M2))
    assert u["mach"] == pytest.approx(M2, rel=0.01)
    assert r.streams["OUT"]["P"] == pytest.approx(P2, rel=0.005)
    # The stream carries the total enthalpy: the flowsheet balance closes without a
    # kinetic-energy term, and an adiabatic ideal gas leaves at its inlet temperature.
    assert r.streams["OUT"]["enthalpy_flow"] == pytest.approx(
        r.streams["IN"]["enthalpy_flow"], abs=1e-6)
    assert all(c["status"] == "ok" for c in r.balance_checks)


def test_gas_pipe_chokes_where_fanno_flow_does():
    r = gas_pipe(40)
    assert not r.ok
    assert "chokes" in r.error["message"] and "speed of sound" in r.error["message"]
    # At 40 mol/s the inlet Mach number is 0.286 and f L* / D = 6.0: sonic after 15.8 m.
    assert "after 15." in r.error["message"]
    shorter = gas_pipe(40, length=14.0)
    assert shorter.ok, shorter.error
    assert 0.5 < shorter.units["U"]["mach"] < 1.0


def test_acceleration_is_nothing_for_a_liquid_pipe():
    r = one_unit(["water"], "nrtl",
                 {"T": 298.15, "P": "5 bar", "flows": {"water": 1000.0 / 18.01528}},
                 {"type": "Pipe", "params": {"length": 10, "diameter": 0.025}})
    u = r.units["U"]
    rho, mu = 997.0, 0.00089
    v = 1.0 / (rho * math.pi * 0.025**2 / 4)
    f = 0.25 / math.log10(0.046e-3 / 0.025 / 3.7 + 5.74 / (rho * v * 0.025 / mu) ** 0.9) ** 2
    assert u["dP"] == pytest.approx(f * 10 / 0.025 * rho * v * v / 2, rel=0.005)
    assert u["mach"] is None  # the liquid model is incompressible
    assert r.streams["OUT"]["T"] == pytest.approx(298.15, abs=1e-6)


# --- vessel flashes along a trajectory ------------------------------------------------------

def test_vessel_flash_does_not_depend_on_the_calls_before():
    pkg = get_package(("water", "ethanol"), "nrtl")
    zs, P = np.array([0.9, 0.1]), 1.2e5
    start = pkg.flash(list(zs), T=298.15, P=P)
    results = []
    for detour in ((), (5.0, -8.0, 3.0), (-20.0, 20.0)):
        flash = ContentsFlash(pkg, start.T)
        flash(zs, P, start.H)
        for dh in detour:
            flash(zs, P, start.H + dh)
        results.append(flash(zs, P, start.H + 1.0).T)
    assert results[0] == results[1] == results[2]  # to the last bit


def test_feed_step_into_an_uncontrolled_tank_integrates():
    # At an exact steady state the stiff integrator stalled on bit-level noise of the
    # vessel flash before the event was even reached.
    model = {
        "components": ["water", "ethanol"], "property_method": "nrtl",
        "units": [
            {"id": "FEED", "type": "Feed", "params": {
                "T": "25 degC", "P": "1.5 bar",
                "flows": {"water": "50 kmol/h", "ethanol": "5 kmol/h"}}},
            {"id": "TANK", "type": "MixingTank", "params": {"P": "1.2 bar"}},
            {"id": "PROD", "type": "Product"},
        ],
        "streams": [{"id": "S1", "from": "FEED.out", "to": "TANK.in"},
                    {"id": "S2", "from": "TANK.liquid", "to": "PROD.in"}],
        "dynamics": {"t_end": "3 h", "points": 13,
                     "units": {"TANK": {"diameter": "1 m", "height": "2 m"}},
                     "events": [{"time": "1 min", "unit": "FEED", "param": "flows",
                                 "component": "water", "value": "60 kmol/h"}]},
    }
    out = simulate(model)
    assert out["ok"], out.get("error")
    level = next(s["values"] for s in out["series"] if s["key"] == "unit.TANK.level")
    assert level[-1] > level[0]  # 10 kmol/h more in, the same out
