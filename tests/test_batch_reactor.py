"""Batch reactor: one batch cycle checked against closed-form results and independent
balances."""

import math
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from fluxsheet import FlowsheetModel, solve
from fluxsheet.flowsheet import ReactionModel, UtilityModel
from fluxsheet.reactions import R, parse_reactions
from fluxsheet.stream import Stream
from fluxsheet.thermo_package import get_package
from fluxsheet.units.base import empty_stream
from fluxsheet.units.batch import BatchReactor
from fluxsheet.utilities import parse_utilities

EXAMPLES = Path(__file__).parent.parent / "examples"

ISOMER = ["n-butane", "isobutane"]
K1 = 0.05  # 1/s
FIRST_ORDER = {"id": "ISO", "stoichiometry": {"n-butane": -1, "isobutane": 1},
               "kinetics": {"k0": K1, "Ea": "0 kJ/mol"}}
SYNGAS = ["carbon monoxide", "water", "carbon dioxide", "hydrogen"]
STEAM = {"id": "LPS", "type": "steam", "P": "3 bar"}
CW = {"id": "CW", "type": "fluid", "T_in": "15 degC", "T_out": "25 degC"}


def run(params, *, components=ISOMER, flows=None, T="400 K", P="1 bar", reactions=None,
        utilities=(), ok=True):
    r = solve({
        "components": list(components),
        "property_method": "ideal",
        "reactions": list(reactions if reactions is not None else [FIRST_ORDER]),
        "utilities": list(utilities),
        "units": [
            {"id": "F", "type": "Feed", "params": {"T": T, "P": P,
                                                   "flows": flows or {"n-butane": "1 mol/s"}}},
            {"id": "B", "type": "BatchReactor", "params": {"reactions": ["ISO"], **params}},
        ],
        "streams": [{"id": "IN", "from": "F.out", "to": "B.in"},
                    {"id": "OUT", "from": "B.out"}, {"id": "VENT", "from": "B.vent"}],
    })
    if ok:
        assert r.ok, r.error
    return r


def series(unit, profile, label):
    prof = unit["profiles"][profile]
    return np.array(prof["x"]), np.array(next(s["values"] for s in prof["series"]
                                              if s["label"] == label))


def h_flow(components, stream):
    """Enthalpy flow with formation enthalpies, W."""
    pkg = get_package(tuple(components), "ideal")
    flows = [stream["component_flows"][n] for n in pkg.names]
    return stream["enthalpy_flow"] + sum(f * h for f, h in zip(flows, pkg.consts.Hfgs))


def mass(stream):
    return stream["mass_flow"]


# --- kinetics against closed forms ----------------------------------------------------------

@pytest.mark.parametrize("t", [5.0, 20.0, 60.0])
def test_first_order_isothermal_follows_exp(t):
    # Ideal gas at constant T and P: dN_A/dt = -k c_A V = -k N_A, so X = 1 - exp(-k t).
    r = run({"phase": "gas", "batch_time": t, "T": "400 K"})
    u = r.units["B"]
    assert u["conversion"] == pytest.approx(1.0 - math.exp(-K1 * t), rel=1e-6)
    x, X = series(u, "conditions", "Conversion of n-butane")
    assert X == pytest.approx(1.0 - np.exp(-K1 * x * 60.0), abs=1e-6)  # minutes on the axis


def test_first_order_in_the_liquid():
    r = run({"phase": "liquid", "batch_time": "30 s", "T": "280 K"}, T="280 K", P="5 bar")
    assert r.units["B"]["conversion"] == pytest.approx(1.0 - math.exp(-K1 * 30.0), rel=1e-6)


def test_second_order_against_closed_form():
    # CO + H2O -> CO2 + H2, equimolar ideal gas at constant T, P: moles and volume stay
    # constant, c_CO = c_H2O = c0 (1 - X), so X = k c0 t / (1 + k c0 t).
    k, T, P, t = 2e-3, 600.0, 2e5, 300.0
    shift = {"id": "ISO", "stoichiometry": {"carbon monoxide": -1, "water": -1,
                                            "carbon dioxide": 1, "hydrogen": 1},
             "kinetics": {"k0": k, "Ea": "0 kJ/mol"}}
    r = run({"phase": "gas", "batch_time": t, "T": T}, components=SYNGAS,
            flows={"carbon monoxide": "1 mol/s", "water": "1 mol/s"}, T=T, P=P,
            reactions=[shift])
    c0 = 0.5 * P / (R * T)
    assert r.units["B"]["conversion"] == pytest.approx(k * c0 * t / (1 + k * c0 * t), rel=1e-6)


# --- energy -----------------------------------------------------------------------------------

def test_adiabatic_rise_closes_the_energy_balance():
    r = run({"phase": "gas", "batch_time": "60 s"}, T="400 K")
    u = r.units["B"]
    assert u["T_out"] > 400.0 + 5.0  # isomerization to isobutane is exothermic
    assert u["heat_per_batch"] == pytest.approx(0.0, abs=1e-6 * abs(u["charge"]) * 1e3)
    h_in, h_out = h_flow(ISOMER, r.streams["IN"]), h_flow(ISOMER, r.streams["OUT"])
    assert h_out == pytest.approx(h_in, rel=1e-7, abs=1e-3)
    # The final temperature is where the products have the feed's enthalpy.
    pkg = get_package(tuple(ISOMER), "ideal")
    out = r.streams["OUT"]
    zs = [out["mole_fractions"][n] for n in pkg.names]
    h_target = (h_in / out["molar_flow"]) - sum(z * h for z, h in zip(zs, pkg.consts.Hfgs))
    assert pkg.flash(zs, P=1e5, H=h_target).T == pytest.approx(u["T_out"], abs=1e-4)


def test_isothermal_duty_is_the_enthalpy_change():
    r = run({"phase": "gas", "batch_time": "60 s", "T": "420 K"}, T="400 K")
    u = r.units["B"]
    dh = h_flow(ISOMER, r.streams["OUT"]) - h_flow(ISOMER, r.streams["IN"])
    assert u["duty"] == pytest.approx(dh, rel=1e-7)  # includes heating the charge to 420 K
    assert u["peak_duty"] < 0.0  # the reaction gives off heat at 420 K


def test_steam_jacket_heats_like_an_independent_integration():
    # No reaction (k0 = 0): the steam jacket heats the gas, N dh/dt = UA (Ts - T).
    still = {**FIRST_ORDER, "kinetics": {"k0": 0.0, "Ea": "0 kJ/mol"}}
    UA, t_end, F = 50.0, 600.0, 1.0
    r = run({"phase": "gas", "batch_time": t_end, "area": "1 m2", "U": UA, "utility": "LPS",
             "cycle_time": t_end}, T="350 K", reactions=[still], utilities=[STEAM])
    u = r.units["B"]
    pkg = get_package(tuple(ISOMER), "ideal")
    steam = next(iter(parse_utilities([UtilityModel(**STEAM)]).values()))
    Ts = steam.T_in
    N = F * t_end

    def dT(t, T):
        h = pkg.phase_props("vapor", T[0], 1e5, [1.0, 0.0]).H
        cp = (pkg.phase_props("vapor", T[0] + 1e-3, 1e5, [1.0, 0.0]).H - h) / 1e-3
        return [UA * (Ts - T[0]) / (N * cp)]

    ref = solve_ivp(dT, (0, t_end), [350.0], rtol=1e-10, atol=1e-10).y[0, -1]
    assert u["T_out"] == pytest.approx(ref, abs=1e-3)
    dh = h_flow(ISOMER, r.streams["OUT"]) - h_flow(ISOMER, r.streams["IN"])
    assert u["duty"] == pytest.approx(dh, rel=1e-6)
    assert u["peak_duty"] == pytest.approx(UA * (Ts - 350.0), rel=1e-3)
    assert u["utility_flow"] == pytest.approx(u["duty"] / steam.latent, rel=1e-6)


def test_scale_invariance_with_a_cooling_jacket():
    # Twice the feed with twice the area and medium flow: the same batch, twice the flows.
    base = {"phase": "gas", "batch_time": "60 s", "U": 200.0, "utility": "CW",
            "utility_flow": "0.5 kg/s", "cycle_time": "60 s"}
    one = run({**base, "area": "1 m2"}, utilities=[CW])
    two = run({**base, "area": "2 m2", "utility_flow": "1 kg/s"},
              flows={"n-butane": "2 mol/s"}, utilities=[CW])
    a, b = one.units["B"], two.units["B"]
    adiabatic = run({"phase": "gas", "batch_time": "60 s"}).units["B"]["T_out"]
    assert a["T_out"] < adiabatic - 5.0  # cooled
    assert b["T_out"] == pytest.approx(a["T_out"], abs=1e-5)
    assert b["conversion"] == pytest.approx(a["conversion"], rel=1e-7)
    assert b["duty"] == pytest.approx(2 * a["duty"], rel=1e-6)
    assert two.streams["OUT"]["molar_flow"] == pytest.approx(
        2 * one.streams["OUT"]["molar_flow"], rel=1e-9)


def test_temperature_program_hits_its_breakpoints():
    prog = [{"time": "0 s", "T": "380 K"}, {"time": "20 s", "T": "420 K"},
            {"time": "40 s", "T": "420 K"}]
    r = run({"phase": "gas", "batch_time": "50 s", "T_profile": prog})
    u = r.units["B"]
    x, T = series(u, "conditions", "Temperature")
    assert T[0] == pytest.approx(380.0)
    assert np.interp(10.0 / 60.0, x, T) == pytest.approx(400.0, abs=1e-6)
    assert u["T_out"] == pytest.approx(420.0)
    dh = h_flow(ISOMER, r.streams["OUT"]) - h_flow(ISOMER, r.streams["IN"])
    assert u["duty"] == pytest.approx(dh, rel=1e-7)


def test_constant_volume_pressure_rises_with_the_moles():
    # Ethane -> ethylene + hydrogen in a closed vessel at constant T: P/P0 = 1 + X for an
    # ideal gas, and the heat is the internal energy of reaction, dH - dn R T per extent.
    comps = ["ethane", "ethylene", "hydrogen"]
    crack = {"id": "ISO", "stoichiometry": {"ethane": -1, "ethylene": 1, "hydrogen": 1},
             "kinetics": {"k0": 0.01, "Ea": "0 kJ/mol"}}
    T = 900.0
    r = run({"phase": "gas", "batch_time": "100 s", "T": T, "vessel": "constant_volume",
             "cycle_time": "100 s"}, components=comps, flows={"ethane": "1 mol/s"}, T=T,
            reactions=[crack])
    u = r.units["B"]
    X = u["conversion"]
    assert X == pytest.approx(1.0 - math.exp(-1.0), rel=1e-6)
    assert u["P_final"] == pytest.approx(1e5 * (1.0 + X), rel=1e-7)
    pkg = get_package(tuple(comps), "ideal")
    rx = parse_reactions([ReactionModel(**crack)], pkg).by_id["ISO"]
    xi = X * u["charge"]
    assert u["heat_per_batch"] == pytest.approx(xi * (rx.heat_of_reaction(pkg, T) - R * T),
                                                rel=1e-6)


# --- stop criteria ----------------------------------------------------------------------------

def test_stop_at_conversion_hits_it():
    r = run({"phase": "gas", "batch_time": "10 min", "stop_conversion": 0.8, "T": "400 K"})
    u = r.units["B"]
    assert u["conversion"] == pytest.approx(0.8, abs=1e-8)
    assert u["batch_time"] == pytest.approx(-math.log(0.2) / K1, rel=1e-6)


def test_stop_at_a_falling_mole_fraction():
    r = run({"phase": "gas", "batch_time": "10 min", "stop_mole_fraction": 0.5,
             "stop_component": "n-butane", "T": "400 K"})
    assert r.units["B"]["batch_time"] == pytest.approx(math.log(2.0) / K1, rel=1e-6)


def test_stop_at_temperature_in_an_adiabatic_batch():
    r = run({"phase": "gas", "batch_time": "10 min", "stop_T": "410 K"})
    u = r.units["B"]
    assert u["T_out"] == pytest.approx(410.0, abs=1e-6)
    assert 0.0 < u["conversion"] < 0.9


# --- vent -------------------------------------------------------------------------------------

def vented():
    slow = {**FIRST_ORDER, "kinetics": {"k0": 1e-3, "Ea": "0 kJ/mol"}}
    return run({"phase": "liquid", "batch_time": "10 min", "duty": "5 kW", "vent": "yes",
                "cycle_time": "10 min"}, T="300 K", P="4 bar", reactions=[slow])


@pytest.mark.parametrize("params", [
    {"phase": "gas", "batch_time": "60 s"},
    {"phase": "gas", "batch_time": "60 s", "vessel": "constant_volume"},
])
def test_batch_reactor_balance_diagnostics(params):
    result = run(params, T="400 K")
    for scope in ("B", "flowsheet"):
        checks = [check for check in result.balance_checks if check["scope"] == scope]
        assert any(check["quantity"] == "element" and check["component"] == "C"
                   for check in checks)
        assert next(check for check in checks if check["quantity"] == "energy")["status"] == "ok"
        assert not any(check["status"] == "warning" for check in checks)


def test_vented_batch_balance_diagnostics():
    result = vented()
    for scope in ("B", "flowsheet"):
        energy = next(check for check in result.balance_checks
                      if check["scope"] == scope and check["quantity"] == "energy")
        assert energy["status"] == "ok"


def test_constant_volume_batch_balance_includes_pressure_volume_change():
    cracking = {"id": "ISO", "stoichiometry": {"ethane": -1, "ethylene": 1, "hydrogen": 1},
                "kinetics": {"k0": 0.01, "Ea": "0 kJ/mol"}}
    result = run({"phase": "gas", "batch_time": "100 s", "T": "900 K",
                  "vessel": "constant_volume", "cycle_time": "100 s"},
                 components=["ethane", "ethylene", "hydrogen"],
                 flows={"ethane": "1 mol/s"}, T="900 K", reactions=[cracking])
    for scope in ("B", "flowsheet"):
        energy = next(check for check in result.balance_checks
                      if check["scope"] == scope and check["quantity"] == "energy")
        assert energy["status"] == "ok"
    assert result.units["B"]["P_final"] > result.streams["IN"]["P"] * 1.1


def test_constant_volume_batch_with_different_fill_pressure_is_not_checkable(
        expected_balance_exceptions):
    for scope in ("B", "flowsheet"):
        expected_balance_exceptions[(scope, "energy", None, "not_checkable")] = (
            "Filling work is unavailable when vessel and feed pressures differ")
    result = run({"phase": "gas", "batch_time": "20 s", "vessel": "constant_volume",
                  "P": "2 bar"}, T="400 K", P="1 bar")
    for scope in ("B", "flowsheet"):
        energy = next(check for check in result.balance_checks
                      if check["scope"] == scope and check["quantity"] == "energy")
        assert energy["status"] == "not_checkable"


def test_zero_flow_batch_has_valid_zero_energy_balance():
    model = {
        "components": ISOMER,
        "property_method": "ideal",
        "reactions": [FIRST_ORDER],
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "T": "400 K", "P": "1 bar", "flows": {"n-butane": "1 mol/s"}}},
            {"id": "S", "type": "Splitter", "params": {"split": 0}},
            {"id": "B", "type": "BatchReactor", "params": {
                "reactions": ["ISO"], "phase": "gas", "batch_time": "20 s",
                "vessel": "constant_volume"}},
        ],
        "streams": [
            {"id": "feed", "from": "F.out", "to": "S.in"},
            {"id": "empty", "from": "S.out1", "to": "B.in"},
            {"id": "bypass", "from": "S.out2"},
            {"id": "out", "from": "B.out"},
        ],
    }
    result = solve(model)
    assert result.ok, result.error
    for scope in ("B", "flowsheet"):
        energy = next(check for check in result.balance_checks
                      if check["scope"] == scope and check["quantity"] == "energy")
        assert energy["status"] == "ok"


def test_vent_closes_the_mass_and_energy_balance():
    r = vented()
    u = r.units["B"]
    feed, out, vent = r.streams["IN"], r.streams["OUT"], r.streams["VENT"]
    assert vent["molar_flow"] > 0.1 * feed["molar_flow"]
    assert mass(out) + mass(vent) == pytest.approx(mass(feed), rel=1e-7)
    assert u["vent_amount"] == pytest.approx(vent["molar_flow"] * 600.0, rel=1e-9)
    dh = h_flow(ISOMER, out) + h_flow(ISOMER, vent) - h_flow(ISOMER, feed)
    assert dh == pytest.approx(5e3 * 600.0 / u["cycle_time"], rel=1e-4)
    # The boiling liquid stays at its bubble point at the vessel pressure.
    pkg = get_package(tuple(ISOMER), "ideal")
    zs = [out["mole_fractions"][n] for n in pkg.names]
    assert pkg.flash(zs, P=4e5, VF=0.0).T == pytest.approx(u["T_out"], abs=0.05)


def test_vent_follows_rayleigh_distillation():
    # No reaction: heating boils the subcooled binary off at its bubble point, so the liquid
    # follows Rayleigh's equation, ln(N/N0) = integral from x0 to x of dx / (y(x) - x).
    still = {**FIRST_ORDER, "kinetics": {"k0": 0.0, "Ea": "0 kJ/mol"}}
    r = run({"phase": "liquid", "batch_time": "10 min", "duty": "10 kW", "vent": "yes",
             "cycle_time": "10 min"}, flows={"n-butane": "0.5 mol/s", "isobutane": "0.5 mol/s"},
            T="300 K", P="4 bar", reactions=[still])
    out = r.streams["OUT"]
    x_end, N_end = out["mole_fractions"]["n-butane"], out["molar_flow"] * 600.0
    pkg = get_package(tuple(ISOMER), "ideal")

    def dlnN(x):
        y = pkg.saturation_point([x, 1.0 - x], "bubble", P=4e5).y[0]
        return 1.0 / (y - x)

    xs = np.linspace(0.5, x_end, 41)
    ln_ratio = np.trapezoid([dlnN(x) for x in xs], xs)
    assert x_end > 0.52  # the lighter isobutane leaves first
    assert N_end / 600.0 == pytest.approx(math.exp(ln_ratio), rel=1e-3)
    vent = r.streams["VENT"]
    assert mass(out) + mass(vent) == pytest.approx(mass(r.streams["IN"]), rel=1e-9)


def energy_gap(r):
    """Enthalpy out (product + vent) minus in, minus the reported average duty, W."""
    out = h_flow(ISOMER, r.streams["OUT"]) + h_flow(ISOMER, r.streams["VENT"])
    return out - h_flow(ISOMER, r.streams["IN"]) - r.units["B"]["duty"]


@pytest.mark.parametrize("duty", ["0 kW", "-5 kW"])
def test_vented_charge_above_its_bubble_point_flashes_off(duty):
    # Feed at 10 bar and 325 K let down into a 4 bar vessel: part of it flashes to the vent
    # at once. Without heat the rest stays at its bubble point; cooled, it ends subcooled.
    still = {**FIRST_ORDER, "kinetics": {"k0": 0.0, "Ea": "0 kJ/mol"}}
    r = run({"phase": "liquid", "batch_time": "10 min", "duty": duty, "vent": "yes",
             "cycle_time": "10 min", "P": "4 bar"},
            flows={"n-butane": "0.5 mol/s", "isobutane": "0.5 mol/s"}, T="325 K", P="10 bar",
            reactions=[still])
    u = r.units["B"]
    assert u["vent_amount"] > 0.05 * u["charge"]
    assert energy_gap(r) == pytest.approx(0.0, abs=1e-6 * abs(h_flow(ISOMER, r.streams["IN"])))
    pkg = get_package(tuple(ISOMER), "ideal")
    out = r.streams["OUT"]
    zs = [out["mole_fractions"][n] for n in pkg.names]
    T_bubble = pkg.saturation_point(zs, "bubble", P=4e5).T
    if duty == "0 kW":
        assert u["T_out"] == pytest.approx(T_bubble, abs=0.02)
    else:
        assert u["T_out"] < T_bubble - 1.0


def test_superheated_pure_liquid_charge_boils_down_to_its_boiling_point():
    still = {**FIRST_ORDER, "kinetics": {"k0": 0.0, "Ea": "0 kJ/mol"}}
    r = run({"phase": "liquid", "batch_time": "5 min", "vent": "yes", "P": "4 bar"},
            T="325 K", P="10 bar", reactions=[still])
    u = r.units["B"]
    pkg = get_package(tuple(ISOMER), "ideal")
    assert u["T_out"] == pytest.approx(pkg.flash([1.0, 0.0], P=4e5, VF=0.0).T, abs=0.02)
    assert u["vent_amount"] > 0.0
    assert energy_gap(r) == pytest.approx(0.0, abs=1e-6 * abs(h_flow(ISOMER, r.streams["IN"])))


def test_constant_volume_adiabatic_conserves_internal_energy():
    # Endothermic cracking in a closed, insulated vessel: U stays, so
    # H_out - P_end V = H_in - P0 V per batch (ideal gas: P V = N R T).
    comps = ["ethane", "ethylene", "hydrogen"]
    crack = {"id": "ISO", "stoichiometry": {"ethane": -1, "ethylene": 1, "hydrogen": 1},
             "kinetics": {"k0": 0.05, "Ea": "150 kJ/mol", "T_ref": "1000 K"}}
    T0, t_c = 1000.0, 100.0
    r = run({"phase": "gas", "batch_time": "60 s", "vessel": "constant_volume",
             "cycle_time": t_c}, components=comps, flows={"ethane": "1 mol/s"}, T=T0,
            reactions=[crack])
    u = r.units["B"]
    assert u["T_out"] < T0 - 20.0
    out = r.streams["OUT"]
    n_end = out["molar_flow"] * t_c
    pv_end = n_end * R * u["T_out"]
    pv_0 = u["charge"] * R * T0
    assert u["P_final"] == pytest.approx(1e5 * n_end * u["T_out"] / (u["charge"] * T0),
                                         rel=1e-7)
    assert h_flow(comps, out) * t_c - pv_end == pytest.approx(
        h_flow(comps, r.streams["IN"]) * t_c - pv_0, rel=1e-8, abs=1e-3)


def test_cooling_water_flow_follows_the_heat_removed():
    r = run({"phase": "gas", "batch_time": "60 s", "U": 200.0, "area": "1 m2",
             "utility": "CW", "cycle_time": "100 s"}, T="420 K", utilities=[CW])
    u = r.units["B"]
    assert u["duty"] < 0.0
    assert u["utility_flow"] == pytest.approx(-u["heat_per_batch"] / (4180.0 * 10.0) / 100.0,
                                              rel=1e-9)


def test_stop_at_a_rising_mole_fraction():
    r = run({"phase": "gas", "batch_time": "10 min", "stop_mole_fraction": 0.75,
             "stop_component": "isobutane", "T": "400 K"})
    assert r.units["B"]["batch_time"] == pytest.approx(math.log(4.0) / K1, rel=1e-6)


def test_changed_parameters_are_recalculated():
    pkg = get_package(tuple(ISOMER), "ideal")
    unit = BatchReactor("B", {"reactions": ["ISO"], "phase": "gas", "batch_time": 20.0,
                              "T": 400.0})
    unit.reactions = parse_reactions([ReactionModel(**FIRST_ORDER)], pkg)
    feed = {"in": [Stream.from_flash(pkg, [1.0, 0.0], T=400.0, P=1e5)]}
    unit.calculate(pkg, feed)
    first = unit.results["conversion"]
    unit.p["batch_time"] = 40.0  # as a dynamic event or controller would
    unit.calculate(pkg, feed)
    assert unit.results["conversion"] == pytest.approx(1.0 - math.exp(-K1 * 40.0), rel=1e-6)
    assert first == pytest.approx(1.0 - math.exp(-K1 * 20.0), rel=1e-6)


# --- cycle and averaging ----------------------------------------------------------------------

def test_cycle_time_averages_the_product():
    r = run({"phase": "gas", "batch_time": "20 s", "T": "400 K", "cycle_time": "100 s",
             "charge_time": "10 s", "down_time": "20 s"})
    u = r.units["B"]
    X = 1.0 - math.exp(-K1 * 20.0)
    assert u["charge"] == pytest.approx(100.0)
    assert u["vessels_needed"] == pytest.approx(0.5)
    out = r.streams["OUT"]
    assert out["component_flows"]["isobutane"] == pytest.approx(X, rel=1e-6)
    assert mass(out) == pytest.approx(mass(r.streams["IN"]), rel=1e-9)
    assert u["productivity"] == pytest.approx(X * 58.12e-3, rel=1e-3)


def test_back_to_back_and_given_charge():
    r = run({"phase": "gas", "batch_time": "20 s", "T": "400 K", "charge_time": "10 s",
             "down_time": "30 s"})
    assert r.units["B"]["cycle_time"] == pytest.approx(60.0)
    assert r.units["B"]["charge"] == pytest.approx(60.0)
    r = run({"phase": "gas", "batch_time": "20 s", "T": "400 K", "charge": "15 mol",
             "down_time": "10 s"})
    u = r.units["B"]
    assert u["cycle_time"] == pytest.approx(15.0)
    assert u["vessels_needed"] == pytest.approx(2.0)
    assert any("2 vessels" in w for w in r.warnings)
    assert mass(r.streams["OUT"]) == pytest.approx(mass(r.streams["IN"]), rel=1e-9)


def test_zero_feed_gives_empty_outlets():
    pkg = get_package(tuple(ISOMER), "ideal")
    unit = BatchReactor("B", {"reactions": ["ISO"], "phase": "gas", "batch_time": 20.0})
    unit.reactions = parse_reactions([ReactionModel(**FIRST_ORDER)], pkg)
    out = unit.calculate(pkg, {"in": [empty_stream(pkg, 400.0, 1e5)]})
    assert out["out"].total_flow == 0.0 and out["vent"].total_flow == 0.0


# --- errors -----------------------------------------------------------------------------------

def error(r):
    assert not r.ok
    return r.error["message"]


def test_missing_kinetics_is_explained():
    plain = {"id": "ISO", "stoichiometry": FIRST_ORDER["stoichiometry"]}
    msg = error(run({"batch_time": "1 h"}, reactions=[plain], ok=False))
    assert "no kinetics" in msg


def test_unreachable_conversion_is_explained():
    msg = error(run({"phase": "gas", "batch_time": "10 s", "stop_conversion": 0.99,
                     "T": "400 K"}, ok=False))
    assert "not reached within the batch time" in msg and "0.3935" in msg


def test_isothermal_with_stop_temperature_is_rejected():
    msg = error(run({"phase": "gas", "batch_time": "10 s", "T": "400 K", "stop_T": "450 K"},
                    ok=False))
    assert "never reaches a stop temperature" in msg


def test_vent_needs_a_liquid_phase_and_an_energy_balance():
    msg = error(run({"phase": "gas", "batch_time": "10 s", "vent": "yes"}, ok=False))
    assert "vent needs" in msg
    msg = error(run({"batch_time": "10 s", "vent": "yes", "T": "300 K"}, ok=False))
    assert "cannot be prescribed" in msg


def test_jacket_with_unknown_charge_and_stop_criterion_is_rejected():
    msg = error(run({"phase": "gas", "batch_time": "10 min", "stop_conversion": 0.5,
                     "area": "1 m2", "U": 50, "utility": "CW"}, utilities=[CW], ok=False))
    assert "charge must be known" in msg


def test_utility_without_jacket_is_rejected():
    msg = error(run({"phase": "gas", "batch_time": "10 s", "utility": "CW"}, utilities=[CW],
                    ok=False))
    assert "jacket" in msg


def test_example_flowsheet_solves():
    model = FlowsheetModel.load(EXAMPLES / "ethyl_acetate_batch.json")
    r = solve(model)
    assert r.ok, r.error
    u = r.units["BATCH"]
    assert u["conversion"] == pytest.approx(0.6, abs=1e-6)
    assert {"conditions", "composition"} <= set(u["profiles"])
