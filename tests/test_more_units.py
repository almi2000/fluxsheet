"""Thin-film evaporator, mixing tank, turbine, multistage compressor, pipe and extraction
column, checked against balances and closed-form results."""

import itertools
import math

import pytest

from fluxsheet import solve
from fluxsheet.thermo_package import get_package
from fluxsheet.units.pressure_pipe import Pipe, friction_factor

GLYCOL = ["water", "ethylene glycol"]
GLYCOL_FEED = {"water": "60 kmol/h", "ethylene glycol": "40 kmol/h"}


def flowsheet(components, feeds, unit, inlets, outlets, method="peng-robinson", reactions=()):
    """Feeds -> one unit 'U'. ``inlets``: feed id -> port; outlets: port names."""
    return solve({
        "components": list(components),
        "property_method": method,
        "units": [{"id": fid, "type": "Feed", "params": params} for fid, params in feeds.items()]
        + [{"id": "U", **unit}],
        "streams": [{"id": fid, "from": f"{fid}.out", "to": f"U.{port}"}
                    for fid, port in inlets.items()]
        + [{"id": port.upper(), "from": f"U.{port}"} for port in outlets],
    })


def closure(r, inlets, outlets, key):
    return (sum(r.streams[i][key] for i in inlets) - sum(r.streams[o][key] for o in outlets))


# --- thin-film evaporator ---------------------------------------------------------------

def evaporator(params, unit="ThinFilmEvaporator", outlets=("vapor", "residue")):
    return flowsheet(GLYCOL, {"F": {"T": "40 degC", "P": "1 bar", "flows": GLYCOL_FEED}},
                     {"type": unit, "params": params}, {"F": "in"}, outlets, method="unifac")


def test_thin_film_evaporator_balances_and_meets_its_fraction():
    r = evaporator({"VF": 0.6, "P": "0.1 bar"})
    assert r.ok, r.error
    u = r.units["U"]
    assert u["vapor_fraction"] == pytest.approx(0.6, abs=1e-6)
    assert abs(closure(r, ["F"], ["VAPOR", "RESIDUE"], "molar_flow")) < 1e-9
    h = closure(r, ["F"], ["VAPOR", "RESIDUE"], "enthalpy_flow") + u["duty"]
    assert abs(h) < 1e-6 * u["duty"]
    T = u["profiles"]["film"]["series"][0]["values"]
    assert all(b >= a - 1e-6 for a, b in itertools.pairwise(T))  # the film gets hotter downwards


def test_film_evaporation_separates_better_than_one_flash():
    film = evaporator({"VF": 0.6, "P": "0.1 bar"})
    flash = evaporator({"duty": film.units["U"]["duty"], "P": "0.1 bar"}, unit="Flash",
                       outlets=("vapor", "liquid"))
    x_film = film.streams["RESIDUE"]["mole_fractions"]["ethylene glycol"]
    x_flash = flash.streams["LIQUID"]["mole_fractions"]["ethylene glycol"]
    assert x_film > x_flash + 0.05


def test_feed_flashing_beyond_the_target_is_explained():
    r = flowsheet(GLYCOL, {"F": {"T": "95 degC", "P": "1 bar", "flows": GLYCOL_FEED}},
                  {"type": "ThinFilmEvaporator", "params": {"VF": 0.01, "P": "0.05 bar"}},
                  {"F": "in"}, ("vapor", "residue"), method="unifac")
    assert not r.ok and "already flashes" in r.error["message"]


# --- mixing tank ------------------------------------------------------------------------

def test_mixing_tank_blends_like_a_mixer_and_reports_residence_time():
    feeds = {"A": {"T": "20 degC", "P": "2 bar", "flows": {"water": "50 kmol/h"}},
             "B": {"T": "80 degC", "P": "2 bar", "flows": {"water": "30 kmol/h"}}}
    tank = flowsheet(["water"], feeds, {"type": "MixingTank", "params": {"V": "2 m3"}},
                     {"A": "in", "B": "in"}, ("liquid",), method="ideal")
    mixer = flowsheet(["water"], feeds, {"type": "Mixer"}, {"A": "in", "B": "in"}, ("out",),
                      method="ideal")
    assert tank.ok, tank.error
    assert not tank.warnings  # the unused vent is fine while no vapor forms
    out = tank.streams["LIQUID"]
    assert out["T"] == pytest.approx(mixer.streams["OUT"]["T"], abs=1e-6)
    pkg = get_package(("water",), "ideal")
    V_liquid = out["molar_flow"] * pkg.flash([1.0], T=out["T"], P=out["P"]).V
    assert tank.units["U"]["residence_time"] == pytest.approx(2.0 / V_liquid, rel=1e-9)


# --- turbine and multistage compressor ------------------------------------------------------

GAS = {"F": {"T": "200 degC", "P": "20 bar", "flows": {"nitrogen": "100 kmol/h"}}}


def test_turbine_produces_efficiency_times_isentropic_work():
    def turbine(eff):
        r = flowsheet(["nitrogen"], GAS, {"type": "Turbine",
                                          "params": {"P": "2 bar", "efficiency": eff}},
                      {"F": "in"}, ("out",))
        assert r.ok, r.error
        return r

    ideal, real = turbine(1.0), turbine(0.8)
    pkg = get_package(("nitrogen",), "peng-robinson")
    s_in = pkg.flash([1.0], T=473.15, P=20e5).S
    out = ideal.streams["OUT"]
    assert pkg.flash([1.0], T=out["T"], P=out["P"]).S == pytest.approx(s_in, abs=1e-6)
    assert ideal.units["U"]["power"] < 0
    assert real.units["U"]["power"] == pytest.approx(0.8 * ideal.units["U"]["power"], rel=1e-6)
    assert real.streams["OUT"]["T"] > out["T"]


def test_intercooled_compressor_needs_less_power_than_one_stage():
    gas = {"F": {"T": "30 degC", "P": "1 bar", "flows": {"nitrogen": "100 kmol/h"}}}
    single = flowsheet(["nitrogen"], gas, {"type": "Compressor", "params": {"P": "27 bar"}},
                       {"F": "in"}, ("out",))
    one = flowsheet(["nitrogen"], gas, {"type": "MultistageCompressor",
                                        "params": {"P": "27 bar", "stages": 1,
                                                   "efficiency": 0.72}},
                    {"F": "in"}, ("out",))
    three = flowsheet(["nitrogen"], gas, {"type": "MultistageCompressor",
                                          "params": {"P": "27 bar", "stages": 3,
                                                     "efficiency": 0.72,
                                                     "T_cooler": "30 degC"}},
                      {"F": "in"}, ("out",))
    assert one.units["U"]["power"] == pytest.approx(single.units["U"]["power"], rel=1e-9)
    u = three.units["U"]
    assert u["pressure_ratio"] == pytest.approx(3.0, rel=1e-9)
    assert u["power"] < 0.8 * single.units["U"]["power"]
    assert u["T_max"] < single.streams["OUT"]["T"]
    assert three.streams["OUT"]["P"] == pytest.approx(27e5)


def test_intercoolers_drain_condensed_water():
    gas = {"F": {"T": "40 degC", "P": "1 bar",
                 "flows": {"methane": "95 kmol/h", "water": "5 kmol/h"}}}
    r = flowsheet(["methane", "water"], gas,
                  {"type": "MultistageCompressor",
                   "params": {"P": "30 bar", "stages": 3, "T_cooler": "35 degC",
                              "aftercooler": "yes"}},
                  {"F": "in"}, ("out", "condensate"))
    assert r.ok, r.error
    water = r.streams["CONDENSATE"]["component_flows"]["water"] * 3.6
    assert 4.0 < water < 5.0  # most of it; the rest stays in the gas as vapor
    assert abs(closure(r, ["F"], ["OUT", "CONDENSATE"], "molar_flow")) < 1e-9


# --- pipe ----------------------------------------------------------------------------------

def test_laminar_pipe_follows_hagen_poiseuille():
    feed = {"F": {"T": "25 degC", "P": "3 bar", "flows": {"water": "0.05 kmol/h"}}}
    r = flowsheet(["water"], feed, {"type": "Pipe", "params": {"length": "50 m",
                                                                "diameter": "10 mm"}},
                  {"F": "in"}, ("out",), method="ideal")
    assert r.ok, r.error
    u = r.units["U"]
    assert u["reynolds"] < 2300
    pkg = get_package(("water",), "ideal")
    st = pkg.flash([1.0], T=298.15, P=3e5)
    mu = Pipe._viscosity(pkg, st, [1.0])
    Q = 0.05 / 3.6 * st.V  # m3/s
    assert u["dP"] == pytest.approx(128 * mu * 50.0 * Q / (math.pi * 0.01 ** 4), rel=1e-3)


def test_swamee_jain_is_close_to_colebrook():
    for Re, eps in ((1e4, 0.0), (1e5, 1e-4), (1e6, 1e-3)):
        f = friction_factor(Re, eps)
        colebrook = 1 / f ** 0.5 + 2 * math.log10(eps / 3.7 + 2.51 / (Re * f ** 0.5))
        assert abs(colebrook) < 0.03 / f ** 0.5


def test_rising_pipe_loses_the_static_head():
    feed = {"F": {"T": "25 degC", "P": "3 bar", "flows": {"water": "0.01 kmol/h"}}}
    r = flowsheet(["water"], feed, {"type": "Pipe", "params": {"length": "20 m",
                                                                "diameter": "50 mm",
                                                                "elevation": "10 m"}},
                  {"F": "in"}, ("out",), method="ideal")
    rho = 0.018015 / get_package(("water",), "ideal").flash([1.0], T=298.15, P=3e5).V
    assert r.units["U"]["dP"] == pytest.approx(rho * 9.80665 * 10.0, rel=1e-3)


# --- extraction column -----------------------------------------------------------------------

def extraction(stages, method="unifac"):
    return flowsheet(
        ["water", "acetic acid", "ethyl acetate"],
        {"F": {"T": "25 degC", "P": "1 atm",
               "flows": {"water": "80 kmol/h", "acetic acid": "10 kmol/h"}},
         "S": {"T": "25 degC", "P": "1 atm", "flows": {"ethyl acetate": "60 kmol/h"}}},
        {"type": "ExtractionColumn", "params": {"stages": stages}},
        {"F": "feed", "S": "solvent"}, ("extract", "raffinate"), method=method)


def test_extraction_column_balances_and_more_stages_extract_more():
    one, five = extraction(1), extraction(5)
    assert five.ok, five.error
    for r in (one, five):
        for c in ("water", "acetic acid", "ethyl acetate"):
            fin = sum(r.streams[i]["component_flows"][c] for i in ("F", "S"))
            out = sum(r.streams[o]["component_flows"][c] for o in ("EXTRACT", "RAFFINATE"))
            assert out == pytest.approx(fin, rel=1e-9)
    acid = [r.streams["RAFFINATE"]["component_flows"]["acetic acid"] for r in (one, five)]
    assert acid[1] < 0.2 * acid[0]
    # The extract is the ethyl acetate phase.
    assert five.streams["EXTRACT"]["mole_fractions"]["ethyl acetate"] > 0.5


def test_extraction_needs_a_method_with_two_liquids():
    r = extraction(3, method="peng-robinson")
    assert not r.ok and "two liquid phases" in r.error["message"]


@pytest.mark.parametrize("flowing_side", ["hot", "cold"])
def test_exchanger_with_empty_side_throttles_flowing_side_isenthalpically(flowing_side):
    empty_side = "cold" if flowing_side == "hot" else "hot"
    result = solve({
        "components": ["propane"], "property_method": "peng-robinson",
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "T": 300, "P": "20 bar", "flows": {"propane": 1}}},
            {"id": "S", "type": "Splitter", "params": {"split": 1}},
            {"id": "H", "type": "HeatExchanger", "params": {
                "UA": 0, f"dP_{flowing_side}": "19 bar"}},
        ],
        "streams": [
            {"id": "feed", "from": "F.out", "to": "S.in"},
            {"id": "in", "from": "S.out1", "to": f"H.{flowing_side}_in"},
            {"id": "empty", "from": "S.out2", "to": f"H.{empty_side}_in"},
            {"id": "out", "from": f"H.{flowing_side}_out"},
            {"id": "empty_out", "from": f"H.{empty_side}_out"},
        ],
    })
    assert result.ok, result.error
    assert result.units["H"]["duty"] == 0
    assert result.streams["out"]["enthalpy_flow"] == pytest.approx(
        result.streams["in"]["enthalpy_flow"], abs=1e-5)
    reference = flowsheet(
        ["propane"], {"F": {"T": 300, "P": "20 bar", "flows": {"propane": 1}}},
        {"type": "Valve", "params": {"P": "1 bar"}}, {"F": "in"}, ("out",),
    )
    assert reference.ok, reference.error
    assert result.streams["out"]["vapor_fraction"] == pytest.approx(
        reference.streams["OUT"]["vapor_fraction"], abs=1e-8)


@pytest.mark.parametrize("temperature, ambient", [(320, 300), (300, 320)])
def test_pipe_heat_transfer_stays_between_inlet_and_ambient(temperature, ambient):
    results = []
    for segments in (1, 20, 100):
        result = flowsheet(
            ["nitrogen"], {"F": {"T": temperature, "P": "2 bar", "flows": {"nitrogen": 1}}},
            {"type": "Pipe", "params": {"length": 3, "diameter": 0.1,
             "segments": segments, "U": 50, "T_ambient": ambient}},
            {"F": "in"}, ("out",))
        assert result.ok, result.error
        outlet = result.streams["OUT"]
        assert min(temperature, ambient) <= outlet["T"] <= max(temperature, ambient)
        assert outlet["enthalpy_flow"] == pytest.approx(
            result.streams["F"]["enthalpy_flow"] + result.units["U"]["duty"], abs=1e-5)
        results.append(outlet["T"])
    assert abs(results[2] - results[1]) < abs(results[1] - results[0])


def test_extraction_stages_reconstruct_component_and_phase_equilibria():
    import numpy as np

    result = extraction(5)
    assert result.ok, result.error
    components = ("water", "acetic acid", "ethyl acetate")
    profile = {s["label"]: s["values"]
               for s in result.units["U"]["profiles"]["stages"]["series"]}
    raffinate = np.array([[profile[f"Raffinate x {name}"][j]
                            * profile["Raffinate flow"][j] for name in components]
                          for j in range(5)])
    feed = np.array([result.streams["F"]["component_flows"][name] for name in components])
    solvent = np.array([result.streams["S"]["component_flows"][name] for name in components])
    extract = np.array([result.streams["EXTRACT"]["component_flows"][name]
                        for name in components])
    pkg = get_package(components, "unifac")
    previous_raffinate = feed
    for j in range(5):
        assert extract.sum() == pytest.approx(profile["Extract flow"][j], rel=1e-8)
        x_extract = extract / extract.sum()
        x_raffinate = raffinate[j] / raffinate[j].sum()
        phi_extract = pkg.phase_props("liquid", 298.15, 101325, list(x_extract)).lnphi
        phi_raffinate = pkg.phase_props("liquid", 298.15, 101325, list(x_raffinate)).lnphi
        assert np.allclose(x_extract * np.exp(phi_extract), x_raffinate * np.exp(phi_raffinate),
                           rtol=1e-5, atol=1e-8)
        extract = extract + raffinate[j] - previous_raffinate
        previous_raffinate = raffinate[j]
    assert np.allclose(extract, solvent, rtol=1e-7, atol=1e-8)


@pytest.mark.parametrize("arrangement", ["countercurrent", "cocurrent"])
def test_gas_exchanger_matches_constant_heat_capacity_solution(arrangement):
    import math

    pkg = get_package(("nitrogen",), "peng-robinson")
    cp = (pkg.phase_props("vapor", 311.0, 2e5, [1.0]).H
          - pkg.phase_props("vapor", 309.0, 2e5, [1.0]).H) / 2
    ua = 25.0
    ntu = ua / cp
    effectiveness = (ntu / (1 + ntu) if arrangement == "countercurrent"
                     else (1 - math.exp(-2 * ntu)) / 2)
    expected = effectiveness * cp * 20.0
    result = flowsheet(
        ["nitrogen"],
        {"H": {"T": 320, "P": "2 bar", "flows": {"nitrogen": 1}},
         "C": {"T": 300, "P": "2 bar", "flows": {"nitrogen": 1}}},
        {"type": "HeatExchanger", "params": {"UA": ua, "arrangement": arrangement}},
        {"H": "hot_in", "C": "cold_in"}, ("hot_out", "cold_out"),
    )
    assert result.ok, result.error
    assert result.units["U"]["duty"] == pytest.approx(expected, rel=0.01)


def test_pipe_grid_converges_to_constant_heat_capacity_solution():
    import math

    pkg = get_package(("nitrogen",), "peng-robinson")
    cp = (pkg.phase_props("vapor", 311.0, 2e5, [1.0]).H
          - pkg.phase_props("vapor", 309.0, 2e5, [1.0]).H) / 2
    expected = 300.0 + 20.0 * math.exp(-50 * math.pi * 0.1 * 3 / cp)
    errors = []
    for segments in (1, 20, 100):
        result = flowsheet(
            ["nitrogen"], {"F": {"T": 320, "P": "2 bar", "flows": {"nitrogen": 1}}},
            {"type": "Pipe", "params": {"length": 3, "diameter": 0.1,
             "segments": segments, "U": 50, "T_ambient": 300}},
            {"F": "in"}, ("out",))
        assert result.ok, result.error
        errors.append(abs(result.streams["OUT"]["T"] - expected))
    assert errors[2] < errors[1] < errors[0]
    assert errors[2] < 0.1


@pytest.mark.parametrize("unit,component,temperature", [
    ("Valve", "nitrogen", 320),
    ("Pump", "water", 300),
    ("Compressor", "nitrogen", 320),
    ("Turbine", "nitrogen", 320),
])
def test_pressure_machine_small_pressure_change_conserves_energy(unit, component, temperature):
    params = {"dP": 1.0}
    result = flowsheet(
        [component], {"F": {"T": temperature, "P": "2 bar", "flows": {component: 1}}},
        {"type": unit, "params": params}, {"F": "in"}, ("out",))
    assert result.ok, result.error
    power = result.units["U"].get("power", 0.0)
    assert result.streams["OUT"]["enthalpy_flow"] == pytest.approx(
        result.streams["F"]["enthalpy_flow"] + power, abs=1e-3)
    assert abs(power) < 0.03
    if unit in ("Pump", "Compressor"):
        assert power >= -1e-6
    elif unit == "Turbine":
        assert power <= 1e-6
