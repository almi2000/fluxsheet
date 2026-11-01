"""Utilities and jacket/coil heat transfer, checked against closed-form results."""

import math

import pytest

from fluxsheet import solve
from fluxsheet.flowsheet import UtilityModel
from fluxsheet.utilities import lmtd, parse_utilities

STEAM = {"id": "LPS", "type": "steam", "P": "3 bar"}
CW = {"id": "CW", "type": "fluid", "T_in": "25 degC", "T_out": "35 degC", "cp": "4.18 kJ/kgK"}
OIL = {"id": "OIL", "type": "fluid", "T_in": "250 degC", "T_out": "200 degC",
       "cp": "2.3 kJ/kgK"}
WATER = {"water": "20 kmol/h"}
SHIFT = {"id": "WGS", "stoichiometry": {"carbon monoxide": -1, "water": -1,
                                        "carbon dioxide": 1, "hydrogen": 1}}
SYNGAS = ["carbon monoxide", "water", "carbon dioxide", "hydrogen"]


def one_unit(unit, *, flows=WATER, T="25 degC", P="3 bar", components=("water",),
             utilities=(STEAM, CW, OIL), reactions=(), outlets=("out",)):
    r = solve({
        "components": list(components),
        "property_method": "ideal",
        "utilities": list(utilities),
        "reactions": list(reactions),
        "units": [{"id": "F", "type": "Feed", "params": {"T": T, "P": P, "flows": flows}},
                  {"id": "U", **unit}],
        "streams": [{"id": "IN", "from": "F.out", "to": "U.in"}]
        + [{"id": port.upper(), "from": f"U.{port}"} for port in outlets],
    })
    return r


def utility(model):
    return next(iter(parse_utilities([UtilityModel(**model)]).values()))


def test_steam_saturation_and_latent_heat_match_the_steam_tables():
    s = utility(STEAM)
    # IAPWS steam tables at 3 bar: 133.52 °C, h_fg = 2163.5 kJ/kg.
    assert s.T_in - 273.15 == pytest.approx(133.52, abs=0.01)
    assert s.latent / 1e3 == pytest.approx(2163.5, abs=0.5)


def test_steam_heated_heater_follows_the_exponential_approach():
    UA = 500.0 * 2.0
    r = one_unit({"type": "Heater", "params": {"utility": "LPS", "U": "500 W/m2K",
                                                "area": "2 m2"}})
    assert r.ok, r.error
    u, T_in, T_out = r.units["U"], r.streams["IN"]["T"], r.streams["OUT"]["T"]
    Ts = utility(STEAM).T_in
    # Condensing steam is at one temperature, so Q = UA * LMTD whatever the flow direction.
    assert u["duty"] == pytest.approx(UA * lmtd(Ts - T_in, Ts - T_out), rel=1e-6)
    # ... and with a constant cp the outlet approaches the steam exponentially in NTU.
    C = u["duty"] / (T_out - T_in)
    assert T_out == pytest.approx(Ts - (Ts - T_in) * math.exp(-UA / C), abs=1e-6)
    assert u["utility_flow"] == pytest.approx(u["duty"] / utility(STEAM).latent, rel=1e-9)
    assert u["UA"] == pytest.approx(UA) and u["heat_transfer_area"] == pytest.approx(2.0)


def test_cooling_water_consumption_and_required_area_for_a_given_temperature():
    r = one_unit({"type": "Flash", "params": {"T": "60 degC", "utility": "CW",
                                               "U": "400 W/m2K"}},
                 T="90 degC", outlets=("vapor", "liquid"))
    assert r.ok, r.error
    u = r.units["U"]
    Q = u["duty"]
    assert Q < 0
    assert u["utility_flow"] == pytest.approx(-Q / (4180.0 * 10.0), rel=1e-9)
    # A well-mixed drum is at 60 °C throughout; the water warms from 25 to 35 °C.
    mean = lmtd(60.0 - 25.0, 60.0 - 35.0)
    assert u["LMTD"] == pytest.approx(mean, rel=1e-9)
    assert u["UA"] == pytest.approx(-Q / mean, rel=1e-9)
    assert u["heat_transfer_area"] == pytest.approx(-Q / mean / 400.0, rel=1e-9)


def test_counter_current_needs_less_area_than_cocurrent():
    def area(arrangement):
        r = one_unit({"type": "Heater", "params": {"T": "150 degC", "utility": "OIL",
                                                    "U": "300 W/m2K",
                                                    "arrangement": arrangement}},
                     P="10 bar")
        assert r.ok, r.error
        return r.units["U"]["heat_transfer_area"]

    # Hot oil 250 -> 200 °C heats water 25 -> 150 °C.
    counter = lmtd(200.0 - 25.0, 250.0 - 150.0)
    co = lmtd(250.0 - 25.0, 200.0 - 150.0)
    assert area("cocurrent") / area("countercurrent") == pytest.approx(counter / co, rel=1e-9)


def test_cooled_reactor_settles_where_jacket_and_energy_balance_agree():
    reactor = {"type": "StoichiometricReactor",
               "params": {"reactions": ["WGS"], "conversions": ["0.6"]}}
    feed = {"carbon monoxide": "50 kmol/h", "water": "60 kmol/h"}
    common = {"flows": feed, "T": "350 degC", "P": "20 bar", "components": SYNGAS,
              "reactions": [SHIFT]}
    adiabatic = one_unit(reactor, **common)
    cooled = one_unit({**reactor, "params": {**reactor["params"], "utility": "CW",
                                              "U": "150 W/m2K", "area": "20 m2"}}, **common)
    assert adiabatic.ok and cooled.ok, cooled.error
    u, T = cooled.units["U"], cooled.units["U"]["T_out"]
    assert u["duty"] < 0 and T < adiabatic.units["U"]["T_out"]
    assert u["duty"] == pytest.approx(-3000.0 * lmtd(T - 308.15, T - 298.15), rel=1e-6)


def test_fixed_coolant_flow_through_a_mixed_vessel_matches_ntu():
    reactor = {"type": "StoichiometricReactor",
               "params": {"reactions": ["WGS"], "conversions": ["0.6"], "utility": "CW",
                          "utility_flow": "2 kg/s", "U": "150 W/m2K", "area": "20 m2"}}
    r = one_unit(reactor, flows={"carbon monoxide": "50 kmol/h", "water": "60 kmol/h"},
                 T="350 degC", P="20 bar", components=SYNGAS, reactions=[SHIFT])
    assert r.ok, r.error
    u = r.units["U"]
    C, T = 2.0 * 4180.0, u["T_out"]
    T_cw = T + (298.15 - T) * math.exp(-3000.0 / C)
    assert u["utility_T_out"] == pytest.approx(T_cw, abs=1e-6)
    assert u["duty"] == pytest.approx(-C * (T_cw - 298.15), rel=1e-6)


@pytest.mark.parametrize("unit, message", [
    ({"type": "Flash", "params": {"T": "30 degC", "utility": "CW"}}, "temperature cross"),
    ({"type": "Heater", "params": {"T": "10 degC", "utility": "LPS"}}, "can only heat"),
    ({"type": "Heater", "params": {"T": "80 degC", "utility": "HPS"}}, "unknown utility"),
    ({"type": "Heater", "params": {"area": "1 m2", "utility": "LPS"}}, "not U"),
    ({"type": "Heater", "params": {"area": "1 m2", "U": 100}}, "needs a heating or cooling"),
])
def test_utility_mistakes_are_explained(unit, message):
    outlets = ("vapor", "liquid") if unit["type"] == "Flash" else ("out",)
    r = one_unit(unit, T="50 degC", outlets=outlets)
    assert not r.ok and message in r.error["message"]
    assert r.error["where"] == "U"
