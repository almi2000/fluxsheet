import math

import numpy as np
import pytest

from fluxsheet import solve
from fluxsheet.columns import (
    ColumnSetup,
    fenske_min_stages,
    gilliland_molokanov,
    kirkbride_ratio,
    solve_column,
    underwood_theta,
)
from fluxsheet.errors import SimulationError
from fluxsheet.stream import Stream
from fluxsheet.thermo_package import get_package
from fluxsheet.units.columns import _feed

BT = {"benzene": "50 kmol/h", "toluene": "50 kmol/h"}
KMOLH = 1 / 3.6  # mol/s


def column(params, *, unit="DistillationColumn", flows=BT, T="90 degC", P="1.5 bar",
           components=("benzene", "toluene"), method="peng-robinson"):
    outlets = {"ShortcutColumn": ("distillate", "bottoms"),
               "DistillationColumn": ("distillate", "bottoms", "vapor_distillate")}[unit]
    return solve({
        "components": list(components),
        "property_method": method,
        "units": [{"id": "F", "type": "Feed", "params": {"T": T, "P": P, "flows": flows}},
                  {"id": "C", "type": unit, "params": params}],
        "streams": [{"id": "FEED", "from": "F.out", "to": "C.in"}]
        + [{"id": port, "from": f"C.{port}"} for port in outlets],
    })


def balances(r, inlets, outlets):
    """Relative closure of every component balance and of the energy balance."""
    s, u = r.streams, next(iter(u for u in r.units.values() if "profiles" in u))
    names = s[inlets[0]]["component_flows"]
    total_in = sum(s[i]["molar_flow"] for i in inlets)
    comp = max(abs(sum(s[i]["component_flows"][c] for i in inlets)
                   - sum(s[o]["component_flows"][c] for o in outlets)) for c in names) / total_in
    duties = u.get("condenser_duty", 0.0) + u.get("reboiler_duty", 0.0)
    h_in = sum(s[i]["enthalpy_flow"] for i in inlets)
    h_out = sum(s[o]["enthalpy_flow"] for o in outlets)
    scale = max(abs(u.get("condenser_duty", 0.0)), abs(u.get("reboiler_duty", 0.0)),
                sum(abs(s[k]["enthalpy_flow"]) for k in [*inlets, *outlets]))
    return comp, abs(h_in + duties - h_out) / scale


# --- shortcut correlations, against hand calculations -------------------------------------

def test_fenske_hand_calculation():
    # 99 % of benzene up, 99 % of toluene down, constant relative volatility 2.5:
    # Nmin = ln[(0.99/0.01)(0.99/0.01)] / ln 2.5 = 9.19024 / 0.91629 = 10.0298
    assert fenske_min_stages(0.99, 0.01, 0.01, 0.99, 2.5) == pytest.approx(10.0298, abs=1e-4)


def test_underwood_binary_saturated_liquid_matches_closed_form():
    # For a binary saturated-liquid feed, Underwood reduces to
    # Rmin = [xD/z - alpha (1 - xD)/(1 - z)] / (alpha - 1).
    alpha, z, xd = np.array([2.5, 1.0]), np.array([0.4, 0.6]), np.array([0.95, 0.05])
    theta = underwood_theta(alpha, z, 1.0, 1.0, 2.5)
    r_min = float(np.sum(alpha * xd / (alpha - theta))) - 1.0
    closed = (xd[0] / z[0] - alpha[0] * xd[1] / z[1]) / (alpha[0] - 1.0)
    assert r_min == pytest.approx(closed, rel=1e-9)


def test_gilliland_limits():
    assert gilliland_molokanov(10.0, 1e7, 1.0) == pytest.approx(10.0, rel=1e-3)  # total reflux
    near_min = [gilliland_molokanov(10.0, r, 1.0) for r in (1.5, 1.1, 1.01)]
    assert near_min == sorted(near_min) and near_min[-1] > 4 * 10.0  # N grows toward Rmin
    # A textbook point: (R - Rmin)/(R + 1) = 0.5 gives (N - Nmin)/(N + 1) = 0.249.
    n = gilliland_molokanov(10.0, 3.0, 1.0)
    assert (n - 10.0) / (n + 1.0) == pytest.approx(0.2489, abs=1e-3)


def test_kirkbride_symmetric_split_puts_feed_in_the_middle():
    assert kirkbride_ratio(0.5, 0.5, 0.01, 0.01, 50.0, 50.0) == pytest.approx(1.0)


# --- shortcut column ---------------------------------------------------------------------

SHORTCUT = {"light_key": "benzene", "heavy_key": "toluene", "lk_recovery": 0.99,
            "hk_recovery": 0.01, "reflux_factor": 1.3, "P": "1 atm"}


@pytest.fixture(scope="module")
def shortcut():
    r = column(SHORTCUT, unit="ShortcutColumn")
    assert r.ok, r.error
    return r


def test_shortcut_benzene_toluene(shortcut):
    u = shortcut.units["C"]
    assert 2.3 < u["alpha"] < 2.6  # benzene/toluene at 1 atm
    # Fenske with the reported (geometric-mean) volatility, by hand.
    assert u["min_stages"] == pytest.approx(math.log(99 * 99) / math.log(u["alpha"]), rel=1e-9)
    assert u["reflux_ratio"] == pytest.approx(1.3 * u["min_reflux"])
    assert 1.2 < u["min_reflux"] < 1.6
    assert u["min_stages"] < u["theoretical_stages"] < 3 * u["min_stages"]
    assert u["stages"] == math.ceil(u["theoretical_stages"]) + 1  # plus the total condenser
    assert abs(u["feed_stage"] - u["stages"] / 2) <= 2  # symmetric split: feed mid-column
    s = shortcut.streams
    assert s["distillate"]["component_flows"]["benzene"] == pytest.approx(0.99 * 50 * KMOLH)
    assert s["distillate"]["phase"] == "liquid" and s["bottoms"]["phase"] == "liquid"
    assert u["condenser_duty"] < 0 < u["reboiler_duty"]


@pytest.mark.parametrize("condenser", ["total", "partial"])
def test_rigorous_column_reproduces_shortcut_design(condenser):
    short = column({**SHORTCUT, "condenser": condenser}, unit="ShortcutColumn")
    assert short.ok, short.error
    u, s = short.units["C"], short.streams
    r = column({"stages": u["stages"], "feed_stage": u["feed_stage"], "P": "1 atm",
                "condenser": condenser, "reflux_ratio": u["reflux_ratio"],
                "distillate_rate": s["distillate"]["molar_flow"]})
    assert r.ok, r.error
    top = "distillate" if condenser == "total" else "vapor_distillate"
    d = r.streams[top]["component_flows"]
    # Gilliland is a correlation; a rigorous column at its N and R should land near the spec.
    assert d["benzene"] / (50 * KMOLH) == pytest.approx(0.99, abs=0.005)
    assert d["toluene"] / (50 * KMOLH) == pytest.approx(0.01, abs=0.005)
    rig = r.units["C"]
    assert rig["condenser_duty"] == pytest.approx(u["condenser_duty"], rel=0.02)
    assert rig["reboiler_duty"] == pytest.approx(u["reboiler_duty"], rel=0.02)


@pytest.mark.parametrize("params, message", [
    ({**SHORTCUT, "reflux_factor": None, "reflux_ratio": 1.0}, "below the minimum"),
    ({**SHORTCUT, "lk_recovery": 0.4, "hk_recovery": 0.6}, "infeasible key recoveries"),
    ({**SHORTCUT, "light_key": "toluene", "heavy_key": "benzene"}, "not more volatile"),
    ({**SHORTCUT, "heavy_key": "water"}, "not in the flowsheet"),
])
def test_shortcut_errors_explain_themselves(params, message):
    r = column(params, unit="ShortcutColumn")
    assert not r.ok
    assert message in r.error["message"]
    assert r.error["where"] == "C" and r.error["hint"]


# --- rigorous column ---------------------------------------------------------------------

BASE = {"stages": 20, "feed_stage": 10, "P": "1 atm"}


@pytest.mark.parametrize("params", [
    {"reflux_ratio": 2, "distillate_rate": "50 kmol/h"},
    {"reflux_ratio": 2, "bottoms_rate": "50 kmol/h", "condenser": "partial"},
    {"boilup_ratio": 3, "distillate_rate": "45 kmol/h", "dP_stage": "0.7 kPa"},
    {"reflux_ratio": 2, "distillate_rate": "50 kmol/h", "efficiency": 0.7},
    {"condenser": "none", "feed_stage": 1, "boilup_ratio": 1.5},
], ids=["total", "partial", "boilup-dP", "murphree", "reboiled-stripper"])
def test_column_balances_close(params):
    r = column({**BASE, **params})
    assert r.ok, r.error
    assert not r.warnings
    outlets = ["distillate", "bottoms", "vapor_distillate"]
    comp, energy = balances(r, ["FEED"], [o for o in outlets if o in r.streams])
    assert comp < 1e-6 and energy < 1e-6
    u = r.units["C"]
    T = u["profiles"]["temperature"]["series"][0]["values"]
    assert len(T) == 20 and T[0] < T[-1]  # cold top, hot bottom
    if params.get("condenser") == "partial":
        assert r.streams["vapor_distillate"]["phase"] == "vapor"
        assert r.streams["distillate"]["molar_flow"] == 0.0


def test_total_condenser_hides_its_empty_vapor_outlet():
    r = solve({
        "components": ["benzene", "toluene"], "property_method": "peng-robinson",
        "units": [{"id": "F", "type": "Feed",
                   "params": {"T": "90 degC", "P": "1.5 bar", "flows": BT}},
                  {"id": "C", "type": "DistillationColumn",
                   "params": {**BASE, "reflux_ratio": 2, "distillate_rate": "50 kmol/h"}}],
        "streams": [{"id": "FEED", "from": "F.out", "to": "C.in"},
                    {"id": "D", "from": "C.distillate"}, {"id": "B", "from": "C.bottoms"}],
    })
    assert r.ok and not r.warnings
    assert set(r.streams) == {"FEED", "D", "B"}


def test_murphree_efficiency_lowers_separation():
    def benzene_up(eff):
        r = column({**BASE, "reflux_ratio": 2, "distillate_rate": "50 kmol/h",
                    "efficiency": eff})
        return r.streams["distillate"]["component_flows"]["benzene"]

    assert benzene_up(0.6) < benzene_up(0.8) < benzene_up(1.0)


def test_nonideal_column_with_unifac():
    r = column({**BASE, "reflux_ratio": 1.5, "distillate_rate": "40 kmol/h"},
               components=("water", "methanol"), method="unifac", T="70 degC",
               flows={"water": "60 kmol/h", "methanol": "40 kmol/h"})
    assert r.ok, r.error
    d = r.streams["distillate"]["mole_fractions"]
    assert d["methanol"] > 0.99
    comp, energy = balances(r, ["FEED"], ["distillate", "bottoms"])
    assert comp < 1e-6 and energy < 1e-6


def test_profiles_report_every_stage_and_component():
    r = column({**BASE, "reflux_ratio": 2, "distillate_rate": "50 kmol/h"})
    profiles = r.units["C"]["profiles"]
    assert set(profiles) == {"temperature", "flows", "liquid_composition", "vapor_composition"}
    x = {s["label"]: s["values"] for s in profiles["liquid_composition"]["series"]}
    assert set(x) == {"benzene", "toluene"}
    assert all(abs(a + b - 1.0) < 1e-9 for a, b in zip(x["benzene"], x["toluene"]))
    assert profiles["temperature"]["x"] == list(range(1, 21))


@pytest.mark.parametrize("params, message, hint", [
    ({"reflux_ratio": 2}, "under-specified", "distillate rate"),
    ({"reflux_ratio": 2, "boilup_ratio": 3}, "over-specified", "only one"),
    ({"distillate_rate": "40 kmol/h", "bottoms_rate": "60 kmol/h"}, "over-specified",
     "only one"),
    ({"reflux_ratio": 2, "distillate_rate": "40 kmol/h", "condenser": "none"},
     "needs a condenser", "distillate or bottoms rate"),
    ({"boilup_ratio": 2, "distillate_rate": "40 kmol/h", "reboiler": "none"},
     "needs a reboiler", "kettle"),
    ({"reflux_ratio": 2, "distillate_rate": "40 kmol/h", "feed_stage": 1},
     "outside the column", "condenser"),
    ({"reflux_ratio": 2, "distillate_rate": "120 kmol/h"}, "not below the feed", "both"),
])
def test_bad_column_specs_explain_themselves(params, message, hint):
    r = column({**BASE, **params})
    assert not r.ok
    assert message in r.error["message"]
    assert hint in r.error["hint"]
    assert r.error["where"] == "C"


def test_newton_failure_explains_itself():
    pkg = get_package(("benzene", "toluene"), "peng-robinson")
    feed = Stream.from_flash(pkg, [50 * KMOLH, 50 * KMOLH], T=363.15, P=1.5e5)
    setup = ColumnSetup(stages=20, P=np.full(20, 101325.0), feeds=[_feed(9, feed)],
                        specs={"reflux_ratio": 2.0, "distillate_rate": 50 * KMOLH},
                        max_iterations=1)
    with pytest.raises(SimulationError, match="did not converge") as info:
        solve_column(pkg, setup)
    assert "reflux" in info.value.hint


# --- absorber / stripper -----------------------------------------------------------------

def absorber(gas, liquid, T_gas, T_liq, P, stages, components):
    return solve({
        "components": components, "property_method": "peng-robinson",
        "units": [
            {"id": "G", "type": "Feed", "params": {"T": T_gas, "P": P, "flows": gas}},
            {"id": "L", "type": "Feed", "params": {"T": T_liq, "P": P, "flows": liquid}},
            {"id": "A", "type": "Absorber", "params": {"stages": stages}},
        ],
        "streams": [{"id": "GAS", "from": "G.out", "to": "A.gas_in"},
                    {"id": "LIQ", "from": "L.out", "to": "A.liquid_in"},
                    {"id": "GAS_OUT", "from": "A.gas_out"},
                    {"id": "LIQ_OUT", "from": "A.liquid_out"}],
    })


def test_absorber_recovers_propane_into_oil():
    r = absorber({"methane": "85 kmol/h", "ethane": "10 kmol/h", "propane": "5 kmol/h"},
                 {"n-decane": "40 kmol/h"}, "30 degC", "30 degC", "30 bar", 8,
                 ["methane", "ethane", "propane", "n-decane"])
    assert r.ok, r.error
    s = r.streams
    absorbed = {c: s["LIQ_OUT"]["component_flows"][c] / s["GAS"]["component_flows"][c]
                for c in ("methane", "ethane", "propane")}
    assert 0.8 < absorbed["propane"] < 0.97
    assert absorbed["methane"] < absorbed["ethane"] < absorbed["propane"]
    assert absorbed["methane"] < 0.15
    assert s["GAS_OUT"]["component_flows"]["n-decane"] < 0.01 * s["LIQ"]["molar_flow"]
    u = r.units["A"]
    assert u["T_bottom"] > 303.15  # heat of absorption warms the liquid
    comp, energy = balances(r, ["GAS", "LIQ"], ["GAS_OUT", "LIQ_OUT"])
    assert comp < 1e-6 and energy < 1e-6


def test_stripper_removes_light_component_from_liquid():
    # Same unit as a stripper: nitrogen strips propane out of n-decane at low pressure.
    r = absorber({"nitrogen": "30 kmol/h"}, {"n-decane": "40 kmol/h", "propane": "4 kmol/h"},
                 "60 degC", "60 degC", "1.5 bar", 6, ["nitrogen", "propane", "n-decane"])
    assert r.ok, r.error
    s = r.streams
    stripped = s["GAS_OUT"]["component_flows"]["propane"] / s["LIQ"]["component_flows"]["propane"]
    assert stripped > 0.9
    comp, energy = balances(r, ["GAS", "LIQ"], ["GAS_OUT", "LIQ_OUT"])
    assert comp < 1e-6 and energy < 1e-6


def two_feed_column(params, top="FEED1", bottom="FEED2"):
    """Benzene-rich and toluene-rich feeds into one column."""
    def feed(uid, flows):
        return {"id": uid, "type": "Feed",
                "params": {"T": "90 degC", "P": "1.5 bar", "flows": flows}}

    return solve({
        "components": ["benzene", "toluene"],
        "property_method": "peng-robinson",
        "units": [feed("A", {"benzene": "35 kmol/h", "toluene": "15 kmol/h"}),
                  feed("B", {"benzene": "15 kmol/h", "toluene": "35 kmol/h"}),
                  {"id": "C", "type": "DistillationColumn", "params": params}],
        "streams": [{"id": top, "from": "A.out", "to": "C.in"},
                    {"id": bottom, "from": "B.out", "to": "C.in"}]
        + [{"id": p, "from": f"C.{p}"} for p in ("distillate", "bottoms", "vapor_distillate")],
    })


TWO_FEEDS = {"stages": 16, "feed_stage": 8, "reflux_ratio": 2, "distillate_rate": "50 kmol/h"}


def test_column_with_two_feeds_on_their_own_stages_balances():
    r = two_feed_column({**TWO_FEEDS, "feed_stages": {"FEED1": 6, "FEED2": 11}})
    assert r.ok, r.error
    comp, energy = balances(r, ["FEED1", "FEED2"], ["distillate", "bottoms"])
    assert comp < 1e-6 and energy < 1e-6
    # Feeding each stream where its composition fits beats mixing them on one stage.
    same = two_feed_column(TWO_FEEDS)
    purity = [res.streams["distillate"]["mole_fractions"]["benzene"] for res in (r, same)]
    assert purity[0] > purity[1]


def test_column_rejects_a_feed_stage_outside_the_column():
    r = two_feed_column({**TWO_FEEDS, "feed_stages": {"FEED2": 30}})
    assert not r.ok
    assert "FEED2" in r.error["message"] and "outside" in r.error["message"]
