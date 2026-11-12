"""Pump and compressor performance curves, and the logic blocks Mult, Dupl and Selector.

Checked against hand calculations: the pump head is rho g H, its power rho g H Q / eta, the
affinity laws scale flow by N and head by N^2, NPSH available follows from the steam table,
the isentropic curve reproduces the pressure-specified compressor and the polytropic head of
argon (an ideal gas with constant cp) has a closed form.
"""

import itertools
import math

import pytest

from fluxsheet import solve
from fluxsheet.curves import G, parse_curve
from fluxsheet.dynamics import simulate
from fluxsheet.quantities import UnitError
from fluxsheet.units import mixers

WATER = {"T": "25 degC", "P": "1 bar", "flows": {"water": "2000 kmol/h"}}  # about 36 m3/h
ARGON = {"T": "300 K", "P": "1 bar", "flows": {"argon": "100 kmol/h"}}  # about 2490 m3/h


def machine(unit_type, params, feed=WATER, components=("water",), method="nrtl"):
    return solve({
        "components": list(components), "property_method": method,
        "units": [{"id": "F", "type": "Feed", "params": feed},
                  {"id": "U", "type": unit_type, "params": params},
                  {"id": "PR", "type": "Product"}],
        "streams": [{"id": "S1", "from": "F.out", "to": "U.in"},
                    {"id": "S2", "from": "U.out", "to": "PR.in"}],
    })


def ok(r):
    assert r.ok, r.error
    return r


# --- curves ------------------------------------------------------------------------------

def test_two_point_curve_is_linear_and_extrapolates_without_negative_head():
    c = parse_curve({"points": [[0, 50], [100, 30]]})
    assert c.head(50 / 3600) == pytest.approx(40 * G)
    assert c.head(120 / 3600) == pytest.approx(26 * G)  # straight on beyond the last point
    assert c.head(1000 / 3600) == 0.0  # but never below zero
    assert c.efficiency(0.01) is None


def test_pchip_curve_passes_its_points_without_overshoot():
    pts = [[0, 50], [10, 50], [20, 49], [40, 40], [60, 20]]
    c = parse_curve({"points": pts})
    for q, h in pts:
        assert c.head(q / 3600) == pytest.approx(h * G)
    heads = [c.head(q / 3600) / G for q in range(61)]
    assert max(heads) <= 50 + 1e-9  # the flat start does not bulge above the points
    assert all(b <= a + 1e-9 for a, b in itertools.pairwise(heads))  # monotone like the data


def test_curve_units_and_efficiency_column():
    c = parse_curve({"flow_unit": "l/s", "head_unit": "ft", "efficiency_unit": "%",
                     "points": [[0, 100, 0], [10, 90, 60], [20, 70, 75]]})
    assert c.flows[1] == pytest.approx(0.01)
    assert c.heads[1] == pytest.approx(90 * 0.3048 * G)
    assert c.eff_flows == (0.01, 0.02)  # the zero at shut-off is dropped
    assert c.efficiency(0.015) == pytest.approx(0.675)
    gas = parse_curve({"head_unit": "kJ/kg", "points": [[1000, 120], [2000, 100]]})
    assert gas.heads == (120e3, 100e3)


@pytest.mark.parametrize("curve, message", [
    ({"points": [[10, 50]]}, "at least 2 points"),
    ({"points": []}, "at least 2 points"),
    ({"points": [[0, 50], [20, 48], [15, 40]]}, "must increase"),
    ({"points": [[0, 50], [0, 48]]}, "must increase"),
    ({"points": [[0, 50], [20, 48, 0.7]]}, "every point needs 2"),
    ({"points": [[0, 50, 0.5], [20, 48, 75]]}, "between 0 and 1"),
    ({"points": [[0, 50], [20, -1]]}, "head is negative"),
    ({"points": [[0, 50], [20, "x"]]}, "not a number"),
    ({"head_unit": "psi", "points": [[0, 50], [20, 40]]}, "unknown head unit"),
    ({"flow_unit": "bbl/d", "points": [[0, 50], [20, 40]]}, "unknown flow unit"),
    ({"points": [[0, 50], [20, 40]], "colour": "red"}, "unknown curve entries"),
    ([[0, 50], [20, 40]], "expected a table"),
])
def test_bad_curves_are_explained(curve, message):
    with pytest.raises(UnitError, match=message):
        parse_curve(curve)
    r = machine("Pump", {"curve": curve})
    assert not r.ok and message in r.error["message"]
    assert r.error["message"].startswith("Performance curve")


def test_speed_needs_rated_speed_and_a_curve():
    r = machine("Pump", {"curve": {"points": [[0, 50], [60, 30]]}, "speed": 1450})
    assert not r.ok and "rated speed" in r.error["message"]
    r = machine("Pump", {"dP": "2 bar", "speed": 1450, "rated_speed": 1450})
    assert not r.ok and "needs a performance curve" in r.error["message"]
    r = machine("Pump", {"dP": "2 bar", "curve": {"points": [[0, 50], [60, 30]]}})
    assert not r.ok and "over-specified" in r.error["message"]


# --- pump --------------------------------------------------------------------------------

PUMP_CURVE = {"points": [[0, 60, 0], [20, 57, 0.55], [40, 48, 0.72], [60, 33, 0.68]]}


def test_pump_on_its_curve_gives_rho_g_h_and_power():
    r = ok(machine("Pump", {"curve": PUMP_CURVE}))
    u, s1, s2 = r.units["U"], r.streams["S1"], r.streams["S2"]
    curve = parse_curve(PUMP_CURVE)
    Q = u["volumetric_flow"]
    rho = s1["mass_flow"] / Q
    assert rho == pytest.approx(997, rel=0.01)  # water at 25 °C
    head = curve.head(Q)  # J/kg
    eta = curve.efficiency(Q)
    assert u["head"] == pytest.approx(head / G, rel=1e-9)
    assert s2["P"] - s1["P"] == pytest.approx(rho * head, rel=1e-9)
    assert u["power"] == pytest.approx(rho * head * Q / eta, rel=1e-9)
    assert u["efficiency"] == pytest.approx(eta)
    assert 45 < u["head"] < 52 and 0.68 < eta < 0.72  # Q is about 36 m3/h
    # Same power as a pressure-specified pump to that pressure at that efficiency.
    same = ok(machine("Pump", {"P": s2["P"], "efficiency": eta})).units["U"]
    assert same["power"] == pytest.approx(u["power"], rel=1e-6)
    assert r.streams["S2"]["T"] > s1["T"]  # the losses warm the water


def test_affinity_laws_scale_flow_by_speed_and_head_by_its_square():
    rated = ok(machine("Pump", {"curve": PUMP_CURVE}))
    Q = rated.units["U"]["volumetric_flow"]
    curve = parse_curve(PUMP_CURVE)
    r = 0.8
    slow = ok(machine("Pump", {"curve": PUMP_CURVE, "speed": 1160, "rated_speed": 1450}))
    u = slow.units["U"]
    assert u["head"] == pytest.approx(curve.head(Q / r) * r**2 / G, rel=1e-9)
    assert u["efficiency"] == pytest.approx(curve.efficiency(Q / r), rel=1e-9)
    # The operating point Q / r on the rated curve, scaled: Q/r at rated speed gives H; at
    # speed r the flow Q gives r^2 H.
    flow = {"water": f"{2000 / r} kmol/h"}
    at_rated = ok(machine("Pump", {"curve": PUMP_CURVE}, feed={**WATER, "flows": flow}))
    assert u["head"] == pytest.approx(at_rated.units["U"]["head"] * r**2, rel=1e-6)


def test_pump_without_flow_reaches_shut_off_head():
    r = ok(solve({  # the splitter sends nothing to the pump
        "components": ["water"], "property_method": "nrtl",
        "units": [{"id": "F", "type": "Feed", "params": WATER},
                  {"id": "SP", "type": "Splitter", "params": {"split": 1.0}},
                  {"id": "U", "type": "Pump", "params": {"curve": PUMP_CURVE}},
                  {"id": "PR", "type": "Product"}, {"id": "PR2", "type": "Product"}],
        "streams": [{"id": "S0", "from": "F.out", "to": "SP.in"},
                    {"id": "A", "from": "SP.out1", "to": "PR2.in"},
                    {"id": "S1", "from": "SP.out2", "to": "U.in"},
                    {"id": "S2", "from": "U.out", "to": "PR.in"}],
    }))
    assert r.streams["S1"]["molar_flow"] == 0.0
    assert r.units["U"]["power"] == 0.0
    dP = r.streams["S2"]["P"] - r.streams["S1"]["P"]
    assert dP / (997 * G) == pytest.approx(60, rel=0.01)


def test_pump_beyond_its_curve_warns_and_extrapolates():
    r = ok(machine("Pump", {"curve": {"points": [[0, 50], [10, 45], [20, 38]]}}))
    assert any("outside the pump curve" in w for w in r.warnings)
    assert r.units["U"]["head"] < 38
    r = machine("Pump", {"curve": {"points": [[0, 50], [10, 30], [20, 5]]}})
    assert any("no head left" in w for w in r.warnings)


def test_npsh_available_matches_the_steam_table():
    hot = {"T": "80 degC", "P": "1 bar", "flows": {"water": "2000 kmol/h"}}
    r = ok(machine("Pump", {"dP": "3 bar"}, feed=hot))
    # Steam table, 80 °C: p_sat = 47.39 kPa, rho = 971.8 kg/m3.
    expected = (1e5 - 47.39e3) / (971.8 * G)
    assert r.units["U"]["NPSH_available"] == pytest.approx(expected, rel=0.02)
    assert not any("NPSH" in w for w in r.warnings)
    low = ok(machine("Pump", {"dP": "3 bar", "NPSH_required": "4 m", "suction_height": "-2 m"},
                     feed=hot))
    assert low.units["U"]["NPSH_available"] == pytest.approx(
        r.units["U"]["NPSH_available"] - 2.0, rel=1e-9)
    assert sum("below NPSH required" in w for w in low.warnings) == 1


def test_warnings_are_reported_once_after_a_recycle_converges():
    hot = {"T": "80 degC", "P": "1 bar", "flows": {"water": "100 kmol/h"}}
    r = ok(solve({
        "components": ["water"], "property_method": "nrtl",
        "units": [{"id": "F", "type": "Feed", "params": hot},
                  {"id": "M", "type": "Mixer", "params": {"P": "1 bar"}},
                  {"id": "P1", "type": "Pump",
                   "params": {"curve": PUMP_CURVE, "NPSH_required": "8 m"}},
                  {"id": "V", "type": "Valve", "params": {"P": "1 bar"}},
                  {"id": "SP", "type": "Splitter", "params": {"split": 0.5}},
                  {"id": "PR", "type": "Product"}],
        "streams": [{"id": "S1", "from": "F.out", "to": "M.in"},
                    {"id": "S2", "from": "M.out", "to": "P1.in"},
                    {"id": "S3", "from": "P1.out", "to": "V.in"},
                    {"id": "S4", "from": "V.out", "to": "SP.in"},
                    {"id": "S5", "from": "SP.out1", "to": "M.in"},
                    {"id": "S6", "from": "SP.out2", "to": "PR.in"}],
    }))
    assert r.iterations > 1
    assert sum("P1: NPSH available" in w for w in r.warnings) == 1


# --- compressor --------------------------------------------------------------------------

GAS_CURVE = {"head_unit": "kJ/kg",
             "points": [[1500, 125, 0.70], [2500, 115, 0.80], [3500, 85, 0.74]]}


def test_isentropic_curve_reproduces_the_pressure_specified_compressor():
    r = ok(machine("Compressor", {"curve": GAS_CURVE, "curve_type": "isentropic"},
                   feed=ARGON, components=("argon",), method="peng-robinson"))
    u = r.units["U"]
    curve = parse_curve(GAS_CURVE)
    assert u["head"] == pytest.approx(curve.head(u["volumetric_flow"]))
    spec = ok(machine("Compressor", {"P": r.streams["S2"]["P"], "efficiency": u["efficiency"]},
                      feed=ARGON, components=("argon",), method="peng-robinson"))
    assert spec.units["U"]["power"] == pytest.approx(u["power"], rel=1e-6)
    assert spec.streams["S2"]["T"] == pytest.approx(r.streams["S2"]["T"], abs=1e-3)
    mass = r.streams["S1"]["mass_flow"]
    assert u["power"] == pytest.approx(mass * u["head"] / u["efficiency"], rel=1e-9)
    assert not r.warnings


def test_polytropic_head_of_argon_matches_the_ideal_gas_formula():
    r = ok(machine("Compressor", {"curve": GAS_CURVE}, feed=ARGON, components=("argon",),
                   method="ideal"))
    u = r.units["U"]
    k, M, T1 = 5.0 / 3.0, 39.948e-3, 300.0
    m = (k - 1.0) / k / u["efficiency"]  # (n - 1) / n
    ratio = (1.0 + u["head"] * M * m / (8.314462618 * T1)) ** (1.0 / m)
    assert u["pressure_ratio"] == pytest.approx(ratio, rel=2e-3)
    # The discharge temperature of the polytropic path: T2 = T1 ratio^((n-1)/n).
    assert r.streams["S2"]["T"] == pytest.approx(T1 * ratio**m, rel=2e-3)


def test_compressor_warns_about_surge_and_stonewall():
    def run(flow, **extra):
        feed = {**ARGON, "flows": {"argon": flow}}
        return ok(machine("Compressor", {"curve": GAS_CURVE, **extra}, feed=feed,
                          components=("argon",), method="ideal"))
    assert any("surge" in w for w in run("50 kmol/h").warnings)  # about 1250 m3/h
    assert not run("50 kmol/h", surge_flow="1000 m3/h").warnings
    assert any("stonewall" in w for w in run("150 kmol/h").warnings)  # about 3740 m3/h
    assert not run("100 kmol/h").warnings
    # At 120 % speed the surge line moves from 1500 to 1800 m3/h.
    assert not run("70 kmol/h").warnings  # about 1750 m3/h
    assert any("surge" in w for w in run("70 kmol/h", speed=1.2, rated_speed=1).warnings)


def test_compressor_curve_refuses_liquid():
    r = machine("Compressor", {"curve": GAS_CURVE}, feed=WATER)
    assert not r.ok and "liquid" in r.error["message"]


# --- logic blocks ------------------------------------------------------------------------

def logic(units, streams, components=("water", "methanol")):
    feeds = [{"id": "F1", "type": "Feed", "params": {
        "T": "30 degC", "P": "2 bar", "flows": {"water": "10 kmol/h", "methanol": "5 kmol/h"}}},
        {"id": "F2", "type": "Feed", "params": {
            "T": "50 degC", "P": "3 bar", "flows": {"water": "1 kmol/h"}}}]
    return solve({"components": list(components), "property_method": "nrtl",
                  "units": feeds + units, "streams": streams})


def balance_warnings(r):
    return [w for w in r.warnings if "mass balance" in w]


def test_mult_scales_flows_at_the_same_state():
    r = ok(logic([{"id": "M", "type": "Mult", "params": {"factor": 2.5}},
                  {"id": "P", "type": "Product"}, {"id": "P2", "type": "Product"}],
                 [{"id": "S1", "from": "F1.out", "to": "M.in"},
                  {"id": "S2", "from": "M.out", "to": "P.in"},
                  {"id": "S3", "from": "F2.out", "to": "P2.in"}]))
    a, b = r.streams["S1"], r.streams["S2"]
    assert b["molar_flow"] == pytest.approx(2.5 * a["molar_flow"])
    assert (b["T"], b["P"]) == (a["T"], a["P"])
    assert b["mole_fractions"] == pytest.approx(a["mole_fractions"], rel=1e-12)
    assert not balance_warnings(r)


def test_dupl_copies_and_drops_unconnected_copies_silently():
    r = ok(logic([{"id": "D", "type": "Dupl"}, {"id": "H", "type": "Heater",
                                                  "params": {"T": "60 degC"}},
                  {"id": "P", "type": "Product"}, {"id": "P2", "type": "Product"},
                  {"id": "P3", "type": "Product"}],
                 [{"id": "S1", "from": "F1.out", "to": "D.in"},
                  {"id": "C1", "from": "D.out1", "to": "H.in"},
                  {"id": "C2", "from": "D.out2", "to": "P2.in"},
                  {"id": "S3", "from": "H.out", "to": "P.in"},
                  {"id": "S4", "from": "F2.out", "to": "P3.in"}]))
    for c in ("C1", "C2"):
        assert r.streams[c]["component_flows"] == r.streams["S1"]["component_flows"]
    assert not [k for k in r.streams if k.startswith("D.")]  # no phantom out3, out4
    assert not r.warnings


def test_selector_passes_the_selected_stream_by_id():
    units = [{"id": "SEL", "type": "Selector", "params": {"selected": "S2"}},
             {"id": "P", "type": "Product"}]
    streams = [{"id": "S1", "from": "F1.out", "to": "SEL.in"},
               {"id": "S2", "from": "F2.out", "to": "SEL.in"},
               {"id": "S3", "from": "SEL.out", "to": "P.in"}]
    r = ok(logic(units, streams))
    assert r.streams["S3"]["component_flows"] == r.streams["S2"]["component_flows"]
    assert not r.warnings
    # Delete S2 from the flowsheet: the selection does not silently fall to another stream.
    r = logic(units, [streams[0], streams[2]])
    assert not r.ok and "'S2' does not enter the selector" in r.error["message"]
    assert "S1" in r.error["hint"]
    r = logic([{"id": "SEL", "type": "Selector"}, units[1]], streams)
    assert not r.ok and "under-specified" in r.error["message"]


def test_balance_check_still_catches_a_real_leak(monkeypatch, expected_balance_exceptions):
    for scope in ("SP", "flowsheet"):
        for quantity, component in (("mass", None), ("energy", None),
                                    ("component", "water"), ("component", "methanol")):
            expected_balance_exceptions[(scope, quantity, component, "warning")] = (
                "Intentional leakage in a splitter downstream of a multiplier")
    def leaky(self, pkg, inlets):
        out = original(self, pkg, inlets)
        out["out2"].flows = [0.9 * f for f in out["out2"].flows]
        return out
    original = mixers.Splitter.calculate
    units = [{"id": "M", "type": "Mult", "params": {"factor": 3}},
             {"id": "SP", "type": "Splitter"},
             {"id": "P", "type": "Product"}, {"id": "P2", "type": "Product"},
             {"id": "P3", "type": "Product"}]
    streams = [{"id": "S1", "from": "F1.out", "to": "M.in"},
               {"id": "S2", "from": "M.out", "to": "SP.in"},
               {"id": "S3", "from": "SP.out1", "to": "P.in"},
               {"id": "S4", "from": "SP.out2", "to": "P2.in"},
               {"id": "S5", "from": "F2.out", "to": "P3.in"}]
    assert not balance_warnings(ok(logic(units, streams)))
    monkeypatch.setattr(mixers.Splitter, "calculate", leaky)
    assert balance_warnings(ok(logic(units, streams)))


# --- dynamics ----------------------------------------------------------------------------

def test_curve_pump_and_logic_blocks_stay_algebraic_in_dynamics():
    model = {
        "components": ["water"], "property_method": "nrtl",
        "units": [
            {"id": "F", "type": "Feed", "params": {"T": "25 degC", "P": "1 bar",
                                                   "flows": {"water": "1000 kmol/h"}}},
            {"id": "PU", "type": "Pump", "params": {"curve": PUMP_CURVE, "speed": 1450,
                                                    "rated_speed": 1450}},
            {"id": "X2", "type": "Mult", "params": {"factor": 2}},
            {"id": "TANK", "type": "MixingTank", "params": {"P": "1.2 bar"}},
            {"id": "PROD", "type": "Product"}],
        "streams": [{"id": "S1", "from": "F.out", "to": "PU.in"},
                    {"id": "S2", "from": "PU.out", "to": "X2.in"},
                    {"id": "S3", "from": "X2.out", "to": "TANK.in"},
                    {"id": "S4", "from": "TANK.liquid", "to": "PROD.in"}],
        "dynamics": {
            "t_end": "20 min", "points": 11,
            "units": {"TANK": {"diameter": "2 m", "height": "3 m"}},
            "controllers": [{
                "id": "PC", "measure": {"stream": "S2", "property": "P"},
                "manipulate": {"unit": "PU", "param": "speed"}, "action": "reverse",
                "gain": 1, "integral_time": "1 min", "pv_range": ["1 bar", "9 bar"],
                "op_range": [0, 3000]}],
            "events": [{"time": "1 min", "controller": "PC", "value": "5 bar"}]},
    }
    out = simulate(model)
    assert out["ok"], out.get("error")
    s = {x["key"]: x["values"] for x in out["series"]}
    assert s["controller.PC.pv"][-1] == pytest.approx(5e5, rel=2e-3)
    assert s["stream.S3.molar_flow"][-1] == pytest.approx(2 * s["stream.S1.molar_flow"][-1])
    # The speed that gives 5 bar, from the affinity laws on the steady curve.
    steady = ok(solve(model))
    head = (5e5 - 1e5) / (steady.streams["S1"]["mass_flow"] / steady.units["PU"]
                          ["volumetric_flow"] * G)
    Q = steady.units["PU"]["volumetric_flow"]
    curve = parse_curve(PUMP_CURVE)
    ratio = s["controller.PC.op"][-1] / 1450
    assert curve.head(Q / ratio) * ratio**2 / G == pytest.approx(head, rel=5e-3)
    assert not math.isnan(ratio)


# --- review follow-ups -------------------------------------------------------------------

def recycle_pump(params):
    return solve({
        "components": ["water"], "property_method": "nrtl",
        "units": [{"id": "F", "type": "Feed", "params": {
                      "T": "25 degC", "P": "1 bar", "flows": {"water": "100 kmol/h"}}},
                  {"id": "M", "type": "Mixer", "params": {"P": "1 bar"}},
                  {"id": "P1", "type": "Pump", "params": params},
                  {"id": "V", "type": "Valve", "params": {"P": "1 bar"}},
                  {"id": "SP", "type": "Splitter", "params": {"split": 0.5}},
                  {"id": "PR", "type": "Product"}],
        "streams": [{"id": "S1", "from": "F.out", "to": "M.in"},
                    {"id": "S2", "from": "M.out", "to": "P1.in"},
                    {"id": "S3", "from": "P1.out", "to": "V.in"},
                    {"id": "S4", "from": "V.out", "to": "SP.in"},
                    {"id": "S5", "from": "SP.out1", "to": "M.in"},
                    {"id": "S6", "from": "SP.out2", "to": "PR.in"}],
    })


@pytest.mark.parametrize("params", [{"dP": "3 bar"}, {"curve": PUMP_CURVE}])
def test_pump_results_are_those_of_the_converged_recycle(params):
    r = ok(recycle_pump(params))
    assert r.iterations > 1
    s2, u = r.streams["S2"], r.units["P1"]
    assert s2["molar_flow"] == pytest.approx(2 * 100 / 3.6, rel=1e-4)  # doubled by the loop
    rho = 997.0
    assert u["volumetric_flow"] == pytest.approx(s2["mass_flow"] / rho, rel=0.01)
    assert u["head"] == pytest.approx((r.streams["S3"]["P"] - s2["P"]) / (
        s2["mass_flow"] / u["volumetric_flow"] * G), rel=1e-9)


def test_curve_compressor_with_no_flow_passes_through():
    r = ok(solve({  # the splitter sends nothing to the compressor; its inlet state is liquid
        "components": ["water"], "property_method": "nrtl",
        "units": [{"id": "F", "type": "Feed", "params": WATER},
                  {"id": "SP", "type": "Splitter", "params": {"split": 1.0}},
                  {"id": "U", "type": "Compressor", "params": {"curve": GAS_CURVE}},
                  {"id": "PR", "type": "Product"}, {"id": "PR2", "type": "Product"}],
        "streams": [{"id": "S0", "from": "F.out", "to": "SP.in"},
                    {"id": "A", "from": "SP.out1", "to": "PR2.in"},
                    {"id": "S1", "from": "SP.out2", "to": "U.in"},
                    {"id": "S2", "from": "U.out", "to": "PR.in"}],
    }))
    assert r.units["U"]["power"] == 0.0
    assert r.streams["S2"]["molar_flow"] == 0.0


def test_two_phase_pump_inlet_moves_its_liquid_on_the_curve():
    feed = {"P": "1 bar", "VF": 0.02, "flows": {"water": "2000 kmol/h"}}
    r = ok(machine("Pump", {"curve": PUMP_CURVE}, feed=feed))
    u, s1 = r.units["U"], r.streams["S1"]
    assert 0.0 < s1["vapor_fraction"] < 0.05
    liquid_mass = s1["mass_flow"] * (1.0 - s1["vapor_fraction"])  # one component
    rho = liquid_mass / u["volumetric_flow"]
    assert rho == pytest.approx(958, rel=0.02)  # saturated water at 100 °C
    assert u["power"] == pytest.approx(liquid_mass * u["head"] * G / u["efficiency"], rel=1e-6)
    assert any("cavitate" in w for w in r.warnings)


@pytest.mark.parametrize("key", ["flow_unit", "head_unit", "efficiency_unit"])
def test_curve_unit_names_must_be_text(key):
    with pytest.raises(UnitError, match=f"{key} must be a unit name"):
        parse_curve({key: 5, "points": [[0, 50], [20, 40]]})


def test_polytropic_compressor_follows_the_fan_laws():
    def run(flow, **extra):
        feed = {**ARGON, "flows": {"argon": flow}}
        return ok(machine("Compressor", {"curve": GAS_CURVE, **extra}, feed=feed,
                          components=("argon",), method="ideal")).units["U"]
    r = 1.1
    fast = run(f"{100 * r} kmol/h", speed=r, rated_speed=1.0)
    rated = run("100 kmol/h")
    assert fast["head"] == pytest.approx(rated["head"] * r**2, rel=1e-6)
    assert fast["efficiency"] == pytest.approx(rated["efficiency"], rel=1e-6)
    assert fast["pressure_ratio"] > rated["pressure_ratio"]


def test_dupl_with_one_copy_connected_warns_once():
    r = ok(logic([{"id": "D", "type": "Dupl"}, {"id": "P", "type": "Product"},
                  {"id": "P2", "type": "Product"}],
                 [{"id": "S1", "from": "F1.out", "to": "D.in"},
                  {"id": "C1", "from": "D.out1", "to": "P.in"},
                  {"id": "S4", "from": "F2.out", "to": "P2.in"}]))
    assert r.warnings == ["D: outlet 'out2' is not connected"]
