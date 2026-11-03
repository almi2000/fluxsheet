"""Optimization around the flowsheet (optimization.py) and equipment costing (costing.py)."""

import json
import math
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from fluxsheet import FlowsheetModel, limits, solve
from fluxsheet.costing import (
    COMPRESSOR,
    EXCHANGER,
    VESSEL,
    annualization_factor,
    exchanger_pressure_factor,
    vessel_pressure_factor,
)
from fluxsheet.optimization import run_optimization
from fluxsheet.sensitivity import run_sensitivity
from fluxsheet.server import app
from fluxsheet.study_models import CostingModel

EXAMPLES = Path(__file__).parent.parent / "examples"


def water_heater(heater=None, feed_T="300 K", **extra):
    """Feed -> Heater -> Product with liquid water, 3.6 kmol/h (1 mol/s)."""
    return {
        "components": ["water"], "property_method": "ideal",
        "units": [
            {"id": "FEED", "type": "Feed",
             "params": {"T": feed_T, "P": "2 bar", "flows": {"water": 1.0}}},
            {"id": "H1", "type": "Heater", "params": heater or {"T": "330 K"}},
            {"id": "OUT", "type": "Product"},
        ],
        "streams": [{"id": "S1", "from": "FEED.out", "to": "H1.in"},
                    {"id": "S2", "from": "H1.out", "to": "OUT.in"}],
        **extra,
    }


def optimize(model: dict, **opt) -> dict:
    model = {**model, "optimizations": [{"id": "OPT", **opt}]}
    return run_optimization(FlowsheetModel.model_validate(model))


# ---- analytic optima ------------------------------------------------------------------------

FEED_T = {"unit": "FEED", "param": "T"}
FEED_P = {"unit": "FEED", "param": "P"}


def test_quadratic_of_two_parameters_finds_the_analytic_minimum():
    out = optimize(water_heater(), variables={"T": FEED_T, "P": FEED_P},
                   objective="(T - 312.5)**2 + 4 * (P / 1e5 - 2.6)**2",
                   vary=[{**FEED_T, "lower": 290, "upper": 320},
                         {**FEED_P, "lower": "1.5 bar", "upper": "4 bar"}])
    assert out["ok"], out.get("error")
    assert out["converged"] and out["method"] == "SLSQP"
    T, P = (v["value"] for v in out["variables"])
    assert T == pytest.approx(312.5, abs=0.05)  # ftol 1e-5 on f/f(start)
    assert P == pytest.approx(2.6e5, rel=2e-3)
    assert out["objective"]["value"] == pytest.approx(0.0, abs=1e-3)
    # The history starts at the model's values and the objective falls to the optimum.
    hist = out["history"]
    assert hist[0]["x"] == pytest.approx([300.0, 2e5])
    assert hist[-1]["objective"] <= hist[0]["objective"]
    assert out["apply"][0] == {"unit": "FEED", "param": "T", "component": None,
                               "value": pytest.approx(T)}


def test_maximize_turns_the_objective_around():
    out = optimize(water_heater(), variables={"T": FEED_T}, sense="maximize",
                   objective="-(T - 305)**2", vary=[{**FEED_T, "lower": 290, "upper": 320}])
    assert out["ok"]
    assert out["variables"][0]["value"] == pytest.approx(305.0, abs=0.02)


def test_an_active_inequality_constraint_holds_the_optimum_on_its_limit():
    out = optimize(water_heater(), variables={"T": FEED_T}, objective="(T - 300)**2",
                   vary=[{**FEED_T, "lower": 290, "upper": 330}],
                   constraints=[{"expression": "T", "relation": ">=", "value": 310}])
    assert out["ok"]
    assert out["variables"][0]["value"] == pytest.approx(310.0, abs=1e-3)
    (c,) = out["constraints"]
    assert c["active"] and c["satisfied"]
    assert c["achieved"] == pytest.approx(310.0, abs=1e-3)
    assert not out["variables"][0]["at_bound"]


def test_an_inactive_constraint_is_reported_as_such():
    out = optimize(water_heater(), variables={"T": FEED_T}, objective="(T - 300)**2",
                   vary=[{**FEED_T, "lower": 290, "upper": 330}],
                   constraints=[{"expression": "T", "relation": "<=", "value": 320}])
    assert out["variables"][0]["value"] == pytest.approx(300.0, abs=0.02)
    assert out["constraints"][0]["satisfied"] and not out["constraints"][0]["active"]


def test_equality_constraint_on_a_line():
    # min T^2 + P^2 (scaled) on T + P = 2: the point (1, 1).
    out = optimize(water_heater(), variables={"T": FEED_T, "P": FEED_P},
                   objective="(T / 300)**2 + (P / 3e5)**2",
                   vary=[{**FEED_T, "lower": 250, "upper": 400},
                         {**FEED_P, "lower": "1 bar", "upper": "5 bar"}],
                   constraints=[{"expression": "T / 300 + P / 3e5", "relation": "=",
                                 "value": 2.0}])
    assert out["ok"], out.get("error")
    T, P = (v["value"] for v in out["variables"])
    assert T == pytest.approx(300.0, rel=1e-3)
    assert P == pytest.approx(3e5, rel=1e-3)


def test_a_result_objective_with_a_design_spec_inside():
    # The heater duty follows the feed temperature through the design spec on its outlet:
    # Q = n cp (T_out - T_feed). Choose the feed temperature that makes the duty 1 kW.
    model = water_heater({"duty": "2 kW"}, design_specs=[
        {"id": "TOUT", "target": {"stream": "S2", "property": "T"}, "value": "330 K",
         "vary": {"unit": "H1", "param": "duty", "lower": "-10 kW", "upper": "20 kW"}}])
    out = optimize(model, variables={"Q": {"unit": "H1", "property": "duty"}},
                   objective="(Q / 1e3 - 1)**2",
                   vary=[{**FEED_T, "lower": 290, "upper": 329}])
    assert out["ok"], out.get("error")
    res = out["result"]
    assert res["units"]["H1"]["duty"] == pytest.approx(1e3, rel=1e-3)
    assert res["streams"]["S2"]["T"] == pytest.approx(330.0, abs=1e-4)
    assert 300.0 < out["variables"][0]["value"] < 329.0


# ---- the costed example against a sweep ------------------------------------------------------

def test_reactor_optimum_matches_a_sensitivity_sweep():
    data = json.loads((EXAMPLES / "ethyl_acetate_optimization.json").read_text("utf-8"))
    model = FlowsheetModel.model_validate(data)
    out = run_optimization(model)
    assert out["ok"] and out["converged"], out.get("message")
    V, T = (v["value"] for v in out["variables"])
    assert T == pytest.approx(358.15, abs=0.01)  # the T <= 85 degC constraint is active
    assert out["constraints"][0]["active"]
    assert out["objective"]["value"] > out["objective"]["start"]

    # Sweep the volume at the optimal temperature: the margin peaks where the optimizer is.
    for u in data["units"]:
        if u["id"] == "CSTR":
            u["params"]["T"] = T
    data["sensitivity"] = [{"id": "V", "vary": [{"unit": "CSTR", "param": "V", "lower": 2,
                                                 "upper": 9, "points": 15}],
                            "outputs": [{"stream": "S3", "property": "component_flow",
                                         "component": "ethyl acetate"}, {"property": "tac"}]}]
    table = run_sensitivity(FlowsheetModel.model_validate(data))
    margins = [0.25 * r["outputs"][0] * 0.08811 * 3600 * 8000 - r["outputs"][1]
               for r in table["rows"]]
    k = max(range(len(margins)), key=margins.__getitem__)
    assert 0 < k < len(margins) - 1  # an interior optimum
    assert abs(V - table["rows"][k]["inputs"][0]) <= 0.5  # within one grid step
    assert out["objective"]["value"] >= max(margins) - 1e-6 * abs(max(margins))


# ---- failures -------------------------------------------------------------------------------

def test_unknown_name_in_the_objective():
    out = optimize(water_heater(), variables={"T": FEED_T}, objective="T + oops",
                   vary=[{**FEED_T, "lower": 290, "upper": 320}])
    assert not out["ok"]
    assert "oops" in out["error"]["message"] and out["error"]["where"] == "OPT"


def test_no_optimization_and_unknown_id():
    model = FlowsheetModel.model_validate(water_heater())
    assert "no optimization" in run_optimization(model)["error"]["message"]
    model = FlowsheetModel.model_validate({**water_heater(), "optimizations": [
        {"id": "A", "objective": "1", "vary": [{**FEED_T, "lower": 290, "upper": 320}]}]})
    assert "B" in run_optimization(model, "B")["error"]["message"]


def test_varying_a_parameter_a_design_spec_adjusts_is_refused():
    model = water_heater({"duty": "2 kW"}, design_specs=[
        {"id": "TOUT", "target": {"stream": "S2", "property": "T"}, "value": "330 K",
         "vary": {"unit": "H1", "param": "duty", "lower": "-10 kW", "upper": "20 kW"}}])
    out = optimize(model, objective="1",
                   vary=[{"unit": "H1", "param": "duty", "lower": 0, "upper": 5000}])
    assert not out["ok"] and "design spec TOUT" in out["error"]["message"]


def test_infeasible_constraints_are_reported_with_the_history():
    out = optimize(water_heater(), variables={"T": FEED_T}, objective="T",
                   vary=[{**FEED_T, "lower": 290, "upper": 320}],
                   constraints=[{"expression": "T", "relation": ">=", "value": 400}])
    assert not out["ok"]
    assert "no feasible point" in out["error"]["message"]
    assert out["history"] and out["variables"][0]["value"] == pytest.approx(320.0)


def pump(P_lower, P_upper, start):
    return {
        "components": ["water"], "property_method": "ideal",
        "units": [
            {"id": "FEED", "type": "Feed",
             "params": {"T": "300 K", "P": "2 bar", "flows": {"water": 10.0}}},
            {"id": "PUMP", "type": "Pump", "params": {"P": start}},
            {"id": "OUT", "type": "Product"},
        ],
        "streams": [{"id": "S1", "from": "FEED.out", "to": "PUMP.in"},
                    {"id": "S2", "from": "PUMP.out", "to": "OUT.in"}],
        "optimizations": [{
            "id": "OPT", "variables": {"W": {"unit": "PUMP", "property": "power"},
                                       "P": {"unit": "PUMP", "param": "P"}},
            # Power falls with the discharge pressure down to the inlet pressure; below it
            # the pump cannot run.
            "objective": "W",
            "vary": [{"unit": "PUMP", "param": "P", "lower": P_lower, "upper": P_upper}]}],
    }


def test_failed_trials_make_the_search_back_off_with_cobyla():
    out = run_optimization(FlowsheetModel.model_validate(pump("1.2 bar", "5 bar", "4 bar")))
    assert out["ok"], out.get("error")
    assert out["failed_evaluations"] > 0
    assert "COBYLA" in out["method"]
    # The best feasible point is the inlet pressure, where the pump still runs.
    assert out["variables"][0]["value"] == pytest.approx(2e5, rel=2e-3)


def test_a_start_where_the_flowsheet_fails_everywhere():
    out = run_optimization(FlowsheetModel.model_validate(pump("1.2 bar", "1.8 bar", "1.5 bar")))
    assert not out["ok"]
    assert "fails at the start point" in out["error"]["message"]
    assert "below inlet pressure" in out["error"]["message"]


def test_demo_time_limit_stops_the_run(monkeypatch):
    monkeypatch.setattr(limits, "MAX_SECONDS", 0.0)
    with limits.calculation(enabled=True):
        out = optimize(water_heater(), variables={"T": FEED_T}, objective="(T - 305)**2",
                       vary=[{**FEED_T, "lower": 290, "upper": 320}])
    assert not out["ok"] and "longer than" in out["error"]["message"]


def test_demo_caps_the_iterations():
    model = FlowsheetModel.model_validate({**water_heater(), "optimizations": [
        {"id": "A", "objective": "1", "max_iterations": 200,
         "vary": [{**FEED_T, "lower": 290, "upper": 320}]}]})
    with pytest.raises(Exception, match="iterations"):
        limits.check_model_size(model, enabled=True)


def test_optimize_endpoint():
    data = {**water_heater(), "optimizations": [
        {"id": "Q", "variables": {"T": FEED_T}, "objective": "(T - 305)**2",
         "vary": [{**FEED_T, "lower": 290, "upper": 320}]}]}
    body = TestClient(app).post("/api/optimize", json={"model": data, "optimization": "Q"}).json()
    assert body["ok"], body.get("error")
    assert body["variables"][0]["value"] == pytest.approx(305.0, abs=0.02)
    assert body["result"]["ok"]


# ---- costing ------------------------------------------------------------------------------------

def test_correlation_values_by_hand():
    # Turton: vertical vessel, 10 m3: log10 Cp = 3.4974 + 0.4485 + 0.1074.
    assert VESSEL.cost(10.0) == (pytest.approx(10 ** 4.0533), 1)
    # Above the range: identical units in parallel.
    cost, n = EXCHANGER.cost(2500.0)
    assert n == 3 and cost == pytest.approx(3 * EXCHANGER.cost(2500.0 / 3)[0])
    # Below the range: the six-tenths rule, continuous at the lower limit.
    assert COMPRESSOR.cost(450.0 * (1 - 1e-9))[0] == pytest.approx(COMPRESSOR.cost(450.0)[0])
    assert COMPRESSOR.cost(45.0)[0] == pytest.approx(COMPRESSOR.cost(450.0)[0] * 0.1 ** 0.6)


def test_pressure_factors():
    assert vessel_pressure_factor(101325.0, 2.0) == 1.0  # ambient: the minimum wall
    # 20 barg, D 2 m: ((21 * 2) / (2 (850 - 12.6)) + 0.00315) / 0.0063
    P = 101325.0 + 20e5
    assert vessel_pressure_factor(P, 2.0) == pytest.approx((42 / 1674.8 + 0.00315) / 0.0063)
    assert vessel_pressure_factor(0.3e5, 1.0) == 1.25  # vacuum
    assert exchanger_pressure_factor(4e5) == 1.0
    lg = math.log10(30.0)
    assert exchanger_pressure_factor(101325.0 + 30e5) == pytest.approx(
        10 ** (0.03881 - 0.11272 * lg + 0.08183 * lg * lg))


def test_annualization_factor():
    assert annualization_factor(CostingModel()) == pytest.approx(0.16275, abs=1e-5)
    assert annualization_factor(CostingModel(interest_rate=0.0, plant_life=5)) == 0.2
    assert annualization_factor(CostingModel(annualization_factor=0.25)) == 0.25


def costed(model: dict, **costing) -> dict:
    result = solve(FlowsheetModel.model_validate({**model, "costing": costing}))
    assert result.ok, result.error
    return result


def test_heater_utility_and_exchanger():
    # 1 mol/s of water from 300 to 330 K needs about 2.26 kW of LP steam.
    r = costed(water_heater())
    unit = r.costing["units"]["H1"]
    (use,) = unit["utilities"]
    Q = r.units["H1"]["duty"]
    assert use["utility"].startswith("LP steam") and use["price"] == 4.54
    assert use["cost"] == pytest.approx(Q * 8000 * 3600 * 4.54e-9)
    (hx,) = unit["items"]
    # Countercurrent against condensing steam at 160 degC, U = 500 W/m2K.
    lmtd = ((433.15 - 300) - (433.15 - 330)) / math.log((433.15 - 300) / (433.15 - 330))
    assert hx["size"] == pytest.approx(Q / (500 * lmtd), rel=1e-6)
    assert r.units["H1"]["bare_module_cost"] == pytest.approx(hx["bare_module_cost"])
    c = r.costing
    assert c["capital_cost"] == pytest.approx(1.18 * c["bare_module_cost"])
    assert c["tac"] == pytest.approx(c["operating_cost"] + c["annualization_factor"]
                                     * c["capital_cost"])


def test_utilities_follow_the_temperature():
    hot = costed(water_heater({"T": "200 degC"}))
    assert hot.costing["units"]["H1"]["utilities"][0]["utility"].startswith("HP steam")
    cool = costed(water_heater({"T": "45 degC"}, feed_T="90 degC"))
    assert cool.costing["units"]["H1"]["utilities"][0]["utility"].startswith("cooling water")
    cold = costed(water_heater({"T": "10 degC"}, feed_T="40 degC"))
    assert cold.costing["units"]["H1"]["utilities"][0]["utility"].startswith("refrigerant")
    own = costed(water_heater({"T": "45 degC"}, feed_T="90 degC"), cooling_price=1.0)
    assert own.costing["units"]["H1"]["utilities"][0]["price"] == 1.0


def test_costs_scale_with_the_cost_index_and_material():
    base = costed(water_heater())
    double = costed(water_heater(), cepci=2 * 397.0)
    assert double.costing["purchased_cost"] == pytest.approx(
        2 * base.costing["purchased_cost"] * 397.0 / 800.0)
    ss = costed(water_heater(), material="SS")
    # B1 + B2 FM: 1.63 + 1.66 * 2.73 against 1.63 + 1.66
    assert ss.costing["bare_module_cost"] / base.costing["bare_module_cost"] == pytest.approx(
        (1.63 + 1.66 * 2.73) / (1.63 + 1.66))


def test_column_size_and_costs_are_plausible():
    data = json.loads((EXAMPLES / "benzene_toluene_column.json").read_text("utf-8"))
    data.pop("design_specs", None)
    r = costed(data)
    cols = [u for u in r.costing["units"].values() if u["type"] == "DistillationColumn"]
    assert cols
    items = {i["item"]: i for i in cols[0]["items"]}
    assert {"shell", "trays", "condenser", "reboiler"} <= set(items)
    D = math.sqrt(4 * items["trays"]["size"] / math.pi)
    assert 0.5 < D < 2.5  # m, for about 1 MW of reboiler duty at 1 atm
    kinds = {u["item"]: u["utility"] for u in cols[0]["utilities"]}
    assert kinds["condenser"].startswith("cooling water")
    assert kinds["reboiler"].startswith("LP steam")


def test_tac_as_a_study_output_and_its_errors():
    data = {**water_heater(), "costing": {}, "sensitivity": [
        {"id": "S", "vary": [{"unit": "H1", "param": "T", "lower": 310, "upper": 350,
                              "points": 3}],
         "outputs": [{"property": "tac"}, {"unit": "H1", "property": "operating_cost"}]}]}
    table = run_sensitivity(FlowsheetModel.model_validate(data))
    tac = [r["outputs"][0] for r in table["rows"]]
    assert tac[0] < tac[1] < tac[2]  # more duty, more steam and area
    data.pop("costing")
    table = run_sensitivity(FlowsheetModel.model_validate(data))
    assert not table["ok"] and "needs costing" in table["error"]["message"]
    with pytest.raises(ValueError, match="costing total"):
        FlowsheetModel.model_validate({**data, "sensitivity": [
            {"id": "S", "vary": data["sensitivity"][0]["vary"],
             "outputs": [{"property": "T"}]}]})


def test_costing_names_unknown_units():
    result = solve(FlowsheetModel.model_validate(
        {**water_heater(), "costing": {"units": {"NOPE": {"material": "SS"}}}}))
    assert not result.ok and "NOPE" in result.error["message"]
