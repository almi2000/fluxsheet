"""Two-sided heat exchanger, design specs and sensitivity studies."""

import itertools
import json
import math
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from fluxsheet import FlowsheetModel, solve
from fluxsheet.sensitivity import run_sensitivity, to_csv
from fluxsheet.server import app

EXAMPLES = Path(__file__).parent.parent / "examples"


# --- heat exchanger ---------------------------------------------------------------------

def exchanger(spec, arrangement="countercurrent", hot=("200 degC", {"nitrogen": "100 kmol/h"}),
              cold=("30 degC", {"nitrogen": "150 kmol/h"}), components=("nitrogen",),
              method="ideal", P="2 bar", P_cold=None):
    return {
        "components": list(components),
        "property_method": method,
        "units": [
            {"id": "HOT", "type": "Feed",
             "params": {"T": hot[0], "P": P, "flows": hot[1]}},
            {"id": "COLD", "type": "Feed",
             "params": {"T": cold[0], "P": P_cold or P, "flows": cold[1]}},
            {"id": "X", "type": "HeatExchanger", "params": {**spec, "arrangement": arrangement}},
        ],
        "streams": [
            {"id": "H1", "from": "HOT.out", "to": "X.hot_in"},
            {"id": "C1", "from": "COLD.out", "to": "X.cold_in"},
            {"id": "H2", "from": "X.hot_out"},
            {"id": "C2", "from": "X.cold_out"},
        ],
    }


def terminal_differences(s, countercurrent=True):
    if countercurrent:
        return s["H1"]["T"] - s["C2"]["T"], s["H2"]["T"] - s["C1"]["T"]
    return s["H1"]["T"] - s["C1"]["T"], s["H2"]["T"] - s["C2"]["T"]


def test_countercurrent_exchanger_matches_lmtd_hand_calculation():
    r = solve(exchanger({"T_hot_out": "80 degC"}))
    assert r.ok, r.error
    s, x = r.streams, r.units["X"]
    released = s["H1"]["enthalpy_flow"] - s["H2"]["enthalpy_flow"]
    absorbed = s["C2"]["enthalpy_flow"] - s["C1"]["enthalpy_flow"]
    assert released == pytest.approx(absorbed, rel=1e-6)
    assert x["duty"] == pytest.approx(released, rel=1e-6)
    assert s["H2"]["T"] == pytest.approx(353.15, abs=1e-6)

    # Hand calculation: nitrogen's cp is nearly constant, so the classic LMTD holds.
    a, b = terminal_differences(s)
    lmtd = (a - b) / math.log(a / b)
    assert x["LMTD"] == pytest.approx(lmtd, rel=5e-3)
    assert x["UA"] == pytest.approx(x["duty"] / lmtd, rel=5e-3)
    assert x["min_approach"] == pytest.approx(min(a, b), abs=0.05)
    # Cold side heated by Q / (F cp): 150 kmol/h of nitrogen takes 2/3 of the hot side's drop.
    assert s["C2"]["T"] - s["C1"]["T"] == pytest.approx(120.0 * 100 / 150, rel=0.01)


def test_ua_rating_reproduces_the_design_duty():
    design = solve(exchanger({"min_approach": "15 K"})).units["X"]
    assert design["min_approach"] == pytest.approx(15.0, abs=1e-6)
    rating = solve(exchanger({"UA": design["UA"]})).units["X"]
    assert rating["duty"] == pytest.approx(design["duty"], rel=1e-5)


def test_cocurrent_outlets_approach_but_never_cross():
    r = solve(exchanger({"min_approach": "5 K"}, arrangement="cocurrent"))
    assert r.ok, r.error
    s = r.streams
    assert s["H2"]["T"] - s["C2"]["T"] == pytest.approx(5.0, abs=1e-3)
    counter = solve(exchanger({"min_approach": "5 K"})).units["X"]
    assert counter["duty"] > r.units["X"]["duty"]  # countercurrent recovers more heat


def test_temperature_cross_raises_a_clear_error():
    # Cooling the hot side to 20 degC needs a cold side below 20 degC; it enters at 30 degC.
    r = solve(exchanger({"T_hot_out": "20 degC"}))
    assert not r.ok
    assert r.error["where"] == "X"
    assert "temperature cross" in r.error["message"]
    assert "cannot be cooled below 30.0 °C" in r.error["hint"]

    r = solve(exchanger({"duty": "5000 kW"}, arrangement="cocurrent"))
    assert "temperature cross" in r.error["message"]
    assert "countercurrent" in r.error["hint"]


def test_infeasible_exchanger_specs_explain_themselves():
    r = solve(exchanger({"min_approach": "200 K"}))
    assert "less than the minimum approach" in r.error["message"]
    r = solve(exchanger({"T_hot_out": "250 degC"}))
    assert "not below the hot inlet" in r.error["message"]
    r = solve(exchanger({"duty": "10 kW"}, hot=("20 degC", {"nitrogen": "100 kmol/h"})))
    assert "swap the streams" in r.error["hint"]


def test_condensing_hot_side_has_an_internal_pinch():
    """A vapor that starts to condense pinches inside the exchanger, not at its ends."""
    r = solve(exchanger(
        {"min_approach": "5 K"}, components=("benzene", "toluene"), method="peng-robinson", P="1 atm",
        P_cold="3 bar",
        hot=("115 degC", {"benzene": "50 kmol/h", "toluene": "50 kmol/h"}),
        cold=("25 degC", {"benzene": "80 kmol/h", "toluene": "20 kmol/h"}),
    ))
    assert r.ok, r.error
    x = r.units["X"]
    assert x["min_approach"] == pytest.approx(5.0, abs=1e-3)
    assert min(terminal_differences(r.streams)) > 6.0
    assert r.streams["H1"]["vapor_fraction"] == 1.0
    assert r.streams["H2"]["vapor_fraction"] < 1.0  # partly condensed
    profile = x["profiles"]["temperatures"]
    hot, cold = (s["values"] for s in profile["series"])
    assert min(h - c for h, c in zip(hot, cold)) == pytest.approx(5.0, abs=0.2)


def test_feed_effluent_example_closes_the_loop():
    r = solve(FlowsheetModel.load(EXAMPLES / "feed_effluent_exchanger.json"))
    assert r.ok, r.error
    assert r.tear_streams == ["S4"]
    fehe = r.units["FEHE"]
    assert fehe["min_approach"] == pytest.approx(20.0, abs=1e-3)
    assert fehe["duty"] > 10 * abs(r.units["TRIM"]["duty"])  # most of the preheat is recovered


# --- design specs -----------------------------------------------------------------------

def heater_with_spec(value="60 degC", upper="2000 kW", **spec):
    return {
        "components": ["water", "ethanol"],
        "property_method": "unifac",
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "T": "25 degC", "P": "1 atm",
                "flows": {"water": "80 kmol/h", "ethanol": "20 kmol/h"}}},
            {"id": "H1", "type": "Heater", "params": {"T": "30 degC"}},
        ],
        "streams": [{"id": "S1", "from": "F.out", "to": "H1.in"},
                    {"id": "S2", "from": "H1.out"}],
        "design_specs": [{
            "id": "TOUT", "target": {"stream": "S2", "property": "T"}, "value": value,
            "vary": {"unit": "H1", "param": "duty", "lower": "0 kW", "upper": upper}, **spec,
        }],
    }


def test_design_spec_on_heater_hits_its_target():
    r = solve(heater_with_spec())
    assert r.ok, r.error
    assert r.streams["S2"]["T"] == pytest.approx(333.15, abs=333.15e-5)
    (spec,) = r.design_specs
    assert spec["converged"] and spec["achieved"] == pytest.approx(333.15, abs=333.15e-5)
    # The varied duty replaced the heater's own T specification, and the unit reports it.
    assert spec["vary_value"] == pytest.approx(r.units["H1"]["duty"], rel=1e-9)
    assert 2 <= spec["iterations"] <= 10


def test_design_spec_tolerance_is_a_temperature_difference():
    r = solve(heater_with_spec(tolerance="0.5 degC"))
    assert abs(r.streams["S2"]["T"] - 333.15) <= 0.5


def test_design_spec_that_cannot_be_bracketed_says_so():
    r = solve(heater_with_spec(value="95 degC", upper="200 kW"))
    assert not r.ok
    assert r.error["where"] == "TOUT"
    assert "cannot reach 95 °C" in r.error["message"]
    assert "do not bracket" in r.error["hint"]


@pytest.mark.parametrize("change, message", [
    ({"target": {"stream": "S9", "property": "T"}}, "no stream called 'S9'"),
    ({"target": {"stream": "S2", "property": "colour"}}, "no property 'colour'"),
    ({"target": {"stream": "S2", "property": "mole_fraction"}}, "no component given"),
    ({"target": {"unit": "H1", "property": "power"}}, "reports no result 'power'"),
    ({"vary": {"unit": "H1", "param": "speed", "lower": 0, "upper": 1}}, "no parameter 'speed'"),
    ({"vary": {"unit": "H1", "param": "duty", "lower": 5, "upper": 1}}, "not below the upper"),
])
def test_design_spec_references_are_checked(change, message):
    model = heater_with_spec()
    model["design_specs"][0].update(change)
    r = solve(model)
    assert not r.ok and r.error["where"] == "TOUT"
    assert message in r.error["message"]


def test_design_spec_restarts_in_the_middle_when_the_start_fails():
    # 190 kW crosses temperatures (at most ~137 kW fit); the search restarts at 100 kW.
    model = exchanger({"duty": "190 kW"})
    model["design_specs"] = [{
        "id": "COLD", "target": {"stream": "C2", "property": "T"}, "value": "80 degC",
        "vary": {"unit": "X", "param": "duty", "lower": "0 kW", "upper": "200 kW"},
    }]
    r = solve(model)
    assert r.ok, r.error
    assert r.streams["C2"]["T"] == pytest.approx(353.15, abs=353.15e-5)


def test_column_purity_spec_converges():
    r = solve(FlowsheetModel.load(EXAMPLES / "benzene_toluene_purity.json"))
    assert r.ok, r.error
    (spec,) = r.design_specs
    assert r.streams["DIST"]["mole_fractions"]["benzene"] == pytest.approx(0.995, abs=1e-5)
    assert 1.2 < spec["vary_value"] < 5.0
    assert r.units["COLUMN"]["reflux_ratio"] == pytest.approx(spec["vary_value"], rel=1e-6)


def test_inactive_design_spec_is_ignored():
    model = heater_with_spec()
    model["design_specs"][0]["active"] = False
    r = solve(model)
    assert r.ok and r.design_specs == []
    assert r.streams["S2"]["T"] == pytest.approx(303.15)


def test_studies_survive_a_file_round_trip(tmp_path):
    model = FlowsheetModel.load(EXAMPLES / "benzene_toluene_purity.json")
    model.dump(tmp_path / "copy.json")
    data = json.loads((tmp_path / "copy.json").read_text(encoding="utf-8"))
    assert data["design_specs"][0]["vary"]["param"] == "reflux_ratio"
    assert data["sensitivity"][0]["outputs"][1] == {"unit": "COLUMN",
                                                    "property": "reboiler_duty"}
    assert FlowsheetModel.model_validate(data) == model


# --- sensitivity ------------------------------------------------------------------------

def heater_study(points=6):
    model = heater_with_spec()
    del model["design_specs"]
    model["sensitivity"] = [{
        "id": "DUTY",
        "vary": [{"unit": "H1", "param": "duty", "lower": "0 kW", "upper": "300 kW",
                  "points": points}],
        "outputs": [{"stream": "S2", "property": "T"},
                    {"stream": "S2", "property": "vapor_fraction"}],
    }]
    return model


def test_sensitivity_outlet_temperature_rises_with_duty():
    table = run_sensitivity(FlowsheetModel.model_validate(heater_study()))
    assert table["ok"], table
    assert [c["label"] for c in table["outputs"]] == ["S2 temperature", "S2 vapor fraction"]
    duties = [row["inputs"][0] for row in table["rows"]]
    T = [row["outputs"][0] for row in table["rows"]]
    assert duties == pytest.approx([0, 60e3, 120e3, 180e3, 240e3, 300e3])
    assert T[0] == pytest.approx(298.15, abs=1e-6)
    assert all(b > a for a, b in itertools.pairwise(T))
    csv = to_csv(table).splitlines()
    assert csv[0] == "H1 duty [W],S2 temperature [K],S2 vapor fraction,error"
    assert len(csv) == 7


def test_sensitivity_column_purity_rises_with_reflux():
    model = FlowsheetModel.load(EXAMPLES / "benzene_toluene_purity.json")
    model.design_specs = None  # the spec would hold purity constant
    model.sensitivity[0].vary[0].points = 4
    table = run_sensitivity(model)
    purity = [row["outputs"][0] for row in table["rows"]]
    duty = [row["outputs"][1] for row in table["rows"]]
    assert all(b > a for a, b in itertools.pairwise(purity))
    assert all(b > a for a, b in itertools.pairwise(duty))


def test_two_parameter_sensitivity_covers_the_grid():
    model = heater_study(points=3)
    model["units"][0]["params"]["T"] = "25 degC"
    model["sensitivity"][0]["vary"].append(
        {"unit": "F", "param": "flows", "component": "water", "lower": "40 kmol/h",
         "upper": "80 kmol/h", "points": 2})
    table = run_sensitivity(FlowsheetModel.model_validate(model))
    assert len(table["rows"]) == 6
    assert table["inputs"][1]["label"] == "F component flows of water"
    # At equal duty, less feed means more of it evaporates.
    by_point = {tuple(r["inputs"]): r["outputs"][1] for r in table["rows"]}
    assert by_point[(150e3, 40 / 3.6)] > by_point[(150e3, 80 / 3.6)]


def test_sensitivity_api_runs_a_study_and_reports_setup_errors():
    client = TestClient(app)
    body = client.post("/api/sensitivity", json={"model": heater_study(points=3)}).json()
    assert body["ok"] and len(body["rows"]) == 3
    body = client.post("/api/sensitivity",
                       json={"model": heater_study(), "study": "NOPE"}).json()
    assert not body["ok"]
    assert "choose one of: DUTY" in body["error"]["hint"]
