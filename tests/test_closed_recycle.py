"""Steady circulation targets and independent initial guesses."""

from pathlib import Path

import pytest

from fluxsheet.flowsheet import FlowsheetModel
from fluxsheet.solver import solve

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "r744_closed_recycle.json"
RETURN = "GASKUEHLERAUSGANG"


def r744_model():
    return FlowsheetModel.load(EXAMPLE).model_dump(by_alias=True, exclude_none=True)


def cut_reference():
    model = r744_model()
    model["units"] += [
        {"id": "FEED", "type": "Feed", "params": {
            "T": "36.5 degC", "P": "90.9 bar", "flows": {"carbon dioxide": "2.272236676 mol/s"},
        }},
        {"id": "PRODUCT", "type": "Product"},
    ]
    # Open the loop at the gas cooler. The specified feed and the calculated return
    # must agree independently of the closed solver's tear and flow constraint.
    model["streams"][-1:] = [
        {"id": "CUT_IN", "from": "FEED.out", "to": "HP_VENTIL.in"},
        {"id": RETURN, "from": "GASKUEHLER.out", "to": "PRODUCT.in"},
    ]
    return solve(model)


@pytest.mark.parametrize("acceleration", ["none", "wegstein", "broyden"])
def test_closed_r744_matches_balanced_cut_with_an_independent_guess(acceleration):
    model = r744_model()
    model["solver"]["acceleration"] = acceleration
    model["streams"][-1]["recycle_guess"] = {
        "T": "35.5 degC", "P": "88 bar", "flows": {"carbon dioxide": "2.1 mol/s"},
    }
    actual = solve(FlowsheetModel.model_validate_json(FlowsheetModel.model_validate(model).to_json()))
    reference = cut_reference()
    assert actual.ok, actual.error
    assert reference.ok, reference.error
    assert actual.tear_streams == [RETURN]
    assert actual.iterations > 1
    assert actual.streams[RETURN]["mass_flow"] == pytest.approx(0.10, rel=1e-8)
    for sid, stream in actual.streams.items():
        for key in ("T", "P", "mass_flow", "vapor_fraction", "enthalpy_flow"):
            assert stream[key] == pytest.approx(reference.streams[sid][key], rel=1e-7, abs=1e-5)
    for uid, key in (("VERDAMPFER", "duty"), ("GASKUEHLER", "duty"), ("VERDICHTER", "power")):
        assert actual.units[uid][key] == pytest.approx(reference.units[uid][key], rel=1e-7)
    assert actual.balance_status == reference.balance_status == "ok"


def test_r744_can_start_without_an_explicit_guess():
    model = r744_model()
    del model["streams"][-1]["recycle_guess"]
    result = solve(model)
    assert result.ok, result.error
    assert result.streams[RETURN]["mass_flow"] == pytest.approx(0.10)
    assert result.balance_status == "ok"


@pytest.mark.parametrize("unrelated_feed", [False, True])
def test_guess_alone_cannot_set_circulation_even_with_an_unrelated_feed(unrelated_feed):
    model = r744_model()
    del model["streams"][-1]["circulation"]
    model["streams"][-1]["recycle_guess"]["flows"] = {"carbon dioxide": 2.1}
    if unrelated_feed:
        model["units"] += [
            {"id": "FEED", "type": "Feed", "params": {
                "T": "25 degC", "P": "1 bar", "flows": {"carbon dioxide": 1},
            }},
            {"id": "PRODUCT", "type": "Product"},
        ]
        model["streams"].append({"id": "FEED_OUT", "from": "FEED.out", "to": "PRODUCT.in"})
    result = solve(model)
    assert not result.ok
    assert result.error["where"] == RETURN
    assert "circulation mass-flow specification" in result.error["message"]
    assert "circulation.mass_flow" in result.error["hint"]


def test_circulation_can_choose_another_internal_stream_as_tear():
    model = r744_model()
    spec = model["streams"][-1].pop("circulation")
    del model["streams"][-1]["recycle_guess"]
    discharge = next(s for s in model["streams"] if s["id"] == "DRUCKGAS")
    discharge["circulation"] = spec
    discharge["recycle_guess"] = {"T": "115 degC", "P": "90.9 bar"}
    result = solve(model)
    assert result.ok, result.error
    assert result.tear_streams == ["DRUCKGAS"]
    assert result.streams[RETURN]["mass_flow"] == pytest.approx(0.10)
    assert result.balance_status == "ok"


def test_a_feed_on_the_other_exchanger_side_does_not_define_circulation():
    model = {
        "components": ["water"],
        "property_method": "ideal",
        "units": [
            {"id": "FEED", "type": "Feed", "params": {
                "T": "20 degC", "P": "1 bar", "flows": {"water": 10},
            }},
            {"id": "HX", "type": "HeatExchanger", "params": {"T_cold_out": "30 degC"}},
            {"id": "HEATER", "type": "Heater", "params": {"T": "60 degC", "P": "1 bar"}},
            {"id": "PRODUCT", "type": "Product"},
        ],
        "streams": [
            {"id": "COLD_IN", "from": "FEED.out", "to": "HX.cold_in"},
            {"id": "COLD_OUT", "from": "HX.cold_out", "to": "PRODUCT.in"},
            {"id": "RETURN", "from": "HX.hot_out", "to": "HEATER.in"},
            {"id": "HOT", "from": "HEATER.out", "to": "HX.hot_in"},
        ],
    }
    result = solve(model)
    assert not result.ok
    assert result.error["where"] == "HOT"
    assert "circulation mass-flow specification" in result.error["message"]
    model["streams"][-1]["circulation"] = {"mass_flow": 0.5}
    result = solve(model)
    assert result.ok, result.error
    assert result.streams["COLD_OUT"]["molar_flow"] == pytest.approx(10.0)
    assert result.streams["RETURN"]["mass_flow"] == pytest.approx(0.5)
    assert result.units["HEATER"]["duty"] == pytest.approx(result.units["HX"]["duty"])
    assert result.balance_status == "ok"


def test_specification_cannot_make_up_a_net_material_loss():
    result = solve({
        "components": ["water"],
        "property_method": "ideal",
        "units": [
            {"id": "SPLIT", "type": "Splitter", "params": {"split": 0.8}},
            {"id": "HEATER", "type": "Heater", "params": {"T": "25 degC", "P": "1 bar"}},
            {"id": "LOSS", "type": "Product"},
        ],
        "streams": [
            {"id": "RETURN", "from": "SPLIT.out1", "to": "HEATER.in",
             "circulation": {"mass_flow": "360 kg/h"}},
            {"id": "HEATED", "from": "HEATER.out", "to": "SPLIT.in"},
            {"id": "LOSS", "from": "SPLIT.out2", "to": "LOSS.in"},
        ],
    })
    assert not result.ok
    assert result.error["where"] == "RETURN"
    assert "target 0.1 kg/s, returned 0.08 kg/s" in result.error["message"]


def test_feedless_mixture_needs_a_composition_specification():
    model = {
        "components": ["water", "ethanol"],
        "property_method": "ideal",
        "units": [
            {"id": "H", "type": "Heater", "params": {"T": "35 degC", "P": "1 bar"}},
            {"id": "C", "type": "Heater", "params": {"T": "20 degC"}},
        ],
        "streams": [
            {"id": "WARM", "from": "H.out", "to": "C.in"},
            {"id": "RETURN", "from": "C.out", "to": "H.in", "circulation": {"mass_flow": 0.10},
             "recycle_guess": {"flows": {"water": 1, "ethanol": 1}}},
        ],
    }
    result = solve(model)
    assert not result.ok
    assert "circulating composition specification" in result.error["message"]
    model["streams"][-1]["circulation"]["mole_fractions"] = {"water": 0.75, "ethanol": 0.25}
    result = solve(model)
    assert result.ok, result.error
    for stream in result.streams.values():
        assert stream["mass_flow"] == pytest.approx(0.10)
        assert stream["mole_fractions"] == pytest.approx({"water": 0.75, "ethanol": 0.25})
    assert sum(unit["duty"] for unit in result.units.values()) == pytest.approx(0.0, abs=1e-4)
    assert result.balance_status == "ok"


@pytest.mark.parametrize(("first", "second", "missing", "anchor"), [
    ({"duty": "100 W"}, {"duty": "-100 W", "P": "1 bar"},
     "temperature (or enthalpy)", {"T": "40 degC"}),
    ({"T": "40 degC"}, {"T": "20 degC"}, "absolute pressure", {"P": "1 bar"}),
])
def test_an_initial_guess_cannot_anchor_a_closed_recycle_state(first, second, missing, anchor):
    model = {
        "components": ["water"],
        "property_method": "ideal",
        "units": [
            {"id": "H", "type": "Heater", "params": first},
            {"id": "C", "type": "Heater", "params": second},
        ],
        "streams": [
            {"id": "HOT", "from": "H.out", "to": "C.in"},
            {"id": "RETURN", "from": "C.out", "to": "H.in", "circulation": {"mass_flow": 0.1},
             "recycle_guess": {"T": "25 degC", "P": "3 bar"}},
        ],
    }
    result = solve(model)
    assert not result.ok
    assert result.error["where"] == "RETURN"
    assert missing in result.error["message"]
    assert "recycle_guess" in result.error["hint"]

    model["units"][0]["params"] = {k: v for k, v in first.items() if k != "duty"} | anchor
    result = solve(model)
    assert result.ok, result.error
    assert result.streams["RETURN"]["mass_flow"] == pytest.approx(0.1)
    assert result.balance_status == "ok"


def test_thermal_return_accepts_a_guess_but_not_an_independent_mass_flow():
    model = {
        "components": ["water"],
        "property_method": "ideal",
        "units": [
            {"id": "FEED", "type": "Feed", "params": {"T": "80 degC", "P": "1 bar",
                                                     "flows": {"water": 10}}},
            {"id": "HX", "type": "HeatExchanger", "params": {"T_cold_out": "30 degC"}},
            {"id": "COOLER", "type": "Heater", "params": {"T": "20 degC"}},
            {"id": "PRODUCT", "type": "Product"},
        ],
        "streams": [
            {"id": "FEED_IN", "from": "FEED.out", "to": "HX.hot_in"},
            {"id": "HOT_OUT", "from": "HX.hot_out", "to": "COOLER.in"},
            {"id": "COLD_IN", "from": "COOLER.out", "to": "HX.cold_in",
             "recycle_guess": {"T": "20 degC", "P": "1 bar", "flows": {"water": 10}}},
            {"id": "PRODUCT_OUT", "from": "HX.cold_out", "to": "PRODUCT.in"},
        ],
    }
    result = solve(model)
    assert result.ok, result.error
    assert result.streams["PRODUCT_OUT"]["molar_flow"] == pytest.approx(10)
    model["streams"][2]["circulation"] = {"mass_flow": result.streams["COLD_IN"]["mass_flow"]}
    result = solve(model)
    assert not result.ok
    assert result.error["where"] == "COLD_IN"
    assert "material recycle" in result.error["message"]
