"""Costing of the batch reactor (costing.py): vessels sized from the batch, the batch cycle,
heating and cooling per batch, and the batch reactor in the TAC seen by studies."""

import json
import math
from pathlib import Path

import pytest

from fluxsheet import FlowsheetModel, solve
from fluxsheet.costing import (
    BASE_CEPCI,
    BATCH_ASPECT,
    BATCH_FILL,
    FBM_STIRRED,
    FM_VESSEL,
    STIRRED,
    TOTAL_MODULE,
    annualization_factor,
    vessel_pressure_factor,
)
from fluxsheet.optimization import run_optimization
from fluxsheet.sensitivity import run_sensitivity
from fluxsheet.study_models import CostingModel

EXAMPLES = Path(__file__).parent.parent / "examples"
HOUR = 3600.0
ISO = {"id": "ISO", "stoichiometry": {"n-butane": -1, "isobutane": 1},
       "kinetics": {"k0": 1e-3, "Ea": "0 kJ/mol"}}


def batch(params=None, costing=None, *, T="300 K", P="4 bar", flow="1 mol/s", **extra):
    """Feed -> BatchReactor with liquid n-butane isomerizing, 1 h isothermal batches every
    3 h (10.8 kmol, about 1.1 m3 of liquid per batch)."""
    return {
        "components": ["n-butane", "isobutane"], "property_method": "ideal",
        "reactions": [ISO],
        "units": [
            {"id": "F", "type": "Feed",
             "params": {"T": T, "P": P, "flows": {"n-butane": flow}}},
            {"id": "B", "type": "BatchReactor",
             "params": {"reactions": ["ISO"], "phase": "liquid", "batch_time": "1 h",
                        "T": "300 K", "cycle_time": "3 h", **(params or {})}},
        ],
        "streams": [{"id": "IN", "from": "F.out", "to": "B.in"},
                    {"id": "OUT", "from": "B.out"}, {"id": "VENT", "from": "B.vent"}],
        "costing": costing if costing is not None else {},
        **extra,
    }


def costed(data):
    r = solve(FlowsheetModel.model_validate(data))
    assert r.ok, r.error
    return r, r.units["B"], r.costing["units"]["B"]


# ---- sizing ---------------------------------------------------------------------------------

def test_vessel_is_a_jacketed_agitated_reactor_holding_the_batch():
    r, u, c = costed(batch())
    assert 0.8 < u["max_volume"] < 1.5  # some 11 kmol of liquid butane
    (item,) = c["items"]
    V = u["max_volume"] / BATCH_FILL
    assert item["equipment"] == STIRRED.label and item["size"] == pytest.approx(V)
    assert item["parallel"] == 1
    index = r.costing["cepci"] / BASE_CEPCI
    assert item["purchased_cost"] == pytest.approx(STIRRED.cost(V)[0] * index)
    # Carbon steel near ambient pressure: exactly the CSTR's bare-module factor.
    assert vessel_pressure_factor(4e5, item["size"] ** (1 / 3)) == 1.0
    assert item["bare_module_cost"] == pytest.approx(item["purchased_cost"] * FBM_STIRRED)
    assert u["bare_module_cost"] == pytest.approx(item["bare_module_cost"])
    assert "B" not in r.costing["not_costed"]


def test_the_same_volume_costs_as_much_as_a_cstr():
    _, u, c = costed(batch())
    V = c["items"][0]["size"]
    cstr = {
        "components": ["n-butane", "isobutane"], "property_method": "ideal",
        "reactions": [ISO],
        "units": [{"id": "F", "type": "Feed",
                   "params": {"T": "300 K", "P": "4 bar", "flows": {"n-butane": 1.0}}},
                  {"id": "R", "type": "CSTR",
                   "params": {"reactions": ["ISO"], "phase": "liquid", "V": V,
                              "T": "300 K"}}],
        "streams": [{"id": "IN", "from": "F.out", "to": "R.in"},
                    {"id": "OUT", "from": "R.out"}],
        "costing": {},
    }
    r = solve(FlowsheetModel.model_validate(cstr))
    assert r.ok, r.error
    assert r.units["R"]["bare_module_cost"] == pytest.approx(u["bare_module_cost"])


def test_material_and_pressure_scale_the_vessel_like_the_other_vessels():
    _, cs, _ = costed(batch())
    _, ss, _ = costed(batch(costing={"material": "SS"}))
    factor = (2.25 + 1.82 * FM_VESSEL["SS"]) / (2.25 + 1.82)
    assert ss["bare_module_cost"] == pytest.approx(cs["bare_module_cost"] * factor)
    assert ss["purchased_cost"] == pytest.approx(cs["purchased_cost"])

    _, hp, c_hp = costed(batch({"P": "40 bar"}, P="40 bar"))
    V = c_hp["items"][0]["size"]
    D = (4.0 * V / (BATCH_ASPECT * math.pi)) ** (1 / 3)
    fp = vessel_pressure_factor(40e5, D)
    assert fp > 1.0
    assert hp["bare_module_cost"] == pytest.approx(
        hp["purchased_cost"] * FBM_STIRRED * (2.25 + 1.82 * fp) / (2.25 + 1.82))


def test_a_large_charge_needs_vessels_in_parallel():
    # 50 mol/s over 3 h: some 55 m3 of liquid per batch; Turton covers up to 35 m3.
    _, _, c = costed(batch(flow="50 mol/s"))
    item = c["items"][0]
    assert item["size"] > STIRRED.hi
    assert item["parallel"] == math.ceil(item["size"] / STIRRED.hi)


# ---- the batch cycle ------------------------------------------------------------------------

def test_cycle_gives_batches_per_year_and_vessels():
    # Reaction 1 h + turnaround 1 h (the default) fits into the 3 h cycle: one vessel.
    _, u, c = costed(batch())
    cyc = c["batch"]
    assert cyc["batches_per_year"] == pytest.approx(8000 / 3)
    assert u["batches_per_year"] == pytest.approx(8000 / 3)
    assert cyc["turnaround"] == HOUR and cyc["turnaround_source"] == "turnaround"
    assert cyc["occupied_time"] == pytest.approx(2 * HOUR)
    assert cyc["vessels"] == 1 and u["vessels_costed"] == 1.0
    assert any("batches/yr" in n for n in c["notes"])

    # A 1.5 h turnaround and 1 h charging: 3.5 h in the vessel per 3 h cycle -> 2 vessels,
    # each as large; the costing's hours per year set the batches.
    _, u2, c2 = costed(batch({"charge_time": "1 h"},
                             {"batch_turnaround": "1.5 h", "hours_per_year": 6000}))
    assert c2["batch"]["occupied_time"] == pytest.approx(3.5 * HOUR)
    assert c2["batch"]["vessels"] == 2 and c2["items"][0]["parallel"] == 2
    assert u2["bare_module_cost"] == pytest.approx(2 * u["bare_module_cost"])
    assert u2["batches_per_year"] == pytest.approx(2000)


def test_the_units_own_down_time_replaces_the_turnaround():
    _, _, c = costed(batch({"down_time": "2.5 h"}, {"batch_turnaround": "0 h"}))
    assert c["batch"]["turnaround"] == pytest.approx(2.5 * HOUR)
    assert c["batch"]["turnaround_source"] == "down time"
    assert c["batch"]["vessels"] == 2  # 1 h + 2.5 h > 3 h
    _, _, c = costed(batch(costing={"batch_turnaround": 0}))
    assert c["batch"]["turnaround"] == 0.0 and c["batch"]["vessels"] == 1


def test_utility_cost_per_batch_times_batches_is_the_operating_cost():
    _, u, c = costed(batch(costing={"cooling_price": 2.0}))
    (cool,) = c["utilities"]
    assert cool["kind"] == "cooling" and cool["price"] == 2.0
    # The exothermic isomerization held at 300 K: all of the heat is cooling.
    assert u["heat_per_batch"] < 0
    assert u["cooling_per_batch"] == pytest.approx(u["heat_per_batch"])
    assert u["heating_per_batch"] == 0.0
    assert cool["per_batch"] == pytest.approx(u["heat_per_batch"])
    per_batch = abs(u["heat_per_batch"]) * 1e-9 * 2.0
    assert cool["cost_per_batch"] == pytest.approx(per_batch)
    assert u["utility_cost_per_batch"] == pytest.approx(per_batch)
    assert u["operating_cost"] == pytest.approx(per_batch * u["batches_per_year"])
    # With capital: the batch reactor's share of the TAC per batch.
    f = annualization_factor(CostingModel())
    share = u["operating_cost"] + f * TOTAL_MODULE * u["bare_module_cost"]
    assert u["cost_per_batch"] == pytest.approx(share / u["batches_per_year"])


def test_heating_and_cooling_are_priced_apart():
    # Feed at 280 K brought to 300 K (heating at the start), then the reaction's heat is
    # removed (cooling): the net heat would hide both.
    _, u, c = costed(batch(T="280 K", costing={"heating_price": 5.0,
                                                "cooling_price": 1.0}))
    assert u["heating_per_batch"] > 0 > u["cooling_per_batch"]
    assert u["heating_per_batch"] + u["cooling_per_batch"] == pytest.approx(
        u["heat_per_batch"], rel=1e-9, abs=1.0)
    kinds = {x["kind"]: x for x in c["utilities"]}
    assert set(kinds) == {"heating", "cooling"}
    assert kinds["heating"]["per_batch"] == pytest.approx(u["heating_per_batch"])
    assert kinds["cooling"]["per_batch"] == pytest.approx(u["cooling_per_batch"])
    expected = (u["heating_per_batch"] * 5.0 - u["cooling_per_batch"] * 1.0) * 1e-9
    assert u["utility_cost_per_batch"] == pytest.approx(expected)


def test_default_utility_follows_the_batch_temperatures():
    # Held at 300 K, cooling water (30 -> 40 degC) is too warm: chilled water.
    _, _, c = costed(batch())
    assert "chilled water" in c["utilities"][0]["utility"]
    _, _, c = costed(batch({"T": "340 K"}, T="340 K", P="20 bar"))
    cool = next(x for x in c["utilities"] if x["kind"] == "cooling")
    assert "cooling water" in cool["utility"]


def test_bad_turnaround():
    r = solve(FlowsheetModel.model_validate(batch(costing={"batch_turnaround": "1 kg"})))
    assert not r.ok and "turnaround" in r.error["message"]


# ---- the TAC with a batch reactor -------------------------------------------------------------

def test_tac_includes_the_batch_reactor():
    r, u, _ = costed(batch())
    cost = r.costing
    assert cost["bare_module_cost"] == pytest.approx(u["bare_module_cost"])
    assert cost["operating_cost"] == pytest.approx(u["operating_cost"])
    f = cost["annualization_factor"]
    assert cost["tac"] == pytest.approx(u["operating_cost"]
                                        + f * TOTAL_MODULE * u["bare_module_cost"])
    assert cost["tac"] > 0


def test_example_batch_flowsheet_is_costed():
    data = json.loads((EXAMPLES / "ethyl_acetate_batch.json").read_text("utf-8"))
    data["costing"] = {}
    r = solve(FlowsheetModel.model_validate(data))
    assert r.ok, r.error
    c = r.costing["units"]["BATCH"]
    assert c["batch"]["turnaround_source"] == "down time"  # the example's 1 h
    assert c["batch"]["occupied_time"] == pytest.approx(
        0.5 * HOUR + r.units["BATCH"]["batch_time"] + HOUR)
    assert c["batch"]["vessels"] == 1
    assert c["batch"]["batches_per_year"] == pytest.approx(8000 / 3)
    assert r.units["BATCH"]["cost_per_batch"] > r.units["BATCH"]["utility_cost_per_batch"]


def test_studies_read_the_batch_costs():
    data = batch(sensitivity=[{
        "id": "S", "vary": [{"unit": "B", "param": "cycle_time", "lower": "2.5 h",
                             "upper": "4 h", "points": 3}],
        "outputs": [{"property": "tac"}, {"unit": "B", "property": "cost_per_batch"},
                    {"unit": "B", "property": "batches_per_year"},
                    {"unit": "B", "property": "bare_module_cost"}]}])
    table = run_sensitivity(FlowsheetModel.model_validate(data))
    assert table["ok"], table.get("error")
    rows = [r["outputs"] for r in table["rows"]]
    # A longer cycle charges more per batch: a larger vessel, fewer and dearer batches.
    assert rows[0][0] < rows[1][0] < rows[2][0]
    assert rows[0][1] < rows[1][1] < rows[2][1]
    assert [r[2] for r in rows] == pytest.approx([8000 / 2.5, 8000 / 3.25, 8000 / 4])


def test_optimization_minimizes_the_tac_of_a_batch_flowsheet():
    data = batch(optimizations=[{
        "id": "OPT", "variables": {"tac": {"property": "tac"}}, "objective": "tac",
        "vary": [{"unit": "B", "param": "cycle_time", "lower": "2.5 h", "upper": "5 h"}]}])
    data["units"][1]["params"]["cycle_time"] = "4 h"
    out = run_optimization(FlowsheetModel.model_validate(data))
    assert out["ok"], out.get("error")
    # The TAC rises with the batch size (one vessel throughout): the shortest cycle wins.
    (var,) = out["variables"]
    assert var["value"] == pytest.approx(2.5 * HOUR, rel=1e-3) and var["at_bound"]
    assert out["objective"]["value"] < out["objective"]["start"]
    data["units"][1]["params"]["cycle_time"] = var["value"]
    data.pop("optimizations")
    r = solve(FlowsheetModel.model_validate(data))
    assert out["objective"]["value"] == pytest.approx(r.costing["tac"], rel=1e-6)


def test_a_design_spec_holds_the_cost_per_batch():
    # The cost per batch of a 3.25 h cycle, found again from a 3 h start.
    _, u, _ = costed(batch({"cycle_time": "3.25 h"}))
    target = u["cost_per_batch"]
    r = solve(FlowsheetModel.model_validate(batch(design_specs=[
        {"id": "COST", "target": {"unit": "B", "property": "cost_per_batch"},
         "value": target,
         "vary": {"unit": "B", "param": "cycle_time", "lower": "2.5 h", "upper": "4 h"}}])))
    assert r.ok, r.error
    assert r.units["B"]["cost_per_batch"] == pytest.approx(target, rel=1e-4)
    assert r.units["B"]["cycle_time"] == pytest.approx(3.25 * HOUR, rel=1e-3)


def test_a_program_that_heats_and_cools_is_priced_both_ways():
    prog = [{"time": "0 h", "T": "300 K"}, {"time": "0.5 h", "T": "330 K"},
            {"time": "1 h", "T": "300 K"}]
    data = batch({"T_profile": prog}, {"heating_price": 5.0, "cooling_price": 1.0},
                 P="20 bar")
    data["units"][1]["params"].pop("T")
    _, u, c = costed(data)
    assert u["heating_per_batch"] > 0 > u["cooling_per_batch"]
    assert u["heating_per_batch"] + u["cooling_per_batch"] == pytest.approx(
        u["heat_per_batch"], rel=1e-9, abs=1.0)
    assert {x["kind"] for x in c["utilities"]} == {"heating", "cooling"}
