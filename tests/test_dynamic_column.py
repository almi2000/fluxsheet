"""Dynamic distillation column (dynamics/column.py): start, balances, control."""

from __future__ import annotations

import copy

import numpy as np
import pytest

from fluxsheet import solve
from fluxsheet.dynamics import simulate
from fluxsheet.dynamics.simulator import Simulator
from fluxsheet.flowsheet import FlowsheetModel

LEVELS = [
    {"id": "LC_DRUM", "measure": {"unit": "COL", "property": "drum_level"},
     "manipulate": {"unit": "COL", "param": "distillate_flow"}, "action": "direct",
     "gain": 2, "integral_time": "10 min", "pv_range": [0, 1],
     "op_range": ["0 kmol/h", "150 kmol/h"]},
    {"id": "LC_SUMP", "measure": {"unit": "COL", "property": "sump_level"},
     "manipulate": {"unit": "COL", "param": "bottoms_flow"}, "action": "direct",
     "gain": 2, "integral_time": "10 min", "pv_range": [0, 1],
     "op_range": ["0 kmol/h", "150 kmol/h"]},
]


def column(dynamics, column_params=None, benzene="50 kmol/h", toluene="50 kmol/h"):
    params = {"stages": 15, "feed_stage": 8, "P": "1 atm", "condenser": "total",
              "reboiler": "kettle", "reflux_ratio": 2.5, "distillate_rate": "50 kmol/h"}
    params.update(column_params or {})
    return {
        "components": ["benzene", "toluene"], "property_method": "peng-robinson",
        "units": [
            {"id": "F", "type": "Feed", "params": {"T": "90 degC", "P": "1.5 bar",
             "flows": {"benzene": benzene, "toluene": toluene}}},
            {"id": "COL", "type": "DistillationColumn", "params": params},
            {"id": "D", "type": "Product"}, {"id": "B", "type": "Product"}],
        "streams": [{"id": "S1", "from": "F.out", "to": "COL.in"},
                    {"id": "DIST", "from": "COL.distillate", "to": "D.in"},
                    {"id": "BOTT", "from": "COL.bottoms", "to": "B.in"}],
        "dynamics": dynamics,
    }


def series(out):
    assert out["ok"], out.get("error")
    return {x["key"]: x["values"] for x in out["series"]}


def start_rates(model):
    sim = Simulator(FlowsheetModel.model_validate(model))
    assert sim.solve().ok
    sim._setup()
    return sim, sim._rhs(0.0, sim.x0 / sim.xscale)


@pytest.mark.parametrize("params", [{}, {"condenser": "partial", "reflux_ratio": 2.5,
                                          "distillate_rate": "50 kmol/h"}])
def test_column_starts_at_its_steady_state(params):
    sim, rates = start_rates(column({"t_end": "1 min", "units": {"COL": {}}}, params))
    m = sim.dynamic["COL"]
    # Scaled states change by less than 1e-6 per second: the start is the steady column.
    assert np.abs(rates).max() < 1e-6
    assert m.results["drum_level"] == pytest.approx(0.5)
    assert m.results["sump_level"] == pytest.approx(0.5)
    assert m.s["diameter"] > 0.5  # sized from the vapor load
    # Tridiagonal stages: the Jacobian needs a few groups, not one per state.
    assert len(sim._column_groups()) <= 3 * m.n


def test_undisturbed_column_stays_put():
    s = series(simulate(column({"t_end": "30 min", "points": 4, "units": {"COL": {}}})))
    for key in ("stream.DIST.x.benzene", "stream.BOTT.x.toluene", "unit.COL.control_T",
                "unit.COL.drum_level"):
        assert s[key][-1] == pytest.approx(s[key][0], rel=1e-6, abs=1e-8)


def test_feed_step_settles_on_the_steady_column_with_the_same_flows():
    model = column({"t_end": "8 h", "points": 5, "units": {"COL": {}}, "controllers": LEVELS,
                    "events": [{"time": "5 min", "unit": "F", "param": "flows",
                                "component": "benzene", "value": "55 kmol/h"}]})
    s = series(simulate(model))
    x_B = s["stream.BOTT.x.toluene"]
    assert x_B[-1] == pytest.approx(x_B[-2], abs=1e-5)  # at rest again
    assert x_B[-1] < x_B[0] - 0.05  # after a real change
    D = s["unit.COL.distillate_flow"][-1]
    B = s["unit.COL.bottoms_flow"][-1]
    R = s["unit.COL.reflux_flow"][-1]
    assert (D + B) * 3.6 == pytest.approx(105.0, rel=1e-5)  # nothing vents, nothing piles up
    assert s["unit.COL.vapor_distillate_flow"][-1] == 0.0
    # At rest the hold-up balances are the MESH equations: the dynamic column ends where
    # the steady column with the same reflux, distillate and feed is.
    steady = copy.deepcopy(model)
    steady.pop("dynamics")
    steady["units"][0]["params"]["flows"]["benzene"] = "55 kmol/h"
    steady["units"][1]["params"].update(reflux_ratio=R / D, distillate_rate=f"{D * 3.6} kmol/h")
    r = solve(steady)
    assert r.ok, r.error
    for stream, comp in (("DIST", "benzene"), ("BOTT", "toluene")):
        assert s[f"stream.{stream}.x.{comp}"][-1] == pytest.approx(
            r.streams[stream]["mole_fractions"][comp], abs=2e-5)
    assert s["unit.COL.reboiler_duty"][-1] == pytest.approx(
        r.units["COL"]["reboiler_duty"], rel=1e-3)
    assert s["unit.COL.condenser_duty"][-1] == pytest.approx(
        r.units["COL"]["condenser_duty"], rel=1e-3)


def test_temperature_control_holds_the_control_stage():
    tc = {"id": "TC", "measure": {"unit": "COL", "property": "control_T"},
          "manipulate": {"unit": "COL", "param": "reboiler_duty"}, "action": "reverse",
          "gain": 1, "integral_time": "10 min", "pv_range": ["80 degC", "120 degC"],
          "op_range": ["0 kW", "3000 kW"]}
    step = [{"time": "5 min", "unit": "F", "param": "flows", "component": "benzene",
             "value": "55 kmol/h"},
            {"time": "5 min", "unit": "F", "param": "flows", "component": "toluene",
             "value": "45 kmol/h"}]
    dyn = {"t_end": "3 h", "points": 7, "units": {"COL": {}}, "events": step}
    open_loop = series(simulate(column({**dyn, "controllers": LEVELS})))
    closed = series(simulate(column({**dyn, "controllers": [*LEVELS, tc]})))
    T0 = open_loop["unit.COL.control_T"][0]
    assert abs(open_loop["unit.COL.control_T"][-1] - T0) > 1.0  # the step moves it
    assert closed["unit.COL.control_T"][-1] == pytest.approx(T0, abs=0.05)
    assert closed["controller.TC.op"][-1] != pytest.approx(closed["controller.TC.op"][0])


def test_stage_temperatures_are_series_and_measurable():
    tc = {"id": "T5", "measure": {"unit": "COL", "property": "T_stage_5"},
          "manipulate": {"unit": "COL", "param": "reflux_flow"}, "action": "direct",
          "gain": 0.5, "integral_time": "20 min", "pv_range": ["80 degC", "120 degC"],
          "op_range": ["0 kmol/h", "300 kmol/h"]}
    s = series(simulate(column({"t_end": "10 min", "points": 3, "units": {"COL": {}},
                                "controllers": [*LEVELS, tc]})))
    assert s["controller.T5.pv"][0] == pytest.approx(s["unit.COL.T_stage_5"][0])
    assert len([k for k in s if k.startswith("unit.COL.T_stage_")]) == 15


def test_duty_condenser_vents_what_it_does_not_condense():
    model = column({"t_end": "1 h", "points": 3, "controllers": LEVELS,
                    "units": {"COL": {"condenser_control": "duty"}},
                    "events": [{"time": "5 min", "unit": "COL", "param": "reboiler_duty",
                                "value": "2000 kW"}]})
    s = series(simulate(model))
    assert s["unit.COL.condenser_duty"][-1] == pytest.approx(s["unit.COL.condenser_duty"][0])
    assert s["unit.COL.vapor_distillate_flow"][-1] > 0.0


@pytest.mark.parametrize("extra, message", [
    ({"side_draws": [{"stage": 4, "phase": "liquid", "rate": "5 kmol/h"}]}, "side draws"),
    ({"efficiency": 0.7}, "Murphree efficiency"),
])
def test_unsupported_column_features_are_named(extra, message):
    out = simulate(column({"t_end": "1 min", "units": {"COL": {}}}, extra))
    assert not out["ok"]
    assert message in out["error"]["message"]
    assert out["error"]["where"] == "COL"


def test_partial_condenser_cannot_condense_everything():
    out = simulate(column({"t_end": "1 min", "units": {"COL": {"condenser_control": "total"}}},
                          {"condenser": "partial"}))
    assert not out["ok"]
    assert "cannot condense everything" in out["error"]["message"]


def test_a_stage_temperature_controller_keeps_the_jacobian_sparse():
    tc = {"id": "TC", "measure": {"unit": "COL", "property": "control_T"},
          "manipulate": {"unit": "COL", "param": "reboiler_duty"},
          "gain": 1, "integral_time": "10 min", "op_range": ["0 kW", "3000 kW"]}
    sim, _ = start_rates(column({"t_end": "1 min", "units": {"COL": {}},
                                 "controllers": [*LEVELS, tc]}))
    # The controllers read one stage each, not the whole column: the column's states still
    # share groups three stages apart, and only the controllers' own states need groups of
    # their own (without this every column state would need one: 45 + 9).
    col = sim.dynamic["COL"]
    controller_states = len(sim.x0) - col.N * col.n
    assert len(sim._column_groups()) <= 3 * col.n + controller_states
