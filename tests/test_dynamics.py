"""Dynamic simulation, calculator blocks and Broyden recycle convergence.

The dynamic tests check physics with known answers: a hold-up at steady state stays put, a
well-mixed tank answers a feed step with exp(-t/τ), an uncontrolled level integrates the
flow imbalance, a PI controller leaves no offset, a reactor settles where the steady-state
solver puts it, and a first-order lag reaches 63 % after one time constant.
"""

import json
import math
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from fluxsheet import solve
from fluxsheet.dynamics import simulate
from fluxsheet.dynamics.simulator import to_csv
from fluxsheet.server import app
from fluxsheet.thermo_package import get_package

EXAMPLES = Path(__file__).parent.parent / "examples"
KMOL_H = 1000.0 / 3600.0


def tank(events=(), controllers=(), t_end="30 min", points=31, settings=None):
    return {
        "components": ["water", "ethanol"], "property_method": "nrtl",
        "units": [
            {"id": "FEED", "type": "Feed",
             "params": {"T": "25 degC", "P": "1.5 bar",
                        "flows": {"water": "50 kmol/h", "ethanol": "5 kmol/h"}}},
            {"id": "TANK", "type": "MixingTank", "params": {"P": "1.2 bar"}},
            {"id": "PROD", "type": "Product"},
        ],
        "streams": [{"id": "S1", "from": "FEED.out", "to": "TANK.in"},
                    {"id": "S2", "from": "TANK.liquid", "to": "PROD.in"}],
        "dynamics": {"t_end": t_end, "points": points,
                     "units": {"TANK": settings or {"diameter": "1 m", "height": "2 m"}},
                     "controllers": list(controllers), "events": list(events)},
    }


def series(out):
    assert out["ok"], out.get("error")
    return {s["key"]: s["values"] for s in out["series"]}


LEVEL_CONTROL = {
    "id": "LC", "measure": {"unit": "TANK", "property": "level"},
    "manipulate": {"unit": "TANK", "param": "liquid_flow"}, "action": "direct",
    "gain": 2, "integral_time": "5 min", "pv_range": ["0 m", "2 m"],
    "op_range": ["0 kmol/h", "120 kmol/h"],
}


def test_steady_state_stays_put():
    s = series(simulate(tank(controllers=[LEVEL_CONTROL])))
    for key in ("unit.TANK.level", "unit.TANK.T", "stream.S2.molar_flow", "controller.LC.op"):
        values = s[key]
        assert max(values) - min(values) <= 1e-6 * max(abs(values[0]), 1e-9), key


def test_mixed_tank_answers_a_composition_step_with_its_residence_time():
    # Same total flow, more ethanol: the hold-up (in mol) stays constant, so the outlet
    # composition follows x_new + (x0 - x_new) exp(-t F / N) exactly.
    model = tank(events=[
        {"time": 0, "unit": "FEED", "param": "flows", "component": "water", "value": "45 kmol/h"},
        {"time": 0, "unit": "FEED", "param": "flows", "component": "ethanol",
         "value": "10 kmol/h"},
    ], t_end="40 min", points=9)
    out = simulate(model)
    s = series(out)
    N = s["unit.TANK.holdup"][0]
    F = 55 * KMOL_H
    assert s["unit.TANK.holdup"][-1] == pytest.approx(N, rel=1e-5)
    x0, x_new = 5 / 55, 10 / 55
    for t, x in zip(out["time"], s["stream.S2.x.ethanol"]):
        assert x == pytest.approx(x_new + (x0 - x_new) * math.exp(-t * F / N), abs=2e-5)


def test_uncontrolled_level_integrates_the_imbalance():
    model = tank(events=[{"time": 0, "unit": "FEED", "param": "flows", "component": "water",
                          "value": "60 kmol/h"}])
    out = simulate(model)
    s = series(out)
    # The moles add up exactly: 10 kmol/h more in than out, linearly in time.
    for t, n in zip(out["time"], s["unit.TANK.holdup"]):
        assert n - s["unit.TANK.holdup"][0] == pytest.approx(10 * KMOL_H * t, rel=1e-4, abs=1e-6)
    # The level is the liquid volume over the cross section and rises steadily.
    area = math.pi / 4
    for level, volume in zip(s["unit.TANK.level"], s["unit.TANK.liquid_volume"]):
        assert level == pytest.approx(volume / area, rel=1e-12)
    # It rises a little faster as the dilution's volume contraction fades out.
    rise = [b - a for a, b in zip(s["unit.TANK.level"], s["unit.TANK.level"][1:])]
    assert min(rise) > 0 and max(rise) / min(rise) < 1.15


def test_pi_level_control_leaves_no_offset():
    model = tank(controllers=[LEVEL_CONTROL], t_end="3 h", points=13, events=[
        {"time": "1 min", "unit": "FEED", "param": "flows", "component": "water",
         "value": "60 kmol/h"}])
    s = series(simulate(model))
    assert s["unit.TANK.level"][-1] == pytest.approx(1.0, abs=2e-3)
    assert s["stream.S2.molar_flow"][-1] == pytest.approx(65 * KMOL_H, rel=2e-3)
    assert max(s["unit.TANK.level"]) > 1.01  # it did move before the controller caught it


def test_setpoint_ramp_is_followed():
    model = tank(controllers=[LEVEL_CONTROL], t_end="3 h", points=13, events=[
        {"time": "10 min", "controller": "LC", "value": "1.4 m", "ramp": "10 min"}])
    out = simulate(model)
    s = series(out)
    sp = dict(zip(out["time"], s["controller.LC.sp"]))
    assert sp[0.0] == pytest.approx(1.0)
    assert sp[900.0] == pytest.approx(1.2)  # half way along the ramp
    assert s["unit.TANK.level"][-1] == pytest.approx(1.4, abs=5e-3)


# --- reactor and lag --------------------------------------------------------------------

def isomerization(F=1.0, energy=None):
    T, P, tau = 400.0, 1e5, 20.0
    V = tau * 1.0 * 8.314462618 * T / P
    return {
        "components": ["n-butane", "isobutane"], "property_method": "ideal",
        "reactions": [{"id": "ISO", "stoichiometry": {"n-butane": -1, "isobutane": 1},
                       "kinetics": {"k0": 0.05, "Ea": "0 kJ/mol"}}],
        "units": [
            {"id": "F", "type": "Feed", "params": {"T": T, "P": P, "flows": {"n-butane": F}}},
            {"id": "R", "type": "CSTR", "params": {"reactions": ["ISO"], "V": V,
                                                    "phase": "gas", **(energy or {"T": T})}},
            {"id": "PR", "type": "Product"}],
        "streams": [{"id": "IN", "from": "F.out", "to": "R.in"},
                    {"id": "OUT", "from": "R.out", "to": "PR.in"}],
    }


def test_cstr_settles_at_the_new_steady_state():
    start = solve(isomerization())
    duty = start.units["R"]["duty"]
    model = isomerization()
    model["dynamics"] = {"t_end": 300, "points": 7, "units": {"R": {}}, "events": [
        {"time": 0, "unit": "F", "param": "flows", "component": "n-butane", "value": 2.0}]}
    s = series(simulate(model))
    # The duty stays at its first value, so compare with a steady solve at that duty.
    end = solve(isomerization(F=2.0, energy={"duty": duty}))
    assert end.ok, end.error
    assert s["stream.OUT.x.n-butane"][-1] == pytest.approx(
        end.streams["OUT"]["mole_fractions"]["n-butane"], abs=1e-5)
    assert s["unit.R.T"][-1] == pytest.approx(end.streams["OUT"]["T"], abs=0.02)
    assert s["stream.OUT.molar_flow"][-1] == pytest.approx(2.0, rel=1e-5)


def test_heater_lag_reaches_63_percent_after_one_time_constant():
    model = tank()
    model["units"][1] = {"id": "H", "type": "Heater", "params": {"duty": "50 kW"}}
    model["streams"] = [{"id": "S1", "from": "FEED.out", "to": "H.in"},
                        {"id": "S2", "from": "H.out", "to": "PROD.in"}]
    model["dynamics"] = {"t_end": "2 min", "points": 3, "units": {"H": {"time_constant": "1 min"}},
                         "events": [{"time": 0, "unit": "H", "param": "duty", "value": "80 kW"}]}
    s = series(simulate(model))
    # The outlet enthalpy is first order: h(t) = h_new + (h0 - h_new) exp(-t/τ).
    pkg = get_package(("water", "ethanol"), "nrtl")
    z, P = [50 / 55, 5 / 55], 1.5e5

    def h(T):
        return pkg.flash(z, T=T, P=P).H

    h0, h1, h2 = (h(T) for T in s["stream.S2.T"])
    feed = solve({**model, "dynamics": None})
    dh = 30e3 / (55 * KMOL_H)  # 30 kW more on 55 kmol/h
    assert h0 == pytest.approx(feed.streams["S2"]["enthalpy"], abs=1e-3)
    assert (h1 - h0) / dh == pytest.approx(1 - math.exp(-1), abs=5e-4)
    assert (h2 - h0) / dh == pytest.approx(1 - math.exp(-2), abs=5e-4)


# --- errors and interfaces --------------------------------------------------------------

def test_unit_without_dynamic_model_is_explained():
    model = tank()
    model["dynamics"]["units"] = {"FEED": {}}
    out = simulate(model)
    assert not out["ok"]
    assert "has no dynamic model" in out["error"]["message"]
    assert "MixingTank" in out["error"]["hint"]


def test_unknown_dynamic_setting_is_explained():
    out = simulate(tank(settings={"diameter": "1 m", "colour": "red"}))
    assert not out["ok"]
    assert "unknown dynamic setting 'colour'" in out["error"]["message"]
    assert "diameter" in out["error"]["hint"]


def test_example_and_api():
    model = json.loads((EXAMPLES / "esterification_control.json").read_text(encoding="utf-8"))
    model["dynamics"]["t_end"] = "15 min"
    model["dynamics"]["points"] = 4
    client = TestClient(app)
    body = client.post("/api/dynamic", json=model).json()
    assert body["ok"], body.get("error")
    assert body["time"][-1] == pytest.approx(900.0)
    csv_text = to_csv(body)
    assert csv_text.splitlines()[0].startswith("time [s],")
    types = {t["type"]: t for t in client.get("/api/unit-types").json()}
    assert "level" in [f["name"] for f in types["MixingTank"]["dynamic"]["result_fields"]]
    assert "dynamic" not in types["Mixer"]


# --- calculator blocks ------------------------------------------------------------------

def gas_recycle():
    return json.loads((EXAMPLES / "gas_recycle.json").read_text(encoding="utf-8"))


def test_calculator_sets_parameters_from_inputs_and_results():
    model = gas_recycle()
    model["calculators"] = [
        {"id": "C1", "variables": {"tf": {"unit": "FEED", "param": "T"}},
         "assign": [{"unit": "CHILLER", "param": "T", "expression": "tf - 70"}]},
        {"id": "C2", "variables": {"ngl": {"stream": "NGL-OUT", "property": "molar_flow"}},
         "assign": [{"unit": "SPLIT", "param": "split",
                     "expression": "min(0.9, 0.3 + ngl * 3.6 / 100)"}]},
    ]
    r = solve(model)
    assert r.ok, r.error
    assert r.streams["S3"]["T"] == pytest.approx(303.15 - 70)
    (assigned,) = r.calculators[1]["assigned"]
    ngl = r.streams["NGL-OUT"]["molar_flow"]
    assert assigned["value"] == pytest.approx(0.3 + ngl * 3.6 / 100, rel=1e-6)


@pytest.mark.parametrize("expression, message", [
    ("__import__('os')", "only these functions can be called"),
    ("tf.real", "Attribute is not allowed"),
    ("tf + nope", "unknown name 'nope'"),
    ("'text'", "is not a number"),
])
def test_calculator_rejects_anything_but_arithmetic(expression, message):
    model = gas_recycle()
    model["calculators"] = [{"id": "C1", "variables": {"tf": {"unit": "FEED", "param": "T"}},
                             "assign": [{"unit": "CHILLER", "param": "T",
                                         "expression": expression}]}]
    r = solve(model)
    assert not r.ok
    assert message in r.error["message"]
    assert r.error["where"] == "C1"


# --- Broyden ----------------------------------------------------------------------------

def test_broyden_matches_wegstein():
    results = {}
    for acc in ("wegstein", "broyden"):
        model = gas_recycle()
        model["solver"] = {"acceleration": acc, "tolerance": 1e-8}
        r = solve(model)
        assert r.ok and r.converged, r.error
        results[acc] = r
    for sid, s in results["wegstein"].streams.items():
        b = results["broyden"].streams[sid]
        assert b["molar_flow"] == pytest.approx(s["molar_flow"], rel=1e-6)
        assert b["T"] == pytest.approx(s["T"], abs=1e-4)
