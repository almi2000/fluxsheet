"""Pressure-driven dynamics: flows that follow from pressures (dynamics/pressure.py).

Checked against independent results: the valve equations against their US-unit forms, an
adiabatic gas blowdown against a separate integration of the ideal-gas balances (and the
isentrope it must follow), valve and pump flows in a run against the sizing equations and
the pump curve solved by hand, undisturbed runs against the flow-driven ones, and pressure
and level controllers on valve openings against the openings the sizing equations need.
"""

import itertools
import json
import math
from pathlib import Path

import numpy as np
import pytest
from fastapi.testclient import TestClient
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from fluxsheet.dynamics import simulate
from fluxsheet.dynamics.pressure import (
    KV_PER_CV,
    characteristic,
    gas_valve_flow,
    liquid_valve_flow,
    opening_for,
)
from fluxsheet.server import app
from fluxsheet.thermo_package import get_package

EXAMPLES = Path(__file__).parent.parent / "examples"
KMOL_H = 1000.0 / 3600.0
R = 8.314462618
G = 9.80665
GPM = 3.785411784e-3 / 60.0  # m3/s
PSI = 6894.757293168  # Pa
LB = 0.45359237  # kg
FT3 = 0.3048 ** 3  # m3


def series(out):
    assert out["ok"], out.get("error")
    return {s["key"]: s["values"] for s in out["series"]}


# --- valve equations ----------------------------------------------------------------------

def test_liquid_valve_equation_matches_the_us_form():
    # Cv: US gpm of water at 1 psi; Q = Cv sqrt(dP / SG) with SG against water at 60 F.
    Cv, dp, rho = 40.0, 2.5e5, 870.0
    q_us = Cv * math.sqrt((dp / PSI) / (rho / 999.0)) * GPM
    assert liquid_valve_flow(Cv * KV_PER_CV, dp, rho) == pytest.approx(q_us, rel=1e-3)
    # Kv: 1 m3/h of water at 1 bar.
    assert liquid_valve_flow(1.0, 1e5, 1000.0) * 3600.0 == pytest.approx(1.0, rel=1e-12)
    assert liquid_valve_flow(1.0, -1e5, 1000.0) == 0.0  # no reverse flow


def test_gas_valve_equation_matches_the_us_form_and_chokes():
    # ISA/IEC in US units: W [lb/h] = 63.3 Cv Y sqrt(x p1 [psia] rho1 [lb/ft3]).
    Cv, p1, rho1, gamma, xT = 25.0, 8e5, 9.0, 1.3, 0.72
    choke = gamma / 1.4 * xT
    for p2 in (7e5, 5e5, 3e5, 1e5):
        x = min((p1 - p2) / p1, choke)
        Y = 1.0 - x / (3.0 * choke)
        w_us = 63.3 * Cv * Y * math.sqrt(x * p1 / PSI * rho1 / LB * FT3) * LB / 3600.0
        w = gas_valve_flow(Cv * KV_PER_CV, p1, p2, rho1, gamma, xT)
        assert w == pytest.approx(w_us, rel=2e-3), p2
    # Choked: no more flow below p2 = p1 (1 - x_choke), and the flow joins smoothly there.
    p_choke = p1 * (1.0 - choke)
    w_choke = gas_valve_flow(10.0, p1, p_choke, rho1, gamma, xT)
    assert gas_valve_flow(10.0, p1, 0.5 * p_choke, rho1, gamma, xT) == w_choke
    slope = (w_choke - gas_valve_flow(10.0, p1, p_choke + 1.0, rho1, gamma, xT)) / 1.0
    assert abs(slope) < 1e-6 * w_choke


def test_characteristics_and_their_inverse():
    for kind in ("linear", "equal_percentage", "quick_opening"):
        assert characteristic(kind, 0.0, 50.0) == 0.0
        assert characteristic(kind, 1.0, 50.0) == pytest.approx(1.0)
        for x in (0.1, 0.37, 0.8):
            assert opening_for(kind, characteristic(kind, x, 50.0), 50.0) == pytest.approx(
                x, abs=1e-12)
    # Equal percentage: each 10 % of travel multiplies the flow by R^0.1 (above ~20 %).
    ratio = characteristic("equal_percentage", 0.9, 50.0) / characteristic(
        "equal_percentage", 0.8, 50.0)
    assert ratio == pytest.approx(50.0 ** 0.1, rel=0.01)


# --- gas drum ---------------------------------------------------------------------------

def drum(events=(), controllers=(), t_end="10 min", points=11, pressure_driven=True,
         vin=None, vout=None, settings=None):
    """Nitrogen from a 20 bar header through VIN into a 10 bar drum, out through VOUT to 1 bar."""
    return {
        "components": ["nitrogen"], "property_method": "ideal",
        "units": [
            {"id": "HDR", "type": "Feed",
             "params": {"T": 300, "P": "20 bar", "flows": {"nitrogen": "100 kmol/h"}}},
            {"id": "VIN", "type": "Valve", "params": {"P": "10 bar", "Kv": 20, **(vin or {})}},
            {"id": "DRUM", "type": "Flash", "params": {"P": "10 bar", "duty": 0}},
            {"id": "VOUT", "type": "Valve", "params": {"P": "1 bar", "Kv": 30, **(vout or {})}},
            {"id": "FLARE", "type": "Product"}, {"id": "LIQ", "type": "Product"}],
        "streams": [{"id": "S1", "from": "HDR.out", "to": "VIN.in"},
                    {"id": "S2", "from": "VIN.out", "to": "DRUM.in"},
                    {"id": "S3", "from": "DRUM.vapor", "to": "VOUT.in"},
                    {"id": "S4", "from": "VOUT.out", "to": "FLARE.in"},
                    {"id": "S5", "from": "DRUM.liquid", "to": "LIQ.in"}],
        "dynamics": {"t_end": t_end, "points": points, "pressure_driven": pressure_driven,
                     "units": {"DRUM": {"diameter": "1 m", "height": "3 m", **(settings or {})}},
                     "events": list(events), "controllers": list(controllers)},
    }


N2 = get_package(("nitrogen",), "ideal")
M_N2 = N2.mws[0] / 1000.0


def h_n2(T):
    return N2.phase_props("vapor", T, 1e5, [1.0]).H


def cp_n2(T):
    return (h_n2(T + 0.01) - h_n2(T - 0.01)) / 0.02


GAMMA = cp_n2(300.0) / (cp_n2(300.0) - R)  # of the drum's steady-state gas


def iec_gas(kv, p1, p2, T):
    """Molar flow of nitrogen, written out from IEC 60534-2-1 (xT = 0.7)."""
    rho1 = p1 * M_N2 / (R * T)
    choke = GAMMA / 1.4 * 0.7
    x = min((p1 - p2) / p1, choke)
    if x <= 0.0:
        return 0.0
    Y = 1.0 - x / (3.0 * choke)
    return 31.6 * kv * Y * math.sqrt(x * p1 / 1e5 * rho1) / 3600.0 / M_N2  # N6 for Kv


def test_gas_blowdown_matches_an_independent_integration():
    model = drum(events=[{"time": 0, "unit": "VIN", "param": "opening", "value": 0}],
                 t_end="100 s", points=21)
    model["dynamics"]["rtol"] = 1e-8  # so the comparison sees the model, not the tolerance
    out = simulate(model)
    s = series(out)
    opening = s["unit.VOUT.opening"][0]
    assert opening == pytest.approx(s["unit.VOUT.opening"][-1])
    # The opening that passes 100 kmol/h at the steady state, from the sizing equation.
    assert opening * iec_gas(30.0, 10e5, 1e5, 300.0) == pytest.approx(100 * KMOL_H, rel=1e-7)
    assert s["stream.S1.molar_flow"][1:] == [0.0] * 20  # the inlet valve is shut

    # Ideal gas in a rigid, adiabatic vessel: dN/dt = -F, N cv dT/dt = -F R T, P = N R T / V.
    V = math.pi / 4.0 * 3.0
    N0 = 10e5 * V / (R * 300.0)

    def rhs(t, y):
        N, T = y
        F = iec_gas(30.0 * opening, N * R * T / V, 1e5, T)
        return [-F, -F * R * T / (N * (cp_n2(T) - R))]

    ref = solve_ivp(rhs, (0.0, 100.0), [N0, 300.0], t_eval=out["time"], method="LSODA",
                    rtol=1e-10, atol=[1e-8, 1e-9])
    checked = 0
    for i, t in enumerate(out["time"]):
        N, T = ref.y[:, i]
        P = N * R * T / V
        if P < 1.2e5:  # close to the back pressure the square root is rounded off
            break
        assert s["unit.DRUM.P"][i] == pytest.approx(P, rel=1e-6), t
        assert s["unit.DRUM.T"][i] == pytest.approx(T, rel=1e-6), t
        assert s["stream.S3.molar_flow"][i] == pytest.approx(
            iec_gas(30.0 * opening, P, 1e5, T), rel=1e-5), t
        checked += 1
    assert checked >= 5
    # The gas left in the drum expands along an isentrope: int cp/T dT = R ln(P/P0).
    for P, T in zip(s["unit.DRUM.P"], s["unit.DRUM.T"]):
        Ts = np.linspace(300.0, T, 201)
        ds = np.trapezoid([cp_n2(x) / x for x in Ts], Ts)
        assert ds == pytest.approx(R * math.log(P / 10e5), rel=1e-5, abs=1e-6)
    assert s["unit.DRUM.P"][-1] == pytest.approx(1e5, rel=1e-4)  # it ends at the flare


def test_undisturbed_runs_equal_the_flow_driven_ones():
    # (A gas-only drum in a flow-driven run needs initial_level 0; a pressure-driven one
    # starts without liquid unless asked to.)
    for model in (drum(settings={"initial_level": 0}), tank()):
        flow = series(simulate({**model, "dynamics": {**model["dynamics"],
                                                      "pressure_driven": False}}))
        press = series(simulate(model))
        for key, values in flow.items():
            if key.startswith(("stream.", "unit.")):
                for a, b in zip(values, press[key]):
                    assert b == pytest.approx(a, rel=1e-9, abs=1e-12), key
        for key, values in press.items():  # and it stays at rest
            assert max(values) - min(values) <= 1e-9 * max(abs(values[0]), 1e-9), key


PC = {"id": "PC", "measure": {"unit": "DRUM", "property": "P"},
      "manipulate": {"unit": "VOUT", "param": "opening"}, "action": "direct",
      "gain": 1.0, "integral_time": "30 s", "op_range": [0, 1],
      "pv_range": ["5 bar", "15 bar"]}


def test_pressure_controller_on_the_outlet_valve_leaves_no_offset():
    out = simulate(drum(controllers=[PC], t_end="20 min", points=21, events=[
        {"time": "1 min", "unit": "VIN", "param": "opening", "value": 0.7}]))
    s = series(out)
    assert max(s["unit.DRUM.P"]) > 10.2e5  # the extra inflow raised the pressure
    assert s["unit.DRUM.P"][-1] == pytest.approx(10e5, rel=1e-4)
    # The inflow follows the inlet valve's sizing equation from the header to the drum.
    for P, T_in, F in zip(s["unit.DRUM.P"][2:], s["stream.S1.T"][2:],
                          s["stream.S1.molar_flow"][2:]):
        assert F == pytest.approx(iec_gas(20.0 * 0.7, 20e5, P, T_in), rel=1e-6)
    # At the end out = in, through the opening the sizing equation needs at 10 bar.
    F_in = s["stream.S1.molar_flow"][-1]
    assert s["stream.S3.molar_flow"][-1] == pytest.approx(F_in, rel=1e-4)
    T = s["unit.DRUM.T"][-1]
    assert s["unit.VOUT.opening"][-1] == pytest.approx(
        F_in / iec_gas(30.0, 10e5, 1e5, T), rel=2e-4)
    assert s["controller.PC.op"][-1] == pytest.approx(s["unit.VOUT.opening"][-1])


def test_closing_the_outlet_builds_up_pressure_towards_the_header():
    out = simulate(drum(t_end="30 min", points=7, events=[
        {"time": 0, "unit": "VOUT", "param": "opening", "value": 0}]))
    s = series(out)
    P = s["unit.DRUM.P"]
    assert all(b > a for a, b in itertools.pairwise(P) if a < 19.9e5)
    assert P[-1] == pytest.approx(20e5, rel=1e-3)  # until the inflow stops
    assert s["stream.S1.molar_flow"][-1] < 1e-3 * s["stream.S1.molar_flow"][0]
    # Moles in = increase of the hold-up: N = P V / (R T) in the drum.
    V = math.pi / 4.0 * 3.0
    n = [p * V / (R * T) for p, T in zip(P, s["unit.DRUM.T"])]
    flows = s["stream.S1.molar_flow"]
    fine = simulate(drum(t_end="30 min", points=361, events=[
        {"time": 0, "unit": "VOUT", "param": "opening", "value": 0}]))
    f = series(fine)
    added = np.trapezoid(f["stream.S1.molar_flow"], fine["time"])
    assert n[-1] - n[0] == pytest.approx(added, rel=2e-3)
    assert flows[0] == pytest.approx(100 * KMOL_H)


# --- liquid tank and pump -----------------------------------------------------------------

def tank(events=(), controllers=(), t_end="30 min", points=7, pump=False):
    """Water from a 3 bar main through VIN into an open tank at 1.2 bar; out through VOUT,
    driven by the level, to a product at 1.2 bar (or through a pump with a curve and VOUT to
    4 bar)."""
    units = [
        {"id": "MAIN", "type": "Feed",
         "params": {"T": "25 degC", "P": "3 bar", "flows": {"water": "50 kmol/h"}}},
        {"id": "VIN", "type": "Valve", "params": {"P": "1.2 bar", "Kv": 3}},
        {"id": "TANK", "type": "MixingTank", "params": {"P": "1.2 bar"}},
        {"id": "VOUT", "type": "Valve",
         "params": {"P": "4 bar" if pump else "1.2 bar", "Kv": 4 if pump else 10,
                    "characteristic": "equal_percentage"}},
        {"id": "PROD", "type": "Product"},
    ]
    streams = [{"id": "S1", "from": "MAIN.out", "to": "VIN.in"},
               {"id": "S2", "from": "VIN.out", "to": "TANK.in"},
               {"id": "S4", "from": "VOUT.out", "to": "PROD.in"}]
    if pump:
        units.append({"id": "PUMP", "type": "Pump", "params": {"curve": PUMP_CURVE}})
        streams += [{"id": "S3", "from": "TANK.liquid", "to": "PUMP.in"},
                    {"id": "S3b", "from": "PUMP.out", "to": "VOUT.in"}]
    else:
        streams.append({"id": "S3", "from": "TANK.liquid", "to": "VOUT.in"})
    return {
        "components": ["water"], "property_method": "ideal",
        "units": units, "streams": streams,
        "dynamics": {"t_end": t_end, "points": points, "pressure_driven": True,
                     "units": {"TANK": {"diameter": "1 m", "height": "2 m"}},
                     "events": list(events), "controllers": list(controllers)},
    }


PUMP_CURVE = {"flow_unit": "m3/h", "head_unit": "m",
              "points": [[0, 60], [0.5, 58], [1.0, 52], [1.5, 42], [2.0, 28]]}
WATER = get_package(("water",), "ideal")
M_W = WATER.mws[0] / 1000.0


def rho_w(T, P):
    return M_W / WATER.phase_props("liquid", T, P, [1.0]).V


def test_valve_flows_in_a_run_follow_the_sizing_equation():
    out = simulate(tank(t_end="8 h", points=17, events=[
        {"time": "2 min", "unit": "VIN", "param": "opening", "value": 0.25}]))
    s = series(out)
    rho_in = rho_w(298.15, 3e5)
    f_out = characteristic("equal_percentage", s["unit.VOUT.opening"][0], 50.0)
    for i, t in enumerate(out["time"]):
        x_in = s["unit.VIN.opening"][i]
        # Into the tank: 3 bar main to the tank's 1.2 bar (the inlet is above the liquid).
        q_in = 3.0 * x_in * math.sqrt((3e5 - 1.2e5) / 1e5 / (rho_in / 1000.0)) / 3600.0
        assert s["stream.S1.molar_flow"][i] == pytest.approx(q_in * rho_in / M_W, rel=1e-6)
        # Out: the static head of the level, to the product at the tank's pressure.
        T = s["unit.TANK.T"][i]
        rho = rho_w(T, 1.2e5)
        dp = rho * G * s["unit.TANK.level"][i]
        q_out = 10.0 * f_out * math.sqrt(dp / 1e5 / (rho / 1000.0)) / 3600.0
        assert s["stream.S3.molar_flow"][i] == pytest.approx(q_out * rho / M_W, rel=1e-6), t
    # Self-regulating: the level rises until the outflow matches the larger inflow, at
    # L = L0 (F / F0)^2 since the outflow goes with the square root of the level.
    F0, F = s["stream.S1.molar_flow"][0], s["stream.S1.molar_flow"][-1]
    assert s["stream.S3.molar_flow"][-1] == pytest.approx(F, rel=2e-3)
    assert s["unit.TANK.level"][-1] == pytest.approx(1.0 * (F / F0) ** 2, rel=5e-3)


def test_pump_curve_sets_the_flow_against_the_valve():
    out = simulate(tank(pump=True, t_end="20 min", points=11, events=[
        {"time": "1 min", "unit": "VOUT", "param": "opening", "value": 0.9}]))
    s = series(out)
    f0 = characteristic("equal_percentage", s["unit.VOUT.opening"][0], 50.0)
    flows = [q for q, _h in PUMP_CURVE["points"]]
    heads = [h for _q, h in PUMP_CURVE["points"]]
    head = PchipInterpolator(np.array(flows) / 3600.0, np.array(heads) * G)
    for i, t in enumerate(out["time"]):
        f = f0 if t < 60.0 else characteristic("equal_percentage", 0.9, 50.0)
        T = s["unit.TANK.T"][i]
        rho = rho_w(T, 1.2e5)
        p_suction = 1.2e5 + rho * G * s["unit.TANK.level"][i]

        def excess(Q, f=f, rho=rho, p_suction=p_suction):
            dp = p_suction + rho * float(head(Q)) - 4e5  # the 4 bar product
            return 4.0 * f * math.sqrt(max(dp, 0.0) / 1e5 / (rho / 1000.0)) / 3600.0 - Q

        Q = brentq(excess, 0.0, 2.0 / 3600.0, xtol=1e-15)
        assert s["stream.S3.molar_flow"][i] == pytest.approx(Q * rho / M_W, rel=1e-5), t
    assert s["stream.S3.molar_flow"][-1] > 1.2 * s["stream.S3.molar_flow"][0]


LC = {"id": "LC", "measure": {"unit": "TANK", "property": "level"},
      "manipulate": {"unit": "VOUT", "param": "opening"}, "action": "direct",
      "gain": 3, "integral_time": "3 min", "pv_range": ["0 m", "2 m"], "op_range": [0, 1]}


def test_level_controller_on_the_outlet_valve():
    out = simulate(tank(controllers=[LC], t_end="2 h", points=13, events=[
        {"time": "1 min", "unit": "VIN", "param": "opening", "value": 0.25}]))
    s = series(out)
    assert max(s["unit.TANK.level"]) > 1.005
    assert s["unit.TANK.level"][-1] == pytest.approx(1.0, abs=2e-3)
    F = s["stream.S1.molar_flow"][-1]
    assert s["stream.S3.molar_flow"][-1] == pytest.approx(F, rel=1e-3)
    # The opening the equal-percentage valve needs for that flow at a 1 m level.
    T = s["unit.TANK.T"][-1]
    rho = rho_w(T, 1.2e5)
    q_full = 10.0 * math.sqrt(rho * G * 1.0 / 1e5 / (rho / 1000.0)) / 3600.0
    need = F * M_W / rho / q_full
    assert s["unit.VOUT.opening"][-1] == pytest.approx(
        opening_for("equal_percentage", need, 50.0), abs=2e-3)


# --- two-phase drum -----------------------------------------------------------------------

def flash_drum():
    return json.loads((EXAMPLES / "flash_drum_pressure_control.json").read_text("utf-8"))


def test_two_phase_drum_with_pressure_and_level_control():
    model = flash_drum()
    model["dynamics"]["t_end"] = "40 min"
    model["dynamics"]["points"] = 9
    out = simulate(model)
    s = series(out)
    assert s["unit.DRUM.P"][0] == pytest.approx(8e5)
    assert max(s["unit.DRUM.P"]) > 8.1e5
    assert s["unit.DRUM.P"][-1] == pytest.approx(8e5, rel=2e-3)
    assert s["unit.DRUM.level"][-1] == pytest.approx(1.5, abs=0.02)
    # The feed valve passes what its sizing equation gives from the 15 bar header.
    pkg = get_package(("propane", "n-butane"), "peng-robinson")
    feed = pkg.flash([0.5, 0.5], T=330.0, P=15e5)
    M = sum(0.5 * mw for mw in pkg.mws) / 1000.0
    rho = M / feed.V
    for i, t in enumerate(out["time"]):
        x = s["unit.VIN.opening"][i]
        q = 10.0 * x * math.sqrt((15e5 - s["unit.DRUM.P"][i]) / 1e5 / (rho / 1000.0)) / 3600.0
        assert s["stream.S1.molar_flow"][i] == pytest.approx(q * rho / M, rel=1e-6), t
    # Close to the new steady state, what comes in leaves as vapor and liquid.
    out_flow = s["stream.S3.molar_flow"][-1] + s["stream.S5.molar_flow"][-1]
    assert out_flow == pytest.approx(s["stream.S1.molar_flow"][-1], rel=2e-2)


def test_example_through_the_api():
    model = flash_drum()
    model["dynamics"]["t_end"] = "3 min"
    model["dynamics"]["points"] = 4
    body = TestClient(app).post("/api/dynamic", json=model).json()
    assert body["ok"], body.get("error")
    keys = {s["key"] for s in body["series"]}
    assert {"unit.DRUM.P", "unit.VV.opening", "unit.VL.dP"} <= keys


# --- what cannot be pressure-driven -------------------------------------------------------

def test_a_valve_too_small_is_explained():
    out = simulate(drum(vout={"Kv": 1}))
    assert not out["ok"]
    assert "too small" in out["error"]["message"]
    assert out["error"]["where"] == "VOUT"


def test_vapor_line_of_a_vessel_without_gas_is_explained():
    model = tank()
    model["units"].append({"id": "VENT", "type": "Valve", "params": {"P": "1 bar", "Kv": 1}})
    model["units"].append({"id": "AIR", "type": "Product"})
    model["streams"] += [{"id": "S5", "from": "TANK.vapor", "to": "VENT.in"},
                         {"id": "S6", "from": "VENT.out", "to": "AIR.in"}]
    out = simulate(model)
    assert not out["ok"]
    assert "holds no gas" in out["error"]["message"]


def test_controller_on_an_outflow_a_valve_sets_is_explained():
    lc = {**LC, "manipulate": {"unit": "TANK", "param": "liquid_flow"},
          "op_range": ["0 kmol/h", "100 kmol/h"]}
    out = simulate(tank(controllers=[lc]))
    assert not out["ok"]
    assert "follows from the pressures" in out["error"]["message"]
    assert "opening of VOUT" in out["error"]["hint"]


def test_valve_off_a_line_keeps_its_flow_with_a_warning():
    # A valve between a mixer and the tank: the mixer is no pressure node.
    model = tank()
    model["units"][1]["params"] = {"P": "1.2 bar"}  # VIN: no flow coefficient
    model["units"] += [{"id": "MIX", "type": "Mixer"},
                       {"id": "V2", "type": "Valve", "params": {"P": "1.2 bar", "Kv": 3}}]
    model["streams"] = [s for s in model["streams"] if s["id"] != "S2"] + [
        {"id": "S2a", "from": "VIN.out", "to": "MIX.in"},
        {"id": "S2b", "from": "MIX.out", "to": "V2.in"},
        {"id": "S2", "from": "V2.out", "to": "TANK.in"}]
    out = simulate(model)
    s = series(out)
    assert any(w.startswith("V2:") for w in out["warnings"])
    assert s["stream.S1.molar_flow"][-1] == pytest.approx(50 * KMOL_H)
