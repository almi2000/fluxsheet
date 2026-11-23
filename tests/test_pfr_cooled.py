"""Plug flow reactor cooled or heated along its length, co- and countercurrent.

Checked against the closed-form profiles of a double-pipe exchanger, the effectiveness-NTU
relation, the condensing-steam limit, and an independent solution of the same
balances (an enthalpy formulation integrated by scipy, the countercurrent case by
collocation with ``solve_bvp`` instead of shooting).
"""

import functools
import math

import numpy as np
import pytest
from scipy.integrate import solve_bvp, solve_ivp
from scipy.interpolate import CubicSpline

from fluxsheet import solve
from fluxsheet.flowsheet import UtilityModel
from fluxsheet.reactions import R
from fluxsheet.thermo_package import get_package
from fluxsheet.units.reactors import PFR
from fluxsheet.utilities import parse_utilities

# Argon is a monatomic ideal gas, its heat capacity 5/2 R at every temperature. A PFR that
# carries only argon is a double-pipe exchanger with constant heat capacities on both
# sides. The reactor needs a reaction; no butane is fed, so its rate is zero.
INERT = ["argon", "n-butane", "isobutane"]
NO_RATE = {"id": "ISO", "stoichiometry": {"n-butane": -1, "isobutane": 1},
           "kinetics": {"k0": 1.0}}
CP_WATER = 4000.0  # J/(kg K) of the fluid utility
CW = {"id": "CW", "type": "fluid", "T_in": "300 K", "T_out": "320 K", "cp": "4000 J/kgK"}
LPS = {"id": "LPS", "type": "steam", "P": "3 bar"}
L, D, U = 10.0, 0.05, 100.0
UA = U * math.pi * D * L
F_AR, T_FEED = 10.0, 600.0  # mol/s, K


def cp_argon() -> float:
    pkg = get_package(tuple(INERT), "ideal")
    h = [pkg.phase_props("vapor", T, 2e5, [1.0, 0.0, 0.0]).H for T in (300.0, 800.0)]
    return (h[1] - h[0]) / 500.0


def run(pfr, *, components=INERT, feed=None, T=T_FEED, P="2 bar", reactions=(NO_RATE,),
        utilities=(CW, LPS), medium=None):
    """Feed -> PFR -> product, with an optional second feed on the medium side."""
    units = [{"id": "F", "type": "Feed",
              "params": {"T": T, "P": P, "flows": feed or {"argon": F_AR}}},
             {"id": "R", "type": "PFR",
              "params": {"length": L, "diameter": D, "phase": "gas",
                         "reactions": [r["id"] for r in reactions], **pfr}}]
    streams = [{"id": "IN", "from": "F.out", "to": "R.in"}, {"id": "OUT", "from": "R.out"}]
    if medium is not None:
        units.append({"id": "M", "type": "Feed", "params": medium})
        streams += [{"id": "MIN", "from": "M.out", "to": "R.coolant_in"},
                    {"id": "MOUT", "from": "R.coolant_out"}]
    r = solve({"components": list(components), "property_method": "ideal",
               "reactions": list(reactions), "utilities": list(utilities),
               "units": units, "streams": streams})
    assert r.ok, r.error
    return r


def series(unit, label):
    conditions = unit["profiles"]["conditions"]
    values = next(s["values"] for s in conditions["series"] if s["label"] == label)
    return np.array(conditions["x"]) / L, np.array(values)


def double_pipe(x, C_p, C_m, T_p0, T_m_in, counter):
    """Closed-form temperatures of both sides at the fractions x of the length.

    The difference theta = T - T_m decays as exp(-k x) with k = UA (1/C_p + 1/C_m)
    cocurrent and UA (1/C_p - 1/C_m) countercurrent; the energy balance
    C_p (T - T_p0) = -/+ C_m (T_m - T_m0) gives each side from it. Countercurrent, the
    medium temperature at the reactor inlet, T_m0, follows from T_m(1) = T_m_in.
    """
    if not counter:
        k = UA * (1.0 / C_p + 1.0 / C_m)
        theta = (T_p0 - T_m_in) * np.exp(-k * x)
        T = (C_p * T_p0 + C_m * (theta + T_m_in)) / (C_p + C_m)
        return T, T - theta
    k = UA * (1.0 / C_p - 1.0 / C_m)
    e = math.exp(-k)
    T_m0 = ((C_p - C_m) * T_m_in - C_p * T_p0 * (1.0 - e)) / (C_p * e - C_m)
    theta = (T_p0 - T_m0) * np.exp(-k * x)
    T_m = (C_p * T_p0 - C_m * T_m0 - C_p * theta) / (C_p - C_m)
    return T_m + theta, T_m


def effectiveness(NTU: float, Cr: float, counter: bool) -> float:
    """Effectiveness of a double-pipe exchanger (e.g. Incropera, Fundamentals of Heat and
    Mass Transfer, Table 11.3)."""
    if not counter:
        return (1.0 - math.exp(-NTU * (1.0 + Cr))) / (1.0 + Cr)
    return (1.0 - math.exp(-NTU * (1.0 - Cr))) / (1.0 - Cr * math.exp(-NTU * (1.0 - Cr)))


# ---- double-pipe exchanger: no reaction, constant heat capacities ------------------------

@pytest.mark.parametrize("arrangement", ["cocurrent", "countercurrent"])
@pytest.mark.parametrize("flow", [0.03, 0.1])  # medium heat capacity flow below / above
def test_double_pipe_profiles_and_duty(arrangement, flow):
    counter = arrangement == "countercurrent"
    r = run({"utility": "CW", "U": U, "utility_flow": flow, "arrangement": arrangement})
    u = r.units["R"]
    C_p, C_m = F_AR * cp_argon(), flow * CP_WATER
    x, T = series(u, "Temperature")
    _, T_m = series(u, "Medium temperature")
    T_ref, T_m_ref = double_pipe(x, C_p, C_m, T_FEED, 300.0, counter)
    assert np.allclose(T, T_ref, rtol=0, atol=1e-5)
    assert np.allclose(T_m, T_m_ref, rtol=0, atol=1e-5)
    # The duty from the effectiveness-NTU relation, and the outlets from it.
    C_min, C_max = min(C_p, C_m), max(C_p, C_m)
    Q = -effectiveness(UA / C_min, C_min / C_max, counter) * C_min * (T_FEED - 300.0)
    assert u["wall_duty"] == pytest.approx(Q, rel=1e-7)
    assert u["duty"] == pytest.approx(Q, rel=1e-7)
    assert u["T_out"] == pytest.approx(T_FEED + Q / C_p, abs=1e-5)
    assert u["utility_T_out"] == pytest.approx(300.0 - Q / C_m, abs=1e-5)
    assert u["UA"] == pytest.approx(UA) and u["heat_transfer_area"] == pytest.approx(UA / U)
    assert u["LMTD"] == pytest.approx(-Q / UA)
    if counter:  # the boundary condition the shooting has to meet
        assert T_m[-1] == pytest.approx(300.0, abs=1e-6)
        assert T_m[0] == pytest.approx(u["utility_T_out"], abs=1e-6)
    else:
        assert T_m[0] == 300.0
    # Cooling only: the hottest point is the inlet.
    assert u["T_hot_spot"] == pytest.approx(T_FEED) and u["hot_spot_position"] == 0.0


@pytest.mark.parametrize("arrangement", ["cocurrent", "countercurrent"])
def test_utility_flow_follows_from_its_outlet_temperature(arrangement):
    counter = arrangement == "countercurrent"
    u = run({"utility": "CW", "U": U, "arrangement": arrangement}).units["R"]
    assert u["utility_T_out"] == pytest.approx(320.0, abs=1e-6)
    C_m = u["utility_flow"] * CP_WATER
    assert u["wall_duty"] == pytest.approx(-C_m * 20.0, rel=1e-9)
    # With that flow the closed form brings the medium to 320 K, and the argon to T_out.
    x, T = series(u, "Temperature")
    T_ref, T_m_ref = double_pipe(x, F_AR * cp_argon(), C_m, T_FEED, 300.0, counter)
    assert np.allclose(T, T_ref, rtol=0, atol=1e-5)
    assert (T_m_ref[0] if counter else T_m_ref[-1]) == pytest.approx(320.0, abs=1e-5)


@pytest.mark.parametrize("arrangement", ["cocurrent", "countercurrent"])
def test_steam_heats_at_its_saturation_temperature(arrangement):
    ut = parse_utilities([UtilityModel(**LPS)])["LPS"]
    Ts = ut.T_in
    u = run({"utility": "LPS", "U": U, "arrangement": arrangement}, T=300.0).units["R"]
    x, T = series(u, "Temperature")
    C_p = F_AR * cp_argon()
    assert np.allclose(T, Ts - (Ts - 300.0) * np.exp(-UA / C_p * x), rtol=0, atol=1e-5)
    assert np.allclose(series(u, "Medium temperature")[1], Ts)
    Q = C_p * (Ts - 300.0) * (1.0 - math.exp(-UA / C_p))
    assert u["wall_duty"] == pytest.approx(Q, rel=1e-7)
    assert u["utility_flow"] == pytest.approx(Q / ut.latent, rel=1e-7)
    assert u["utility_T_out"] == Ts


# ---- a process stream as the medium -------------------------------------------------------

@pytest.mark.parametrize("arrangement", ["cocurrent", "countercurrent"])
def test_process_stream_cools_the_reactor(arrangement):
    counter = arrangement == "countercurrent"
    medium = {"T": 300.0, "P": "3 bar", "flows": {"argon": 6.0}}
    r = run({"medium": "stream", "U": U, "arrangement": arrangement}, medium=medium)
    u = r.units["R"]
    C_p, C_m = F_AR * cp_argon(), 6.0 * cp_argon()
    x, T = series(u, "Temperature")
    _, T_m = series(u, "Medium temperature")
    T_ref, T_m_ref = double_pipe(x, C_p, C_m, T_FEED, 300.0, counter)
    assert np.allclose(T, T_ref, rtol=0, atol=1e-5)
    assert np.allclose(T_m, T_m_ref, rtol=0, atol=1e-5)
    out, m_in, m_out = r.streams["OUT"], r.streams["MIN"], r.streams["MOUT"]
    assert m_out["T"] == pytest.approx(T_m_ref[0] if counter else T_m_ref[-1], abs=1e-5)
    assert m_out["P"] == m_in["P"] and m_out["component_flows"] == m_in["component_flows"]
    # The heat stays inside the unit: what the argon loses, the medium gains.
    Q = out["enthalpy_flow"] - r.streams["IN"]["enthalpy_flow"]
    assert u["wall_duty"] == pytest.approx(Q, rel=1e-9)
    assert m_out["enthalpy_flow"] - m_in["enthalpy_flow"] == pytest.approx(-Q, rel=1e-9)
    assert u["duty"] == 0.0
    energy = next(c for c in r.balance_checks if c["scope"] == "R" and c["quantity"] == "energy")
    assert energy["status"] == "ok" and abs(energy["residual"]) < 1e-6 * abs(Q)


# ---- exothermic reaction: an independent solution of the same balances -------------------

REACTIVE = ["cyclopropane", "propylene", "argon"]
# Cyclopropane isomerizes to propylene, first order, releasing about 33 kJ/mol (from the
# formation enthalpies). The rate constant is made up to give a hot spot in this tube.
K0, EA = 3e11, 150e3
ISOMERIZATION = {"id": "CP", "stoichiometry": {"cyclopropane": -1, "propylene": 1},
                 "kinetics": {"k0": K0, "Ea": "150 kJ/mol"}}
OIL = {"id": "OIL", "type": "fluid", "T_in": "620 K", "T_out": "660 K", "cp": "2300 J/kgK"}
U_R, OIL_FLOW, P_R, T_R = 40.0, 0.05, 2e5, 680.0
N0 = np.array([0.1, 0.0, 0.4])


def reactive(arrangement="countercurrent", **pfr):
    return run({"utility": "OIL", "U": U_R, "utility_flow": OIL_FLOW,
                "arrangement": arrangement, **pfr},
               components=REACTIVE, reactions=(ISOMERIZATION,), utilities=(OIL,),
               feed=dict(zip(REACTIVE, N0)), T=T_R, P=P_R)


@functools.cache
def _enthalpies():
    """Each component's ideal-gas enthalpy with its formation enthalpy, J/mol, as a cubic
    spline in T (0.5 K steps), tabulated once from the pure components."""
    pkg = get_package(tuple(REACTIVE), "ideal")
    T = np.linspace(250.0, 2000.0, 3501)
    splines = []
    for i, hf in enumerate(pkg.consts.Hfgs):
        pure = [1.0 if j == i else 0.0 for j in range(len(REACTIVE))]
        splines.append(CubicSpline(T, [pkg.phase_props("vapor", Tk, P_R, pure).H + hf
                                       for Tk in T]))
    return splines


@functools.cache
def reference(arrangement, tol=1e-8):
    """The balances in another form: the enthalpy flow H (with formation enthalpies) as the
    state, the temperature found from it by Newton on the pure components' tabulated
    enthalpies (an ideal-gas mixture), concentrations from the ideal-gas law. Cocurrent is
    integrated by Radau; countercurrent is solved as the boundary-value problem it is, by
    collocation (solve_bvp) instead of shooting.

    Returns a function of the length fraction (array) giving conversion, T and T_medium.
    """
    h = _enthalpies()
    V = math.pi / 4.0 * D**2 * L
    UA_v = U_R * math.pi * D * L / V
    C = OIL_FLOW * 2300.0

    def state(F_A, H):
        F_A, H = np.asarray(F_A, dtype=float), np.asarray(H, dtype=float)
        n = [F_A, N0[0] - F_A, np.full_like(F_A, N0[2])]
        T = np.full_like(F_A, T_R)
        for _ in range(50):
            step = (sum(ni * hi(T) for ni, hi in zip(n, h)) - H) / sum(
                ni * hi(T, 1) for ni, hi in zip(n, h))
            T = T - step
            if np.max(np.abs(step)) < 1e-10:
                break
        c_A = F_A / N0.sum() * P_R / (R * T)  # total moles do not change
        return T, K0 * np.exp(-EA / (R * T)) * c_A

    H0 = sum(ni * hi(T_R) for ni, hi in zip(N0, h))
    if arrangement == "cocurrent":
        def rhs(v, y):  # F_A, H, heat through the wall so far
            T, r = state(y[0], y[1])
            q = UA_v * (620.0 - y[2] / C - T)
            return [-r, q, q]

        sol = solve_ivp(rhs, (0.0, V), [N0[0], H0, 0.0], method="Radau", rtol=tol,
                        atol=[1e-13, 1e-8, 1e-8], dense_output=True)
        assert sol.success, sol.message

        def at_co(x):
            F_A, H, Q = sol.sol(np.asarray(x) * V)
            return 1.0 - F_A / N0[0], state(F_A, H)[0], 620.0 - Q / C
        return at_co

    def fun(v, y):  # F_A, H, medium temperature; the medium flows towards v = 0
        T, r = state(y[0], y[1])
        q = UA_v * (y[2] - T)
        return np.vstack([-r, q, q / C])

    def bc(ya, yb):
        return np.array([ya[0] - N0[0], ya[1] - H0, yb[2] - 620.0])

    # Start the collocation from the tube with the medium at its inlet temperature all along.
    def rhs_guess(v, y):
        T, r = state(y[0], y[1])
        return [-r, UA_v * (620.0 - T)]

    mesh = np.linspace(0.0, V, 60)
    start = solve_ivp(rhs_guess, (0.0, V), [N0[0], H0], method="Radau", t_eval=mesh,
                      rtol=1e-6, atol=[1e-12, 1e-6])
    guess = np.vstack([start.y, np.full(len(mesh), 620.0)])
    bvp = solve_bvp(fun, bc, mesh, guess, tol=tol, max_nodes=100000)
    assert bvp.success, bvp.message

    def at_counter(x):
        F_A, H, T_m = bvp.sol(np.asarray(x) * V)
        return 1.0 - F_A / N0[0], state(F_A, H)[0], T_m
    return at_counter


@pytest.mark.parametrize("arrangement", ["cocurrent", "countercurrent"])
def test_exothermic_reaction_matches_an_independent_solution(arrangement):
    u = reactive(arrangement).units["R"]
    ref = reference(arrangement)
    x, T = series(u, "Temperature")
    _, T_m = series(u, "Medium temperature")
    _, X = series(u, "Conversion of cyclopropane")
    X_ref, T_ref, T_m_ref = ref(x)
    assert np.allclose(T, T_ref, rtol=0, atol=2e-4)
    assert np.allclose(T_m, T_m_ref, rtol=0, atol=2e-4)
    assert np.allclose(X, X_ref, rtol=0, atol=1e-6)
    assert u["conversion"] == pytest.approx(X_ref[-1], abs=1e-6)
    # A hot spot inside the tube, where the reference has its maximum too.
    fine = np.linspace(0.0, 1.0, 4001)
    T_fine = ref(fine)[1]
    k = int(np.argmax(T_fine))
    assert 0 < k < len(fine) - 1 and T_fine[k] > T_R + 1.0
    assert u["T_hot_spot"] == pytest.approx(T_fine[k], abs=2e-4)
    assert u["hot_spot_position"] == pytest.approx(fine[k] * L, abs=2 * L / 4000)
    if arrangement == "countercurrent":
        assert T_m[-1] == pytest.approx(620.0, abs=1e-6)  # the medium inlet condition
    # Overall energy balance: the oil takes what the reacting gas gives off.
    assert u["wall_duty"] == pytest.approx(-OIL_FLOW * 2300.0 * (u["utility_T_out"] - 620.0),
                                           rel=1e-9)


def test_tolerance_refinement_converges(monkeypatch):
    """The integration is adaptive: tightening its tolerance converges, and the solution
    does not depend on how many points the profile reports."""
    results, profile = [], None
    for rtol in (1e-4, 1e-6, 1e-8, 1e-10):
        monkeypatch.setattr(PFR, "RTOL", rtol)
        u = reactive().units["R"]
        results.append(np.array([u["T_out"], u["T_hot_spot"], u["utility_T_out"]]))
        if rtol == 1e-8:
            profile = series(u, "Temperature")[1]
    errors = [np.max(np.abs(r - results[-1])) for r in results[:-1]]
    assert errors[0] > errors[1] > errors[2] and errors[2] < 1e-5
    monkeypatch.setattr(PFR, "RTOL", 1e-8)
    monkeypatch.setattr(PFR, "N_POINTS", 4 * (PFR.N_POINTS - 1) + 1)
    fine = reactive().units["R"]
    assert fine["T_out"] == pytest.approx(results[2][0], abs=1e-9)
    assert series(fine, "Temperature")[1][::4] == pytest.approx(profile, abs=1e-6)


def test_collocation_refinement_converges_to_the_shooting_solution():
    u = reactive().units["R"]
    misses = [abs(reference("countercurrent", tol=tol)(1.0)[1] - u["T_out"])
              for tol in (1e-3, 1e-5, 1e-8)]
    assert misses[0] > misses[2] and misses[2] < 1e-4


# ---- the other modes, and mistakes --------------------------------------------------------

def test_with_a_temperature_the_utility_reports_its_consumption():
    """Isothermal with a utility: nothing is rated along the wall; the utility's flow
    follows from the duty, as for any unit, and the reactor has no medium ports."""
    r = run({"utility": "CW", "T": 500.0})
    u = r.units["R"]
    assert u["duty"] == pytest.approx(F_AR * cp_argon() * (500.0 - T_FEED), rel=1e-9)
    assert u["utility_flow"] == pytest.approx(-u["duty"] / (CP_WATER * 20.0), rel=1e-9)
    assert "wall_duty" not in u and "MOUT" not in r.streams


def test_with_a_duty_the_medium_stream_takes_it_up():
    r = run({"medium": "stream", "duty": -5000.0},
            medium={"T": 300.0, "P": "3 bar", "flows": {"argon": 6.0}})
    u = r.units["R"]
    assert u["wall_duty"] == pytest.approx(-5000.0) and u["duty"] == 0.0
    assert r.streams["MOUT"]["T"] == pytest.approx(300.0 + 5000.0 / (6.0 * cp_argon()))


@pytest.mark.parametrize("pfr, message", [
    ({"medium": "stream", "utility": "CW", "U": U}, "utility is chosen too"),
    ({"medium": "stream"}, "which needs U"),
    ({"utility": "CW", "U": U, "V": 0.02, "length": None}, "needs its length"),
    ({"utility": "LPS", "U": U}, "can only heat"),  # argon at 600 K is hotter than the steam
    ({"utility": "CW", "U": U, "utility_flow": 0.0}, "flow must be positive"),
])
def test_mistakes_are_explained(pfr, message):
    pfr = {"length": L, **pfr}
    units = [{"id": "F", "type": "Feed", "params": {"T": T_FEED, "P": "2 bar",
                                                    "flows": {"argon": F_AR}}},
             {"id": "R", "type": "PFR", "params": {
                 "diameter": D, "phase": "gas", "reactions": ["ISO"],
                 **{k: v for k, v in pfr.items() if v is not None}}}]
    r = solve({"components": INERT, "property_method": "ideal", "reactions": [NO_RATE],
               "utilities": [CW, LPS], "units": units,
               "streams": [{"id": "IN", "from": "F.out", "to": "R.in"},
                           {"id": "OUT", "from": "R.out"}]})
    assert not r.ok and message in r.error["message"]


def test_a_cooling_utility_colder_than_nothing_is_a_temperature_cross():
    r = solve({"components": INERT, "property_method": "ideal", "reactions": [NO_RATE],
               "utilities": [CW],
               "units": [{"id": "F", "type": "Feed",
                          "params": {"T": 290.0, "P": "2 bar", "flows": {"argon": F_AR}}},
                         {"id": "R", "type": "PFR",
                          "params": {"length": L, "diameter": D, "phase": "gas",
                                     "reactions": ["ISO"], "utility": "CW", "U": U}}],
               "streams": [{"id": "IN", "from": "F.out", "to": "R.in"},
                           {"id": "OUT", "from": "R.out"}]})
    assert not r.ok and "temperature cross" in r.error["message"]


def test_a_medium_stream_needs_its_ports_connected():
    r = solve({"components": INERT, "property_method": "ideal", "reactions": [NO_RATE],
               "units": [{"id": "F", "type": "Feed",
                          "params": {"T": T_FEED, "P": "2 bar", "flows": {"argon": F_AR}}},
                         {"id": "R", "type": "PFR",
                          "params": {"length": L, "diameter": D, "phase": "gas",
                                     "reactions": ["ISO"], "medium": "stream", "U": U}}],
               "streams": [{"id": "IN", "from": "F.out", "to": "R.in"},
                           {"id": "OUT", "from": "R.out"}]})
    assert not r.ok and "coolant_in" in r.error["message"]
    schema = PFR.schema()
    assert [p["name"] for p in schema["inlets"]] == ["in", "coolant_in"]
    assert [p["name"] for p in schema["outlets"]] == ["out", "coolant_out"]
