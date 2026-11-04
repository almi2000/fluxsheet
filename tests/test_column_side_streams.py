"""Side draws, pumparounds, side heaters and reactive stages of the rigorous column."""

import json
from dataclasses import replace
from itertools import pairwise
from pathlib import Path

import numpy as np
import pytest

from fluxsheet import solve
from fluxsheet.columns import ColumnSetup, Pumparound, SideDraw, solve_column, stages
from fluxsheet.reactions import Kinetics, RateLaw, Reaction
from fluxsheet.stream import Stream
from fluxsheet.thermo_package import get_package
from fluxsheet.units.columns import _feed

KMOLH = 1 / 3.6  # mol/s
EXAMPLES = Path(__file__).parent.parent / "examples"
SIDE = ["side1", "side2", "side3", "side4"]
BTX = {"benzene": "40 kmol/h", "toluene": "30 kmol/h", "o-xylene": "30 kmol/h"}
BASE = {"stages": 20, "feed_stage": 10, "reflux_ratio": 2, "distillate_rate": "40 kmol/h"}


def column(params, *, flows=BTX, components=("benzene", "toluene", "o-xylene"),
           method="peng-robinson", reactions=None, T="107 degC"):
    n_draws = len(params.get("side_draws", []))
    return solve({
        "components": list(components),
        "property_method": method,
        "reactions": reactions or [],
        "units": [{"id": "F", "type": "Feed", "params": {"T": T, "P": "1.5 bar", "flows": flows}},
                  {"id": "C", "type": "DistillationColumn", "params": params}],
        "streams": [{"id": "FEED", "from": "F.out", "to": "C.in"},
                    {"id": "distillate", "from": "C.distillate"},
                    {"id": "bottoms", "from": "C.bottoms"}]
        + [{"id": port, "from": f"C.{port}"} for port in SIDE[:n_draws]],
    })


def closure(r, *, formation=None):
    """Largest relative error of the component balances and of the energy balance, with all
    duties of the column (condenser, reboiler, pumparounds, side heaters). With reactions,
    ``formation`` (J/mol per component) adds the formation enthalpies on both sides."""
    s, u = r.streams, r.units["C"]
    outlets = [k for k in s if k != "FEED"]
    names = list(s["FEED"]["component_flows"])
    flow_in = np.array([s["FEED"]["component_flows"][c] for c in names])
    flow_out = np.array([sum(s[o]["component_flows"][c] for o in outlets) for c in names])
    duties = (u["condenser_duty"] + u["reboiler_duty"] + (u.get("pumparound_duty") or 0.0)
              + sum(d["duty"] for d in u["params_stage_duties"]))
    h_in = s["FEED"]["enthalpy_flow"]
    h_out = sum(s[o]["enthalpy_flow"] for o in outlets)
    if formation is not None:
        h_in += flow_in @ formation
        h_out += flow_out @ formation
    scale = max(abs(u["condenser_duty"]), abs(u["reboiler_duty"]), abs(h_in), abs(h_out))
    comp = np.max(np.abs(flow_in - flow_out)) / flow_in.sum()
    return comp, abs(h_in + duties - h_out) / scale


def with_duties(r, params):
    """Attach the side heater duties (W) to the unit results for closure()."""
    from fluxsheet.quantities import to_si

    r.units["C"]["params_stage_duties"] = [
        {"duty": to_si(d["duty"], "power")} for d in params.get("stage_duties", [])]
    return r


# --- side draws ---------------------------------------------------------------------------

def test_side_draws_close_the_balances_and_leave_in_their_phase():
    params = {**BASE, "side_draws": [
        {"stage": 6, "phase": "liquid", "rate": "8 kmol/h"},
        {"stage": 15, "phase": "vapor", "fraction": 0.1}]}
    r = with_duties(column(params), params)
    assert r.ok, r.error
    assert not r.warnings
    comp, energy = closure(r)
    assert comp < 1e-8 and energy < 1e-7
    s = r.streams
    assert s["side1"]["molar_flow"] == pytest.approx(8 * KMOLH, rel=1e-8)
    assert s["side1"]["phase"] == "liquid" and s["side2"]["phase"] == "vapor"
    # The vapor draw takes a tenth of the vapor leaving stage 15 upward and sideways.
    V = r.units["C"]["profiles"]["flows"]["series"][1]["values"]
    assert s["side2"]["molar_flow"] == pytest.approx(0.1 * (V[14] + s["side2"]["molar_flow"]),
                                                     rel=1e-8)
    # A middle-stage draw is richest in the middle component.
    x = s["side1"]["mole_fractions"]
    assert x["toluene"] > s["distillate"]["mole_fractions"]["toluene"]
    assert x["toluene"] > s["bottoms"]["mole_fractions"]["toluene"]
    # The distillate rate holds; the draws come out of the bottoms.
    assert s["distillate"]["molar_flow"] == pytest.approx(40 * KMOLH, rel=1e-8)


@pytest.mark.parametrize("draw", [{"rate": 0.0}, {"fraction": 0.0}])
def test_a_zero_side_draw_reproduces_the_plain_column(draw):
    plain = column(BASE)
    r = column({**BASE, "side_draws": [{"stage": 8, "phase": "liquid", **draw}]})
    assert r.ok, r.error
    assert r.streams["side1"]["molar_flow"] == 0.0
    for port in ("distillate", "bottoms"):
        for c, f in plain.streams[port]["component_flows"].items():
            assert r.streams[port]["component_flows"][c] == pytest.approx(f, rel=1e-8, abs=1e-12)
    assert r.units["C"]["reboiler_duty"] == pytest.approx(plain.units["C"]["reboiler_duty"],
                                                          rel=1e-8)


def test_unused_side_draw_outlets_raise_no_warning():
    r = column(BASE)
    assert r.ok and not r.warnings


# --- pumparounds and side heaters ---------------------------------------------------------

@pytest.fixture(scope="module")
def btx():
    pkg = get_package(("benzene", "toluene", "o-xylene"), "peng-robinson")
    F = 100 * KMOLH
    feed = Stream.from_flash(pkg, [0.4 * F, 0.3 * F, 0.3 * F], T=380.0, P=1.5e5)
    setup = ColumnSetup(stages=20, P=np.full(20, 1.5e5), feeds=[_feed(9, feed)],
                        specs={"reflux_ratio": 2.0, "distillate_rate": 40 * KMOLH})
    return pkg, setup, solve_column(pkg, setup)


def test_an_adiabatic_pumparound_to_the_stage_below_changes_nothing(btx):
    """Liquid pumped from stage a to a + 1 without heat exchange only joins the liquid that
    flows there anyway: the steady state is that of the plain column."""
    pkg, setup, plain = btx
    sol = solve_column(pkg, replace(setup, pumparounds=[Pumparound(5, 6, rate=20 * KMOLH)]))
    assert np.max(np.abs(sol.T - plain.T)) < 1e-6
    assert np.allclose(sol.distillate, plain.distillate, rtol=1e-8)
    assert sol.reboiler_duty == pytest.approx(plain.reboiler_duty, rel=1e-8)
    assert sol.pumparound_duties == [0.0]
    assert sol.L[5] == pytest.approx(plain.L[5] - 20 * KMOLH, rel=1e-6)


def test_a_pumparound_duty_acts_like_a_side_cooler_on_its_return_stage(btx):
    pkg, setup, _ = btx
    pa = solve_column(pkg, replace(setup, pumparounds=[
        Pumparound(5, 6, fraction=0.3, duty=-1e5)]))
    cooler = solve_column(pkg, replace(setup, stage_duties={6: -1e5}))
    assert np.max(np.abs(pa.T - cooler.T)) < 1e-6
    assert pa.reboiler_duty == pytest.approx(cooler.reboiler_duty, rel=1e-8)
    assert pa.pumparound_duties[0] == pytest.approx(-1e5)


@pytest.mark.parametrize("pa", [
    {"from_stage": 13, "to_stage": 4, "rate": "30 kmol/h", "duty": "-300 kW"},
    {"from_stage": 13, "to_stage": 4, "fraction": 0.3, "T": "60 degC"},
    {"from_stage": 20, "to_stage": 1, "rate": "5 kmol/h"},
], ids=["duty", "temperature", "reboiler-to-condenser"])
def test_pumparounds_close_the_balances(pa):
    params = {**BASE, "pumparounds": [pa],
              "stage_duties": [{"stage": 16, "duty": "50 kW"}],
              "side_draws": [{"stage": 8, "phase": "liquid", "rate": "5 kmol/h"}]}
    r = with_duties(column(params), params)
    assert r.ok, r.error
    comp, energy = closure(r)
    assert comp < 1e-8 and energy < 1e-7
    res = r.units["C"]["pumparounds"][0]
    assert res["from_stage"] == pa["from_stage"] and res["to_stage"] == pa["to_stage"]
    if "T" in pa:
        assert res["T_return"] == pytest.approx(333.15)
        assert res["duty"] < 0.0  # cooled from the stage temperature
    if "duty" in pa:
        assert res["duty"] == pytest.approx(-3e5)
        T = r.units["C"]["profiles"]["temperature"]["series"][0]["values"]
        assert res["T_return"] < T[12]


def test_a_pumparound_cooler_lowers_the_condenser_duty():
    """Heat taken out by a pumparound in the upper column is heat the condenser no longer
    has to remove (at fixed reflux the internal reflux below the return grows instead)."""
    plain = column(BASE)
    r = column({**BASE, "pumparounds": [
        {"from_stage": 6, "to_stage": 3, "rate": "40 kmol/h", "duty": "-200 kW"}]})
    assert r.ok, r.error
    assert r.units["C"]["reboiler_duty"] > plain.units["C"]["reboiler_duty"]
    L_plain = plain.units["C"]["profiles"]["flows"]["series"][0]["values"]
    L = r.units["C"]["profiles"]["flows"]["series"][0]["values"]
    assert L[3] > L_plain[3]


# --- Jacobian -------------------------------------------------------------------------------

def _ester_reaction(pkg):
    nu = np.array([-1.0, -1.0, 1.0, 1.0])
    law = RateLaw(k0=2e-7, Ea=5e4, orders=(1.0, 1.0, 0.0, 0.0), T_ref=343.15)
    return Reaction("EST", nu, tuple(pkg.names), Kinetics(law, "concentration",
                                                          equilibrium_limited=True),
                    keq=(1.386, 0.0, 0.0, 0.0))


@pytest.fixture(scope="module")
def ester():
    pkg = get_package(("acetic acid", "ethanol", "ethyl acetate", "water"), "nrtl")
    feed = Stream.from_flash(pkg, [10 * KMOLH, 10 * KMOLH, 0.0, 0.0], T=343.15, P=1.1e5)
    hf = np.array(pkg.consts.Hfgs, dtype=float)
    return pkg, feed, hf


def test_colored_jacobian_equals_the_dense_one(ester):
    """The coloring must see every coupling: a pumparound from the reboiler to the condenser,
    one from the middle, a rate-specified draw, a stage duty and reactive stages."""
    pkg, feed, hf = ester
    N = 7
    setup = ColumnSetup(
        stages=N, P=np.full(N, 1.1e5), feeds=[_feed(3, feed)],
        specs={"reflux_ratio": 2.0, "distillate_rate": 10 * KMOLH},
        draws=[SideDraw(2, "liquid", rate=1 * KMOLH), SideDraw(4, "vapor", fraction=0.1)],
        pumparounds=[Pumparound(N - 1, 0, rate=2 * KMOLH),
                     Pumparound(4, 1, fraction=0.2, duty=-2e4),
                     Pumparound(3, 5, rate=1 * KMOLH, T=320.0)],
        stage_duties={2: 1e4}, reactions=[_ester_reaction(pkg)],
        holdup={1: 0.05, 2: 0.05, 3: 0.05}, formation_enthalpies=hf)
    m = stages._Model(pkg, setup)
    x = stages._initial(m)
    props = m.all_props(x)
    r0 = m.residual(x, props)
    J = m.jacobian(x, props, r0).toarray()
    steps = m._steps(x)
    dense = np.empty_like(J)
    for col in range(m.size):
        xp = x.copy()
        xp[col] += steps[col]
        dense[:, col] = (m.residual(xp, m.all_props(xp)) - r0) / steps[col]
    assert np.max(np.abs(J - dense)) < 1e-9 * np.max(np.abs(dense))


# --- reactive stages --------------------------------------------------------------------------

REACTIONS = json.loads((EXAMPLES / "ethyl_acetate.json").read_text())["reactions"]
ESTER = {"acetic acid": "10 kmol/h", "ethanol": "10 kmol/h"}
ESTER_COMPONENTS = ("acetic acid", "ethanol", "ethyl acetate", "water")


def ester_column(holdup: float, **extra):
    params = {"stages": 14, "feed_stage": 4, "reflux_ratio": 2, "distillate_rate": "10 kmol/h",
              **extra}
    if holdup:
        params |= {"reactions": ["EST"],
                   "reactive_stages": [{"stage": j, "holdup": f"{holdup} m3"}
                                       for j in range(3, 11)]}
    return column(params, flows=ESTER, components=ESTER_COMPONENTS, method="nrtl",
                  reactions=REACTIONS, T="70 degC")


def test_reactive_stages_conserve_elements_and_energy(ester):
    pkg, _, hf = ester
    r = with_duties(ester_column(0.05), {})
    assert r.ok, r.error
    s = r.streams
    names = list(s["FEED"]["component_flows"])
    outlets = ("distillate", "bottoms")
    for el in ("C", "H", "O"):
        atoms = np.array([pkg.consts.atomss[pkg.names.index(c)].get(el, 0) for c in names])
        fed = atoms @ [s["FEED"]["component_flows"][c] for c in names]
        out = sum(atoms @ [s[o]["component_flows"][c] for c in names] for o in outlets)
        assert out == pytest.approx(fed, rel=1e-9)
    formed = np.sum(r.units["C"]["profiles"]["reaction"]["series"][2]["values"])
    assert formed == pytest.approx(s["distillate"]["component_flows"]["ethyl acetate"]
                                   + s["bottoms"]["component_flows"]["ethyl acetate"], rel=1e-8)
    assert formed > 1 * KMOLH
    _, energy = closure(r, formation=hf)
    assert energy < 1e-7
    # Without the heat of reaction the sensible enthalpies do not close.
    assert closure(r)[1] > 10 * energy


def test_conversion_rises_with_the_hold_up():
    conversion = []
    for holdup in (0.0, 0.01, 0.05, 0.2):
        r = ester_column(holdup)
        assert r.ok, r.error
        s = r.streams
        acid = sum(s[o]["component_flows"]["acetic acid"] for o in ("distillate", "bottoms"))
        conversion.append(1.0 - acid / s["FEED"]["component_flows"]["acetic acid"])
    assert conversion[0] == pytest.approx(0.0, abs=1e-9)
    assert all(b > a + 0.02 for a, b in pairwise(conversion))
    # The lightest product, ethyl acetate, goes overhead; acetic acid, the heaviest, down.
    assert s["distillate"]["mole_fractions"]["ethyl acetate"] > 0.4
    assert s["bottoms"]["mole_fractions"]["acetic acid"] > s["distillate"]["mole_fractions"][
        "acetic acid"]


def test_single_reactive_stage_is_a_cstr_at_stage_conditions(ester):
    """One stage without condenser and reboiler: a flash drum that reacts its liquid like a
    CSTR of the hold-up volume at the stage's temperature and liquid composition."""
    pkg, feed, hf = ester
    rxn = _ester_reaction(pkg)
    hot = Stream.from_flash(pkg, feed.flows, P=1.1e5, VF=0.4)
    setup = ColumnSetup(stages=1, P=np.array([1.1e5]), feeds=[_feed(0, hot)],
                        condenser="none", reboiler=False, reactions=[rxn], holdup={0: 0.3},
                        formation_enthalpies=hf)
    sol = solve_column(pkg, setup)
    x, T = sol.x[0], sol.T[0]
    Vm = pkg.phase_props("liquid", T, 1.1e5, list(x)).V
    rate = rxn.rate(pkg, T, 1.1e5, x, Vm, "liquid")
    assert rate > 0.0
    formed = sol.bottoms + sol.vapor_distillate - np.asarray(hot.flows)
    assert np.allclose(formed, 0.3 * rate * rxn.nu, rtol=1e-8, atol=1e-12)
    # Adiabatic: sensible enthalpy out minus in is the heat of reaction released.
    H_out = sol.L[0] * sol.h[0] + sol.V[0] * sol.Hv[0]
    assert H_out - hot.enthalpy_flow == pytest.approx(-formed @ hf, rel=1e-7)


# --- input checks and schema ------------------------------------------------------------------

@pytest.mark.parametrize("params, message", [
    ({"side_draws": [{"stage": 1, "rate": "1 kmol/h"}]}, "on the condenser"),
    ({"side_draws": [{"stage": 20, "rate": "1 kmol/h"}]}, "on the reboiler"),
    ({"side_draws": [{"stage": 6}]}, "give either rate or fraction"),
    ({"side_draws": [{"stage": 6, "rate": "1 kmol/h", "fraction": 0.1}]}, "either rate or"),
    ({"side_draws": [{"stage": 30, "rate": "1 kmol/h"}]}, "outside the column"),
    ({"side_draws": [{"stage": 6, "rate": "1 kmol/h"}] * 5}, "outlets for 4"),
    ({"side_draws": [{"stage": 6, "flow": 1}]}, "unknown field 'flow'"),
    ({"side_draws": [{"stage": 6, "fraction": 0.6}, {"stage": 6, "fraction": 0.5}]},
     "adds up to a fraction"),
    ({"pumparounds": [{"from_stage": 5, "to_stage": 5, "rate": "1 kmol/h"}]},
     "returns to the stage it is drawn from"),
    ({"pumparounds": [{"from_stage": 8, "to_stage": 5, "rate": "1 kmol/h", "duty": "-1 kW",
                       "T": "50 degC"}]}, "at most one of duty or T"),
    ({"reactive_stages": [{"stage": 5, "holdup": "0.1 m3"}]}, "no reactions"),
    ({"reactions": ["R1"]}, "no reactive stages"),
    ({"side_draws": [{"stage": 6, "rate": "200 kmol/h"}]}, "cannot be met"),
])
def test_bad_side_streams_explain_themselves(params, message):
    r = column({**BASE, **params})
    assert not r.ok
    assert message in r.error["message"]
    assert r.error["hint"] or "field" in message  # unit conversion errors carry no hint


def test_records_parameters_describe_their_fields():
    from fluxsheet.units.columns import DistillationColumn

    schema = DistillationColumn.schema()
    params = {p["name"]: p for p in schema["params"]}
    fields = {f["name"]: f for f in params["side_draws"]["fields"]}
    assert params["side_draws"]["quantity"] == "records"
    assert fields["phase"]["options"] == ["liquid", "vapor"]
    assert fields["rate"]["quantity"] == "molar_flow"
    assert params["reflux_ratio"]["fields"] is None
    assert [o["name"] for o in schema["outlets"]][-4:] == SIDE


def test_a_liquid_fraction_draw_takes_its_share_of_the_stage_outflow():
    r = column({**BASE, "side_draws": [{"stage": 6, "phase": "liquid", "fraction": 0.25}]})
    assert r.ok, r.error
    L = r.units["C"]["profiles"]["flows"]["series"][0]["values"]
    W = r.streams["side1"]["molar_flow"]
    assert W == pytest.approx(0.25 * (L[5] + W), rel=1e-8)


def test_warm_start_with_side_streams_gives_the_same_column(btx):
    pkg, setup, _ = btx
    side = replace(setup, draws=[SideDraw(14, "liquid", rate=10 * KMOLH)],
                   pumparounds=[Pumparound(8, 4, fraction=0.3, duty=-5e4)])
    cold = solve_column(pkg, side)
    warm = solve_column(pkg, side, start=cold)
    assert warm.iterations < cold.iterations
    assert np.max(np.abs(warm.T - cold.T)) < 1e-6
    assert np.allclose(warm.side_draws[0], cold.side_draws[0], rtol=1e-7)
    # A draw that changes phase is another column: no warm start from the old one.
    vapor = replace(side, draws=[SideDraw(14, "vapor", rate=10 * KMOLH)])
    assert stages._restart_key(stages._Model(pkg, vapor)) != cold.restart[0]


def test_a_pumparound_to_the_condenser_keeps_the_reflux_spec():
    r = column({**BASE, "pumparounds": [{"from_stage": 8, "to_stage": 1, "rate": "20 kmol/h"}]})
    assert r.ok, r.error
    assert r.units["C"]["reflux_ratio"] == pytest.approx(2.0, rel=1e-8)
    assert closure(with_duties(r, {}))[1] < 1e-7


@pytest.mark.parametrize("params, message", [
    ({"stage_duties": [{"stage": 1, "duty": "10 kW"}]}, "on the condenser"),
    ({"stage_duties": [{"stage": 20, "duty": "10 kW"}]}, "on the reboiler"),
    ({"pumparounds": [{"from_stage": 1, "to_stage": 4, "fraction": 0.2}]},
     "fraction of the condenser's outflow"),
    ({"pumparounds": [{"from_stage": 8, "to_stage": 4, "rate": 0, "duty": "-10 kW"}]},
     "no flow to heat or cool"),
])
def test_more_bad_side_streams_explain_themselves(params, message):
    r = column({**BASE, **params})
    assert not r.ok
    assert message in r.error["message"] and r.error["hint"]
