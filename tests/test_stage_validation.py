"""Stage balances reconstructed from the public column solution."""

from dataclasses import replace

import numpy as np
import pytest

from fluxsheet.columns import ColumnSetup, Pumparound, SideDraw, solve_column
from fluxsheet.stream import Stream
from fluxsheet.thermo_package import get_package
from fluxsheet.units.columns import _feed


def stage_balances(pkg, setup, sol):
    """Use reported flows and properties, without the solver residual function."""
    n = setup.stages
    component = np.zeros((n, pkg.n))
    energy = np.zeros(n)
    for j in range(n):
        component[j] -= sol.L[j] * sol.x[j] + sol.V[j] * sol.y[j]
        energy[j] -= sol.L[j] * sol.h[j] + sol.V[j] * sol.Hv[j]
        if j:
            component[j] += sol.L[j - 1] * sol.x[j - 1]
            energy[j] += sol.L[j - 1] * sol.h[j - 1]
        if j + 1 < n:
            component[j] += sol.V[j + 1] * sol.y[j + 1]
            energy[j] += sol.V[j + 1] * sol.Hv[j + 1]
    for feed in setup.feeds:
        component[feed.stage] += feed.flows
        energy[feed.stage] += feed.H
    if setup.condenser == "total":
        component[0] -= sol.distillate
        energy[0] -= sol.distillate.sum() * sol.h[0]
    energy[0] += sol.condenser_duty
    energy[-1] += sol.reboiler_duty
    for draw, flows in zip(setup.draws, sol.side_draws):
        component[draw.stage] -= flows
        enthalpy = sol.h[draw.stage] if draw.phase == "liquid" else sol.Hv[draw.stage]
        energy[draw.stage] -= flows.sum() * enthalpy
    for pa, flows, returned_H in zip(
            setup.pumparounds, sol.pumparound_flows, sol.pumparound_return_H):
        component[pa.source] -= flows
        energy[pa.source] -= flows.sum() * sol.h[pa.source]
        component[pa.target] += flows
        energy[pa.target] += returned_H
    for stage, duty in setup.stage_duties.items():
        energy[stage] += duty
    if sol.generation is not None:
        component += sol.generation
        energy -= sol.generation @ setup.formation_enthalpies
    return component, energy


@pytest.fixture(scope="module")
def reference():
    pkg = get_package(("benzene", "toluene", "o-xylene"), "peng-robinson")
    feed = Stream.from_flash(pkg, [11.111111, 8.333333, 8.333333], T=380.0, P=1.5e5)
    setup = ColumnSetup(stages=20, P=np.full(20, 1.5e5), feeds=[_feed(9, feed)],
                        specs={"reflux_ratio": 2.0, "distillate_rate": 40 / 3.6})
    return pkg, setup


@pytest.mark.parametrize("condenser,efficiency", [("total", 1.0), ("partial", 1.0),
                                                 ("total", 0.7)])
def test_independent_stage_balances_and_equilibrium(reference, condenser, efficiency):
    pkg, base = reference
    setup = replace(base, condenser=condenser, efficiency=efficiency)
    sol = solve_column(pkg, setup)
    component, energy = stage_balances(pkg, setup, sol)
    assert np.max(np.abs(component)) < 1e-7
    assert np.max(np.abs(energy)) < 1e-2
    for j in range(setup.stages):
        K = np.array(pkg.k_values(sol.T[j], sol.P[j], list(sol.x[j]), list(sol.y[j])))
        if j == 0 and condenser == "total":
            assert sum(K * sol.x[j]) == pytest.approx(1.0, abs=1e-6)
            continue
        eta = 1.0 if j in (0, setup.stages - 1) else efficiency
        incoming_y = sol.y[j + 1] if j + 1 < setup.stages else np.zeros(pkg.n)
        assert np.allclose(sol.y[j], eta * K * sol.x[j] + (1 - eta) * incoming_y,
                           rtol=1e-6, atol=1e-7), j


def test_independent_stage_balances_with_side_streams_and_multiple_feeds(reference):
    pkg, base = reference
    second = Stream.from_flash(pkg, [0.3, 0.2, 0.1], T=390.0, P=1.5e5)
    setup = replace(base, feeds=[base.feeds[0], _feed(13, second)],
                    draws=[SideDraw(5, "liquid", rate=1.0),
                           SideDraw(14, "vapor", fraction=0.1)],
                    pumparounds=[Pumparound(7, 3, rate=3.0, duty=-1e4)],
                    stage_duties={15: 5000.0})
    sol = solve_column(pkg, setup)
    component, energy = stage_balances(pkg, setup, sol)
    assert np.max(np.abs(component)) < 1e-7
    assert np.max(np.abs(energy)) < 1e-2


def test_reactive_stage_independent_component_and_energy_balances():
    from fluxsheet.reactions import Kinetics, RateLaw, Reaction

    pkg = get_package(("acetic acid", "ethanol", "ethyl acetate", "water"), "nrtl")
    flows = [10 / 3.6, 10 / 3.6, 0.0, 0.0]
    feed = Stream.from_flash(pkg, flows, P=1.1e5, VF=0.4)
    law = RateLaw(k0=2e-7, Ea=5e4, orders=(1.0, 1.0, 0.0, 0.0), T_ref=343.15)
    reaction = Reaction("EST", np.array([-1.0, -1.0, 1.0, 1.0]), tuple(pkg.names),
                        Kinetics(law, "concentration", equilibrium_limited=True),
                        keq=(1.386, 0.0, 0.0, 0.0))
    setup = ColumnSetup(stages=1, P=np.array([1.1e5]), feeds=[_feed(0, feed)],
                        condenser="none", reboiler=False, reactions=[reaction],
                        holdup={0: 0.3}, formation_enthalpies=np.array(pkg.consts.Hfgs))
    sol = solve_column(pkg, setup)
    component, energy = stage_balances(pkg, setup, sol)
    assert np.max(np.abs(component)) < 1e-7
    assert np.max(np.abs(energy)) < 1e-2
    assert sol.generation[0, 2] > 0


def test_absorber_stage_balances_and_equilibrium():
    pkg = get_package(("methane", "ethane", "propane", "n-decane"), "peng-robinson")
    gas = Stream.from_flash(pkg, [85 / 3.6, 10 / 3.6, 5 / 3.6, 0.0],
                            T=303.15, P=30e5)
    oil = Stream.from_flash(pkg, [0.0, 0.0, 0.0, 40 / 3.6], T=303.15, P=30e5)
    setup = ColumnSetup(stages=8, P=np.full(8, 30e5),
                        feeds=[_feed(0, oil), _feed(7, gas)],
                        condenser="none", reboiler=False)
    sol = solve_column(pkg, setup)
    component, energy = stage_balances(pkg, setup, sol)
    assert np.max(np.abs(component)) < 1e-7
    assert np.max(np.abs(energy)) < 1e-2
    for j in range(setup.stages):
        K = np.array(pkg.k_values(sol.T[j], sol.P[j], list(sol.x[j]), list(sol.y[j])))
        assert np.allclose(sol.y[j], K * sol.x[j], rtol=1e-6, atol=1e-7)


def test_ideal_binary_single_stage_matches_raoult_rachford_rice():
    from scipy.optimize import brentq

    pkg = get_package(("benzene", "toluene"), "ideal")
    T, P = 365.0, 1e5
    composition = np.array([0.5, 0.5])
    K = np.array([correlation(T) / P for correlation in pkg.props.VaporPressures])
    vapor_fraction = brentq(
        lambda beta: sum(composition * (K - 1) / (1 + beta * (K - 1))), 0.0, 1.0)
    liquid_x = composition / (1 + vapor_fraction * (K - 1))
    vapor_y = K * liquid_x
    feed = Stream.from_flash(pkg, composition.tolist(), T=T, P=P)
    setup = ColumnSetup(stages=1, P=np.array([P]), feeds=[_feed(0, feed)],
                        condenser="none", reboiler=False)
    sol = solve_column(pkg, setup)
    assert sol.T[0] == pytest.approx(T, abs=1e-4)
    assert sol.V[0] == pytest.approx(vapor_fraction, abs=1e-6)
    assert np.allclose(sol.x[0], liquid_x, rtol=1e-6, atol=1e-7)
    assert np.allclose(sol.y[0], vapor_y, rtol=1e-6, atol=1e-7)


def test_multistage_reactive_column_stage_balances():
    from fluxsheet.reactions import Kinetics, RateLaw, Reaction

    pkg = get_package(("acetic acid", "ethanol", "ethyl acetate", "water"), "nrtl")
    feed = Stream.from_flash(pkg, [10 / 3.6, 10 / 3.6, 0.0, 0.0], T=343.15, P=1.1e5)
    law = RateLaw(k0=2e-7, Ea=5e4, orders=(1.0, 1.0, 0.0, 0.0), T_ref=343.15)
    reaction = Reaction("EST", np.array([-1.0, -1.0, 1.0, 1.0]), tuple(pkg.names),
                        Kinetics(law, "concentration", equilibrium_limited=True),
                        keq=(1.386, 0.0, 0.0, 0.0))
    setup = ColumnSetup(
        stages=14, P=np.full(14, 1.1e5), feeds=[_feed(3, feed)],
        specs={"reflux_ratio": 2.0, "distillate_rate": 10 / 3.6},
        reactions=[reaction], holdup={j: 0.05 for j in range(2, 10)},
        formation_enthalpies=np.array(pkg.consts.Hfgs))
    sol = solve_column(pkg, setup)
    component, energy = stage_balances(pkg, setup, sol)
    assert np.max(np.abs(component)) < 1e-7
    assert np.max(np.abs(energy)) < 1e-2
    assert sum(sol.generation[:, 2]) > 0.1
    for stage, volume in setup.holdup.items():
        molar_volume = pkg.phase_props("liquid", sol.T[stage], sol.P[stage],
                                       list(sol.x[stage])).V
        rate = reaction.rate(pkg, sol.T[stage], sol.P[stage], sol.x[stage],
                             molar_volume, "liquid")
        assert np.allclose(sol.generation[stage], volume * rate * reaction.nu,
                           rtol=1e-8, atol=1e-10)
