"""A gas-phase kinetic reactor at high pressure, checked against an independent solution.

Toluene hydrodealkylation (HDA) style chemistry: toluene + H2 -> benzene + CH4 with a
fractional order in hydrogen, and benzene <-> biphenyl + H2 as a forward/reverse pair, all on
partial pressures, in an isothermal plug flow reactor at 30 to 60 bar. The kinetic constants
are illustrative, chosen for this test, not taken from literature. The reference integrates the
same rate laws with scipy, independently of Fluxsheet's reactor code; the flash downstream
checks that the hydrogen-rich gas at high pressure is recognised as gas. The last test puts
both into a complete, equally illustrative plant with gas and toluene recycles and three
distillation columns.
"""

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from fluxsheet import solve

R_GAS = 8.314462618
T_REF = 923.15
COMPONENTS = ["hydrogen", "methane", "benzene", "toluene", "biphenyl"]
# id: (stoichiometry, k0 at T_REF in mol/(m3 s Pa^n), Ea in J/mol, orders)
REACTIONS = {
    "HDA": ({"toluene": -1, "hydrogen": -1, "benzene": 1, "methane": 1}, 3e-9, 200e3,
            {"toluene": 1, "hydrogen": 0.5}),
    "DIPH": ({"benzene": -2, "biphenyl": 1, "hydrogen": 1}, 1e-12, 130e3, {"benzene": 2}),
    "DIPH_REV": ({"biphenyl": -1, "hydrogen": -1, "benzene": 2}, 5e-12, 100e3,
                 {"biphenyl": 1, "hydrogen": 1}),
}
FEED = {"toluene": 100.0, "hydrogen": 400.0, "methane": 20.0}  # kmol/h
DIAMETER = 1.2


def reactor_flowsheet(T_C, P_bar, length, spec_conversion=None):
    model = {
        "components": COMPONENTS,
        "property_method": "peng-robinson",
        "reactions": [{"id": rid, "stoichiometry": nu,
                       "kinetics": {"basis": "partial_pressure", "k0": k0, "Ea": f"{Ea} J/mol",
                                    "T_ref": f"{T_REF} K", "orders": orders}}
                      for rid, (nu, k0, Ea, orders) in REACTIONS.items()],
        "units": [
            {"id": "FEED", "type": "Feed", "params": {
                "T": "25 degC", "P": f"{P_bar} bar",
                "flows": {c: f"{f} kmol/h" for c, f in FEED.items()}}},
            {"id": "FURNACE", "type": "Heater", "params": {"T": f"{T_C} degC"}},
            {"id": "REACTOR", "type": "PFR", "params": {
                "reactions": list(REACTIONS), "length": f"{length} m",
                "diameter": f"{DIAMETER} m", "phase": "gas", "T": f"{T_C} degC",
                "key": "toluene"}},
            {"id": "COOLER", "type": "Heater", "params": {"T": "40 degC", "P": f"{P_bar - 2} bar"}},
            {"id": "FLASH", "type": "Flash"},
        ],
        "streams": [
            {"id": "S1", "from": "FEED.out", "to": "FURNACE.in"},
            {"id": "S2", "from": "FURNACE.out", "to": "REACTOR.in"},
            {"id": "S3", "from": "REACTOR.out", "to": "COOLER.in"},
            {"id": "S4", "from": "COOLER.out", "to": "FLASH.in"},
            {"id": "VAP", "from": "FLASH.vapor"},
            {"id": "LIQ", "from": "FLASH.liquid"},
        ],
    }
    if spec_conversion is not None:
        model["design_specs"] = [{
            "id": "CONVERSION", "target": {"unit": "REACTOR", "property": "conversion"},
            "value": spec_conversion,
            "vary": {"unit": "REACTOR", "param": "length", "lower": 0.5, "upper": 30}}]
    return model


def reference_outlet(T_C, P_bar, length):
    """Outlet flows in mol/s from integrating the rate laws along the tube (ideal gas)."""
    T, P, area = T_C + 273.15, P_bar * 1e5, math.pi * DIAMETER**2 / 4
    idx = {c: i for i, c in enumerate(COMPONENTS)}
    nus = [np.array([nu.get(c, 0.0) for c in COMPONENTS]) for nu, *_ in REACTIONS.values()]

    def rhs(_z, F):
        p = np.maximum(F, 0.0) / F.sum() * P
        rates = []
        for _nu, k0, Ea, orders in REACTIONS.values():
            r = k0 * math.exp(-Ea / R_GAS * (1 / T - 1 / T_REF))
            for c, a in orders.items():
                r *= p[idx[c]] ** a
            rates.append(r)
        return area * sum(r * nu for r, nu in zip(rates, nus))

    F0 = np.array([FEED.get(c, 0.0) / 3.6 for c in COMPONENTS])
    sol = solve_ivp(rhs, (0.0, length), F0, method="LSODA", rtol=1e-10, atol=1e-12)
    return dict(zip(COMPONENTS, sol.y[:, -1]))


@pytest.mark.parametrize("T_C, P_bar, length", [(650, 60, 9.0), (600, 30, 9.0), (700, 60, 3.0)])
def test_kinetic_reactor_matches_independent_integration(T_C, P_bar, length):
    r = solve(reactor_flowsheet(T_C, P_bar, length))
    assert r.ok, r.error
    expected = reference_outlet(T_C, P_bar, length)
    got = r.streams["S3"]["component_flows"]
    for c in COMPONENTS:
        assert got[c] == pytest.approx(expected[c], rel=1e-4, abs=1e-6), c
    tol0 = FEED["toluene"] / 3.6
    assert r.units["REACTOR"]["conversion"] == pytest.approx(
        1 - expected["toluene"] / tol0, abs=1e-5)


def test_reactor_length_for_a_conversion():
    """The design spec finds the tube length for 60 % conversion, like a hand iteration on
    the reference solution."""
    r = solve(reactor_flowsheet(650, 60, 5.0, spec_conversion=0.6))
    assert r.ok, r.error
    tol0 = FEED["toluene"] / 3.6
    L_ref = brentq(lambda L: 1 - reference_outlet(650, 60, L)["toluene"] / tol0 - 0.6, 0.5, 30)
    (spec,) = r.design_specs
    assert spec["converged"]
    assert spec["vary_value"] == pytest.approx(L_ref, rel=1e-3)
    assert r.units["REACTOR"]["conversion"] == pytest.approx(0.6, abs=1e-4)


def test_high_pressure_separator_splits_gas_from_aromatics():
    r = solve(reactor_flowsheet(650, 60, 9.0))
    assert r.ok, r.error
    assert r.streams["S1"]["vapor_fraction"] == pytest.approx(0.81, abs=0.02)  # 25 °C, 60 bar
    vap, liq = r.streams["VAP"]["component_flows"], r.streams["LIQ"]["component_flows"]
    assert vap["hydrogen"] > 0.98 * (vap["hydrogen"] + liq["hydrogen"])
    assert vap["methane"] > 0.9 * (vap["methane"] + liq["methane"])
    assert liq["toluene"] > 0.98 * (vap["toluene"] + liq["toluene"])
    assert liq["biphenyl"] > 0.999 * (vap["biphenyl"] + liq["biphenyl"])


# --- the whole plant ------------------------------------------------------------------------

def plant_flowsheet():
    """The reactor and separator above in a complete HDA plant: gas recycle with purge,
    stabilizer, benzene column, and a column that recycles toluene. Flowsheet and specs are
    illustrative, made up for this test; the distillate of the toluene column is set so
    high that some biphenyl recycles too, which the reverse reaction turns back to benzene."""
    model = reactor_flowsheet(650, 60, 9.0)
    units = {u["id"]: u for u in model["units"]}
    units["COOLER"]["params"]["P"] = "58 bar"
    feeds = [
        {"id": "TOLUENE", "type": "Feed", "params": {
            "T": "25 degC", "P": "60 bar", "flows": {"toluene": "100 kmol/h"}}},
        {"id": "HYDROGEN", "type": "Feed", "params": {
            "T": "25 degC", "P": "60 bar",
            "flows": {"hydrogen": "150 kmol/h", "methane": "8 kmol/h"}}},
        {"id": "MIX", "type": "Mixer"},
    ]

    def column(stages, feed_stage, P, condenser, R, D_to_F):
        return {"type": "DistillationColumn", "params": {
            "stages": stages, "feed_stage": feed_stage, "P": P, "condenser": condenser,
            "reflux_ratio": R, "distillate_to_feed": D_to_F}}

    rest = [
        {"id": "PURGE", "type": "Splitter", "params": {"split": 0.8}},
        {"id": "RECYCLE_COMP", "type": "Compressor", "params": {"P": "60 bar"}},
        {"id": "LETDOWN", "type": "Valve", "params": {"P": "3 bar"}},
        {"id": "STABILIZER", **column(8, 3, "3 bar", "partial", 1.0, 0.12)},
        {"id": "PRODUCT_COL", **column(22, 11, "1.2 bar", "total", 1.5, 0.62)},
        {"id": "HEAVY_COL", **column(10, 5, "1.2 bar", "total", 1.0, 0.95)},
        {"id": "TOL_PUMP", "type": "Pump", "params": {"P": "60 bar"}},
    ]
    model["units"] = feeds + [u for u in model["units"] if u["id"] != "FEED"] + rest
    model["streams"] = [
        {"id": "TOL-FRESH", "from": "TOLUENE.out", "to": "MIX.in"},
        {"id": "H2-FRESH", "from": "HYDROGEN.out", "to": "MIX.in"},
        {"id": "S1", "from": "MIX.out", "to": "FURNACE.in"},
        *[s for s in model["streams"] if s["id"] in ("S2", "S3", "S4")],
        {"id": "VAP", "from": "FLASH.vapor", "to": "PURGE.in"},
        {"id": "GAS-REC", "from": "PURGE.out1", "to": "RECYCLE_COMP.in"},
        {"id": "H2-REC", "from": "RECYCLE_COMP.out", "to": "MIX.in"},
        {"id": "PURGE-GAS", "from": "PURGE.out2"},
        {"id": "LIQ", "from": "FLASH.liquid", "to": "LETDOWN.in"},
        {"id": "STAB-IN", "from": "LETDOWN.out", "to": "STABILIZER.in"},
        {"id": "LIGHTS", "from": "STABILIZER.vapor_distillate"},
        {"id": "STAB-LIQ", "from": "STABILIZER.bottoms", "to": "PRODUCT_COL.in"},
        {"id": "BENZENE", "from": "PRODUCT_COL.distillate"},
        {"id": "PROD-BOT", "from": "PRODUCT_COL.bottoms", "to": "HEAVY_COL.in"},
        {"id": "TOL-REC-LP", "from": "HEAVY_COL.distillate", "to": "TOL_PUMP.in"},
        {"id": "TOL-REC", "from": "TOL_PUMP.out", "to": "MIX.in"},
        {"id": "BIPHENYL", "from": "HEAVY_COL.bottoms"},
    ]
    return model


def test_hda_plant_converges_from_scratch():
    # Early recycle passes feed the columns compositions far from the converged plant; the
    # toluene column used to fail there. Each column starts from its previous pass.
    r = solve(plant_flowsheet())
    assert r.ok, r.error
    assert r.converged
    s = r.streams
    # Aromatic rings: every toluene fed leaves as benzene, toluene or half a biphenyl.
    out = [s[k]["component_flows"] for k in ("BENZENE", "BIPHENYL", "PURGE-GAS", "LIGHTS")]
    rings = sum(f["benzene"] + f["toluene"] + 2 * f["biphenyl"] for f in out)
    assert rings == pytest.approx(100 / 3.6, rel=1e-5)
    assert s["BENZENE"]["mole_fractions"]["benzene"] > 0.95
    assert s["TOL-REC"]["component_flows"]["biphenyl"] > 0.1 / 3.6  # pushed overhead
    for column in ("STABILIZER", "PRODUCT_COL", "HEAVY_COL"):
        assert r.units[column]["iterations"] <= 3  # warm start from the previous pass
