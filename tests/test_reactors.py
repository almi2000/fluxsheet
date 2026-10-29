"""Reactions and reactors, checked against closed-form results and independent balances."""

import itertools
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

from fluxsheet import FlowsheetModel, solve
from fluxsheet.reactions import T0, R, parse_reactions
from fluxsheet.thermo_package import get_package

EXAMPLES = Path(__file__).parent.parent / "examples"

SHIFT = {"id": "WGS", "stoichiometry": {"carbon monoxide": -1, "water": -1,
                                        "carbon dioxide": 1, "hydrogen": 1}}
REFORMING = {"id": "SMR", "stoichiometry": {"methane": -1, "water": -1,
                                            "carbon monoxide": 1, "hydrogen": 3}}
SYNGAS = ["methane", "water", "carbon monoxide", "carbon dioxide", "hydrogen"]


def run(components, feed, reactor, reactions=(), T="400 K", P="1 bar", method="ideal"):
    """Feed -> one reactor -> product; returns (result, inlet stream, outlet stream)."""
    r = solve({
        "components": components,
        "property_method": method,
        "reactions": list(reactions),
        "units": [
            {"id": "F", "type": "Feed", "params": {"T": T, "P": P, "flows": feed}},
            {"id": "R", **reactor},
        ],
        "streams": [{"id": "IN", "from": "F.out", "to": "R.in"}, {"id": "OUT", "from": "R.out"}],
    })
    assert r.ok, r.error
    return r, r.streams["IN"], r.streams["OUT"]


def element_flows(components, flows):
    pkg = get_package(tuple(components), "ideal")
    out: dict[str, float] = {}
    for name, atoms in zip(pkg.names, pkg.consts.atomss):
        for el, n in atoms.items():
            out[el] = out.get(el, 0.0) + n * flows[name]
    return out


# --- first-order kinetics: closed-form conversions --------------------------------------
# Gas-phase isomerization A -> B (equal moles, ideal gas, isothermal, no pressure drop) has
# a constant volumetric flow, so tau = V/Q exactly and r = k c_A gives the textbook results.

ISOMER = ["n-butane", "isobutane"]
K1 = 0.05  # 1/s
FIRST_ORDER = {"id": "ISO", "stoichiometry": {"n-butane": -1, "isobutane": 1},
               "kinetics": {"k0": K1, "Ea": "0 kJ/mol"}}


def _first_order(reactor_type, tau):
    T, P, F = 400.0, 1e5, 1.0
    V = tau * F * R * T / P
    r, _, out = run(ISOMER, {"n-butane": F}, {
        "type": reactor_type,
        "params": {"reactions": ["ISO"], "V": V, "phase": "gas", "T": T},
    }, [FIRST_ORDER], T=T, P=P)
    unit = r.units["R"]
    assert unit["residence_time"] == pytest.approx(tau, rel=1e-9)
    X = 1.0 - out["component_flows"]["n-butane"] / F
    assert unit["conversion"] == pytest.approx(X, rel=1e-12)
    return X, unit


@pytest.mark.parametrize("tau", [2.0, 20.0, 100.0])
def test_first_order_cstr(tau):
    X, _ = _first_order("CSTR", tau)
    k_tau = K1 * tau
    assert X == pytest.approx(k_tau / (1.0 + k_tau), rel=1e-7)


@pytest.mark.parametrize("tau", [2.0, 20.0, 100.0])
def test_first_order_pfr(tau):
    X, unit = _first_order("PFR", tau)
    assert X == pytest.approx(1.0 - math.exp(-K1 * tau), rel=1e-6)
    conditions = unit["profiles"]["conditions"]
    conv = next(s for s in conditions["series"] if s["label"].startswith("Conversion"))
    V = np.array(conditions["x"])
    expected = 1.0 - np.exp(-K1 * tau * V / V[-1])
    assert np.allclose(conv["values"], expected, rtol=1e-5, atol=1e-9)


def test_pfr_length_and_diameter_give_the_volume():
    T, P, F, tau = 400.0, 1e5, 1.0, 20.0
    V = tau * F * R * T / P
    D = 0.5
    L = V / (math.pi / 4 * D**2)
    r, _, out = run(ISOMER, {"n-butane": F}, {
        "type": "PFR",
        "params": {"reactions": ["ISO"], "length": L, "diameter": D, "phase": "gas", "T": T},
    }, [FIRST_ORDER], T=T, P=P)
    assert 1.0 - out["component_flows"]["n-butane"] == pytest.approx(1 - math.exp(-1.0), 1e-6)
    assert r.units["R"]["profiles"]["composition"]["x"][-1] == pytest.approx(L)


# --- equilibrium constants --------------------------------------------------------------

def test_water_gas_shift_equilibrium_constant():
    pkg = get_package(tuple(SYNGAS), "ideal")
    (wgs,) = parse_reactions([FlowsheetModel.model_validate(
        {"components": SYNGAS, "units": [], "reactions": [SHIFT]}).reactions[0]], pkg
    ).by_id.values()
    # NIST-JANAF: ΔfG°(1000 K) of CO2, CO, H2O(g) = -395.886, -200.275, -192.590 kJ/mol,
    # so ln K = 3.021 kJ/mol / (R * 1000 K) and K = 1.44.
    assert math.exp(wgs.ln_keq(pkg, 1000.0)) == pytest.approx(1.44, rel=0.03)
    # Exothermic by about 41 kJ/mol, so K falls with temperature (van 't Hoff).
    assert wgs.heat_of_reaction(pkg) == pytest.approx(-41.2e3, rel=0.01)
    assert wgs.ln_keq(pkg, 800.0) > wgs.ln_keq(pkg, 1000.0) > wgs.ln_keq(pkg, 1200.0)


def test_user_keq_correlation_is_used():
    pkg = get_package(tuple(SYNGAS), "ideal")
    model = FlowsheetModel.model_validate({"components": SYNGAS, "units": [], "reactions": [
        {**SHIFT, "keq": {"A": -4.0, "B": 4600.0}}]})
    (wgs,) = parse_reactions(model.reactions, pkg).by_id.values()
    assert wgs.ln_keq(pkg, 920.0) == pytest.approx(-4.0 + 4600.0 / 920.0)


# --- equilibrium and Gibbs reactors -----------------------------------------------------

STEAM_REFORMER = {"methane": 1.0, "water": 3.0}


@pytest.mark.parametrize("T, P", [("1000 K", "1 bar"), ("1100 K", "20 bar")])
def test_gibbs_and_equilibrium_reactor_agree_on_steam_reforming(T, P):
    _, _, gibbs = run(SYNGAS, STEAM_REFORMER, {"type": "GibbsReactor", "params": {"T": T}},
                      T=T, P=P)
    r, _, eq = run(SYNGAS, STEAM_REFORMER, {
        "type": "EquilibriumReactor", "params": {"reactions": ["SMR", "WGS"], "T": T},
    }, [REFORMING, SHIFT], T=T, P=P)
    for c in SYNGAS:
        assert eq["component_flows"][c] == pytest.approx(gibbs["component_flows"][c],
                                                         rel=1e-6, abs=1e-12)
    # Both reactions are at equilibrium: Q = K for each, from the outlet composition.
    pkg = get_package(tuple(SYNGAS), "ideal")
    reactions = parse_reactions(FlowsheetModel.model_validate(
        {"components": SYNGAS, "units": [], "reactions": [REFORMING, SHIFT]}).reactions, pkg)
    T_K, P_bar = eq["T"], eq["P"] / 1e5
    y = eq["mole_fractions"]
    for rx in reactions.by_id.values():
        ln_q = sum(v * math.log(y[c] * P_bar) for c, v in zip(SYNGAS, rx.nu) if v)
        assert ln_q == pytest.approx(rx.ln_keq(pkg, T_K), abs=1e-6)
    assert r.units["R"]["conversion"] > 0.5  # methane is mostly reformed at these conditions


def test_gibbs_reactor_pressure_shifts_reforming_back():
    def methane_left(P):
        _, _, out = run(SYNGAS, STEAM_REFORMER, {"type": "GibbsReactor", "params": {"T": "1000 K"}},
                        T="1000 K", P=P)
        return out["component_flows"]["methane"]

    # Le Chatelier: reforming makes moles, so pressure leaves more methane unconverted.
    assert methane_left("30 bar") > methane_left("1 bar")


def test_gibbs_reactor_conserves_elements_and_keeps_inerts():
    components = [*SYNGAS, "nitrogen"]
    feed = {**STEAM_REFORMER, "nitrogen": 0.5}
    _, inlet, out = run(components, feed, {
        "type": "GibbsReactor", "params": {"T": "900 K", "products": SYNGAS},
    }, T="900 K")
    a, b = element_flows(components, inlet["component_flows"]), element_flows(
        components, out["component_flows"])
    for el in a:
        assert b[el] == pytest.approx(a[el], rel=1e-9)
    assert out["component_flows"]["nitrogen"] == pytest.approx(0.5)


def test_adiabatic_equilibrium_reactor_closes_energy_balance():
    r, inlet, out = run(SYNGAS, {"carbon monoxide": 1.0, "water": 2.0, "hydrogen": 0.2},
                        {"type": "EquilibriumReactor", "params": {"reactions": ["WGS"]}},
                        [SHIFT], T="600 K")
    unit = r.units["R"]
    assert abs(unit["duty"]) < 1e-3  # W, with ~100 kW of reaction heat
    assert out["T"] > inlet["T"] + 50  # exothermic shift heats the gas
    pkg = get_package(tuple(SYNGAS), "ideal")
    (wgs,) = parse_reactions(FlowsheetModel.model_validate(
        {"components": SYNGAS, "units": [], "reactions": [SHIFT]}).reactions, pkg).by_id.values()
    y = out["mole_fractions"]
    ln_q = math.log(y["carbon dioxide"] * y["hydrogen"] / (y["carbon monoxide"] * y["water"]))
    assert ln_q == pytest.approx(wgs.ln_keq(pkg, out["T"]), abs=1e-5)


# --- energy balance: adiabatic temperature rise -----------------------------------------

def _independent_adiabatic_T(components, feed, products, T_in):
    """Solve sum_out F_i H_i(T) = sum_in F_i H_i(T_in), H_i = ΔHf_i + ∫Cp dT (ideal gas)."""
    pkg = get_package(tuple(components), "ideal")

    def H(flows, T):
        total = 0.0
        for i, name in enumerate(components):
            cp = pkg.props.HeatCapacityGases[i]
            total += flows.get(name, 0.0) * (pkg.consts.Hfgs[i]
                                             + cp.T_dependent_property_integral(T0, T))
        return total

    h_in = H(feed, T_in)
    return brentq(lambda T: H(products, T) - h_in, 300.0, 3000.0, xtol=1e-9)


def test_adiabatic_temperature_rise_matches_enthalpy_balance():
    feed = {"carbon monoxide": 1.0, "water": 2.0}
    r, _, out = run(SYNGAS, feed, {
        "type": "StoichiometricReactor",
        "params": {"reactions": ["WGS"], "conversions": ["0.6"]},
    }, [SHIFT], T="550 K")
    products = {"carbon monoxide": 0.4, "water": 1.4, "carbon dioxide": 0.6, "hydrogen": 0.6}
    expected = _independent_adiabatic_T(SYNGAS, feed, products, 550.0)
    assert out["T"] == pytest.approx(expected, abs=0.01)
    assert out["T"] - 550.0 > 50.0
    assert r.units["R"]["conversion"] == pytest.approx(0.6)


def test_adiabatic_pfr_follows_the_adiabatic_line():
    """The ODE energy balance and a direct enthalpy balance agree on outlet temperature."""
    shift = {**SHIFT, "kinetics": {"k0": 1000.0, "Ea": "60 kJ/mol",
                                   "orders": {"carbon monoxide": 1, "water": 1},
                                   "equilibrium_limited": True}}
    feed = {"carbon monoxide": 1.0, "water": 2.0}
    r, _, out = run(SYNGAS, feed, {
        "type": "PFR", "params": {"reactions": ["WGS"], "V": "0.5 m3", "phase": "gas"},
    }, [shift], T="600 K")
    xi = 1.0 - out["component_flows"]["carbon monoxide"]
    assert 0.1 < xi < 0.99
    products = {"carbon monoxide": 1 - xi, "water": 2 - xi, "carbon dioxide": xi,
                "hydrogen": xi}
    expected = _independent_adiabatic_T(SYNGAS, feed, products, 600.0)
    assert out["T"] == pytest.approx(expected, abs=0.01)
    T_profile = next(s for s in r.units["R"]["profiles"]["conditions"]["series"]
                     if s["label"] == "Temperature")["values"]
    assert T_profile[0] == pytest.approx(600.0)
    assert T_profile[-1] == pytest.approx(expected, abs=0.05)
    assert all(b >= a - 1e-6 for a, b in itertools.pairwise(T_profile))  # heats monotonically
    assert T_profile[len(T_profile) // 2] < T_profile[-1] - 1.0  # still reacting midway


def test_adiabatic_cstr_matches_its_energy_balance():
    shift = {**SHIFT, "kinetics": {"k0": 2e5, "Ea": "60 kJ/mol",
                                   "orders": {"carbon monoxide": 1, "water": 1},
                                   "equilibrium_limited": True}}
    feed = {"carbon monoxide": 1.0, "water": 2.0}
    _, _, out = run(SYNGAS, feed, {
        "type": "CSTR", "params": {"reactions": ["WGS"], "V": "0.5 m3", "phase": "gas"},
    }, [shift], T="600 K")
    xi = 1.0 - out["component_flows"]["carbon monoxide"]
    products = {"carbon monoxide": 1 - xi, "water": 2 - xi, "carbon dioxide": xi,
                "hydrogen": xi}
    assert out["T"] == pytest.approx(_independent_adiabatic_T(SYNGAS, feed, products, 600.0),
                                     abs=0.01)
    # The rate at outlet conditions times the volume equals the extent.
    pkg = get_package(tuple(SYNGAS), "ideal")
    (wgs,) = parse_reactions(FlowsheetModel.model_validate(
        {"components": SYNGAS, "units": [], "reactions": [shift]}).reactions, pkg).by_id.values()
    z = np.array([out["mole_fractions"][c] for c in SYNGAS])
    Vm = R * out["T"] / out["P"]
    assert 0.5 * wgs.rate(pkg, out["T"], out["P"], z, Vm, "gas") == pytest.approx(xi, rel=1e-6)


# --- stoichiometric and yield reactors --------------------------------------------------

@pytest.mark.parametrize("mode", ["sequential", "parallel"])
def test_stoichiometric_reactor_conserves_mass_and_elements(mode):
    combustion = {"id": "COMB", "stoichiometry": {"methane": -1, "oxygen": -2,
                                                  "carbon dioxide": 1, "water": 2}}
    partial = {"id": "POX", "stoichiometry": {"methane": -2, "oxygen": -1,
                                              "carbon monoxide": 2, "hydrogen": 4}}
    components = [*SYNGAS, "oxygen"]
    feed = {"methane": 1.0, "oxygen": 1.2}
    r, inlet, out = run(components, feed, {
        "type": "StoichiometricReactor",
        "params": {"reactions": ["COMB", "POX"], "conversions": ["40%", "0.5"], "mode": mode,
                   "T": "900 K"},
    }, [combustion, partial], T="500 K")
    assert out["mass_flow"] == pytest.approx(inlet["mass_flow"], rel=1e-12)
    a = element_flows(components, inlet["component_flows"])
    b = element_flows(components, out["component_flows"])
    for el in a:
        assert b[el] == pytest.approx(a[el], rel=1e-12)
    methane_after_comb = 0.6
    pox_basis = methane_after_comb if mode == "sequential" else 1.0
    xi = r.units["R"]["extents"]
    assert xi["COMB"] == pytest.approx(0.4)
    assert xi["POX"] == pytest.approx(0.5 * pox_basis / 2)
    assert out["T"] == pytest.approx(900.0)
    assert r.units["R"]["duty"] < 0  # combustion heat must be removed to hold 900 K


def test_yield_reactor_conserves_mass():
    components = ["n-hexane", "methane", "ethylene", "propylene", "hydrogen"]
    r, inlet, out = run(components, {"n-hexane": 1.0, "hydrogen": 0.1}, {
        "type": "YieldReactor",
        "params": {"yields": {"methane": 0.2, "ethylene": 0.5, "propylene": 0.3},
                   "inerts": ["hydrogen"], "T": "1100 K"},
    }, T="700 K")
    assert out["mass_flow"] == pytest.approx(inlet["mass_flow"], rel=1e-12)
    assert out["component_flows"]["n-hexane"] == pytest.approx(0.0, abs=1e-15)
    assert out["component_flows"]["hydrogen"] == pytest.approx(0.1)
    pkg = get_package(tuple(components), "ideal")
    mw = dict(zip(pkg.names, pkg.mws))
    ethylene_mass = out["component_flows"]["ethylene"] * mw["ethylene"]
    assert ethylene_mass == pytest.approx(0.5 * mw["n-hexane"])
    assert r.units["R"]["duty"] > 0  # cracking is endothermic


# --- liquid-phase kinetics and the examples ---------------------------------------------

def test_liquid_esterification_approaches_user_equilibrium():
    components = ["acetic acid", "ethanol", "ethyl acetate", "water"]
    rxn = {"id": "EST",
           "stoichiometry": {"acetic acid": -1, "ethanol": -1, "ethyl acetate": 1, "water": 1},
           "kinetics": {"k0": 1e-4, "Ea": "0 kJ/mol", "orders": {"acetic acid": 1, "ethanol": 1},
                        "equilibrium_limited": True},
           "keq": {"A": math.log(4.0)}}
    feed = {"acetic acid": 1.0, "ethanol": 1.0}
    _, _, out = run(components, feed, {
        "type": "PFR", "params": {"reactions": ["EST"], "V": "50 m3", "phase": "liquid",
                                  "T": "60 degC"},
    }, [rxn], T="60 degC", method="ideal")
    x = out["mole_fractions"]
    Q = x["ethyl acetate"] * x["water"] / (x["acetic acid"] * x["ethanol"])
    assert Q == pytest.approx(4.0, rel=1e-3)  # K = 4: two thirds conversion of 1:1 feed
    assert 1.0 - out["component_flows"]["acetic acid"] == pytest.approx(2.0 / 3.0, rel=1e-3)


@pytest.mark.parametrize("name", ["steam_reforming", "ethyl_acetate"])
def test_reactor_examples_solve(name):
    r = solve(FlowsheetModel.load(EXAMPLES / f"{name}.json"))
    assert r.ok, r.error
    assert not [w for w in r.warnings if "balance" in w]
    assert any(u.get("profiles") for u in r.units.values())


# --- errors -----------------------------------------------------------------------------

def _error(reactions, reactor, components=SYNGAS):
    r = solve({
        "components": components, "property_method": "ideal", "reactions": reactions,
        "units": [{"id": "F", "type": "Feed",
                   "params": {"T": "600 K", "P": "1 bar", "flows": {components[0]: 1.0}}},
                  {"id": "R", **reactor}],
        "streams": [{"id": "IN", "from": "F.out", "to": "R.in"}, {"id": "OUT", "from": "R.out"}],
    })
    assert not r.ok
    return r.error


def test_unbalanced_reaction_is_explained():
    bad = {"id": "BAD", "stoichiometry": {"methane": -1, "water": -1, "carbon monoxide": 1,
                                          "hydrogen": 2}}
    err = _error([bad], {"type": "GibbsReactor"})
    assert err["where"] == "reaction BAD"
    assert "H: 6 in reactants, 4 in products" in err["message"]
    assert "coefficients" in err["hint"]


def test_unknown_component_in_reaction():
    bad = {"id": "R1", "stoichiometry": {"methane": -1, "unobtainium": 1}}
    err = _error([bad], {"type": "GibbsReactor"})
    assert "'unobtainium' is not in the flowsheet's component list" in err["message"]
    assert "methane" in err["hint"]


def test_cstr_without_kinetics_is_explained():
    err = _error([REFORMING], {"type": "CSTR", "params": {"reactions": ["SMR"], "V": 1.0}})
    assert err["where"] == "R"
    assert "SMR has no kinetics" in err["message"]
    assert "rate law" in err["hint"]


def test_unknown_reaction_id_lists_the_available_ones():
    err = _error([REFORMING], {"type": "EquilibriumReactor", "params": {"reactions": ["SMR2"]}})
    assert "unknown reaction 'SMR2'" in err["message"]
    assert "SMR" in err["hint"]


def test_missing_formation_enthalpy_names_the_component():
    # The database has no ΔHf° for HMF, so no reactor energy balance can include it.
    err = _error([], {"type": "YieldReactor",
                      "params": {"yields": {"5-hydroxymethylfurfural": 1.0}, "T": "400 K"}},
                 components=["water", "5-hydroxymethylfurfural"])
    assert "no formation enthalpy in the database for 5-hydroxymethylfurfural" in err["message"]


def test_missing_gibbs_energy_names_the_component():
    # Limonene has ΔHf° but no ΔGf°, so its equilibrium is unknown.
    err = _error([], {"type": "GibbsReactor", "params": {"T": "600 K"}},
                 components=["methane", "hydrogen", "limonene"])
    assert "no ideal-gas formation data" in err["message"]
    assert "limonene" in err["message"]
    assert "leave these compounds out" in err["hint"]


def test_reaction_listed_twice_is_explained():
    err = _error([SHIFT], {"type": "StoichiometricReactor",
                           "params": {"reactions": ["WGS", "WGS"], "conversions": ["1", "1"]}},
                 components=["carbon monoxide", "water", "carbon dioxide", "hydrogen"])
    assert "listed twice" in err["message"]
