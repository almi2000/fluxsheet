"""Henry components in activity-coefficient methods, and user-defined (estimated) compounds.

Reference values: gas solubilities in water at 25 °C and 1 atm partial pressure (mole
fractions) from Wilhelm, Battino and Wilcock, Chem. Rev. 77 (1977) 219 and, for CO2, Carroll,
Slupsky and Mather, J. Phys. Chem. Ref. Data 20 (1991) 1201; heat of solution of O2 in water
from the same review (about -12 kJ/mol). Database values for ethanol from the chemicals
library. Joback's method has typical errors of a few percent for Tb and Tc and up to about
10 % for Pc (Joback and Reid 1987; Poling et al., The Properties of Gases and Liquids).
"""

import math

import pytest
from fastapi.testclient import TestClient

from fluxsheet import solve
from fluxsheet.custom_compounds import CustomCompound, estimate
from fluxsheet.errors import SimulationError
from fluxsheet.henry import (
    IAPWS_G704,
    WATER,
    HenryDilutedModel,
    henry_constant,
    iapws_henry,
)
from fluxsheet.server import app
from fluxsheet.thermo_package import PropertyPackage, get_package, suggest_method

ATM = 101325.0
T25 = 298.15
ACTIVITY = ["nrtl", "uniquac", "unifac", "wilson"]

# Mole fraction of the gas in water at 25 °C and 1 atm partial pressure of the gas.
SOLUBILITY_25C = {
    "oxygen": 2.29e-5,
    "nitrogen": 1.18e-5,
    "carbon dioxide": 6.1e-4,
    "hydrogen": 1.41e-5,
    "methane": 2.48e-5,
    "argon": 2.52e-5,
}

ETHANOL_GROUPS = {"joback_groups": {"-CH3": 1, "-CH2-": 1, "-OH (alcohol)": 1},
                  "unifac_groups": {"CH3": 1, "CH2": 1, "OH(P)": 1}}


def custom(**kw) -> CustomCompound:
    return CustomCompound.from_dict({"name": "my-ethanol", "formula": "C2H6O", **kw})


# --- Henry's law: data and solubilities -------------------------------------------------

@pytest.mark.parametrize("gas", list(SOLUBILITY_25C))
def test_gas_solubility_in_water_matches_literature(gas):
    """Flash gas over water at 25 °C; the liquid's gas content per atm of partial pressure."""
    pkg = get_package(("water", gas), "nrtl")
    assert pkg.henry_components == [gas]
    st = pkg.flash([0.5, 0.5], T=T25, P=ATM + 3170.0)  # gas partial pressure about 1 atm
    x_per_atm = st.x[1] / (st.y[1] * st.P / ATM)
    assert x_per_atm == pytest.approx(SOLUBILITY_25C[gas], rel=0.05)


def test_iapws_fit_and_temperature_dependence():
    """The ln H form thermo evaluates reproduces IAPWS G7-04 within 1 % up to 250 °C, and
    H of O2 passes through its maximum near 100 °C, as measured."""
    for cas in IAPWS_G704:
        k = henry_constant(cas, WATER)
        for T in (275.0, 300.0, 350.0, 400.0, 450.0, 520.0):
            assert k.value(T) == pytest.approx(iapws_henry(cas, T), rel=0.01)
    o2 = henry_constant("7782-44-7", WATER)
    assert o2.value(353.15) > o2.value(T25)
    assert o2.value(473.15) < o2.value(373.15)


def test_heat_of_solution_of_oxygen():
    """The liquid enthalpy of a Henry solute is the ideal-gas enthalpy plus the heat of
    solution, -R T^2 dlnH/dT: dissolving O2 in water releases about 12 kJ/mol."""
    pkg = get_package(("water", "oxygen"), "nrtl")
    x = [1.0 - 1e-5, 1e-5]
    h_liq = pkg.phase_props("liquid", T25, ATM, x).H
    h_water = pkg.phase_props("liquid", T25, ATM, [1.0, 0.0]).H
    h_gas = pkg.phase_props("vapor", T25, ATM, [0.0, 1.0]).H
    dissolution = (h_liq - x[0] * h_water) / x[1] - h_gas
    assert dissolution == pytest.approx(-12e3, abs=2e3)


@pytest.mark.parametrize("method", ACTIVITY)
def test_air_water_flash_with_activity_methods(method):
    """Air over water: dissolved N2 and O2 at their Henry solubility, no liquid nitrogen,
    and a PH flash returns the temperature (enthalpies consistent)."""
    pkg = get_package(("water", "nitrogen", "oxygen"), method)
    assert pkg.henry_components == ["nitrogen", "oxygen"]
    st = pkg.flash([0.5, 0.39, 0.11], T=T25, P=ATM)
    assert 0.4 < st.VF < 0.6 and not st.liquid_split
    x_water, x_n2, x_o2 = st.x
    assert x_water > 0.9999
    assert x_n2 == pytest.approx(SOLUBILITY_25C["nitrogen"] * st.y[1] * st.P / ATM, rel=0.05)
    assert x_o2 == pytest.approx(SOLUBILITY_25C["oxygen"] * st.y[2] * st.P / ATM, rel=0.05)
    assert pkg.flash([0.5, 0.39, 0.11], P=ATM, H=st.H).T == pytest.approx(T25, abs=0.01)


def test_pure_gas_and_gas_only_mixtures():
    """A stream of only the dissolved gas is plain vapor (no crash in the mixing rule), and a
    gas-only flowsheet has no Henry components, as before."""
    pkg = get_package(("water", "nitrogen"), "nrtl")
    assert pkg.flash([0.0, 1.0], T=T25, P=ATM).phase == "vapor"
    assert pkg.phase_props("liquid", T25, ATM, [0.0, 1.0]).H < 0.0
    gases = PropertyPackage(["nitrogen", "oxygen"], "nrtl")
    assert gases.henry_components == []


def test_bubble_point_of_gas_saturated_water():
    """Water saturated with air at 25 °C and 1 atm is at its bubble point there: the first
    bubble is mostly air. With a tenth of that gas it boils between 25 and 100 °C."""
    pkg = get_package(("water", "nitrogen", "oxygen"), "nrtl")
    x_n2, x_o2 = 0.78 * SOLUBILITY_25C["nitrogen"], 0.21 * SOLUBILITY_25C["oxygen"]
    sat = pkg.saturation_point([1 - x_n2 - x_o2, x_n2, x_o2], "bubble", P=ATM)
    assert sat.T == pytest.approx(T25, abs=3.0)
    assert sat.y[0] < 0.05
    lean = pkg.saturation_point([1 - (x_n2 + x_o2) / 10, x_n2 / 10, x_o2 / 10], "bubble", P=ATM)
    assert 330.0 < lean.T < 373.0


def test_mixed_solvents_use_their_own_constants():
    """N2 is far more soluble in toluene than in water; with both liquids present each
    layer gets its own solubility, from the constant of its main solvent."""
    pkg = get_package(("water", "toluene", "nitrogen"), "nrtl")
    pairs = {s["pair"]: s["source"] for s in pkg.parameter_sources}
    assert "IAPWS" in pairs["nitrogen in water"]
    assert "Peng-Robinson" in pairs["nitrogen in toluene"]
    st = pkg.flash([0.45, 0.45, 0.1], T=T25, P=ATM)
    assert st.liquid_split
    assert st.x_light[2] > 50 * st.x_heavy[2]  # toluene-rich layer holds most N2


def test_explicit_henry_components_and_errors():
    """``henry_components`` overrides the automatic choice: H2S becomes one, [] turns it off;
    unknown names and a list without any solvent are refused."""
    pkg = PropertyPackage(["water", "hydrogen sulfide"], "nrtl",
                          henry_components=("hydrogen sulfide",))
    assert pkg.henry_components == ["hydrogen sulfide"]
    off = PropertyPackage(["water", "nitrogen"], "nrtl", henry_components=())
    assert off.henry_components == []
    eos = PropertyPackage(["water", "nitrogen"], "srk", henry_components=("nitrogen",))
    assert any("henry_components is ignored" in w for w in eos.warnings)
    with pytest.raises(SimulationError, match="not among the components"):
        PropertyPackage(["water", "nitrogen"], "nrtl", henry_components=("argon",))
    with pytest.raises(SimulationError, match="every component"):
        PropertyPackage(["water", "nitrogen"], "nrtl", henry_components=("water", "nitrogen"))


def test_missing_henry_constant_is_warned():
    """Without a constant for a solvent, the other solvents' constants are used, with a
    warning (an extrapolated vapor pressure would be orders of magnitude off)."""
    pkg = PropertyPackage(["water", "my-ethanol", "oxygen"], "nrtl",
                          custom=(custom(**ETHANOL_GROUPS),))
    pairs = {s["pair"]: s["source"] for s in pkg.parameter_sources}
    assert "none" in pairs["oxygen in my-ethanol"]
    x = pkg.phase_props("liquid", T25, ATM, [0.9, 0.0999, 1e-4])  # still finite
    assert all(math.isfinite(v) for v in x.lnphi)
    assert any("no Henry constant for oxygen in my-ethanol" in w for w in pkg.warnings)


def test_eos_methods_ignore_henry():
    pkg = get_package(("water", "nitrogen"), "peng-robinson")
    assert pkg.henry_components == []


def test_henry_diluted_model_derivatives():
    """The wrapper's G^E derivatives agree with finite differences; solutes have gamma = 1
    and the solvents those of their solute-free mixture."""
    from thermo import NRTL

    inner = NRTL(T=300.0, xs=[0.5, 0.5], tau_as=[[0, 0.3], [-0.2, 0]],
                 tau_bs=[[0, 400.0], [200.0, 0]], alpha_cs=[[0, 0.3], [0.3, 0]])
    xs = [0.4, 0.1, 0.45, 0.05]
    m = HenryDilutedModel(inner, [0, 2], 4, T=330.0, xs=xs)
    gammas = m.gammas()
    ref = inner.to_T_xs(330.0, [0.4 / 0.85, 0.45 / 0.85]).gammas()
    assert gammas[1] == gammas[3] == 1.0
    assert gammas[0] == pytest.approx(ref[0]) and gammas[2] == pytest.approx(ref[1])
    h = 1e-6
    for i in range(4):
        up, down = list(xs), list(xs)
        up[i] += h
        down[i] -= h
        fd = (m.to_T_xs(330.0, up).GE() - m.to_T_xs(330.0, down).GE()) / (2 * h)
        assert m.dGE_dxs()[i] == pytest.approx(fd, rel=1e-5, abs=1e-6)
        for j in range(4):
            fd2 = (m.to_T_xs(330.0, up).dGE_dxs()[j] - m.to_T_xs(330.0, down).dGE_dxs()[j]) / (2 * h)
            assert m.d2GE_dxixjs()[i][j] == pytest.approx(fd2, rel=1e-4, abs=1e-3)
    dT = (m.to_T_xs(330.001, xs).GE() - m.to_T_xs(329.999, xs).GE()) / 0.002
    assert m.dGE_dT() == pytest.approx(dT, rel=1e-6)
    assert m.HE() == pytest.approx(m.GE() - 330.0 * m.dGE_dT())


def gas_water_flowsheet(P="1 atm", **extra):
    return {
        "components": ["water", "nitrogen"], "property_method": "nrtl", **extra,
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "T": "25 degC", "P": P, "flows": {"water": "90 kmol/h", "nitrogen": "10 kmol/h"}}},
            {"id": "D", "type": "Flash", "params": {"T": "25 degC", "P": P}},
        ],
        "streams": [{"id": "S1", "from": "F.out", "to": "D.in"},
                    {"id": "V", "from": "D.vapor"}, {"id": "L", "from": "D.liquid"}],
    }


def test_flowsheet_reports_henry_components_and_pressure_warning():
    r = solve(gas_water_flowsheet())
    assert r.ok, r.error
    assert any("nitrogen is a Henry component" in w for w in r.warnings)
    assert any(s["pair"] == "nitrogen in water" for s in r.parameter_sources)
    assert r.streams["L"]["mole_fractions"]["nitrogen"] < 1e-4
    high = solve(gas_water_flowsheet(P="60 bar"))
    assert high.ok, high.error
    assert any("no pressure correction" in w for w in high.warnings)
    off = solve(gas_water_flowsheet(henry_components=[]))
    assert off.ok and not any("Henry component" in w for w in off.warnings)


# --- user-defined compounds ---------------------------------------------------------------

def test_joback_estimates_for_ethanol_within_typical_errors():
    db = get_package(("ethanol",), "peng-robinson").consts
    est = estimate(custom(**ETHANOL_GROUPS))
    v = est.values
    assert v["MW"] == pytest.approx(db.MWs[0], rel=1e-3)
    assert v["Tb"] == pytest.approx(db.Tbs[0], rel=0.05)
    assert v["Tc"] == pytest.approx(db.Tcs[0], rel=0.05)
    assert v["Pc"] == pytest.approx(db.Pcs[0], rel=0.12)
    assert v["Hf"] == pytest.approx(db.Hfgs[0], abs=5e3)
    assert v["Gf"] == pytest.approx(db.Gfgs[0], abs=5e3)
    assert v["r"] == pytest.approx(db.UNIFAC_Rs[0], rel=1e-3)
    assert {"Tb", "Tc", "Pc", "omega", "Hf"} <= set(est.estimated)
    pkg = PropertyPackage(["my-ethanol"], "ideal", custom=(custom(**ETHANOL_GROUPS),))
    ref = get_package(("ethanol",), "ideal")
    for T in (300.0, 400.0, 600.0):
        assert pkg.props.HeatCapacityGases[0](T) == pytest.approx(
            ref.props.HeatCapacityGases[0](T), rel=0.03)


def test_joback_groups_from_smiles_with_rdkit():
    pytest.importorskip("rdkit")
    est = estimate(custom(smiles="CCO"))
    assert est.joback == {1: 1, 2: 1, 20: 1}  # -CH3, -CH2-, -OH (alcohol)
    assert est.estimated["Joback groups"].startswith("RDKit")


def test_joback_group_without_a_contribution():
    """Joback has no Gf or Cp value for -N= (nonring): Gf stays unknown, Cp falls back to
    Lastovka-Shaw; =NH has no Tc, which must then be given."""
    est = estimate(CustomCompound.from_dict({
        "name": "imine", "formula": "C3H7N", "joback_groups": {"-CH3": 2, "-N= (nonring)": 1,
                                                                 "=CH-": 1}}))
    assert est.values["Gf"] is None and est.cp_coeffs is None
    with pytest.raises(SimulationError) as err:
        estimate(CustomCompound.from_dict({
            "name": "x", "formula": "C2H5N", "joback_groups": {"-CH3": 1, "=CH-": 1, "=NH": 1}}))
    assert "=NH" in (err.value.hint or "")


def test_flowsheet_file_with_antoine_constants():
    r = solve(custom_flowsheet(antoine={"A": 10.2370, "B": 1592.864, "C": -46.976},
                               constants={"Tc": 514.7, "Pc": "62.68 bar"}, **ETHANOL_GROUPS))
    assert r.ok, r.error


def test_given_constants_are_used_and_not_flagged():
    """With the database's Tb, Tc, Pc and omega given, the Lee-Kesler vapor pressure is
    close to ethanol's own correlation, and only what was missing is reported."""
    db = get_package(("ethanol",), "peng-robinson").consts
    c = custom(constants={"Tb": db.Tbs[0], "Tc": f"{db.Tcs[0] - 273.15} degC",
                          "Pc": f"{db.Pcs[0] / 1e5} bar", "omega": db.omegas[0]},
               **ETHANOL_GROUPS)
    est = estimate(c)
    assert est.values["Tc"] == pytest.approx(db.Tcs[0])
    assert est.values["Pc"] == pytest.approx(db.Pcs[0])
    assert not {"Tb", "Tc", "Pc", "omega"} & set(est.estimated)
    assert est.notes == []
    mine = PropertyPackage(["water", "my-ethanol"], "nrtl", custom=(c,))
    ref = get_package(("water", "ethanol"), "nrtl")
    t_mine = mine.saturation_point([0.0, 1.0], "bubble", P=ATM).T
    t_ref = ref.saturation_point([0.0, 1.0], "bubble", P=ATM).T
    assert t_mine == pytest.approx(t_ref, abs=0.5)


def test_antoine_constants_set_vapor_pressure_and_boiling_point():
    # Ethanol, log10(P/Pa) = 10.33 - 1642.89 / (T/K - 42.85) (NIST WebBook form, converted)
    c = custom(antoine={"A": 10.2370, "B": 1592.864, "C": -46.976},
               constants={"Tc": 514.7, "Pc": "62.68 bar"}, **ETHANOL_GROUPS)
    est = estimate(c)
    assert est.values["Tb"] == pytest.approx(351.4, abs=1.0)
    assert "vapor pressure" not in est.estimated
    pkg = PropertyPackage(["my-ethanol"], "ideal", custom=(c,))
    assert pkg.props.VaporPressures[0](est.values["Tb"]) == pytest.approx(ATM, rel=1e-3)


@pytest.mark.parametrize("method", ["peng-robinson", "srk", "ideal", "nrtl", "uniquac",
                                    "unifac", "wilson"])
def test_custom_compound_in_every_method(method):
    pkg = PropertyPackage(["water", "my-ethanol"], method, custom=(custom(**ETHANOL_GROUPS),))
    sat = pkg.saturation_point([0.5, 0.5], "bubble", P=ATM)
    assert 330.0 < sat.T < 360.0
    assert sat.y[1] > 0.5  # the alcohol is the lighter component
    assert any("my-ethanol is a user-defined compound" in w for w in pkg.warnings)


def custom_flowsheet(**compound):
    return {
        "components": ["water", "my-ethanol"], "property_method": "auto",
        "custom_compounds": [{"name": "my-ethanol", "formula": "C2H6O", **compound}],
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "T": "20 degC", "P": "1 atm",
                "flows": {"water": "60 kmol/h", "my-ethanol": "40 kmol/h"}}},
            {"id": "H", "type": "Heater", "params": {"T": "85 degC"}},
            {"id": "D", "type": "Flash", "params": {"duty": 0}},
        ],
        "streams": [{"id": "S1", "from": "F.out", "to": "H.in"},
                    {"id": "S2", "from": "H.out", "to": "D.in"},
                    {"id": "V", "from": "D.vapor"}, {"id": "L", "from": "D.liquid"}],
    }


def test_flowsheet_with_custom_compound_solves():
    r = solve(custom_flowsheet(**ETHANOL_GROUPS))
    assert r.ok, r.error
    assert r.property_method == "unifac"
    assert 0.0 < r.streams["V"]["molar_flow"] < r.streams["S1"]["molar_flow"]
    assert any("estimated" in w and "my-ethanol" in w for w in r.warnings)
    sources = [s for s in r.parameter_sources if s["pair"] == "my-ethanol"]
    assert any(s["source"].startswith("Tc: Joback") for s in sources)


def test_editing_a_definition_rebuilds_the_package():
    a = get_package(("water", "my-ethanol"), "peng-robinson",
                    custom=(custom(**ETHANOL_GROUPS),))
    b = get_package(("water", "my-ethanol"), "peng-robinson",
                    custom=(custom(constants={"Tc": 520.0}, **ETHANOL_GROUPS),))
    assert a is not b and b.consts.Tcs[1] == 520.0


@pytest.mark.parametrize(("definition", "message"), [
    ({"name": "x"}, "has no formula"),
    ({"name": "x", "formula": "C2H6Q"}, "formula"),
    ({"name": "x", "formula": "C2H6O", "joback_groups": {"-CH4": 1}}, "unknown Joback group"),
    ({"name": "x", "formula": "C3H8O", **ETHANOL_GROUPS}, "do not match the formula"),
    ({"name": "x", "formula": "C2H6O"}, "no normal boiling point"),
    ({"name": "x", "formula": "C2H6O", "constants": {"Tb": 350, "Tc": 340, "Pc": 6e6}},
     "must lie below"),
    ({"name": "x", "formula": "C2H6O", "constants": {"Tcrit": 500}}, "unknown constant"),
    ({"name": "x", "formula": "C2H6O", "joback_groups": {"-CH3": 1.5}}, "whole number"),
])
def test_incomplete_or_wrong_definitions_are_refused(definition, message):
    with pytest.raises(SimulationError, match=message):
        c = CustomCompound.from_dict(definition)
        PropertyPackage(["water", c.name], "peng-robinson", custom=(c,))


def test_unifac_needs_groups_for_a_custom_compound():
    c = custom(joback_groups=ETHANOL_GROUPS["joback_groups"])
    with pytest.raises(SimulationError, match="unifac_groups"):
        PropertyPackage(["water", "my-ethanol"], "unifac", custom=(c,))
    with pytest.raises(SimulationError, match="UNIQUAC"):
        PropertyPackage(["water", "my-ethanol"], "uniquac", custom=(c,))


def test_server_accepts_custom_compounds():
    client = TestClient(app)
    body = {"components": ["water", "my-ethanol"],
            "custom_compounds": [{"name": "my-ethanol", "formula": "C2H6O", **ETHANOL_GROUPS}]}
    r = client.post("/api/methods/suggest", json=body)
    assert r.status_code == 200 and r.json()["method"] == "unifac"
    d = client.post("/api/diagrams/binary", json={**body, "method": "unifac", "P": "1 atm",
                                                  "points": 7}).json()
    assert d["ok"], d
    assert suggest_method(["water", "ethanol"]).method == "unifac"  # unchanged without them


def test_custom_reaction_uses_joback_formation_enthalpy():
    pkg = PropertyPackage(["my-ethanol"], "ideal", custom=(custom(**ETHANOL_GROUPS),))
    assert pkg.formation_enthalpies[0] == pytest.approx(-236.84e3, rel=1e-3)
    assert math.isfinite(pkg.flash([1.0], T=400.0, P=ATM).H_reactive)
