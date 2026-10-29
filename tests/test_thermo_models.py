"""Activity-coefficient models, vapor-liquid-liquid equilibrium, decanter and phase diagrams.

Reference values: azeotrope and heteroazeotrope data as tabulated in standard references
(e.g. the CRC Handbook azeotrope tables, Gmehling et al.), mutual solubilities of water and
1-butanol from IUPAC-NIST Solubility Data Series. Tolerances reflect the models, not the data.
"""

import pytest
from fastapi.testclient import TestClient

from fluxsheet import solve
from fluxsheet.diagrams import binary_diagram, phase_envelope
from fluxsheet.server import app
from fluxsheet.thermo_package import METHODS, get_package, suggest_method

ATM = 101325.0
C = 273.15


def decanter_flowsheet(method="nrtl", decanter_params=None, T="25 degC"):
    return {
        "components": ["water", "1-butanol"],
        "property_method": method,
        "units": [
            {"id": "F", "type": "Feed",
             "params": {"T": T, "P": "1 atm",
                        "flows": {"water": "70 kmol/h", "1-butanol": "30 kmol/h"}}},
            {"id": "D", "type": "Decanter", "params": decanter_params or {}},
        ],
        "streams": [
            {"id": "S1", "from": "F.out", "to": "D.in"},
            {"id": "L", "from": "D.light"},
            {"id": "H", "from": "D.heavy"},
            {"id": "V", "from": "D.vapor"},
        ],
    }


# --- methods and parameter provenance ---------------------------------------------------

def test_activity_methods_are_listed():
    assert {"nrtl", "uniquac", "wilson"} <= set(METHODS)
    assert "liquid-liquid" in METHODS["wilson"].description  # says what it cannot do


@pytest.mark.parametrize("method", ["nrtl", "uniquac", "wilson"])
def test_measured_parameters_are_used_without_warning(method):
    pkg = get_package(("water", "ethanol"), method)
    assert pkg.parameter_sources == [{"pair": "water/ethanol", "source": "ChemSep"}]
    assert pkg.warnings == []


@pytest.mark.parametrize("method", ["nrtl", "uniquac", "wilson"])
def test_missing_pair_is_estimated_from_unifac_and_flagged(method):
    # ChemSep has no water/n-hexane set for any of the three models.
    pkg = get_package(("water", "hexane"), method)
    assert pkg.parameter_sources == [{"pair": "water/hexane", "source": "UNIFAC estimate"}]
    assert len(pkg.warnings) == 1 and "estimated" in pkg.warnings[0]


def test_estimate_warning_reaches_the_solve_result():
    r = solve({
        "components": ["water", "hexane"], "property_method": "nrtl",
        "units": [{"id": "F", "type": "Feed",
                   "params": {"T": "25 degC", "P": "1 atm", "flows": {"water": 1, "hexane": 1}}}],
        "streams": [{"id": "S1", "from": "F.out"}],
    })
    assert r.ok, r.error
    assert any("water/hexane" in w and "UNIFAC" in w for w in r.warnings)


# --- vapor-liquid equilibrium -----------------------------------------------------------

def test_nrtl_water_ethanol_azeotrope():
    d = binary_diagram(["ethanol", "water"], "nrtl", P=ATM)
    az = d["azeotrope"]
    assert az["type"] == "homogeneous"
    assert az["x"] == pytest.approx(0.89, abs=0.025)  # measured 0.894
    assert az["T"] - C == pytest.approx(78.2, abs=0.2)  # measured 78.15 degC
    # Minimum-boiling: below both pure boiling points.
    assert az["T"] < min(d["T"][0], d["T"][-1])


@pytest.mark.parametrize("method", ["nrtl", "uniquac", "wilson", "unifac", "ideal"])
def test_txy_endpoints_are_pure_boiling_points(method):
    d = binary_diagram(["ethanol", "water"], method, P=ATM)
    assert d["x"][0] == 0.0 and d["x"][-1] == 1.0
    assert d["T"][0] == pytest.approx(373.124, abs=0.05)  # water
    assert d["T"][-1] == pytest.approx(351.57, abs=0.05)  # ethanol
    assert d["y"][0] == pytest.approx(0.0, abs=1e-6)
    assert d["y"][-1] == pytest.approx(1.0, abs=1e-6)


def test_pxy_endpoints_are_vapor_pressures():
    d = binary_diagram(["ethanol", "water"], "nrtl", T=C + 60)
    pkg = get_package(("ethanol", "water"), "nrtl")
    psat_ethanol, psat_water = (vp(C + 60) for vp in pkg.props.VaporPressures)
    assert d["P"][0] == pytest.approx(psat_water, rel=1e-4)
    assert d["P"][-1] == pytest.approx(psat_ethanol, rel=1e-4)
    assert psat_ethanol == pytest.approx(47.0e3, rel=0.03)  # steam-table class check


# --- liquid-liquid equilibrium ----------------------------------------------------------

def test_water_butanol_splits_at_25C_with_plausible_solubilities():
    pkg = get_package(("water", "1-butanol"), "nrtl")
    st = pkg.flash([0.7, 0.3], T=C + 25, P=ATM)
    assert st.liquid_split and st.VF == 0.0
    water_rich, butanol_rich = st.x_heavy, st.x_light
    assert st.rho_heavy > st.rho_light
    # Measured: about 1.9 mol% butanol in the water layer, 51 mol% water in the butanol layer.
    assert 0.003 < water_rich[1] < 0.03
    assert 0.45 < butanol_rich[0] < 0.65
    # The combined liquid is still the feed, as before the split was modelled.
    assert st.x == pytest.approx((0.7, 0.3))
    lever = st.light_fraction * butanol_rich[0] + st.heavy_fraction * water_rich[0]
    assert lever == pytest.approx(0.7, rel=1e-6)


def test_miscible_liquid_does_not_split():
    st = get_package(("water", "ethanol"), "nrtl").flash([0.5, 0.5], T=C + 25, P=ATM)
    assert not st.liquid_split and st.x_heavy is None and st.heavy_fraction == 0.0


@pytest.mark.parametrize("pair, measured_C", [(("water", "1-butanol"), 92.7),
                                              (("hexane", "water"), 61.6)])
def test_heteroazeotrope_temperatures(pair, measured_C):
    d = binary_diagram(list(pair), "nrtl", P=ATM)
    assert d["azeotrope"]["type"] == "heterogeneous"
    assert d["azeotrope"]["T"] - C == pytest.approx(measured_C, abs=1.0)
    assert d["liquid_liquid"]["curve"]  # solubility curve below the three-phase line


def test_stream_reports_second_liquid():
    r = solve(decanter_flowsheet())
    assert r.ok, r.error
    split = r.streams["S1"]["liquid_split"]
    assert set(split) == {"light", "heavy"}
    assert split["light"]["density"] < split["heavy"]["density"]
    assert split["light"]["fraction"] + split["heavy"]["fraction"] == pytest.approx(1.0)
    assert r.streams["L"]["liquid_split"] is None  # each layer is one phase


# --- decanter ---------------------------------------------------------------------------

def mass(stream):
    return stream["mass_flow"]


def test_decanter_conserves_mass_and_separates_layers():
    r = solve(decanter_flowsheet(decanter_params={"T": "25 degC"}))
    assert r.ok, r.error
    s = r.streams
    assert mass(s["L"]) + mass(s["H"]) + mass(s["V"]) == pytest.approx(mass(s["S1"]), rel=1e-9)
    for c in ("water", "1-butanol"):
        out = sum(s[k]["component_flows"][c] for k in ("L", "H", "V"))
        assert out == pytest.approx(s["S1"]["component_flows"][c], rel=1e-9)
    assert s["H"]["mole_fractions"]["water"] > 0.98
    assert s["L"]["mole_fractions"]["1-butanol"] > 0.35
    assert s["V"]["molar_flow"] == 0.0
    assert 0.0 < r.units["D"]["heavy_fraction"] < 1.0
    assert r.units["D"]["duty"] == pytest.approx(0.0, abs=1.0)  # feed already at 25 degC
    assert not [w for w in r.warnings if "balance" in w]


def test_adiabatic_decanter_conserves_energy():
    r = solve(decanter_flowsheet(T="60 degC"))
    assert r.ok, r.error
    s = r.streams
    h_out = sum(s[k]["enthalpy_flow"] for k in ("L", "H", "V"))
    assert h_out == pytest.approx(s["S1"]["enthalpy_flow"], rel=1e-6)
    assert s["L"]["T"] == pytest.approx(C + 60, abs=0.05)


def test_hot_decanter_vents_vapor():
    # Above the 92.7 degC heteroazeotrope at 1 atm, part of the feed boils off.
    r = solve(decanter_flowsheet(decanter_params={"T": "94 degC"}))
    assert r.ok, r.error
    assert r.units["D"]["vapor_fraction"] > 0 and r.streams["V"]["molar_flow"] > 0
    assert r.units["D"]["duty"] > 0


def test_decanter_needs_a_method_that_predicts_two_liquids():
    r = solve(decanter_flowsheet(method="peng-robinson"))
    assert not r.ok and r.error["where"] == "D"
    assert "NRTL" in r.error["hint"]


# --- method assistant -------------------------------------------------------------------

def test_suggestion_detects_liquid_liquid_split():
    s = suggest_method(["water", "1-butanol"])
    assert s.method == "nrtl"
    assert "two liquid phases" in s.reason and "ChemSep" in s.reason


def test_suggestion_without_measured_parameters_falls_back_to_unifac():
    s = suggest_method(["water", "hexane"])
    assert s.method == "unifac" and "water/hexane" in s.reason


def test_suggestion_keeps_miscible_polar_mixtures_on_unifac():
    assert suggest_method(["water", "ethanol"]).method == "unifac"
    assert suggest_method(["methane", "propane"]).method == "peng-robinson"


# --- phase envelope ---------------------------------------------------------------------

def test_natural_gas_phase_envelope():
    e = phase_envelope(["methane", "ethane", "propane", "n-butane"], [0.6, 0.2, 0.12, 0.08],
                       "peng-robinson")
    assert len(e["bubble"]) > 10 and len(e["dew"]) > 10
    # At any pressure the dew point lies above the bubble point.
    dew = {round(p["P"]): p["T"] for p in e["dew"]}
    for p in e["bubble"][:10]:
        assert dew[round(p["P"])] > p["T"]
    assert e["cricondentherm"]["T"] > max(p["T"] for p in e["bubble"])
    assert 50e5 < e["cricondenbar"]["P"] < 150e5


def test_envelope_rejects_activity_models():
    with pytest.raises(Exception, match="equation of state"):
        phase_envelope(["water", "ethanol"], [0.5, 0.5], "nrtl")


# --- API --------------------------------------------------------------------------------

client = TestClient(app)


def test_binary_diagram_endpoint_accepts_units_and_auto():
    body = client.post("/api/diagrams/binary", json={
        "components": ["ethanol", "water"], "method": "nrtl", "P": "1 atm", "points": 11,
    }).json()
    assert body["ok"] and body["kind"] == "txy" and len(body["T"]) == 11
    body = client.post("/api/diagrams/binary", json={
        "components": ["water", "1-butanol"], "T": "80 degC", "points": 11,
    }).json()
    assert body["ok"] and body["kind"] == "pxy" and body["method"] == "nrtl"


def test_diagram_endpoint_explains_errors():
    body = client.post("/api/diagrams/binary", json={
        "components": ["ethanol", "water", "methanol"], "method": "nrtl", "P": "1 atm",
    }).json()
    assert not body["ok"] and "two components" in body["error"]["message"]
    assert body["error"]["hint"]
    body = client.post("/api/diagrams/envelope", json={
        "components": ["methane", "propane"], "fractions": [0.5, 0.5], "method": "unifac",
    }).json()
    assert not body["ok"] and "Peng-Robinson" in body["error"]["hint"]
