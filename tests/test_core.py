from pathlib import Path

import pytest

from fluxsheet import FlowsheetModel, solve
from fluxsheet.compounds import lookup, search
from fluxsheet.quantities import UnitError, to_si
from fluxsheet.thermo_package import suggest_method

EXAMPLES = Path(__file__).parent.parent / "examples"


def flowsheet(units, streams, components=("water", "ethanol"), method="unifac"):
    return {
        "components": list(components),
        "property_method": method,
        "units": units,
        "streams": streams,
    }


FEED = {"id": "F", "type": "Feed",
        "params": {"T": "25 degC", "P": "1 atm",
                   "flows": {"water": "80 kmol/h", "ethanol": "20 kmol/h"}}}


# --- quantities -------------------------------------------------------------------------

@pytest.mark.parametrize("value, quantity, expected", [
    ("25 degC", "temperature", 298.15),
    ("77 F", "temperature", 298.15),
    (300, "temperature", 300.0),
    ("1 atm", "pressure", 101325.0),
    ("1.5 bar", "pressure", 150000.0),
    ("36 kmol/h", "molar_flow", 10.0),
    ("2 kW", "power", 2000.0),
    ("40 %", "fraction", 0.4),
])
def test_to_si(value, quantity, expected):
    assert to_si(value, quantity) == pytest.approx(expected)


def test_to_si_rejects_wrong_unit():
    with pytest.raises(UnitError, match="unknown unit"):
        to_si("3 bar", "temperature")


# --- compounds and methods --------------------------------------------------------------

def test_lookup_by_name_and_cas():
    assert lookup("ethanol").cas == "64-17-5"
    assert lookup("64-17-5").name == "ethanol"


def test_search_puts_exact_match_first():
    results = search("methanol")
    assert results[0].cas == "67-56-1"
    assert results[0].simulation_ready


def test_method_suggestion():
    assert suggest_method(["water", "ethanol"]).method == "unifac"
    assert suggest_method(["methane", "propane"]).method == "peng-robinson"
    assert suggest_method(["water"]).method == "ideal"


# --- unit operations --------------------------------------------------------------------

def test_flash_conserves_mass_and_enriches_ethanol():
    r = solve(flowsheet(
        [FEED,
         {"id": "H", "type": "Heater", "params": {"T": "92 degC"}},
         {"id": "D", "type": "Flash"}],
        [{"id": "S1", "from": "F.out", "to": "H.in"},
         {"id": "S2", "from": "H.out", "to": "D.in"},
         {"id": "V", "from": "D.vapor"},
         {"id": "L", "from": "D.liquid"}],
    ))
    assert r.ok, r.error
    s = r.streams
    assert s["V"]["molar_flow"] + s["L"]["molar_flow"] == pytest.approx(s["S1"]["molar_flow"])
    assert s["V"]["mole_fractions"]["ethanol"] > s["S1"]["mole_fractions"]["ethanol"]
    assert s["V"]["T"] == pytest.approx(92 + 273.15)
    assert r.units["H"]["duty"] > 0


def test_adiabatic_mixer_conserves_energy():
    hot = {"id": "HOT", "type": "Feed",
           "params": {"T": "80 degC", "P": "1 atm", "flows": {"water": "10 kmol/h"}}}
    cold = {"id": "COLD", "type": "Feed",
            "params": {"T": "20 degC", "P": "1 atm", "flows": {"water": "10 kmol/h"}}}
    r = solve(flowsheet(
        [hot, cold, {"id": "M", "type": "Mixer"}],
        [{"id": "A", "from": "HOT.out", "to": "M.in"},
         {"id": "B", "from": "COLD.out", "to": "M.in"},
         {"id": "C", "from": "M.out"}],
        components=["water"], method="ideal",
    ))
    assert r.ok, r.error
    s = r.streams
    assert s["C"]["enthalpy_flow"] == pytest.approx(
        s["A"]["enthalpy_flow"] + s["B"]["enthalpy_flow"], rel=1e-6)
    assert s["C"]["T"] == pytest.approx(50 + 273.15, abs=1.5)  # Cp varies with T


def test_valve_is_isenthalpic():
    r = solve(flowsheet(
        [{"id": "F", "type": "Feed",
          "params": {"T": "20 degC", "P": "40 bar", "flows": {"methane": 1, "propane": 1}}},
         {"id": "V", "type": "Valve", "params": {"P": "5 bar"}}],
        [{"id": "S1", "from": "F.out", "to": "V.in"}, {"id": "S2", "from": "V.out"}],
        components=["methane", "propane"], method="peng-robinson",
    ))
    assert r.ok, r.error
    assert r.streams["S2"]["enthalpy"] == pytest.approx(r.streams["S1"]["enthalpy"], rel=1e-6)
    assert r.streams["S2"]["T"] < r.streams["S1"]["T"]  # Joule-Thomson cooling


def test_pump_needs_power():
    r = solve(flowsheet(
        [{"id": "F", "type": "Feed",
          "params": {"T": "25 degC", "P": "1 bar", "flows": {"water": "36 kmol/h"}}},
         {"id": "P", "type": "Pump", "params": {"P": "11 bar", "efficiency": 0.5}}],
        [{"id": "S1", "from": "F.out", "to": "P.in"}, {"id": "S2", "from": "P.out"}],
        components=["water"], method="ideal",
    ))
    assert r.ok, r.error
    # Ideal pump work for water: V * dP = 18e-6 m3/mol * 1e6 Pa * 10 mol/s = 180 W
    assert r.units["P"]["power"] == pytest.approx(180 / 0.5, rel=0.1)


def test_component_splitter():
    r = solve(flowsheet(
        [FEED, {"id": "C", "type": "ComponentSplitter",
                "params": {"fractions": {"ethanol": 0.95, "water": 0.05}}}],
        [{"id": "S1", "from": "F.out", "to": "C.in"},
         {"id": "TOP", "from": "C.top"}, {"id": "BOT", "from": "C.bottom"}],
    ))
    assert r.ok, r.error
    assert r.streams["TOP"]["component_flows"]["ethanol"] == pytest.approx(
        0.95 * r.streams["S1"]["component_flows"]["ethanol"])


# --- recycle ----------------------------------------------------------------------------

def test_recycle_example_converges_and_balances():
    r = solve(FlowsheetModel.load(EXAMPLES / "gas_recycle.json"))
    assert r.ok, r.error
    assert r.converged and r.tear_streams
    assert not [w for w in r.warnings if "mass balance" in w]
    s = r.streams
    out = s["PURGE"]["mass_flow"] + s["NGL-OUT"]["mass_flow"]
    assert out == pytest.approx(s["S1"]["mass_flow"], rel=1e-5)


@pytest.mark.parametrize("path", sorted(EXAMPLES.glob("*.json")), ids=lambda p: p.stem)
def test_examples_solve(path):
    r = solve(FlowsheetModel.load(path))
    assert r.ok, r.error


# --- helpful errors ---------------------------------------------------------------------

def test_overspecified_heater_explains_itself():
    r = solve(flowsheet(
        [FEED, {"id": "H", "type": "Heater", "params": {"T": 350, "duty": 1000}}],
        [{"id": "S1", "from": "F.out", "to": "H.in"}],
    ))
    assert not r.ok
    assert r.error["where"] == "H"
    assert "over-specified" in r.error["message"]
    assert r.error["hint"]


def test_unconnected_inlet_is_reported():
    r = solve(flowsheet([FEED, {"id": "H", "type": "Heater", "params": {"T": 350}}], []))
    assert not r.ok and r.error["where"] == "H"


def test_unknown_unit_in_feed_flows():
    feed = {**FEED, "params": {**FEED["params"], "flows": {"methanol": 1}}}
    r = solve(flowsheet([feed], [{"id": "S1", "from": "F.out"}]))
    assert not r.ok and "component list" in r.error["message"]


def test_valve_cannot_raise_pressure():
    r = solve(flowsheet(
        [FEED, {"id": "V", "type": "Valve", "params": {"P": "5 bar"}}],
        [{"id": "S1", "from": "F.out", "to": "V.in"}],
    ))
    assert not r.ok and "pump or compressor" in r.error["hint"]


def test_units_downstream_of_a_recycle_run_once_after_convergence(monkeypatch):
    # Units that do not feed back into the loop must never see unconverged loop states
    # (a column behind a reactor loop failed on its early, too-small feed before this).
    model = FlowsheetModel.load(EXAMPLES / "gas_recycle.json")
    model.units.append(model.units[0].model_copy(
        update={"id": "NGL_HEATER", "type": "Heater", "params": {"T": "20 degC"}, "ui": None}))
    for s in model.streams:
        if s.id == "NGL-OUT":
            s.target = "NGL_HEATER.in"
    model.streams.append(type(model.streams[0])(id="NGL-WARM", source="NGL_HEATER.out",
                                                target="NGL.in"))
    from fluxsheet.units.heat import Heater

    calls = []
    original = Heater.calculate

    def counting(self, pkg, inlets):
        calls.append(self.id)
        return original(self, pkg, inlets)

    monkeypatch.setattr(Heater, "calculate", counting)
    r = solve(model)
    assert r.ok, r.error
    assert r.iterations > 1
    assert calls.count("NGL_HEATER") == 1
    assert calls.count("CHILLER") == r.iterations


def test_methanol_example_meets_its_product_purities():
    # Guards the combined example: a fixed distillate rate that no longer matches the
    # column feed silently pushes methanol into the bottoms.
    r = solve(FlowsheetModel.load(EXAMPLES / "methanol_plant.json"))
    assert r.ok, r.error
    assert r.streams["MEOH"]["mole_fractions"]["methanol"] > 0.997
    assert r.streams["BOTTOMS"]["mole_fractions"]["water"] > 0.997
    sources = {s["pair"]: s["source"] for s in r.parameter_sources}
    assert sources["methanol/water"] == "ThermoML fit"
