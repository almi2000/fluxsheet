"""Different property methods in different parts of a flowsheet (property sections)."""

from pathlib import Path

import pytest

from fluxsheet import FlowsheetModel, solve
from fluxsheet.property_sections import _grouped
from fluxsheet.thermo_package import get_package

EXAMPLES = Path(__file__).parent.parent / "examples"


def ethanol_flash(heater_method=None, drum_method=None, T="88 degC", default="unifac"):
    """Feed -> Heater -> adiabatic Flash -> two products, methods per unit."""
    def unit(uid, type_, method, **params):
        u = {"id": uid, "type": type_, "params": params}
        if method:
            u["property_method"] = method
        return u

    return {
        "components": ["water", "ethanol"],
        "property_method": default,
        "units": [
            unit("F", "Feed", None, T="25 degC", P="1 atm",
                 flows={"water": "80 kmol/h", "ethanol": "20 kmol/h"}),
            unit("H", "Heater", heater_method, T=T),
            unit("DRUM", "Flash", drum_method),
            unit("TOP", "Product", None),
            unit("BOTTOM", "Product", None),
        ],
        "streams": [
            {"id": "S1", "from": "F.out", "to": "H.in"},
            {"id": "S2", "from": "H.out", "to": "DRUM.in"},
            {"id": "V", "from": "DRUM.vapor", "to": "TOP.in"},
            {"id": "L", "from": "DRUM.liquid", "to": "BOTTOM.in"},
        ],
    }


def strip_methods(d: dict) -> dict:
    return {k: v for k, v in d.items() if not k.startswith("property_method")}


def test_single_method_gives_identical_results():
    model = FlowsheetModel.load(EXAMPLES / "gas_recycle.json")
    base = solve(model)
    explicit = model.model_copy(deep=True)
    for u in explicit.units:
        u.property_method = "peng-robinson"  # the same as the flowsheet default
    again = solve(explicit)
    assert base.ok and again.ok
    assert again.property_crossings == []
    assert again.streams == base.streams
    assert {k: strip_methods(v) for k, v in again.units.items()} == \
        {k: strip_methods(v) for k, v in base.units.items()}
    assert again.iterations == base.iterations
    assert all(u["property_method"] == "peng-robinson" for u in base.units.values())
    assert base.units["SEP"]["property_method_source"] == "flowsheet"
    assert again.units["SEP"]["property_method_source"] == "unit"


def test_crossing_keeps_flows_T_P_and_reflashes_with_the_consumer_method():
    r = solve(ethanol_flash(heater_method="peng-robinson", drum_method="nrtl", T="85 degC"))
    assert r.ok, r.error
    assert r.units["H"]["property_method"] == "peng-robinson"
    assert r.units["DRUM"]["property_method"] == "nrtl"
    assert r.units["F"]["property_method"] == "unifac"  # flowsheet default
    assert r.units["TOP"]["property_method"] == "nrtl"  # products follow their feed
    assert r.units["TOP"]["property_method_source"] == "inlet"
    assert r.streams["S2"]["property_method"] == "peng-robinson"
    assert r.streams["V"]["property_method"] == "nrtl"

    crossings = {c["stream"]: c for c in r.property_crossings}
    assert set(crossings) == {"S1", "S2"}  # unifac -> PR, PR -> NRTL; products never cross
    c = crossings["S2"]
    s2 = r.streams["S2"]
    assert (c["from_method"], c["to_method"]) == ("peng-robinson", "nrtl")
    assert c["T"] == s2["T"] and c["P"] == s2["P"]
    assert c["vapor_fraction_before"] == pytest.approx(s2["vapor_fraction"])
    expected = get_package(("water", "ethanol"), "nrtl").flash(
        list(s2["mole_fractions"].values()), T=s2["T"], P=s2["P"])
    assert c["vapor_fraction_after"] == pytest.approx(expected.VF)
    assert c["enthalpy_change"] == pytest.approx((expected.H - s2["enthalpy"]) * s2["molar_flow"])

    # Mass is conserved component by component across the boundary and the drum.
    for comp in ("water", "ethanol"):
        out = r.streams["V"]["component_flows"][comp] + r.streams["L"]["component_flows"][comp]
        assert out == pytest.approx(r.streams["S1"]["component_flows"][comp], rel=1e-9)
    assert not [w for w in r.warnings if "mass balance" in w]

    # The adiabatic drum's energy balance uses NRTL: its products carry the re-flashed enthalpy.
    h_out = r.streams["V"]["enthalpy_flow"] + r.streams["L"]["enthalpy_flow"]
    assert h_out == pytest.approx(expected.H * s2["molar_flow"], rel=1e-6, abs=1e-3)


def test_crossing_that_changes_the_phase_is_warned_about():
    # Ideal (Raoult) says 88 °C is below the bubble point; UNIFAC flashes 40 % of it.
    r = solve(ethanol_flash(heater_method="ideal"))
    assert r.ok, r.error
    c = next(c for c in r.property_crossings if c["stream"] == "S2")
    assert c["vapor_fraction_before"] == pytest.approx(0.0, abs=1e-9)
    assert c["vapor_fraction_after"] > 0.3
    assert any(w.startswith("S2: switching from ideal (H) to unifac (DRUM)") for w in r.warnings)


def test_recycle_spanning_two_methods_converges():
    model = FlowsheetModel.load(EXAMPLES / "gas_recycle.json")
    for u in model.units:
        if u.id in ("SEP", "SPLIT", "COMP"):
            u.property_method = "srk"
    r = solve(model)
    assert r.ok, r.error
    assert r.converged and r.tear_streams
    crossed = {c["stream"] for c in r.property_crossings}
    assert {"S3", "R2"} <= crossed  # into the SRK section and back into the PR mixer
    assert not [w for w in r.warnings if "balance" in w]
    s = r.streams
    out = s["PURGE"]["mass_flow"] + s["NGL-OUT"]["mass_flow"]
    assert out == pytest.approx(s["S1"]["mass_flow"], rel=1e-5)
    # The loop really runs both methods: SRK and PR give a different NGL split than PR alone.
    base = solve(FlowsheetModel.load(EXAMPLES / "gas_recycle.json"))
    assert s["NGL-OUT"]["molar_flow"] != pytest.approx(base.streams["NGL-OUT"]["molar_flow"],
                                                       rel=1e-6)


@pytest.mark.parametrize("method, hint", [
    ("Peng Robinson", "did you mean 'peng-robinson'"),
    ("NRTL", "did you mean 'nrtl'"),
    ("auto", "only for the flowsheet default"),
    ("magic", "use one of: peng-robinson"),
])
def test_unknown_unit_method_is_a_clear_error(method, hint):
    r = solve(ethanol_flash(drum_method=method))
    assert not r.ok
    assert r.error["where"] == "DRUM"
    assert f"unknown property method {method!r}" in r.error["message"]
    assert hint in r.error["hint"]
    assert "flowsheet default (unifac)" in r.error["hint"]


def test_unit_method_round_trips_through_the_file_format():
    model = FlowsheetModel.model_validate(ethanol_flash(drum_method="nrtl"))
    text = model.to_json()
    assert text.count('"property_method": "nrtl"') == 1
    assert FlowsheetModel.model_validate_json(text) == model


def test_separation_units_get_a_suggestion_from_their_own_feed():
    r = solve(ethanol_flash())
    assert r.ok, r.error
    s = r.units["DRUM"]["property_method_suggestion"]
    assert s["method"] == "nrtl"  # water/ethanol have measured ChemSep NRTL parameters
    assert "water and ethanol" in s["reason"]
    assert "property_method_suggestion" not in r.units["H"]


def test_methanol_example_uses_nrtl_for_the_column_and_pr_in_the_loop():
    r = solve(FlowsheetModel.load(EXAMPLES / "methanol_plant.json"))
    assert r.ok, r.error
    methods = {uid: u["property_method"] for uid, u in r.units.items()}
    assert methods["REACTOR"] == methods["HP_SEP"] == "peng-robinson"
    assert methods["DEGASSER"] == methods["COLUMN"] == "nrtl"
    assert [c["stream"] for c in r.property_crossings] == ["S6"]
    assert r.units["HP_SEP"]["property_method_suggestion"]["method"] == "peng-robinson"
    assert r.units["COLUMN"]["property_method_suggestion"]["method"] == "nrtl"
    assert not [w for w in r.warnings if "balance" in w]


def test_pair_warnings_are_grouped():
    w = [f"NRTL: no binary parameters for {p}, and UNIFAC cannot estimate them; this pair is "
         "treated as an ideal solution. Use a method with data for it." for p in ("a/b", "a/c")]
    assert _grouped(w) == [(
        "NRTL: no binary parameters for 2 pairs (a/b, a/c), and UNIFAC cannot estimate them; "
        "these pairs are treated as ideal solutions. Use a method with data for them."
    )]
    assert _grouped(w[:1]) == w[:1]
    assert _grouped(["something else"]) == ["something else"]
