"""The extension points that column, reactor and thermo work builds on."""

import math

import pytest

from fluxsheet import solve
from fluxsheet.errors import SimulationError
from fluxsheet.thermo_package import PropertyPackage
from fluxsheet.units import UNIT_TYPES, UnitOperation
from fluxsheet.units.base import IN, OUT, Param, ResultField, single_inlet

# --- registry -------------------------------------------------------------------------------


def test_registry_finds_all_builtin_units_in_palette_order():
    types = list(UNIT_TYPES)
    assert types[:2] == ["Feed", "Product"]
    assert {"Mixer", "Heater", "Flash", "Valve", "Pump", "Compressor"} <= set(types)
    assert "_Machine" not in types


def test_schema_includes_result_fields():
    schema = UNIT_TYPES["Flash"].schema()
    assert [f["name"] for f in schema["result_fields"]] == ["duty", "vapor_fraction"]
    assert UNIT_TYPES["Pump"].schema()["result_fields"][0]["quantity"] == "power"


# --- parameter kinds ------------------------------------------------------------------------


class _Demo(UnitOperation):
    # Not registered: defined after the registry was built. Used to exercise base features.
    type = "_Demo"
    label = "Demo"
    category = "Test"
    description = "test unit"
    inlets = (IN,)
    outlets = (OUT,)
    params = (
        Param("stages", "Stages", "integer", default=10, min=2),
        Param("condenser", "Condenser", "choice", default="total",
              options=("total", "partial", "none")),
        Param("key", "Key component", "component"),
        Param("reactions", "Reactions", "string_list"),
        Param("volume", "Volume", "volume"),
    )
    result_fields = (ResultField("x", "X"),)


def test_integer_choice_component_and_list_params():
    u = _Demo("D", {"stages": 12, "condenser": "partial", "key": "water",
                    "reactions": ["R1", "R2"], "volume": "250 L"})
    assert u.p == {"stages": 12, "condenser": "partial", "key": "water",
                   "reactions": ["R1", "R2"], "volume": pytest.approx(0.25)}
    assert isinstance(u.p["stages"], int)


@pytest.mark.parametrize("raw, message", [
    ({"stages": 2.5}, "whole number"),
    ({"stages": 1}, "at least 2"),
    ({"condenser": "magic"}, "not an option"),
    ({"reactions": "R1"}, "list of names"),
])
def test_bad_params_explain_themselves(raw, message):
    with pytest.raises(SimulationError, match=message):
        _Demo("D", raw)


def test_component_index_and_profiles():
    pkg = PropertyPackage(["water", "ethanol"], "ideal")
    u = _Demo("D", {})
    assert u.component_index(pkg, "ethanol") == 1
    with pytest.raises(SimulationError, match="not in the flowsheet"):
        u.component_index(pkg, "methanol")
    u.add_profile("temp", "Temperature profile", "Stage", [1, 2],
                  {"T": ("temperature", [350.0, 360.0])})
    assert u.profiles["temp"]["series"][0] == {
        "label": "T", "quantity": "temperature", "values": [350.0, 360.0]}
    s = single_inlet({"in": ["stream"]})
    assert s == "stream"


# --- property package -----------------------------------------------------------------------


def test_k_values_match_flash_equilibrium():
    pkg = PropertyPackage(["water", "ethanol"], "unifac")
    st = pkg.flash([0.8, 0.2], P=101325.0, VF=0.5)
    K = pkg.k_values(st.T, st.P, list(st.x), list(st.y))
    for k, x, y in zip(K, st.x, st.y):
        assert k == pytest.approx(y / x, rel=1e-4)


def test_phase_props_and_reactive_enthalpy():
    pkg = PropertyPackage(["methane", "water", "carbon dioxide", "hydrogen"], "peng-robinson")
    zs = [0.25, 0.25, 0.25, 0.25]
    vap = pkg.phase_props("vapor", 600.0, 1e5, zs)
    assert len(vap.lnphi) == 4 and all(math.isfinite(v) for v in vap.lnphi)
    hf = pkg.formation_enthalpies
    assert vap.H_reactive == pytest.approx(vap.H + sum(z * h for z, h in zip(zs, hf)))
    st = pkg.flash(zs, T=600.0, P=1e5)
    assert st.H_reactive == pytest.approx(st.H + sum(z * h for z, h in zip(zs, hf)))
    with pytest.raises(ValueError):
        pkg.phase_props("solid", 300.0, 1e5, zs)


def test_unit_results_carry_profiles_only_when_present():
    r = solve({
        "components": ["water"], "property_method": "ideal",
        "units": [
            {"id": "F", "type": "Feed",
             "params": {"T": "25 degC", "P": "1 bar", "flows": {"water": 1}}},
            {"id": "H", "type": "Heater", "params": {"T": "50 degC"}},
        ],
        "streams": [{"id": "S1", "from": "F.out", "to": "H.in"}, {"id": "S2", "from": "H.out"}],
    })
    assert r.ok, r.error
    assert "profiles" not in r.units["H"]
