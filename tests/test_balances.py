"""Checks use physical stream boundaries and identify conservation errors."""

from fluxsheet import solve
from fluxsheet.units.mixers import Splitter


def flowsheet():
    return {
        "components": ["propane", "ethane"],
        "property_method": "peng-robinson",
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "T": 300, "P": "2 bar", "flows": {"propane": 1, "ethane": 0.001}}},
            {"id": "S", "type": "Splitter", "params": {"split": 0.5}},
        ],
        "streams": [
            {"id": "feed", "from": "F.out", "to": "S.in"},
            {"id": "A", "from": "S.out1"},
            {"id": "B", "from": "S.out2"},
        ],
    }


def test_balances_report_components_and_energy():
    result = solve(flowsheet())
    assert result.ok, result.error
    for scope in ("S", "flowsheet"):
        checks = [c for c in result.balance_checks if c["scope"] == scope]
        assert {c["quantity"] for c in checks} == {"mass", "component", "energy"}
        assert {c["component"] for c in checks if c["quantity"] == "component"} == {
            "propane", "ethane"}
        assert all(c["status"] == "ok" for c in checks)
        assert all(c["unit"] and c["tolerance"] is not None for c in checks)
    assert result.to_dict()["balance_checks"] == result.balance_checks


def test_balance_error_identifies_unit_component_and_flowsheet(monkeypatch):
    original = Splitter.calculate

    def broken(self, pkg, inlets):
        outputs = original(self, pkg, inlets)
        outputs["out1"].flows[1] *= 0.5
        return outputs

    monkeypatch.setattr(Splitter, "calculate", broken)
    result = solve(flowsheet())
    assert result.ok
    failures = [c for c in result.balance_checks if c["status"] == "warning"]
    assert {c["scope"] for c in failures if c["component"] == "ethane"} == {
        "S", "flowsheet"}
    assert any("ethane" in warning for warning in result.warnings)


def test_zero_flow_branch_keeps_balances_valid():
    model = flowsheet()
    model["units"][1]["params"]["split"] = 0
    result = solve(model)
    assert result.ok, result.error
    assert all(c["status"] == "ok" for c in result.balance_checks)


def test_energy_error_is_reported_locally_and_globally(monkeypatch):
    from fluxsheet.stream import Stream

    original = Splitter.calculate

    def broken(self, pkg, inlets):
        outputs = original(self, pkg, inlets)
        stream = outputs["out1"]
        outputs["out1"] = Stream.from_flash(pkg, stream.flows, T=stream.T + 10, P=stream.P)
        return outputs

    monkeypatch.setattr(Splitter, "calculate", broken)
    result = solve(flowsheet())
    failures = [c for c in result.balance_checks
                if c["quantity"] == "energy" and c["status"] == "warning"]
    assert {c["scope"] for c in failures} == {"S", "flowsheet"}


def test_reactive_example_checks_elements_and_energy():
    import json
    from pathlib import Path

    model = json.loads((Path(__file__).parent.parent / "examples" / "ethyl_acetate.json")
                       .read_text(encoding="utf-8"))
    result = solve(model)
    assert result.ok, result.error
    for scope in ("CSTR", "PFR", "flowsheet"):
        checks = [c for c in result.balance_checks if c["scope"] == scope]
        assert any(c["quantity"] == "element" and c["component"] == "C" for c in checks)
        assert any(c["quantity"] == "energy" and c["status"] == "ok" for c in checks)
        assert all(c["status"] != "warning" for c in checks)


def test_property_crossing_is_accounted_for_in_global_energy():
    model = flowsheet()
    model["units"][1]["property_method"] = "ideal"
    result = solve(model)
    assert result.ok, result.error
    assert result.property_crossings
    energy = next(c for c in result.balance_checks
                  if c["scope"] == "flowsheet" and c["quantity"] == "energy")
    assert energy["status"] == "ok"


def test_property_crossing_into_virtual_boundary_does_not_change_global_energy():
    model = {
        "components": ["water", "ethanol"], "property_method": "nrtl",
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "T": "85 degC", "P": "2 bar", "flows": {"water": 4, "ethanol": 6}}},
            {"id": "M", "type": "Mult", "property_method": "ideal",
             "params": {"factor": 2}},
        ],
        "streams": [{"id": "in", "from": "F.out", "to": "M.in"},
                    {"id": "out", "from": "M.out"}],
    }
    result = solve(model)
    assert result.ok, result.error
    assert abs(result.property_crossings[0]["enthalpy_change"]) > 100
    energy = next(c for c in result.balance_checks
                  if c["scope"] == "flowsheet" and c["quantity"] == "energy")
    assert energy["status"] == "ok"


def test_property_crossing_from_virtual_boundary_is_included_in_global_energy():
    model = {
        "components": ["water", "ethanol"], "property_method": "nrtl",
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "T": "85 degC", "P": "2 bar", "flows": {"water": 4, "ethanol": 6}}},
            {"id": "M", "type": "Mult", "params": {"factor": 2}},
            {"id": "S", "type": "Splitter", "property_method": "ideal"},
        ],
        "streams": [{"id": "in", "from": "F.out", "to": "M.in"},
                    {"id": "cross", "from": "M.out", "to": "S.in"},
                    {"id": "A", "from": "S.out1"}, {"id": "B", "from": "S.out2"}],
    }
    result = solve(model)
    assert result.ok, result.error
    assert abs(result.property_crossings[0]["enthalpy_change"]) > 100
    energy = next(c for c in result.balance_checks
                  if c["scope"] == "flowsheet" and c["quantity"] == "energy")
    assert energy["status"] == "ok"


def test_property_crossing_into_product_is_outside_global_energy_boundary():
    model = {
        "components": ["water", "ethanol"], "property_method": "nrtl",
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "T": "85 degC", "P": "2 bar", "flows": {"water": 4, "ethanol": 6}}},
            {"id": "P", "type": "Product", "property_method": "ideal"},
        ],
        "streams": [{"id": "out", "from": "F.out", "to": "P.in"}],
    }
    result = solve(model)
    assert result.ok, result.error
    assert abs(result.property_crossings[0]["enthalpy_change"]) > 100
    energy = next(c for c in result.balance_checks
                  if c["scope"] == "flowsheet" and c["quantity"] == "energy")
    assert energy["status"] == "ok"


def test_missing_formation_enthalpy_is_marked_not_checkable():
    from types import SimpleNamespace

    from fluxsheet.balances import _energy, _record

    stream = SimpleNamespace(flows=[1.0], enthalpy_flow=10.0)
    checks = []
    _record(checks, "reactor", "energy", None, [_energy(stream, [None])], [0.0], "W", 0.1)
    assert checks[0]["status"] == "not_checkable"
    assert checks[0]["residual"] is None
