"""Fixed thermal duties are checked independently of equipment energy balances."""

import pytest

from fluxsheet import FlowsheetModel
from fluxsheet.balances import check_specifications
from fluxsheet.solver import Solver
from fluxsheet.units.heat import HeatExchanger


def _feed(uid, T):
    return {"id": uid, "type": "Feed", "params": {"T": T, "P": "1 atm",
                                                "flows": {"water": 2.0}}}


def _heater(params, method=None):
    unit = {"id": "H", "type": "Heater", "params": params}
    if method is not None:
        unit["property_method"] = method
    return {
        "components": ["water"], "property_method": "ideal",
        "units": [_feed("F", "55 degC"), unit],
        "streams": [{"id": "feed", "from": "F.out", "to": "H.in"}],
    }


def _exchanger(params):
    return {
        "components": ["water"], "property_method": "ideal",
        "units": [_feed("F1", "80 degC"), _feed("F2", "20 degC"),
                  {"id": "HX", "type": "HeatExchanger", "params": params}],
        "streams": [{"id": "cold", "from": "F2.out", "to": "HX.cold_in"},
                    {"id": "hot", "from": "F1.out", "to": "HX.hot_in"},
                    {"id": "chilled", "from": "HX.hot_out"}],
    }


def _solve(model):
    solver = Solver(FlowsheetModel.model_validate(model))
    result = solver.solve()
    assert result.ok, result.error
    return result, check_specifications(solver, solver.pkg, solver.streams)


@pytest.mark.parametrize("target", [3000.0, -2000.0, 0.0])
def test_heater_uses_the_local_inlet_enthalpy_across_property_methods(target):
    result, (check,) = _solve(_heater({"duty": target}, method="peng-robinson"))
    assert result.property_crossings[0]["stream"] == "feed"
    assert abs(result.property_crossings[0]["enthalpy_change"]) > check["tolerance"]
    assert check["scope"] == "H"
    assert check["quantity"] == "duty"
    assert check["side"] == "process"
    assert check["unit"] == "W"
    assert check["target"] == target
    assert check["actual"] == pytest.approx(target, abs=0.01)
    assert check["status"] == "ok"


def test_exchanger_detects_wrong_duty_even_with_a_closed_energy_balance(monkeypatch):
    original = HeatExchanger.calculate

    def miscalculate(unit, pkg, inlets):
        target = unit.p["duty"]
        try:
            unit.p["duty"] = target / 2
            return original(unit, pkg, inlets)
        finally:
            unit.p["duty"] = target

    monkeypatch.setattr(HeatExchanger, "calculate", miscalculate)
    result, checks = _solve(_exchanger({"duty": "3 kW"}))
    assert all(c["status"] == "ok" for c in result.balance_checks)
    assert result.units["HX"]["duty"] == pytest.approx(1500.0)
    assert {check["side"] for check in checks} == {"hot", "cold"}
    for check in checks:
        assert check["target"] == 3000.0
        assert check["actual"] == pytest.approx(1500.0)
        assert check["residual"] == pytest.approx(-1500.0)
        assert check["status"] == "warning"


def test_exchanger_checks_both_sides_with_a_positive_transfer_sign():
    _, checks = _solve(_exchanger({"duty": "3 kW"}))
    assert {check["side"] for check in checks} == {"hot", "cold"}
    for check in checks:
        assert check["scope"] == "HX"
        assert check["target"] == 3000.0
        assert check["actual"] == pytest.approx(3000.0, abs=0.01)
        assert check["status"] == "ok"


@pytest.mark.parametrize("make_model", [_heater, _exchanger])
def test_zero_flow_reports_zero_actual_duty_and_an_unmet_target(make_model):
    model = make_model({"duty": "1 kW"})
    model["units"].append({"id": "S", "type": "Splitter", "params": {"split": 0.0}})
    conn = model["streams"][0]
    model["streams"].append({"id": "supply", "from": conn["from"], "to": "S.in"})
    conn["from"] = "S.out1"
    _, checks = _solve(model)
    assert checks
    for check in checks:
        assert check["actual"] == pytest.approx(0.0, abs=0.01)
        assert check["residual"] == pytest.approx(-1000.0, abs=0.01)
        assert check["status"] == "warning"


@pytest.mark.parametrize("model", [
    _heater({"T": "60 degC"}),
    _exchanger({"UA": "50 W/K"}),
])
def test_calculated_duties_are_not_treated_as_fixed_targets(model):
    result, checks = _solve(model)
    assert any(unit.get("duty", 0) > 0 for unit in result.units.values())
    assert checks == []
