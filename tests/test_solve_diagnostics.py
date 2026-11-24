import json

import pytest

from fluxsheet import solve
from fluxsheet.cli import main
from fluxsheet.units.heat import Heater


@pytest.fixture
def heating_flowsheet():
    return {
        "components": ["water"],
        "property_method": "ideal",
        "units": [
            {"id": "F", "type": "Feed", "params": {
                "T": "25 degC", "P": "1 atm", "flows": {"water": "36 kmol/h"},
            }},
            {"id": "H", "type": "Heater", "params": {"duty": "20 kW"}},
        ],
        "streams": [
            {"id": "inlet", "from": "F.out", "to": "H.in"},
            {"id": "outlet", "from": "H.out"},
        ],
    }


@pytest.mark.parametrize(("parameters", "expected"), [
    ({"duty": "20 kW"}, "ok"),
    ({"T": "50 degC"}, "not_applicable"),
])
def test_solve_reports_which_checks_apply(heating_flowsheet, parameters, expected):
    heating_flowsheet["units"][1]["params"] = parameters
    result = solve(heating_flowsheet).to_dict()

    assert result["ok"] and result["converged"]
    assert result["balance_status"] == "ok"
    assert result["specification_status"] == expected
    if expected == "not_applicable":
        assert result["specification_checks"] == []
    else:
        check, = result["specification_checks"]
        assert check["target"] == 20_000
        assert check["actual"] == pytest.approx(20_000, abs=0.1)


def test_failed_solve_does_not_claim_checks_ran(heating_flowsheet):
    heating_flowsheet["streams"] = []
    result = solve(heating_flowsheet).to_dict()

    assert not result["ok"] and not result["converged"]
    assert result["error"]["where"] == "H"
    assert result["specification_status"] == result["balance_status"] == "not_run"
    assert result["specification_checks"] == result["balance_checks"] == []


def test_cli_keeps_a_duty_deviation_distinct_from_the_balance(
    heating_flowsheet, tmp_path, monkeypatch, capsys,
):
    calculate = Heater.calculate

    def calculate_half_duty(unit, pkg, inlets):
        target = unit.p["duty"]
        unit.p["duty"] = target / 2
        try:
            return calculate(unit, pkg, inlets)
        finally:
            unit.p["duty"] = target

    monkeypatch.setattr(Heater, "calculate", calculate_half_duty)
    path = tmp_path / "heater.json"
    path.write_text(json.dumps(heating_flowsheet), encoding="utf-8")

    assert main(["run", str(path), "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["ok"] and result["converged"]
    assert result["balance_status"] == "ok"
    assert result["specification_status"] == "warning"
    check, = result["specification_checks"]
    assert check["residual"] == pytest.approx(-10_000, abs=0.1)

    assert main(["run", str(path)]) == 0
    output = capsys.readouterr().out
    assert "Calculation: converged | Duty specifications: warning | Balances: ok" in output
    assert "H: process duty specification target 20000 W, actual 10000 W, residual -10000 W" in output
