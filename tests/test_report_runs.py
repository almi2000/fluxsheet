"""Reports with costing and the results of sensitivity studies, optimizations and dynamic
runs (report.ReportRuns)."""

import io
import json

import pytest
from fastapi.testclient import TestClient

from fluxsheet import report
from fluxsheet.cli import main
from fluxsheet.flowsheet import FlowsheetModel
from fluxsheet.server import app

openpyxl = pytest.importorskip("openpyxl")
client = TestClient(app)

TINY = {
    "format": "fluxsheet/1", "name": "Tiny heater", "components": ["water"],
    "property_method": "ideal",
    "units": [{"id": "F", "type": "Feed",
               "params": {"T": "25 degC", "P": "1 bar", "flows": {"water": "10 kmol/h"}}},
              {"id": "H", "type": "Heater", "params": {"T": "60 degC"}}],
    "streams": [{"id": "S1", "from": "F.out", "to": "H.in"}, {"id": "S2", "from": "H.out"}],
    "sensitivity": [{"id": "HOT", "vary": [{"unit": "H", "param": "T", "lower": "40 degC",
                                            "upper": "80 degC", "points": 3}],
                     "outputs": [{"unit": "H", "property": "duty"}]}],
}

SENSITIVITY = {
    "ok": True, "id": "HOT",
    "inputs": [{"label": "H T", "quantity": "temperature"}],
    "outputs": [{"label": "H duty", "quantity": "power"}, {"label": "x", "quantity": "fraction"}],
    "rows": [{"inputs": [313.15], "outputs": [12000.0, 0.5], "error": None},
             {"inputs": [353.15], "outputs": [None, None],
              "error": {"message": "flash failed", "where": "H", "hint": None}}],
}

OPTIMIZATION = {
    "ok": True, "id": "OPT", "sense": "minimize", "method": "SLSQP", "converged": True,
    "message": "Optimization terminated successfully", "iterations": 2, "evaluations": 7,
    "failed_evaluations": 1,
    "history": [{"iteration": 0, "objective": 120.0, "x": [323.15], "max_violation": 0.0},
                {"iteration": 1, "objective": 100.0, "x": [333.15], "max_violation": 0.2},
                {"iteration": 2, "objective": 95.5, "x": [343.15], "max_violation": 0.0}],
    "objective": {"expression": "tac", "value": 95.5, "start": 120.0},
    "variables": [{"label": "H T", "quantity": "temperature", "unit": "H", "param": "T",
                   "component": None, "lower": 313.15, "upper": 343.15, "start": 323.15,
                   "value": 343.15, "at_bound": True}],
    "constraints": [{"expression": "duty", "relation": "<=", "value": 5e4, "achieved": 5e4,
                     "satisfied": True, "active": True},
                    {"expression": "x", "relation": ">=", "value": 0.9, "achieved": 0.95,
                     "satisfied": True, "active": False}],
    "apply": [], "result": {"anything": "ignored"}, "new_field": 1,
}

DYNAMICS = {
    "ok": True, "time": [0.0, 60.0, 120.0], "warnings": ["controller LC saturated"],
    "series": [{"key": "S2.T", "group": "Streams", "owner": "S2", "label": "T",
                "quantity": "temperature", "values": [333.15, 334.15, 335.15]},
               {"key": "TANK.level", "group": "Units", "owner": "TANK", "label": "Level",
                "quantity": "length", "values": [1.0, 1.1, None]}],
}


def _runs(**kw) -> report.ReportRuns:
    return report.ReportRuns.model_validate(kw)


@pytest.fixture(scope="module")
def tiny():
    model = FlowsheetModel.model_validate(TINY)
    result = report.solve_for_report(model)
    assert result.ok, result.error
    return model, result


@pytest.fixture(scope="module")
def built(tiny):
    model, result = tiny
    return report.build_report(model, result, runs=_runs(
        sensitivity=[SENSITIVITY], optimizations=[OPTIMIZATION], dynamics=DYNAMICS))


@pytest.fixture(scope="module")
def workbook(built):
    return openpyxl.load_workbook(io.BytesIO(report.to_xlsx(built)))


def _rows(ws) -> list[list]:
    return [list(r) for r in ws.iter_rows(values_only=True)]


def _find(ws, first: str) -> list:
    for row in _rows(ws):
        if row[0] == first:
            return row
    raise AssertionError(f"no row {first!r} in {ws.title}")


def test_sheets_and_summary(workbook):
    assert workbook.sheetnames[-3:] == ["Sensitivity", "Optimization", "Dynamics"]
    assert _find(workbook["Summary"], "Included runs")[1] == \
        "sensitivity study HOT, optimization OPT, dynamic run"


def test_sensitivity_sheet_in_display_units(workbook):
    ws = workbook["Sensitivity"]
    assert _find(ws, "H T [°C]") == ["H T [°C]", "H duty [kW]", "x", "Error"]
    first = _find(ws, pytest.approx(40.0))
    assert first[1] == pytest.approx(12.0) and first[2] == 0.5 and first[3] is None
    failed = _find(ws, pytest.approx(80.0))
    assert failed[1] is None and failed[3] == "flash failed"


def test_optimization_sheet(workbook):
    ws = workbook["Optimization"]
    assert _find(ws, "Status")[1] == "converged"
    assert _find(ws, "Objective at optimum")[1] == 95.5
    assert _find(ws, "Failed solves")[1] == 1
    var = _find(ws, "H T")
    assert var[1:7] == [pytest.approx(50.0), pytest.approx(70.0), pytest.approx(40.0),
                        pytest.approx(70.0), "°C", "yes"]
    assert _find(ws, "duty")[4] == "active"
    assert _find(ws, "x")[4] == "met"
    assert _find(ws, "Iteration") == ["Iteration", "Objective", "Max. violation", "H T [°C]",
                                      None, None, None]
    assert _find(ws, 1)[:4] == [1, 100.0, 0.2, pytest.approx(60.0)]


def test_dynamics_sheet_in_excel_only(workbook, built):
    rows = _rows(workbook["Dynamics"])
    head = rows.index(_find(workbook["Dynamics"], "Time [s]"))
    assert rows[head][:3] == ["Time [s]", "S2 T [°C]", "TANK Level [m]"]
    assert rows[head + 2][:3] == [60.0, pytest.approx(61.0), 1.1]
    assert rows[head + 3][2] is None
    ws = workbook["Dynamics"]
    assert _find(ws, "controller LC saturated")
    page = report.to_html(built)
    assert "3 time points of 2 variables, in the Excel report" in page
    assert "334.15" not in page and "61.00" not in page


def test_html_sections(built):
    page = report.to_html(built)
    for text in ("<h2>Sensitivity</h2>", "<h2>Optimization</h2>", "Optimization OPT: constraints",
                 "Optimization OPT: iteration history", "flash failed", "H duty [kW]"):
        assert text in page


def test_failed_runs_show_their_error(tiny):
    model, result = tiny
    built = report.build_report(model, result, runs=_runs(
        sensitivity=[{"ok": False, "id": "BAD", "error": {"message": "no such unit X"}}],
        optimizations=[{**OPTIMIZATION, "ok": False, "converged": False,
                        "error": {"message": "no feasible point found", "where": "OPT"}}],
        dynamics={"ok": False, "error": {"message": "no dynamics"}}))
    wb = openpyxl.load_workbook(io.BytesIO(report.to_xlsx(built)))
    assert "no such unit X" in json.dumps(_rows(wb["Sensitivity"]))
    assert _find(wb["Optimization"], "Status")[1] == "failed"
    assert _find(wb["Optimization"], "Objective at least violation")[1] == 95.5
    assert "OPT: no feasible point found" in _find(wb["Optimization"], "Error")[1]
    assert "no dynamics" in json.dumps(_rows(wb["Dynamics"]))


def test_costing_sheet():
    model = FlowsheetModel.model_validate({**TINY, "costing": {}})
    result = report.solve_for_report(model)
    assert result.ok and result.costing, result.error
    built = report.build_report(model, result)
    ws = openpyxl.load_workbook(io.BytesIO(report.to_xlsx(built)))["Costing"]
    assert _find(ws, "Total annualized cost (TAC)")[1:3] == [pytest.approx(result.costing["tac"]),
                                                            "$/yr"]
    assert _find(ws, "Unit")[:4] == ["Unit", "Equipment", "Size", "Purchased [$]"]
    assert "Cost summary" in report.to_html(built)


def test_runs_validation_caps_and_ignores_extras():
    runs = _runs(optimizations=[OPTIMIZATION])
    assert not hasattr(runs.optimizations[0], "result")
    with pytest.raises(ValueError):
        _runs(sensitivity=[SENSITIVITY] * 21)
    # As many points as a study and a dynamic run can have are accepted.
    _runs(dynamics={**DYNAMICS, "time": [0.0] * 2001})


def test_api_takes_the_model_or_the_model_with_runs():
    old = client.post("/api/report?format=html", json=TINY)
    assert old.status_code == 200 and "Sensitivity" not in old.text
    new = client.post("/api/report?format=xlsx", json={
        "model": TINY, "runs": {"sensitivity": [SENSITIVITY], "dynamics": DYNAMICS}})
    assert new.status_code == 200
    assert openpyxl.load_workbook(io.BytesIO(new.content)).sheetnames[-2:] == \
        ["Sensitivity", "Dynamics"]
    assert client.post("/api/report", json={"model": TINY}).status_code == 200


@pytest.mark.parametrize("body", [
    {"model": {**TINY, "units": "nope"}},
    {"model": TINY, "runs": {"sensitivity": [{"rows": []}]}},  # "ok" missing
    {"name": "no units"},
])
def test_api_rejects_invalid_bodies(body):
    r = client.post("/api/report", json=body)
    assert r.status_code == 422
    assert isinstance(r.json()["detail"], list)


def test_cli_with_runs(tmp_path):
    src, xlsx = tmp_path / "tiny.json", tmp_path / "tiny.xlsx"
    src.write_text(json.dumps(TINY))
    assert main(["report", str(src), "--xlsx", str(xlsx), "--with-runs"]) == 0
    ws = openpyxl.load_workbook(xlsx)["Sensitivity"]
    assert [r[0] for r in _rows(ws)[1:]] == [pytest.approx(40.0), pytest.approx(60.0),
                                             pytest.approx(80.0)]
    assert main(["report", str(src), "--xlsx", str(xlsx)]) == 0
    assert "Sensitivity" not in openpyxl.load_workbook(xlsx).sheetnames
