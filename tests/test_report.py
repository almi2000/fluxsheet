import io
import json
import re
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from fluxsheet import report
from fluxsheet.cli import main
from fluxsheet.flowsheet import FlowsheetModel
from fluxsheet.server import app

openpyxl = pytest.importorskip("openpyxl")

ROOT = Path(__file__).parent.parent
FLASH = ROOT / "examples" / "ethanol_flash.json"
COLUMN = ROOT / "examples" / "benzene_toluene_purity.json"  # column profiles + a design spec


@pytest.fixture(scope="module")
def column():
    model = FlowsheetModel.load(COLUMN)
    result = report.solve_for_report(model)
    assert result.ok, result.error
    return model, result


@pytest.fixture(scope="module")
def workbook(column):
    model, result = column
    return openpyxl.load_workbook(io.BytesIO(report.to_xlsx(report.build_report(model, result))))


def _row(ws, label: str) -> list:
    for row in ws.iter_rows(values_only=True):
        if row[0] == label:
            return list(row)
    raise AssertionError(f"no row {label!r}")


def test_workbook_sheets(workbook):
    assert workbook.sheetnames == ["Summary", "Streams", "Units", "Profiles", "Design specs"]


def test_summary_names_flowsheet_and_method(workbook, column):
    model, result = column
    rows = {r[0]: r[1] for r in workbook["Summary"].iter_rows(values_only=True) if r[0]}
    assert rows["Flowsheet"] == model.name
    assert rows["Property method"] == result.property_method
    assert rows["Components"] == ", ".join(model.components)
    assert rows["Status"] == "solved"


def test_stream_table_in_display_units(workbook, column):
    _, result = column
    ws = workbook["Streams"]
    ids = list(result.streams)
    assert list(next(ws.iter_rows(values_only=True))) == ["Stream", "Unit", *ids]
    assert ws.freeze_panes == "C2"
    T = _row(ws, "Temperature")
    assert T[1] == "°C"
    assert T[2] == pytest.approx(result.streams[ids[0]]["T"] - 273.15)
    P = _row(ws, "Pressure")
    assert P[1] == "bar" and P[2] == pytest.approx(result.streams[ids[0]]["P"] / 1e5)
    flow = _row(ws, "Molar flow")
    assert flow[1] == "kmol/h"
    assert flow[2] == pytest.approx(result.streams[ids[0]]["molar_flow"] * 3.6)
    # Mole fractions of every stream with flow add up to one.
    comps = list(result.streams[ids[0]]["mole_fractions"])
    rows = list(ws.iter_rows(values_only=True))
    start = next(i for i, r in enumerate(rows) if r[0] == "Mole fractions")
    fracs = rows[start + 1:start + 1 + len(comps)]
    for col in range(2, 2 + len(ids)):
        if result.streams[ids[col - 2]]["molar_flow"] > 0:
            assert sum(r[col] for r in fracs) == pytest.approx(1.0)
    assert ws.cell(row=rows.index(tuple(T)) + 1, column=3).number_format == "0.00"


def test_unit_results_carry_units(workbook, column):
    _, result = column
    ws = workbook["Units"]
    rows = [r for r in ws.iter_rows(values_only=True)][1:]
    duty = next(r for r in rows if r[0] == "COLUMN" and r[2].lower().startswith("condenser"))
    assert duty[4] == "kW"
    assert duty[3] == pytest.approx(result.units["COLUMN"]["condenser_duty"] / 1e3)
    assert not any(r[2] == "Profiles" for r in rows)


def test_profiles_sheet(workbook, column):
    _, result = column
    ws = workbook["Profiles"]
    rows = list(ws.iter_rows(values_only=True))
    profile = next(iter(result.units["COLUMN"]["profiles"].values()))
    assert rows[0][0] == f"COLUMN — {profile['label']}"
    header = rows[1]
    assert header[0] == profile["x_label"]
    temps = [s for s in profile["series"] if s["quantity"] == "temperature"]
    if temps:
        col = header.index(f"{temps[0]['label']} [°C]")
        assert rows[2][col] == pytest.approx(temps[0]["values"][0] - 273.15)


def test_design_spec_sheet(workbook, column):
    _, result = column
    row = list(workbook["Design specs"].iter_rows(values_only=True))[1]
    spec = result.design_specs[0]
    assert row[0] == spec["id"] and row[8] == "yes"


def test_html_contains_stream_table(column):
    model, result = column
    page = report.to_html(report.build_report(model, result))
    assert page.startswith("<!doctype html>")
    assert "@media print" in page and "<script" not in page
    assert "Temperature" in page and "°C" in page
    for sid in result.streams:
        assert f"<th>{sid}</th>" in page


def test_user_text_is_escaped_and_not_a_formula():
    model = json.loads(FLASH.read_text())
    model["name"] = "=HYPERLINK(\"http://x\")<script>\x07"
    model = FlowsheetModel.model_validate(model)
    built = report.build_report(model, report.solve_for_report(model))
    assert "<script>" not in report.to_html(built)
    wb = openpyxl.load_workbook(io.BytesIO(report.to_xlsx(built)))
    cell = wb["Summary"]["B3"]
    assert cell.data_type == "s" and cell.value.startswith("=HYPERLINK")


def test_failed_solve_still_reports_the_error():
    model = json.loads(FLASH.read_text())
    model["components"] = ["water", "unobtainium"]
    model = FlowsheetModel.model_validate(model)
    result = report.solve_for_report(model)
    assert not result.ok
    built = report.build_report(model, result)
    assert list(built.sheets) == ["Summary"]
    assert "unobtainium" in report.to_html(built)
    wb = openpyxl.load_workbook(io.BytesIO(report.to_xlsx(built)))
    assert wb.sheetnames == ["Summary"]


def test_display_units_match_the_editor():
    """report.DISPLAY mirrors DISPLAY in web/src/quantities.ts."""
    source = (ROOT / "web" / "src" / "quantities.ts").read_text(encoding="utf-8")
    block = source[source.index("const DISPLAY"):]
    block = block[:block.index("};")]
    editor = {m[0]: (float(m[1]), float(m[2]), m[3])
              for m in re.findall(r"(\w+): \[([-\d.e]+), ([-\d.e]+), \"([^\"]*)\"\]", block)}
    assert editor, "could not read DISPLAY from quantities.ts"
    for quantity, (factor, offset, unit) in editor.items():
        assert report.DISPLAY[quantity] == (pytest.approx(factor), pytest.approx(offset), unit)


def test_filename_is_ascii():
    assert report.filename("Ethanol / Wasser — Flash", "xlsx") == "ethanol_wasser_flash_report.xlsx"
    assert report.filename("äöü", "html") == "flowsheet_report.html"


client = TestClient(app)


@pytest.mark.parametrize("fmt, kind", [
    ("xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
    ("html", "text/html; charset=utf-8"),
])
def test_api_report(fmt, kind):
    r = client.post(f"/api/report?format={fmt}", json=json.loads(FLASH.read_text()))
    assert r.status_code == 200
    assert r.headers["content-type"] == kind
    assert f'filename="ethanol-water_flash_report.{fmt}"' in r.headers["content-disposition"]
    if fmt == "xlsx":
        assert "Streams" in openpyxl.load_workbook(io.BytesIO(r.content)).sheetnames
    else:
        assert "Stream table" in r.text


def test_api_report_rejects_unknown_format():
    r = client.post("/api/report?format=pdf", json=json.loads(FLASH.read_text()))
    assert r.status_code == 422


def test_cli_report(tmp_path):
    xlsx, page = tmp_path / "out.xlsx", tmp_path / "out.html"
    assert main(["report", str(FLASH), "--xlsx", str(xlsx), "--html", str(page)]) == 0
    assert "Streams" in openpyxl.load_workbook(xlsx).sheetnames
    assert "Stream table" in page.read_text(encoding="utf-8")
    assert main(["report", str(FLASH)]) == 2


def test_small_flowsheet_formats_each_unit_result():
    """Few unit results (fewer rows than columns), each with its own number style."""
    model = FlowsheetModel.model_validate({
        "format": "fluxsheet/1", "name": "Tiny", "components": ["water"],
        "property_method": "ideal",
        "units": [{"id": "F", "type": "Feed",
                   "params": {"T": "25 degC", "P": "1 bar", "flows": {"water": "10 kmol/h"}}},
                  {"id": "H", "type": "Heater", "params": {"T": "60 degC"}}],
        "streams": [{"id": "S1", "from": "F.out", "to": "H.in"}, {"id": "S2", "from": "H.out"}],
    })
    result = report.solve_for_report(model)
    assert result.ok, result.error
    built = report.build_report(model, result)
    assert "S2" in report.to_html(built)
    ws = openpyxl.load_workbook(io.BytesIO(report.to_xlsx(built)))["Units"]
    duty = next(r for r in ws.iter_rows() if r[0].value == "H" and r[4].value == "kW")
    assert duty[3].number_format in ("#,##0.0000", "0.000E+00")
