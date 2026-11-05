"""Reports of a solved flowsheet: an Excel workbook and a printable HTML page.

``build_report`` turns a model and its ``SolveResult`` into a format-neutral ``Report``
(summary rows plus titled tables); ``to_xlsx`` and ``to_html`` write it, so both formats
always show the same numbers. Values are converted from SI to the units the editor shows
(°C, bar, kmol/h, kg/h, kW, ...; see ``DISPLAY``), and every column or row names its unit.

The HTML page is self-contained (no scripts, no external files) and has print styles, so
the browser's "Print to PDF" gives a PDF report. Excel output needs ``openpyxl``
(``pip install fluxsheet[report]``).

A model with a ``costing`` section gets a Costing sheet from the report's own solve. Results of
sensitivity studies, optimizations and a dynamic run are added when given as ``runs``
(``ReportRuns``): the editor passes the results of its last runs that still match the
flowsheet, so the report shows what the user saw without repeating long calculations on the
server. ``compute_runs`` runs them all instead (the CLI's ``--with-runs``). Dynamic time
series go to the Excel workbook only; the HTML report names the run and points there.
"""

from __future__ import annotations

import html
import io
import math
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any

from pydantic import BaseModel, ConfigDict, Field

from . import __version__
from .compounds import UnknownCompoundError
from .errors import SimulationError
from .flowsheet import FlowsheetModel
from .solver import SolveResult, solve
from .units import UNIT_TYPES

# How SI values are shown: quantity -> (factor, offset, unit label); mirrors DISPLAY in
# web/src/quantities.ts. Quantities not listed are shown in SI.
DISPLAY: dict[str, tuple[float, float, str]] = {
    "temperature": (1.0, -273.15, "°C"),
    "pressure": (1e-5, 0.0, "bar"),
    "pressure_difference": (1e-5, 0.0, "bar"),
    "molar_flow": (3.6, 0.0, "kmol/h"),
    "mass_flow": (3600.0, 0.0, "kg/h"),
    "power": (1e-3, 0.0, "kW"),
    "volume": (1.0, 0.0, "m³"),
    "length": (1.0, 0.0, "m"),
    "area": (1.0, 0.0, "m²"),
    "time": (1.0, 0.0, "s"),
    "amount": (1e-3, 0.0, "kmol"),
    "molar_energy": (1e-3, 0.0, "kJ/mol"),
    "temperature_difference": (1.0, 0.0, "K"),
    "thermal_conductance": (1e-3, 0.0, "kW/K"),
    "heat_transfer_coefficient": (1.0, 0.0, "W/m²K"),
    "specific_heat": (1e-3, 0.0, "kJ/kg·K"),
    "velocity": (1.0, 0.0, "m/s"),
    "volumetric_flow": (3600.0, 0.0, "m³/h"),
    "specific_energy": (1e-3, 0.0, "kJ/kg"),
}

MAX_HTML_STREAMS = 8  # stream columns per printed table; wider tables are split


def display_unit(quantity: str | None) -> str:
    return DISPLAY[quantity][2] if quantity in DISPLAY else ""


def to_display(value: Any, quantity: str | None) -> Any:
    """An SI number in its display unit; anything that is not a finite number as it is."""
    if not _is_number(value) or quantity not in DISPLAY:
        return value
    factor, offset, _ = DISPLAY[quantity]
    return value * factor + offset


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _label(text: str, unit: str) -> str:
    return f"{text} [{unit}]" if unit else text


# ---- the report structure -----------------------------------------------------------------


@dataclass
class Table:
    """A titled table. ``kinds`` gives each column's number style: "T" (temperature, two
    decimals), "frac" (four decimals), "int", or "num" (four significant digits)."""

    title: str
    columns: list[str]
    rows: list[list[Any]]
    kinds: list[str] | None = None
    # Per-row number style instead (Units sheet: one value per row, of varying quantity).
    row_kinds: list[str] | None = None
    # Streams-style table: the first column holds row labels, values run across.
    label_column: bool = False
    # Table of the Streams sheet split in chunks for print; see to_html.
    chunk: bool = False
    # Shown in the HTML report too; long time series are in the Excel workbook only.
    html: bool = True


@dataclass
class Report:
    name: str
    description: str
    ok: bool
    summary: list[tuple[str, Any]]
    warnings: list[str]
    error: dict | None
    # Sheet name -> tables shown on it, in order.
    sheets: dict[str, list[Table]] = field(default_factory=dict)


# ---- results of runs the report includes ----------------------------------------------------
# They mirror the results of run_sensitivity, run_optimization and simulate (as the editor
# gets them from the API); unknown fields are ignored and sizes are capped, as a report has to
# stay a sensible size.

_Num = float | None
# As many as the runs can produce: study points (sensitivity.MAX_POINTS), dynamic output
# points (DynamicsModel.points).
MAX_STUDY_ROWS = 101 * 101
MAX_RUN_POINTS = 2001


class _Loose(BaseModel):
    model_config = ConfigDict(extra="ignore")


class RunError(_Loose):
    message: str = ""
    where: str | None = None
    hint: str | None = None


class SensitivityColumnRun(_Loose):
    label: str
    quantity: str | None = None


class SensitivityRowRun(_Loose):
    inputs: list[_Num] = Field(max_length=10)
    outputs: list[_Num] = Field(max_length=200)
    error: RunError | None = None


class SensitivityRun(_Loose):
    ok: bool
    id: str | None = None
    inputs: list[SensitivityColumnRun] = Field(default_factory=list, max_length=10)
    outputs: list[SensitivityColumnRun] = Field(default_factory=list, max_length=200)
    rows: list[SensitivityRowRun] = Field(default_factory=list, max_length=MAX_STUDY_ROWS)
    error: RunError | None = None


class OptimizationVariableRun(_Loose):
    label: str
    quantity: str | None = None
    lower: _Num = None
    upper: _Num = None
    start: _Num = None
    value: _Num = None
    at_bound: bool = False


class OptimizationConstraintRun(_Loose):
    expression: str
    relation: str
    value: _Num = None
    achieved: _Num = None
    satisfied: bool = True
    active: bool = False


class OptimizationObjectiveRun(_Loose):
    expression: str = ""
    value: _Num = None
    start: _Num = None


class OptimizationHistoryRun(_Loose):
    iteration: int
    objective: _Num = None
    x: list[_Num] = Field(default_factory=list, max_length=200)
    max_violation: _Num = None


class OptimizationRun(_Loose):
    ok: bool
    id: str | None = None
    sense: str | None = None
    method: str | None = None
    converged: bool | None = None
    message: str | None = None
    iterations: int | None = None
    evaluations: int | None = None
    failed_evaluations: int | None = None
    history: list[OptimizationHistoryRun] = Field(default_factory=list, max_length=5_000)
    objective: OptimizationObjectiveRun | None = None
    variables: list[OptimizationVariableRun] = Field(default_factory=list, max_length=200)
    constraints: list[OptimizationConstraintRun] = Field(default_factory=list, max_length=200)
    error: RunError | None = None


class DynamicSeriesRun(_Loose):
    key: str = ""
    owner: str = ""
    label: str
    quantity: str | None = None
    values: list[_Num] = Field(max_length=MAX_RUN_POINTS)


class DynamicRun(_Loose):
    ok: bool
    time: list[_Num] = Field(default_factory=list, max_length=MAX_RUN_POINTS)
    series: list[DynamicSeriesRun] = Field(default_factory=list, max_length=2_000)
    warnings: list[str] = Field(default_factory=list, max_length=200)
    error: RunError | None = None


class ReportRuns(_Loose):
    sensitivity: list[SensitivityRun] = Field(default_factory=list, max_length=20)
    optimizations: list[OptimizationRun] = Field(default_factory=list, max_length=20)
    dynamics: DynamicRun | None = None


def compute_runs(model: FlowsheetModel, dynamics: bool = True) -> ReportRuns:
    """Run every sensitivity study and optimization of ``model``, and its dynamic run, for a
    report. This can take long: each study point and iteration is a flowsheet solve."""
    from .dynamics.simulator import simulate
    from .optimization import run_optimization
    from .sensitivity import run_sensitivity

    def safe(run, *args) -> dict:
        try:
            return run(*args)
        except SimulationError as exc:
            return {"ok": False, "error": exc.to_dict()}
        except UnknownCompoundError as exc:
            return {"ok": False, "error": {"message": str(exc)}}

    return ReportRuns.model_validate({
        "sensitivity": [safe(run_sensitivity, model, s.id) for s in model.sensitivity or []],
        "optimizations": [safe(run_optimization, model, o.id) for o in model.optimizations or []],
        "dynamics": safe(simulate, model) if dynamics and model.dynamics is not None else None,
    })


def solve_for_report(model: FlowsheetModel) -> SolveResult:
    """Solve ``model``; errors become a failed result (as /api/solve reports them), so a
    report can still be written for a flowsheet that does not solve."""
    try:
        return solve(model)
    except SimulationError as exc:
        return SolveResult(ok=False, error=exc.to_dict())
    except UnknownCompoundError as exc:
        return SolveResult(ok=False, error={
            "message": str(exc), "where": None,
            "hint": "search for the compound and pick it from the list"})


def build_report(model: FlowsheetModel, result: SolveResult, now: datetime | None = None,
                 runs: ReportRuns | None = None) -> Report:
    now = now or datetime.now(timezone.utc)
    runs = runs or ReportRuns()
    summary: list[tuple[str, Any]] = [
        ("Flowsheet", model.name),
        ("Description", model.description or ""),
        ("Components", ", ".join(model.components)),
        ("Property method", result.property_method or model.property_method),
    ]
    if result.property_method_reason:
        summary.append(("Method choice", result.property_method_reason))
    summary.append(("Status", "solved" if result.ok else "failed"))
    if result.ok:
        recycle = (f"tear {', '.join(result.tear_streams)}; converged in {result.iterations} "
                   "iterations" if result.tear_streams else "no recycle loops")
        summary.append(("Recycle", recycle))
    included = [*(f"sensitivity study {r.id}" for r in runs.sensitivity),
                *(f"optimization {r.id}" for r in runs.optimizations),
                *(["dynamic run"] if runs.dynamics else [])]
    if included:
        summary.append(("Included runs", ", ".join(included)))
    summary += [("Generated", now.strftime("%Y-%m-%d %H:%M UTC")),
                ("Fluxsheet version", __version__)]

    report = Report(name=model.name, description=model.description or "", ok=result.ok,
                    summary=summary, warnings=list(result.warnings), error=result.error)
    sources = [[s.get("pair", ""), s.get("source", ""), s.get("method", "")]
               for s in result.parameter_sources]
    summary_tables = []
    if sources:
        summary_tables.append(Table("Binary parameter sources", ["Pair", "Source", "Method"],
                                    sources))
    if result.calculators:
        summary_tables.append(Table(
            "Calculator blocks", ["Calculator", "Assigned (SI units)", "Flowsheet solves"],
            [[c["id"], ", ".join(f"{k} = {_fmt_num(v)}" for k, v in c["assigned"].items())
              if isinstance(c.get("assigned"), dict) else str(c.get("assigned")), c["passes"]]
             for c in result.calculators], kinds=["", "", "int"]))
    report.sheets["Summary"] = summary_tables
    if not result.ok:
        return report

    report.sheets["Streams"] = [_stream_table(model, result)]
    report.sheets["Units"] = [_unit_table(result)]
    profiles = _profile_tables(result)
    if profiles:
        report.sheets["Profiles"] = profiles
    if result.design_specs:
        report.sheets["Design specs"] = [_spec_table(result)]
    if result.costing:
        report.sheets["Costing"] = _costing_tables(result.costing)
    if runs.sensitivity:
        report.sheets["Sensitivity"] = [t for r in runs.sensitivity for t in _sensitivity_tables(r)]
    if runs.optimizations:
        report.sheets["Optimization"] = [t for r in runs.optimizations
                                         for t in _optimization_tables(r)]
    if runs.dynamics:
        report.sheets["Dynamics"] = _dynamic_tables(runs.dynamics)
    return report


def _stream_table(model: FlowsheetModel, result: SolveResult) -> Table:
    ids = list(result.streams)
    ends = {s.id: (s.source, s.target or "") for s in model.streams}
    streams = [result.streams[i] for i in ids]
    components = list(streams[0]["component_flows"]) if streams else []

    def row(label: str, unit: str, kind: str, values: list[Any]) -> tuple[list[Any], str]:
        return [label, unit, *values], kind

    def conv(key: str, quantity: str | None) -> list[Any]:
        return [to_display(s.get(key), quantity) for s in streams]

    rows = [
        row("From", "", "", [ends.get(i, ("", ""))[0] for i in ids]),
        row("To", "", "", [ends.get(i, ("", ""))[1] for i in ids]),
        row("Phase", "", "", [s.get("phase", "") for s in streams]),
    ]
    if len({s.get("property_method") for s in streams}) > 1:
        rows.append(row("Property method", "", "", [s.get("property_method") for s in streams]))
    rows += [
        row("Temperature", "°C", "T", conv("T", "temperature")),
        row("Pressure", "bar", "num", conv("P", "pressure")),
        row("Vapor fraction", "", "frac", conv("vapor_fraction", None)),
        row("Molar flow", "kmol/h", "num", conv("molar_flow", "molar_flow")),
        row("Mass flow", "kg/h", "num", conv("mass_flow", "mass_flow")),
        row("Enthalpy flow", "kW", "num", conv("enthalpy_flow", "power")),
        row("Component flows", "", "", [""] * len(ids)),
    ]
    rows += [row(f"  {c}", "kmol/h", "num",
                 [to_display(s["component_flows"].get(c), "molar_flow") for s in streams])
             for c in components]
    rows.append(row("Mole fractions", "", "", [""] * len(ids)))
    rows += [row(f"  {c}", "", "frac",
                 [s["mole_fractions"].get(c) if (s.get("molar_flow") or 0) > 0 else None
                  for s in streams])
             for c in components]
    return Table("Stream table", ["Stream", "Unit", *ids], [r for r, _ in rows],
                 kinds=[k for _, k in rows], label_column=True, chunk=True)


# Results kept out of the Units sheet: the type is a column, profiles have their own sheet.
_SKIP_RESULTS = {"type", "profiles"}


def _unit_table(result: SolveResult) -> Table:
    rows, kinds = [], []
    for uid, unit in result.units.items():
        cls = UNIT_TYPES.get(unit.get("type"))
        declared = {f.name: f for f in cls.result_fields} if cls else {}
        names = [n for n in declared if n in unit]
        names += [n for n in unit if n not in declared and n not in _SKIP_RESULTS]
        for name in names:
            value = unit[name]
            if not (_is_number(value) or isinstance(value, (str, bool))):
                continue  # nested results (e.g. a method suggestion) are not scalars
            f = declared.get(name)
            quantity = f.quantity if f else None
            label = f.label if f else name.replace("_", " ").capitalize()
            rows.append([uid, unit.get("type", ""), label, to_display(value, quantity),
                         display_unit(quantity)])
            kinds.append("int" if isinstance(value, int) and not isinstance(value, bool) else
                         "T" if quantity == "temperature" else
                         "frac" if quantity == "fraction" else "num")
    return Table("Unit results", ["Unit", "Type", "Result", "Value", "Unit of measure"], rows,
                 row_kinds=kinds)


def _profile_tables(result: SolveResult) -> list[Table]:
    tables = []
    for uid, unit in result.units.items():
        for profile in (unit.get("profiles") or {}).values():
            series = profile.get("series", [])
            columns = [profile.get("x_label", "x")]
            columns += [_label(s["label"], display_unit(s.get("quantity"))) for s in series]
            rows = [[x, *(to_display(s["values"][i], s.get("quantity"))
                         if i < len(s["values"]) else None for s in series)]
                    for i, x in enumerate(profile.get("x", []))]
            kinds = ["num"] + ["T" if s.get("quantity") == "temperature" else
                               "frac" if s.get("quantity") in (None, "fraction") else "num"
                               for s in series]
            tables.append(Table(f"{uid} — {profile.get('label', '')}", columns, rows, kinds))
    return tables


def _spec_table(result: SolveResult) -> Table:
    rows = []
    for d in result.design_specs:
        q, vq = d.get("quantity"), d.get("vary_quantity")
        rows.append([d.get("id"), d.get("target"), to_display(d.get("value"), q),
                     to_display(d.get("achieved"), q), display_unit(q), d.get("vary"),
                     to_display(d.get("vary_value"), vq), display_unit(vq),
                     "yes" if d.get("converged") else "no", d.get("iterations")])
    return Table("Design specifications",
                 ["Spec", "Target", "Target value", "Achieved", "Unit", "Varied", "Value",
                  "Unit", "Converged", "Flowsheet solves"], rows,
                 kinds=["", "", "num", "num", "", "", "num", "", "", "int"])


def _costing_tables(costing: dict) -> list[Table]:
    summary = Table("Cost summary", ["Item", "Value", "Unit"], [
        ["Purchased equipment", costing.get("purchased_cost"), "$"],
        ["Bare module", costing.get("bare_module_cost"), "$"],
        ["Capital (bare module + 18 %)", costing.get("capital_cost"), "$"],
        ["Operating", costing.get("operating_cost"), "$/yr"],
        ["Annualization factor", costing.get("annualization_factor"), "1/yr"],
        ["Annualized capital", costing.get("annualized_capital"), "$/yr"],
        ["Total annualized cost (TAC)", costing.get("tac"), "$/yr"],
        ["CEPCI", costing.get("cepci"), ""],
        ["Operating hours", costing.get("hours_per_year"), "h/yr"],
    ], kinds=["", "num", ""])
    rows = []
    for uid, u in (costing.get("units") or {}).items():
        items = u.get("items") or []
        rows.append([
            uid, ", ".join(str(i.get("item", "")) for i in items),
            ", ".join(f"{_fmt_num(i.get('size'))} {i.get('size_unit', '')}"
                      + (f" ×{i['parallel']}" if (i.get("parallel") or 1) > 1 else "")
                      for i in items),
            u.get("purchased_cost"), u.get("bare_module_cost"),
            ", ".join(f"{x.get('utility')} "
                      f"{_fmt_num(to_display(abs(x.get('duty') or 0), 'power'))} kW"
                      for x in u.get("utilities") or []),
            u.get("operating_cost"), "; ".join(u.get("notes") or []),
        ])
    tables = [summary, Table("Equipment costs", [
        "Unit", "Equipment", "Size", "Purchased [$]", "Bare module [$]", "Utilities",
        "Operating [$/yr]", "Notes"], rows, kinds=["", "", "", "num", "num", "", "num", ""])]
    if costing.get("not_costed"):
        tables.append(Table("Not costed", ["Unit"], [[u] for u in costing["not_costed"]]))
    return tables


def _kind_of(quantity: str | None) -> str:
    return "T" if quantity == "temperature" else "frac" if quantity == "fraction" else "num"


def _error_table(title: str, error: RunError | None) -> Table:
    return Table(title, ["Error"], [[_error_text(error.model_dump()) if error else "failed"]])


def _sensitivity_tables(run: SensitivityRun) -> list[Table]:
    title = f"Sensitivity study {run.id or ''}".strip()
    if not run.ok:
        return [_error_table(title, run.error)]
    cols = [*run.inputs, *run.outputs]
    rows = []
    for row in run.rows:
        values = [*row.inputs, *row.outputs]
        rows.append([to_display(values[i], c.quantity) if i < len(values) else None
                     for i, c in enumerate(cols)]
                    + [row.error.message if row.error else ""])
    return [Table(title, [_label(c.label, display_unit(c.quantity)) for c in cols] + ["Error"],
                  rows, kinds=[_kind_of(c.quantity) for c in cols] + [""])]


def _optimization_tables(run: OptimizationRun) -> list[Table]:
    name = f"Optimization {run.id or ''}".strip()
    if not run.ok:
        status = "failed"
    elif run.converged:
        status = "converged"
    else:
        status = f"not converged: {run.message or ''}".strip()
    summary: list[list[Any]] = [
        ["Status", status], ["Method", run.method], ["Sense", run.sense],
        ["Iterations", run.iterations], ["Flowsheet solves", run.evaluations],
        ["Failed solves", run.failed_evaluations]]
    if run.objective:
        at = "optimum" if run.ok and run.converged else "best point" if run.ok else             "least violation"
        summary += [["Objective", run.objective.expression],
                    ["Objective at start", run.objective.start],
                    [f"Objective at {at}", run.objective.value]]
    if not run.ok and run.error:
        summary.append(["Error", _error_text(run.error.model_dump())])
    tables = [Table(name, ["Item", "Value"], [r for r in summary if r[1] is not None],
                    kinds=["", "num"])]
    if run.variables:
        tables.append(Table(f"{name}: variables", [
            "Variable", "Start", "Result", "Lower bound", "Upper bound", "Unit", "At bound"],
            [[v.label, to_display(v.start, v.quantity), to_display(v.value, v.quantity),
              to_display(v.lower, v.quantity), to_display(v.upper, v.quantity),
              display_unit(v.quantity), "yes" if v.at_bound else ""] for v in run.variables],
            kinds=["", "num", "num", "num", "num", "", ""]))
    if run.constraints:
        tables.append(Table(f"{name}: constraints", [
            "Constraint", "Relation", "Limit", "Achieved", "Status"],
            [[c.expression, c.relation, c.value, c.achieved,
              "violated" if not c.satisfied else "active" if c.active else "met"]
             for c in run.constraints], kinds=["", "", "num", "num", ""]))
    if run.history:
        variables = run.variables
        n = max(len(h.x) for h in run.history)
        labels = [_label(variables[k].label, display_unit(variables[k].quantity))
                  if k < len(variables) else f"x{k + 1}" for k in range(n)]
        quantity = [variables[k].quantity if k < len(variables) else None for k in range(n)]
        tables.append(Table(f"{name}: iteration history",
                            ["Iteration", "Objective", "Max. violation", *labels],
                            [[h.iteration, h.objective, h.max_violation,
                              *(to_display(h.x[k], quantity[k]) if k < len(h.x) else None
                                for k in range(n))] for h in run.history],
                            kinds=["int", "num", "num", *(_kind_of(q) for q in quantity)]))
    return tables


# Excel's last column is XFD (16 384); the time takes the first.
MAX_XLSX_SERIES = 16_383


def _dynamic_tables(run: DynamicRun) -> list[Table]:
    if not run.ok and not run.series:
        return [_error_table("Dynamic run", run.error)]
    series = run.series[:MAX_XLSX_SERIES]
    columns = ["Time [s]"] + [_label(f"{s.owner} {s.label}".strip(), display_unit(s.quantity))
                              for s in series]
    rows = [[t, *(to_display(s.values[i], s.quantity) if i < len(s.values) else None
                  for s in series)] for i, t in enumerate(run.time)]
    tables = [Table("Dynamic run", columns, rows,
                    kinds=["num", *(_kind_of(s.quantity) for s in series)], html=False)]
    notes = [[w] for w in run.warnings]
    if run.error:
        notes.insert(0, [f"stopped: {_error_text(run.error.model_dump())}"])
    if notes:
        tables.append(Table("Dynamic run: notes", ["Note"], notes))
    return tables


# ---- Excel ---------------------------------------------------------------------------------

_XLSX_FORMATS = {"T": "0.00", "frac": "0.0000", "int": "0", "num": "#,##0.0000"}
_FORMULA_START = ("=", "+", "-", "@", "\t", "\r")


def to_xlsx(report: Report) -> bytes:
    """The report as an .xlsx workbook (bytes)."""
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font, PatternFill
        from openpyxl.utils import get_column_letter
    except ImportError as exc:  # optional dependency
        raise RuntimeError("Excel reports need openpyxl: pip install fluxsheet[report]") from exc

    bold = Font(bold=True)
    head_fill = PatternFill("solid", fgColor="DDE5F0")
    title_font = Font(bold=True, size=13)
    wb = Workbook()
    wb.remove(wb.active)

    def put(ws, r: int, c: int, value: Any, kind: str = ""):
        cell = ws.cell(row=r, column=c)
        if _is_number(value):
            if math.isfinite(value):
                cell.value = value
                cell.number_format = _number_format(value, kind)
        elif isinstance(value, bool):
            cell.value = "yes" if value else "no"
        elif value is not None:
            text = _clean_text(str(value))
            cell.value = text
            if text.startswith(_FORMULA_START):
                cell.data_type = "s"  # user text, never a formula
        return cell

    def header(ws, r: int, columns: list[str]):
        for c, text in enumerate(columns, 1):
            cell = put(ws, r, c, text)
            cell.font, cell.fill = bold, head_fill
            cell.alignment = Alignment(wrap_text=True, vertical="top")

    # Summary
    ws = wb.create_sheet("Summary")
    put(ws, 1, 1, "Fluxsheet report").font = title_font
    r = 3
    for key, value in report.summary:
        put(ws, r, 1, key).font = bold
        put(ws, r, 2, value).alignment = Alignment(wrap_text=True, vertical="top")
        r += 1
    if report.error:
        r += 1
        put(ws, r, 1, "Error").font = bold
        put(ws, r, 2, _error_text(report.error))
        r += 1
    if report.warnings:
        r += 1
        put(ws, r, 1, "Warnings").font = bold
        for w in report.warnings:
            put(ws, r, 2, w)
            r += 1
    for table in report.sheets.get("Summary", []):
        r += 1
        put(ws, r, 1, table.title).font = bold
        r += 1
        header(ws, r, table.columns)
        for i, row in enumerate(table.rows):
            r += 1
            for c, value in enumerate(row, 1):
                put(ws, r, c, value, _kind(table, i, c - 1))
    _fit_columns(ws, get_column_letter, cap={2: 80})

    for sheet, tables in report.sheets.items():
        if sheet == "Summary":
            continue
        ws = wb.create_sheet(sheet)
        r = 1
        for table in tables:
            if len(tables) > 1:  # several blocks on one sheet get a title each
                put(ws, r, 1, table.title).font = bold
                r += 1
            header(ws, r, table.columns)
            if len(tables) == 1:
                ws.freeze_panes = ws.cell(row=r + 1, column=3 if table.label_column else 1)
            for i, row in enumerate(table.rows):
                kind = (table.kinds or [""] * len(table.rows))[i] if table.label_column else ""
                for c, value in enumerate(row, 1):
                    k = kind if table.label_column else _kind(table, i, c - 1)
                    cell = put(ws, r + 1 + i, c, value, k)
                    if table.label_column and c == 1:
                        cell.font = Font(bold=not str(value).startswith("  "))
            r += len(table.rows) + 2
        _fit_columns(ws, get_column_letter)

    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()


def _kind(table: Table, row: int, column: int) -> str:
    """Number style of a cell of a table without a label column."""
    if table.row_kinds is not None:
        return table.row_kinds[row] if row < len(table.row_kinds) else ""
    return table.kinds[column] if table.kinds and column < len(table.kinds) else ""


def _number_format(value: float, kind: str) -> str:
    if kind == "int" and float(value).is_integer():
        return "0"
    if kind in ("num", "") and value != 0 and not 1e-3 <= abs(value) < 1e7:
        return "0.000E+00"
    return _XLSX_FORMATS.get(kind, "General")


_CONTROL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")


def _clean_text(text: str) -> str:
    """Excel refuses control characters in cells."""
    return _CONTROL.sub("", text)


def _fit_columns(ws, letter, cap: dict[int, int] | None = None) -> None:
    widths: dict[int, int] = {}
    for row in ws.iter_rows():
        for cell in row:
            if cell.value is not None:
                n = len(str(cell.value)) if not _is_number(cell.value) else 12
                widths[cell.column] = max(widths.get(cell.column, 0), n)
    for col, n in widths.items():
        ws.column_dimensions[letter(col)].width = min(max(n + 2, 8), (cap or {}).get(col, 50))


def _error_text(error: dict) -> str:
    text = error.get("message", "")
    if error.get("where"):
        text = f"{error['where']}: {text}"
    if error.get("hint"):
        text += f" (hint: {error['hint']})"
    return text


# ---- HTML ----------------------------------------------------------------------------------

_CSS = """
:root { color-scheme: light; }
body { font: 13px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; color: #1c2430;
       margin: 24px auto; max-width: 1100px; padding: 0 16px; background: #fff; }
h1 { font-size: 22px; margin: 0 0 4px; }
h2 { font-size: 16px; margin: 28px 0 8px; border-bottom: 2px solid #2f6fb3; padding-bottom: 3px; }
h3 { font-size: 13px; margin: 16px 0 6px; }
.muted { color: #5b6675; }
table { border-collapse: collapse; margin: 0 0 12px; font-variant-numeric: tabular-nums; }
th, td { border: 1px solid #c9d1dc; padding: 3px 8px; vertical-align: top; }
thead th { background: #dde5f0; text-align: left; }
td.n { text-align: right; white-space: nowrap; }
th.row { text-align: left; font-weight: 600; background: #f3f6fa; white-space: nowrap; }
th.sub { font-weight: normal; padding-left: 20px; }
td.u { color: #5b6675; white-space: nowrap; }
dl.summary { display: grid; grid-template-columns: max-content 1fr; gap: 3px 16px; margin: 0; }
dl.summary dt { font-weight: 600; }
dl.summary dd { margin: 0; }
.error { border-left: 4px solid #c0392b; background: #fbeeee; padding: 8px 12px; }
.warnings li { color: #8a5a00; }
@page { size: A4 landscape; margin: 12mm; }
@media print {
  body { margin: 0; max-width: none; font-size: 10.5px; }
  h2 { break-before: page; }
  h2.first { break-before: avoid; }
  table, tr { break-inside: avoid; }
  thead { display: table-header-group; }
}
"""


def to_html(report: Report) -> str:
    """The report as a self-contained, printable HTML page."""
    e = html.escape
    out = [
        "<!doctype html>",
        '<html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        f"<title>{e(report.name)} — Fluxsheet report</title>",
        f"<style>{_CSS}</style></head><body>",
        f"<h1>{e(report.name)}</h1>",
    ]
    if report.description:
        out.append(f'<p class="muted">{e(report.description)}</p>')
    out.append('<h2 class="first">Summary</h2><dl class="summary">')
    for key, value in report.summary:
        out.append(f"<dt>{e(key)}</dt><dd>{e(str(value))}</dd>")
    out.append("</dl>")
    if report.error:
        out.append(f'<p class="error"><b>Error.</b> {e(_error_text(report.error))}</p>')
    if report.warnings:
        out.append('<h3>Warnings</h3><ul class="warnings">')
        out += [f"<li>{e(w)}</li>" for w in report.warnings]
        out.append("</ul>")
    for table in report.sheets.get("Summary", []):
        out.append(_html_table(table, heading="h3"))
    for sheet, tables in report.sheets.items():
        if sheet == "Summary":
            continue
        out.append(f"<h2>{e(_HTML_SECTION.get(sheet, sheet))}</h2>")
        for table in tables:
            if not table.html:
                out.append(f'<p class="muted">{e(table.title)}: {len(table.rows)} time points of '
                           f'{len(table.columns) - 1} variables, in the Excel report.</p>')
                continue
            out.append(_html_table(table, heading="h3" if len(tables) > 1 else None))
    out.append("</body></html>")
    return "\n".join(out)


_HTML_SECTION = {"Streams": "Stream table", "Units": "Unit results"}


def _html_table(table: Table, heading: str | None) -> str:
    e = html.escape
    parts = [f"<{heading}>{e(table.title)}</{heading}>"] if heading else []
    if table.label_column:
        # Properties down, streams across; split wide tables so each fits a printed page.
        n = len(table.columns) - 2
        step = MAX_HTML_STREAMS if table.chunk else max(n, 1)
        for start in range(0, max(n, 1), step):
            cols = list(range(2 + start, 2 + min(start + step, n)))
            parts.append("<table><thead><tr><th></th><th></th>")
            parts += [f"<th>{e(str(table.columns[c]))}</th>" for c in cols]
            parts.append("</tr></thead><tbody>")
            for i, row in enumerate(table.rows):
                label = str(row[0])
                cls = "row sub" if label.startswith("  ") else "row"
                kind = (table.kinds or [""] * len(table.rows))[i]
                parts.append(f'<tr><th class="{cls}">{e(label.strip())}</th>'
                             f'<td class="u">{e(str(row[1] or ""))}</td>')
                parts += [_html_cell(row[c], kind) for c in cols]
                parts.append("</tr>")
            parts.append("</tbody></table>")
        return "".join(parts)
    parts.append("<table><thead><tr>")
    parts += [f"<th>{e(c)}</th>" for c in table.columns]
    parts.append("</tr></thead><tbody>")
    for i, row in enumerate(table.rows):
        parts.append("<tr>" + "".join(_html_cell(v, _kind(table, i, c)) for c, v in enumerate(row))
                     + "</tr>")
    parts.append("</tbody></table>")
    return "".join(parts)


def _html_cell(value: Any, kind: str) -> str:
    if _is_number(value):
        return f'<td class="n">{html.escape(_fmt_num(value, kind))}</td>'
    if value is None:
        return '<td class="n">–</td>'
    return f"<td>{html.escape(str(value))}</td>"


def _fmt_num(value: Any, kind: str = "num") -> str:
    """Numbers as the editor shows them: °C with one decimal... here two, fractions with
    four decimals, otherwise four significant digits."""
    if not _is_number(value):
        return str(value)
    if not math.isfinite(value):
        return "–"
    if kind == "int" or isinstance(value, int):
        return str(int(value)) if float(value).is_integer() else f"{value:.4g}"
    if kind == "T":
        return f"{value:.2f}"
    if kind == "frac":
        return f"{value:.4f}"
    if value == 0:
        return "0"
    if not 1e-3 <= abs(value) < 1e7:
        return f"{value:.3e}"
    return f"{value:,.{max(0, 3 - math.floor(math.log10(abs(value))))}f}"


def filename(name: str, extension: str) -> str:
    """An ASCII file name for a report of the flowsheet ``name``."""
    stem = re.sub(r"[^A-Za-z0-9_-]+", "_", name).strip("_").lower()[:60] or "flowsheet"
    return f"{stem}_report.{extension}"
