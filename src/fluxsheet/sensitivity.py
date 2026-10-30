"""Sensitivity studies: solve the flowsheet over a range of one or two parameters.

Each point is a complete solve (with the design specs, if any), so a study answers questions
like "how do purity and reboiler duty change with the reflux ratio?". A point that fails is
recorded with its error and the study goes on.
"""

from __future__ import annotations

import csv
import io
import itertools

import numpy as np

from .errors import SimulationError
from .flowsheet import FlowsheetModel
from .solver import solve
from .specs import Handle, Variable

MAX_POINTS = 400  # all combinations of a two-parameter study


def run_sensitivity(model: FlowsheetModel, study_id: str | None = None) -> dict:
    """Run one study of ``model`` (the first if ``study_id`` is None).

    Returns ``{"ok": True, "id", "inputs", "outputs", "rows"}``: ``inputs`` and ``outputs`` hold
    a label and quantity per column, each row ``{"inputs", "outputs", "error"}`` holds SI
    values (``None`` where a point failed). On a setup error: ``{"ok": False, "error"}``.
    """
    try:
        study = _pick(model, study_id)
        handles = [Handle.resolve(v, model, study.id) for v in study.vary]
        outputs = [Variable.resolve(o, model, study.id) for o in study.outputs]
        grids = [np.linspace(h.lower, h.upper, v.points).tolist()
                 for h, v in zip(handles, study.vary)]
        total = int(np.prod([len(g) for g in grids]))
        if total > MAX_POINTS:
            raise SimulationError(f"the study has {total} points, more than {MAX_POINTS}",
                                  where=study.id, hint="use fewer points per parameter")
    except SimulationError as exc:
        return {"ok": False, "error": exc.to_dict()}

    rows = []
    for combo in itertools.product(*grids):
        trial = model
        for h, x in zip(handles, combo):
            trial = h.apply(trial, x)
        result = solve(trial)
        values, error = [None] * len(outputs), result.error
        if result.ok:
            try:
                values = [o.read(result) for o in outputs]
            except SimulationError as exc:
                error = exc.to_dict()
        rows.append({"inputs": list(combo), "outputs": values, "error": error})
    return {
        "ok": True,
        "id": study.id,
        "inputs": [{"label": h.label, "quantity": h.quantity} for h in handles],
        "outputs": [{"label": o.label, "quantity": o.quantity} for o in outputs],
        "rows": rows,
    }


def _pick(model: FlowsheetModel, study_id: str | None):
    studies = model.sensitivity or []
    if not studies:
        raise SimulationError("the flowsheet has no sensitivity study",
                              hint="add one under 'sensitivity' (in the editor: flowsheet "
                              "settings, Sensitivity)")
    if study_id is None:
        return studies[0]
    for s in studies:
        if s.id == study_id:
            return s
    raise SimulationError(f"there is no sensitivity study {study_id!r}",
                          hint=f"choose one of: {', '.join(s.id for s in studies)}")


def to_csv(table: dict) -> str:
    """The rows of a finished study as CSV, SI units, one column per input and output."""
    out = io.StringIO()
    writer = csv.writer(out, lineterminator="\n")
    cols = table["inputs"] + table["outputs"]
    writer.writerow([_header(c) for c in cols] + ["error"])
    for row in table["rows"]:
        values = row["inputs"] + row["outputs"]
        writer.writerow(["" if v is None else f"{v:.10g}" for v in values]
                        + [(row["error"] or {}).get("message", "")])
    return out.getvalue()


_SI = {"temperature": "K", "pressure": "Pa", "pressure_difference": "Pa", "molar_flow": "mol/s",
       "mass_flow": "kg/s", "power": "W", "temperature_difference": "K",
       "thermal_conductance": "W/K", "volume": "m3", "length": "m", "area": "m2", "time": "s"}


def _header(col: dict) -> str:
    unit = _SI.get(col["quantity"] or "")
    return f"{col['label']} [{unit}]" if unit else col["label"]
