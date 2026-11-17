"""Compare a fixed-pressure, single-liquid stirred vessel with registered step data."""

import csv
import hashlib
import io
import math
from itertools import pairwise
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, FiniteFloat

from fluxsheet import FlowsheetModel
from fluxsheet.dynamics import simulate
from fluxsheet.fitting.thermoml import MAX_BYTES
from fluxsheet.quantities import to_si

if __package__:
    from .equipment_metrics import dynamic_step, interpolate
else:
    from equipment_metrics import dynamic_step, interpolate


class Record(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)


class CacheFile(Record):
    path: str = Field(min_length=1)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")

    def read(self, cache):
        relative = Path(self.path)
        cache = cache.resolve()
        path = (cache / relative).resolve()
        if relative.anchor or not path.is_relative_to(cache):
            raise ValueError(f"cache path must be relative and stay inside the cache: {self.path}")
        with path.open("rb") as stream:
            data = stream.read(MAX_BYTES + 1)
        if len(data) > MAX_BYTES or hashlib.sha256(data).hexdigest() != self.sha256:
            raise ValueError(f"reference size or SHA-256 mismatch: {self.path}")
        return data


class Column(Record):
    column: str = Field(min_length=1)
    unit: str = Field(min_length=1)


class Input(Column):
    interpolation: Literal["step", "linear"]


class Measurement(Column):
    uncertainty: FiniteFloat | None = Field(ge=0)


class Feed(Record):
    temperature: Input
    flow: Input


class InitialState(Record):
    temperature_K: FiniteFloat = Field(gt=0)
    temperature_tolerance_K: FiniteFloat = Field(ge=0)
    inventory_quantity: Literal["level", "holdup", "liquid_volume"]
    inventory_si: FiniteFloat = Field(gt=0)
    inventory_tolerance_si: FiniteFloat = Field(ge=0)
    evidence: str = Field(min_length=1)


class Run(Record):
    id: str = Field(min_length=1)
    flowsheet: CacheFile
    measurements: CacheFile
    time: Column
    feeds: dict[str, Feed] = Field(min_length=1)
    liquid_flow: Input
    duty: Input
    outputs: dict[Literal["temperature", "level"], Measurement] = Field(min_length=1)
    initial: InitialState
    baseline_s: tuple[FiniteFloat, FiniteFloat]
    comparison_s: tuple[FiniteFloat, FiniteFloat]
    end_start_s: FiniteFloat


def _read_columns(data, columns):
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig")), strict=True)
    headers = reader.fieldnames or []
    required = {col.column for col, _ in columns}
    if not required.issubset(headers) or len(headers) != len(set(headers)):
        raise ValueError("CSV columns are missing or duplicate")
    values = [[] for _ in columns]
    for line, row in enumerate(reader, 2):
        if None in row or None in row.values():
            raise ValueError(f"CSV row {line} has an inconsistent number of fields")
        for series, (column, quantity) in zip(values, columns):
            try:
                number = float(row[column.column])
                value = to_si(f"{number} {column.unit}", quantity)
                if not math.isfinite(value):
                    raise ValueError("nonfinite value")
            except ValueError as exc:
                raise ValueError(f"CSV row {line}, column {column.column}: {exc}") from exc
            series.append(value)
    return values


def _configure(data, run):
    model = FlowsheetModel.model_validate_json(data)
    if model.components != ["water"] or model.custom_compounds:
        raise ValueError("the vessel adapter requires the built-in pure water component")
    vessels = [u for u in model.units if u.type == "MixingTank"]
    if len(vessels) != 1 or any(u.type not in {"Feed", "MixingTank", "Product"} for u in model.units):
        raise ValueError("expected one MixingTank, direct feeds and a liquid outlet")
    tank = vessels[0]
    feed_ids = {unit.id for unit in model.units if unit.type == "Feed"}
    if feed_ids != run.feeds.keys():
        raise ValueError("input mappings must specify every feed")
    inlet_sources = {f"{uid}.out" for uid in feed_ids}
    outlet = [s for s in model.streams if s.source == f"{tank.id}.liquid"]
    products = {f"{u.id}.in" for u in model.units if u.type == "Product"}
    if (len(outlet) != 1 or (outlet[0].target is not None and outlet[0].target not in products)
        or len(model.streams) != len(feed_ids) + 1
        or {s.source for s in model.streams if s.target == f"{tank.id}.in"} != inlet_sources):
        raise ValueError("expected direct feeds and one liquid outlet")
    if tank.params.get("P") is None or any(tank.params.get(k) is not None for k in ("T", "area", "dP")):
        raise ValueError("the tank needs explicit pressure and measured net duty as its heat input")
    dyn = model.dynamics
    if not dyn or set(dyn.units) != {tank.id} or dyn.events or dyn.controllers:
        raise ValueError("register one dynamic tank with no preconfigured events or controllers")
    settings = dyn.units[tank.id]
    if (not {"diameter", "height", "initial_level"}.issubset(settings)
        or settings.get("liquid_outflow", "fixed") != "fixed"):
        raise ValueError("explicit geometry, initial level and imposed liquid outflow are required")
    return model.model_dump(by_alias=True, exclude_none=True), tank.id


def _events(times, values, target, interpolation):
    for i, (left, right) in enumerate(pairwise(values)):
        if left == right:
            continue
        if interpolation == "step":
            yield {**target, "time": times[i + 1], "value": right}
        else:
            yield {**target, "time": times[i], "ramp": times[i + 1] - times[i], "value": right}


def _prepare(model, tank, run, data):
    outputs = {"temperature": ("T", "temperature", "K"), "level": ("level", "length", "m")}
    channels = [(run.duty, "power", tank, "duty"),
                (run.liquid_flow, "molar_flow", tank, "liquid_flow")]
    for uid, feed in run.feeds.items():
        channels.extend([(feed.temperature, "temperature", uid, "T"),
                         (feed.flow, "molar_flow", uid, "flows")])
    columns = [(col, quantity) for col, quantity, _, _ in channels]
    columns.extend((col, outputs[name][1]) for name, col in run.outputs.items())
    times, *data = _read_columns(data, [(run.time, "time"), *columns])
    if len(times) < 2 or times[0] != 0 or any(a >= b for a, b in pairwise(times)):
        raise ValueError("CSV times must start at zero and strictly increase")
    b0, b1 = run.baseline_s
    c0, c1 = run.comparison_s
    if not 0 <= b0 < b1 < c0 <= run.end_start_s < c1 <= times[-1]:
        raise ValueError("invalid baseline, comparison or end window")
    baseline = [i for i, t in enumerate(times) if b0 <= t <= b1]
    compare = [i for i, t in enumerate(times) if c0 <= t <= c1]
    if len(baseline) < 3 or sum(times[i] >= run.end_start_s for i in compare) < 2:
        raise ValueError("need three baseline samples and two samples in the end window")
    dyn = model["dynamics"]
    units = {u["id"]: u["params"] for u in model["units"]}
    for (column, quantity, uid, parameter), values in zip(channels, data):
        if (quantity == "temperature" or parameter == "flows") and min(values) <= 0:
            raise ValueError(f"{column.column} must remain positive")
        if parameter == "liquid_flow" and min(values) < 0:
            raise ValueError(f"{column.column} must be nonnegative")
        target = {"unit": uid, "param": parameter}
        if parameter == "flows":
            target["component"] = "water"
        dyn["events"].extend(event for event in _events(times, values, target, column.interpolation)
                             if event["time"] <= c1)
        if uid == tank:
            dyn["units"][tank][parameter] = values[0]
            if parameter == "duty":
                units[tank]["duty"] = values[0]
        elif parameter == "flows":
            units[uid]["flows"] = {"water": values[0]}
        else:
            units[uid].pop("VF", None)
            units[uid]["T"] = values[0]
    if not any(b1 <= event["time"] <= c0 and (event["time"] > b1 or "ramp" in event)
               for event in dyn["events"]):
        raise ValueError("no input change between the baseline and start of comparison")
    dyn["t_end"] = c1
    return times, baseline, compare, outputs, data[len(channels):]


def _simulate(model, tank, run, points, tolerance):
    model["dynamics"].update(points=points, rtol=tolerance)
    result = simulate(model)
    if not result["ok"]:
        error = result["error"]
        raise ValueError(f"vessel simulation failed: {error.get('message', error)}")
    series = {s["key"]: s["values"] for s in result["series"]}
    values = {key.removeprefix(f"unit.{tank}."): val for key, val in series.items()
              if key.startswith(f"unit.{tank}.")}
    if any(not math.isfinite(value) for samples in values.values() for value in samples):
        raise ValueError("simulation produced nonfinite vessel results")
    initial = run.initial
    for field, expected, tol in [
        ("T", initial.temperature_K, initial.temperature_tolerance_K),
        (initial.inventory_quantity, initial.inventory_si, initial.inventory_tolerance_si)
    ]:
        if abs(values[field][0] - expected) > tol:
            raise ValueError(f"initial {field} differs from reference: {values[field][0]:g} "
                             f"versus {expected:g} (tolerance {tol:g})")
    if any(v > 0 for v in values["vapor_flow"]):
        raise ValueError("vessel entered the vapor region")
    if min(values["level"]) < 0.01 or max(values["level_fraction"]) > 1:
        raise ValueError("vessel level left the range for an imposed liquid outflow")
    return result["time"], values, result["warnings"]


def run_case(record, cache):
    run = Run.model_validate(record)
    model, tank = _configure(run.flowsheet.read(cache), run)
    times, baseline, compare, outputs, references = _prepare(
        model, tank, run, run.measurements.read(cache))
    coarse_t, coarse, _ = _simulate(model, tank, run, 501, 1e-6)
    fine_t, fine, warnings = _simulate(model, tank, run, 1001, 1e-7)
    targets = [times[i] for i in compare]
    end_start = next(i for i, t in enumerate(targets) if t >= run.end_start_s)
    metrics = []
    for (name, column), reference in zip(run.outputs.items(), references):
        key, quantity, unit = outputs[name]
        prediction = interpolate(fine_t, fine[key], targets)
        previous = interpolate(coarse_t, coarse[key], targets)
        metric = dynamic_step([times[i] for i in baseline], [reference[i] for i in baseline],
                              [reference[i] for i in compare], prediction, end_start)
        refinement = max(abs(p - q) for p, q in zip(prediction, previous))
        if refinement > 0.01 * metric["step_amplitude"]:
            raise ValueError(f"{name}: numerical refinement exceeds 1% of step amplitude")
        worst = max(range(len(targets)), key=lambda i: abs(metric["residuals"][i]))
        uncertainty = column.uncertainty
        if uncertainty is not None:
            difference = "temperature_difference" if quantity == "temperature" else quantity
            uncertainty = to_si(f"{uncertainty} {column.unit}", difference)
        metrics.append(dict(metric, quantity=name, unit=unit, time_s=targets,
                            uncertainty=uncertainty, maximum_error_time_s=targets[worst],
                            maximum_absolute_error=abs(metric["residuals"][worst]),
                            numerical_difference=refinement))
    return {"id": run.id, "passed": all(metric["passed"] for metric in metrics),
            "metrics": metrics, "warnings": warnings, "inputs": run.model_dump(mode="json"),
            "property_method": model["property_method"],
            "initial_prediction": {key: fine[key][0] for key in ("T", run.initial.inventory_quantity)}}
