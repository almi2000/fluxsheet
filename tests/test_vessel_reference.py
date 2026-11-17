"""Offline CLI tests using a closed-form stirred vessel enthalpy balance."""

import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest
from scipy.optimize import brentq
from scripts import validate_vessel_reference as vessel

from fluxsheet.thermo_package import get_package

ROOT = Path(__file__).resolve().parents[1]


def _file(path):
    return {"path": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


@pytest.fixture
def campaign(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    pressure, t_initial, feed_rate = 101325.0, 293.15, 1.0
    pkg = get_package(("water",), "ideal")

    def props(temperature):
        return pkg.phase_props("liquid", temperature, pressure, [1.0])

    area = math.pi * 0.1**2 / 4
    inventory = area * 0.1 / props(t_initial).V
    model = {"components": ["water"], "property_method": "ideal",
             "units": [{"id": "F", "type": "Feed", "params": {"P": pressure}},
                       {"id": "T", "type": "MixingTank", "params": {"P": pressure}}],
             "streams": [{"id": "feed", "from": "F.out", "to": "T.in"},
                         {"id": "out", "from": "T.liquid"}],
             "dynamics": {"t_end": 600, "units": {
                 "T": {"diameter": 0.1, "height": 0.2, "initial_level": 0.5}}}}
    flowsheet = cache / "tank.json"
    flowsheet.write_text(json.dumps(model), encoding="utf-8")
    runs = []
    for name in ("heater-step", "feed-temperature-step"):
        data = cache / f"{name}.csv"
        with data.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["time", "inlet", "flow", "duty", "temperature", "level"])
            for t in (0, 5, 10, 15, 23, 41, 69, 100, 160, 240, 340, 450, 520, 600):
                inlet = t_initial + (20 if t >= 15 and name == "feed-temperature-step" else 0)
                duty = 1000 if t >= 15 and name == "heater-step" else 0
                # N is constant because the two flows are equal. Solve N dH/dt = F(Hin-H)+Q
                # with fixed post-step inputs, then invert the liquid property equation.
                delta = (props(inlet).H - props(t_initial).H + duty / feed_rate) * (
                    1 - math.exp(-max(t-15, 0)*feed_rate/inventory))
                h = props(t_initial).H + delta
                temperature = brentq(lambda T, enthalpy=h: props(T).H - enthalpy, 282, 335)
                level = inventory * props(temperature).V / area
                writer.writerow([t / 60, inlet - 273.15, feed_rate * 3.6, duty / 1000,
                                 temperature - 273.15, 100 * level])
        def input_column(name, unit):
            return {"column": name, "unit": unit, "interpolation": "step"}
        runs.append({
            "id": name, "flowsheet": _file(flowsheet), "measurements": _file(data),
            "time": {"column": "time", "unit": "min"},
            "feeds": {"F": {"temperature": input_column("inlet", "degC"),
                            "flow": input_column("flow", "kmol/h")}},
            "liquid_flow": input_column("flow", "kmol/h"), "duty": input_column("duty", "kW"),
            "outputs": {"temperature": {"column": "temperature", "unit": "degC", "uncertainty": 0.1},
                        "level": {"column": "level", "unit": "cm", "uncertainty": None}},
            "initial": {"temperature_K": t_initial, "temperature_tolerance_K": 1e-4,
                        "inventory_quantity": "level", "inventory_si": 0.1,
                        "inventory_tolerance_si": 1e-5, "evidence": "synthetic initial condition"},
            "baseline_s": [0, 10.1], "comparison_s": [15, 600], "end_start_s": 500})
    manifest = {"version": 2, "sources": [], "equipment_sources": [{
        "id": "vessel-fixture", "family": "dynamics", "admission": "qualified",
        "url": "urn:fluxsheet:synthetic-vessel", "evidence": "closed-form enthalpy balance",
        "data_kind": "synthetic", "independence": "separate",
        "parameter_provenance": "fixed water properties and synthetic apparatus geometry",
        "runs": runs}]}
    return cache, manifest, runs


def _cli(tmp_path, cache, manifest):
    spec = tmp_path / "manifest.json"
    spec.write_text(json.dumps(manifest), encoding="utf-8")
    output = tmp_path / "output"
    process = subprocess.run([sys.executable, str(ROOT / "scripts" / "validate_models.py"),
        "--manifest", str(spec), "--offline", "--cache-dir", str(cache),
        "--output-dir", str(output)], capture_output=True, text=True, timeout=90, check=False)
    assert not process.stderr, process.stderr
    return process.returncode, json.loads((output / "validation.json").read_text()), output


def test_two_registered_steps_through_cli(tmp_path, campaign):
    cache, manifest, _ = campaign
    code, report, output = _cli(tmp_path, cache, manifest)
    assert code == 0, report["errors"]
    source = report["equipment_sources"][0]
    assert source["experimental_validation"] == "open"
    assert source["screening"] == "pass"
    assert len(source["run_results"]) == 2
    for run in source["run_results"]:
        temp, level = run["metrics"]
        assert temp["uncertainty"] == pytest.approx(0.1)
        assert level["uncertainty"] is None
        assert temp["time_s"][1] == pytest.approx(23)
        assert temp["rmse"] < 0.01 * temp["step_amplitude"]
        assert level["rmse"] < 0.01 * level["step_amplitude"]
    assert "| temperature | pass |" in (output / "validation.md").read_text()


@pytest.mark.parametrize("scenario,code", [("biased_measurement", 1), ("hash", 2), ("initial", 2),
                                          ("weak_signal", 2), ("missing_sample", 2)])
def test_rejected_reference_and_retained_comparisons(tmp_path, campaign, scenario, code):
    cache, manifest, runs = campaign
    run = runs[1]
    if scenario == "hash":
        run["measurements"]["sha256"] = "0" * 64
    elif scenario == "initial":
        run["initial"]["temperature_K"] = 295
    else:
        path = cache / run["measurements"]["path"]
        with path.open(newline="", encoding="utf-8") as stream:
            header, *rows = list(csv.reader(stream))
        if scenario == "weak_signal":
            for row in rows:
                row[4] = "20"
        elif scenario == "missing_sample":
            rows[4][4] = ""
        else:
            for row in rows[3:]:
                row[4] = str(float(row[4]) + 5)
        with path.open("w", newline="", encoding="utf-8") as stream:
            csv.writer(stream).writerows([header, *rows])
        run["measurements"] = _file(path)
    actual, report, _ = _cli(tmp_path, cache, manifest)
    assert actual == code, report["errors"]
    source = report["equipment_sources"][0]
    assert source["run_results"][0]["passed"]
    if code == 2:
        assert source["screening"] == "not_checkable"
        assert report["errors"][0]["run"] == run["id"]
    else:
        assert source["screening"] == "fail"
        assert not source["run_results"][1]["metrics"][0]["passed"]


def test_refinement_gate(campaign, monkeypatch):
    cache, _, runs = campaign
    real_simulate = vessel._simulate

    def perturbed(model, tank, run, points, tolerance):
        times, values, warnings = real_simulate(model, tank, run, points, tolerance)
        if points == 501:
            values["T"] = [value + 1 for value in values["T"]]
        return times, values, warnings

    monkeypatch.setattr(vessel, "_simulate", perturbed)
    with pytest.raises(ValueError, match="numerical refinement"):
        vessel.run_case(runs[0], cache)
