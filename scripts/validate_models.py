"""Run the registered experimental benchmarks with a verified offline cache."""

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import re
import subprocess
import sys
from pathlib import Path

from fluxsheet.fitting.archive import ArchiveError, fetch
from fluxsheet.fitting.thermoml import MAX_BYTES
from fluxsheet.thermo_package import get_package


def _sibling(name):
    """A module next to this one, whether this file runs as a script or as ``scripts.*``."""
    return importlib.import_module(f"{__package__}.{name}" if __package__ else name)


_registry = _sibling("equipment_registry")
equipment_markdown, qualify = _registry.markdown, _registry.qualify
_btx = _sibling("validate_btx_reference")
_ethanol_vle = _sibling("validate_ethanol_water_reference")
_ethanol_holdout = _sibling("validate_ethanol_water_holdout")
_butanol_2010 = _sibling("validate_butanol_water_reference")
_chen_screening = json.loads(
    _sibling("validate_chen_lle").MANIFEST.read_text(encoding="utf-8"))["screening"]

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = {
    "btx": ("validate_btx_reference.py", {"ideal"}, 2, 11),
    "ethanol_vle": ("validate_ethanol_water_reference.py", {"unifac"}, 2, 16),
    "ethanol_holdout": ("validate_ethanol_water_holdout.py", {"nrtl", "uniquac"}, 1, 12),
    "butanol_2010": ("validate_butanol_water_reference.py", {"unifac", "nrtl", "uniquac"}, 1, 4),
    "butanol_2007": ("validate_chen_lle.py", {"unifac", "nrtl", "uniquac"}, 1, 4),
}
# The limits the validators themselves apply, one [mean, max] pair per reported metric. A
# campaign manifest must register exactly these.
THRESHOLDS = {
    "btx": [[_btx.MEAN_LIMIT_PERCENT, _btx.MAX_LIMIT_PERCENT]] * 2,
    "ethanol_vle": [[_ethanol_vle.MEAN_T_K, _ethanol_vle.MAX_T_K],
                    [_ethanol_vle.MEAN_Y, _ethanol_vle.MAX_Y]],
    "ethanol_holdout": [[_ethanol_holdout.MEAN_LIMIT_PERCENT,
                         _ethanol_holdout.MAX_LIMIT_PERCENT]],
    "butanol_2010": [[_butanol_2010.MEAN_LIMIT, _butanol_2010.MAX_LIMIT]],
    "butanol_2007": [[_chen_screening["mean_absolute"], _chen_screening["max_absolute"]]],
}
# Quantity and unit of each summary line in a legacy validator's text output.
QUANTITIES = {
    "btx": [("bubble pressure", "%")] * 2,
    "ethanol_vle": [("bubble temperature", "K"), ("vapor mole fraction", "fraction")],
    "ethanol_holdout": [("bubble pressure", "%")],
    "butanol_2010": [("butanol mass fraction", "fraction")],
}
METRIC = re.compile(r"AAD=([0-9.]+).*?max=([0-9.]+).*?at\s+(.+?)"
                    r"(?:;\s*limits\b.*|\s+(?:PASS|FAIL))?$")
SERIES = re.compile(r"^(T=[0-9.]+ K) n=")
PHASE_ROW = re.compile(r"^T=([0-9.]+) K (aqueous|organic): .*absolute error=([0-9.]+)$")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def versions():
    result = {"python": sys.version.split()[0]}
    for name in ("fluxsheet", "thermo", "chemicals", "scipy"):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = None
    return result


def parse_legacy(source_id, output, expected_count, points):
    lines = [line.strip() for line in output.splitlines() if "AAD=" in line and "max=" in line]
    if len(lines) != expected_count:
        raise ValueError(f"{source_id}: unexpected validator output")
    metrics = []
    for line, (quantity, unit) in zip(lines, QUANTITIES[source_id]):
        match = METRIC.search(line)
        if not match:
            raise ValueError(f"{source_id}: incomplete validator metric")
        series = SERIES.match(line)  # one summary line per isotherm
        if series:
            quantity = f"{quantity}, {series[1]}"
        metrics.append({"quantity": quantity, "unit": unit,
                        "points": points, "mean_absolute_error": float(match[1]),
                        "maximum_absolute_error": float(match[2]),
                        "maximum_error_point": match[3]})
    phase_errors = {}
    for line in output.splitlines():
        match = PHASE_ROW.match(line)
        if match:
            phase_errors.setdefault(match[2], []).append(float(match[3]))
    phases = {name: {"points": len(values), "mean_absolute_error": sum(values) / len(values),
                     "maximum_absolute_error": max(values)}
              for name, values in phase_errors.items()}
    if source_id == "butanol_2010" and set(phases) != {"aqueous", "organic"}:
        raise ValueError("missing liquid-phase summaries")
    return metrics, phases


def run_case(source, method, path):
    script, _, count, points = SCRIPTS[source["id"]]
    command = [sys.executable, str(ROOT / "scripts" / script), str(path)]
    if source["id"] in ("ethanol_holdout", "butanol_2010", "butanol_2007"):
        command.extend(["--method", method])
    result = subprocess.run(command, text=True, capture_output=True, timeout=180, check=False)
    if result.returncode not in (0, 1):
        raise ValueError(f"validator failed: {result.stderr.strip()}")
    if source["id"] == "butanol_2007":
        value = json.loads(result.stdout)
        if value["doi"] != source["doi"] or value["sha256"] != source["sha256"]:
            raise ValueError("validator returned an unexpected source")
        worst = value["maximum_error_point"]
        metrics = [{"quantity": "butanol mass fraction", "unit": "fraction",
                    "points": len(value["measurements"]),
                    "mean_absolute_error": value["mean_absolute_error"],
                    "maximum_absolute_error": worst["absolute_error"],
                    "maximum_error_point": worst}]
        phases = {}
        for phase in ("aqueous", "organic"):
            entries = [v["absolute_error"] for v in value["measurements"] if v["phase"] == phase]
            if len(entries) != 2:
                raise ValueError("unexpected phase measurement count")
            phases[phase] = {"points": len(entries), "mean_absolute_error": sum(entries) / len(entries),
                             "maximum_absolute_error": max(entries)}
        if value["passed"] != (result.returncode == 0):
            raise ValueError("inconsistent validator exit status")
        screening_passed = all(metric["mean_absolute_error"] <= limits[0] and
                               metric["maximum_absolute_error"] <= limits[1]
                               for metric, limits in zip(metrics, source["thresholds"]))
        if screening_passed != (result.returncode == 0):
            raise ValueError("validator result disagrees with registered thresholds")
    else:
        if "AAD=" not in result.stdout:  # the validator stopped before its summary
            reason = result.stderr.strip().splitlines()[-1:] or ["no output"]
            raise ValueError(f"validator failed: {reason[0]}")
        metrics, phases = parse_legacy(source["id"], result.stdout, count, points)
        # The text summary is rounded, so the limits are not applied to it again: the
        # validator's exit status decides, and THRESHOLDS are its own limits.
        if ("FAIL" in result.stdout) != (result.returncode == 1):
            raise ValueError("inconsistent validator exit status")
    package = get_package(tuple(source["components"]), method)
    return {"source": source["id"], "doi": source["doi"], "sha256": source["sha256"],
            "method": method, "passed": result.returncode == 0, "metrics": metrics,
            "phases": phases, "independence": "unconfirmed", "limits": source["limits"],
            "uncertainty": source["uncertainty"], "parameter_sources": package.parameter_sources,
            "point_details": value["measurements"] if source["id"] == "butanol_2007" else None,
            "validator_sha256": digest(ROOT / "scripts" / script)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--offline", action="store_true", help="Use cached files only (default)")
    mode.add_argument("--download", action="store_true", help="Fetch missing files from NIST")
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        if manifest["version"] not in (1, 2) or not isinstance(manifest["sources"], list):
            raise ValueError("unsupported or empty campaign manifest")
        if manifest["version"] == 1 and "equipment_sources" in manifest:
            raise ValueError("equipment_sources requires manifest version 2")
        if not manifest["sources"] and not manifest.get("equipment_sources"):
            raise ValueError("empty campaign manifest")
        if args.cache_dir.resolve().is_relative_to(ROOT):
            raise ValueError("raw-data cache must be outside the repository")
        seen = set()
        for source in manifest["sources"]:
            name = source["id"]
            if name not in SCRIPTS or name in seen:
                raise ValueError(f"unknown or duplicate source {name}")
            seen.add(name)
            if not all(isinstance(source.get(key), str) for key in ("doi", "sha256")):
                raise ValueError(f"missing doi or sha256 for {name}")
            if not source["methods"] or not set(source["methods"]).issubset(SCRIPTS[name][1]):
                raise ValueError(f"unsupported method for {name}")
            if source["thresholds"] != THRESHOLDS[name]:
                raise ValueError(f"registered thresholds differ from validator limits for {name}")
        equipment = qualify(manifest.get("equipment_sources", []), seen)
    except (OSError, KeyError, ValueError, TypeError) as exc:
        parser.exit(2, f"invalid campaign: {exc}\n")

    results = []
    errors = [{"source": r["id"], "error": r["reason"], "kind": "equipment_data_gap"}
              for r in equipment]
    for source in manifest["sources"]:
        try:
            path = fetch(source["doi"], args.cache_dir, offline=not args.download, pause=0)
            if path.stat().st_size > MAX_BYTES or digest(path) != source["sha256"]:
                raise ValueError("reference size or SHA-256 mismatch")
            for method in source["methods"]:
                try:
                    results.append(run_case(source, method, path))
                except (ValueError, KeyError, subprocess.TimeoutExpired) as exc:
                    errors.append({"source": source["id"], "method": method, "error": str(exc)})
        except (ArchiveError, OSError, ValueError) as exc:
            errors.append({"source": source["id"], "error": str(exc)})
    parameter_hashes = {model: digest(ROOT / "src" / "fluxsheet" / "data" /
                                      "binary_parameters" / f"{model}.json")
                        for model in ("nrtl", "uniquac")}
    report = {"manifest_sha256": digest(args.manifest), "versions": versions(),
              "parameter_file_sha256": parameter_hashes, "results": results, "errors": errors,
              "equipment_sources": equipment}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "validation.json").write_text(json.dumps(report, indent=2, allow_nan=False)
                                                     + "\n", encoding="utf-8")
    rows = ["# Experimental validation campaign", "", ("| Source | Method | Status | Quantity | "
            "Points | Mean | Maximum | Maximum-error point |"), ("| --- | --- | --- | --- | ---: | "
            "---: | ---: | --- |")]
    for result in results:
        for metric in result["metrics"]:
            point = metric["maximum_error_point"]
            if isinstance(point, dict):
                point = f"{point['temperature_K']} K, {point['phase']}"
            rows.append(f"| {result['source']} | {result['method']} | "
                        f"{'pass' if result['passed'] else 'fail'} | {metric['quantity']} | "
                        f"{metric['points']} | {metric['mean_absolute_error']:.6g} {metric['unit']} | "
                        f"{metric['maximum_absolute_error']:.6g} {metric['unit']} | {point} |")
    rows.extend(["", "Source details, phase summaries and errors are in validation.json.", ""])
    if equipment:
        rows.extend(equipment_markdown(equipment))
    if errors:
        rows.extend(["## Unavailable or invalid", ""])
        rows.extend(f"- {error['source']}: {error['error']}" for error in errors)
        rows.append("")
    (args.output_dir / "validation.md").write_text("\n".join(rows), encoding="utf-8")
    print(f"{len(results)} comparisons; {sum(not r['passed'] for r in results)} failed; "
          f"{len(errors)} unavailable or invalid")
    return 2 if errors else (1 if any(not r["passed"] for r in results) else 0)


if __name__ == "__main__":
    raise SystemExit(main())
