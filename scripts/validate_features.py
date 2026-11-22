"""Compare the equilibrium reactive column and the Gibbs reactor on fugacities with literature.

The reference values are small published tables, kept in the campaign manifest. Nothing is
downloaded. Exit status: 0 when every registered limit holds, 1 when a screening limit fails,
2 when the manifest is invalid or a case cannot be run.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from fluxsheet import solve

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "docs" / "validation-campaign-features.json"


def metric(quantity, unit, errors, thresholds):
    """Mean and maximum absolute error over (label, error) pairs against [mean, max] limits."""
    labels, values = zip(*errors, strict=True)
    values = np.abs(values)
    worst = int(np.argmax(values))
    mean, maximum = float(values.mean()), float(values[worst])
    return {"quantity": quantity, "unit": unit, "points": len(values),
            "mean_absolute_error": mean, "maximum_absolute_error": maximum,
            "maximum_error_point": labels[worst], "thresholds": list(thresholds),
            "passed": mean <= thresholds[0] and maximum <= thresholds[1]}


def _gibbs(case, T_C, P_atm, fugacities):
    model = {
        "components": ["nitrogen", "hydrogen", "ammonia"],
        "property_method": case["property_method"],
        "units": [{"id": "F", "type": "Feed",
                   "params": {"T": f"{T_C} degC", "P": f"{P_atm} atm", "flows": case["feed"]}},
                  {"id": "R", "type": "GibbsReactor",
                   "params": {"fugacities": fugacities, "T": f"{T_C} degC"}},
                  {"id": "P", "type": "Product"}],
        "streams": [{"id": "S1", "from": "F.out", "to": "R.in"},
                    {"id": "S2", "from": "R.out", "to": "P.in"}],
    }
    r = solve(model)
    if not r.ok:
        raise ValueError(f"{T_C} degC, {P_atm} atm: {r.error['message']}")
    return 100.0 * r.streams["S2"]["mole_fractions"]["ammonia"]


def run_ammonia(case):
    metrics, points = [], []
    for series in case["series"]:
        errors, ideal = [], []
        for T, row in case["measured"].items():
            for P in series["pressures_atm"]:
                measured = row[str(P)]
                real = _gibbs(case, T, P, "property_method")
                gas = _gibbs(case, T, P, "ideal_gas")
                at = f"{T} degC, {P} atm"
                errors.append((f"{at} (predicted {real:.2f}, measured {measured:.2f})",
                               real - measured))
                ideal.append((f"{at} (predicted {gas:.2f}, measured {measured:.2f})",
                              gas - measured))
                points.append({"T_C": float(T), "P_atm": P, "measured": measured,
                               "property_method": real, "ideal_gas": gas})
        metrics.append(metric(f"{case['quantity']}, {series['name']}", case["unit"], errors,
                              series["thresholds"]))
        # The ideal gas is reported beside it for contrast; it is not screened.
        reference = metric(f"{case['quantity']}, {series['name']}, ideal gas", case["unit"],
                           ideal, series["thresholds"])
        reference["screened"] = False
        metrics.append(reference)
    return metrics, points


def run_column(case):
    model = json.loads((ROOT / case["example"]).read_text(encoding="utf-8"))
    r = solve(model)
    if not r.ok:
        raise ValueError(r.error["message"])
    profiles = r.units[case["column"]]["profiles"]
    T = profiles["temperature"]["series"][0]["values"]
    x = {s["label"]: s["values"] for s in profiles["liquid_composition"]["series"]}
    limits = case["thresholds"]
    t_errors = [(f"stage {j} (predicted {T[int(j) - 1]:.2f} K, measured {m} K)",
                 T[int(j) - 1] - m) for j, m in case["measured_temperature_K"].items()]
    x_errors = [(f"stage {j}, {c} (predicted {x[c][int(j) - 1]:.4f}, measured {m})",
                 x[c][int(j) - 1] - m)
                for j, row in case["measured_liquid_mole_fraction"].items()
                for c, m in row.items()]
    metrics = [metric("stage temperature", "K", t_errors, limits["temperature"]),
               metric("liquid mole fraction", "fraction", x_errors, limits["composition"])]
    points = {"temperature_K": T, "liquid_mole_fraction": x}
    return metrics, points


RUNNERS = {"ammonia_larson": run_ammonia, "methyl_acetate_mekala": run_column}


def validate(manifest):
    if manifest.get("version") != 1 or not manifest.get("cases"):
        raise ValueError("unsupported or empty feature campaign")
    seen = set()
    for case in manifest["cases"]:
        if case.get("id") not in RUNNERS or case["id"] in seen:
            raise ValueError(f"unknown or duplicate case {case.get('id')}")
        if not case.get("references"):
            raise ValueError(f"case {case['id']} names no reference")
        seen.add(case["id"])


def run(manifest):
    results, errors = [], []
    for case in manifest["cases"]:
        try:
            metrics, points = RUNNERS[case["id"]](case)
        except (ValueError, KeyError, IndexError) as exc:
            errors.append({"case": case["id"], "error": str(exc)})
            continue
        screened = [m for m in metrics if m.get("screened", True)]
        results.append({"case": case["id"], "feature": case["feature"],
                        "references": case["references"],
                        "passed": all(m["passed"] for m in screened),
                        "metrics": metrics, "points": points})
    return results, errors


def markdown(results, errors):
    rows = ["# Feature reference campaign", "",
            "| Case | Quantity | Points | Mean | Maximum | Limits | Status | Maximum-error point |",
            "| --- | --- | ---: | ---: | ---: | --- | --- | --- |"]
    for result in results:
        for m in result["metrics"]:
            status = ("pass" if m["passed"] else "fail") if m.get("screened", True) else "—"
            rows.append(f"| {result['case']} | {m['quantity']} | {m['points']} | "
                        f"{m['mean_absolute_error']:.4g} {m['unit']} | "
                        f"{m['maximum_absolute_error']:.4g} {m['unit']} | "
                        f"{m['thresholds'][0]:g}/{m['thresholds'][1]:g} | {status} | "
                        f"{m['maximum_error_point']} |")
    if errors:
        rows.extend(["", "## Not run", ""])
        rows.extend(f"- {e['case']}: {e['error']}" for e in errors)
    return "\n".join(rows) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        validate(manifest)
    except (OSError, ValueError, TypeError, AttributeError) as exc:
        parser.exit(2, f"invalid campaign: {exc}\n")
    results, errors = run(manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = {"results": results, "errors": errors}
    (args.output_dir / "features.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (args.output_dir / "features.md").write_text(markdown(results, errors), encoding="utf-8")
    failed = sum(not r["passed"] for r in results)
    print(f"{len(results)} cases; {failed} failed screening; {len(errors)} not run")
    return 2 if errors else (1 if failed else 0)


if __name__ == "__main__":
    sys.exit(main())
