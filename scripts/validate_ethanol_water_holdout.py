"""Compare fitted activity models with a withheld ThermoML pressure series."""

import argparse
import hashlib
from pathlib import Path

from fluxsheet.fitting.store import lookup
from fluxsheet.fitting.thermoml import parse_thermoml
from fluxsheet.thermo_package import get_package

DOI = "10.1016/j.fluid.2012.12.014"
MEAN_LIMIT_PERCENT = 5.0
MAX_LIMIT_PERCENT = 10.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xml", type=Path)
    parser.add_argument("--method", choices=("nrtl", "uniquac"), default="nrtl")
    args = parser.parse_args()
    raw = args.xml.read_bytes()
    report = parse_thermoml(raw)
    if report.doi != DOI:
        raise ValueError(f"expected DOI {DOI}, got {report.doi}")
    data = [d for d in report.datasets if d.number == 2 and d.kind == "isothermal"
            and [c.name for c in d.compounds] == ["water", "ethanol"]
            and abs(d.points[0].T - 363.3) < 0.01]
    if len(data) != 1 or len(data[0].points) != 12 or data[0].has_vapor:
        raise ValueError("expected the 12-point isothermal water/ethanol data set")
    fitted = lookup(args.method, "7732-18-5", "64-17-5")
    if fitted is None or DOI in fitted.dois:
        raise ValueError("the reference DOI is part of the current fit or fit metadata is missing")
    excluded = fitted.record.get("rejected_datasets", [])
    if not any(row.get("doi") == DOI and row.get("conditions") == "isothermal 363.30 K"
               and not row.get("used", True) for row in excluded):
        raise ValueError("the reference series is not recorded as excluded from the fit")
    package = get_package(("water", "ethanol"), args.method)
    expected_source = "ThermoML fit"
    if not any(p["source"] == expected_source for p in package.parameter_sources):
        raise ValueError("the expected fitted parameter set is unavailable")
    errors = []
    for point in data[0].points:
        predicted = package.saturation_point(
            [point.x1, 1 - point.x1], "bubble", T=point.T).P
        errors.append((abs(predicted - point.P) / point.P * 100, point.x1,
                       predicted, point.P))
    mean = sum(e[0] for e in errors) / len(errors)
    maximum, x, predicted, measured = max(errors)
    passed = mean <= MEAN_LIMIT_PERCENT and maximum <= MAX_LIMIT_PERCENT
    print(f"DOI {DOI} | SHA-256 {hashlib.sha256(raw).hexdigest()}")
    print(f"{args.method} ThermoML fit, T=363.3 K, n={len(errors)}")
    print(f"Bubble pressure AAD={mean:.4f}%, max={maximum:.4f}% at x_water={x:g} "
          f"(predicted={predicted:.1f} Pa, measured={measured:.1f} Pa)")
    print("Screening limits: mean/max 5%/10%; " + ("PASS" if passed else "FAIL"))
    return int(not passed)


if __name__ == "__main__":
    raise SystemExit(main())
