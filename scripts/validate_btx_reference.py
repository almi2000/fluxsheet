"""Compare the ideal benzene/toluene bubble pressure with one ThermoML article.

Pass the downloaded XML path. Raw ThermoML records are kept outside this repository.
"""

import argparse
import hashlib
from pathlib import Path

from fluxsheet.fitting.thermoml import parse_thermoml
from fluxsheet.thermo_package import get_package

DOI = "10.1016/j.fluid.2009.03.012"
MEAN_LIMIT_PERCENT = 5.0
MAX_LIMIT_PERCENT = 10.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xml", type=Path)
    args = parser.parse_args()
    raw = args.xml.read_bytes()
    report = parse_thermoml(raw)
    if report.doi != DOI:
        raise ValueError(f"expected DOI {DOI}, got {report.doi}")
    series = [d for d in report.datasets
              if [c.name for c in d.compounds] == ["benzene", "toluene"]
              and d.kind == "isothermal"]
    if len(series) != 2 or any(len(d.points) != 11 for d in series):
        raise ValueError("expected two benzene/toluene series of 11 points")
    package = get_package(("benzene", "toluene"), "ideal")
    failed = False
    print(f"DOI {DOI} | SHA-256 {hashlib.sha256(raw).hexdigest()}")
    print("Method ideal | mean limit 5%, maximum limit 10%")
    for dataset in series:
        errors = []
        for point in dataset.points:
            predicted = package.flash([point.x1, 1 - point.x1], T=point.T, VF=0.0).P
            errors.append((abs(predicted - point.P) / point.P * 100, point.x1,
                           predicted, point.P))
        mean = sum(e[0] for e in errors) / len(errors)
        maximum, x, predicted, measured = max(errors)
        passed = mean <= MEAN_LIMIT_PERCENT and maximum <= MAX_LIMIT_PERCENT
        failed |= not passed
        print(f"T={dataset.points[0].T:g} K n={len(errors)} AAD={mean:.4f}% "
              f"max={maximum:.4f}% at x_benzene={x:g} "
              f"(predicted={predicted:.1f} Pa, measured={measured:.1f} Pa) "
              f"{'PASS' if passed else 'FAIL'}")
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
