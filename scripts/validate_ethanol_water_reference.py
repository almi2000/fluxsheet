"""Compare UNIFAC ethanol/water VLE predictions with one ThermoML article."""

import argparse
import hashlib
from pathlib import Path

from fluxsheet.fitting.thermoml import parse_thermoml
from fluxsheet.thermo_package import get_package

DOI = "10.1016/j.jct.2013.08.020"
# Screening limits set before the comparison, not experimental uncertainty.
MEAN_T_K, MAX_T_K = 2.0, 5.0
MEAN_Y, MAX_Y = 0.03, 0.08


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xml", type=Path)
    args = parser.parse_args()
    raw = args.xml.read_bytes()
    report = parse_thermoml(raw)
    if report.doi != DOI:
        raise ValueError(f"expected DOI {DOI}, got {report.doi}")
    datasets = [d for d in report.datasets if d.number == 9 and d.kind == "isobaric"
                and [c.name for c in d.compounds] == ["ethanol", "water"]]
    if len(datasets) != 1 or len(datasets[0].points) != 16 or not datasets[0].has_vapor:
        raise ValueError("expected the 16-point ethanol/water T-x-y series")
    dataset = datasets[0]
    package = get_package(("ethanol", "water"), "unifac")
    temperature_errors, vapor_errors = [], []
    for point in dataset.points:
        state = package.saturation_point([point.x1, 1 - point.x1], "bubble", P=point.P)
        temperature_errors.append((abs(state.T - point.T), point.x1, state.T, point.T))
        vapor_errors.append((abs(state.y[0] - point.y1), point.x1, state.y[0], point.y1))
    mean_T = sum(e[0] for e in temperature_errors) / len(temperature_errors)
    mean_y = sum(e[0] for e in vapor_errors) / len(vapor_errors)
    max_T, x_T, pred_T, obs_T = max(temperature_errors)
    max_y, x_y, pred_y, obs_y = max(vapor_errors)
    passed = (mean_T <= MEAN_T_K and max_T <= MAX_T_K
              and mean_y <= MEAN_Y and max_y <= MAX_Y)
    print(f"DOI {DOI} | SHA-256 {hashlib.sha256(raw).hexdigest()}")
    print(f"UNIFAC Dortmund, P={dataset.points[0].P:.0f} Pa, n={len(dataset.points)}")
    print(f"Temperature: AAD={mean_T:.4f} K, max={max_T:.4f} K at x={x_T:g} "
          f"(predicted={pred_T:.4f} K, measured={obs_T:.4f} K)")
    print(f"Vapor composition: AAD={mean_y:.5f}, max={max_y:.5f} at x={x_y:g} "
          f"(predicted={pred_y:.5f}, measured={obs_y:.5f})")
    print("Screening limits: mean/max 2/5 K and 0.03/0.08 mole fraction; "
          + ("PASS" if passed else "FAIL"))
    return int(not passed)


if __name__ == "__main__":
    raise SystemExit(main())
