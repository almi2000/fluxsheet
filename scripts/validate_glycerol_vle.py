"""Refit local NRTL parameters from binary ThermoML data and check a ternary holdout."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from thermo import NRTL

from fluxsheet.flowsheet import BinaryParameterModel, FlowsheetModel
from fluxsheet.thermo_package import get_package, package_options

COELHO = "10.1016/j.jct.2011.06.016"
ZAOUI = "10.1016/j.jct.2013.09.046"
ZHANG = "10.1021/acs.jced.5b00116"
HASHES = {
    COELHO: "958a077c28caf5ee1b0c4c8267fda330e258a6b0a563f561117945f98eeae875",
    ZAOUI: "b34691b5195ff97d1507095a16aca72153e9ebc2601ff59b8cbd4e48da222ceb",
    ZHANG: "7c1df83928f553b787f4db345d23dae44e51943552dce2fbea0ebc245e102aae",
}
COMPONENTS = ("ethanol", "water", "glycerol")
ALPHA = 0.2  # Coelho et al., Table 6
# Screening thresholds applied to the entire held-out ternary series.
LIMITS = {"T_K": (2.0, 5.0), "y_ethanol": (0.03, 0.08)}


def load_tables(path: Path, doi: str) -> list[dict]:
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != HASHES[doi]:
        raise ValueError(f"unexpected SHA-256 for {doi}: {path}")
    return json.loads(raw)["PureOrMixtureData"]


def rows(table: dict) -> list[tuple[float, float, float]]:
    return [(values[1], values[2], row["PropertyValue"][0]["nPropValue"])
            for row in table["NumValues"]
            for values in [{v["nVarNumber"]: v["nVarValue"] for v in row["VariableValue"]}]]


def fit_pair(component: str, points: list[tuple[float, float, float]]) -> BinaryParameterModel:
    pkg = get_package((component, "glycerol"), "ideal")
    points = [(T, P, x) for T, P, x in points if 0 < x < 1]
    rows = [(T, P, [x, 1 - x], np.array(pkg.psats(T, [x, 1 - x]))) for T, P, x in points]

    def residual(coeffs):
        a12, a21, b12, b21 = coeffs
        model = NRTL(T=350.0, xs=[0.5, 0.5], tau_as=[[0.0, a12], [a21, 0.0]],
                     tau_bs=[[0.0, b12], [b21, 0.0]], alpha_cs=[[0.0, ALPHA], [ALPHA, 0.0]])
        return [np.log(np.dot(np.multiply(x, model.to_T_xs(T, x).gammas()), psat) / P)
                for T, P, x, psat in rows]

    fit = least_squares(residual, [0.0] * 4, x_scale=[1.0, 1.0, 350.0, 350.0])
    if not fit.success:
        raise ValueError(f"NRTL fit failed for {component}/glycerol: {fit.message}")
    temperatures, pressures, _ = zip(*points)
    return BinaryParameterModel(
        method="nrtl", components=(component, "glycerol"), convention="thermo", alpha=ALPHA,
        source=f"Fluxsheet regression of binary NIST TRC ThermoML data (DOIs {COELHO}; {ZAOUI})",
        a=tuple(round(value, 9) for value in fit.x[:2]),
        b=tuple(round(value, 9) for value in fit.x[2:]),
        T_range_K=(min(temperatures), max(temperatures)),
        P_range_Pa=(min(pressures), max(pressures)),
    )


def holdout(model: FlowsheetModel, tables: list[dict]) -> dict:
    pkg = get_package(COMPONENTS, "nrtl", **package_options(model))
    mwe, mww, mwg = pkg.mws
    errors = {"T_K": [], "y_ethanol": []}
    for row, (temperature, vapor) in enumerate(zip(rows(tables[2]), rows(tables[1]), strict=True), 1):
        xe_solvent, wg, measured_T = temperature
        if temperature[:2] != vapor[:2]:
            raise ValueError(f"unmatched ternary temperature/vapor data at row {row}")
        liquid_moles = (1 - wg) / (xe_solvent * mwe + (1 - xe_solvent) * mww)
        amounts = np.array([xe_solvent * liquid_moles, (1 - xe_solvent) * liquid_moles, wg / mwg])
        sat = pkg.saturation_point((amounts / amounts.sum()).tolist(), "bubble", P=100000.0)
        errors["T_K"].append((abs(sat.T - measured_T), row))
        errors["y_ethanol"].append((abs(sat.y[0] - vapor[2]), row))
    return {key: {"points": len(values), "mean": float(np.mean([v for v, _ in values])),
                  "max": max(values)[0], "max_row": max(values)[1]}
            for key, values in errors.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coelho", type=Path, help=f"ThermoML JSON for {COELHO}")
    parser.add_argument("zaoui", type=Path, help=f"ThermoML JSON for {ZAOUI}")
    parser.add_argument("zhang", type=Path, help=f"ThermoML JSON for {ZHANG}")
    args = parser.parse_args()
    coelho, zaoui, zhang = (load_tables(path, doi) for path, doi in
                            ((args.coelho, COELHO), (args.zaoui, ZAOUI), (args.zhang, ZHANG)))
    pairs = [fit_pair(component, [(T, 1000 * P, x) for x, P, T in rows(coelho[i])]
                      + [(T, 1000 * P, x) for T, x, P in rows(zaoui[j])])
             for component, i, j in (("ethanol", 6, 9), ("water", 5, 10))]
    model = FlowsheetModel(components=list(COMPONENTS), property_method="nrtl", units=[],
                           property_data={"binary_parameters": pairs})
    prediction = holdout(model, zhang)
    print(json.dumps({"parameters": [p.model_dump(mode="json") for p in pairs],
                      "holdout": prediction,
                      "baseline": holdout(model.model_copy(update={"property_data": None}), zhang),
                      "limits_mean_max": LIMITS}, indent=2))
    return int(any(prediction[k]["mean"] > mean or prediction[k]["max"] > maximum
                   for k, (mean, maximum) in LIMITS.items()))


if __name__ == "__main__":
    raise SystemExit(main())
