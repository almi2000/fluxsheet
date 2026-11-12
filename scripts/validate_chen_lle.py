"""Evaluate a pre-registered binary LLE reference without redistributing raw data."""

import argparse
import hashlib
import json
import math
from pathlib import Path

from defusedxml import ElementTree

from fluxsheet.fitting.thermoml import MAX_BYTES, parse_thermoml
from fluxsheet.thermo_package import get_package

NS = {"t": "http://www.iupac.org/namespaces/ThermoML"}
MANIFEST = Path(__file__).resolve().parents[1] / "docs" / "chen2007-lle-manifest.json"


def mass_fraction(x, mw_water, mw_butanol):
    return x * mw_butanol / (x * mw_butanol + (1 - x) * mw_water)


def extract(raw, manifest):
    if len(raw) > MAX_BYTES or hashlib.sha256(raw).hexdigest() != manifest["sha256"]:
        raise ValueError("Reference file exceeds size limit or SHA-256 does not match manifest")
    report = parse_thermoml(raw)
    if report.doi != manifest["doi"]:
        raise ValueError("Unexpected DOI")
    if "water" not in report.compounds[4].names or "butan-1-ol" not in report.compounds[5].names:
        raise ValueError("Unexpected compound identities")
    root = ElementTree.fromstring(raw, forbid_dtd=True)
    measurements = []
    for number, expected_phase, expected_pressure in zip(
            manifest["datasets"], ("Liquid mixture 2", "Liquid mixture 1"),
            manifest["reported_pressures_Pa"]):
        tables = [d for d in root.findall("t:PureOrMixtureData", NS)
                  if d.findtext("t:nPureOrMixtureDataNumber", namespaces=NS) == str(number)]
        if len(tables) != 1:
            raise ValueError(f"Missing or duplicate dataset {number}")
        table = tables[0]
        compounds = {c.text for c in table.findall("t:Component/t:RegNum/t:nOrgNum", NS)}
        prop = table.find("t:Property", NS)
        if (compounds != {"4", "5"} or prop is None
                or prop.findtext(".//t:ePropName", namespaces=NS) != "Mole fraction"
                or prop.findtext(".//t:RegNum/t:nOrgNum", namespaces=NS) != "5"
                or prop.findtext("t:PropPhaseID/t:ePropPhase", namespaces=NS) != expected_phase):
            raise ValueError(f"Unexpected component, phase or property in dataset {number}")
        pressure_type = table.findtext("t:Constraint/t:ConstraintID/t:ConstraintType/t:ePressure",
                                       namespaces=NS)
        pressure = float(table.findtext("t:Constraint/t:nConstraintValue", namespaces=NS)) * 1000
        if pressure_type != "Pressure, kPa" or pressure != expected_pressure:
            raise ValueError(f"Unexpected pressure in dataset {number}")
        rows = table.findall("t:NumValues", NS)
        if len(rows) != len(manifest["temperatures_K"]):
            raise ValueError(f"Unexpected row count in dataset {number}")
        temperatures = []
        for row in rows:
            temperature = float(row.findtext("t:VariableValue/t:nVarValue", namespaces=NS))
            temperatures.append(temperature)
            x = float(row.findtext("t:PropertyValue/t:nPropValue", namespaces=NS))
            uncertainty = row.findtext(
                "t:PropertyValue/t:CombinedUncertainty/t:nCombExpandUncertValue", namespaces=NS)
            uncertainty = float(uncertainty) if uncertainty is not None else None
            if (temperature not in manifest["temperatures_K"] or not math.isfinite(x)
                    or not 0 <= x <= 1 or (uncertainty is not None and
                                            (not math.isfinite(uncertainty) or uncertainty < 0))):
                raise ValueError(f"Invalid measurement in dataset {number}")
            measurements.append((number, temperature, pressure, x, uncertainty))
        if sorted(temperatures) != sorted(manifest["temperatures_K"]):
            raise ValueError(f"Missing or duplicate temperatures in dataset {number}")
    return measurements


def evaluate(raw, manifest, method):
    data = extract(raw, manifest)
    pkg = get_package(("water", "1-butanol"), method)
    water, butanol = pkg.mws
    results = []
    for number, temperature, pressure, x, uncertainty in data:
        state = pkg.flash([0.8, 0.2], T=temperature, P=pressure)
        if not state.liquid_split:
            raise ValueError(f"No two-liquid split at {temperature} K, {pressure} Pa")
        predicted = sorted(mass_fraction(phase[1], water, butanol)
                           for phase in (state.x_light, state.x_heavy))
        measured = mass_fraction(x, water, butanol)
        derivative = water * butanol / (x * butanol + (1 - x) * water) ** 2
        deviation = abs(predicted[0 if number == 3 else 1] - measured)
        results.append({"dataset": number, "temperature_K": temperature,
                        "pressure_Pa": pressure,
                        "phase": "aqueous" if number == 3 else "organic",
                        "predicted_mass_fraction": predicted[0 if number == 3 else 1],
                        "measured_mass_fraction": measured,
                        "expanded_uncertainty_mass_fraction": (
                            uncertainty * derivative if uncertainty is not None else None),
                        "absolute_error": deviation})
    mean = sum(r["absolute_error"] for r in results) / len(results)
    worst = max(results, key=lambda r: r["absolute_error"])
    limits = manifest["screening"]
    return {"doi": manifest["doi"], "sha256": manifest["sha256"], "method": method,
            "parameter_sources": pkg.parameter_sources, "independence": manifest["independence"],
            "measurements": results, "mean_absolute_error": mean,
            "maximum_error_point": worst,
            "passed": mean <= limits["mean_absolute"] and
            worst["absolute_error"] <= limits["max_absolute"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xml", type=Path)
    parser.add_argument("--method", choices=("unifac", "nrtl", "uniquac"), default="unifac")
    args = parser.parse_args()
    try:
        if args.xml.stat().st_size > MAX_BYTES:
            raise ValueError("Reference file exceeds size limit")
        result = evaluate(args.xml.read_bytes(), json.loads(MANIFEST.read_text()), args.method)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"invalid reference data: {exc}\n")
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
