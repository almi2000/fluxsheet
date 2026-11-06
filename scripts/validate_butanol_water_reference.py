"""Compare UNIFAC LLE with the zero-salt rows of a ThermoML study."""

import argparse
import hashlib
from pathlib import Path

from defusedxml import ElementTree

from fluxsheet.fitting.thermoml import MAX_BYTES, parse_thermoml
from fluxsheet.thermo_package import get_package

DOI = "10.1016/j.fluid.2010.09.013"
NS = {"t": "http://www.iupac.org/namespaces/ThermoML"}
MEAN_LIMIT, MAX_LIMIT = 0.05, 0.10  # absolute 1-butanol mass-fraction error


def values(row, tag, number, value):
    return {e.find(f"t:{number}", NS).text: float(e.find(f"t:{value}", NS).text)
            for e in row.findall(f"t:{tag}", NS)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xml", type=Path)
    parser.add_argument("--method", choices=("unifac", "nrtl", "uniquac"), default="unifac")
    args = parser.parse_args()
    if args.xml.stat().st_size > MAX_BYTES:
        raise ValueError("ThermoML file is too large")
    raw = args.xml.read_bytes()
    report = parse_thermoml(raw)
    if report.doi != DOI:
        raise ValueError(f"expected DOI {DOI}, got {report.doi}")
    for index, expected in ((1, "butan-1-ol"), (2, "sodium nitrate"), (3, "water")):
        if expected not in report.compounds[index].names:
            raise ValueError(f"unexpected component {index}")
    root = ElementTree.fromstring(raw, forbid_dtd=True)
    tables = [d for d in root.findall("t:PureOrMixtureData", NS)
              if d.findtext("t:nPureOrMixtureDataNumber", namespaces=NS) == "1"]
    if len(tables) != 1:
        raise ValueError("expected a single table 1")
    table = tables[0]
    properties = table.findall("t:Property", NS)
    for prop, expected in zip(properties, (("1", "1"), ("1", "2"), ("3", "2"))):
        comp = prop.findtext(".//t:RegNum/t:nOrgNum", namespaces=NS)
        phase = prop.findtext("t:PropPhaseID/t:ePropPhase", namespaces=NS)
        name = prop.findtext(".//t:ePropName", namespaces=NS)
        if (comp, phase, name) != (expected[0], f"Liquid mixture {expected[1]}", "Mass fraction"):
            raise ValueError("unexpected phase or component in property table")
    if len(properties) != 3:
        raise ValueError("expected three properties")
    rows = []
    for row in table.findall("t:NumValues", NS):
        variables = values(row, "VariableValue", "nVarNumber", "nVarValue")
        props = values(row, "PropertyValue", "nPropNumber", "nPropValue")
        # Water and butanol mass fractions sum to one in both phases.
        if (abs(props["1"] + variables["2"] - 1) <= 1e-8
                and abs(props["2"] + props["3"] - 1) <= 1e-8):
            rows.append((variables["1"], props["1"], props["2"]))
    if len(rows) != 2 or sorted(t for t, _, _ in rows) != [298.15, 308.15]:
        raise ValueError("expected exactly two zero-salt measurements")
    package = get_package(("water", "1-butanol"), args.method)
    errors = []
    print(f"DOI {DOI} | SHA-256 {hashlib.sha256(raw).hexdigest()}")
    print(f"{args.method} | P=101000 Pa | 1-butanol mass fraction")
    for T, aqueous, organic in rows:
        state = package.flash([0.8, 0.2], T=T, P=101000.0)
        if not state.liquid_split:
            raise ValueError(f"no two-liquid split at {T:g} K")
        fractions = []
        for phase in (state.x_light, state.x_heavy):
            masses = [x * mw for x, mw in zip(phase, package.mws)]
            fractions.append(masses[1] / sum(masses))
        calculated = sorted(fractions)
        for label, predicted, measured in zip(("aqueous", "organic"), calculated,
                                               (aqueous, organic)):
            error = abs(predicted - measured)
            errors.append((error, T, label, predicted, measured))
            print(f"T={T:g} K {label}: predicted={predicted:.6f} measured={measured:.6f} "
                  f"absolute error={error:.6f}")
    mean = sum(e[0] for e in errors) / len(errors)
    maximum, temp, phase, _, _ = max(errors)
    passed = mean <= MEAN_LIMIT and maximum <= MAX_LIMIT
    print(f"AAD={mean:.6f}, max={maximum:.6f} at {temp:g} K ({phase}); "
          f"limits {MEAN_LIMIT:g}/{MAX_LIMIT:g}: {'PASS' if passed else 'FAIL'}")
    return int(not passed)


if __name__ == "__main__":
    raise SystemExit(main())
