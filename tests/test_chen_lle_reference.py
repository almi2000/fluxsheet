"""Synthetic LLE data exercise source validation without external measurements."""

import hashlib
from types import SimpleNamespace

import pytest
from scripts import validate_chen_lle as reference


def synthetic_xml():
    tables = []
    for number, phase, pressure, x in [(3, 2, 101, 0.1), (4, 1, 101.325, 0.6)]:
        tables.append(f'''<PureOrMixtureData><nPureOrMixtureDataNumber>{number}</nPureOrMixtureDataNumber>
          <Component><RegNum><nOrgNum>4</nOrgNum></RegNum></Component>
          <Component><RegNum><nOrgNum>5</nOrgNum></RegNum></Component>
          <Property><Property-MethodID><PropertyGroup><CompositionAtPhaseEquilibrium>
            <ePropName>Mole fraction</ePropName></CompositionAtPhaseEquilibrium></PropertyGroup>
            <RegNum><nOrgNum>5</nOrgNum></RegNum></Property-MethodID>
            <PropPhaseID><ePropPhase>Liquid mixture {phase}</ePropPhase></PropPhaseID></Property>
          <Constraint><ConstraintID><ConstraintType><ePressure>Pressure, kPa</ePressure>
          </ConstraintType></ConstraintID><nConstraintValue>{pressure}</nConstraintValue></Constraint>
          <NumValues><VariableValue><nVarValue>300</nVarValue></VariableValue>
          <PropertyValue><nPropValue>{x}</nPropValue><CombinedUncertainty>
          <nCombExpandUncertValue>0.01</nCombExpandUncertValue></CombinedUncertainty>
          </PropertyValue></NumValues></PureOrMixtureData>''')
    return ('<ThermoML xmlns="http://www.iupac.org/namespaces/ThermoML">'
            + ''.join(tables) + '</ThermoML>').encode()


def test_extracts_both_phases_without_erasing_pressure_difference(monkeypatch):
    raw = synthetic_xml()
    manifest = {"sha256": hashlib.sha256(raw).hexdigest(), "doi": "example",
                "datasets": [3, 4], "reported_pressures_Pa": [101000, 101325],
                "temperatures_K": [300]}
    monkeypatch.setattr(reference, "parse_thermoml", lambda _: SimpleNamespace(
        doi="example", compounds={4: SimpleNamespace(names=("water",)),
                                  5: SimpleNamespace(names=("butan-1-ol",))}))
    data = reference.extract(raw, manifest)
    assert data == [(3, 300, 101000, 0.1, 0.01), (4, 300, 101325, 0.6, 0.01)]
    with pytest.raises(ValueError, match="SHA-256"):
        reference.extract(raw + b" ", manifest)


def test_mass_fraction_conversion_at_endpoints():
    assert reference.mass_fraction(0, 18, 74) == 0
    assert reference.mass_fraction(1, 18, 74) == 1
    assert reference.mass_fraction(0.5, 18, 74) == pytest.approx(74 / 92)
