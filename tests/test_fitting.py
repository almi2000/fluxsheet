"""Binary parameter regression from ThermoML data (no network: data are made up or synthetic).

Synthetic data are generated with the same phase model the fit uses (modified Raoult's law
with thermo's vapor pressures), so a correct fit recovers the generating parameters.
"""

import json

import numpy as np
import pytest
from thermo import NRTL, UNIQUAC

from fluxsheet import solve
from fluxsheet.cli import main
from fluxsheet.diagrams import binary_diagram
from fluxsheet.fitting import pipeline, store
from fluxsheet.fitting.consistency import area_test, point_test
from fluxsheet.fitting.models import Model, Pair, nrtl_ln_gammas, uniquac_ln_gammas
from fluxsheet.fitting.regression import FitOptions, fit_binary
from fluxsheet.fitting.thermoml import (
    ThermoMLCompound,
    ThermoMLError,
    VLEDataset,
    VLEPoint,
    parse_thermoml,
)
from fluxsheet.thermo_package import get_package

ATM = 101325.0
WATER, ETHANOL = "7732-18-5", "64-17-5"

# A ThermoML file shaped like the archive's: the T-x and y-x halves of one isobaric T-P-x-y
# set are separate blocks, plus a density set and a liquid-liquid set that must be ignored.
SNIPPET = """<?xml version="1.0" encoding="UTF-8"?>
<DataReport xmlns="http://www.iupac.org/namespaces/ThermoML">
  <Citation>
    <sAuthor>Doe, J.[Jane]</sAuthor><sAuthor>Roe, R.[Rick]</sAuthor>
    <sPubName>J. Chem. Eng. Data</sPubName><yrPubYr>2001</yrPubYr><sVol>46</sVol>
    <sPage>1-5</sPage><sTitle>Made-up VLE of ethanol + water</sTitle><sDOI>10.9999/test01</sDOI>
  </Citation>
  <Compound><RegNum><nOrgNum>1</nOrgNum></RegNum>
    <sStandardInChIKey>LFQSCWFLJHTTHZ-UHFFFAOYSA-N</sStandardInChIKey>
    <sCommonName>ethanol</sCommonName><sFormulaMolec>C2H6O</sFormulaMolec></Compound>
  <Compound><RegNum><nOrgNum>2</nOrgNum></RegNum>
    <sStandardInChIKey>XLYOFNOQVPJJNP-UHFFFAOYSA-N</sStandardInChIKey>
    <sCommonName>water</sCommonName><sFormulaMolec>H2O</sFormulaMolec></Compound>
  <PureOrMixtureData>
    <nPureOrMixtureDataNumber>1</nPureOrMixtureDataNumber>
    <Component><RegNum><nOrgNum>1</nOrgNum></RegNum></Component>
    <Component><RegNum><nOrgNum>2</nOrgNum></RegNum></Component>
    <Property><nPropNumber>1</nPropNumber><Property-MethodID><PropertyGroup>
      <VaporPBoilingTAzeotropTandP><ePropName>Boiling temperature at pressure P, K</ePropName>
      </VaporPBoilingTAzeotropTandP></PropertyGroup></Property-MethodID>
      <PropPhaseID><ePropPhase>Liquid</ePropPhase></PropPhaseID></Property>
    <PhaseID><ePhase>Liquid</ePhase></PhaseID><PhaseID><ePhase>Gas</ePhase></PhaseID>
    <Constraint><ConstraintID><ConstraintType><ePressure>Pressure, kPa</ePressure>
      </ConstraintType></ConstraintID><nConstraintValue>101.325</nConstraintValue></Constraint>
    <Variable><nVarNumber>1</nVarNumber><VariableID><VariableType>
      <eComponentComposition>Mole fraction</eComponentComposition></VariableType>
      <RegNum><nOrgNum>1</nOrgNum></RegNum></VariableID>
      <VarPhaseID><eVarPhase>Liquid</eVarPhase></VarPhaseID></Variable>
    {rows_T}
  </PureOrMixtureData>
  <PureOrMixtureData>
    <nPureOrMixtureDataNumber>2</nPureOrMixtureDataNumber>
    <Component><RegNum><nOrgNum>1</nOrgNum></RegNum></Component>
    <Component><RegNum><nOrgNum>2</nOrgNum></RegNum></Component>
    <Property><nPropNumber>1</nPropNumber><Property-MethodID><PropertyGroup>
      <CompositionAtPhaseEquilibrium><ePropName>Mole fraction</ePropName>
      </CompositionAtPhaseEquilibrium></PropertyGroup><RegNum><nOrgNum>1</nOrgNum></RegNum>
      </Property-MethodID><PropPhaseID><ePropPhase>Gas</ePropPhase></PropPhaseID></Property>
    <PhaseID><ePhase>Gas</ePhase></PhaseID><PhaseID><ePhase>Liquid</ePhase></PhaseID>
    <Constraint><ConstraintID><ConstraintType><ePressure>Pressure, kPa</ePressure>
      </ConstraintType></ConstraintID><nConstraintValue>101.325</nConstraintValue></Constraint>
    <Variable><nVarNumber>1</nVarNumber><VariableID><VariableType>
      <eComponentComposition>Mole fraction</eComponentComposition></VariableType>
      <RegNum><nOrgNum>1</nOrgNum></RegNum></VariableID>
      <VarPhaseID><eVarPhase>Liquid</eVarPhase></VarPhaseID></Variable>
    {rows_y}
  </PureOrMixtureData>
  <PureOrMixtureData>
    <nPureOrMixtureDataNumber>3</nPureOrMixtureDataNumber>
    <Component><RegNum><nOrgNum>1</nOrgNum></RegNum></Component>
    <Component><RegNum><nOrgNum>2</nOrgNum></RegNum></Component>
    <Property><nPropNumber>1</nPropNumber><Property-MethodID><PropertyGroup>
      <VolumetricProp><ePropName>Mass density, kg/m3</ePropName></VolumetricProp>
      </PropertyGroup></Property-MethodID></Property>
    <PhaseID><ePhase>Liquid</ePhase></PhaseID>
    <NumValues><PropertyValue><nPropNumber>1</nPropNumber><nPropValue>900</nPropValue>
    </PropertyValue></NumValues>
  </PureOrMixtureData>
  <PureOrMixtureData>
    <nPureOrMixtureDataNumber>4</nPureOrMixtureDataNumber>
    <Component><RegNum><nOrgNum>1</nOrgNum></RegNum></Component>
    <Component><RegNum><nOrgNum>2</nOrgNum></RegNum></Component>
    <PhaseID><ePhase>Liquid mixture 1</ePhase></PhaseID>
    <PhaseID><ePhase>Liquid mixture 2</ePhase></PhaseID>
  </PureOrMixtureData>
</DataReport>
"""

XS = (0.1, 0.3, 0.5, 0.7, 0.9)
TS = (360.15, 354.95, 352.85, 351.85, 351.35)
YS = (0.44, 0.58, 0.65, 0.73, 0.88)


def snippet() -> str:
    def rows(values):
        return "".join(
            f"<NumValues><VariableValue><nVarNumber>1</nVarNumber><nVarValue>{x}</nVarValue>"
            f"</VariableValue><PropertyValue><nPropNumber>1</nPropNumber>"
            f"<nPropValue>{v}</nPropValue></PropertyValue></NumValues>"
            for x, v in zip(XS, values))
    return SNIPPET.replace("{rows_T}", rows(TS)).replace("{rows_y}", rows(YS))


# --- ThermoML parser --------------------------------------------------------------------

def test_parser_joins_split_tables_and_skips_other_data():
    report = parse_thermoml(snippet())
    assert report.doi == "10.9999/test01"
    assert report.citation == "Doe and Roe, J. Chem. Eng. Data 46 (2001) 1-5"
    assert len(report.datasets) == 1
    ds = report.datasets[0]
    assert [c.name for c in ds.compounds] == ["ethanol", "water"]
    assert ds.kind == "isobaric" and ds.has_vapor
    assert [p.x1 for p in ds.points] == list(XS)
    assert [p.T for p in ds.points] == list(TS)
    assert [p.y1 for p in ds.points] == list(YS)
    assert all(p.P == pytest.approx(ATM) for p in ds.points)


def test_parser_refuses_entities_and_foreign_documents():
    bomb = ('<?xml version="1.0"?><!DOCTYPE r [<!ENTITY a "aaaa">]>'
            '<DataReport xmlns="http://www.iupac.org/namespaces/ThermoML">&a;</DataReport>')
    with pytest.raises(ThermoMLError):
        parse_thermoml(bomb)
    with pytest.raises(ThermoMLError):
        parse_thermoml("<html><body>not ThermoML</body></html>")
    with pytest.raises(ThermoMLError):
        parse_thermoml("<DataReport><unclosed></DataReport>")


def test_compounds_map_to_cas_and_data_are_oriented(tmp_path):
    path = tmp_path / "test01.xml"
    path.write_text(snippet(), encoding="utf-8")
    pair = Pair.from_cas(WATER, ETHANOL)
    found = pipeline.datasets_for_pair(pair, [path], skipped := [])
    assert not skipped and len(found) == 1
    ds, citation = found[0]
    assert "Doe" in citation
    assert [c.name for c in ds.compounds] == ["water", "ethanol"]  # swapped to the pair
    assert ds.points[0].x1 == pytest.approx(0.9) and ds.points[0].y1 == pytest.approx(0.56)


# --- models -----------------------------------------------------------------------------

def test_activity_models_match_thermo():
    pair = Pair.from_cas(WATER, ETHANOL)
    p, T, x = (0.5, -1.2, 300.0, 150.0), 340.0, 0.3
    ref = NRTL(T=T, xs=[x, 1 - x], tau_as=[[0, p[0]], [p[1], 0]], tau_bs=[[0, p[2]], [p[3], 0]],
               alpha_cs=[[0, 0.3], [0.3, 0]]).gammas()
    assert np.exp(nrtl_ln_gammas(p, 0.3, T, x)) == pytest.approx(ref, rel=1e-10)
    (r1, r2), (q1, q2) = pair.rq
    p = (0.5, -1.2, -300.0, -150.0)
    ref = UNIQUAC(T=T, xs=[x, 1 - x], rs=[r1, r2], qs=[q1, q2], tau_as=[[0, p[0]], [p[1], 0]],
                  tau_bs=[[0, p[2]], [p[3], 0]]).gammas()
    assert np.exp(uniquac_ln_gammas(p, pair.rq, T, x)) == pytest.approx(ref, rel=1e-10)


# --- regression on synthetic data -------------------------------------------------------

COMPOUNDS = (ThermoMLCompound(1, ("water",)), ThermoMLCompound(2, ("ethanol",)))
XGRID = np.linspace(0.04, 0.96, 13)


def synthetic(model: Model, p, temperatures=(), pressures=(ATM,)) -> list:
    sets = []
    for i, T in enumerate(temperatures):
        P, y = model.bubble_P(p, np.full_like(XGRID, T), XGRID)
        sets.append(VLEDataset("10.9999/synthetic", i, COMPOUNDS, "isothermal", tuple(
            VLEPoint(T, float(pp), float(x), float(yy)) for pp, x, yy in zip(P, XGRID, y))))
    for j, P in enumerate(pressures):
        T, y = model.bubble_T(p, np.full_like(XGRID, P), XGRID)
        sets.append(VLEDataset("10.9999/synthetic", 100 + j, COMPOUNDS, "isobaric", tuple(
            VLEPoint(float(t), P, float(x), float(yy)) for t, x, yy in zip(T, XGRID, y))))
    return [(ds, "synthetic") for ds in sets]


@pytest.mark.parametrize(("kind", "true"), [
    ("nrtl", (0.0, 0.0, 620.0, -30.0)),  # ChemSep-like water/ethanol
    ("uniquac", (0.0, 0.0, -60.0, -90.0)),
])
def test_fit_recovers_parameters_from_isobaric_data(kind, true):
    pair = Pair.from_cas(WATER, ETHANOL)
    data = synthetic(Model(kind, pair, 0.3), true)
    result = fit_binary(pair, data, FitOptions(model=kind, alpha=0.3))
    assert result.free == ("b12", "b21")  # 22 K of data: no temperature dependence
    assert result.params == pytest.approx(true, abs=1e-3)
    assert result.n_rejected == 0
    assert result.stats["AAD_T_K"] < 1e-6 and result.stats["AAD_y"] < 1e-7


def test_fit_recovers_temperature_dependence_from_isotherms():
    pair = Pair.from_cas(WATER, ETHANOL)
    true = (1.2, -0.8, 250.0, 260.0)
    data = synthetic(Model("nrtl", pair, 0.3), true, temperatures=(300.0, 330.0, 360.0),
                     pressures=())
    result = fit_binary(pair, data, FitOptions(alpha=0.3))
    assert result.free == ("a12", "a21", "b12", "b21")
    assert result.params == pytest.approx(true, rel=1e-4)
    assert result.stats["AAD_P_pct"] < 1e-5


def test_outlier_is_rejected_and_reported():
    pair = Pair.from_cas(WATER, ETHANOL)
    true = (0.0, 0.0, 620.0, -30.0)
    (ds, cite), = synthetic(Model("nrtl", pair, 0.3), true)
    points = list(ds.points)
    p = points[6]
    points[6] = VLEPoint(p.T + 3.0, p.P, p.x1, p.y1 + 0.05)  # a bad point
    bad = VLEDataset(ds.doi, ds.number, ds.compounds, ds.kind, tuple(points))
    result = fit_binary(pair, [(bad, cite)], FitOptions(alpha=0.3, require_consistency=False))
    assert result.n_rejected == 1 and result.datasets[0].n_rejected == 1
    assert result.params == pytest.approx(true, abs=0.01)


def test_consistency_tests_pass_model_data_and_catch_shifted_vapor():
    pair = Pair.from_cas(WATER, ETHANOL)
    (ds, _), = synthetic(Model("nrtl", pair, 0.3), (0.0, 0.0, 620.0, -30.0))
    assert point_test(pair, ds).passed
    assert area_test(pair, ds).passed
    shifted = VLEDataset(ds.doi, ds.number, ds.compounds, ds.kind, tuple(
        VLEPoint(p.T, p.P, p.x1, min(p.y1 + 0.04, 0.999)) for p in ds.points))
    assert not point_test(pair, shifted).passed
    report = fit_binary(pair, [(ds, "ok"), (shifted, "bad")], FitOptions(alpha=0.3)).datasets
    assert [r.used for r in report] == [True, False]
    assert "consistency" in report[1].reason


# --- parameter store and property packages ----------------------------------------------

def test_store_round_trip_in_either_order(tmp_path):
    entry = {"components": [{"name": "water", "cas": WATER}, {"name": "ethanol", "cas": ETHANOL}],
             "a12": 0.1, "a21": 0.2, "b12": 300.0, "b21": -40.0, "alpha": 0.3,
             "sources": [{"doi": "10.9999/x"}]}
    path = store.save_entry("nrtl", entry, tmp_path)
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["format_version"] == store.FORMAT_VERSION and len(data["pairs"]) == 1
    store.save_entry("nrtl", entry, tmp_path)  # replaces, does not duplicate
    assert len(json.loads(path.read_text(encoding="utf-8"))["pairs"]) == 1
    fwd = store.lookup("nrtl", WATER, ETHANOL, tmp_path)
    rev = store.lookup("nrtl", ETHANOL, WATER, tmp_path)
    assert (fwd.a_ij, fwd.b_ij, fwd.b_ji) == (0.1, 300.0, -40.0)
    assert (rev.a_ij, rev.b_ij, rev.b_ji) == (0.2, -40.0, 300.0)
    assert fwd.dois == ["10.9999/x"]


def test_seed_parameters_are_committed_with_provenance():
    for model in ("nrtl", "uniquac"):
        data = store.load(model)
        assert data["pairs"], model
        for e in data["pairs"]:
            assert e["sources"] and all(s["doi"].startswith("10.") for s in e["sources"])
            assert e["n_points"] >= 10 and e["T_range_K"][0] <= e["T_range_K"][1]
            assert "AAD_y" in e["fit"]


@pytest.mark.parametrize("method", ["nrtl", "uniquac"])
def test_package_prefers_thermoml_fit_and_reports_it(method):
    pkg = get_package(("water", "ethanol"), method)
    (src,) = pkg.parameter_sources
    assert src["source"] == "ThermoML fit" and src["references"]
    assert pkg.warnings == []
    chemsep = get_package(("water", "ethanol"), method, prefer_fitted=False)
    assert chemsep.parameter_sources == [{"pair": "water/ethanol", "source": "ChemSep"}]


def test_fitted_water_ethanol_reproduces_the_azeotrope():
    d = binary_diagram(["ethanol", "water"], "nrtl", P=ATM)
    assert d["parameter_sources"][0]["source"] == "ThermoML fit"
    az = d["azeotrope"]
    assert az["type"] == "homogeneous"
    assert az["x"] == pytest.approx(0.89, abs=0.02)  # measured 0.894
    assert az["T"] - 273.15 == pytest.approx(78.2, abs=0.15)  # measured 78.15 degC


def test_solve_result_lists_parameter_sources():
    r = solve({
        "components": ["water", "ethanol"], "property_method": "nrtl",
        "units": [{"id": "F", "type": "Feed", "params": {
            "T": "25 degC", "P": "1 atm", "flows": {"water": "1 kmol/h", "ethanol": "1 kmol/h"}}}],
        "streams": [{"id": "S", "from": "F.out"}],
    })
    assert r.ok, r.error
    assert r.parameter_sources[0]["source"] == "ThermoML fit"
    assert r.to_dict()["parameter_sources"] == r.parameter_sources


def test_cli_fit_binary_offline_from_cache(tmp_path, capsys):
    # The search result and file of the made-up article stand in for the archive.
    cache = tmp_path / "cache"
    keys = sorted(["LFQSCWFLJHTTHZ-UHFFFAOYSA-N", "XLYOFNOQVPJJNP-UHFFFAOYSA-N"])
    (cache / "search").mkdir(parents=True)
    (cache / "search" / f"{'_'.join(keys)}.json").write_text('[["10.9999/test01", 10]]')
    (cache / "10.9999").mkdir()
    (cache / "10.9999" / "test01.xml").write_text(snippet(), encoding="utf-8")
    code = main(["fit-binary", "water", "ethanol", "--offline", "--cache", str(cache),
                 "--no-consistency", "--write", "--output", str(tmp_path), "--quiet"])
    assert code == 0
    out = capsys.readouterr().out
    assert "1 VLE data sets in 1 ThermoML files" in out and "ChemSep" in out
    entry = store.load("nrtl", tmp_path)["pairs"][0]
    assert entry["sources"][0]["doi"] == "10.9999/test01"
    assert entry["n_points"] + entry["n_rejected_points"] == 5
