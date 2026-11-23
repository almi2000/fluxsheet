"""DWSIM import (importers/dwsim.py): sample files in tests/data/dwsim, see the README there."""

import hashlib
import io
import json
import os
import zipfile
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from fluxsheet import limits
from fluxsheet.cli import main
from fluxsheet.flowsheet import FlowsheetModel
from fluxsheet.importers import DwsimImportError, dwsim, import_dwsim
from fluxsheet.server import app
from fluxsheet.solver import solve

DATA = Path(__file__).parent / "data" / "dwsim"


def _units(model):
    return {u["id"]: u for u in model["units"]}


def _streams(model):
    return {s["id"]: s for s in model["streams"]}


def _solve(model):
    result = solve(FlowsheetModel.model_validate(model))
    assert result.ok, result.error
    return result


def _zipped(xml: bytes, name="simulation.xml") -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(name, xml)
    return buf.getvalue()


# --- the sample flowsheets -----------------------------------------------------------------


def test_heater_flash_imports_and_solves():
    model, warnings = import_dwsim(DATA / "heater_flash.dwxml")
    assert model["format"] == "fluxsheet/1"
    assert model["name"] == "heater_flash"
    assert model["components"] == ["water", "ethanol"]
    assert model["property_method"] == "nrtl"
    units = _units(model)
    assert {u: units[u]["type"] for u in units} == {
        "FEED": "Feed", "E-101": "Heater", "V-101": "Flash", "VAPOR": "Product",
        "LIQUID": "Product"}
    assert units["FEED"]["params"] == {"T": "25 degC", "P": "1.2 bar",
                                       "flows": {"water": "80 kmol/h", "ethanol": "20 kmol/h"}}
    assert units["E-101"]["params"] == {"T": "92 degC", "dP": "0.1 bar"}
    assert _streams(model)["S2"] == {"id": "S2", "from": "E-101.out", "to": "V-101.in"}
    assert _streams(model)["VAPOR"]["from"] == "V-101.vapor"
    # DWSIM's layout, moved to the origin.
    assert units["FEED"]["ui"] == {"x": 40, "y": 120}
    assert units["E-101"]["ui"] == {"x": 200, "y": 120}
    # Unsupported parts are listed, not dropped silently.
    assert any("energy streams left out (Q-101)" in w for w in warnings)
    assert any("annotations" in w for w in warnings)

    result = _solve(model)
    vapor = result.streams["VAPOR"]
    assert vapor["T"] == pytest.approx(365.15)
    assert vapor["P"] == pytest.approx(1.1e5)
    assert vapor["mole_fractions"]["ethanol"] > 0.3  # the flash enriches ethanol
    feed, liquid = result.streams["FEED"], result.streams["LIQUID"]
    assert vapor["molar_flow"] + liquid["molar_flow"] == pytest.approx(feed["molar_flow"])


def test_gas_recycle_loop_closes_without_the_recycle_block():
    model, warnings = import_dwsim(DATA / "gas_recycle.dwxml")
    assert model["property_method"] == "peng-robinson"
    assert model["components"] == ["methane", "ethane", "propane", "butane"]
    units, streams = _units(model), _streams(model)
    assert units["SPL-100"]["params"] == {"split": 0.6}
    assert units["K-100"]["params"] == {"P": "30 bar", "efficiency": 0.72}
    assert units["E-100"]["params"] == {"T": "-40 degC", "dP": "0.5 bar"}
    assert "REC-1" not in units and "ADJ-1" not in units
    # The recycle block's outlet stream R2 is merged into S6, which now closes the loop.
    assert "R2" not in streams
    assert streams["S6"] == {"id": "S6", "from": "K-100.out", "to": "MIX-100.in"}
    assert streams["S5"]["from"] == "SPL-100.out1"
    assert streams["FUELGAS"]["from"] == "SPL-100.out2"
    assert any("recycle block REC-1" in w for w in warnings)
    assert any("ADJ-1 left out" in w for w in warnings)

    result = _solve(model)
    assert result.tear_streams
    feed = result.streams["FEED"]["mass_flow"]
    out = result.streams["FUELGAS"]["mass_flow"] + result.streams["NGL"]["mass_flow"]
    assert out == pytest.approx(feed, rel=1e-4)
    assert result.streams["S4"]["T"] == pytest.approx(233.15)


def test_unsupported_unit_becomes_feed_and_product():
    model, warnings = import_dwsim(DATA / "shortcut_column.dwxml")
    assert model["property_method"] == "ideal"
    units, streams = _units(model), _streams(model)
    assert "PIPE-1" not in units
    assert units["B"]["type"] == "Product"  # the pipe's inlet
    assert units["B2"]["type"] == "Feed"  # the pipe's outlet, at its stored state
    assert units["B2"]["params"]["P"] == "0.95 bar"
    assert streams["B2"]["to"] == "V-200.in"
    assert units["P-100"]["params"] == {"dP": "1 bar", "efficiency": 0.75}
    col = units["T-100"]["params"]
    assert col["light_key"] == "benzene" and col["reflux_ratio"] == 2.5
    assert col["lk_recovery"] == pytest.approx(0.9798, abs=1e-3)
    assert col["hk_recovery"] == pytest.approx(0.0099, abs=1e-3)
    assert any("PIPE-1 left out" in w and "pipe" in w for w in warnings)

    result = _solve(model)
    assert result.streams["D"]["mole_fractions"]["toluene"] == pytest.approx(0.01, abs=1e-3)
    assert result.streams["B"]["mole_fractions"]["benzene"] == pytest.approx(0.02, abs=1e-3)
    assert result.streams["B3"]["P"] == pytest.approx(8e4)


def test_zipped_file_with_exchanger_reactor_and_separator(tmp_path):
    xml = (DATA / "hydrogenation.dwxml").read_bytes()
    path = tmp_path / "hydrogenation.dwxmz"
    path.write_bytes(_zipped(xml, "hydrogenation.xml"))
    model, warnings = import_dwsim(path)
    assert model == import_dwsim(xml, name="hydrogenation")[0]  # bytes work the same
    units, streams = _units(model), _streams(model)
    # The hot side is on connector 1 here: assigned by the stored inlet temperatures.
    hx = units["E-1"]
    assert hx["params"] == {"min_approach": "20 K", "dP_hot": "0.3 bar", "dP_cold": "0.3 bar"}
    assert streams["FEED"]["to"] == "E-1.cold_in"
    assert streams["S4"] == {"id": "S4", "from": "R-1.out", "to": "E-1.hot_in"}
    assert streams["S5"]["from"] == "E-1.hot_out"
    assert model["reactions"] == [{"id": "Hydrogenation", "stoichiometry": {
        "benzene": -1.0, "hydrogen": -3.0, "cyclohexane": 1.0}}]
    assert units["R-1"]["params"] == {"reactions": ["Hydrogenation"], "conversions": ["0.99"],
                                      "keys": ["benzene"], "T": "240 degC"}
    assert units["X-1"]["params"] == {"fractions": {"cyclohexane": 0.99, "benzene": 0.5}}
    assert streams["RESIDUE"]["from"] == "X-1.bottom"
    assert any("X-1: no separation spec for hydrogen" in w for w in warnings)

    result = _solve(model)
    assert result.streams["S3"]["T"] == pytest.approx(473.15)
    assert result.streams["S4"]["T"] == pytest.approx(513.15)
    product = result.streams["PRODUCT"]["mole_fractions"]
    assert product["cyclohexane"] > 0.95


# --- mapping details, on small files built here --------------------------------------------


def _stream(name, T=300.0, P=2e5, flow=10.0, x=(("Water", 1.0),), vf=0.0):
    comps = "".join(f"<Compound><Name>{n}</Name><MoleFraction>{v}</MoleFraction></Compound>"
                    for n, v in x)
    return (f"<SimulationObject><Type>DWSIM.Thermodynamics.Streams.MaterialStream</Type>"
            f"<Name>{name}</Name><Phases><Phase><ID>0</ID><Compounds>{comps}</Compounds>"
            f"<Properties><temperature>{T}</temperature><pressure>{P}</pressure>"
            f"<molarflow>{flow}</molarflow></Properties></Phase>"
            f"<Phase><ID>2</ID><Properties><molarfraction>{vf}</molarfraction></Properties>"
            f"</Phase></Phases></SimulationObject>")


def _graphic(name, otype, ins=(), outs=(), x=0, y=0):
    """ins: (index, from object, its outlet index); outs: (index, to object, inlet index)."""
    def conns(items, kind):
        attr = "From" if kind == "In" else "To"
        n = max([i for i, *_ in items], default=-1) + 1
        by = {i: (o, j) for i, o, j in items}
        out = []
        for i in range(n):
            if i in by:
                o, j = by[i]
                out.append(f'<Connector IsAttached="true" Attached{attr}ObjID="{o}" '
                           f'Attached{attr}ConnIndex="{j}"/>')
            else:
                out.append('<Connector IsAttached="false"/>')
        return "".join(out)
    return (f"<GraphicObject><Name>{name}</Name><Tag>{name}</Tag><ObjectType>{otype}"
            f"</ObjectType><X>{x}</X><Y>{y}</Y><InputConnectors>{conns(ins, 'In')}"
            f"</InputConnectors><OutputConnectors>{conns(outs, 'Out')}</OutputConnectors>"
            f"</GraphicObject>")


def _file(unit_cls, otype, props, outlets=1, states=None, package="RaoultPropertyPackage",
          compounds=(("Water", "7732-18-5"),), inlets=1, extra=""):
    """FEED(s) -> unit -> OUT0..OUTn, with the unit's spec elements ``props``."""
    states = states or {}
    sims = [_stream(f"FEED{i}", **states.get(f"FEED{i}", {})) for i in range(inlets)]
    sims += [_stream(f"OUT{i}", **states.get(f"OUT{i}", {})) for i in range(outlets)]
    sims.append(f"<SimulationObject><Type>DWSIM.UnitOperations.UnitOperations.{unit_cls}"
                f"</Type><Name>U</Name>{props}</SimulationObject>")
    gos = [_graphic(f"FEED{i}", "MaterialStream", outs=[(0, "U", i)]) for i in range(inlets)]
    gos += [_graphic(f"OUT{i}", "MaterialStream", ins=[(0, "U", i)]) for i in range(outlets)]
    gos.append(_graphic("U", otype, ins=[(i, f"FEED{i}", 0) for i in range(inlets)],
                        outs=[(i, f"OUT{i}", 0) for i in range(outlets)]))
    comps = "".join(f"<Compound><Name>{n}</Name><CAS_Number>{c}</CAS_Number></Compound>"
                    for n, c in compounds)
    return (f"<?xml version='1.0' encoding='utf-8'?><DWSIM_Simulation_Data>"
            f"<SimulationObjects>{''.join(sims)}</SimulationObjects>"
            f"<GraphicObjects>{''.join(gos)}</GraphicObjects>"
            f"<PropertyPackages><PropertyPackage><ID>PP</ID>"
            f"<Type>DWSIM.Thermodynamics.PropertyPackages.{package}</Type>"
            f"<ComponentName>{package}</ComponentName></PropertyPackage></PropertyPackages>"
            f"<Compounds>{comps}</Compounds>{extra}</DWSIM_Simulation_Data>").encode()


def test_cooler_duty_is_negative():
    model, _ = import_dwsim(_file("Cooler", "Cooler",
                                  "<CalcMode>HeatRemoved</CalcMode><DeltaQ>50</DeltaQ>"))
    assert _units(model)["U"]["params"] == {"duty": "-50 kW"}


def test_spec_without_counterpart_uses_the_stored_outlet_state():
    data = _file("Heater", "Heater", "<CalcMode>EnergyStream</CalcMode>",
                 states={"OUT0": {"T": 350.15}})
    model, warnings = import_dwsim(data)
    assert _units(model)["U"]["params"] == {"T": "77 degC"}
    assert any("calculation mode EnergyStream" in w and "77 degC" in w for w in warnings)
    data = _file("Valve", "Valve", "<CalcMode>Kv_Liquid</CalcMode>",
                 states={"OUT0": {"P": 1.5e5}})
    assert _units(import_dwsim(data)[0])["U"]["params"] == {"P": "1.5 bar"}


def test_three_way_splitter_becomes_two_splitters():
    data = _file("Splitter", "NodeOut", "<SplitRatios><SplitRatio>0.5</SplitRatio>"
                 "<SplitRatio>0.3</SplitRatio><SplitRatio>0.2</SplitRatio></SplitRatios>",
                 outlets=3)
    model, warnings = import_dwsim(data)
    units, streams = _units(model), _streams(model)
    assert units["U"]["params"] == {"split": 0.5}
    assert units["U_2"]["params"] == {"split": 0.6}
    assert streams["OUT1"]["from"] == "U_2.out1" and streams["OUT2"]["from"] == "U_2.out2"
    result = _solve(model)
    assert result.streams["OUT2"]["molar_flow"] == pytest.approx(2.0)
    assert any("two splitters" in w for w in warnings)


def test_mixer_and_vessel_with_two_liquids():
    data = _file("Mixer", "NodeIn", "", inlets=2)
    model, _ = import_dwsim(data)
    assert [s["to"] for s in model["streams"] if s["id"].startswith("FEED")] == ["U.in", "U.in"]
    data = _file("Vessel", "Vessel", "<OverrideT>true</OverrideT>"
                 "<FlashTemperature>310.15</FlashTemperature>", outlets=3,
                 package="NRTLPropertyPackage")
    units = _units(import_dwsim(data)[0])
    assert units["U"]["type"] == "Decanter" and units["U"]["params"] == {"T": "37 degC"}


def test_unmapped_property_package_falls_back_to_auto():
    model, warnings = import_dwsim(_file("Valve", "Valve", "<DeltaP>1e5</DeltaP>",
                                         package="SteamTablesPropertyPackage"))
    assert model["property_method"] == "auto"
    assert any("SteamTablesPropertyPackage has no Fluxsheet counterpart" in w for w in warnings)
    assert _units(model)["U"]["params"] == {"dP": "1 bar"}


def test_file_with_only_unsupported_units_is_an_error():
    # (An unsupported unit inside a flowsheet: see test_unsupported_unit_becomes_feed_...)
    with pytest.raises(DwsimImportError, match="nothing to import"):
        import_dwsim(_file("OrificePlate", "OrificePlate", ""))


# --- bad files -----------------------------------------------------------------------------


@pytest.mark.parametrize("data, message", [
    (b"this is not xml", "not a valid XML file"),
    (b"<?xml version='1.0'?><Flowsheet/>", "not a DWSIM simulation"),
    (b"<DWSIM_Simulation_Data><SimulationObjects>", "not a valid XML file"),
    (b"<DWSIM_Simulation_Data/>", "has no compounds"),
    (b"PK\x03\x04garbage", "zip archive is damaged"),
])
def test_malformed_files_give_clear_errors(data, message):
    with pytest.raises(DwsimImportError, match=message):
        import_dwsim(data)


def test_entity_expansion_is_refused():
    bomb = (b'<?xml version="1.0"?><!DOCTYPE d [<!ENTITY a "aaaaaaaaaa">'
            b'<!ENTITY b "&a;&a;&a;&a;&a;&a;&a;&a;&a;&a;">]><DWSIM_Simulation_Data>&b;'
            b'</DWSIM_Simulation_Data>')
    with pytest.raises(DwsimImportError, match="DTD or entities"):
        import_dwsim(bomb)
    external = (b'<?xml version="1.0"?><!DOCTYPE d [<!ENTITY x SYSTEM "file:///etc/passwd">]>'
                b'<DWSIM_Simulation_Data>&x;</DWSIM_Simulation_Data>')
    with pytest.raises(DwsimImportError, match="DTD or entities"):
        import_dwsim(external)


def test_zip_without_simulation_and_zip_bomb(monkeypatch):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("readme.txt", "hello")
    with pytest.raises(DwsimImportError, match="no .xml simulation file"):
        import_dwsim(buf.getvalue())
    # 8 MB of zeros compress to a few kB; with a 1 MB limit the entry is refused.
    monkeypatch.setattr(dwsim, "MAX_XML_BYTES", 1024 * 1024)
    bomb = _zipped(b"<DWSIM_Simulation_Data><!--" + b"0" * (8 * 1024 * 1024) + b"-->"
                   b"</DWSIM_Simulation_Data>")
    assert len(bomb) < 100_000
    with pytest.raises(DwsimImportError, match="too large"):
        import_dwsim(bomb)
    with pytest.raises(DwsimImportError, match="entries"):
        import_dwsim(_many_entries(dwsim.MAX_ZIP_ENTRIES + 1))


def _many_entries(n: int) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for i in range(n):
            zf.writestr(f"f{i}.xml", "<DWSIM_Simulation_Data/>")
    return buf.getvalue()


def test_oversized_and_missing_files(monkeypatch, tmp_path):
    monkeypatch.setattr(dwsim, "MAX_FILE_BYTES", 1000)
    with pytest.raises(DwsimImportError, match="too large"):
        import_dwsim(DATA / "heater_flash.dwxml")
    with pytest.raises(DwsimImportError, match="too large"):
        import_dwsim(b"<" * 2000)
    with pytest.raises(DwsimImportError, match="cannot read"):
        import_dwsim(tmp_path / "missing.dwxml")


def test_unknown_compound_is_an_error():
    data = _file("Valve", "Valve", "", compounds=(("Unobtainium", "0000-00-0"),))
    with pytest.raises(DwsimImportError, match="Unobtainium") as exc:
        import_dwsim(data)
    assert "database" in str(exc.value)


# --- entry points --------------------------------------------------------------------------


def test_cli_import(tmp_path, capsys):
    out = tmp_path / "flash.json"
    assert main(["import-dwsim", str(DATA / "heater_flash.dwxml"), "-o", str(out)]) == 0
    model = FlowsheetModel.load(out)
    assert {u.type for u in model.units} >= {"Feed", "Heater", "Flash"}
    err = capsys.readouterr().err
    assert "written to" in err and "warning: energy streams left out" in err
    assert main(["import-dwsim", str(DATA / "heater_flash.dwxml")]) == 0
    assert json.loads(capsys.readouterr().out)["components"] == ["water", "ethanol"]
    bad = tmp_path / "bad.dwxml"
    bad.write_text("<nope/>")
    assert main(["import-dwsim", str(bad)]) == 1
    assert "error: not a DWSIM simulation" in capsys.readouterr().err


def test_api_import(monkeypatch):
    client = TestClient(app)
    data = _zipped((DATA / "heater_flash.dwxml").read_bytes())
    r = client.post("/api/import/dwsim", content=data)
    assert r.status_code == 200
    body = r.json()
    assert body["ok"] and body["model"]["property_method"] == "nrtl"
    assert any("Q-101" in w for w in body["warnings"])
    r = client.post("/api/import/dwsim", content=b"<nope/>")
    assert r.json()["ok"] is False and "not a DWSIM" in r.json()["error"]["message"]
    assert client.post("/api/import/dwsim", content=b"").status_code == 400
    monkeypatch.setattr(dwsim, "MAX_FILE_BYTES", 1000)
    assert client.post("/api/import/dwsim", content=data).status_code == 413
    monkeypatch.setattr(limits, "DEMO", True)
    monkeypatch.setattr(limits, "MAX_BODY_BYTES", 500)
    assert client.post("/api/import/dwsim", content=data).status_code == 413


# ---- DWSIM's own layout and real DWSIM files ---------------------------------------------------


def test_members_written_twice_are_read_from_the_one_with_data():
    # DWSIM writes a phase's Properties and Compounds twice: a type name first, the data after.
    xml = (DATA / "heater_flash.dwxml").read_text(encoding="utf-8")
    doubled = (xml.replace("<Properties>", "<Properties>DWSIM.Thermodynamics.BaseClasses."
                           "PhaseProperties</Properties>\n<Properties>")
               .replace("<Compounds>\n            <Compound>",
                        "<Compounds>System.Collections.Generic.Dictionary</Compounds>\n"
                        "<Compounds>\n            <Compound>"))
    assert doubled.count("PhaseProperties</Properties>") == xml.count("<Properties>")
    assert import_dwsim(doubled.encode("utf-8")) == import_dwsim(xml.encode("utf-8"))


REAL = json.loads((DATA / "real-samples.json").read_text(encoding="utf-8"))
REAL_DIR = os.environ.get("FLUXSHEET_DWSIM_SAMPLES")
# Fluxsheet against the states DWSIM stored: temperature (K), relative molar flow, mole
# fraction. Peng-Robinson's interaction parameters differ between the two programs.
REAL_LIMITS = {
    "Cavett's Problem.dwxmz": (1.5, 0.10, 0.01),
    "Three Phase Separator.dwxmz": (0.1, 0.03, 0.01),
    "Gibbs and Equilibrium Reactors.dwxml": (0.01, 1e-3, 1e-3),
}


@pytest.mark.skipif(not REAL_DIR, reason="set FLUXSHEET_DWSIM_SAMPLES to the folder "
                    "scripts/fetch_dwsim_samples.py filled")
@pytest.mark.parametrize("sample", REAL["samples"], ids=lambda s: s["file"])
def test_a_real_dwsim_file_solves_like_dwsim(sample):
    path = Path(REAL_DIR) / sample["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == sample["sha256"]
    model, warnings = import_dwsim(path)
    stored = dwsim.stored_streams(path)
    result = _solve(model)
    dT, dF, dx = REAL_LIMITS[sample["file"]]
    compared = 0
    for sid, st in result.streams.items():
        ref = stored.get(sid)
        if ref is None or ref["T"] is None:  # streams the import adds
            continue
        compared += 1
        assert st["T"] == pytest.approx(ref["T"], abs=dT), sid
        assert st["P"] == pytest.approx(ref["P"], rel=1e-6), sid
        if not ref["molar_flow"]:
            assert st["molar_flow"] == pytest.approx(0.0, abs=1e-9), sid
            continue
        assert st["molar_flow"] == pytest.approx(ref["molar_flow"], rel=dF), sid
        for c in model["components"]:
            assert st["mole_fractions"][c] == pytest.approx(
                ref["mole_fractions"].get(c, 0.0), abs=dx), (sid, c)
    assert compared >= 3
    assert all(c["status"] == "ok" for c in result.balance_checks)
    if sample["file"].startswith("Three Phase"):
        assert any("UNIFAC for its phase split" in w for w in warnings)
    if sample["file"].startswith("Gibbs"):
        assert _units(model)["RG-000"]["params"]["fugacities"] == "property_method"
