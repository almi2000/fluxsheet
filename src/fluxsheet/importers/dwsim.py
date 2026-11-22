"""Import DWSIM simulations (``.dwxml``, and ``.dwxmz``, which is that XML zipped).

    model, warnings = import_dwsim("plant.dwxmz")

``model`` is a Fluxsheet flowsheet (the JSON dict of flowsheet.py) and ``warnings`` lists
everything of the DWSIM file that has no Fluxsheet counterpart and was left out or
approximated: nothing is dropped silently.

The DWSIM file format, as DWSIM writes it (DWSIM.FlowsheetBase SaveToXML; objects are
serialized by reflection, so element names are .NET property names, numbers use the
invariant culture, values are SI except heat duties, which are kW)::

    <DWSIM_Simulation_Data>
      <GeneralInfo><BuildVersion>...</BuildVersion>...</GeneralInfo>
      <SimulationObjects>
        <SimulationObject>
          <Type>DWSIM.UnitOperations.UnitOperations.Heater</Type>
          <Name>HEAT-1a2b</Name>                  (the id the drawing refers to)
          <CalcMode>OutletTemperature</CalcMode><OutletTemperature>365.15</OutletTemperature> ...
        </SimulationObject>
        <SimulationObject>                        (a material stream: its solved state)
          <Type>DWSIM.Thermodynamics.Streams.MaterialStream</Type> <Name>MAT-..</Name>
          <Phases><Phase><ID>0</ID>              (0 = overall, 1 liquid, 2 vapor, ...)
            <Compounds><Compound><Name>Water</Name><MoleFraction>0.8</MoleFraction>...
            <Properties><temperature>298.15</temperature><pressure>101325</pressure>
                        <molarflow>27.8</molarflow> ...
      <GraphicObjects>
        <GraphicObject>                           (same Name; Tag is the user's label)
          <Name>HEAT-1a2b</Name><Tag>E-101</Tag><ObjectType>Heater</ObjectType><X>..</X><Y>..</Y>
          <InputConnectors><Connector IsAttached="true" AttachedFromObjID="MAT-.."
                                      AttachedFromConnIndex="0" .../></InputConnectors>
          <OutputConnectors><Connector IsAttached="true" AttachedToObjID=".." AttachedToConnIndex="0"/>
          <EnergyConnector>...</EnergyConnector>
      <PropertyPackages><PropertyPackage><ID>..</ID><Type>...PengRobinsonPropertyPackage</Type>
                        <ComponentName>Peng-Robinson (PR)</ComponentName>...
      <Compounds><Compound><Name>Water</Name><CAS_Number>7732-18-5</CAS_Number>...
      <ReactionSets>, <Reactions>, <OptimizationCases>, <SensitivityAnalysis>, ...

In DWSIM, streams are objects between units; in Fluxsheet they are connections. A stream
without an upstream unit becomes a Feed (at its state in the file), one without a
downstream unit a Product. A DWSIM unit without a Fluxsheet counterpart is left out; the
streams around it then become products (its inlets) and feeds at the state DWSIM last
computed (its outlets), so the rest of the flowsheet still solves. Specs that Fluxsheet
cannot express (energy-stream or curve modes, ...) are replaced by the outlet state DWSIM
computed, with a warning.

Safety: the XML is parsed with defusedxml (a dependency of Fluxsheet), which refuses DTDs,
entity declarations and external references, so neither "billion laughs" nor external
entity tricks work. Files and the XML inside a zip are size-limited, the zip entry is read
with a hard cap rather than trusting its header, so zip bombs are refused.
"""

from __future__ import annotations

import io
import math
import re
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
from xml.etree.ElementTree import Element

from defusedxml import DefusedXmlException
from defusedxml import ElementTree as SafeET

from ..compounds import UnknownCompoundError, lookup
from ..quantities import UnitError, to_si

# Size limits. DWSIM files with many compounds reach a few MB (each compound carries its
# full property data); these leave ample room while refusing absurd or malicious input.
MAX_FILE_BYTES = 64 * 1024 * 1024  # the file as given (zip or XML)
MAX_XML_BYTES = 64 * 1024 * 1024  # the XML, after unzipping
MAX_ZIP_ENTRIES = 64

ROOT = "DWSIM_Simulation_Data"


class DwsimImportError(ValueError):
    """The file cannot be imported (not DWSIM, malformed, too large, unknown compounds)."""

    def __init__(self, message: str, hint: str = ""):
        super().__init__(message)
        self.message = message
        self.hint = hint

    def to_dict(self) -> dict:
        return {"message": self.message, "where": None, "hint": self.hint}


class _Unsupported(Exception):
    """A unit (or its spec) cannot be expressed in Fluxsheet; it is left out."""


# --- reading -------------------------------------------------------------------------------


def _read_source(source: str | Path | bytes | bytearray) -> tuple[bytes, str | None]:
    if isinstance(source, (bytes, bytearray, memoryview)):
        data = bytes(source)
        if len(data) > MAX_FILE_BYTES:
            raise DwsimImportError(_too_large("the file", len(data), MAX_FILE_BYTES))
        return data, None
    path = Path(source)
    try:
        size = path.stat().st_size
        if size > MAX_FILE_BYTES:
            raise DwsimImportError(_too_large(path.name, size, MAX_FILE_BYTES))
        with path.open("rb") as f:
            data = f.read(MAX_FILE_BYTES + 1)
    except OSError as exc:
        raise DwsimImportError(f"cannot read {path}: {exc.strerror or exc}") from exc
    if len(data) > MAX_FILE_BYTES:
        raise DwsimImportError(_too_large(path.name, len(data), MAX_FILE_BYTES))
    return data, path.name


def _too_large(what: str, size: int, limit: int) -> str:
    return f"{what} is too large ({size / 2**20:.1f} MB; the limit is {limit / 2**20:.0f} MB)"


def _unzip(data: bytes) -> tuple[bytes, list[str]]:
    """The XML document inside a .dwxmz (DWSIM's loader takes the entry ending in .xml)."""
    notes = []
    try:
        zf = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile as exc:
        raise DwsimImportError(f"not a valid .dwxmz file: the zip archive is damaged ({exc})") \
            from exc
    with zf:
        infos = zf.infolist()
        if len(infos) > MAX_ZIP_ENTRIES:
            raise DwsimImportError(f"the zip archive has {len(infos)} entries; a DWSIM file has "
                                   f"one or two (at most {MAX_ZIP_ENTRIES} are accepted)")
        xmls = [i for i in infos if not i.is_dir() and i.filename.lower().endswith(".xml")]
        if not xmls:
            raise DwsimImportError("the zip archive contains no .xml simulation file",
                                   hint="a .dwxmz file is a zip with one .xml document inside")
        info = xmls[0]
        if len(xmls) > 1:
            notes.append(f"the archive holds {len(xmls)} XML files; imported {info.filename}")
        if info.file_size > MAX_XML_BYTES:
            raise DwsimImportError(_too_large(f"{info.filename} (unzipped)", info.file_size,
                                              MAX_XML_BYTES))
        try:
            with zf.open(info) as f:  # read with a hard cap: the header's size may lie
                xml = f.read(MAX_XML_BYTES + 1)
        except (RuntimeError, NotImplementedError, zipfile.BadZipFile, OSError, EOFError,
                ValueError) as exc:  # encrypted, unsupported compression, corrupt data
            raise DwsimImportError(f"cannot unpack {info.filename}: {exc}") from exc
        if len(xml) > MAX_XML_BYTES:
            raise DwsimImportError(f"{info.filename} unzips to more than "
                                   f"{MAX_XML_BYTES / 2**20:.0f} MB, the limit")
    return xml, notes


def _parse(xml: bytes) -> Element:
    try:
        # defusedxml: no DTDs, no entity declarations, no external references.
        root = SafeET.fromstring(xml, forbid_dtd=True)
    except DefusedXmlException as exc:
        raise DwsimImportError(
            "refused: the XML declares a DTD or entities, which DWSIM files never contain",
            hint="open the file in DWSIM and save it again") from exc
    except ET.ParseError as exc:
        raise DwsimImportError(f"not a valid XML file ({exc})",
                               hint="the file may be truncated or not a DWSIM file") from exc
    if _local(root.tag) != ROOT:
        raise DwsimImportError(
            f"not a DWSIM simulation: the XML root is <{_local(root.tag)}>, not <{ROOT}>",
            hint="export the simulation from DWSIM as .dwxml or .dwxmz")
    return root


# --- XML helpers ---------------------------------------------------------------------------


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _child(el: Element | None, name: str) -> Element | None:
    if el is None:
        return None
    for c in el:
        if _local(c.tag) == name:
            return c
    return None


def _children(el: Element | None, name: str) -> list[Element]:
    return [] if el is None else [c for c in el if _local(c.tag) == name]


def _s(el: Element | None, name: str) -> str | None:
    c = _child(el, name)
    if c is None or c.text is None:
        return None
    text = c.text.strip()
    return text or None


def _num(text: str | None) -> float | None:
    if text is None:
        return None
    try:  # DWSIM writes numbers in the invariant culture ("R" format)
        v = float(text.strip())
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def _f(el: Element | None, name: str) -> float | None:
    return _num(_s(el, name))


def _b(el: Element | None, name: str) -> bool:
    return (_s(el, name) or "").lower() in ("true", "1")


def _enum(el: Element | None, name: str, members: tuple[str, ...], default: str) -> str:
    """An enum property: DWSIM writes the member name; tolerate the integer value too."""
    text = _s(el, name)
    if text is None:
        return default
    if text.isdigit():
        i = int(text)
        return members[i] if i < len(members) else text
    return text


def _attr_bool(el: Element, name: str) -> bool:
    return (el.get(name) or "").lower() in ("true", "1")


def _class_name(type_name: str | None) -> str:
    """``DWSIM.UnitOperations.UnitOperations.Heater`` -> ``Heater``."""
    if not type_name:
        return ""
    return type_name.split("+")[0].split(",")[0].strip().rsplit(".", 1)[-1]


# --- values with units ---------------------------------------------------------------------


def _g(x: float) -> str:
    return f"{x:.8g}"


def _T(K: float) -> str:
    return f"{_g(K - 273.15)} degC"


def _P(Pa: float) -> str:
    return f"{_g(Pa / 1e5)} bar"


def _kW(kW: float) -> str:
    return f"{_g(kW)} kW"


def _kmolh(mol_s: float) -> str:
    return f"{_g(mol_s * 3.6)} kmol/h"


# --- DWSIM vocabulary ----------------------------------------------------------------------

# ObjectType of a drawing object -> class name of its simulation object, for files whose
# SimulationObject <Type> is missing or unfamiliar.
_OBJECT_TYPES = {
    "NodeIn": "Mixer", "NodeOut": "Splitter", "RCT_Conversion": "Reactor_Conversion",
    "RCT_Equilibrium": "Reactor_Equilibrium", "RCT_Gibbs": "Reactor_Gibbs",
    "RCT_CSTR": "Reactor_CSTR", "RCT_PFR": "Reactor_PFR", "OT_Recycle": "Recycle",
    "OT_EnergyRecycle": "EnergyRecycle", "OT_Adjust": "Adjust", "OT_Spec": "Spec",
    "TPVessel": "Vessel", "NodeEn": "EnergyMixer",
}

# DWSIM objects without a Fluxsheet counterpart, by class name, with what they are.
_UNSUPPORTED = {
    "DistillationColumn": "rigorous distillation column (its specs do not translate; a "
                          "ShortcutColumn or Fluxsheet's DistillationColumn can be set up by hand)",
    "AbsorptionColumn": "absorption column", "RefluxedAbsorber": "refluxed absorber",
    "ReboiledAbsorber": "reboiled absorber",
    "Reactor_Equilibrium": "equilibrium reactor (DWSIM's equilibrium reactions are not imported)",
    "Reactor_CSTR": "CSTR (DWSIM's kinetic reactions are not imported)",
    "Reactor_PFR": "PFR (DWSIM's kinetic reactions are not imported)",
    "Reactor_ReaktoroGibbs": "Reaktoro Gibbs reactor",
    "Pipe": "pipe segment (DWSIM's pipe profile is not imported)",
    "OrificePlate": "orifice plate", "CustomUO": "script unit operation",
    "ExcelUO": "spreadsheet unit operation", "CapeOpenUO": "CAPE-OPEN unit operation",
    "Flowsheet": "sub-flowsheet", "FlowsheetUO": "sub-flowsheet",
    "SolidsSeparator": "solids separator", "Filter": "filter", "AirCooler2": "air cooler",
    "WindTurbine": "wind turbine", "HydroelectricTurbine": "hydroelectric turbine",
    "SolarPanel": "solar panel", "PEMFuelCell": "fuel cell", "WaterElectrolyzer": "electrolyzer",
}

_LOGICAL = {
    "Adjust": "adjust (controller) block", "Spec": "spec block",
    "EnergyRecycle": "energy recycle block", "PIDController": "PID controller",
    "PythonController": "Python controller", "Switch": "switch", "Input": "input box",
    "AnalogGauge": "gauge", "DigitalGauge": "gauge", "LevelGauge": "gauge",
    "EnergyMixer": "energy mixer",
}

# DWSIM property package class -> (Fluxsheet method, note when only approximated).
_PR = "is imported as Peng-Robinson"
_UNIFAC = "is imported as UNIFAC (Dortmund)"
_PACKAGES: dict[str, tuple[str, str | None]] = {
    "PengRobinsonPropertyPackage": ("peng-robinson", None),
    "PengRobinson1978PropertyPackage": ("peng-robinson", f"Peng-Robinson 1978 {_PR}"),
    "PengRobinson1978AdvancedPropertyPackage": ("peng-robinson", f"Peng-Robinson 1978 {_PR}"),
    "PRSV2PropertyPackage": ("peng-robinson", f"PRSV2 {_PR}"),
    "PRSV2VLPropertyPackage": ("peng-robinson", f"PRSV2-VL {_PR}"),
    "PengRobinsonLKPropertyPackage": ("peng-robinson", f"Peng-Robinson/Lee-Kesler {_PR}"),
    "SRKPropertyPackage": ("srk", None),
    "SoaveRedlichKwongAdvancedPropertyPackage": ("srk", None),
    "NRTLPropertyPackage": ("nrtl", None),
    "UNIQUACPropertyPackage": ("uniquac", None),
    "WilsonPropertyPackage": ("wilson", None),
    "MODFACPropertyPackage": ("unifac", None),  # Fluxsheet's UNIFAC is the Dortmund version
    "UNIFACPropertyPackage": ("unifac", f"original UNIFAC {_UNIFAC}"),
    "UNIFACLLPropertyPackage": ("unifac", f"UNIFAC-LL {_UNIFAC}"),
    "NISTMFACPropertyPackage": ("unifac", f"Modified UNIFAC (NIST) {_UNIFAC}"),
    "RaoultPropertyPackage": ("ideal", None),
}
# The same by DWSIM's display name, for files whose <Type> is unfamiliar.
_PACKAGE_NAMES = {
    "peng-robinson (pr)": "PengRobinsonPropertyPackage",
    "peng-robinson 1978 (pr78)": "PengRobinson1978PropertyPackage",
    "soave-redlich-kwong (srk)": "SRKPropertyPackage", "nrtl": "NRTLPropertyPackage",
    "uniquac": "UNIQUACPropertyPackage", "wilson": "WilsonPropertyPackage",
    "unifac": "UNIFACPropertyPackage", "unifac-ll": "UNIFACLLPropertyPackage",
    "modified unifac (dortmund)": "MODFACPropertyPackage",
    "modified unifac (nist)": "NISTMFACPropertyPackage", "raoult's law": "RaoultPropertyPackage",
}
_WITH_PARAMETERS = ("peng-robinson", "srk", "nrtl", "uniquac", "wilson")

_HEATER_MODES = ("HeatAdded", "OutletTemperature", "EnergyStream", "OutletVaporFraction",
                 "TemperatureChange", "HeatAddedRemoved")
_COOLER_MODES = ("HeatRemoved", "OutletTemperature", "OutletVaporFraction", "TemperatureChange",
                 "EnergyStream")
_PUMP_MODES = ("Delta_P", "OutletPressure", "EnergyStream", "Curves", "Power")
_COMPRESSOR_MODES = ("OutletPressure", "Delta_P", "EnergyStream", "PowerRequired", "Head",
                     "Curves", "PressureRatio")
_EXPANDER_MODES = ("OutletPressure", "Delta_P", "PowerGenerated", "Head", "Curves",
                   "PressureRatio")
_VALVE_MODES = ("DeltaP", "OutletPressure", "Kv_Liquid", "Kv_Gas", "Kv_Steam", "Kv_General")
_HX_MODES = ("CalcTempHotOut", "CalcTempColdOut", "CalcBothTemp", "CalcBothTemp_UA", "CalcArea",
             "ShellandTube_Rating", "ShellandTube_CalcFoulingFactor", "PinchPoint",
             "ThermalEfficiency", "OutletVaporFraction1", "OutletVaporFraction2")
_REACTOR_MODES = ("Isothermic", "Adiabatic", "OutletTemperature", "NonIsothermalNonAdiabatic")


# --- the model of the file -----------------------------------------------------------------


@dataclass
class _Conn:
    index: int
    other: str
    other_index: int
    energy: bool


@dataclass
class _Graphic:
    name: str
    tag: str
    object_type: str
    x: float
    y: float
    inputs: list[_Conn]
    outputs: list[_Conn]


@dataclass
class _Obj:
    name: str
    cls: str
    el: Element | None
    graphic: _Graphic | None
    id: str = ""

    @property
    def label(self) -> str:
        return self.graphic.tag if self.graphic and self.graphic.tag else self.name


@dataclass
class _State:
    T: float | None
    P: float | None
    flow: float | None  # mol/s
    x: dict[str, float]  # Fluxsheet component -> mole fraction
    vf: float | None
    spec: str


@dataclass
class _Unit:
    """A converted unit: Fluxsheet unit dicts and where DWSIM's connectors land."""

    units: list[dict]
    inlets: dict[int, str]  # DWSIM inlet connector index -> "UNIT.port"
    outlets: dict[int, str]  # DWSIM outlet connector index -> "UNIT.port"
    streams: list[dict] = field(default_factory=list)  # internal connections


def _graphic(el: Element) -> _Graphic | None:
    name = _s(el, "Name")
    if not name:
        return None

    def conns(section: str, prefix: str) -> list[_Conn]:
        out = []
        for i, c in enumerate(_children(_child(el, section), "Connector")):
            other = c.get(f"Attached{prefix}ObjID")
            if _attr_bool(c, "IsAttached") and other:
                idx = _num(c.get(f"Attached{prefix}ConnIndex"))
                out.append(_Conn(i, other, int(idx) if idx is not None else 0,
                                 _attr_bool(c, f"Attached{prefix}EnergyConn")))
        return out

    return _Graphic(name=name, tag=(_s(el, "Tag") or "").strip(),
                    object_type=_s(el, "ObjectType") or "", x=_f(el, "X") or 0.0,
                    y=_f(el, "Y") or 0.0, inputs=conns("InputConnectors", "From"),
                    outputs=conns("OutputConnectors", "To"))


class _Importer:
    def __init__(self, root: Element, name: str, notes: list[str]):
        self.root = root
        self.name = name
        self.warnings: list[str] = list(notes)
        self.components: list[str] = []
        self.dwsim_names: dict[str, str] = {}  # DWSIM compound name (lower) -> Fluxsheet name
        self.mw: dict[str, float] = {}  # Fluxsheet name -> kg/mol
        self.objects: dict[str, _Obj] = {}
        self.states: dict[str, _State] = {}
        self.src: dict[str, tuple[str, int]] = {}  # stream -> (unit, outlet index)
        self.dst: dict[str, tuple[str, int]] = {}  # stream -> (unit, inlet index)
        self.reactions: list[dict] = []
        self.reaction_ids: dict[str, str] = {}  # DWSIM reaction id -> Fluxsheet id
        self.unit_ids: set[str] = set()
        self.stream_ids: set[str] = set()

    def warn(self, text: str) -> None:
        if text not in self.warnings:
            self.warnings.append(text)

    # -- the steps --

    def run(self) -> tuple[dict, list[str]]:
        self._compounds()
        self._objects()
        self._connections()
        self._collapse_recycles()
        methods = self._property_packages()
        converted = self._convert_units()
        model = self._assemble(converted, methods)
        from ..flowsheet import FlowsheetModel  # validate what we built

        try:
            FlowsheetModel.model_validate(model)
        except Exception as exc:  # a bug here, not in the user's file
            raise DwsimImportError(f"the imported flowsheet is not valid: {exc}") from exc
        return model, self.warnings

    def _compounds(self) -> None:
        entries = [(_s(c, "Name"), _s(c, "CAS_Number"))
                   for c in _children(_child(self.root, "Compounds"), "Compound")]
        if not entries:  # fall back on the compounds named in the streams
            seen = []
            for so in _children(_child(self.root, "SimulationObjects"), "SimulationObject"):
                for ph in _children(_child(so, "Phases"), "Phase"):
                    for c in _children(_child(ph, "Compounds"), "Compound"):
                        n = _s(c, "Name")
                        if n and n not in seen:
                            seen.append(n)
            entries = [(n, None) for n in seen]
        if not entries:
            raise DwsimImportError("the simulation has no compounds",
                                   hint="add compounds in DWSIM and save the file again")
        unknown = []
        for dname, cas in entries:
            compound = None
            for ident in (cas, dname):
                if not ident:
                    continue
                try:
                    compound = lookup(ident)
                    break
                except UnknownCompoundError:
                    continue
            if compound is None:
                unknown.append(f"{dname or '?'}" + (f" (CAS {cas})" if cas else ""))
                continue
            fname = compound.name
            try:  # the name must resolve back to the same compound, else use the CAS number
                if lookup(fname).cas != compound.cas:
                    fname = compound.cas
            except UnknownCompoundError:
                fname = compound.cas
            if not compound.simulation_ready:
                self.warn(f"compound {dname} lacks critical data in Fluxsheet's database; "
                          "results involving it may fail")
            if fname not in self.components:
                self.components.append(fname)
                self.mw[fname] = compound.mw / 1000.0
            if dname:
                self.dwsim_names[dname.lower()] = fname
        if unknown:
            raise DwsimImportError(
                f"compounds not found in Fluxsheet's database: {', '.join(unknown)}",
                hint="user-defined or petroleum compounds of DWSIM cannot be imported; replace "
                "them in DWSIM by database compounds, or define them in Fluxsheet")

    def comp(self, dwsim_name: str | None) -> str | None:
        return self.dwsim_names.get((dwsim_name or "").strip().lower())

    def _objects(self) -> None:
        graphics = {}
        annotations = Counter()
        for el in _children(_child(self.root, "GraphicObjects"), "GraphicObject"):
            g = _graphic(el)
            if g is not None:
                graphics[g.name] = g
        for so in _children(_child(self.root, "SimulationObjects"), "SimulationObject"):
            name = _s(so, "Name") or _s(so, "ComponentName")
            if not name:
                continue
            g = graphics.get(name)
            cls = _class_name(_s(so, "Type"))
            if g is not None and cls not in _KNOWN_CLASSES:
                cls = _OBJECT_TYPES.get(g.object_type, g.object_type) or cls
            self.objects[name] = _Obj(name, cls, so, g)
        for name, g in graphics.items():
            if name not in self.objects:
                if g.object_type.startswith("GO_") or not g.object_type:
                    annotations[g.object_type.removeprefix("GO_").lower() or "drawing"] += 1
                else:
                    self.warn(f"{g.tag or name} ({g.object_type}) has no simulation data in the "
                              "file; left out")
        if annotations:
            what = ", ".join(f"{n} {k}" for k, n in sorted(annotations.items()))
            self.warn(f"drawing annotations left out ({what})")
        for o in self.objects.values():
            if o.cls == "MaterialStream":
                self.states[o.name] = self._state(o)

    def _state(self, o: _Obj) -> _State:
        phases = _children(_child(o.el, "Phases"), "Phase")

        def phase(pid: int) -> Element | None:
            names = {0: "mixture", 2: "vapor"}
            for ph in phases:
                if _f(ph, "ID") == pid or (_s(ph, "Name") or "").lower() == names.get(pid):
                    return ph
            return None

        mix = phase(0)
        props = _child(mix, "Properties")
        T, P = _f(props, "temperature"), _f(props, "pressure")
        flow = _f(props, "molarflow")
        x: dict[str, float] = {}
        w: dict[str, float] = {}
        for c in _children(_child(mix, "Compounds"), "Compound"):
            name = self.comp(_s(c, "Name"))
            if name is None:
                continue
            if (v := _f(c, "MoleFraction")) is not None:
                x[name] = x.get(name, 0.0) + v
            if (v := _f(c, "MassFraction")) is not None:
                w[name] = w.get(name, 0.0) + v
        if sum(x.values()) <= 0.0 and sum(w.values()) > 0.0:  # from mass fractions
            moles = {k: v / self.mw[k] for k, v in w.items() if self.mw.get(k)}
            total = sum(moles.values())
            x = {k: v / total for k, v in moles.items()} if total > 0 else {}
        if flow is None and (m := _f(props, "massflow")) is not None and x:
            mw = sum(x[k] * self.mw.get(k, 0.0) for k in x)
            flow = m / mw if mw > 0 else None
        total = sum(x.values())
        if total > 0:
            x = {k: v / total for k, v in x.items()}
        vf = _f(_child(phase(2), "Properties"), "molarfraction")
        return _State(T, P, flow, x, vf, _s(o.el, "SpecType") or "Temperature_and_Pressure")

    def _connections(self) -> None:
        edges = set()
        for o in self.objects.values():
            g = o.graphic
            if g is None:
                continue
            for c in g.outputs:
                edges.add((o.name, c.index, c.other, c.other_index, c.energy))
            for c in g.inputs:
                edges.add((c.other, c.other_index, o.name, c.index, c.energy))
        for src, si, dst, di, energy in sorted(edges, key=str):
            a, b = self.objects.get(src), self.objects.get(dst)
            if a is None or b is None or energy:
                continue
            if b.cls == "MaterialStream" and a.cls != "MaterialStream":
                if self.src.setdefault(dst, (src, si)) != (src, si):
                    self.warn(f"stream {b.label} has two sources in the file; kept the first")
            elif (a.cls == "MaterialStream" and b.cls != "MaterialStream"
                  and self.dst.setdefault(src, (dst, di)) != (dst, di)):
                self.warn(f"stream {a.label} has two destinations in the file; kept the first")

    def _collapse_recycles(self) -> None:
        """DWSIM marks tear streams with a Recycle block; Fluxsheet finds them itself."""
        for r in [o for o in self.objects.values() if o.cls == "Recycle"]:
            ins = [s for s, (u, _) in self.dst.items() if u == r.name]
            outs = [s for s, (u, _) in self.src.items() if u == r.name]
            if len(ins) == 1 and len(outs) == 1:
                a, b = ins[0], outs[0]
                if b in self.dst:
                    self.dst[a] = self.dst.pop(b)
                else:
                    del self.dst[a]
                del self.src[b]
                self.objects[b].cls = "MergedStream"
                self.warn(f"recycle block {r.label} left out: Fluxsheet finds and converges "
                          f"recycles itself (stream {self.objects[b].label} merged into "
                          f"{self.objects[a].label})")
            else:
                self.warn(f"recycle block {r.label} is not connected to one inlet and one "
                          "outlet stream; left out")

    def _property_packages(self) -> dict[str, str]:
        """Fluxsheet method per DWSIM package (by ID and by display name)."""
        out: dict[str, str] = {}
        for pp in _children(_child(self.root, "PropertyPackages"), "PropertyPackage"):
            cls = _class_name(_s(pp, "Type"))
            label = _s(pp, "ComponentName") or cls or "?"
            if cls not in _PACKAGES:
                cls = _PACKAGE_NAMES.get(label.lower(), cls)
            mapped = _PACKAGES.get(cls)
            if mapped is None:
                method = "auto"
                self.warn(f"property package {label} has no Fluxsheet counterpart; the property "
                          "method is set to auto (Fluxsheet chooses one for the components)")
            else:
                method, note = mapped
                if note:
                    self.warn(f"property package {label}: {note}")
                if method in _WITH_PARAMETERS:
                    self.warn(f"property package {label}: interaction parameters and package "
                              "options stored in the DWSIM file are not imported; Fluxsheet "
                              "uses its own parameter database")
            for key in (_s(pp, "ID"), _s(pp, "ComponentName"), _s(pp, "Tag")):
                if key:
                    out.setdefault(key, method)
        return out

    # -- units --

    def _convert_units(self) -> dict[str, _Unit]:
        converted: dict[str, _Unit] = {}
        energy = []
        for o in self.objects.values():
            if o.cls in ("MaterialStream", "MergedStream", "Recycle"):
                continue
            if o.cls == "EnergyStream":
                energy.append(o.label)
                continue
            if o.cls in _LOGICAL:
                self.warn(f"{_LOGICAL[o.cls]} {o.label} left out")
                continue
            o.id = self._unit_id(o.label)
            convert = _CONVERTERS.get(o.cls)
            try:
                if convert is None:
                    what = _UNSUPPORTED.get(o.cls, f"unit type {o.cls or '?'}")
                    raise _Unsupported(f"{what} has no Fluxsheet counterpart")
                unit = convert(self, o)
                self._check_inlets(o, unit)
            except _Unsupported as exc:
                self.unit_ids.discard(o.id)
                self.warn(f"{o.label} left out: {exc}; the streams around it become feeds "
                          "and products at the state stored in the file")
                continue
            converted[o.name] = unit
        if energy:
            self.warn(f"energy streams left out ({', '.join(energy)}): Fluxsheet units compute "
                      "their duty and power themselves")
        return converted

    def _check_inlets(self, o: _Obj, unit: _Unit) -> None:
        from ..units import UNIT_TYPES

        main = unit.units[0]
        connected = {unit.inlets[i] for s, (u, i) in self.dst.items()
                     if u == o.name and i in unit.inlets}
        connected |= {s["to"] for s in unit.streams}
        for port in UNIT_TYPES[main["type"]].inlets:
            if f"{main['id']}.{port.name}" not in connected:
                raise _Unsupported(f"its inlet {port.name!r} is not connected in the file")

    def _unit_id(self, label: str) -> str:
        base = re.sub(r"\s*\.\s*", "_", label).strip() or "UNIT"
        uid, n = base, 2
        while uid in self.unit_ids:
            uid, n = f"{base}_{n}", n + 1
        self.unit_ids.add(uid)
        return uid

    def _stream_id(self, label: str) -> str:
        base = label.strip() or "S"
        sid, n = base, 2
        while sid in self.stream_ids:
            sid, n = f"{base}_{n}", n + 1
        self.stream_ids.add(sid)
        return sid

    def inlet_state(self, o: _Obj, index: int = 0) -> _State | None:
        for s, (u, i) in self.dst.items():
            if u == o.name and i == index:
                return self.states.get(s)
        return None

    def outlet_state(self, o: _Obj, index: int = 0) -> _State | None:
        for s, (u, i) in self.src.items():
            if u == o.name and i == index:
                return self.states.get(s)
        return None

    def attached_outlets(self, o: _Obj) -> list[int]:
        return sorted(i for u, i in self.src.values() if u == o.name)

    def stored_outlet_T(self, o: _Obj, why: str, index: int = 0) -> str:
        st = self.outlet_state(o, index)
        if st is None or st.T is None:
            raise _Unsupported(f"{why} has no Fluxsheet counterpart and the file stores no "
                               "outlet temperature to use instead")
        self.warn(f"{o.label}: {why} has no Fluxsheet counterpart; the outlet temperature "
                  f"DWSIM computed ({_T(st.T)}) is used as the spec")
        return _T(st.T)

    def stored_outlet_P(self, o: _Obj, why: str, index: int = 0) -> str:
        st = self.outlet_state(o, index)
        if st is None or st.P is None:
            raise _Unsupported(f"{why} has no Fluxsheet counterpart and the file stores no "
                               "outlet pressure to use instead")
        self.warn(f"{o.label}: {why} has no Fluxsheet counterpart; the outlet pressure "
                  f"DWSIM computed ({_P(st.P)}) is used as the spec")
        return _P(st.P)

    # -- assembling the flowsheet --

    def _assemble(self, converted: dict[str, _Unit], methods: dict[str, str]) -> dict:
        used: Counter = Counter()
        for o in self.objects.values():
            ref = _s(o.el, "PropertyPackage")
            if ref and ref in methods and (o.name in converted or o.cls == "MaterialStream"):
                used[methods[ref]] += 1
        method = used.most_common(1)[0][0] if used else next(iter(methods.values()), "auto")
        if len(set(methods.values())) > 1:
            self.warn("the file has several property packages; the flowsheet uses "
                      f"{method}, units on another package get their own property method")

        feeds, middle, products, streams = [], [], [], []
        positions = []
        for o in self.objects.values():
            if o.name not in converted:
                continue
            unit = converted[o.name]
            ref = _s(o.el, "PropertyPackage")
            for k, u in enumerate(unit.units):
                if o.graphic:
                    u["ui"] = {"x": o.graphic.x + 120 * k, "y": o.graphic.y}
                if ref in methods and methods[ref] != method:
                    u["property_method"] = methods[ref]
                middle.append(u)
            streams.extend(unit.streams)
        for o in self.objects.values():
            if o.cls != "MaterialStream":
                continue
            up, down = self.src.get(o.name), self.dst.get(o.name)
            source = target = None
            if up and up[0] in converted:
                source = converted[up[0]].outlets.get(up[1])
                if source is None:
                    self.warn(f"stream {o.label} leaves outlet {up[1]} of "
                              f"{self.objects[up[0]].label}, which has no Fluxsheet "
                              "counterpart; it becomes a feed")
            if down and down[0] in converted:
                target = converted[down[0]].inlets.get(down[1])
                if target is None:
                    self.warn(f"stream {o.label} enters inlet {down[1]} of "
                              f"{self.objects[down[0]].label}, which has no Fluxsheet "
                              "counterpart; it becomes a product")
            if source is None and target is None:
                if not up and not down:
                    self.warn(f"stream {o.label} is not connected to any unit; left out")
                elif (up and up[0] not in converted) and (down and down[0] not in converted):
                    self.warn(f"stream {o.label} runs between left-out units; left out")
                elif up and up[0] not in converted and not down:
                    self.warn(f"stream {o.label}, a product of a left-out unit, is left out")
                elif down and down[0] not in converted and not up:
                    self.warn(f"stream {o.label}, a feed of a left-out unit, is left out")
                continue
            sid = self._stream_id(o.label)
            ui = {"x": o.graphic.x, "y": o.graphic.y} if o.graphic else None
            if source is None:
                feed = self._feed(o)
                if feed is None:
                    continue
                if ui:
                    feed["ui"] = ui
                feeds.append(feed)
                source = f"{feed['id']}.out"
                if up:
                    self.warn(f"stream {o.label}, an outlet of the left-out unit "
                              f"{self.objects[up[0]].label}, becomes a feed at its stored state")
            if target is None:
                product = {"id": self._unit_id(o.label), "type": "Product"}
                if ui:
                    product["ui"] = ui
                products.append(product)
                target = f"{product['id']}.in"
            streams.append({"id": sid, "from": source, "to": target})
        units = feeds + middle + products
        for u in units:
            if "ui" in u:
                positions.append((u["ui"]["x"], u["ui"]["y"]))
        if positions:  # move the drawing to the origin, keep DWSIM's arrangement
            x0 = min(p[0] for p in positions) - 40
            y0 = min(p[1] for p in positions) - 40
            for u in units:
                if "ui" in u:
                    u["ui"] = {"x": round(u["ui"]["x"] - x0), "y": round(u["ui"]["y"] - y0)}
        if not middle and not feeds:
            raise DwsimImportError("nothing to import: the file has no unit operations or "
                                   "feeds that Fluxsheet can represent")
        info = _child(self.root, "GeneralInfo")
        build = _s(info, "BuildVersion")
        model: dict[str, Any] = {
            "format": "fluxsheet/1",
            "name": self.name,
            "description": "Imported from DWSIM" + (f" {build}" if build else "")
                           + " by fluxsheet import-dwsim.",
            "components": self.components,
            "property_method": method,
        }
        if self.reactions:
            model["reactions"] = self.reactions
        model["units"] = units
        model["streams"] = streams
        for section, what in (("OptimizationCases", "optimization cases"),
                              ("SensitivityAnalysis", "sensitivity studies"),
                              ("ScriptItems", "scripts"), ("ChartItems", "charts")):
            n = len(list(_child(self.root, section) or ()))
            if n:
                self.warn(f"{n} DWSIM {what} left out")
        return model

    def _feed(self, o: _Obj) -> dict | None:
        st = self.states.get(o.name)
        if st is None or st.P is None or (st.T is None and st.vf is None):
            self.warn(f"stream {o.label} would be a feed but the file stores no temperature "
                      "and pressure for it; left out")
            return None
        flow = st.flow or 0.0
        flows = {k: _kmolh(flow * v) for k, v in st.x.items() if v > 0.0}
        if flow <= 0.0 or not flows:
            self.warn(f"feed {o.label} has no flow in the file")
            flows = {self.components[0]: "0 kmol/h"}
        params: dict[str, Any] = {}
        spec = st.spec
        if spec in ("Pressure_and_VaporFraction", "Temperature_and_VaporFraction") \
                and st.vf is not None:
            params["VF"] = round(st.vf, 10)
            if spec == "Temperature_and_VaporFraction":
                self.warn(f"feed {o.label}: given by temperature and vapor fraction in DWSIM; "
                          "imported with its vapor fraction and the pressure DWSIM computed")
        else:
            if st.T is None:
                self.warn(f"stream {o.label} has no stored temperature; left out")
                return None
            params["T"] = _T(st.T)
            if spec not in ("Temperature_and_Pressure", ""):
                self.warn(f"feed {o.label}: specified by {spec.replace('_', ' ').lower()} in "
                          "DWSIM; imported at the temperature DWSIM computed")
        params["P"] = _P(st.P)
        params["flows"] = flows
        return {"id": self._unit_id(o.label), "type": "Feed", "params": params}


# --- unit converters -----------------------------------------------------------------------


def _unit(o: _Obj, type_: str, params: dict, inlets: dict[int, str] | None = None,
          outlets: dict[int, str] | None = None) -> _Unit:
    u = {"id": o.id, "type": type_}
    if params:
        u["params"] = params
    return _Unit([u], {i: f"{o.id}.{p}" for i, p in (inlets or {0: "in"}).items()},
                 {i: f"{o.id}.{p}" for i, p in (outlets or {0: "out"}).items()})


def _pressure_drop(params: dict, el: Element | None, name: str = "DeltaP") -> None:
    dp = _f(el, name)
    if dp and dp > 0:
        params["dP"] = _P(dp)
    elif dp and dp < 0:  # a pressure rise: Fluxsheet takes the outlet pressure instead
        raise _Unsupported(f"a negative pressure drop ({_P(dp)}) has no Fluxsheet counterpart")


def _heater(imp: _Importer, o: _Obj, cooler: bool) -> _Unit:
    el = o.el
    modes = _COOLER_MODES if cooler else _HEATER_MODES
    mode = _enum(el, "CalcMode", modes, modes[0])
    p: dict[str, Any] = {}
    if mode in ("HeatAdded", "HeatRemoved", "HeatAddedRemoved"):
        dq = _f(el, "DeltaQ")
        if dq is None:
            raise _Unsupported("no heat duty given")
        eff = _f(el, "Eficiencia")
        if eff is not None and abs(eff - 100.0) > 1e-9:
            dq *= eff / 100.0
            imp.warn(f"{o.label}: efficiency {_g(eff)} % folded into the duty")
        p["duty"] = _kW(-dq if cooler else dq)
    elif mode == "OutletTemperature":
        T = _f(el, "OutletTemperature")
        if T is None:
            raise _Unsupported("no outlet temperature given")
        p["T"] = _T(T)
    elif mode == "OutletVaporFraction":
        vf = _f(el, "OutletVaporFraction")
        if vf is None:
            raise _Unsupported("no outlet vapor fraction given")
        p["VF"] = vf
    else:
        p["T"] = imp.stored_outlet_T(o, f"calculation mode {mode}")
    _pressure_drop(p, el)
    return _unit(o, "Heater", p)


def _pump(imp: _Importer, o: _Obj) -> _Unit:
    el = o.el
    mode = _enum(el, "CalcMode", _PUMP_MODES, "Delta_P")
    p: dict[str, Any] = {}
    if mode == "Delta_P" and _f(el, "DeltaP") is not None:
        p["dP"] = _P(_f(el, "DeltaP"))
    elif mode == "OutletPressure" and _f(el, "Pout"):
        p["P"] = _P(_f(el, "Pout"))
    else:
        p["P"] = imp.stored_outlet_P(o, f"calculation mode {mode}")
    if (eff := _f(el, "Eficiencia")) is not None and 0 < eff <= 100:
        p["efficiency"] = round(eff / 100.0, 10)
    return _unit(o, "Pump", p)


def _machine(imp: _Importer, o: _Obj, expander: bool) -> _Unit:
    el = o.el
    modes = _EXPANDER_MODES if expander else _COMPRESSOR_MODES
    mode = _enum(el, "CalcMode", modes, "OutletPressure")
    p: dict[str, Any] = {}
    if mode == "OutletPressure" and _f(el, "POut"):
        p["P"] = _P(_f(el, "POut"))
    elif mode == "Delta_P" and _f(el, "DeltaP") is not None:
        p["dP"] = _P(abs(_f(el, "DeltaP")))
    else:
        p["P"] = imp.stored_outlet_P(o, f"calculation mode {mode}")
    path = _enum(el, "ProcessPath", ("Adiabatic", "Polytropic"), "Adiabatic")
    eff = _f(el, "AdiabaticEfficiency")
    if eff is None:
        eff = _f(el, "EficienciaAdiabatica")
    if path == "Polytropic":
        poly = _f(el, "PolytropicEfficiency")
        imp.warn(f"{o.label}: polytropic path imported as isentropic"
                 + (f"; the polytropic efficiency {_g(poly)} % is used as the isentropic one"
                    if poly is not None else ""))
        eff = poly if poly is not None else eff
    if eff is not None and 0 < eff <= 100:
        p["efficiency"] = round(eff / 100.0, 10)
    return _unit(o, "Turbine" if expander else "Compressor", p)


def _valve(imp: _Importer, o: _Obj) -> _Unit:
    el = o.el
    mode = _enum(el, "CalcMode", _VALVE_MODES, "DeltaP")
    p: dict[str, Any] = {}
    if mode == "DeltaP" and _f(el, "DeltaP") is not None:
        p["dP"] = _P(_f(el, "DeltaP"))
    elif mode == "OutletPressure" and _f(el, "OutletPressure"):
        p["P"] = _P(_f(el, "OutletPressure"))
    else:
        p["P"] = imp.stored_outlet_P(o, f"valve mode {mode}")
    return _unit(o, "Valve", p)


def _mixer(imp: _Importer, o: _Obj) -> _Unit:
    p: dict[str, Any] = {}
    how = _enum(o.el, "PressureCalculation", ("Average", "Maximum", "Minimum"), "Minimum")
    if how != "Minimum":  # Fluxsheet's mixer takes the lowest inlet pressure by default
        p["P"] = imp.stored_outlet_P(o, f"outlet pressure as the {how.lower()} of the inlets")
    return _unit(o, "Mixer", p, {i: "in" for i in range(6)})


def _tank(imp: _Importer, o: _Obj) -> _Unit:
    imp.warn(f"{o.label}: tank imported as a pass-through mixer (its volume is not used in "
             "a steady state)")
    return _unit(o, "Mixer", {}, {0: "in"})


def _splitter(imp: _Importer, o: _Obj) -> _Unit:
    el = o.el
    ratios = [_num(r.text) for r in _children(_child(el, "SplitRatios"), "SplitRatio")]
    if not ratios:
        ratios = [_num(r.text) for r in _children(_child(el, "Ratios"), "Item")]
    ratios = [r if r is not None else 0.0 for r in ratios]
    outs = imp.attached_outlets(o)
    mode = _enum(el, "OperationMode", ("SplitRatios", "StreamMassFlowSpec", "StreamMoleFlowSpec"),
                 "SplitRatios")
    if mode != "SplitRatios" or len(ratios) < max(outs, default=0) + 1:
        inlet = imp.inlet_state(o)
        flows = [(imp.outlet_state(o, i) or _State(None, None, None, {}, None, "")).flow
                 for i in outs]
        if inlet is None or not inlet.flow or any(f is None for f in flows):
            raise _Unsupported(f"split mode {mode} has no Fluxsheet counterpart and the file "
                               "stores no stream flows to derive split fractions from")
        ratios = [0.0] * (max(outs) + 1)
        for i, f in zip(outs, flows):
            ratios[i] = f / inlet.flow
        imp.warn(f"{o.label}: split mode {mode} imported as the split fractions DWSIM computed")
    if len(outs) <= 1:
        i = outs[0] if outs else 0
        return _unit(o, "Splitter", {"split": 1.0}, outlets={i: "out1"})
    if len(outs) > 3:
        raise _Unsupported("more than three outlets")
    r = [max(ratios[i], 0.0) for i in outs]
    total = sum(r) or 1.0
    r = [v / total for v in r]
    if len(outs) == 2:
        return _unit(o, "Splitter", {"split": round(r[0], 12)},
                     outlets={outs[0]: "out1", outs[1]: "out2"})
    # Three outlets: two splitters in a row.
    second = imp._unit_id(f"{o.id}_2")
    rest = r[1] + r[2]
    unit = _unit(o, "Splitter", {"split": round(r[0], 12)}, outlets={outs[0]: "out1"})
    unit.units.append({"id": second, "type": "Splitter",
                       "params": {"split": round(r[1] / rest, 12) if rest > 0 else 0.5}})
    unit.outlets[outs[1]] = f"{second}.out1"
    unit.outlets[outs[2]] = f"{second}.out2"
    unit.streams.append({"id": imp._stream_id(f"{o.id}_internal"), "from": f"{o.id}.out2",
                         "to": f"{second}.in"})
    imp.warn(f"{o.label}: three-way splitter imported as two splitters in a row ({o.id}, "
             f"{second})")
    return unit


def _vessel(imp: _Importer, o: _Obj) -> _Unit:
    el = o.el
    p: dict[str, Any] = {}
    if _b(el, "OverrideT") and _f(el, "FlashTemperature") is not None:
        p["T"] = _T(_f(el, "FlashTemperature"))
    if _b(el, "OverrideP") and _f(el, "FlashPressure"):
        p["P"] = _P(_f(el, "FlashPressure"))
    mode = _enum(el, "CalculationMode", ("Adiabatic", "Legacy", "HeatingCoolingIsothermic",
                                         "HeatingCoolingIsobaric"), "Legacy")
    if mode.startswith("HeatingCooling") and "T" not in p:
        p["T"] = imp.stored_outlet_T(o, f"vessel mode {mode}")
    if _enum(el, "PressureCalculation", ("Average", "Maximum", "Minimum"), "Minimum") != \
            "Minimum" and "P" not in p:
        p["P"] = imp.stored_outlet_P(o, "vessel pressure from the average or maximum inlet")
    inlets = {i: "in" for i in range(6)}
    if 2 in imp.attached_outlets(o):  # a second liquid outlet: three-phase separator
        imp.warn(f"{o.label}: three-phase vessel imported as a Decanter (liquid 1 -> light, "
                 "liquid 2 -> heavy)")
        return _unit(o, "Decanter", p, inlets, {0: "vapor", 1: "light", 2: "heavy"})
    return _unit(o, "Flash", p, inlets, {0: "vapor", 1: "liquid"})


def _heat_exchanger(imp: _Importer, o: _Obj) -> _Unit:
    el = o.el
    s0, s1 = imp.inlet_state(o, 0), imp.inlet_state(o, 1)
    if s0 and s1 and s0.T is not None and s1.T is not None:
        hot = 0 if s0.T >= s1.T else 1  # DWSIM assigns the sides by temperature, too
    else:
        hot = 0
        imp.warn(f"{o.label}: the file stores no inlet temperatures; inlet 1 is taken as the "
                 "hot side")
    cold = 1 - hot
    p: dict[str, Any] = {}
    mode = _enum(el, "CalculationMode", _HX_MODES, "CalcBothTemp_UA")
    t_hot, t_cold = _f(el, "HotSideOutletTemperature"), _f(el, "ColdSideOutletTemperature")
    q = _f(el, "Q")
    if mode == "CalcTempHotOut" and t_cold is not None:
        p["T_cold_out"] = _T(t_cold)
    elif mode == "CalcTempColdOut" and t_hot is not None:
        p["T_hot_out"] = _T(t_hot)
    elif mode == "CalcBothTemp" and q is not None:
        p["duty"] = _kW(q)
    elif mode == "CalcBothTemp_UA" and _f(el, "OverallCoefficient") and _f(el, "Area"):
        p["UA"] = _f(el, "OverallCoefficient") * _f(el, "Area")
    elif mode == "PinchPoint" and _f(el, "MITA"):
        p["min_approach"] = f"{_g(_f(el, 'MITA'))} K"
    elif mode == "CalcArea" and t_hot is not None:
        p["T_hot_out"] = _T(t_hot)
        imp.warn(f"{o.label}: area calculation imported with the hot outlet temperature as the "
                 "spec (Fluxsheet takes one spec; the cold outlet follows)")
    elif q is not None:
        p["duty"] = _kW(q)
        imp.warn(f"{o.label}: exchanger mode {mode} has no Fluxsheet counterpart; the duty "
                 f"DWSIM computed ({_kW(q)}) is used as the spec")
    else:
        raise _Unsupported(f"exchanger mode {mode} without a stored duty")
    dp = {hot: _f(el, "HotSidePressureDrop"), cold: _f(el, "ColdSidePressureDrop")}
    if dp[hot]:
        p["dP_hot"] = _P(dp[hot])
    if dp[cold]:
        p["dP_cold"] = _P(dp[cold])
    if (_s(el, "FlowDir") or "").lower() == "cocurrent":
        p["arrangement"] = "cocurrent"
    side = {hot: "hot", cold: "cold"}
    return _unit(o, "HeatExchanger", p, {i: f"{side[i]}_in" for i in (0, 1)},
                 {i: f"{side[i]}_out" for i in (0, 1)})


def _shortcut_column(imp: _Importer, o: _Obj) -> _Unit:
    el = o.el
    lk, hk = imp.comp(_s(el, "m_lightkey")), imp.comp(_s(el, "m_heavykey"))
    if lk is None or hk is None:
        raise _Unsupported("its light or heavy key is not set")
    p: dict[str, Any] = {"light_key": lk, "heavy_key": hk}
    if (rr := _f(el, "m_refluxratio")) is not None:
        p["reflux_ratio"] = rr
    if (pc := _f(el, "m_condenserpressure")) is not None:
        p["P"] = _P(pc)
        pb = _f(el, "m_boilerpressure")
        if pb is not None and abs(pb - pc) > 1e-3 * pc:
            imp.warn(f"{o.label}: reboiler pressure {_P(pb)} differs from the condenser "
                     "pressure; Fluxsheet's shortcut column runs at the condenser pressure")
    if _enum(el, "condtype", ("TotalCond", "PartialCond"), "TotalCond") == "PartialCond":
        p["condenser"] = "partial"
    # DWSIM specifies the light key in the bottoms and the heavy key in the distillate (mole
    # fractions); Fluxsheet the key recoveries to the distillate.
    d, b, f = imp.outlet_state(o, 0), imp.outlet_state(o, 1), imp.inlet_state(o, 0)
    rec = None
    if d and b and d.flow and b.flow:
        def recovery(c):
            nd, nb = d.flow * d.x.get(c, 0.0), b.flow * b.x.get(c, 0.0)
            return nd / (nd + nb) if nd + nb > 0 else None
        rec = recovery(lk), recovery(hk)
        how = "from the product streams stored in the file"
    if (rec is None or None in rec) and f and len([v for v in f.x.values() if v > 0]) == 2:
        z = f.x.get(lk, 0.0)
        xlb = _f(el, "m_lightkeymolarfrac") or 0.01
        xhd = _f(el, "m_heavykeymolarfrac") or 0.01
        xld = 1.0 - xhd
        if 0 < z < 1 and xld > xlb:
            dfrac = (z - xlb) / (xld - xlb)
            rec = dfrac * xld / z, dfrac * xhd / (1.0 - z)
            how = "from the key mole fractions and the feed composition"
    if rec is None or None in rec:
        imp.warn(f"{o.label}: key specs (mole fractions in the products) could not be turned "
                 "into recoveries; Fluxsheet's defaults (0.99 and 0.01) are used")
    else:
        p["lk_recovery"] = round(min(max(rec[0], 1e-6), 1 - 1e-6), 8)
        p["hk_recovery"] = round(min(max(rec[1], 1e-6), 1 - 1e-6), 8)
        imp.warn(f"{o.label}: key recoveries derived {how}")
    return _unit(o, "ShortcutColumn", p, {0: "in"}, {0: "distillate", 1: "bottoms"})


def _component_separator(imp: _Importer, o: _Obj) -> _Unit:
    el = o.el
    to_spec = int(_f(el, "SpecifiedStreamIndex") or 0)
    inlet = imp.inlet_state(o)
    fractions: dict[str, float] = {}
    for spec in _children(_child(el, "SeparationSpecs"), "SeparationSpec"):
        c = imp.comp(spec.get("CompID"))
        v = _num(spec.get("SpecValue"))
        kind = spec.get("SepSpec") or "PercentInletMassFlow"
        if c is None or v is None:
            continue
        if kind.startswith("Percent"):
            frac = v / 100.0
        else:
            unit = (spec.get("SpecUnit") or "").strip()
            quantity = "mass_flow" if kind == "MassFlow" else "molar_flow"
            try:
                value = to_si(f"{v} {unit}" if unit else v, quantity)
            except UnitError:
                value = v
            n_in = inlet.flow * inlet.x.get(c, 0.0) if inlet and inlet.flow else None
            if not n_in:
                raise _Unsupported(f"the absolute flow spec of {c} needs the inlet flow, which "
                                   "the file does not store")
            frac = (value / imp.mw[c] if kind == "MassFlow" else value) / n_in
            imp.warn(f"{o.label}: the flow spec of {c} is imported as the fraction "
                     f"{_g(frac)} of the stored inlet flow")
        fractions[c] = min(max(frac, 0.0), 1.0)
    missing = [c for c in imp.components if c not in fractions]
    if missing and fractions:
        imp.warn(f"{o.label}: no separation spec for {', '.join(missing)}; "
                 f"{'it leaves' if len(missing) == 1 else 'they leave'} entirely with outlet "
                 f"{2 if to_spec == 0 else 1}")
    if to_spec != 0:  # the specs refer to the second outlet
        fractions = {c: 1.0 - fractions.get(c, 0.0) for c in imp.components}
    if not fractions:
        raise _Unsupported("it has no separation specs")
    return _unit(o, "ComponentSplitter", {"fractions": fractions}, {0: "in"},
                 {0: "top", 1: "bottom"})


def _reactor_mode(imp: _Importer, o: _Obj, p: dict) -> None:
    el = o.el
    mode = _enum(el, "ReactorOperationMode", _REACTOR_MODES, "Adiabatic")
    if mode == "OutletTemperature" and _f(el, "OutletTemperature") is not None:
        p["T"] = _T(_f(el, "OutletTemperature"))
    elif mode == "Isothermic":
        st = imp.inlet_state(o)
        if st is None or st.T is None:
            raise _Unsupported("isothermal operation needs the inlet temperature, which the "
                               "file does not store")
        p["T"] = _T(st.T)
        imp.warn(f"{o.label}: isothermal reactor imported at its stored inlet temperature "
                 f"({_T(st.T)})")
    elif mode == "NonIsothermalNonAdiabatic":
        dq = _f(el, "DeltaQ")
        if dq is None:
            raise _Unsupported("mode NonIsothermalNonAdiabatic without a stored duty")
        p["duty"] = _kW(dq)
        imp.warn(f"{o.label}: reactor mode with an energy stream imported with the duty DWSIM "
                 f"computed ({_kW(dq)})")
    _pressure_drop(p, el)


def _reactor_outlets(imp: _Importer, o: _Obj, unit: _Unit) -> _Unit:
    """DWSIM reactors have a vapor (0) and a liquid (1) outlet; Fluxsheet's one outlet, so
    with both connected an adiabatic flash drum is added behind the reactor."""
    outs = imp.attached_outlets(o)
    if len(outs) <= 1:
        unit.outlets = {i: f"{o.id}.out" for i in (outs or [0])}
        return unit
    sep = imp._unit_id(f"{o.id}_SEP")
    unit.units.append({"id": sep, "type": "Flash"})
    unit.streams.append({"id": imp._stream_id(f"{o.id}_OUT"), "from": f"{o.id}.out",
                         "to": f"{sep}.in"})
    unit.outlets = {0: f"{sep}.vapor", 1: f"{sep}.liquid"}
    imp.warn(f"{o.label}: its vapor and liquid outlets come from an added flash drum {sep}")
    return unit


def _reaction_set(imp: _Importer, set_id: str | None) -> list[tuple[int, Element]]:
    reactions = {(_s(r, "ID") or ""): r for r in _children(_child(imp.root, "Reactions"),
                                                           "Reaction")}
    for rs in _children(_child(imp.root, "ReactionSets"), "ReactionSet"):
        if _s(rs, "ID") == (set_id or "DefaultSet") or _s(rs, "Name") == set_id:
            out = []
            for item in _children(_child(rs, "Reactions"), "Reaction"):
                rid = item.get("ReactionID") or item.get("Key")
                if rid in reactions and (item.get("IsActive") or "true").lower() != "false":
                    out.append((int(_num(item.get("Rank")) or 0), reactions[rid]))
            return sorted(out, key=lambda t: t[0])
    return []


def _conversion_reactor(imp: _Importer, o: _Obj) -> _Unit:
    el = o.el
    stored = {r.get("ID"): _num(r.text)
              for r in _children(_child(el, "ReactionConversions"), "Reaction")}
    ids, conversions, keys, ranks = [], [], [], []
    for rank, rx in _reaction_set(imp, _s(el, "ReactionSetID")):
        if _s(rx, "ReactionType") not in (None, "Conversion"):
            continue
        rid = _s(rx, "ID") or ""
        label = _s(rx, "Name") or rid
        stoich, base = {}, imp.comp(_s(rx, "BaseReactant"))
        for c in _children(_child(rx, "Compounds"), "Compound"):
            name, nu = imp.comp(c.get("Name")), _num(c.get("StoichCoeff"))
            if name is None or nu is None:
                raise _Unsupported(f"reaction {label} names a compound that is not imported")
            if nu:
                stoich[name] = stoich.get(name, 0.0) + nu
            if base is None and _attr_bool(c, "IsBaseReactant"):
                base = name
        if not stoich:
            raise _Unsupported(f"reaction {label} has no stoichiometry")
        x = _num(_s(rx, "Expression"))
        if x is not None:
            x /= 100.0  # DWSIM's conversion is in percent
        elif stored.get(rid) is not None:
            x = stored[rid] if stored[rid] <= 1.0 else stored[rid] / 100.0
            imp.warn(f"{o.label}: the conversion of {label} is an expression "
                     f"({_s(rx, 'Expression')}); the value DWSIM computed ({_g(x)}) is used")
        else:
            raise _Unsupported(f"the conversion of reaction {label} is not a number")
        if rid not in imp.reaction_ids:
            fid = re.sub(r"\s+", "_", label.strip()) or "R"
            existing = {r["id"] for r in imp.reactions}
            base_id, n = fid, 2
            while fid in existing:
                fid, n = f"{base_id}_{n}", n + 1
            imp.reaction_ids[rid] = fid
            imp.reactions.append({"id": fid, "stoichiometry": stoich})
        ids.append(imp.reaction_ids[rid])
        conversions.append(_g(min(max(x, 0.0), 1.0)))
        keys.append(base or next(c for c, nu in stoich.items() if nu < 0))
        ranks.append(rank)
    if not ids:
        raise _Unsupported("its reaction set holds no active conversion reactions")
    p: dict[str, Any] = {"reactions": ids, "conversions": conversions, "keys": keys}
    if len(ids) > 1:
        if len(set(ranks)) == 1:
            p["mode"] = "parallel"
        elif len(set(ranks)) != len(ranks):
            imp.warn(f"{o.label}: reactions of equal rank run in parallel in DWSIM; Fluxsheet "
                     "runs them one after another in rank order")
    _reactor_mode(imp, o, p)
    return _reactor_outlets(imp, o, _unit(o, "StoichiometricReactor", p))


def _gibbs_reactor(imp: _Importer, o: _Obj) -> _Unit:
    el = o.el
    p: dict[str, Any] = {}
    products = [imp.comp(c.get("ID")) for c in _children(_child(el, "Compounds"), "Compound")]
    products = [c for c in products if c]
    if products and set(products) != set(imp.components):
        p["products"] = products
    phase = _enum(el, "ReactivePhaseBehavior", ("CalculateEquilibrium", "Vapor", "Liquid",
                                                "Solid"), "CalculateEquilibrium")
    if phase != "Vapor":
        imp.warn(f"{o.label}: Fluxsheet's Gibbs reactor treats the reacting mixture as an "
                 f"ideal gas (DWSIM: {phase})")
    _reactor_mode(imp, o, p)
    return _reactor_outlets(imp, o, _unit(o, "GibbsReactor", p))


_CONVERTERS = {
    "Heater": lambda imp, o: _heater(imp, o, cooler=False),
    "Cooler": lambda imp, o: _heater(imp, o, cooler=True),
    "Pump": _pump,
    "Compressor": lambda imp, o: _machine(imp, o, expander=False),
    "Expander": lambda imp, o: _machine(imp, o, expander=True),
    "Valve": _valve,
    "Mixer": _mixer,
    "Splitter": _splitter,
    "Tank": _tank,
    "Vessel": _vessel,
    "HeatExchanger": _heat_exchanger,
    "ShortcutColumn": _shortcut_column,
    "ComponentSeparator": _component_separator,
    "Reactor_Conversion": _conversion_reactor,
    "Reactor_Gibbs": _gibbs_reactor,
}
_KNOWN_CLASSES = {*_CONVERTERS, *_UNSUPPORTED, *_LOGICAL, "MaterialStream", "EnergyStream",
                  "Recycle"}


def import_dwsim(source: str | Path | bytes, *, name: str | None = None
                 ) -> tuple[dict, list[str]]:
    """Read a DWSIM simulation (path to, or bytes of, a .dwxml or .dwxmz file) and return
    the Fluxsheet flowsheet as a JSON dict plus the list of warnings about what was left
    out or approximated. Raises :class:`DwsimImportError` for files that cannot be read."""
    data, filename = _read_source(source)
    notes: list[str] = []
    if data[:4] == b"PK\x03\x04":
        data, notes = _unzip(data)
    elif data[:2] == b"PK":
        raise DwsimImportError("not a valid .dwxmz file: the zip archive is damaged")
    root = _parse(data)
    title = name or (Path(filename).stem if filename else "Imported DWSIM flowsheet")
    return _Importer(root, title, notes).run()
