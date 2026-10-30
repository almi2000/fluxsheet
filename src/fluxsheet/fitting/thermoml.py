"""Read binary vapor-liquid equilibrium data from ThermoML files.

ThermoML (IUPAC) is the XML format of the NIST TRC ThermoML Archive. Downloaded files are
untrusted input: they are parsed with ``defusedxml`` (no DTDs, entities or external
references), size-limited, and every number is checked before use.

Only what a binary VLE regression needs is extracted: data sets with exactly two
components, a liquid and a gas phase, and temperature, pressure and liquid mole fraction
for every point (vapor mole fraction where measured). Temperature and pressure may be a
constraint (isothermal, isobaric), a variable or a measured property.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from pathlib import Path

from defusedxml import ElementTree

MAX_BYTES = 50_000_000  # the largest archive files are a few MB

_UNIT_FACTORS = {"K": 1.0, "Pa": 1.0, "kPa": 1e3, "MPa": 1e6, "bar": 1e5}


class ThermoMLError(ValueError):
    """The file is not a ThermoML document this parser can read."""


@dataclass(frozen=True)
class ThermoMLCompound:
    org_num: int
    names: tuple[str, ...]
    formula: str | None = None
    inchi: str | None = None
    inchikey: str | None = None

    @property
    def name(self) -> str:
        return self.names[0] if self.names else f"compound {self.org_num}"


@dataclass(frozen=True)
class VLEPoint:
    T: float  # K
    P: float  # Pa
    x1: float  # liquid mole fraction of the first component
    y1: float | None = None  # vapor mole fraction of the first component, if measured


@dataclass(frozen=True)
class VLEDataset:
    """One binary VLE data set of an article, in the component order of ``compounds``."""

    doi: str
    number: int  # nPureOrMixtureDataNumber within the file
    compounds: tuple[ThermoMLCompound, ThermoMLCompound]
    kind: str  # "isothermal", "isobaric" or "other" (both T and P vary)
    points: tuple[VLEPoint, ...]

    @property
    def has_vapor(self) -> bool:
        return all(p.y1 is not None for p in self.points)

    def swapped(self) -> VLEDataset:
        return replace(self, compounds=self.compounds[::-1], points=tuple(
            VLEPoint(p.T, p.P, 1.0 - p.x1, None if p.y1 is None else 1.0 - p.y1)
            for p in self.points))

    @property
    def label(self) -> str:
        if self.kind == "isothermal":
            return f"isothermal {self.points[0].T:.2f} K"
        if self.kind == "isobaric":
            return f"isobaric {self.points[0].P / 1e3:.4g} kPa"
        return "T and P varying"


@dataclass(frozen=True)
class ThermoMLReport:
    doi: str
    title: str
    citation: str  # "Author et al., J. Chem. Eng. Data 53 (2008) 820-825"
    compounds: dict[int, ThermoMLCompound]
    datasets: list[VLEDataset]


def _local(tag) -> str:
    return tag.rsplit("}", 1)[-1] if isinstance(tag, str) else ""


def _children(el, name: str) -> list:
    return [c for c in el if _local(c.tag) == name]


def _child(el, *path: str):
    for name in path:
        if el is None:
            return None
        el = next((c for c in el if _local(c.tag) == name), None)
    return el


def _text(el, *path: str) -> str | None:
    node = _child(el, *path)
    if node is None or node.text is None:
        return None
    return node.text.strip() or None


def _number(text: str | None) -> float | None:
    if text is None:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def _org_num(el) -> int | None:
    value = _number(_text(el, "RegNum", "nOrgNum"))
    return int(value) if value is not None else None


def _unit_factor(label: str) -> float | None:
    """'Pressure, kPa' -> 1000.0; None for units this parser does not know."""
    unit = label.rsplit(",", 1)[-1].strip() if "," in label else ""
    return _UNIT_FACTORS.get(unit)


def _column(kind_el, phase: str | None, org: int | None) -> tuple[str, float] | None:
    """Map a variable or constraint type element to (quantity, SI factor)."""
    tag, text = _local(kind_el.tag), (kind_el.text or "").strip()
    if tag in ("eTemperature", "ePressure"):
        factor = _unit_factor(text)
        if factor is None or (tag == "eTemperature") != text.startswith("Temperature"):
            return None
        return ("T" if tag == "eTemperature" else "P"), factor
    if tag == "eComponentComposition" and text == "Mole fraction" and org is not None:
        if phase == "Liquid":
            return f"x{org}", 1.0
        if phase == "Gas":
            return f"y{org}", 1.0
    return None


def _property_column(prop) -> tuple[str, float] | None:
    group = _child(prop, "Property-MethodID", "PropertyGroup")
    if group is None or not len(group):
        return None
    name = _text(group[0], "ePropName") or ""
    phase = _text(prop, "PropPhaseID", "ePropPhase")
    org = _org_num(_child(prop, "Property-MethodID"))
    if name == "Mole fraction" and org is not None:
        return {"Liquid": (f"x{org}", 1.0), "Gas": (f"y{org}", 1.0)}.get(phase)
    if name.startswith(("Vapor or sublimation pressure", "Boiling temperature")):
        factor = _unit_factor(name)
        return None if factor is None else ("P" if name.startswith("Vapor") else "T", factor)
    return None


def _parse_compound(el) -> ThermoMLCompound | None:
    org = _org_num(el)
    if org is None:
        return None
    names = tuple(c.text.strip() for c in _children(el, "sCommonName") if c.text and c.text.strip())
    return ThermoMLCompound(org, names, _text(el, "sFormulaMolec"), _text(el, "sStandardInChI"),
                            _text(el, "sStandardInChIKey"))


def _citation(cit) -> tuple[str, str, str]:
    doi = _text(cit, "sDOI") or ""
    title = _text(cit, "sTitle") or ""
    authors = [c.text.split(",")[0].strip() for c in _children(cit, "sAuthor") if c.text]
    who = (authors[0] + (" et al." if len(authors) > 2 else f" and {authors[1]}"
                         if len(authors) == 2 else "")) if authors else "unknown authors"
    where = " ".join(v for v in (_text(cit, "sPubName"), _text(cit, "sVol"),
                                 f"({_text(cit, 'yrPubYr')})" if _text(cit, "yrPubYr") else None,
                                 _text(cit, "sPage")) if v)
    return doi, title, f"{who}, {where}".strip(", ")


def _parse_table(el, compounds: dict[int, ThermoMLCompound]) -> _Table | None:
    orgs = [_org_num(c) for c in _children(el, "Component")]
    if len(orgs) != 2 or None in orgs or len(set(orgs)) != 2:
        return None
    if any(o not in compounds for o in orgs):
        return None
    phases = {_text(p, "ePhase") for p in _children(el, "PhaseID")}
    if phases != {"Liquid", "Gas"}:
        return None

    columns: dict[tuple[str, int], tuple[str, float]] = {}
    for prop in _children(el, "Property"):
        number = _number(_text(prop, "nPropNumber"))
        col = _property_column(prop)
        if number is None:
            return None
        if col is not None:
            columns[("prop", int(number))] = col
    for var in _children(el, "Variable"):
        number = _number(_text(var, "nVarNumber"))
        kind = _child(var, "VariableID", "VariableType")
        if number is None or kind is None or not len(kind):
            return None
        col = _column(kind[0], _text(var, "VarPhaseID", "eVarPhase"),
                      _org_num(_child(var, "VariableID")))
        if col is None:  # a variable we cannot interpret makes the rows ambiguous
            return None
        columns[("var", int(number))] = col
    fixed: dict[str, float] = {}
    for con in _children(el, "Constraint"):
        kind = _child(con, "ConstraintID", "ConstraintType")
        value = _number(_text(con, "nConstraintValue"))
        if kind is None or not len(kind) or value is None:
            return None
        col = _column(kind[0], _text(con, "ConstraintPhaseID", "eConstraintPhase"),
                      _org_num(_child(con, "ConstraintID")))
        if col is None:
            return None
        fixed[col[0]] = value * col[1]

    first, second = sorted(orgs)  # canonical order, so split tables of one set line up
    rows = []
    for row in _children(el, "NumValues"):
        values = dict(fixed)
        for tag, key, num, val in (("VariableValue", "var", "nVarNumber", "nVarValue"),
                                   ("PropertyValue", "prop", "nPropNumber", "nPropValue")):
            for item in _children(row, tag):
                number, value = _number(_text(item, num)), _number(_text(item, val))
                col = columns.get((key, int(number))) if number is not None else None
                if col is not None and value is not None:
                    values[col[0]] = value * col[1]
        rows.append(_row(values, first, second))
    if not rows:
        return None
    return _Table(int(_number(_text(el, "nPureOrMixtureDataNumber")) or 0), (first, second),
                  tuple(sorted(fixed.items())), rows)


@dataclass
class _Table:
    """One PureOrMixtureData block. ThermoML often stores a T-P-x-y set as two blocks, the
    boiling temperature (or pressure) against x and the vapor composition against x."""

    number: int
    orgs: tuple[int, int]
    fixed: tuple
    rows: list[dict]  # keys T, P, x1, y1 (first component = lower nOrgNum), each optional

    def has(self, key: str) -> bool:
        return all(r.get(key) is not None for r in self.rows)


def _row(values: dict, first: int, second: int) -> dict:
    return {"T": values.get("T"), "P": values.get("P"),
            "x1": _fraction(values, "x", first, second), "y1": _fraction(values, "y", first, second)}


def _same_rows(a: _Table, b: _Table) -> bool:
    """Whether two blocks describe the same experimental points (same x, and same T or P
    wherever both give one)."""
    if a.orgs != b.orgs or a.fixed != b.fixed or len(a.rows) != len(b.rows):
        return False
    for ra, rb in zip(a.rows, b.rows):
        for key in ("x1", "T", "P"):
            va, vb = ra.get(key), rb.get(key)
            if key == "x1" and (va is None or vb is None):
                return False
            if va is not None and vb is not None and abs(va - vb) > 1e-9 * max(1.0, abs(va)):
                return False
    return True


def _merge(tables: list[_Table]) -> list[_Table]:
    """Join vapor-composition blocks to the block with the matching T/P, row by row."""
    out, used = [], set()
    for i, t in enumerate(tables):
        if i in used or not (t.has("T") and t.has("P") and t.has("x1")):
            continue
        rows = [dict(r) for r in t.rows]
        if not t.has("y1"):
            for j, other in enumerate(tables):
                if j != i and j not in used and other.has("y1") and _same_rows(t, other):
                    for r, o in zip(rows, other.rows):
                        r["y1"] = o["y1"]
                    used.add(j)
                    break
        used.add(i)
        out.append(_Table(t.number, t.orgs, t.fixed, rows))
    return out


def _datasets(table: _Table, doi: str, compounds: dict[int, ThermoMLCompound]) -> list[VLEDataset]:
    points = [p for r in table.rows if (p := _point(r)) is not None]
    if any(p.y1 is None for p in points) and any(p.y1 is not None for p in points):
        points = [p for p in points if p.y1 is not None]  # keep the complete T-P-x-y rows
    pair = (compounds[table.orgs[0]], compounds[table.orgs[1]])
    kind = _kind(points, dict(table.fixed)) if points else None
    if kind == "other":  # several isotherms (or isobars) in one table: one data set each
        for key, label in ((lambda p: round(p.T, 2), "isothermal"),
                           (lambda p: round(p.P, -1), "isobaric")):
            groups: dict[float, list[VLEPoint]] = {}
            for p in points:
                groups.setdefault(key(p), []).append(p)
            if len(groups) < len(points) / 3:
                return [VLEDataset(doi, table.number, pair, label, tuple(g))
                        for g in groups.values()]
    return [VLEDataset(doi, table.number, pair, kind, tuple(points))] if points else []


def _fraction(values: dict, phase: str, first: int, second: int) -> float | None:
    if f"{phase}{first}" in values:
        return values[f"{phase}{first}"]
    if f"{phase}{second}" in values:
        return 1.0 - values[f"{phase}{second}"]
    return None


def _point(row: dict) -> VLEPoint | None:
    T, P, x1, y1 = row["T"], row["P"], row["x1"], row["y1"]
    if T is None or P is None or x1 is None:
        return None
    if not (50.0 < T < 2000.0 and 0.0 < P < 1e9 and 0.0 <= x1 <= 1.0):
        return None
    if y1 is not None and not 0.0 <= y1 <= 1.0:
        y1 = None
    return VLEPoint(T, P, x1, y1)


def _kind(points: list[VLEPoint], fixed: dict) -> str:
    Ts, Ps = [p.T for p in points], [p.P for p in points]
    if "T" in fixed or max(Ts) - min(Ts) < 0.01:
        return "isothermal"
    if "P" in fixed or max(Ps) - min(Ps) < 1e-3 * max(Ps):
        return "isobaric"
    return "other"


def parse_thermoml(source: bytes | str | Path) -> ThermoMLReport:
    """Parse a ThermoML document (bytes, XML text or a path) into binary VLE data sets."""
    if isinstance(source, Path):
        if source.stat().st_size > MAX_BYTES:
            raise ThermoMLError(f"{source} is larger than {MAX_BYTES} bytes")
        source = source.read_bytes()
    if isinstance(source, str):
        source = source.encode("utf-8")
    if len(source) > MAX_BYTES:
        raise ThermoMLError(f"document is larger than {MAX_BYTES} bytes")
    try:
        root = ElementTree.fromstring(source, forbid_dtd=True)
    except Exception as exc:  # ParseError and defusedxml's DTD/entity refusals
        raise ThermoMLError(f"not a well-formed, safe XML document: {exc}") from exc
    if _local(root.tag) != "DataReport":
        raise ThermoMLError("not a ThermoML document (root element is not DataReport)")
    doi, title, citation = _citation(_child(root, "Citation"))
    compounds = {}
    for el in _children(root, "Compound"):
        compound = _parse_compound(el)
        if compound is not None:
            compounds[compound.org_num] = compound
    tables = [t for el in _children(root, "PureOrMixtureData")
              if (t := _parse_table(el, compounds)) is not None]
    datasets = [ds for t in _merge(tables) for ds in _datasets(t, doi, compounds)]
    return ThermoMLReport(doi, title, citation, compounds, datasets)
