"""User-defined compounds: name, formula and groups, with estimates for what is not given.

A flowsheet can define compounds the database lacks (or redefine one it has)::

    "custom_compounds": [
      {"name": "my-ester", "formula": "C5H10O2", "smiles": "CCCC(=O)OC",
       "joback_groups": {"-CH3": 2, "-CH2-": 2, "-COO- (ester)": 1},
       "unifac_groups": {"CH3": 2, "CH2": 1, "CH2COO": 1},
       "constants": {"Tb": "102.7 degC", "Pc": "34.7 bar"},
       "antoine": {"A": 9.321, "B": 1389.0, "C": -54.0}}
    ]

Missing constants are estimated: Tb, Tc, Pc, Vc, the ideal-gas formation enthalpy and Gibbs
energy and the ideal-gas heat capacity by Joback's group contributions (Joback and Reid
1987), the acentric factor from Tb, Tc and Pc by the Lee-Kesler vapor-pressure equation, the
vapor pressure by Lee-Kesler (or Ambrose-Walton) from Tc, Pc and omega unless Antoine
constants are given (``log10(P/Pa) = A - B / (T/K + C)``), the liquid volume by the COSTALD
method. Joback groups come from the definition, or from the SMILES when RDKit is
installed (``pip install fluxsheet[estimation]``). Every estimate is recorded, so results can
say which numbers are not measured.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from functools import lru_cache

from chemicals.acentric import LK_omega
from chemicals.elements import molecular_weight, similarity_variable, simple_formula_parser
from thermo import (
    ChemicalConstantsPackage,
    HeatCapacityGas,
    PropertyCorrelationsPackage,
    VaporPressure,
    VolumeLiquid,
)
from thermo.group_contribution.joback import JOBACK_GROUPS, Joback
from thermo.unifac import DOUFSG, UFSG

from .errors import SimulationError
from .quantities import to_si

PREFIX = "custom:"
R = 8.314462618

# Known constants a definition may give: key -> (label, quantity for unit conversion or None
# for plain SI numbers).
CONSTANTS: dict[str, tuple[str, str | None]] = {
    "MW": ("molecular weight [g/mol]", None),
    "Tb": ("normal boiling point", "temperature"),
    "Tm": ("melting point", "temperature"),
    "Tc": ("critical temperature", "temperature"),
    "Pc": ("critical pressure", "pressure"),
    "Vc": ("critical volume [m3/mol]", None),
    "omega": ("acentric factor", None),
    "Hf": ("ideal-gas enthalpy of formation at 25 °C", "molar_energy"),
    "Gf": ("ideal-gas Gibbs energy of formation at 25 °C", "molar_energy"),
    "dipole": ("dipole moment [debye]", None),
    "r": ("UNIQUAC volume parameter r", None),
    "q": ("UNIQUAC surface parameter q", None),
}
VAPOR_PRESSURE_METHODS = {"lee-kesler": "LEE_KESLER_PSAT", "ambrose-walton": "AMBROSE_WALTON"}


@dataclass(frozen=True)
class CustomCompound:
    """A definition as written in the flowsheet, constants already in SI."""

    name: str
    formula: str
    smiles: str = ""
    joback_groups: tuple[tuple[str, int], ...] = ()
    unifac_groups: tuple[tuple[str, int], ...] = ()
    constants: tuple[tuple[str, float], ...] = ()
    vapor_pressure: str = "lee-kesler"
    # log10(P/Pa) = A - B / (T/K + C), with its range in K: (A, B, C, Tmin, Tmax).
    antoine: tuple[float, float, float, float | None, float | None] | None = None

    @classmethod
    def from_dict(cls, d: dict) -> CustomCompound:
        """From the flowsheet JSON; unit strings are converted, everything is checked."""
        name = str(d.get("name") or "").strip()
        if not name:
            raise SimulationError("a custom compound has no name", hint="give it a unique name")
        where = f"custom compound {name!r}"
        formula = str(d.get("formula") or "").strip()
        if not formula:
            raise SimulationError(f"{where} has no formula",
                                  hint="give the molecular formula, e.g. 'C2H6O'")
        constants = []
        for key, value in (d.get("constants") or {}).items():
            if key not in CONSTANTS:
                raise SimulationError(f"{where}: unknown constant {key!r}",
                                      hint=f"known constants: {', '.join(CONSTANTS)}")
            if value is None or value == "":
                continue
            quantity = CONSTANTS[key][1]
            try:
                si = to_si(value, quantity) if quantity else float(value)
            except (ValueError, TypeError) as exc:
                raise SimulationError(f"{where}: {key} = {value!r} is not a valid value",
                                      hint=str(exc)) from exc
            constants.append((key, si))
        vp = str(d.get("vapor_pressure") or "lee-kesler").lower()
        if vp not in VAPOR_PRESSURE_METHODS:
            raise SimulationError(f"{where}: unknown vapor pressure method {vp!r}",
                                  hint=f"use one of: {', '.join(VAPOR_PRESSURE_METHODS)}")

        def groups(key):
            out = []
            for g, n in (d.get(key) or {}).items():
                try:
                    count = int(n)
                except (TypeError, ValueError):
                    count = -1
                if count < 0 or count != n:
                    raise SimulationError(f"{where}: {key} count for {g!r} must be a "
                                          f"non-negative whole number, got {n!r}")
                if count:
                    out.append((str(g), count))
            return tuple(sorted(out))

        antoine = None
        if d.get("antoine"):
            a = d["antoine"]
            try:
                antoine = (float(a["A"]), float(a["B"]), float(a["C"]),
                           *(to_si(a[k], "temperature") if a.get(k) not in (None, "") else None
                             for k in ("Tmin", "Tmax")))
            except (KeyError, TypeError, ValueError) as exc:
                raise SimulationError(
                    f"{where}: the Antoine constants need numbers A, B and C",
                    hint="log10(P/Pa) = A - B / (T/K + C); Tmin and Tmax are optional") from exc
        return cls(name=name, formula=formula, smiles=str(d.get("smiles") or "").strip(),
                   joback_groups=groups("joback_groups"), unifac_groups=groups("unifac_groups"),
                   constants=tuple(sorted(constants)), vapor_pressure=vp, antoine=antoine)

    @property
    def id(self) -> str:
        """Stands in for a CAS number; changes with the definition, so caches stay valid."""
        text = json.dumps([self.name, self.formula, self.smiles, self.joback_groups,
                           self.unifac_groups, self.constants, self.vapor_pressure, self.antoine])
        return PREFIX + hashlib.sha1(text.encode()).hexdigest()[:12]


def is_custom(identifier: str) -> bool:
    return identifier.startswith(PREFIX)


# Definitions by id, so that caches keyed by ids (like CAS numbers) can rebuild packages.
_REGISTRY: dict[str, CustomCompound] = {}


def register(compound: CustomCompound) -> str:
    _REGISTRY[compound.id] = compound
    return compound.id


def definition(identifier: str) -> CustomCompound:
    return _REGISTRY[identifier]


@dataclass
class Estimate:
    """Constants of a custom compound, and which of them were estimated, and how."""

    compound: CustomCompound
    values: dict[str, float | None]
    atoms: dict[str, int]
    joback: dict[int, int]
    unifac: dict[int, int]
    cp_coeffs: list[float] | None  # Joback ideal-gas Cp polynomial, J/mol/K, T in K
    estimated: dict[str, str] = field(default_factory=dict)  # property -> method
    notes: list[str] = field(default_factory=list)  # doubts about the values, for warnings

    def summary(self) -> str:
        by_method: dict[str, list[str]] = {}
        for prop, method in self.estimated.items():
            by_method.setdefault(method, []).append(prop)
        return "; ".join(f"{', '.join(props)} by {method}" for method, props in by_method.items())


def _resolve_groups(pairs, table: dict, kind: str, where: str) -> dict[int, int]:
    """Group names (case-insensitive) or numeric ids to ``{id: count}``."""
    by_name: dict[str, list[int]] = {}
    for gid, g in table.items():
        by_name.setdefault(g.group.lower(), []).append(gid)
    out: dict[int, int] = {}
    for key, count in pairs:
        gid = None
        if key.strip().isdigit() and int(key) in table:
            gid = int(key)
        elif key.strip().lower() in by_name:
            ids = by_name[key.strip().lower()]
            if len(ids) > 1:
                raise SimulationError(f"{where}: {kind} group name {key!r} is ambiguous",
                                      hint=f"use the group number instead: {ids}")
            gid = ids[0]
        if gid is None:
            raise SimulationError(f"{where}: unknown {kind} group {key!r}",
                                  hint=f"see the {kind} group table in the README "
                                  "(names or numbers)")
        out[gid] = out.get(gid, 0) + count
    return out


def _joback_from_smiles(smiles: str, where: str) -> dict[int, int] | None:
    """Joback groups from a SMILES string; None without RDKit."""
    try:
        from rdkit import Chem
    except ImportError:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SimulationError(f"{where}: SMILES {smiles!r} could not be parsed")
    fragments = Joback(mol)
    if not fragments.success:
        raise SimulationError(f"{where}: the SMILES could not be split into Joback groups "
                              f"({fragments.status})", hint="give joback_groups explicitly")
    return dict(fragments.counts)


def rdkit_available() -> bool:
    try:
        import rdkit  # noqa: F401
    except ImportError:
        return False
    return True


@lru_cache(maxsize=256)
def estimate(compound: CustomCompound) -> Estimate:
    """Fill in every constant the property methods need, or say what is missing."""
    where = f"custom compound {compound.name!r}"
    try:
        atoms = simple_formula_parser(compound.formula)
    except Exception as exc:  # chemicals raises plain exceptions for bad formulas
        raise SimulationError(f"{where}: formula {compound.formula!r} could not be read",
                              hint="write it like 'C2H6O' (element symbols and counts)") from exc
    if not atoms:
        raise SimulationError(f"{where}: formula {compound.formula!r} has no atoms")
    try:
        mw_formula = molecular_weight(atoms)
    except Exception as exc:
        raise SimulationError(f"{where}: formula {compound.formula!r} has an unknown element",
                              hint="write it like 'C2H6O' (element symbols and counts)") from exc
    given = dict(compound.constants)
    v: dict[str, float | None] = {k: given.get(k) for k in CONSTANTS}
    estimated: dict[str, str] = {}
    if v["MW"] is None:
        v["MW"] = mw_formula

    joback = _resolve_groups(compound.joback_groups, JOBACK_GROUPS, "Joback", where)
    if not joback and compound.smiles:
        joback = _joback_from_smiles(compound.smiles, where) or {}
        if joback:
            estimated["Joback groups"] = "RDKit fragmentation of the SMILES"
    unifac = _resolve_groups(compound.unifac_groups, DOUFSG, "UNIFAC (Dortmund)", where)
    if joback:
        _check_group_atoms(joback, JOBACK_GROUPS, "Joback", atoms, where)
    if unifac:
        _check_group_atoms(unifac, DOUFSG, "UNIFAC", atoms, where)

    n_atoms = sum(atoms.values())
    if v["Tb"] is None and compound.antoine is not None:
        A, B, C = compound.antoine[:3]
        v["Tb"] = B / (A - math.log10(101325.0)) - C
        estimated["Tb"] = "the given Antoine equation at 1 atm"
    if joback:
        jb = "Joback group contribution"
        if v["Tb"] is None:
            v["Tb"], estimated["Tb"] = Joback.Tb(joback), jb
        if v["Tc"] is None:
            v["Tc"], estimated["Tc"] = Joback.Tc(joback, v["Tb"]), jb
        if v["Pc"] is None:
            v["Pc"], estimated["Pc"] = Joback.Pc(joback, n_atoms), jb
        if v["Vc"] is None:
            v["Vc"], estimated["Vc"] = Joback.Vc(joback), jb
        if v["Tm"] is None:
            v["Tm"], estimated["Tm"] = Joback.Tm(joback), jb
        if v["Hf"] is None:
            v["Hf"], estimated["Hf"] = Joback.Hf(joback), jb
        if v["Gf"] is None:
            v["Gf"], estimated["Gf"] = Joback.Gf(joback), jb
    missing = [k for k in ("Tb", "Tc", "Pc") if v[k] is None]
    if missing:
        raise SimulationError(
            f"{where}: no {', '.join(CONSTANTS[k][0] for k in missing)}",
            hint="give the missing constants, or Joback groups (or a SMILES with RDKit "
            "installed) to estimate them",
        )
    if not v["Tb"] < v["Tc"]:
        raise SimulationError(f"{where}: the boiling point ({v['Tb']:.1f} K) must lie below "
                              f"the critical temperature ({v['Tc']:.1f} K)")
    if v["Pc"] <= 0.0:
        raise SimulationError(f"{where}: the critical pressure must be positive")
    if v["omega"] is None:
        v["omega"], estimated["omega"] = (LK_omega(v["Tb"], v["Tc"], v["Pc"]),
                                          "Lee-Kesler, from Tb, Tc and Pc")
    notes = []
    if estimated.get("Tb") == "Joback group contribution":
        notes.append("its boiling point is a Joback estimate, typically 10 to 40 K off, so "
                     "vapor pressures may be off by a factor of 2 or more; give a measured Tb "
                     "or Antoine constants if known")
    if not -0.1 <= v["omega"] <= 1.0:
        notes.append(f"its acentric factor {v['omega']:.2f} lies outside the usual range "
                     "(-0.1 to 1), which points to inconsistent Tb, Tc and Pc")
    if not 1.2 <= v["Tc"] / v["Tb"] <= 1.9:
        notes.append(f"Tc/Tb = {v['Tc'] / v['Tb']:.2f} is unusual (most compounds: 1.3 to "
                     "1.8); check Tb and Tc")
    if v["Vc"] is None:  # no groups: from Pitzer's Zc = 0.291 - 0.080 omega
        v["Vc"] = (0.291 - 0.080 * v["omega"]) * R * v["Tc"] / v["Pc"]
        estimated["Vc"] = "Pitzer correlation Zc = 0.291 - 0.080 omega"
    if v["r"] is None and v["q"] is None and unifac:
        rq = _uniquac_rq(unifac)
        if rq is not None:
            v["r"], v["q"] = rq
            estimated["UNIQUAC r, q"] = "UNIFAC group volumes and areas"
    cp = list(Joback.Cpig_coeffs(joback)) if joback else None
    estimated["ideal-gas heat capacity"] = (
        "Joback group contribution" if cp else "Lastovka-Shaw, from the formula")
    if compound.antoine is None:
        estimated["vapor pressure"] = {"lee-kesler": "Lee-Kesler",
                                       "ambrose-walton": "Ambrose-Walton"}[compound.vapor_pressure]
    estimated["liquid volume"] = "COSTALD"
    return Estimate(compound, v, atoms, joback, unifac, cp, estimated, notes)


def _check_group_atoms(groups: dict[int, int], table: dict, kind: str, atoms: dict[str, int],
                       where: str) -> None:
    """Heavy atoms of the groups must match the formula (hydrogens may be implicit)."""
    counted: dict[str, int] = {}
    for gid, n in groups.items():
        for el, k in (table[gid].atoms or {}).items():
            counted[el] = counted.get(el, 0) + k * n
    heavy = {el: n for el, n in atoms.items() if el != "H"}
    heavy_groups = {el: n for el, n in counted.items() if el != "H"}
    if heavy != heavy_groups:
        raise SimulationError(
            f"{where}: the {kind} groups do not match the formula",
            hint=f"groups give {_formula(heavy_groups)} (without H), the formula "
            f"{_formula(heavy)}",
        )


def _formula(atoms: dict[str, int]) -> str:
    return "".join(f"{el}{n if n > 1 else ''}" for el, n in sorted(atoms.items())) or "nothing"


def _uniquac_rq(unifac: dict[int, int]) -> tuple[float, float] | None:
    """UNIQUAC r and q from the original UNIFAC volumes and areas of the same subgroups."""
    r = q = 0.0
    for gid, n in unifac.items():
        name = DOUFSG[gid].group
        same = gid in UFSG and (UFSG[gid].group == name or name.startswith("OH("))
        if not same:
            return None
        sub = UFSG[gid if not name.startswith("OH(") else 14]
        r += n * sub.R
        q += n * sub.Q
    return r, q


def packages(compounds: tuple[CustomCompound, ...]):
    """``ChemicalConstantsPackage`` and ``PropertyCorrelationsPackage`` for the compounds."""
    ests = [estimate(c) for c in compounds]
    T_ref = 298.15
    cols: dict[str, list] = {k: [] for k in (
        "names", "CASs", "formulas", "smiless", "atomss", "MWs", "Tbs", "Tms", "Tcs", "Pcs",
        "Vcs", "Zcs", "omegas", "Hfgs", "Gfgs", "Sfgs", "dipoles", "UNIFAC_Dortmund_groups",
        "UNIFAC_groups", "UNIFAC_Rs", "UNIFAC_Qs", "similarity_variables")}
    for e in ests:
        c, v = e.compound, e.values
        cols["names"].append(c.name)
        cols["CASs"].append(c.id)
        cols["formulas"].append(c.formula)
        cols["smiless"].append(c.smiles or None)
        cols["atomss"].append(e.atoms)
        cols["MWs"].append(v["MW"])
        cols["Tbs"].append(v["Tb"])
        cols["Tms"].append(v["Tm"])
        cols["Tcs"].append(v["Tc"])
        cols["Pcs"].append(v["Pc"])
        cols["Vcs"].append(v["Vc"])
        cols["Zcs"].append(v["Pc"] * v["Vc"] / (R * v["Tc"]))
        cols["omegas"].append(v["omega"])
        cols["Hfgs"].append(v["Hf"])
        cols["Gfgs"].append(v["Gf"])
        cols["Sfgs"].append(None if v["Hf"] is None or v["Gf"] is None
                            else (v["Hf"] - v["Gf"]) / T_ref)
        cols["dipoles"].append(v["dipole"])
        cols["UNIFAC_Dortmund_groups"].append(e.unifac or None)
        cols["UNIFAC_groups"].append(None)
        cols["UNIFAC_Rs"].append(v["r"])
        cols["UNIFAC_Qs"].append(v["q"])
        cols["similarity_variables"].append(similarity_variable(e.atoms, v["MW"]))
    consts = ChemicalConstantsPackage(**cols)
    vps, cps, vls = [], [], []
    for e in ests:
        v = e.values
        vp = VaporPressure(Tb=v["Tb"], Tc=v["Tc"], Pc=v["Pc"], omega=v["omega"], CASRN="")
        vp.method = VAPOR_PRESSURE_METHODS[e.compound.vapor_pressure]
        if e.compound.antoine is not None:
            A, B, C, t_min, t_max = e.compound.antoine
            vp.add_correlation("given Antoine", "Antoine", Tmin=t_min or 0.5 * v["Tb"],
                               Tmax=t_max or v["Tc"], A=A, B=B, C=C, base=10.0)
        vps.append(vp)
        if e.cp_coeffs:
            a, b, c2, d = e.cp_coeffs
            # Joback's polynomial holds from about 280 to 1000 K; extended linearly beyond.
            cp = HeatCapacityGas(CASRN="", MW=v["MW"], similarity_variable=None,
                                 poly_fit=(280.0, 1000.0, [d, c2, b, a]),
                                 extrapolation="linear")
        else:
            cp = HeatCapacityGas(CASRN="", MW=v["MW"],
                                 similarity_variable=similarity_variable(e.atoms, v["MW"]))
            cp.method = "LASTOVKA_SHAW"
        cps.append(cp)
        Zc = v["Pc"] * v["Vc"] / (R * v["Tc"])
        vl = VolumeLiquid(MW=v["MW"], Tb=v["Tb"], Tc=v["Tc"], Pc=v["Pc"], Vc=v["Vc"], Zc=Zc,
                          omega=v["omega"], CASRN="")
        vl.method = "HTCOSTALD"
        vls.append(vl)
    # thermo fills the other correlations (viscosity, conductivity, ...) by corresponding-
    # states methods from the constants; its database lookups need a CAS number or none.
    blank = ChemicalConstantsPackage(**{**cols, "CASs": [""] * len(ests)})
    props = PropertyCorrelationsPackage(blank, VaporPressures=vps, HeatCapacityGases=cps,
                                        VolumeLiquids=vls)
    props.constants = consts
    return consts, props
