"""Property packages: the thermodynamic model behind every flash calculation.

A property package bundles the compound list with a phase-equilibrium method
(equation of state or activity-coefficient model). It is built once per
flowsheet and cached, because loading the constants takes a moment.
"""

from __future__ import annotations

import math
from contextlib import contextmanager
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations

import numpy as np
from fluids.numerics import OscillationError, UnconvergedError
from scipy.optimize import brentq, least_squares
from thermo import (
    NRTL,
    PR,
    PRMIX,
    SRK,
    SRKMIX,
    UNIFAC,
    UNIQUAC,
    CEOSGas,
    CEOSLiquid,
    ChemicalConstantsPackage,
    FlashPureVLS,
    FlashVL,
    FlashVLN,
    IdealGas,
    Wilson,
)
from thermo.bulk import BulkSettings
from thermo.interaction_parameters import IPDB
from thermo.unifac import DOUFIP2016, DOUFSG

from . import custom_compounds, plugins
from .activity_liquid import ActivityLiquid, ClausiusClapeyron
from .association import (
    DIMERIZATION,
    AssociatedLatentHeat,
    AssociatingGas,
    FlashPureAssociating,
    MonomerVaporPressure,
    vapor_from_monomers,
)
from .compounds import Compound, UnknownCompoundError, lookup
from .custom_compounds import CustomCompound
from .errors import SimulationError
from .fitting import store as fitted_store
from .flowsheet import BinaryParameterModel
from .henry import (
    HENRY_TC,
    HenryDilutedModel,
    HenryLiquid,
    henry_constant,
    ideal_solvent_model,
)

WATER_ICE = "water-ice"
_ICE_MELTING_T = 273.15  # K, rounded near atmospheric pressure
_ICE_MIN_T, _ICE_MAX_T = 263.15, 293.15  # K


@dataclass(frozen=True)
class Method:
    key: str
    label: str
    description: str
    plugin: str | None = None  # the plug-in (distribution) it comes from


METHODS: dict[str, Method] = {
    m.key: m
    for m in (
        Method(
            "peng-robinson",
            "Peng-Robinson",
            "Cubic equation of state. The workhorse for hydrocarbons, light gases and "
            "high pressures. Uses ChemSep binary interaction parameters where available.",
        ),
        Method(
            "srk",
            "Soave-Redlich-Kwong",
            "Cubic equation of state, similar to Peng-Robinson; common in gas processing.",
        ),
        Method(
            "unifac",
            "UNIFAC (Dortmund)",
            "Predictive activity-coefficient model for non-ideal liquids such as "
            "water/alcohol mixtures, at low to moderate pressure. Ideal vapor.",
        ),
        Method(
            "nrtl",
            "NRTL",
            "Activity-coefficient model with binary parameters regressed to measured "
            "phase-equilibrium data (Fluxsheet's NIST ThermoML fits, else the ChemSep "
            "databank). The first choice for strongly non-ideal "
            "liquids, azeotropes and liquid-liquid splits (decanters), at low to moderate "
            "pressure. Pairs without data are estimated from UNIFAC and flagged. Ideal vapor.",
        ),
        Method(
            "uniquac",
            "UNIQUAC",
            "Activity-coefficient model with binary parameters regressed to measured data "
            "(NIST ThermoML fits, else ChemSep); molecular size and "
            "surface from UNIFAC groups. Handles azeotropes and liquid-liquid splits. Pairs "
            "without data are estimated from UNIFAC and flagged. Ideal vapor.",
        ),
        Method(
            "wilson",
            "Wilson",
            "Activity-coefficient model with ChemSep binary parameters. Accurate for miscible "
            "polar mixtures such as alcohols, but by construction it cannot predict a "
            "liquid-liquid split. Pairs without data are estimated from UNIFAC and flagged.",
        ),
        Method(
            "ideal",
            "Ideal (Raoult's law)",
            "Ideal liquid and ideal gas. Fast, and fine for similar molecules at low pressure.",
        ),
        Method(
            WATER_ICE,
            "Water / ice reference",
            "Pure water and ice near atmospheric pressure. Temperature or enthalpy specifications.",
        ),
    )
}


# Property methods and compounds from plug-ins (fluxsheet.plugins), by key and lower-case name.
PLUGIN_METHODS: dict[str, plugins.ActivityMethod | plugins.EquationOfState] = {}
PLUGIN_COMPOUNDS: dict[str, CustomCompound] = {}


def is_eos(method: str) -> bool:
    """Whether a method describes both phases by an equation of state."""
    return method in ("peng-robinson", "srk") or isinstance(
        PLUGIN_METHODS.get(method), plugins.EquationOfState)


def _register_methods(plugin, value) -> list[str]:
    found = plugins.objects(value)
    kinds = (plugins.ActivityMethod, plugins.EquationOfState)
    if not found or not all(isinstance(m, kinds) for m in found):
        raise ValueError("expected an ActivityMethod or EquationOfState, or a list of them")
    keys = [m.key for m in found]
    taken = [k for k in keys if k in METHODS or k == "auto"]
    if taken or len(set(keys)) < len(keys) or not all(isinstance(k, str) and k for k in keys):
        raise ValueError(f"property method already defined: {', '.join(taken) or keys}")
    for m in found:
        if not callable(getattr(m, "gibbs_excess", None) or getattr(m, "eos_mix", None)):
            raise TypeError(f"{m.key}: no gibbs_excess function or eos_mix class")
    for m in found:
        PLUGIN_METHODS[m.key] = m
        METHODS[m.key] = Method(m.key, m.label, m.description,
                                plugin.distribution or plugin.name)
    return keys


def _register_compounds(plugin, value) -> list[str]:
    found = [c if isinstance(c, CustomCompound) else CustomCompound.from_dict(c)
             for c in plugins.objects(value)]
    if not found:
        raise ValueError("provides no compounds")
    names = [c.name.lower() for c in found]
    taken = [c.name for c in found if c.name.lower() in PLUGIN_COMPOUNDS]
    if taken or len(set(names)) < len(names):
        raise ValueError(f"compound already defined by a plug-in: {', '.join(taken) or names}")
    for c in found:
        try:
            lookup(c.name)
        except UnknownCompoundError:
            continue
        raise ValueError(f"{c.name!r} is already in the compound database")
    for c in found:
        PLUGIN_COMPOUNDS[c.name.lower()] = c
    return [c.name for c in found]


@contextmanager
def _plugin_errors(method: str):
    """A plug-in method's own failure, reported like any other flowsheet error."""
    try:
        yield
    except SimulationError:
        raise
    except Exception as exc:
        raise SimulationError(f"{METHODS[method].label}: {exc}",
                              hint=f"raised by the plug-in {METHODS[method].plugin!r}") from exc


def plugin_compound(name: str) -> Compound | None:
    """A compound defined by a plug-in, by name."""
    c = PLUGIN_COMPOUNDS.get(name.strip().lower())
    return None if c is None else _custom_compound(c)


def _custom_compound(c: CustomCompound) -> Compound:
    est = custom_compounds.estimate(c)
    return Compound(cas=custom_compounds.register(c), name=c.name, formula=c.formula,
                    mw=est.values["MW"], smiles=c.smiles, simulation_ready=True)


# UNIFAC subgroups (original and Dortmund numbering agree) of carboxylic acids: COOH, HCOOH.
_CARBOXYL = {42, 43}


@dataclass(frozen=True)
class PhaseState:
    """Result of a flash: everything per mole of mixture, SI units."""

    T: float
    P: float
    VF: float
    H: float  # J/mol, relative to ideal gas at 298.15 K (no formation enthalpy)
    S: float  # J/mol/K
    phase: str  # "vapor", "liquid", "two-phase", "solid" or "liquid-solid"
    x: tuple[float, ...] | None  # liquid mole fractions
    y: tuple[float, ...] | None  # vapor mole fractions
    V: float  # molar volume, m3/mol
    # H plus formation enthalpy: the basis for energy balances with reactions.
    # None when a component has no formation enthalpy in the database.
    H_reactive: float | None = None
    # Liquid-liquid split (vapor-liquid-liquid flash). ``x`` stays the combined liquid; these
    # describe its two phases, "light" and "heavy" by mass density. None without a split.
    x_light: tuple[float, ...] | None = None
    x_heavy: tuple[float, ...] | None = None
    heavy_fraction: float = 0.0  # mol heavy liquid per mol mixture; light = 1 - VF - this
    rho_light: float | None = None  # kg/m3
    rho_heavy: float | None = None  # kg/m3
    solid_fraction: float = 0.0  # mol solid per mol mixture
    x_solid: tuple[float, ...] | None = None

    @property
    def liquid_fraction(self) -> float:
        return max(1.0 - self.VF - self.solid_fraction, 0.0)

    @property
    def liquid_split(self) -> bool:
        return self.x_heavy is not None

    @property
    def light_fraction(self) -> float:
        return max(self.liquid_fraction - self.heavy_fraction, 0.0) if self.liquid_split else 0.0


@dataclass(frozen=True)
class PhaseProps:
    """One phase at fixed T, P and composition (no phase split). For stage-wise models."""

    lnphi: tuple[float, ...]  # log fugacity coefficients
    H: float  # J/mol
    H_reactive: float | None  # J/mol, see PhaseState
    V: float  # m3/mol


@dataclass(frozen=True)
class TransportProps:
    """One phase's properties for heat transfer correlations (film coefficients), SI."""

    rho: float  # density, kg/m3
    cp: float  # heat capacity, J/(kg K)
    mu: float  # dynamic viscosity, Pa s
    k: float  # thermal conductivity, W/(m K)

    @property
    def Pr(self) -> float:
        return self.cp * self.mu / self.k


@dataclass(frozen=True)
class Saturation:
    """A bubble or dew point: liquid composition x and vapor composition y in equilibrium."""

    T: float
    P: float
    x: tuple[float, ...]
    y: tuple[float, ...]


@dataclass(frozen=True)
class Suggestion:
    method: str
    reason: str


def suggest_method(components: list[str], custom: tuple[CustomCompound, ...] = (), *,
                   binary_parameters: tuple[BinaryParameterModel, ...] = (),
                   group_assignments: tuple[GroupAssignment, ...] = ()) -> Suggestion:
    """Pick a sensible default method and say why, the way an experienced engineer would."""
    consts = _constants(tuple(c.cas for c in resolve(components, custom)))[0]
    if len(components) == 1:
        return Suggestion(
            "ideal",
            f"A single compound has no mixing effects, so the ideal method is exact for phase "
            f"equilibrium: it uses {consts.names[0]}'s own vapor-pressure correlation. Switch to "
            "Peng-Robinson for high pressures or near the critical point.",
        )
    consts = _with_groups(consts, components, group_assignments)
    polar = [
        name
        for name, dipole in zip(consts.names, consts.dipoles)
        if dipole is not None and dipole >= 1.0
    ]
    splits = _liquid_liquid_pairs(consts, polar, bool(group_assignments))
    for method, table in _GE_TABLES.items():
        if splits and method not in LLE_METHODS:
            continue
        local = {frozenset(p.components) for p in binary_parameters if p.method == method}
        if local and all(
            frozenset((components[i], components[j])) in local
            or has_ip(table, consts.CASs[i], consts.CASs[j], "bij")
            or fitted_store.lookup(method, consts.CASs[i], consts.CASs[j]) is not None
            for i, j in combinations(range(consts.N), 2)
        ) and (method != "uniquac" or None not in (*consts.UNIFAC_Rs, *consts.UNIFAC_Qs)):
            return Suggestion(method, f"{METHODS[method].label} uses the flowsheet's local "
                              "binary parameters and has data for every pair.")
    if splits:
        pairs = " and ".join(f"{a}/{b}" for a, b, _, _ in splits)
        verb = "forms" if len(splits) == 1 else "form"
        lle = (f"UNIFAC predicts that {pairs} {verb} two liquid phases at 25 °C, so expect a "
               "liquid-liquid split (a Decanter can separate it). ")
        unmeasured = [f"{a}/{b}" for a, b, ca, cb in splits
                      if not has_ip(_GE_TABLES["nrtl"], ca, cb, "bij")
                      and fitted_store.lookup("nrtl", ca, cb) is None]
        if not unmeasured:
            return Suggestion(
                "nrtl",
                lle + "NRTL with binary parameters regressed to measured data (Fluxsheet's "
                "NIST ThermoML fits or the ChemSep databank) describes both liquid phases and "
                "the vapor in equilibrium with them; Wilson cannot represent a split at all.",
            )
        return Suggestion(
            "unifac",
            lle + f"There are no measured NRTL parameters (ChemSep or ThermoML fit) for "
            f"{', '.join(unmeasured)}, so "
            "the predictive UNIFAC model, which also captures liquid-liquid splits, is the "
            "safer choice.",
        )
    if polar and len(components) > 1:
        return Suggestion(
            "unifac",
            f"{', '.join(polar)} {'is' if len(polar) == 1 else 'are'} polar, so the liquid is "
            "likely non-ideal. An activity-coefficient model (UNIFAC) captures that; "
            "an equation of state like Peng-Robinson would not.",
        )
    return Suggestion(
        "peng-robinson",
        "All compounds are non-polar or only weakly polar, so a cubic equation of state "
        "(Peng-Robinson) describes both phases well, up to high pressure.",
    )


def resolve(components: list[str], custom: tuple[CustomCompound, ...] = ()) -> list[Compound]:
    """The compounds behind component names; the flowsheet's own definitions win over a
    plug-in's, and both over the database."""
    defined = {**PLUGIN_COMPOUNDS, **{c.name.lower(): c for c in custom}}
    out = []
    for name in components:
        c = defined.get(name.strip().lower())
        out.append(lookup(name) if c is None else _custom_compound(c))
    return out


def has_ip(table: str, cas_a: str, cas_b: str, key: str) -> bool:
    """Whether a binary parameter table has the pair; never for a user-defined compound."""
    if custom_compounds.is_custom(cas_a) or custom_compounds.is_custom(cas_b):
        return False
    return IPDB.has_ip_specific(table, [cas_a, cas_b], key)


def ip_matrix(table: str, cas_numbers, key: str) -> list[list[float]]:
    """A binary parameter matrix; pairs with a user-defined compound get the table's
    default (no data)."""
    idx = [i for i, cas in enumerate(cas_numbers) if not custom_compounds.is_custom(cas)]
    if len(idx) == len(cas_numbers):
        return IPDB.get_ip_asymmetric_matrix(table, list(cas_numbers), key)
    missing = IPDB.metadata[table]["missing"][key]
    n = len(cas_numbers)
    out = [[0.0 if i == j else missing for j in range(n)] for i in range(n)]
    known = IPDB.get_ip_asymmetric_matrix(table, [cas_numbers[i] for i in idx], key) if idx else []
    for a, i in enumerate(idx):
        for b, j in enumerate(idx):
            out[i][j] = known[a][b]
    return out


@lru_cache(maxsize=32)
def _constants(cas_numbers: tuple[str, ...]):
    """Constants and correlations; ids of custom compounds (see custom_compounds) may stand
    in for CAS numbers."""
    if not any(custom_compounds.is_custom(c) for c in cas_numbers):
        return ChemicalConstantsPackage.from_IDs(list(cas_numbers))
    parts = []
    for cas in cas_numbers:  # one package per compound, joined in the given order
        if custom_compounds.is_custom(cas):
            parts.append(custom_compounds.packages((custom_compounds.definition(cas),)))
        else:
            parts.append(ChemicalConstantsPackage.from_IDs([cas]))
    consts, props = parts[0]
    for c, p in parts[1:]:
        consts, props = consts + c, props + p
    return consts, props


# ---- activity-coefficient models ---------------------------------------------------------

# Binary parameter tables from the ChemSep databank, as shipped with `thermo`.
_GE_TABLES = {"nrtl": "ChemSep NRTL", "uniquac": "ChemSep UNIQUAC", "wilson": "ChemSep Wilson"}
# Methods with parameters regressed by Fluxsheet to NIST ThermoML data (fluxsheet.fitting).
_FITTED_METHODS = ("nrtl", "uniquac")
# Methods whose liquid model can describe two coexisting liquid phases.
LLE_METHODS = ("nrtl", "uniquac", "unifac")
# Activity-coefficient methods, where light gases are Henry components (see henry.py).
HENRY_METHODS = ("nrtl", "uniquac", "unifac", "wilson")
# Beyond these the Henry treatment (ideal dilute solution, no Poynting term) is flagged.
HENRY_X_MAX = 0.03
HENRY_P_MAX = 20e5


@dataclass(frozen=True)
class GroupAssignment:
    component: str
    groups: tuple[tuple[int, int], ...]
    source: str


def _with_groups(consts, names: list[str], assignments: tuple[GroupAssignment, ...]):
    if not assignments:
        return consts
    groups = list(consts.UNIFAC_Dortmund_groups)
    for assignment in assignments:
        for group, _ in assignment.groups:
            if group not in DOUFSG:
                raise SimulationError(f"unknown Dortmund UNIFAC subgroup {group}",
                                      where=assignment.component,
                                      hint="use a subgroup ID from the Dortmund UNIFAC tables")
        groups[names.index(assignment.component)] = dict(assignment.groups)
    return consts.with_new_constants(UNIFAC_Dortmund_groups=groups)


def _unifac_gaps(groups) -> list[tuple[int, int]]:
    main = sorted({DOUFSG[i].main_group_id for counts in groups for i in (counts or {})})
    return [(i, j) for i, j in permutations(main, 2) if j not in DOUFIP2016.get(i, {})]


def _unifac(consts, T: float = 298.15, idx: list[int] | None = None):
    """UNIFAC (Dortmund) model for these compounds (or those at ``idx``), or None if a
    compound lacks groups."""
    groups = [consts.UNIFAC_Dortmund_groups[i] for i in (idx or range(consts.N))]
    if any(not g for g in groups):
        return None
    return UNIFAC.from_subgroups(
        T=T, xs=[1.0 / len(groups)] * len(groups), chemgroups=groups,
        version=1, interaction_data=DOUFIP2016, subgroups=DOUFSG,
    )


def _ge_model(method: str, consts, a, b, alpha, idx: list[int] | None = None):
    """A regressed activity model; ``a``, ``b`` (and ``alpha`` for NRTL) are n x n matrices.
    With ``idx``, the model of only those compounds (rows and columns of the matrices)."""
    idx = list(idx if idx is not None else range(consts.N))

    def sub(m):
        return None if m is None else [[m[i][j] for j in idx] for i in idx]

    a, b, alpha = sub(a), sub(b), sub(alpha)
    xs = [1.0 / len(idx)] * len(idx)
    if method == "nrtl":  # tau_ij = a_ij + b_ij / T
        return NRTL(T=298.15, xs=xs, tau_as=a, tau_bs=b, alpha_cs=alpha)
    if method == "uniquac":  # tau_ij = exp(a_ij + b_ij / T)
        return UNIQUAC(T=298.15, xs=xs, rs=[consts.UNIFAC_Rs[i] for i in idx],
                       qs=[consts.UNIFAC_Qs[i] for i in idx], tau_as=a, tau_bs=b)
    return Wilson(T=298.15, xs=xs, lambda_as=a, lambda_bs=b)  # lambda_ij = exp(a + b / T)


_X_GRID = [i / 60 for i in range(1, 60)]


@lru_cache(maxsize=256)
def _binary_splits(cas_a: str, cas_b: str, T: float = 298.15) -> bool:
    """Whether UNIFAC predicts two liquid phases for this binary at ``T``.

    The Gibbs energy of mixing of a stable liquid is convex in composition; a region where it
    curves downward means the liquid lowers its energy by splitting.
    """
    return _unifac_splits(_unifac(_constants((cas_a, cas_b))[0], T), T)


def _unifac_splits(model, T: float) -> bool:
    if model is None:
        return False
    g = []
    for x in _X_GRID:
        gammas = model.to_T_xs(T, [x, 1.0 - x]).gammas()
        g.append(x * math.log(x * gammas[0]) + (1 - x) * math.log((1 - x) * gammas[1]))
    return any(g[i - 1] - 2 * g[i] + g[i + 1] < -1e-9 for i in range(1, len(g) - 1))


def _liquid_liquid_pairs(consts, polar: list[str], local_groups: bool = False
                         ) -> list[tuple[str, str, str, str]]:
    """Binary pairs (names, CAS numbers) likely to split; only pairs with a polar member."""
    out = []
    for i in range(consts.N):
        for j in range(i + 1, consts.N):
            a, b = consts.names[i], consts.names[j]
            if a not in polar and b not in polar:
                continue
            split = (_unifac_splits(_unifac(consts, idx=[i, j]), 298.15) if local_groups
                     else _binary_splits(consts.CASs[i], consts.CASs[j]))
            if split:
                out.append((a, b, consts.CASs[i], consts.CASs[j]))
    return out


@lru_cache(maxsize=128)
def _estimate_binary(method: str, cas_a: str, cas_b: str):
    """Fit binary parameters to UNIFAC activity coefficients when no measured set exists.

    Returns ``(a12, a21, b12, b21, alpha)`` or None when UNIFAC cannot describe the pair.
    The fit covers the whole composition range at three temperatures spanning the two
    boiling points, so it reproduces UNIFAC's predictions, not measured data.
    """
    return _fit_unifac_binary(method, _constants((cas_a, cas_b))[0])


def _fit_unifac_binary(method: str, consts):
    unifac = _unifac(consts)
    if unifac is None or (method == "uniquac" and None in (*consts.UNIFAC_Rs, *consts.UNIFAC_Qs)):
        return None
    tbs = [tb for tb in consts.Tbs if tb is not None] or [350.0]
    t_lo, t_hi = max(min(tbs) - 40.0, 250.0), min(max(max(tbs), min(tbs) + 20.0), 500.0)
    temps = [t_lo, 0.5 * (t_lo + t_hi), t_hi]
    xs = [0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 0.98]
    target = [math.log(g) for T in temps for x in xs
              for g in unifac.to_T_xs(T, [x, 1.0 - x]).gammas()]
    # Renon and Prausnitz: alpha = 0.2 suits partially miscible pairs, 0.3 the others.
    alpha_value = 0.2 if _unifac_splits(unifac, 298.15) else 0.3
    alpha = [[0.0, alpha_value], [alpha_value, 0.0]]

    def residuals(p):
        a12, a21, b12, b21 = p
        model = _ge_model(method, consts, [[0.0, a12], [a21, 0.0]], [[0.0, b12], [b21, 0.0]],
                          alpha)
        calc = [math.log(g) for T in temps for x in xs
                for g in model.to_T_xs(T, [x, 1.0 - x]).gammas()]
        # A light pull towards a = 0 keeps the correlated a and b from drifting apart.
        return [c - t for c, t in zip(calc, target)] + [0.05 * a12, 0.05 * a21]

    try:
        fit = least_squares(residuals, [0.0, 0.0, 0.0, 0.0], x_scale=[1.0, 1.0, 300.0, 300.0])
    except (ValueError, OverflowError, ZeroDivisionError):
        return None
    if not fit.success:
        return None
    return (*(float(v) for v in fit.x), alpha_value)


def _water_ice_specification_error() -> SimulationError:
    return SimulationError(
        "the water-ice reference method supports only T-P and P-H flashes",
        hint="specify temperature or enthalpy with pressure, or select a fluid property method",
    )


class PropertyPackage:
    def __init__(self, components: list[str], method: str, prefer_fitted: bool = True,
                 henry_components: tuple[str, ...] | None = None,
                 custom: tuple[CustomCompound, ...] = (), *,
                 binary_parameters: tuple[BinaryParameterModel, ...] = (),
                 group_assignments: tuple[GroupAssignment, ...] = ()):
        """``prefer_fitted``: for NRTL and UNIQUAC, use Fluxsheet's ThermoML fits over ChemSep
        parameters when both exist (ThermoML fits always fill pairs ChemSep lacks).
        ``henry_components``: names of the Henry components of activity methods; None picks
        them automatically (critical temperature below ``HENRY_TC``), () means none.
        ``custom``: definitions of user-defined compounds among ``components``."""
        if method not in METHODS:
            raise SimulationError(
                f"unknown property method {method!r}",
                hint=f"choose one of: {', '.join(METHODS)}",
            )
        if not components:
            raise SimulationError(
                "the flowsheet has no components", hint="add at least one compound"
            )
        self.method = method
        self.prefer_fitted = prefer_fitted
        self.custom = tuple(custom)
        self.binary_parameters = binary_parameters
        self.group_assignments = group_assignments
        self.compounds = resolve(list(components), self.custom)
        self.names = list(components)
        cas = tuple(c.cas for c in self.compounds)
        if len(set(cas)) != len(cas):
            raise SimulationError(
                "the component list contains the same compound twice",
                hint="different names can refer to the same compound (e.g. 'ethanol' "
                "and 'ethyl alcohol'); remove one",
            )
        if self.method == WATER_ICE and cas != ("7732-18-5",):
            raise SimulationError(
                "the water-ice property method requires pure water",
                hint="set water as the only component or choose another property method",
            )
        self._database_consts, self.props = _constants(cas)
        # Where each binary parameter set came from, and warnings worth showing with results.
        self.parameter_sources: list[dict] = []
        self.warnings: list[str] = []
        self.assumptions: list[str] = []
        self.consts = _with_groups(self._database_consts, self.names, group_assignments)
        for assignment in group_assignments:
            self.parameter_sources.append({
                "pair": assignment.component, "kind": "unifac", "source": "Flowsheet groups",
                "reference": assignment.source, "convention": "dortmund",
                "groups": dict(assignment.groups),
            })
        self._ge = None  # activity-coefficient model, if the method has one
        self._shifts: list[float] | None = None  # volume translations (equations of state)
        self._assoc: dict[int, tuple[float, float]] = {}  # dimerizing acids (activity methods)
        self._liquid_vps = list(self.props.VaporPressures)  # liquid reference fugacities
        self._report_estimates()
        self._check_data()
        # Henry components (dissolved gases) of activity methods: flags and names.
        self._henry_request = henry_components
        self.henry = [False] * self.n
        self._flasher = self._build()
        self._warn_association()
        if henry_components and not (self.has_henry_model and self.n > 1):
            self.warnings.append(
                f"henry_components is ignored by {METHODS[method].label}: "
                + ("equations of state describe dissolved gases directly." if is_eos(method)
                   else "a single compound has no solvent." if self.has_henry_model
                   else "this plug-in method has no Henry components (ActivityMethod "
                   "henry=False)."))

    def _carboxylic(self) -> list[bool]:
        """Which components are carboxylic acids (by their UNIFAC groups)."""
        return [bool(_CARBOXYL & set(groups or dortmund or ()))
                for groups, dortmund in zip(self._database_consts.UNIFAC_groups,
                                            self._database_consts.UNIFAC_Dortmund_groups)]

    def _warn_association(self) -> None:
        """Carboxylic acids dimerize in the vapor.

        The activity methods describe it for the acids in ``association.DIMERIZATION``
        (chemical theory). An acid without dimerization constants keeps an ideal monomer
        vapor and the liquid enthalpy of Clausius-Clapeyron: half of the enthalpy between
        the liquid and the monomer vapor (51 kJ/mol for acetic acid at 25 °C) is the dimers'
        bond, so with an ideal monomer vapor a liquid enthalpy can match either the measured
        heat of vaporization (23 kJ/mol) or the formation enthalpies, not both, and
        Clausius-Clapeyron (about 40 kJ/mol) lies between them.
        """
        acids = [name for name, acid in zip(self.names, self._carboxylic()) if acid]
        if not is_eos(self.method):
            described = set(self._associating())
            acids = [name for name in acids if self.names.index(name) not in described]
        if not acids:
            return
        label = METHODS[self.method].label
        names = ", ".join(acids)
        if is_eos(self.method):
            self.warnings.append(
                f"{label}: {names} associate{'s' if len(acids) == 1 else ''} in the vapor "
                "(carboxylic acid dimers), which an equation of state without association "
                "cannot describe: the heat of vaporization comes out 20 to 60 % too high, and "
                "with it condenser and reboiler duties. Use NRTL, UNIQUAC or UNIFAC, which "
                "describe the dimers of acetic, formic and propionic acid.")
        else:
            self.warnings.append(
                f"{label}: {names} associate{'s' if len(acids) == 1 else ''} in the vapor "
                "(carboxylic acid dimers), but no dimerization constants are known for "
                + ("it" if len(acids) == 1 else "them")
                + ", so the vapor is treated as an ideal gas of monomers: the heat of "
                "vaporization comes out 40 to 65 % too high, liquid-phase heats of reaction "
                "about 10 kJ per mole of acid off, and vapor densities too low.")

    def _associating(self) -> dict[int, tuple[float, float]]:
        """Acids whose vapor dimerization is described: index -> constants."""
        return {i: DIMERIZATION[cas] for i, (cas, acid) in
                enumerate(zip(self.consts.CASs, self._carboxylic()))
                if acid and cas in DIMERIZATION}

    @property
    def henry_components(self) -> list[str]:
        return self.henry_names()

    def henry_validity(self, states) -> list[str]:
        """Warnings for solved streams where Henry's law as used here is stretched: a
        concentrated solute (ideal dilute solution) or high pressure (no Poynting
        correction). ``states``: stream id -> PhaseState."""
        if not any(self.henry):
            return []
        concentrated: dict[str, list[str]] = {}
        high_p: list[str] = []
        for sid, st in states.items():
            if st.x is None or st.VF >= 1.0:
                continue
            for i in (i for i in range(self.n) if self.henry[i]):
                if st.x[i] > HENRY_X_MAX:
                    concentrated.setdefault(self.names[i], []).append(sid)
            if st.P > HENRY_P_MAX:
                high_p.append(sid)
        out = [f"{name} is a Henry component but makes up more than {HENRY_X_MAX:.0%} of the "
               f"liquid in {', '.join(sids)}; Henry's law assumes a dilute solution there. "
               "Use Peng-Robinson, or leave it out of henry_components if it condenses."
               for name, sids in concentrated.items()]
        if high_p:
            out.append(f"liquid with dissolved gases above {HENRY_P_MAX / 1e5:.0f} bar in "
                       f"{', '.join(high_p)}: the Henry constants have no pressure correction "
                       "(Krichevsky-Kasarnovsky), so solubilities may be some 10 % off; "
                       "consider Peng-Robinson.")
        return out

    def _report_estimates(self) -> None:
        """Say which constants of user-defined compounds are estimates."""
        for comp in self.compounds:
            if not custom_compounds.is_custom(comp.cas):
                continue
            est = custom_compounds.estimate(custom_compounds.definition(comp.cas))
            for prop, how in est.estimated.items():
                value = est.values.get(prop)
                self.parameter_sources.append({
                    "pair": comp.name, "kind": "estimate", "source": f"{prop}: {how}",
                    **({"value": value} if isinstance(value, float) else {})})
            self._fallback(
                f"{comp.name} is a user-defined compound; estimated: {est.summary()}. Treat "
                "results that depend on it as estimates."
                + "".join(f" Note: {note}." for note in est.notes))

    def _fallback(self, message: str) -> None:
        self.assumptions.append(message)
        self.warnings.append(message)

    def require_property_data(self, where: str) -> None:
        if self.assumptions:
            raise SimulationError(
                "strict property data rejects these assumptions: " + " ".join(self.assumptions),
                where=where, hint="supply the missing data or select a method with suitable data",
            )

    @property
    def n(self) -> int:
        return len(self.names)

    @property
    def mws(self) -> list[float]:
        return list(self.consts.MWs)  # g/mol

    @property
    def has_henry_model(self) -> bool:
        """Whether light gases can be Henry components: the activity methods in
        ``HENRY_METHODS`` and plug-in activity methods with ``henry=True``."""
        plugin = PLUGIN_METHODS.get(self.method)
        if isinstance(plugin, plugins.ActivityMethod):
            return plugin.henry
        return self.method in HENRY_METHODS

    @property
    def supports_lle(self) -> bool:
        """Whether flashes can find two liquid phases (vapor-liquid-liquid equilibrium)."""
        plugin = PLUGIN_METHODS.get(self.method)
        lle = plugin.liquid_liquid if isinstance(plugin, plugins.ActivityMethod) else (
            self.method in LLE_METHODS)
        return lle and self.n > 1

    def _check_data(self) -> None:
        for field, label in (("Tcs", "critical temperature"), ("Pcs", "critical pressure"),
                             ("omegas", "acentric factor")):
            for name, value in zip(self.names, getattr(self.consts, field)):
                if value is None:
                    raise SimulationError(
                        f"no {label} available for {name!r}",
                        hint="this compound lacks the data needed for simulation; "
                        "pick a different compound or a pure-component alternative",
                    )

    def _build(self):
        c, p = self.consts, self.props
        cp_gas = p.HeatCapacityGases
        plugin = PLUGIN_METHODS.get(self.method)
        if isinstance(plugin, plugins.EquationOfState):
            n = self.n
            with _plugin_errors(self.method):
                kijs = [list(map(float, row)) for row in plugin.kijs(c)] if plugin.kijs else (
                    [[0.0] * n for _ in range(n)])
                extra = plugin.extra(c) if plugin.extra else {}
            if len(kijs) != n or any(len(row) != n for row in kijs):
                raise SimulationError(f"{plugin.label}: kijs must be a {n} x {n} matrix")
            kwargs = {"Tcs": c.Tcs, "Pcs": c.Pcs, "omegas": c.omegas, "kijs": kijs, **extra}
            gas = CEOSGas(plugin.eos_mix, kwargs, HeatCapacityGases=cp_gas)
            liquid = CEOSLiquid(plugin.eos_mix, kwargs, HeatCapacityGases=cp_gas)
        elif self.method in ("peng-robinson", "srk"):
            pr = self.method == "peng-robinson"
            eos = PRMIX if pr else SRKMIX
            kijs = (
                ip_matrix("ChemSep PR", c.CASs, "kij") if pr
                else [[0.0] * self.n for _ in range(self.n)]
            )
            missing = [f"{self.names[i]}/{self.names[j]}" for i, j in combinations(range(self.n), 2)
                       if not pr or not has_ip("ChemSep PR", c.CASs[i], c.CASs[j], "kij")]
            if missing:
                self.assumptions.append(
                    f"{METHODS[self.method].label}: binary interactions set to zero "
                    f"(kij = 0) for {', '.join(missing)}.")
            kwargs = {"Tcs": c.Tcs, "Pcs": c.Pcs, "omegas": c.omegas, "kijs": kijs}
            self._shifts = self._volume_translations(PR if pr else SRK)
            gas = CEOSGas(eos, kwargs, HeatCapacityGases=cp_gas)
            liquid = CEOSLiquid(eos, kwargs, HeatCapacityGases=cp_gas)
        else:
            self._assoc = self._associating()
            gas = (AssociatingGas(self._assoc, HeatCapacityGases=cp_gas) if self._assoc
                   else IdealGas(HeatCapacityGases=cp_gas))
            henry_kwargs: dict = {}
            if self.has_henry_model and self.n > 1:
                self.henry, henry_kwargs = self._henry_setup()
            solvents = [i for i in range(self.n) if not self.henry[i]]
            ge_model = None
            if isinstance(plugin, plugins.ActivityMethod):
                # Henry components left out: the solvents only (all without henry=True).
                with _plugin_errors(self.method):
                    ge_model = plugin.gibbs_excess(c, solvents)
                if ge_model is not None and getattr(ge_model, "N", len(solvents)) != len(
                        solvents):
                    raise SimulationError(
                        f"{plugin.label}: gibbs_excess returned a model of {ge_model.N} "
                        f"compounds for {len(solvents)}",
                        hint=f"from the plug-in {METHODS[self.method].plugin!r}; it must "
                        "describe exactly the components whose indices it is given")
            elif self.method == "unifac":
                missing = [self.names[i] for i in solvents if not c.UNIFAC_Dortmund_groups[i]]
                if missing:
                    raise SimulationError(
                        f"No Dortmund UNIFAC group assignment for {', '.join(missing)}",
                        hint="define the compound's groups in property_data.unifac_groups "
                             "or select a property method with the necessary data",
                    )
                gaps = _unifac_gaps(c.UNIFAC_Dortmund_groups[i] for i in solvents)
                if gaps:
                    pairs = ', '.join(f'{i}→{j}' for i, j in gaps)
                    self.assumptions.append(
                        "UNIFAC (Dortmund): missing main-group interactions use "
                        f"neutral factors (psi = 1): {pairs}.")
                ge_model = _unifac(c, idx=solvents)
            elif self.method in _GE_TABLES and len(solvents) > 1:
                ge_model = self._regressed_model(solvents)
            if henry_kwargs:
                ge_model = HenryDilutedModel(ge_model or ideal_solvent_model(len(solvents)),
                                             solvents, self.n)
            self._ge = ge_model
            # Fugacities from vapor pressures and activity coefficients; the liquid enthalpy
            # from the heat of vaporization correlations (see activity_liquid.py).
            # A dimerizing acid's liquid is in equilibrium with the monomer pressure of its
            # saturated vapor, and its enthalpy lies below the monomer gas by the heat of
            # vaporization plus the dimers' bond enthalpy (association.py). Acids without
            # dimerization constants keep Clausius-Clapeyron: see _warn_association.
            self._liquid_vps = [MonomerVaporPressure(vp, self._assoc[i]) if i in self._assoc
                                else vp for i, vp in enumerate(p.VaporPressures)]
            latent = []
            for i, (acid, correlation) in enumerate(zip(self._carboxylic(),
                                                        p.EnthalpyVaporizations)):
                if i in self._assoc:
                    latent.append(AssociatedLatentHeat(correlation, p.VaporPressures[i],
                                                       self._assoc[i]))
                else:
                    latent.append(ClausiusClapeyron() if acid else correlation)
            liquid = (HenryLiquid if henry_kwargs else ActivityLiquid)(
                VaporPressures=self._liquid_vps,
                HeatCapacityGases=cp_gas,
                VolumeLiquids=p.VolumeLiquids,
                EnthalpyVaporizations=latent,
                GibbsExcessModel=ge_model,
                equilibrium_basis="Psat",
                caloric_basis="Hvap",
                **henry_kwargs,
            )
        if self.n == 1:
            flasher = FlashPureAssociating if self._assoc else FlashPureVLS
            return flasher(c, p, gas=gas, liquids=[liquid], solids=[], settings=_SETTINGS)
        if self.supports_lle:
            return _FlashVLLE(c, p, liquids=[liquid, liquid], gas=gas, settings=_SETTINGS)
        return FlashVL(c, p, liquid=liquid, gas=gas, settings=_SETTINGS)

    def _volume_translations(self, pure_eos) -> list[float]:
        """Volume translation c_i (m3/mol) of each component for the cubic equation of state.

        Cubic equations give poor liquid volumes: Peng-Robinson is 15 % too large for water
        and 13 % for methanol, Soave-Redlich-Kwong 8 to 24 % for everything. Pump power,
        velocities, levels, residence times and liquid concentrations (kinetics) inherit
        that. Peneloux's constant shift V = V_eos - sum(x_i c_i) corrects it without changing
        any phase equilibrium: every fugacity gets the same factor exp(-c_i P / R T). So the
        equation is solved as it is and the shift is taken off the volumes it reports (see
        ``_volume``); enthalpies stay those of the plain equation.

        c_i makes the saturated liquid volume right at 25 °C, or at Tr = 0.7 for compounds
        whose critical temperature lies below 426 K (the classical Peneloux reference). A
        compound without a liquid volume correlation keeps the plain equation (c_i = 0).
        """
        c = self.consts
        shifts = []
        for i in range(self.n):
            shift = 0.0
            try:
                T_ref = min(298.15, 0.7 * c.Tcs[i])
                measured = self.props.VolumeLiquids[i].T_dependent_property(T_ref)
                if measured:
                    eos = pure_eos(Tc=c.Tcs[i], Pc=c.Pcs[i], omega=c.omegas[i], T=T_ref, P=1e5)
                    shift = float(eos.V_l_sat(T_ref)) - float(measured)
            except Exception:  # noqa: BLE001 - no data or no liquid root: plain equation
                shift = 0.0
            shifts.append(shift if math.isfinite(shift) else 0.0)
        return shifts

    def _volume(self, V: float, zs) -> float:
        """Molar volume with the equation of state's volume translation taken off. For a
        mixture of phases the overall composition gives the same as phase by phase."""
        if self._shifts is None:
            return float(V)
        return float(V) - sum(z * c for z, c in zip(zs, self._shifts))

    def _henry_setup(self) -> tuple[list[bool], dict]:
        """Which components follow Henry's law, and their constants for the liquid.

        A candidate (critical temperature below ``HENRY_TC``, or listed explicitly) becomes
        a Henry component when a constant exists for it in at least one solvent present.
        Solvents without a constant get a rough estimate, so that every solvent weighs in the
        mixing rule with its real mole fraction.
        """
        c, n, request = self.consts, self.n, self._henry_request
        label = METHODS[self.method].label
        if request is None:
            candidate = [not custom_compounds.is_custom(cas) and tc is not None and tc < HENRY_TC
                         for cas, tc in zip(c.CASs, c.Tcs)]
        else:
            wanted = {name.strip().lower() for name in request}
            unknown = wanted - {name.lower() for name in self.names}
            if unknown:
                raise SimulationError(
                    f"henry_components lists {', '.join(sorted(unknown))}, which "
                    f"{'is' if len(unknown) == 1 else 'are'} not among the components",
                    hint="use the component names exactly as in the component list")
            candidate = [name.lower() in wanted for name in self.names]
            if all(candidate):
                raise SimulationError(
                    "every component is listed as a Henry component",
                    hint="Henry's law describes gases dissolved in a solvent; leave at least "
                    "one component (the solvent) out of henry_components")
        if all(candidate):  # gases only: no solvent to dissolve in
            return [False] * n, {}
        solvents = [j for j in range(n) if not candidate[j]]
        coeffs = [[(0.0,) * 6 for _ in range(n)] for _ in range(n)]
        flags = [False] * n
        for i in (i for i in range(n) if candidate[i]):
            found = {j: k for j in solvents
                     if (k := henry_constant(c.CASs[i], c.CASs[j])) is not None}
            missing = [j for j in solvents if j not in found]
            if not found:
                self._fallback(
                    f"{label}: no Henry constant for {self.names[i]} in "
                    f"{', '.join(self.names[j] for j in missing)}; its liquid fugacity uses its "
                    "vapor pressure correlation, extrapolated above the critical point. Its "
                    "solubility is not reliable; use Peng-Robinson for this gas.")
                continue
            flags[i] = True
            for j, k in found.items():
                coeffs[i][j] = k.coeffs
                self.parameter_sources.append({"pair": f"{self.names[i]} in {self.names[j]}",
                                               "kind": "henry",
                                               "source": f"Henry constant: {k.source}"})
                if k.source.startswith("predicted") or "25 °C value only" in k.source:
                    self.assumptions.append(f"{label}: Henry constant for {self.names[i]} in "
                                            f"{self.names[j]} uses {k.source}.")
            if missing:
                # thermo's mixing rule ("solvents_with_parameters") then averages over the
                # solvents with data only; an extrapolated vapor pressure would be far worse.
                for j in missing:
                    self.parameter_sources.append({
                        "pair": f"{self.names[i]} in {self.names[j]}", "kind": "henry",
                        "source": "Henry constant: none (other solvents' constants used)"})
                self._fallback(
                    f"{label}: no Henry constant for {self.names[i]} in "
                    f"{', '.join(self.names[j] for j in missing)}; in liquids rich in "
                    f"{'that solvent' if len(missing) == 1 else 'those solvents'} its "
                    "solubility is taken from the constants of the other solvents "
                    f"({', '.join(self.names[j] for j in found)}) and may be off by a factor "
                    "of several.")
        if not any(flags):
            return flags, {}
        names = ", ".join(self.henry_names(flags))
        self.warnings.append(
            f"{label}: {names} {'is a Henry component' if sum(flags) == 1 else 'are Henry components'}"
            " (dissolved gas: Henry's law, ideal dilute solution, no pressure correction). "
            "Valid for dilute solutions at moderate pressure; set henry_components to change.")
        kwargs = {"henry_components": flags, "henry_mode": "solvents_with_parameters"}
        for k, key in enumerate(("henry_as", "henry_bs", "henry_cs", "henry_ds", "henry_es",
                                 "henry_fs")):
            kwargs[key] = [[coeffs[i][j][k] for j in range(n)] for i in range(n)]
        return flags, kwargs

    def henry_names(self, flags=None) -> list[str]:
        return [name for name, h in zip(self.names, flags or self.henry) if h]

    def _regressed_model(self, idx: list[int]):
        """NRTL, UNIQUAC or Wilson with measured parameters; UNIFAC estimates fill the gaps.

        Per pair: flowsheet data, Fluxsheet's ThermoML fit (first if ``prefer_fitted``),
        the ChemSep set, a ThermoML fit, then a UNIFAC estimate. Only the compounds at ``idx``
        (the solvents: Henry components are left out).
        """
        c, n, table = self.consts, self.n, _GE_TABLES[self.method]
        if self.method == "uniquac":
            for i in idx:
                if c.UNIFAC_Rs[i] is None or c.UNIFAC_Qs[i] is None:
                    custom = custom_compounds.is_custom(c.CASs[i])
                    raise SimulationError(
                        f"{self.names[i]!r} has no UNIQUAC size and surface parameters (r, q)",
                        hint=("give both r and q in its constants, or unifac_groups made of "
                              "original UNIFAC subgroups" if custom else
                              "they come from UNIFAC groups, which this compound lacks; "
                              "use NRTL or Peng-Robinson instead"),
                    )
        b = ip_matrix(table, c.CASs, "bij")
        a = (ip_matrix(table, c.CASs, "aij") if self.method == "wilson"
             else [[0.0] * n for _ in range(n)])
        alpha = (ip_matrix(table, c.CASs, "alphaij")
                 if self.method == "nrtl" else None)
        label = METHODS[self.method].label
        local = {frozenset(p.components): p for p in self.binary_parameters
                 if p.method == self.method}
        local_groups = {g.component for g in self.group_assignments}
        for i in idx:
            for j in (j for j in idx if j > i):
                names = (self.names[i], self.names[j])
                pair = "/".join(names)
                if (params := local.get(frozenset(names))) is not None:
                    forward = names == params.components
                    a[i][j], a[j][i] = params.a if forward else params.a[::-1]
                    b[i][j], b[j][i] = params.b if forward else params.b[::-1]
                    if alpha is not None:
                        alpha[i][j] = alpha[j][i] = params.alpha
                    self.parameter_sources.append({
                        "pair": pair, "source": "Flowsheet", "reference": params.source,
                        "convention": params.convention, "equation": fitted_store.equation(self.method),
                        "T_range_K": params.T_range_K, "P_range_Pa": params.P_range_Pa,
                        "a": (a[i][j], a[j][i]), "b": (b[i][j], b[j][i]), "alpha": params.alpha,
                    })
                    continue
                chemsep = has_ip(table, c.CASs[i], c.CASs[j], "bij")
                fitted = (fitted_store.lookup(self.method, c.CASs[i], c.CASs[j])
                          if self.method in _FITTED_METHODS else None)
                if fitted is not None and (self.prefer_fitted or not chemsep):
                    a[i][j], a[j][i] = fitted.a_ij, fitted.a_ji
                    b[i][j], b[j][i] = fitted.b_ij, fitted.b_ji
                    if alpha is not None:
                        alpha[i][j] = alpha[j][i] = fitted.alpha
                    self.parameter_sources.append({"pair": pair, "source": "ThermoML fit",
                                                   **fitted.provenance()})
                    continue
                if chemsep:
                    self.parameter_sources.append({"pair": pair, "source": "ChemSep"})
                    continue
                est = (_fit_unifac_binary(self.method, c.subset([i, j]))
                       if local_groups.intersection(names)
                       else _estimate_binary(self.method, c.CASs[i], c.CASs[j]))
                if est is None:
                    self.parameter_sources.append({"pair": pair, "source": "none (ideal)"})
                    self._fallback(
                        f"{label}: no binary parameters for {pair}, and UNIFAC cannot estimate "
                        "them; this pair is treated as an ideal solution. Use a method with "
                        "data for it, or Peng-Robinson for non-polar mixtures."
                    )
                    continue
                a[i][j], a[j][i], b[i][j], b[j][i], alpha_value = est
                if alpha is not None:
                    alpha[i][j] = alpha[j][i] = alpha_value
                self.parameter_sources.append({"pair": pair, "source": "UNIFAC estimate"})
                self._fallback(
                    f"{label}: there are no measured binary parameters (ChemSep or ThermoML "
                    f"fit) for {pair}; they were "
                    "estimated by fitting to UNIFAC (Dortmund) predictions. Treat results for "
                    "this pair as predictive, not data-based."
                )
        return _ge_model(self.method, c, a, b, alpha, idx)

    def psats(self, T: float, zs: list[float]) -> list[float]:
        """Reference fugacities of the liquid: vapor pressures, the Henry constant in this
        liquid for Henry components, and for a dimerizing acid the monomer pressure of its
        saturated vapor."""
        if self.method == WATER_ICE:
            raise _water_ice_specification_error()
        if not any(self.henry):
            return [vp(T) for vp in self._liquid_vps]
        return list(self._liquid_template().to(T=T, P=101325.0, zs=list(zs)).Psats())

    def _liquid_template(self):
        f = self._flasher
        return f.liquid if hasattr(f, "liquid") else f.liquids[0]

    @property
    def formation_enthalpies(self) -> list[float] | None:
        """Ideal-gas formation enthalpies at 298.15 K, J/mol; None if any is unknown."""
        hfs = self.consts.Hfgs
        return None if any(h is None for h in hfs) else list(hfs)

    def missing_formation_data(self) -> list[str]:
        return [n for n, h in zip(self.names, self.consts.Hfgs) if h is None]

    def _h_reactive(self, H: float, zs) -> float | None:
        hfs = self.formation_enthalpies
        return None if hfs is None else H + sum(z * h for z, h in zip(zs, hfs))

    def phase_props(self, kind: str, T: float, P: float, zs: list[float]) -> PhaseProps:
        """Properties of a single ``"vapor"`` or ``"liquid"`` phase without a flash.

        Equations of state may return the vapor-like root for the liquid when no liquid
        root exists (e.g. far above the boiling point); callers should expect K close to 1.
        """
        def props(ph):
            H = float(ph.H())
            return PhaseProps(tuple(float(v) for v in ph.lnphis()), H, self._h_reactive(H, zs),
                              self._volume(ph.V(), zs))

        return self._single_phase(kind, T, P, zs, props)

    def phase_state(self, kind: str, T: float, P: float, zs: list[float]) -> PhaseState:
        """A :class:`PhaseState` for one ``"vapor"`` or ``"liquid"`` phase, without a flash.

        For products a unit has already separated (e.g. each layer of a decanter), so the
        product keeps exactly the phase and enthalpy it had inside the unit.
        """
        vapor = kind == "vapor"

        def state(ph):
            H = float(ph.H())
            return PhaseState(
                T=T, P=P, VF=1.0 if vapor else 0.0, H=H, S=float(ph.S()),
                phase="vapor" if vapor else "liquid", x=None if vapor else tuple(zs),
                y=tuple(zs) if vapor else None, V=self._volume(ph.V(), zs),
                H_reactive=self._h_reactive(H, zs),
            )

        return self._single_phase(kind, T, P, zs, state)

    def transport_properties(self, kind: str, T: float, P: float,
                             zs: list[float]) -> TransportProps:
        """Density, heat capacity, viscosity and thermal conductivity of one ``"vapor"`` or
        ``"liquid"`` phase at T, P and mole fractions ``zs``, for film coefficients.

        The density comes from the package's own phase model (the same as its flashes), and
        so does a vapor's heat capacity. A liquid's heat capacity is the mole-fraction
        average of the pure-liquid correlations, which is closer to measurements than an
        activity model's (built from the heat of vaporization), with the phase model as a
        fallback. Viscosity and thermal conductivity come from the pure-component
        correlations of ``thermo`` (DIPPR and others) mixed by fixed rules: liquids by
        logarithmic (molar) viscosity mixing and DIPPR 9H, gases by Brokaw and
        Lindsay-Bromley.
        """
        zs = list(zs)
        mws = self.mws
        mw = sum(z * m for z, m in zip(zs, mws)) / 1e3  # kg/mol
        ws = [z * m / 1e3 / mw for z, m in zip(zs, mws)]

        def basic(ph):
            return (mw / self._volume(ph.V(), zs), float(ph.Cp()) / mw)

        rho, cp = self._single_phase(kind, T, P, zs, basic)
        p = self.props
        vapor = kind == "vapor"
        if not vapor:
            try:
                cp_l = p.HeatCapacityLiquidMixture.calculate(T, P, zs, ws, "LINEAR")
            except (ArithmeticError, ValueError, TypeError):
                cp_l = None
            if cp_l is not None and math.isfinite(cp_l) and cp_l > 0.0:
                cp = float(cp_l) / mw
        mixtures = (((p.ViscosityGasMixture, "BROKAW"),
                     (p.ThermalConductivityGasMixture, "LINDSAY_BROMLEY")) if vapor
                    else ((p.ViscosityLiquidMixture, "Logarithmic mixing, molar"),
                          (p.ThermalConductivityLiquidMixture, "DIPPR_9H")))
        values = []
        for prop, (mixture, method) in zip(("viscosity", "thermal conductivity"), mixtures):
            try:
                v = mixture.calculate(T, P, zs, ws, method)
            except (ArithmeticError, ValueError, TypeError):  # out of range, missing data
                v = None
            if v is None or not math.isfinite(v) or v <= 0.0:
                names = ", ".join(n for n, z in zip(self.names, zs) if z > 0.0)
                raise SimulationError(
                    f"no {kind} {prop} for {names} at {T:.6g} K",
                    hint="the film coefficients need it; give the heat transfer "
                    "coefficient U (or that side's film coefficient) instead")
            values.append(float(v))
        return TransportProps(rho=rho, cp=cp, mu=values[0], k=values[1])

    def _single_phase(self, kind: str, T: float, P: float, zs: list[float], evaluate):
        if self.method == WATER_ICE:
            raise _water_ice_specification_error()
        if kind == "vapor":
            template = self._flasher.gas
        elif kind == "liquid":
            template = self._liquid_template()
        else:
            raise ValueError(f"kind must be 'vapor' or 'liquid', got {kind!r}")
        try:
            return evaluate(template.to(T=T, P=P, zs=list(zs)))
        except Exception as exc:
            raise SimulationError(
                f"property evaluation failed for the {kind} phase at T={T:.6g} K, P={P:.6g} Pa",
                hint="check that temperature and pressure are physically reasonable",
            ) from exc

    def ln_mixture_fugacity_coefficients(self, zs: list[float], T: float, P: float
                                         ) -> tuple[np.ndarray, PhaseState]:
        """ln Phi_i = ln(f_i / (z_i P)) of a mixture of overall composition ``zs`` at T, P.

        The mixture first settles into its phases (a T-P flash); f_i is then the same in each
        of them, so one phase gives it: f_i = y_i phi_i P from the vapor if there is vapor,
        else from the liquid. With an equation of state phi comes from it, with an activity
        model the liquid's phi is gamma_i Psat_i / P (times its corrections), so
        sum_i nu_i ln(z_i Phi_i P / P°) is the reaction's ln Q on fugacities. A component
        absent from the mixture (z_i = 0) gets the fugacity coefficient of that phase at
        infinite dilution. Also returns the flash.
        """
        st = self.flash(list(zs), T=T, P=P)
        if st.VF > 0.0 and st.y is not None:
            kind, comp = "vapor", st.y
        elif st.liquid_split and st.x_light is not None:
            kind, comp = "liquid", st.x_light
        else:
            kind, comp = "liquid", st.x if st.x is not None else tuple(zs)
        lnphi = np.asarray(self.phase_props(kind, T, P, list(comp)).lnphi, dtype=float)
        z = np.asarray(zs, dtype=float)
        c = np.asarray(comp, dtype=float)
        out = lnphi.copy()
        present = z > 0.0
        out[present] = (np.log(np.maximum(c[present], 1e-300)) + lnphi[present]
                        - np.log(z[present]))
        return out, st

    def k_values(self, T: float, P: float, x: list[float], y: list[float]) -> list[float]:
        """Equilibrium ratios K_i = y_i / x_i = phi_i(liquid) / phi_i(vapor)."""
        lnphi_l = self.phase_props("liquid", T, P, x).lnphi
        lnphi_v = self.phase_props("vapor", T, P, y).lnphi
        return [math.exp(a - b) for a, b in zip(lnphi_l, lnphi_v)]

    def flash(self, zs: list[float], **spec: float) -> PhaseState:
        """Flash at two of T, P, H, S, VF (SI, molar basis)."""
        zs = list(zs)
        if self.method == WATER_ICE:
            return self._water_ice_flash(zs, spec)
        if spec.keys() == {"P", "H"}:
            return self._flash_ph(zs, spec["P"], spec["H"])
        if self.n > 1 and spec.keys() == {"P", "S"}:
            boiling = self._pure_boiling(zs, spec)
            if boiling is not None:
                return boiling
        return self._state(self._flash_raw(zs, spec)[0], zs)

    def _water_ice_flash(self, zs: list[float], spec: dict[str, float]) -> PhaseState:
        """Equilibrium between one liquid and one solid, neglecting the melting curve's slope."""
        if spec.keys() not in ({"P", "T"}, {"P", "H"}):
            raise _water_ice_specification_error()
        if zs != [1.0]:
            raise SimulationError(
                "the water-ice reference method requires a pure water mole fraction of 1",
                hint="pass the normalized water composition",
            )
        P = spec["P"]
        if not 90_000.0 <= P <= 110_000.0:
            raise SimulationError(
                f"water-ice pressure {P:.6g} Pa is outside the reference range",
                hint="use a pressure between 0.9 and 1.1 bar",
            )

        T = spec.get("T")
        if T is not None and not _ICE_MIN_T <= T <= _ICE_MAX_T:
            raise SimulationError(
                f"water-ice temperature {T:.6g} K is outside the reference range",
                hint=f"use a temperature between {_ICE_MIN_T} and {_ICE_MAX_T} K",
            )

        tm = _ICE_MELTING_T
        # Anchor both phases to the existing ideal-gas reference via liquid water at melting.
        liquid = self._flasher.flash(T=tm, P=P, zs=zs)
        h_liquid, s_liquid, v_liquid = liquid.H(), liquid.S(), liquid.V()
        fusion = self.consts.Hfus_Tms[0]  # J/mol
        h_solid = h_liquid - fusion
        s_solid = s_liquid - fusion / tm
        v_solid = self.mws[0] / (1e3 * 917.0)  # m3/mol, ice Ih near melting
        cp_solid, cp_liquid = 37.8, 76.0  # J/(mol K), rounded near-melting values

        if T is None:
            H = spec["H"]
            minimum = h_solid + cp_solid * (_ICE_MIN_T - tm)
            maximum = h_liquid + cp_liquid * (_ICE_MAX_T - tm)
            if not minimum <= H <= maximum:
                raise SimulationError(
                    f"water-ice enthalpy {H:.6g} J/mol is outside the reference range",
                    hint=f"use an enthalpy between {minimum:.6g} and {maximum:.6g} J/mol",
                )
            if H <= h_solid:
                solid = 1.0
                T = max(_ICE_MIN_T, tm + (H - h_solid) / cp_solid)
            elif H >= h_liquid:
                solid = 0.0
                T = min(_ICE_MAX_T, tm + (H - h_liquid) / cp_liquid)
            else:
                T = tm
                solid = (h_liquid - H) / fusion
        elif T < tm:
            solid = 1.0
            H = h_solid + cp_solid * (T - tm)
        else:
            # At melting T and P do not fix the phase split. T-P returns liquid.
            solid = 0.0
            H = h_liquid + cp_liquid * (T - tm)

        if solid == 1.0:
            phase = "solid"
        elif solid == 0.0:
            phase = "liquid"
        else:
            phase = "liquid-solid"
        heat_capacity = solid * cp_solid + (1.0 - solid) * cp_liquid
        entropy = solid * s_solid + (1.0 - solid) * s_liquid + heat_capacity * math.log(T / tm)
        return PhaseState(
            T=T, P=P, VF=0.0, H=H, S=entropy, phase=phase, x=tuple(zs) if solid < 1.0 else None,
            y=None, V=solid * v_solid + (1.0 - solid) * v_liquid, H_reactive=self._h_reactive(H, zs),
            solid_fraction=solid, x_solid=tuple(zs) if solid > 0.0 else None,
        )

    def _flash_ph(self, zs: list[float], P: float, H: float) -> PhaseState:
        active_components = sum(z > 0.0 for z in zs)
        if self.n > 1 and active_components == 1:
            boiling = self._pure_boiling(zs, {"P": P, "H": H})
            if boiling is not None and _enthalpy_matches(boiling.H - H, H):
                return boiling

        standard = (
            lambda: self._flasher.flash(zs=zs, P=P, H=H),
            lambda: self._by_temperature(zs, P, H),
        )
        by_fraction = (lambda: self._ph_by_vapor_fraction(zs, P, H),)
        # Narrow boiling bands defeat the default temperature iteration. With
        # composition-independent K values the fraction search is especially fast,
        # so try it first for a dominant component without dropping any traces.
        narrow = (active_components > 1 and max(zs) > 1.0 - 1e-6
                  and getattr(self._flasher, "K_composition_independent", False))
        attempts = by_fraction + standard if narrow else standard + by_fraction

        closest_H = None
        last_error = None
        for flash in attempts:
            try:
                raw = flash()
                actual_H = float(raw.H())
                if _enthalpy_matches(actual_H - H, H):
                    return self._state(raw, zs)
                if math.isfinite(actual_H) and (
                    closest_H is None or abs(actual_H - H) < abs(closest_H - H)
                ):
                    closest_H = actual_H
            except Exception as exc:  # noqa: BLE001 - thermo raises several exception types
                last_error = exc

        details = "no enthalpy candidate available"
        if closest_H is not None:
            details = f"calculated H={closest_H:.9g} J/mol, residual={closest_H - H:+.9g} J/mol"
        raise SimulationError(
            f"pressure-enthalpy flash did not converge at P={P:.6g} Pa: "
            f"requested H={H:.9g} J/mol, {details}",
            hint="check the specified duty or pressure and whether the property method "
            "supports the phases at this condition",
        ) from last_error

    def _ph_by_vapor_fraction(self, zs: list[float], P: float, H: float):
        """Resolve latent heat when the equilibrium is too narrow for a temperature search.

        At fixed P, a vapor-fraction flash uses the actual composition of both phases,
        including traces. This avoids replacing a barely separated mixture with a pure
        compound or interpolating between unrelated T-P states across a phase boundary.
        """
        def at_fraction(vf: float):
            return self._flasher.flash(zs=zs, P=P, VF=vf)

        vf = brentq(lambda vf: float(at_fraction(vf).H()) - H, 0.0, 1.0, xtol=1e-12)
        result = at_fraction(vf)
        # A fixed-VF flash searches for vapor-liquid equilibrium. Reject its candidate
        # if a phase split missed by that search has lower Gibbs energy at the same T/P.
        if self.supports_lle:
            stable = self._flasher.flash(zs=zs, P=P, T=float(result.T))
            if result.G() > stable.G() + 1e-3:
                raise ValueError("vapor-fraction flash found a metastable equilibrium")
        return result

    def _pure_boiling(self, zs: list[float], spec: dict) -> PhaseState | None:
        """A single compound in a package of several (the others at zero flow) that boils at
        the given enthalpy or entropy, or None.

        A pure compound boils at one temperature, so its enthalpy fixes the vapor fraction,
        not the temperature. The multicomponent flash cannot see that: it snaps to all
        liquid or all vapor, losing up to the whole heat of vaporization, and takes seconds.
        Here the saturated liquid and vapor are evaluated at the boiling point and mixed.
        """
        if max(zs) < 1.0 - 1e-12 or sum(1 for z in zs if z > 0.0) == len(zs):
            return None
        P, by = spec["P"], "H" if "H" in spec else "S"
        try:
            T = float(self._flash_raw(zs, {"P": P, "VF": 0.0})[0].T)
            liquid = self.phase_state("liquid", T, P, zs)
            vapor = self.phase_state("vapor", T, P, zs)
        except (SimulationError, ValueError, ArithmeticError):
            return None
        lo, hi = getattr(liquid, by), getattr(vapor, by)
        if not lo < spec[by] < hi:
            return None  # one phase: the ordinary flash finds it
        vf = (spec[by] - lo) / (hi - lo)

        def mixed(a, b):
            return None if a is None or b is None else (1.0 - vf) * a + vf * b

        return PhaseState(
            T=T, P=P, VF=vf, H=mixed(liquid.H, vapor.H), S=mixed(liquid.S, vapor.S),
            phase="two-phase", x=tuple(zs), y=tuple(zs), V=mixed(liquid.V, vapor.V),
            H_reactive=mixed(liquid.H_reactive, vapor.H_reactive))

    def saturation_point(self, zs: list[float], kind: str, *, T: float | None = None,
                         P: float | None = None, guess: float | None = None) -> Saturation:
        """Bubble (``kind="bubble"``) or dew point at fixed T or P, with the incipient phase.

        ``guess`` (the unknown T or P, e.g. the neighbouring point of a diagram) speeds up
        mixtures that may split into two liquids.
        """
        if self.method == WATER_ICE:
            raise _water_ice_specification_error()
        spec = {"VF": 0.0 if kind == "bubble" else 1.0, **({"T": T} if T is not None else {}),
                **({"P": P} if P is not None else {})}
        res, incipient = self._flash_raw(list(zs), spec, guess)
        other = incipient or tuple(zs)
        if kind == "bubble":
            return Saturation(float(res.T), float(res.P), tuple(zs), other)
        return Saturation(float(res.T), float(res.P), other, tuple(zs))

    def _flash_raw(self, zs: list[float], spec: dict, guess: float | None = None):
        """The flash result, plus the incipient phase composition at a VF = 0 or 1 boundary."""
        try:
            if "VF" in spec and self.n > 1 and (self.supports_lle or self._ge is not None):
                res, incipient = self._checked_vf(zs, spec, guess)
            else:
                try:
                    res = self._flasher.flash(zs=zs, **spec)
                except Exception:
                    if spec.keys() != {"P", "S"}:
                        raise
                    res = self._by_temperature(zs, spec["P"], spec["S"], "S")
                incipient = _incipient(res, spec.get("VF"), self.n)
        except Exception as exc:
            spec_text = ", ".join(f"{k}={v:.6g}" for k, v in spec.items())
            raise SimulationError(
                f"phase equilibrium calculation did not converge ({spec_text})",
                hint="check that temperature and pressure are physically reasonable for "
                "these compounds, or try a different property method",
            ) from exc
        if spec.get("VF") in (0.0, 1.0) and self.n > 1:
            self._check_saturation(zs, spec, res, incipient)
        return res, incipient

    def _check_saturation(self, zs: list[float], spec: dict, res, incipient) -> None:
        """Refuse a bubble or dew point that is not one.

        The solvers in thermo can converge to the other boundary (a "bubble point" whose
        vapor has the feed composition is a dew point), which looks like a valid number.
        An incipient phase equal to the feed is no sign of failure: it is right at an
        azeotrope and for a pure component.
        """
        def is_feed(phase_zs) -> bool:
            return phase_zs is not None and max(abs(a - b) for a, b in zip(phase_zs, zs)) < 1e-3

        bubble = spec["VF"] == 0.0
        gas = getattr(res, "gas", None)
        liquids = list(getattr(res, "liquids", None) or [])
        y = tuple(gas.zs) if gas is not None else None
        x = _combined_liquid(res, liquids, self.n) if liquids else None
        # Only the swap is an error: next to a liquid that splits, the liquid phase shown
        # at the bubble point need not have the feed composition.
        requested, other = (x, y) if bubble else (y, x)
        name = "bubble" if bubble else "dew"
        at = ", ".join(f"{k}={v:.6g}" for k, v in spec.items() if k != "VF")
        if is_feed(other) and not is_feed(requested):
            raise SimulationError(
                f"no {name} point found ({at}): the calculation converged to the "
                f"{'dew' if bubble else 'bubble'} point instead",
                hint="a mixture rich in a supercritical gas such as hydrogen may have no "
                f"{name} point at this condition; specify the temperature instead",
            )

    def flash_ph_near(self, zs: list[float], P: float, H: float, T_guess: float) -> PhaseState:
        """Pressure-enthalpy flash from a nearby temperature, by T-P flashes.

        For marching models (an evaporator film, a pipe) where the previous step's T is a
        good start: secant steps on H(T), kept inside a bracket, need a few T-P flashes,
        much less than a full PH flash. Falls back to ``flash`` if they do not settle.
        """
        if self.method == WATER_ICE:
            return self.flash(zs, P=P, H=H)

        def h(T: float):
            res = self._flasher.flash(zs=zs, T=T, P=P)
            return float(res.H()) - H, res

        try:
            T0, (f0, res) = T_guess, h(T_guess)
            if _enthalpy_matches(f0, H):
                return self._state(res, zs)
            T1 = T0 + (1.0 if f0 < 0.0 else -1.0)
            lo = hi = None
            for _ in range(40):
                f1, res = h(T1)
                if _enthalpy_matches(f1, H):
                    return self._state(res, zs)
                if f1 < 0.0:
                    lo = T1 if lo is None or T1 > lo else lo
                else:
                    hi = T1 if hi is None or T1 < hi else hi
                step = -f1 * (T1 - T0) / (f1 - f0) if f1 != f0 else (1.0 if f1 < 0 else -1.0)
                T_next = T1 + max(min(step, 200.0), -200.0)
                if lo is not None and hi is not None and not lo < T_next < hi:
                    T_next = 0.5 * (lo + hi)
                T0, f0, T1 = T1, f1, max(T_next, 20.0)
        except Exception:  # noqa: BLE001, S110 - any failure: use the general flash
            pass
        return self.flash(zs, P=P, H=H)

    def _by_temperature(self, zs: list[float], P: float, value: float, by: str = "H"):
        """Pressure-enthalpy (or pressure-entropy, ``by="S"``) flash by solving for T with
        T-P flashes.

        The fallback when thermo's own flash fails, e.g. for an almost pure component in a
        mixture model (the residue at the end of an evaporator), or an isentropic expansion
        of a mixture into the two-phase region with an activity-coefficient method. H and S
        rise monotonically with T at fixed P, and T-P flashes stay robust there.
        """
        def h(T: float):
            res = self._flasher.flash(zs=zs, T=T, P=P)
            return float(res.H() if by == "H" else res.S()) - value, res

        lo = hi = sum(z * tb for z, tb in zip(zs, self.consts.Tbs) if tb) or 300.0
        f_lo = f_hi = h(lo)[0]
        step = 20.0
        while f_lo > 0.0:
            lo = max(lo - step, 20.0)
            f_lo = h(lo)[0]
            step *= 1.6
            if lo <= 20.0 and f_lo > 0.0:
                raise ValueError("specification below the range of T-P flashes")
        step = 20.0
        while f_hi < 0.0:
            hi += step
            f_hi = h(hi)[0]
            step *= 1.6
            if hi > 5000.0:
                raise ValueError("specification above the range of T-P flashes")
        T = brentq(lambda T: h(T)[0], lo, hi, xtol=1e-7, rtol=1e-12)
        return h(T)[1]

    def _checked_vf(self, zs: list[float], spec: dict, guess: float | None):
        """A vapor-fraction flash that agrees with T-P flashes, which weigh liquid splits.

        The vapor-fraction solvers in thermo look for vapor-liquid solutions only, so for a
        heterogeneous mixture they can return a metastable bubble point (for water/hexane
        about 10 K too low), and they are slow there. The vapor fraction rises monotonically
        with T (falls with P), so a first estimate is checked by T-P flashes on both sides
        and, if it fails, the boundary is bracketed and bisected.
        """
        target, by_T = spec["VF"], "P" in spec
        res = y_bubble = None
        if guess is not None:
            u0 = guess if by_T else -math.log(guess)
        elif target == 0.0 and (hom := self._homogeneous_bubble(zs, spec)) is not None:
            u0, y_bubble = hom
        else:
            try:
                res = self._flasher.flash(zs=zs, **spec)
                u0 = float(res.T) if by_T else -math.log(float(res.P))
            except (UnconvergedError, OscillationError, ValueError, ArithmeticError,
                    UnboundLocalError):  # thermo can fail outright on immiscible pairs
                u0 = (sum(z * tb for z, tb in zip(zs, self.consts.Tbs)) if by_T else
                      -math.log(sum(z * p for z, p in zip(zs, self.psats(spec["T"], zs)))))

        def tp(u):
            T, P = (u, spec["P"]) if by_T else (spec["T"], math.exp(-u))
            return self._flasher.flash(zs=zs, T=T, P=P)

        def above(r):  # past the requested vapor fraction (raising u = T or -ln P adds vapor)
            vf = _vapor_fraction(r)
            return vf > 1e-9 if target == 0.0 else vf >= target - 1e-7

        tol = 1e-3 if by_T else 1e-6
        r_lo = tp(u0 - tol)
        if y_bubble is not None and len(getattr(r_lo, "liquids", None) or []) < 2:
            # The estimate is exact unless the boiling liquid splits; near an azeotrope the
            # two-phase band is too narrow for a T-P flash to be a reliable judge.
            return r_lo, y_bubble
        r_hi = tp(u0 + tol)
        if above(r_hi) and not above(r_lo):
            if res is not None:
                return res, _incipient(res, target, self.n)
            u_lo = u_hi = u0
        elif above(r_lo):  # the boundary lies below u0
            u_hi, step = u0 - tol, 0.5 if by_T else 5e-3
            u_lo = u_hi - step
            while above(r_lo := tp(u_lo)):
                u_hi, r_hi, step = u_lo, r_lo, step * 2
                u_lo = u_hi - step
        else:
            u_lo, step = u0 + tol, 0.5 if by_T else 5e-3
            u_hi = u_lo + step
            while not above(r_hi := tp(u_hi)):
                u_lo, r_lo, step = u_hi, r_hi, step * 2
                u_hi = u_lo + step
        while u_hi - u_lo > tol:
            mid = 0.5 * (u_lo + u_hi)
            r = tp(mid)
            if above(r):
                u_hi, r_hi = mid, r
            else:
                u_lo, r_lo = mid, r
        if target == 0.0:
            return r_lo, _incipient(r_hi, target, self.n)
        if target == 1.0:
            return r_hi, _incipient(r_lo, target, self.n)
        return min((r_lo, r_hi), key=lambda r: abs(_vapor_fraction(r) - target)), None

    def _homogeneous_bubble(self, zs: list[float], spec: dict):
        """Bubble point of the unsplit liquid by modified Raoult's law: (T or -ln P, vapor y).

        Exact for the activity models here (ideal vapor of monomers and acid dimers, no
        Poynting correction) as long as the liquid does not split, and far cheaper than a
        flash.
        """
        ge, vps = self._ge, self.props.VaporPressures
        if ge is None:
            return None

        def vapor(T):
            """Total pressure and vapor composition over the liquid."""
            gammas = ge.to_T_xs(T, list(zs)).gammas()
            pp = [z * g * p for z, g, p in zip(zs, gammas, self.psats(T, zs))]
            if self._assoc:  # the acids' monomers, plus their dimers
                return vapor_from_monomers(self._assoc, T, pp)
            total = sum(pp)
            return total, [p / total for p in pp]

        try:
            if "T" in spec:
                T = spec["T"]
            else:
                P = spec["P"]
                tsat = [vp.solve_property(P) for vp, h in zip(vps, self.henry) if not h]
                T = brentq(lambda T: math.log(vapor(T)[0] / P), min(tsat) - 80, max(tsat) + 80)
            total, y = vapor(T)
            return (T if "P" in spec else -math.log(total)), tuple(y)
        except (ValueError, ArithmeticError):
            return None

    def _state(self, res, zs) -> PhaseState:
        vf = _vapor_fraction(res)
        gas = getattr(res, "gas", None)
        liquids = list(getattr(res, "liquids", None) or [])
        liquid_x = _combined_liquid(res, liquids, self.n) if liquids else None
        split: dict = {}
        if len(liquids) > 1:
            # thermo's betas_liquids are fractions of the liquid; the split is per mole of
            # mixture, so with vapor present they are scaled by the liquid fraction.
            betas = [float(b) * (1.0 - vf) for b in res.betas_liquids]
            by_density = sorted(zip(liquids, betas), key=lambda pb: pb[0].rho_mass())
            (light, _), (heavy, beta_heavy) = by_density[0], by_density[-1]
            split = {"x_light": tuple(light.zs), "x_heavy": tuple(heavy.zs),
                     "heavy_fraction": beta_heavy, "rho_light": float(light.rho_mass()),
                     "rho_heavy": float(heavy.rho_mass())}
        phase = "two-phase" if 0.0 < vf < 1.0 else ("vapor" if vf >= 1.0 else "liquid")
        return PhaseState(
            T=float(res.T),
            P=float(res.P),
            VF=vf,
            H=float(res.H()),
            S=float(res.S()),
            phase=phase,
            x=liquid_x if liquid_x is not None and vf < 1.0 else None,
            y=tuple(gas.zs) if gas is not None and vf > 0.0 else None,
            V=self._volume(res.V(), zs),
            H_reactive=self._h_reactive(float(res.H()), zs),
            **split,
        )


def _enthalpy_matches(residual: float, target: float) -> bool:
    # Same molar tolerance as the nearby-temperature iteration: an absolute floor
    # at the enthalpy reference, plus a relative term for states far from it.
    return (math.isfinite(target) and math.isfinite(residual)
            and abs(residual) <= 1e-3 + 1e-7 * abs(target))


def _vapor_fraction(res) -> float:
    vf = float(res.VF) if res.VF is not None else (1.0 if res.phase == "V" else 0.0)
    return min(max(vf, 0.0), 1.0)


def _incipient(res, vf: float | None, n: int) -> tuple[float, ...] | None:
    """The phase about to appear at a bubble (vf = 0: vapor) or dew point (vf = 1: liquid)."""
    if vf == 0.0:
        gas = getattr(res, "gas", None)
        return tuple(gas.zs) if gas is not None else None
    if vf == 1.0:
        liquids = list(getattr(res, "liquids", None) or [])
        return _combined_liquid(res, liquids, n) if liquids else None
    return None


def _combined_liquid(res, liquids, n: int) -> tuple[float, ...]:
    """Composition of all liquid phases together; with one liquid, simply its composition."""
    if len(liquids) == 1:
        return tuple(liquids[0].zs)
    betas = [float(b) for b in res.betas_liquids]
    return tuple(sum(b * ph.zs[i] for b, ph in zip(betas, liquids)) / sum(betas)
                 for i in range(n))


# Which equation-of-state root counts as vapor: V T^2 above Vpc Tpc^2, with Kay's
# pseudo-critical volume and temperature. thermo's default, the phase identification
# parameter, calls a compressed supercritical gas a liquid (hydrogen at 25 °C and 70 bar),
# so a hydrogen-rich gas over toluene came out as two liquids with no vapor at all. A
# volume-only rule fails the same way at higher pressure (hydrogen at 130 bar, nitrogen at
# 340 bar), a temperature-only rule calls superheated benzene vapor a liquid. This one
# passes the pure-component phase checks in tests/test_phase_physics.py.
_SETTINGS = BulkSettings(VL_ID="Tpc Vpc")


class _FlashVLLE(FlashVLN):
    """FlashVLN that also weighs a liquid-liquid split against a vapor-liquid answer.

    For strongly immiscible pairs (water/hexane) the T-P flash in thermo can accept a
    vapor-liquid solution while two liquids have the lower Gibbs energy. Checking that
    candidate too puts the three-phase temperature where the model says it is.
    """

    def flash_TPV(self, T, P, V, zs=None, solution=None, hot_start=None, solution_target=None):
        sln = super().flash_TPV(T, P, V, zs, solution, hot_start, solution_target)
        if T is None or P is None or sln[0] is None or len(sln[1]) != 1:
            return sln
        try:
            ll = self.flash_TP_stability_test(T, P, zs, self.liquids[0], self.liquids[1], LL=True)
        except (UnconvergedError, OscillationError, ValueError, ArithmeticError):
            return sln
        if ll[0] is None and len(ll[1]) == 2 and self._gibbs(ll) < self._gibbs(sln):
            return ll
        return sln

    def _gibbs(self, sln) -> float:
        phases = ([sln[0]] if sln[0] is not None else []) + list(sln[1])
        return sum(b * (p.G_min_criteria() if self.ideal_gas_basis else p.G())
                   for b, p in zip(sln[3], phases))


@lru_cache(maxsize=16)
def get_package(components: tuple[str, ...], method: str, prefer_fitted: bool = True,
                henry_components: tuple[str, ...] | None = None,
                custom: tuple[CustomCompound, ...] = (), *,
                binary_parameters: tuple[BinaryParameterModel, ...] = (),
                group_assignments: tuple[GroupAssignment, ...] = ()) -> PropertyPackage:
    """A cached package; ``custom`` holds the full definitions, so editing a user-defined
    compound under the same name builds a new one."""
    return PropertyPackage(list(components), method, prefer_fitted, henry_components, custom,
                           binary_parameters=binary_parameters, group_assignments=group_assignments)


def package_options(model) -> dict:
    """The flowsheet's thermodynamic options as keyword arguments for :func:`get_package`."""
    henry = getattr(model, "henry_components", None)
    custom = getattr(model, "custom_compounds", None) or []
    data = getattr(model, "property_data", None)
    return {"henry_components": tuple(henry) if henry is not None else None,
            "custom": tuple(CustomCompound.from_dict(c.model_dump(exclude_none=True))
                            for c in custom),
            **({"binary_parameters": tuple(p for p in data.binary_parameters
                                           if all(c in model.components for c in p.components)),
                "group_assignments": tuple(GroupAssignment(g.component,
                                                          tuple(sorted(g.groups.items())), g.source)
                                           for g in data.unifac_groups
                                           if g.component in model.components)} if data else {})}


def default_package(model) -> tuple[PropertyPackage, str | None]:
    options = package_options(model)
    method, reason = model.property_method, None
    if method == "auto":
        suggestion = suggest_method(model.components, **{k: v for k, v in options.items()
                                                        if k != "henry_components"})
        method, reason = suggestion.method, suggestion.reason
    elif method not in METHODS:
        raise SimulationError(f"unknown property method {method!r}",
                              hint=f"use 'auto' or one of: {', '.join(METHODS)}")
    return get_package(tuple(model.components), method, **options), reason


plugins.load("property_methods", _register_methods)
plugins.load("components", _register_compounds)
