"""Property packages: the thermodynamic model behind every flash calculation.

A property package bundles the compound list with a phase-equilibrium method
(equation of state or activity-coefficient model). It is built once per
flowsheet and cached, because loading the constants takes a moment.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from fluids.numerics import OscillationError, UnconvergedError
from scipy.optimize import brentq, least_squares
from thermo import (
    NRTL,
    PRMIX,
    SRKMIX,
    UNIFAC,
    UNIQUAC,
    CEOSGas,
    CEOSLiquid,
    ChemicalConstantsPackage,
    FlashPureVLS,
    FlashVL,
    FlashVLN,
    GibbsExcessLiquid,
    IdealGas,
    Wilson,
)
from thermo.bulk import BulkSettings
from thermo.interaction_parameters import IPDB
from thermo.unifac import DOUFIP2016, DOUFSG

from . import custom_compounds
from .compounds import Compound, lookup
from .custom_compounds import CustomCompound
from .errors import SimulationError
from .fitting import store as fitted_store
from .henry import (
    HENRY_TC,
    HenryDilutedModel,
    HenryLiquid,
    henry_constant,
    ideal_solvent_model,
    raoult_estimate,
)


@dataclass(frozen=True)
class Method:
    key: str
    label: str
    description: str


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
    )
}


@dataclass(frozen=True)
class PhaseState:
    """Result of a flash: everything per mole of mixture, SI units."""

    T: float
    P: float
    VF: float
    H: float  # J/mol, relative to ideal gas at 298.15 K (no formation enthalpy)
    S: float  # J/mol/K
    phase: str  # "vapor", "liquid" or "two-phase"
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

    @property
    def liquid_split(self) -> bool:
        return self.x_heavy is not None

    @property
    def light_fraction(self) -> float:
        return max(1.0 - self.VF - self.heavy_fraction, 0.0) if self.liquid_split else 0.0


@dataclass(frozen=True)
class PhaseProps:
    """One phase at fixed T, P and composition (no phase split). For stage-wise models."""

    lnphi: tuple[float, ...]  # log fugacity coefficients
    H: float  # J/mol
    H_reactive: float | None  # J/mol, see PhaseState
    V: float  # m3/mol


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


def suggest_method(components: list[str],
                   custom: tuple[CustomCompound, ...] = ()) -> Suggestion:
    """Pick a sensible default method and say why, the way an experienced engineer would."""
    consts = _constants(tuple(c.cas for c in resolve(components, custom)))[0]
    if len(components) == 1:
        return Suggestion(
            "ideal",
            f"A single compound has no mixing effects, so the ideal method is exact for phase "
            f"equilibrium: it uses {consts.names[0]}'s own vapor-pressure correlation. Switch to "
            "Peng-Robinson for high pressures or near the critical point.",
        )
    polar = [
        name
        for name, dipole in zip(consts.names, consts.dipoles)
        if dipole is not None and dipole >= 1.0
    ]
    splits = _liquid_liquid_pairs(consts, polar)
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
    """The compounds behind component names; a custom definition wins over the database."""
    defined = {c.name.lower(): c for c in custom}
    out = []
    for name in components:
        c = defined.get(name.strip().lower())
        if c is None:
            out.append(lookup(name))
            continue
        est = custom_compounds.estimate(c)
        out.append(Compound(cas=custom_compounds.register(c), name=c.name, formula=c.formula,
                            mw=est.values["MW"], smiles=c.smiles, simulation_ready=True))
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
    model = _unifac(_constants((cas_a, cas_b))[0], T)
    if model is None:
        return False
    g = []
    for x in _X_GRID:
        gammas = model.to_T_xs(T, [x, 1.0 - x]).gammas()
        g.append(x * math.log(x * gammas[0]) + (1 - x) * math.log((1 - x) * gammas[1]))
    return any(g[i - 1] - 2 * g[i] + g[i + 1] < -1e-9 for i in range(1, len(g) - 1))


def _liquid_liquid_pairs(consts, polar: list[str]) -> list[tuple[str, str, str, str]]:
    """Binary pairs (names, CAS numbers) likely to split; only pairs with a polar member."""
    out = []
    for i in range(consts.N):
        for j in range(i + 1, consts.N):
            a, b = consts.names[i], consts.names[j]
            if (a in polar or b in polar) and _binary_splits(consts.CASs[i], consts.CASs[j]):
                out.append((a, b, consts.CASs[i], consts.CASs[j]))
    return out


@lru_cache(maxsize=128)
def _estimate_binary(method: str, cas_a: str, cas_b: str):
    """Fit binary parameters to UNIFAC activity coefficients when no measured set exists.

    Returns ``(a12, a21, b12, b21, alpha)`` or None when UNIFAC cannot describe the pair.
    The fit covers the whole composition range at three temperatures spanning the two
    boiling points, so it reproduces UNIFAC's predictions, not measured data.
    """
    consts = _constants((cas_a, cas_b))[0]
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
    alpha_value = 0.2 if _binary_splits(cas_a, cas_b) else 0.3
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


class PropertyPackage:
    def __init__(self, components: list[str], method: str, prefer_fitted: bool = True,
                 henry_components: tuple[str, ...] | None = None,
                 custom: tuple[CustomCompound, ...] = ()):
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
        self.compounds = resolve(list(components), self.custom)
        self.names = list(components)
        cas = tuple(c.cas for c in self.compounds)
        if len(set(cas)) != len(cas):
            raise SimulationError(
                "the component list contains the same compound twice",
                hint="different names can refer to the same compound (e.g. 'ethanol' "
                "and 'ethyl alcohol'); remove one",
            )
        self.consts, self.props = _constants(cas)
        # Where each binary parameter set came from, and warnings worth showing with results.
        self.parameter_sources: list[dict] = []
        self.warnings: list[str] = []
        self._ge = None  # activity-coefficient model, if the method has one
        self._report_estimates()
        self._check_data()
        # Henry components (dissolved gases) of activity methods: flags and names.
        self._henry_request = henry_components
        self.henry = [False] * self.n
        self._flasher = self._build()

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
                    "pair": comp.name, "source": f"{prop}: {how}",
                    **({"value": value} if isinstance(value, float) else {})})
            self.warnings.append(
                f"{comp.name} is a user-defined compound; estimated: {est.summary()}. Treat "
                "results that depend on it as estimates."
                + "".join(f" Note: {note}." for note in est.notes))

    @property
    def n(self) -> int:
        return len(self.names)

    @property
    def mws(self) -> list[float]:
        return list(self.consts.MWs)  # g/mol

    @property
    def supports_lle(self) -> bool:
        """Whether flashes can find two liquid phases (vapor-liquid-liquid equilibrium)."""
        return self.method in LLE_METHODS and self.n > 1

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
        if self.method in ("peng-robinson", "srk"):
            eos = PRMIX if self.method == "peng-robinson" else SRKMIX
            kijs = (
                ip_matrix("ChemSep PR", c.CASs, "kij")
                if self.method == "peng-robinson"
                else [[0.0] * self.n for _ in range(self.n)]
            )
            kwargs = {"Tcs": c.Tcs, "Pcs": c.Pcs, "omegas": c.omegas, "kijs": kijs}
            gas = CEOSGas(eos, kwargs, HeatCapacityGases=cp_gas)
            liquid = CEOSLiquid(eos, kwargs, HeatCapacityGases=cp_gas)
        else:
            gas = IdealGas(HeatCapacityGases=cp_gas)
            henry_kwargs: dict = {}
            if self.method in HENRY_METHODS and self.n > 1:
                self.henry, henry_kwargs = self._henry_setup()
            solvents = [i for i in range(self.n) if not self.henry[i]]
            ge_model = None
            if self.method == "unifac":
                for i in solvents:
                    if not c.UNIFAC_Dortmund_groups[i]:
                        raise SimulationError(
                            f"{self.names[i]!r} has no UNIFAC group assignment",
                            hint="use Peng-Robinson for this mixture instead"
                            + (" (a user-defined compound needs unifac_groups)"
                               if custom_compounds.is_custom(c.CASs[i]) else ""),
                        )
                ge_model = _unifac(c, idx=solvents)
            elif self.method in _GE_TABLES and len(solvents) > 1:
                ge_model = self._regressed_model(solvents)
            if henry_kwargs:
                ge_model = HenryDilutedModel(ge_model or ideal_solvent_model(len(solvents)),
                                             solvents, self.n)
            self._ge = ge_model
            liquid = (HenryLiquid if henry_kwargs else GibbsExcessLiquid)(
                VaporPressures=p.VaporPressures,
                HeatCapacityGases=cp_gas,
                VolumeLiquids=p.VolumeLiquids,
                GibbsExcessModel=ge_model,
                equilibrium_basis="Psat",
                caloric_basis="Psat",
                **henry_kwargs,
            )
        if self.n == 1:
            return FlashPureVLS(c, p, gas=gas, liquids=[liquid], solids=[], settings=_SETTINGS)
        if self.supports_lle:
            return _FlashVLLE(c, p, liquids=[liquid, liquid], gas=gas, settings=_SETTINGS)
        return FlashVL(c, p, liquid=liquid, gas=gas, settings=_SETTINGS)

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
                self.warnings.append(
                    f"{label}: no Henry constant for {self.names[i]} in "
                    f"{', '.join(self.names[j] for j in missing)}; its liquid fugacity uses its "
                    "vapor pressure correlation, extrapolated above the critical point. Its "
                    "solubility is not reliable; use Peng-Robinson for this gas.")
                continue
            flags[i] = True
            for j, k in found.items():
                coeffs[i][j] = k.coeffs
                self.parameter_sources.append({"pair": f"{self.names[i]} in {self.names[j]}",
                                               "source": f"Henry constant: {k.source}"})
            if missing:
                guess = raoult_estimate(self.props.VaporPressures[i])
                for j in missing:
                    coeffs[i][j] = guess.coeffs
                    self.parameter_sources.append({"pair": f"{self.names[i]} in {self.names[j]}",
                                                   "source": f"Henry constant: {guess.source}"})
                self.warnings.append(
                    f"{label}: no Henry constant for {self.names[i]} in "
                    f"{', '.join(self.names[j] for j in missing)}; it is estimated from "
                    f"{self.names[i]}'s extrapolated vapor pressure (ideal solution), which "
                    "can be off by a factor of several. Its solubility in liquids rich in "
                    f"{'that solvent' if len(missing) == 1 else 'those solvents'} is rough.")
        if not any(flags):
            return flags, {}
        names = ", ".join(self.henry_names(flags))
        self.warnings.append(
            f"{label}: {names} {'is a Henry component' if sum(flags) == 1 else 'are Henry components'}"
            " (dissolved gas: Henry's law, ideal dilute solution, no pressure correction). "
            "Valid for dilute solutions at moderate pressure; set henry_components to change.")
        kwargs = {"henry_components": flags}
        for k, key in enumerate(("henry_as", "henry_bs", "henry_cs", "henry_ds", "henry_es",
                                 "henry_fs")):
            kwargs[key] = [[coeffs[i][j][k] for j in range(n)] for i in range(n)]
        return flags, kwargs

    def henry_names(self, flags=None) -> list[str]:
        return [name for name, h in zip(self.names, flags or self.henry) if h]

    def _regressed_model(self, idx: list[int]):
        """NRTL, UNIQUAC or Wilson with measured parameters; UNIFAC estimates fill the gaps.

        Per pair: Fluxsheet's ThermoML fit (NRTL, UNIQUAC; first if ``prefer_fitted``), the
        ChemSep set, a ThermoML fit, then a UNIFAC estimate. Only the compounds at ``idx``
        (the solvents: Henry components are left out).
        """
        c, n, table = self.consts, self.n, _GE_TABLES[self.method]
        if self.method == "uniquac":
            for i in idx:
                if c.UNIFAC_Rs[i] is None or c.UNIFAC_Qs[i] is None:
                    raise SimulationError(
                        f"{self.names[i]!r} has no UNIQUAC size and surface parameters (r, q)",
                        hint="they come from UNIFAC groups, which this compound lacks; "
                        "use NRTL or Peng-Robinson instead",
                    )
        b = ip_matrix(table, c.CASs, "bij")
        a = (ip_matrix(table, c.CASs, "aij") if self.method == "wilson"
             else [[0.0] * n for _ in range(n)])
        alpha = (ip_matrix(table, c.CASs, "alphaij")
                 if self.method == "nrtl" else None)
        label = METHODS[self.method].label
        for i in idx:
            for j in (j for j in idx if j > i):
                pair = f"{self.names[i]}/{self.names[j]}"
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
                est = _estimate_binary(self.method, c.CASs[i], c.CASs[j])
                if est is None:
                    self.parameter_sources.append({"pair": pair, "source": "none (ideal)"})
                    self.warnings.append(
                        f"{label}: no binary parameters for {pair}, and UNIFAC cannot estimate "
                        "them; this pair is treated as an ideal solution. Use a method with "
                        "data for it, or Peng-Robinson for non-polar mixtures."
                    )
                    continue
                a[i][j], a[j][i], b[i][j], b[j][i], alpha_value = est
                if alpha is not None:
                    alpha[i][j] = alpha[j][i] = alpha_value
                self.parameter_sources.append({"pair": pair, "source": "UNIFAC estimate"})
                self.warnings.append(
                    f"{label}: there are no measured binary parameters (ChemSep or ThermoML "
                    f"fit) for {pair}; they were "
                    "estimated by fitting to UNIFAC (Dortmund) predictions. Treat results for "
                    "this pair as predictive, not data-based."
                )
        return _ge_model(self.method, c, a, b, alpha, idx)

    def psats(self, T: float, zs: list[float]) -> list[float]:
        """Vapor pressures, with the Henry constant in this liquid for Henry components."""
        if not any(self.henry):
            return [vp(T) for vp in self.props.VaporPressures]
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
                              float(ph.V()))

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
                y=tuple(zs) if vapor else None, V=float(ph.V()),
                H_reactive=self._h_reactive(H, zs),
            )

        return self._single_phase(kind, T, P, zs, state)

    def _single_phase(self, kind: str, T: float, P: float, zs: list[float], evaluate):
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

    def k_values(self, T: float, P: float, x: list[float], y: list[float]) -> list[float]:
        """Equilibrium ratios K_i = y_i / x_i = phi_i(liquid) / phi_i(vapor)."""
        lnphi_l = self.phase_props("liquid", T, P, x).lnphi
        lnphi_v = self.phase_props("vapor", T, P, y).lnphi
        return [math.exp(a - b) for a, b in zip(lnphi_l, lnphi_v)]

    def flash(self, zs: list[float], **spec: float) -> PhaseState:
        """Flash at two of T, P, H, S, VF (SI, molar basis)."""
        return self._state(self._flash_raw(list(zs), spec)[0], zs)

    def saturation_point(self, zs: list[float], kind: str, *, T: float | None = None,
                         P: float | None = None, guess: float | None = None) -> Saturation:
        """Bubble (``kind="bubble"``) or dew point at fixed T or P, with the incipient phase.

        ``guess`` (the unknown T or P, e.g. the neighbouring point of a diagram) speeds up
        mixtures that may split into two liquids.
        """
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
                    if set(spec) != {"P", "H"}:
                        raise
                    res = self._ph_by_temperature(zs, spec["P"], spec["H"])
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
        def h(T: float):
            res = self._flasher.flash(zs=zs, T=T, P=P)
            return float(res.H()) - H, res

        try:
            T0, (f0, res) = T_guess, h(T_guess)
            if abs(f0) <= 1e-7 * abs(H) + 1e-3:
                return self._state(res, zs)
            T1 = T0 + (1.0 if f0 < 0.0 else -1.0)
            lo = hi = None
            for _ in range(40):
                f1, res = h(T1)
                if abs(f1) <= 1e-7 * abs(H) + 1e-3:
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

    def _ph_by_temperature(self, zs: list[float], P: float, H: float):
        """Pressure-enthalpy flash by solving for T with T-P flashes.

        The fallback when thermo's PH flash fails, e.g. for an almost pure component in a
        mixture model (the residue at the end of an evaporator). H rises monotonically with
        T at fixed P, and T-P flashes stay robust there.
        """
        def h(T: float):
            res = self._flasher.flash(zs=zs, T=T, P=P)
            return float(res.H()) - H, res

        lo = hi = sum(z * tb for z, tb in zip(zs, self.consts.Tbs) if tb) or 300.0
        f_lo = f_hi = h(lo)[0]
        step = 20.0
        while f_lo > 0.0:
            lo = max(lo - step, 20.0)
            f_lo = h(lo)[0]
            step *= 1.6
            if lo <= 20.0 and f_lo > 0.0:
                raise ValueError("enthalpy below the range of T-P flashes")
        step = 20.0
        while f_hi < 0.0:
            hi += step
            f_hi = h(hi)[0]
            step *= 1.6
            if hi > 5000.0:
                raise ValueError("enthalpy above the range of T-P flashes")
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

        Exact for the activity models here (ideal vapor, no Poynting correction) as long as
        the liquid does not split, and far cheaper than a flash.
        """
        ge, vps = self._ge, self.props.VaporPressures
        if ge is None:
            return None

        def partial_pressures(T):
            gammas = ge.to_T_xs(T, list(zs)).gammas()
            return [z * g * p for z, g, p in zip(zs, gammas, self.psats(T, zs))]

        try:
            if "T" in spec:
                T = spec["T"]
            else:
                P = spec["P"]
                tsat = [vp.solve_property(P) for vp, h in zip(vps, self.henry) if not h]
                T = brentq(lambda T: math.log(sum(partial_pressures(T)) / P),
                           min(tsat) - 80, max(tsat) + 80)
            pp = partial_pressures(T)
            total = sum(pp)
            return (T if "P" in spec else -math.log(total)), tuple(p / total for p in pp)
        except (ValueError, ArithmeticError):
            return None

    def _state(self, res, zs) -> PhaseState:
        vf = _vapor_fraction(res)
        gas = getattr(res, "gas", None)
        liquids = list(getattr(res, "liquids", None) or [])
        liquid_x = _combined_liquid(res, liquids, self.n) if liquids else None
        split: dict = {}
        if len(liquids) > 1:
            betas = [float(b) for b in res.betas_liquids]  # mol per mol mixture
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
            V=float(res.V()),
            H_reactive=self._h_reactive(float(res.H()), zs),
            **split,
        )


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
                custom: tuple[CustomCompound, ...] = ()) -> PropertyPackage:
    """A cached package; ``custom`` holds the full definitions, so editing a user-defined
    compound under the same name builds a new one."""
    return PropertyPackage(list(components), method, prefer_fitted, henry_components, custom)


def package_options(model) -> dict:
    """The flowsheet's thermodynamic options as keyword arguments for :func:`get_package`."""
    henry = getattr(model, "henry_components", None)
    custom = getattr(model, "custom_compounds", None) or []
    return {"henry_components": tuple(henry) if henry is not None else None,
            "custom": tuple(CustomCompound.from_dict(c.model_dump(exclude_none=True))
                            for c in custom)}
