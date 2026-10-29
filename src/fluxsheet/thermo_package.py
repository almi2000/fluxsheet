"""Property packages: the thermodynamic model behind every flash calculation.

A property package bundles the compound list with a phase-equilibrium method
(equation of state or activity-coefficient model). It is built once per
flowsheet and cached, because loading the constants takes a moment.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from thermo import (
    PRMIX,
    SRKMIX,
    UNIFAC,
    CEOSGas,
    CEOSLiquid,
    ChemicalConstantsPackage,
    FlashPureVLS,
    FlashVL,
    GibbsExcessLiquid,
    IdealGas,
)
from thermo.interaction_parameters import IPDB
from thermo.unifac import DOUFIP2016, DOUFSG

from .compounds import lookup
from .errors import SimulationError


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
    H: float  # J/mol
    S: float  # J/mol/K
    phase: str  # "vapor", "liquid" or "two-phase"
    x: tuple[float, ...] | None  # liquid mole fractions
    y: tuple[float, ...] | None  # vapor mole fractions
    V: float  # molar volume, m3/mol


@dataclass(frozen=True)
class Suggestion:
    method: str
    reason: str


def suggest_method(components: list[str]) -> Suggestion:
    """Pick a sensible default method and say why, the way an experienced engineer would."""
    consts = _constants(tuple(lookup(c).cas for c in components))[0]
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


@lru_cache(maxsize=32)
def _constants(cas_numbers: tuple[str, ...]):
    return ChemicalConstantsPackage.from_IDs(list(cas_numbers))


class PropertyPackage:
    def __init__(self, components: list[str], method: str):
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
        self.compounds = [lookup(c) for c in components]
        self.names = list(components)
        cas = tuple(c.cas for c in self.compounds)
        if len(set(cas)) != len(cas):
            raise SimulationError(
                "the component list contains the same compound twice",
                hint="different names can refer to the same compound (e.g. 'ethanol' "
                "and 'ethyl alcohol'); remove one",
            )
        self.consts, self.props = _constants(cas)
        self._check_data()
        self._flasher = self._build()

    @property
    def n(self) -> int:
        return len(self.names)

    @property
    def mws(self) -> list[float]:
        return list(self.consts.MWs)  # g/mol

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
                IPDB.get_ip_asymmetric_matrix("ChemSep PR", c.CASs, "kij")
                if self.method == "peng-robinson"
                else [[0.0] * self.n for _ in range(self.n)]
            )
            kwargs = {"Tcs": c.Tcs, "Pcs": c.Pcs, "omegas": c.omegas, "kijs": kijs}
            gas = CEOSGas(eos, kwargs, HeatCapacityGases=cp_gas)
            liquid = CEOSLiquid(eos, kwargs, HeatCapacityGases=cp_gas)
        else:
            gas = IdealGas(HeatCapacityGases=cp_gas)
            ge_model = None
            if self.method == "unifac":
                for name, groups in zip(self.names, c.UNIFAC_Dortmund_groups):
                    if not groups:
                        raise SimulationError(
                            f"{name!r} has no UNIFAC group assignment",
                            hint="use Peng-Robinson for this mixture instead",
                        )
                ge_model = UNIFAC.from_subgroups(
                    T=298.15,
                    xs=[1.0 / self.n] * self.n,
                    chemgroups=c.UNIFAC_Dortmund_groups,
                    version=1,
                    interaction_data=DOUFIP2016,
                    subgroups=DOUFSG,
                )
            liquid = GibbsExcessLiquid(
                VaporPressures=p.VaporPressures,
                HeatCapacityGases=cp_gas,
                VolumeLiquids=p.VolumeLiquids,
                GibbsExcessModel=ge_model,
                equilibrium_basis="Psat",
                caloric_basis="Psat",
            )
        if self.n == 1:
            return FlashPureVLS(c, p, gas=gas, liquids=[liquid], solids=[])
        return FlashVL(c, p, liquid=liquid, gas=gas)

    def flash(self, zs: list[float], **spec: float) -> PhaseState:
        """Flash at two of T, P, H, S, VF (SI, molar basis)."""
        try:
            res = self._flasher.flash(zs=list(zs), **spec)
        except Exception as exc:
            spec_text = ", ".join(f"{k}={v:.6g}" for k, v in spec.items())
            raise SimulationError(
                f"phase equilibrium calculation did not converge ({spec_text})",
                hint="check that temperature and pressure are physically reasonable for "
                "these compounds, or try a different property method",
            ) from exc
        vf = float(res.VF) if res.VF is not None else (1.0 if res.phase == "V" else 0.0)
        vf = min(max(vf, 0.0), 1.0)
        gas = getattr(res, "gas", None)
        liquid = res.liquids[0] if getattr(res, "liquids", None) else None
        phase = "two-phase" if 0.0 < vf < 1.0 else ("vapor" if vf >= 1.0 else "liquid")
        return PhaseState(
            T=float(res.T),
            P=float(res.P),
            VF=vf,
            H=float(res.H()),
            S=float(res.S()),
            phase=phase,
            x=tuple(liquid.zs) if liquid is not None and vf < 1.0 else None,
            y=tuple(gas.zs) if gas is not None and vf > 0.0 else None,
            V=float(res.V()),
        )


@lru_cache(maxsize=16)
def get_package(components: tuple[str, ...], method: str) -> PropertyPackage:
    return PropertyPackage(list(components), method)
