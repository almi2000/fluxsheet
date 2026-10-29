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
            "nrtl",
            "NRTL",
            "Activity-coefficient model with binary parameters regressed to measured "
            "phase-equilibrium data (ChemSep databank). The first choice for strongly non-ideal "
            "liquids, azeotropes and liquid-liquid splits (decanters), at low to moderate "
            "pressure. Pairs without data are estimated from UNIFAC and flagged. Ideal vapor.",
        ),
        Method(
            "uniquac",
            "UNIQUAC",
            "Activity-coefficient model with ChemSep binary parameters; molecular size and "
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
    splits = _liquid_liquid_pairs(consts, polar)
    if splits:
        pairs = " and ".join(f"{a}/{b}" for a, b, _, _ in splits)
        verb = "forms" if len(splits) == 1 else "form"
        lle = (f"UNIFAC predicts that {pairs} {verb} two liquid phases at 25 °C, so expect a "
               "liquid-liquid split (a Decanter can separate it). ")
        unmeasured = [f"{a}/{b}" for a, b, ca, cb in splits
                      if not IPDB.has_ip_specific(_GE_TABLES["nrtl"], [ca, cb], "bij")]
        if not unmeasured:
            return Suggestion(
                "nrtl",
                lle + "NRTL with ChemSep binary parameters, regressed to measured data, "
                "describes both liquid phases and the vapor in equilibrium with them; Wilson "
                "cannot represent a split at all.",
            )
        return Suggestion(
            "unifac",
            lle + f"ChemSep has no measured NRTL parameters for {', '.join(unmeasured)}, so "
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


@lru_cache(maxsize=32)
def _constants(cas_numbers: tuple[str, ...]):
    return ChemicalConstantsPackage.from_IDs(list(cas_numbers))


# ---- activity-coefficient models ---------------------------------------------------------

# Binary parameter tables from the ChemSep databank, as shipped with `thermo`.
_GE_TABLES = {"nrtl": "ChemSep NRTL", "uniquac": "ChemSep UNIQUAC", "wilson": "ChemSep Wilson"}
# Methods whose liquid model can describe two coexisting liquid phases.
LLE_METHODS = ("nrtl", "uniquac", "unifac")


def _unifac(consts, T: float = 298.15):
    """UNIFAC (Dortmund) model for these compounds, or None if a compound lacks groups."""
    if any(not groups for groups in consts.UNIFAC_Dortmund_groups):
        return None
    return UNIFAC.from_subgroups(
        T=T, xs=[1.0 / consts.N] * consts.N, chemgroups=consts.UNIFAC_Dortmund_groups,
        version=1, interaction_data=DOUFIP2016, subgroups=DOUFSG,
    )


def _ge_model(method: str, consts, a, b, alpha):
    """A regressed activity model; ``a``, ``b`` (and ``alpha`` for NRTL) are n x n matrices."""
    xs = [1.0 / consts.N] * consts.N
    if method == "nrtl":  # tau_ij = a_ij + b_ij / T
        return NRTL(T=298.15, xs=xs, tau_as=a, tau_bs=b, alpha_cs=alpha)
    if method == "uniquac":  # tau_ij = exp(a_ij + b_ij / T)
        return UNIQUAC(T=298.15, xs=xs, rs=consts.UNIFAC_Rs, qs=consts.UNIFAC_Qs,
                       tau_as=a, tau_bs=b)
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
        # Where each binary parameter set came from, and warnings worth showing with results.
        self.parameter_sources: list[dict] = []
        self.warnings: list[str] = []
        self._ge = None  # activity-coefficient model, if the method has one
        self._check_data()
        self._flasher = self._build()

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
                ge_model = _unifac(c)
            elif self.method in _GE_TABLES and self.n > 1:
                ge_model = self._regressed_model()
            self._ge = ge_model
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
        if self.supports_lle:
            return _FlashVLLE(c, p, liquids=[liquid, liquid], gas=gas)
        return FlashVL(c, p, liquid=liquid, gas=gas)

    def _regressed_model(self):
        """NRTL, UNIQUAC or Wilson with ChemSep parameters; UNIFAC estimates fill the gaps."""
        c, n, table = self.consts, self.n, _GE_TABLES[self.method]
        if self.method == "uniquac":
            for name, r, q in zip(self.names, c.UNIFAC_Rs, c.UNIFAC_Qs):
                if r is None or q is None:
                    raise SimulationError(
                        f"{name!r} has no UNIQUAC size and surface parameters (r, q)",
                        hint="they come from UNIFAC groups, which this compound lacks; "
                        "use NRTL or Peng-Robinson instead",
                    )
        b = IPDB.get_ip_asymmetric_matrix(table, c.CASs, "bij")
        a = (IPDB.get_ip_asymmetric_matrix(table, c.CASs, "aij") if self.method == "wilson"
             else [[0.0] * n for _ in range(n)])
        alpha = (IPDB.get_ip_asymmetric_matrix(table, c.CASs, "alphaij")
                 if self.method == "nrtl" else None)
        label = METHODS[self.method].label
        for i in range(n):
            for j in range(i + 1, n):
                pair = f"{self.names[i]}/{self.names[j]}"
                if IPDB.has_ip_specific(table, [c.CASs[i], c.CASs[j]], "bij"):
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
                    f"{label}: ChemSep has no measured binary parameters for {pair}; they were "
                    "estimated by fitting to UNIFAC (Dortmund) predictions. Treat results for "
                    "this pair as predictive, not data-based."
                )
        return _ge_model(self.method, c, a, b, alpha)

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
        f = self._flasher
        if kind == "vapor":
            template = f.gas
        elif kind == "liquid":
            template = f.liquid if hasattr(f, "liquid") else f.liquids[0]
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
                return self._checked_vf(zs, spec, guess)
            res = self._flasher.flash(zs=zs, **spec)
            return res, _incipient(res, spec.get("VF"), self.n)
        except Exception as exc:
            spec_text = ", ".join(f"{k}={v:.6g}" for k, v in spec.items())
            raise SimulationError(
                f"phase equilibrium calculation did not converge ({spec_text})",
                hint="check that temperature and pressure are physically reasonable for "
                "these compounds, or try a different property method",
            ) from exc

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
                      -math.log(sum(z * vp(spec["T"])
                                    for z, vp in zip(zs, self.props.VaporPressures))))

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
            return [z * g * vp(T) for z, g, vp in zip(zs, gammas, vps)]

        try:
            if "T" in spec:
                T = spec["T"]
            else:
                P = spec["P"]
                tsat = [vp.solve_property(P) for vp in vps]
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
def get_package(components: tuple[str, ...], method: str) -> PropertyPackage:
    return PropertyPackage(list(components), method)
