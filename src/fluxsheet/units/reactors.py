"""Reactors.

All reactors close the energy balance on enthalpies that include formation enthalpies
(``H + sum z_i ΔHf_i``), so the heat of reaction comes out of the balance by itself. Each
runs isothermally at a given outlet temperature, adiabatically (the default), or with a
given duty, with an optional pressure drop. Reactions are defined at flowsheet level (see
reactions.py) and referenced by id; the solver hands every reactor the parsed set.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space
from scipy.optimize import brentq, minimize_scalar, root

from ..errors import SimulationError
from ..reactions import (
    P_STD,
    R,
    Reaction,
    ReactionSet,
    element_matrix,
    equilibrium_extents,
    gibbs_over_rt,
)
from ..stream import Stream
from ..utilities import WallMedium, stream_medium, utility_mass_flow, utility_medium
from .base import (
    DP,
    DUTY,
    IN,
    OUT,
    P_OUT,
    Param,
    Port,
    ResultField,
    UnitOperation,
    empty_stream,
    single_inlet,
)
from .heat_transfer import AREA, HEAT_PARAMS, HEAT_RESULTS, UTILITY, UTILITY_FLOW

T_SPEC = Param("T", "Outlet temperature", "temperature", min=0.0,
               help="Isothermal operation at this temperature.")
DUTY_SPEC = Param("duty", "Duty", "power",
                  help="Heat added (positive) or removed (negative). Leave the temperature and "
                  "duty unset for an adiabatic reactor.")
REACTIONS = Param("reactions", "Reactions", "string_list",
                  help="Ids of flowsheet reactions that run in this reactor.")
KEY = Param("key", "Key component", "component",
            help="Component whose conversion is reported. Defaults to the first reactant of "
            "the first reaction.")
PHASE = Param("phase", "Reacting phase", "choice", default="gas", options=("gas", "liquid"),
              help="Phase the reactions take place in; sets concentrations and activities.")

# Isothermal, fixed duty, adiabatic (none of them), or a jacket/coil with a utility.
ENERGY = (T_SPEC, DUTY_SPEC, AREA, *HEAT_PARAMS)

T_OUT = ResultField("T_out", "Outlet temperature", "temperature")
CONVERSION = ResultField("conversion", "Conversion of key component", "fraction")
RESIDENCE_TIME = ResultField("residence_time", "Residence time (inlet conditions)", "time")
ENERGY_GROUPS = (("T", "duty", "area"), ("P", "dP"))


class Reactor(UnitOperation):
    """Shared energy balance and bookkeeping. Not a unit type itself (no ``type``)."""

    category = "Reactors"
    inlets = (IN,)
    outlets = (OUT,)
    exclusive_groups = ENERGY_GROUPS
    result_fields = (T_OUT, DUTY, CONVERSION, *HEAT_RESULTS)
    reactions = ReactionSet()  # replaced by the solver with the flowsheet's reactions

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P = self.outlet_pressure(s.P)
        if P <= 0.0:
            raise SimulationError(f"the pressure drop leaves no pressure (outlet {P:g} Pa)",
                                  hint="reduce the pressure drop")
        if s.total_flow <= 0.0:
            self.results.update(duty=0.0, T_out=self.p["T"] or s.T)
            return {"out": empty_stream(pkg, self.p["T"] or s.T, P)}
        return {"out": self.react(pkg, s, P)}

    def react(self, pkg, s: Stream, P: float) -> Stream:
        raise NotImplementedError

    # ---- energy balance ----------------------------------------------------------------

    def formation_enthalpies(self, pkg, flows=(), species=()) -> np.ndarray:
        """ΔHf° per component; components that enter, leave or react must have one."""
        needed = {i for f in flows for i, v in enumerate(f) if v > 0.0} | set(species)
        missing = [pkg.names[i] for i in sorted(needed) if pkg.consts.Hfgs[i] is None]
        if missing:
            raise SimulationError(
                f"no formation enthalpy in the database for {', '.join(missing)}; the "
                "reactor's energy balance needs it for every compound that enters, leaves or "
                "reacts",
                where=self.id,
                hint="choose compounds that have formation data, or keep these out of the "
                "reactor",
            )
        return np.array([h if h is not None else 0.0 for h in pkg.consts.Hfgs])

    @staticmethod
    def h_flow(s: Stream, hf: np.ndarray) -> float:
        """Enthalpy flow including formation enthalpies, W."""
        return s.total_flow * (s.state.H + float(np.dot(s.zs, hf)))

    def clean(self, flows: np.ndarray, scale: float) -> list[float]:
        if (flows < -1e-8 * scale).any():
            bad = [name for name, f in zip(self._names, flows) if f < -1e-8 * scale]
            raise SimulationError(
                f"the reactions consume more {', '.join(bad)} than is available",
                where=self.id, hint="lower the conversions, or check the feed",
            )
        return [max(float(f), 0.0) for f in flows]

    def finish(self, pkg, s: Stream, flows, P: float, hf: np.ndarray,
               duty: float | None = None) -> Stream:
        """Outlet stream from the energy specification, or from ``duty`` when the unit found
        it itself; records duty, T_out, conversion."""
        flows = self.clean(np.asarray(flows, dtype=float), s.total_flow)
        if self.p["T"] is not None and duty is None:
            out = Stream.from_flash(pkg, flows, T=self.p["T"], P=P)
        else:
            if duty is None:
                duty = self.p["duty"] or 0.0
            target = self.h_flow(s, hf) + duty
            total = sum(flows)
            h = target / total - sum(f / total * v for f, v in zip(flows, hf))
            out = Stream.from_flash(pkg, flows, H=h, P=P)
        self.results["duty"] = self.h_flow(out, hf) - self.h_flow(s, hf)
        self.results["T_out"] = out.T
        key = self.key_index(pkg)
        if key is not None and s.flows[key] > 0.0:
            self.results["conversion"] = 1.0 - out.flows[key] / s.flows[key]
        return out

    def key_index(self, pkg) -> int | None:
        if self.p.get("key"):
            return self.component_index(pkg, self.p["key"], "key component")
        if self._key is not None:
            return self._key
        for r in self._rxns:
            return r.first_reactant
        return None

    def temperature_for_energy(self, pkg, s: Stream, P: float, hf: np.ndarray, flows_at):
        """Outlet T of an adiabatic / fixed-duty reactor whose outlet flows depend on T."""
        target = self.h_flow(s, hf) + (self.p["duty"] or 0.0)

        def residual(T: float) -> float:
            flows = np.maximum(flows_at(T), 0.0)
            total = flows.sum()
            st = pkg.flash(list(flows / total), T=T, P=P)
            return total * (st.H + float(np.dot(flows / total, hf))) - target

        lo = hi = s.T
        f_lo = f_hi = residual(s.T)
        step = 50.0
        while f_lo > 0.0 or f_hi < 0.0:
            if f_hi < 0.0:
                lo, f_lo = hi, f_hi
                hi += step
                if hi > 4000.0:
                    raise SimulationError("the adiabatic temperature is above 4000 K",
                                          where=self.id, hint="remove more heat (duty)")
                f_hi = residual(hi)
            else:
                hi, f_hi = lo, f_lo
                lo -= step
                if lo < 50.0:
                    raise SimulationError("the outlet temperature would fall below 50 K",
                                          where=self.id, hint="add heat (duty)")
                f_lo = residual(lo)
            step *= 1.5
        return brentq(residual, lo, hi, xtol=1e-7, rtol=1e-12) if f_lo < 0 < f_hi else (
            lo if f_lo == 0.0 else hi)

    def select_reactions(self, pkg) -> list[Reaction]:
        self._rxns = self.reactions.select(self.p["reactions"], self.id)
        self._names = pkg.names
        return self._rxns

    def __post_init__(self) -> None:
        super().__post_init__()
        self._rxns: list[Reaction] = []
        self._names: list[str] = []
        self._key: int | None = None  # key component when not given as a parameter


# ---- conversion / yield --------------------------------------------------------------------

def _fraction(text: str, where: str) -> float:
    t = text.strip()
    try:
        x = float(t[:-1]) / 100.0 if t.endswith("%") else float(t)
    except ValueError:
        raise SimulationError(f"conversion {text!r} is not a number", where=where,
                              hint="write fractions like 0.9 or 90%") from None
    if not 0.0 <= x <= 1.0:
        raise SimulationError(f"conversion {text!r} must be between 0 and 1", where=where)
    return x


class StoichiometricReactor(Reactor):
    type = "StoichiometricReactor"
    label = "Stoichiometric Reactor"
    description = ("Reactions run to a given conversion of a key component each, one after "
                   "another or in parallel. No kinetics or equilibrium needed.")
    params = (
        REACTIONS,
        Param("conversions", "Conversions", "string_list",
              help="Fractional conversion of each reaction's key component, in the order of the "
              "reactions, e.g. 0.9, 50%."),
        Param("keys", "Key components", "string_list",
              help="Key component of each reaction, in order. Defaults to each reaction's first "
              "reactant."),
        Param("mode", "Reactions run", "choice", default="sequential",
              options=("sequential", "parallel"),
              help="Sequential: each conversion applies to what the previous reactions left. "
              "Parallel: all conversions refer to the inlet."),
        *ENERGY, P_OUT, DP,
    )

    def react(self, pkg, s, P):
        rxns = self.select_reactions(pkg)
        conv = [_fraction(x, self.id) for x in (self.p["conversions"] or [])]
        if len(conv) != len(rxns):
            raise SimulationError(
                f"{len(rxns)} reactions but {len(conv)} conversions", where=self.id,
                hint="give one conversion per reaction, in the same order",
            )
        keys = self.p["keys"] or []
        if keys and len(keys) != len(rxns):
            raise SimulationError(f"{len(rxns)} reactions but {len(keys)} key components",
                                  where=self.id, hint="give one key per reaction, or none")
        inlet = np.array(s.flows)
        flows = inlet.copy()
        extents = {}
        for j, (r, x) in enumerate(zip(rxns, conv)):
            k = (self.component_index(pkg, keys[j], "key component") if keys
                 else r.first_reactant)
            if r.nu[k] >= 0:
                raise SimulationError(
                    f"key component {pkg.names[k]!r} is not a reactant of {r.id}", where=self.id,
                    hint=f"pick a reactant of {r.equation}",
                )
            basis = inlet if self.p["mode"] == "parallel" else flows
            xi = x * basis[k] / -r.nu[k]
            flows = flows + xi * r.nu
            extents[r.id] = xi
            if j == 0:
                self._key = k
        hf = self.formation_enthalpies(pkg, [inlet, np.maximum(flows, 0.0)])
        self.results["extents"] = extents
        return self.finish(pkg, s, flows, P, hf)


class YieldReactor(Reactor):
    type = "YieldReactor"
    label = "Yield Reactor"
    description = ("Turns the feed into products with given mass yields, scaled so mass is "
                   "conserved. For reactions whose stoichiometry is unknown, e.g. pyrolysis.")
    result_fields = (T_OUT, DUTY, *HEAT_RESULTS)
    params = (
        Param("yields", "Mass yields", "component_fractions",
              help="Relative mass of each product; normalized so the outlet mass equals the "
              "reacting inlet mass."),
        Param("inerts", "Inert components", "string_list",
              help="Components that pass through unchanged."),
        *ENERGY, P_OUT, DP,
    )

    def react(self, pkg, s, P):
        self._names = pkg.names
        yields = self.p["yields"] or {}
        inerts = {self.component_index(pkg, n, "inert component") for n in self.p["inerts"] or []}
        weights = np.zeros(pkg.n)
        for name, w in yields.items():
            i = self.component_index(pkg, name, "yield component")
            if i in inerts:
                raise SimulationError(f"{name!r} is both inert and a product", where=self.id,
                                      hint="remove it from one of the lists")
            if w < 0:
                raise SimulationError(f"the yield of {name!r} is negative", where=self.id)
            weights[i] = w
        if weights.sum() <= 0.0:
            raise SimulationError("no product yields given", where=self.id,
                                  hint="give a positive mass yield for at least one component")
        mw = np.array(pkg.mws)
        inlet = np.array(s.flows)
        reacting = np.array([0.0 if i in inerts else f for i, f in enumerate(inlet)])
        mass = float(reacting @ mw)
        flows = inlet - reacting + mass * weights / weights.sum() / mw
        hf = self.formation_enthalpies(pkg, [inlet, flows])
        return self.finish(pkg, s, flows, P, hf)


# ---- equilibrium ---------------------------------------------------------------------------

class EquilibriumReactor(Reactor):
    type = "EquilibriumReactor"
    label = "Equilibrium Reactor"
    description = ("Chosen reactions reach chemical equilibrium. K(T) from Gibbs energies of "
                   "formation (van 't Hoff with ideal-gas heat capacities) or from the "
                   "reaction's own correlation.")
    params = (REACTIONS, PHASE, KEY, *ENERGY, P_OUT, DP)

    def react(self, pkg, s, P):
        rxns = self.select_reactions(pkg)
        gas = self.p["phase"] == "gas"
        if not gas:
            for r in rxns:
                if r.keq is None:
                    raise SimulationError(
                        f"reaction {r.id} needs a user equilibrium constant (keq) in the liquid "
                        "phase", where=self.id,
                        hint="the built-in K refers to ideal gases; add ln K = A + B/T + C ln T "
                        "+ D T to the reaction",
                    )
        nu = np.array([r.nu for r in rxns])
        n0 = np.array(s.flows)
        ln_p = math.log(P / P_STD) if gas else 0.0
        species = sorted({i for r in rxns for i in r.involved})
        hf = self.formation_enthalpies(pkg, [n0], species)

        def flows_at(T: float) -> np.ndarray:
            ln_k = np.array([r.ln_keq(pkg, T) for r in rxns])
            self._xi = equilibrium_extents(n0, nu, ln_k, ln_p)
            return n0 + self._xi @ nu

        T = self.p["T"] if self.p["T"] is not None else (
            self.temperature_for_energy(pkg, s, P, hf, flows_at))
        flows = flows_at(T)
        self.results["extents"] = {r.id: float(x) for r, x in zip(rxns, self._xi)}
        return self.finish(pkg, s, flows, P, hf)


class GibbsReactor(Reactor):
    type = "GibbsReactor"
    label = "Gibbs Reactor"
    description = ("Minimizes the Gibbs energy subject to element balances: of an ideal-gas "
                   "mixture, or on fugacities from the property method, with the vapor and "
                   "liquid phases the mixture forms. Needs no reactions, only which compounds "
                   "may form.")
    params = (
        Param("products", "Reacting components", "string_list",
              help="Components that may form or react. Others pass through unchanged. "
              "Default: all components."),
        Param("fugacities", "Fugacities", "choice", default="ideal_gas",
              options=("ideal_gas", "property_method"),
              help="ideal_gas: an ideal-gas mixture (exact at low pressure and high "
              "temperature). property_method: fugacities from the flowsheet's property "
              "method, so high pressures (an equation of state) and condensing products "
              "(vapor and liquid in phase equilibrium) are accounted for."),
        KEY, *ENERGY, P_OUT, DP,
    )

    def react(self, pkg, s, P):
        self._names = pkg.names
        names = self.p["products"] or pkg.names
        species = sorted({self.component_index(pkg, n, "reacting component") for n in names})
        _, A = element_matrix(pkg, species)
        for i in species:
            if not pkg.consts.atomss[i]:
                raise SimulationError(f"the database has no formula for {pkg.names[i]!r}",
                                      where=self.id, hint="leave it out of the reacting "
                                      "components")
        basis = null_space(A).T  # independent reactions among the species
        nu = np.zeros((basis.shape[0], pkg.n))
        nu[:, species] = basis
        n0 = np.array(s.flows)
        ln_p = math.log(P / P_STD)
        hf = self.formation_enthalpies(pkg, [n0], species)

        real = self.p["fugacities"] == "property_method"

        def flows_at(T: float) -> np.ndarray:
            if not len(nu):
                return n0
            g = gibbs_over_rt(pkg, T, species, why="the Gibbs reactor")
            ln_k = -nu[:, species] @ g[species]
            n = n0 + equilibrium_extents(n0, nu, ln_k, ln_p) @ nu
            return self._fugacity_equilibrium(pkg, T, P, n0, nu, ln_k, ln_p, n) if real else n

        T = self.p["T"] if self.p["T"] is not None else (
            self.temperature_for_energy(pkg, s, P, hf, flows_at))
        return self.finish(pkg, s, flows_at(T), P, hf)

    def _fugacity_equilibrium(self, pkg, T, P, n0, nu, ln_k, ln_p, n):
        """Equilibrium on fugacities f_i = z_i Phi_i P of the whole mixture (Phi from its
        phase split, see PropertyPackage.ln_mixture_fugacity_coefficients).

        sum_i nu_i ln(f_i/P°) = ln K reads sum_i nu_i ln(z_i P/P°) = ln K - sum_i nu_i ln Phi_i:
        the ideal-gas problem with a corrected K. Phi follows the composition, so the two are
        iterated (successive substitution, damped when it oscillates); the fixed point is the
        Gibbs minimum with phase and chemical equilibrium together.
        """
        lnphi = None
        weight, last = 1.0, math.inf
        total0 = n0.sum()
        for _ in range(200):
            F = np.maximum(n, 0.0)
            new, _ = pkg.ln_mixture_fugacity_coefficients(list(F / F.sum()), T, P)
            change = math.inf if lnphi is None else float(np.max(np.abs(new - lnphi)))
            if change > last:  # oscillating: take smaller steps
                weight = max(0.5 * weight, 0.05)
            last = change
            lnphi = new if lnphi is None else lnphi + weight * (new - lnphi)
            n_new = n0 + equilibrium_extents(n0, nu, ln_k - nu @ lnphi, ln_p) @ nu
            if change < 1e-10 and np.max(np.abs(n_new - n)) <= 1e-12 * total0:
                return n_new
            n = n_new
        raise SimulationError(
            "the chemical and phase equilibrium did not converge", where=self.id,
            hint="try the ideal-gas fugacities, or a temperature away from where a phase "
            "appears")

    def key_index(self, pkg):
        return super().key_index(pkg) if self.p.get("key") else None


# ---- kinetic reactors ----------------------------------------------------------------------

class KineticReactor(Reactor):
    """CSTR and PFR: rates from the reactions' kinetics."""

    result_fields = (T_OUT, DUTY, CONVERSION, RESIDENCE_TIME, *HEAT_RESULTS)

    def setup(self, pkg) -> list[Reaction]:
        rxns = self.select_reactions(pkg)
        self._phase = self.p["phase"]
        self._kind = "vapor" if self._phase == "gas" else "liquid"
        for r in rxns:
            kin = r.kinetics
            if kin is None:
                raise SimulationError(
                    f"reaction {r.id} has no kinetics, which a {self.label} needs",
                    where=self.id,
                    hint=f"add a rate law to {r.id} in the flowsheet's reactions, or use an "
                    "Equilibrium or Stoichiometric Reactor",
                )
            if kin.basis == "partial_pressure" and self._phase == "liquid":
                raise SimulationError(
                    f"reaction {r.id} has a partial-pressure rate law but the reacting phase "
                    "is liquid", where=self.id,
                    hint="use the concentration basis, or set the reacting phase to gas",
                )
            if kin.equilibrium_limited and self._phase == "liquid" and r.keq is None:
                raise SimulationError(
                    f"reaction {r.id} is equilibrium-limited in the liquid phase but has no "
                    "user equilibrium constant", where=self.id,
                    hint="add keq (ln K = A + B/T + C ln T + D T, K on mole fractions) to the "
                    "reaction; the built-in K refers to ideal gases",
                )
        self._nu = np.array([r.nu for r in rxns])
        return rxns

    def rates(self, pkg, T: float, P: float, flows: np.ndarray) -> np.ndarray:
        F = np.maximum(flows, 0.0)
        z = F / F.sum()
        Vm = pkg.phase_props(self._kind, T, P, list(z)).V
        return np.array([r.rate(pkg, T, P, z, Vm, self._phase) for r in self._rxns])

    def h_phase(self, pkg, T: float, P: float, flows: np.ndarray, hf: np.ndarray) -> float:
        """Enthalpy flow (with formation enthalpies) of the reacting phase, W."""
        F = np.maximum(flows, 0.0)
        total = F.sum()
        z = F / total
        return total * (pkg.phase_props(self._kind, T, P, list(z)).H + float(z @ hf))

    def residence_time(self, pkg, s: Stream, V: float) -> float:
        q = s.total_flow * pkg.phase_props(self._kind, s.T, s.P, s.zs).V
        return V / q


class CSTR(KineticReactor):
    type = "CSTR"
    label = "CSTR"
    description = ("Continuous stirred tank: perfectly mixed, so the whole volume reacts at "
                   "outlet conditions. Rates from the reactions' kinetics.")
    params = (
        REACTIONS,
        Param("V", "Volume", "volume", min=0.0, help="Volume of the reacting phase."),
        PHASE, KEY, *ENERGY, P_OUT, DP,
    )
    exclusive_groups = (("V",), *ENERGY_GROUPS)
    required_groups = (("V",),)

    def react(self, pkg, s, P):
        rxns = self.setup(pkg)
        nu, V = self._nu, self.p["V"]
        n0 = np.array(s.flows)
        scale = s.total_flow
        species = sorted({i for r in rxns for i in r.involved})
        hf = self.formation_enthalpies(pkg, [n0], species)
        fixed_T = self.p["T"]
        target = self.h_flow(s, hf) + (self.p["duty"] or 0.0)
        e_scale = scale * R * s.T

        def residual(x, vol):
            xi = x[: len(rxns)] * scale
            T = fixed_T if fixed_T is not None else x[-1] * s.T
            flows = n0 + xi @ nu
            res = (xi - vol * self.rates(pkg, T, P, flows)) / scale
            if fixed_T is None:
                res = np.r_[res, (self.h_phase(pkg, T, P, flows, hf) - target) / e_scale]
            return res

        # Continuation in volume from an empty tank: robust, and stays on the branch that
        # a start-up from cold would reach when several steady states exist.
        x = np.zeros(len(rxns)) if fixed_T is not None else np.r_[np.zeros(len(rxns)), 1.0]
        done, frac = 0.0, 1e-4
        for _ in range(200):
            trial = min(1.0, done + frac)
            sol = root(residual, x, args=(trial * V,), method="hybr", options={"xtol": 1e-12})
            ok = sol.success and np.max(np.abs(residual(sol.x, trial * V))) < 1e-8
            if ok:
                flows = n0 + sol.x[: len(rxns)] * scale @ nu
                ok = bool((flows >= -1e-9 * scale).all())
            if ok:
                x, done = sol.x, trial
                if done >= 1.0:
                    break
                frac *= 3.0
            else:
                frac /= 4.0
                if frac < 1e-9:
                    break
        if done < 1.0:
            raise SimulationError(
                "the CSTR balances did not converge", where=self.id,
                hint="check the rate constants and units, or give a temperature",
            )
        xi = x[: len(rxns)] * scale
        self.results["extents"] = {r.id: float(v) for r, v in zip(rxns, xi)}
        self.results["residence_time"] = self.residence_time(pkg, s, V)
        return self.finish(pkg, s, n0 + xi @ nu, P, hf)




# ---- plug flow, cooled or heated along the tube --------------------------------------------

COOLANT_IN = Port("coolant_in", "Medium inlet")
COOLANT_OUT = Port("coolant_out", "Medium outlet")
MEDIUM = Param(
    "medium", "Medium side", "choice", default="utility", options=("utility", "stream"),
    help="What flows outside the tube: a utility from the flowsheet's list, or a process "
    "stream connected to the medium inlet and outlet (e.g. the cold feed that cools the "
    "reactor before it enters it).")
WALL_U = Param(
    "U", "Heat transfer coefficient U", "heat_transfer_coefficient", min=0.0,
    help="Overall coefficient through the tube wall, on the heat transfer area. With a "
    "medium and neither temperature nor duty, the heat exchanged along the tube follows "
    "from U, the area and the local temperatures of both sides.")
AREA_PER_LENGTH = Param(
    "area_per_length", "Heat transfer area per length", "length", min=0.0, requires="U",
    help="Wall area per metre of tube in m²/m, entered as a length. Defaults to the tube "
    "circumference π·D.")
WALL_ARRANGEMENT = Param(
    "arrangement", "Flow arrangement", "choice", default="countercurrent",
    options=("countercurrent", "cocurrent"),
    help="Direction of the medium relative to the reacting stream. Countercurrent: the "
    "medium enters at the reactor outlet.")
WALL_DUTY = ResultField("wall_duty", "Heat through the wall", "power")
T_HOT_SPOT = ResultField("T_hot_spot", "Hot-spot temperature", "temperature")
HOT_SPOT_AT = ResultField("hot_spot_position", "Hot-spot position", "length")


class PFR(KineticReactor):
    """Plug flow reactor: adiabatic, isothermal, with a fixed duty, or exchanging heat
    along its length with a medium outside the tube.

    Along the volume V the balances are dξ/dV = r and, for the energy, dH/dV = q with
    q = U a_V (T_m - T): U the overall coefficient, a_V the wall area per reactor volume
    and T_m the medium temperature there. The medium's balance is integrated with them as
    Q(V), the heat passed into the process up to V. Cocurrent, the medium has given up Q(V)
    at V. Countercurrent, it enters at the reactor outlet, so at V it has given up
    Q_tot - Q(V), where Q_tot = Q(V_end) is the unknown total duty: a boundary-value
    problem, solved by shooting on Q_tot (equivalently on the medium outlet temperature)
    with a bracketing root finder. A fluid utility without a given flow takes the flow that
    brings it to its outlet temperature; that flow depends on Q_tot, so it is shot on Q_tot
    as well, in either arrangement.
    """

    type = "PFR"
    label = "Plug Flow Reactor"
    description = ("Tubular reactor without back-mixing: the balances are integrated along the "
                   "volume. Gives conversion, temperature and composition profiles. Can be "
                   "cooled or heated along its length by a utility or a process stream, in "
                   "co- or countercurrent.")
    inlets = (IN, COOLANT_IN)
    outlets = (OUT, COOLANT_OUT)
    params = (
        REACTIONS,
        Param("V", "Volume", "volume", min=0.0, help="Reactor volume."),
        Param("length", "Length", "length", min=0.0, help="Tube length; needs the diameter."),
        Param("diameter", "Diameter", "length", min=0.0,
              help="Tube inner diameter, used with the length."),
        PHASE, KEY, T_SPEC, DUTY_SPEC, MEDIUM, UTILITY, UTILITY_FLOW, WALL_U, AREA_PER_LENGTH,
        WALL_ARRANGEMENT, P_OUT, DP,
    )
    exclusive_groups = (("V", "length"), ("T", "duty"), ("P", "dP"))
    required_groups = (("V", "length"),)
    result_fields = (T_OUT, DUTY, CONVERSION, RESIDENCE_TIME, WALL_DUTY, T_HOT_SPOT,
                     HOT_SPOT_AT, *HEAT_RESULTS)
    process_side = "plug"
    N_POINTS = 41
    RTOL = 1e-8  # relative tolerance of the integration along the tube

    def __post_init__(self) -> None:
        self._utility = None
        self._rated = False
        self._coolant: Stream | None = None
        self._wall_heat: float | None = None
        self._feed_T = 0.0
        super().__post_init__()
        p = self.p
        stream = p["medium"] == "stream"
        if stream and p["utility"] is not None:
            raise SimulationError("the medium is a process stream, but a utility is chosen too",
                                  where=self.id, hint="clear the utility, or set the medium "
                                  "side to utility")
        if stream and p["U"] is None and p["T"] is None and p["duty"] is None:
            raise SimulationError(
                "the medium stream exchanges heat through the tube wall, which needs U",
                where=self.id, hint="give U (with the length and diameter), or a temperature "
                "or duty for the medium stream to take up")
        # The wall is rated along the tube: the duty follows from U, the area and the local
        # temperatures. Otherwise a utility reports its consumption for the given T or duty.
        self._rated = (p["U"] is not None and p["T"] is None and p["duty"] is None
                       and (stream or p["utility"] is not None))
        if not stream:  # no medium stream: the reactor has one inlet and one outlet
            self.inlets, self.outlets = (IN,), (OUT,)

    # The solver's utilities.calculate() finds the duty of a well-mixed jacket or reports a
    # utility's consumption. A wall rated along the tube is integrated here instead, so the
    # wrapper sees no utility then.
    @property
    def utility(self):
        return None if self._rated else self._utility

    @utility.setter
    def utility(self, medium) -> None:
        self._utility = medium

    def calculate(self, pkg, inlets):
        c = inlets["coolant_in"][0] if "coolant_in" in inlets else None
        self._coolant, self._wall_heat = c, None
        outputs = super().calculate(pkg, inlets)
        if c is None:
            return outputs
        # Heat the medium stream passed into the process: rated along the wall, or the duty
        # that keeps the given temperature or duty.
        Q = self._wall_heat if self._wall_heat is not None else self.results.get("duty") or 0.0
        if c.total_flow > 0.0:
            out = Stream.from_flash(pkg, c.flows, P=c.P, H=c.state.H - Q / c.total_flow)
        elif abs(Q) > 1e-9:
            raise SimulationError("the medium stream has no flow, so it cannot take up the "
                                  "reactor's duty", where=self.id,
                                  hint="give the medium stream a flow, or clear the "
                                  "temperature or duty")
        else:
            out = empty_stream(pkg, c.T, c.P)
        # The heat stays inside the unit: no duty crosses its boundary.
        self.results.update(duty=0.0, wall_duty=Q, utility_T_out=out.T,
                            utility_flow=sum(f * mw for f, mw in zip(c.flows, pkg.mws)) / 1e3)
        outputs["coolant_out"] = out
        return outputs

    def _geometry(self) -> tuple[float, float | None]:
        """Volume and, if the tube is given by length and diameter, its length."""
        L = self.p["length"]
        if L is not None:
            if not self.p["diameter"]:
                raise SimulationError("the length is given but not the diameter", where=self.id,
                                      hint="set the diameter, or give the volume instead")
            V = math.pi / 4.0 * self.p["diameter"] ** 2 * L
        else:
            V = self.p["V"]
        if not V or V <= 0.0:
            raise SimulationError("the reactor volume must be positive", where=self.id)
        return V, L

    def react(self, pkg, s, P_out):
        rxns = self.setup(pkg)
        nu, m = self._nu, len(rxns)
        V, L = self._geometry()
        n0 = np.array(s.flows)
        scale = s.total_flow
        species = sorted({i for r in rxns for i in r.involved})
        hf = self.formation_enthalpies(pkg, [n0], species)
        fixed_T = self.p["T"]
        q_v = (self.p["duty"] or 0.0) / V
        P_in = s.P
        dPdV = (P_out - P_in) / V
        e_scale = scale * R * s.T  # W, scales the heat through the wall
        counter = self.p["arrangement"] == "countercurrent"
        wall = self._rated and not (self._coolant is not None and self._coolant.total_flow <= 0)
        if wall and L is None:
            raise SimulationError(
                "the heat exchanged along the tube needs its length", where=self.id,
                hint="give the length and diameter instead of the volume")
        area = (self.p["area_per_length"] or math.pi * self.p["diameter"]) * L if wall else 0.0
        UA_v = self.p["U"] * area / V if wall else 0.0  # W/(m³ K)

        def rhs(vol, y, medium, Q_tot):
            xi = y[:m] * scale
            T = fixed_T if fixed_T is not None else y[m]
            if not 1.0 < T < 6000.0:
                raise SimulationError(
                    f"the temperature along the reactor runs away ({T:.4g} K)", where=self.id,
                    hint="check the rate constants and units, or remove more heat")
            P = P_in + dPdV * vol
            flows = n0 + xi @ nu
            r = self.rates(pkg, T, P, flows)
            if fixed_T is not None:
                return r / scale
            if medium is None:
                q = q_v
            else:
                Q = y[m + 1] * e_scale
                q = UA_v * (medium.T(Q_tot - Q if counter else Q) - T)
            # dH/dV = q: split into dH/dT dT/dV + sum_j dH/dxi_j r_j + dH/dP dP/dV. The heat
            # capacity by a central difference: a forward one is off by about 3e-5, which
            # shifts the temperature profile by a millikelvin.
            h0 = self.h_phase(pkg, T, P, flows, hf)
            dT = 1e-4 * T
            cp = (self.h_phase(pkg, T + dT, P, flows, hf)
                  - self.h_phase(pkg, T - dT, P, flows, hf)) / (2.0 * dT)
            d_xi = 1e-7 * scale
            dh = [(self.h_phase(pkg, T, P, flows + d_xi * nu[j], hf) - h0) / d_xi
                  for j in range(m)]
            dP = 1e-4 * P
            dh_dp = (self.h_phase(pkg, T, P + dP, flows, hf) - h0) / dP if dPdV else 0.0
            dTdV = (q - float(np.dot(dh, r)) - dh_dp * dPdV) / cp
            if medium is None:
                return np.r_[r / scale, dTdV]
            return np.r_[r / scale, dTdV, q / e_scale]

        def integrate(medium: WallMedium | None = None, Q_tot: float = 0.0, t_eval=None,
                      dense: bool = False):
            y0 = np.zeros(m) if fixed_T is not None else np.r_[np.zeros(m), s.T]
            atol = np.r_[np.full(m, 1e-12), [1e-6]][: len(y0)]
            if medium is not None:
                y0, atol = np.r_[y0, 0.0], np.r_[atol, 1e-9]
            sol = solve_ivp(rhs, (0.0, V), y0, method="LSODA", t_eval=t_eval, rtol=self.RTOL,
                            atol=atol, dense_output=dense, args=(medium, Q_tot))
            if not sol.success:
                raise SimulationError(f"integration along the reactor failed: {sol.message}",
                                      where=self.id,
                                      hint="check the rate constants and units, or give a "
                                      "temperature")
            return sol

        grid = np.linspace(0.0, V, self.N_POINTS)
        medium, Q_tot = None, 0.0
        if wall:
            self._feed_T = s.T
            medium, Q_tot = self._shoot(pkg, integrate, e_scale, counter)
        sol = integrate(medium, Q_tot, t_eval=grid, dense=fixed_T is None)
        xis = sol.y[:m].T * scale
        flows_profile = n0 + xis @ nu
        T_profile = np.full(len(grid), fixed_T) if fixed_T is not None else sol.y[m]

        self.results["extents"] = {r.id: float(v) for r, v in zip(rxns, xis[-1])}
        self.results["residence_time"] = self.residence_time(pkg, s, V)
        extra: dict[str, tuple[str | None, list[float]]] = {}
        duty = None
        if medium is not None:
            Q_wall = sol.y[m + 1] * e_scale
            duty = float(Q_wall[-1])
            if counter and not medium.isothermal:
                miss = abs(duty - Q_tot)  # integration noise is of order RTOL
                if miss > max(1e-6, 10.0 * self.RTOL) * max(abs(duty), e_scale):
                    raise SimulationError(
                        "the countercurrent wall balance did not converge: the medium does "
                        f"not meet its inlet condition (off by {miss / 1e3:.3g} kW)",
                        where=self.id,
                        hint="a small medium flow along a long, well-cooled tube makes the "
                        "problem ill-conditioned; raise the medium flow or split the reactor")
            T_medium = np.array([medium.T(Q_tot - q if counter else q) for q in Q_wall])
            extra = {"Medium temperature": ("temperature", [float(t) for t in T_medium]),
                     "Heat through the wall": ("power", [float(q) for q in Q_wall])}
            self._wall_results(medium, duty, area, T_profile, T_medium)
        out = self.finish(pkg, s, flows_profile[-1], P_out, hf, duty)
        if fixed_T is None:
            self._hot_spot(sol, m, V, L)
        self._profiles(pkg, grid, V, flows_profile, T_profile, P_in, dPdV, extra)
        return out

    def _shoot(self, pkg, integrate, e_scale: float, counter: bool):
        """The medium along the wall and the total duty Q_tot it passes into the process."""
        ut, flow = self._utility, self.p["utility_flow"]

        def heat(medium: WallMedium, Q_tot: float) -> float:  # into the process, W
            return float(integrate(medium, Q_tot).y[-1, -1]) * e_scale

        if self._coolant is not None:
            fixed = stream_medium(pkg, self._coolant)
        elif ut.kind == "steam" or flow is not None:
            if ut.kind != "steam" and flow <= 0.0:
                raise SimulationError("the medium flow must be positive", where=self.id,
                                      hint="give a positive medium flow, or leave it empty")
            fixed = utility_medium(ut, flow)
        else:
            fixed = None
        if fixed is not None and (fixed.isothermal or not counter):
            return fixed, 0.0  # a single pass from the reactor inlet

        if fixed is not None:
            # Countercurrent: the medium's state at the reactor inlet is unknown. The search
            # starts where the medium leaves at the feed temperature, the most it could
            # exchange without a reaction. The residual falls with Q_tot (a medium that has
            # given up more is colder and takes more heat), so its sign there says which way
            # the root lies. Starting at Q_tot = 0 instead would ask the medium to come out
            # colder than it went in, where a process stream may condense.
            def medium_for(Q: float) -> WallMedium:
                return fixed

            def residual(Q: float) -> float:
                return heat(fixed, Q) - Q

            Q_ref = fixed.given_up(self._feed_T)
            a = Q_ref if abs(Q_ref) > 1e-6 * e_scale else 0.0
            fa = residual(a)
            if fa == 0.0:
                return fixed, a
            step = math.copysign(abs(Q_ref) / 2.0 if a else e_scale, fa)
            found = _bracket(residual, a, fa, a + step)
            if found is None:
                raise SimulationError(
                    "no medium outlet temperature meets the medium inlet temperature in "
                    "countercurrent", where=self.id,
                    hint="check U, the area and the medium, or try cocurrent flow")
        else:
            # A fluid utility at the flow that brings it to its outlet temperature: that flow
            # is Q_tot / (cp (T_in - T_out)). With an unlimited flow the medium stays at T_in
            # and passes the most heat, Q_inf; the duty has its sign and lies near it or
            # below it in magnitude, so it is searched on ln|Q_tot| from there.
            Q_inf = heat(WallMedium(ut.describe(), ut.T_in), 0.0)
            if Q_inf == 0.0 or (Q_inf > 0.0) != ut.heats:
                raise SimulationError(
                    f"temperature cross: {ut.describe()} cannot "
                    f"{'heat' if ut.heats else 'cool'} the reactor along the tube",
                    where=self.id, hint="choose a utility with more temperature difference")

            def medium_for(Q: float) -> WallMedium:
                return utility_medium(ut, Q=Q)

            def residual(Q: float) -> float:
                return heat(medium_for(Q), Q) - Q

            def in_log(x: float) -> float:
                return residual(math.copysign(math.exp(x), Q_inf))

            x0 = math.log(abs(Q_inf))
            f_inf = in_log(x0)
            step = math.log(2.0)
            found = (_bracket(in_log, x0, f_inf, x0 - step, tries=40, grow=False)
                     or _bracket(in_log, x0, f_inf, x0 + step, tries=20, grow=False))
            if found is None:
                raise SimulationError(
                    f"{ut.describe()} cannot reach its outlet temperature along the reactor",
                    where=self.id, hint="give the medium flow, so its outlet temperature "
                    "follows, or choose another outlet temperature")
            found = tuple((math.copysign(math.exp(x), Q_inf), f) for x, f in found)
        (a, fa), (b, fb) = found
        if fa == 0.0 or fb == 0.0:
            Q = a if fa == 0.0 else b
            return medium_for(Q), Q
        Q = brentq(residual, min(a, b), max(a, b), rtol=1e-12,
                   xtol=1e-10 * max(abs(a), abs(b)) + 1e-9 * e_scale)
        return medium_for(Q), Q

    def _wall_results(self, medium: WallMedium, Q: float, area: float, T, T_medium) -> None:
        UA = self.p["U"] * area
        self.results.update(wall_duty=Q, UA=UA, heat_transfer_area=area,
                            LMTD=abs(Q) / UA if UA > 0.0 else None)
        if self._coolant is not None:
            self._wall_heat = Q
            return
        ut = self._utility
        if medium.heats_only and Q < 0.0:
            raise SimulationError(
                f"{ut.describe()} can only heat, but the reactor gives {-Q / 1e3:.4g} kW to it",
                where=self.id, hint="choose a cooling utility")
        if medium.heats_only and (T > T_medium + 1e-9).any():
            self.notes.append(f"the reactor is hotter than {ut.describe()} over part of its "
                              "length; the jacket takes heat there at the saturation "
                              "temperature")
        self.results.update(utility_flow=utility_mass_flow(ut, medium, Q),
                            utility_T_out=medium.T(Q))

    def _hot_spot(self, sol, m: int, V: float, L: float | None) -> None:
        """Highest temperature along the tube and where it is, from the dense solution."""
        v = np.linspace(0.0, V, 801)
        T = sol.sol(v)[m]
        k = int(np.argmax(T))
        v_hot, T_hot = float(v[k]), float(T[k])
        if 0 < k < len(v) - 1:
            res = minimize_scalar(lambda x: -float(sol.sol(x)[m]), bounds=(v[k - 1], v[k + 1]),
                                  method="bounded", options={"xatol": 1e-10 * V})
            if -res.fun > T_hot:
                v_hot, T_hot = float(res.x), float(-res.fun)
        self.results["T_hot_spot"] = T_hot
        self.results["hot_spot_position"] = v_hot / V * L if L is not None else None

    def _profiles(self, pkg, grid, V, flows, T, P_in, dPdV, extra=None) -> None:
        if self.p["length"] is not None:
            x, x_label = grid / V * self.p["length"], "Length [m]"
        else:
            x, x_label = grid, "Volume [m³]"
        x = [float(v) for v in x]
        F = np.maximum(flows, 0.0)
        z = F / F.sum(axis=1, keepdims=True)
        series: dict[str, tuple[str | None, list[float]]] = {}
        key = self.key_index(pkg)
        if key is not None and flows[0, key] > 0.0:
            series[f"Conversion of {pkg.names[key]}"] = (
                "fraction", [float(1.0 - f / flows[0, key]) for f in flows[:, key]])
        series["Temperature"] = ("temperature", [float(t) for t in T])
        series.update(extra or {})
        if dPdV:
            series["Pressure"] = ("pressure", [float(P_in + dPdV * v) for v in grid])
        self.add_profile("conditions", "Along the reactor", x_label, x, series)
        present = [i for i in range(pkg.n) if F[:, i].max() > 0.0]
        self.add_profile(
            "composition", "Composition along the reactor", x_label, x,
            {f"{pkg.names[i]} (mole fraction)": ("fraction", [float(v) for v in z[:, i]])
             for i in present},
        )


def _bracket(residual, a: float, fa: float, b: float, tries: int = 60, grow: bool = True):
    """Two points around a sign change of ``residual``, walking from ``a`` (where it is
    ``fa``) through ``b`` and on in the same direction; None if there is none within
    ``tries`` steps. Each step doubles if ``grow``. Where the residual cannot be evaluated
    (the reactor runs away, a flash fails), the step is halved instead.
    """
    step, bad = b - a, None
    for _ in range(tries):
        try:
            fb = residual(b)
        except SimulationError:
            bad = b
            step /= 2.0
            b = a + step
            continue
        if fb == 0.0 or (fb > 0.0) != (fa > 0.0):
            return (a, fa), (b, fb)
        a, fa = b, fb
        if bad is not None:
            step = (bad - a) / 2.0
        elif grow:
            step *= 2.0
        b = a + step
    return None
