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
from scipy.optimize import brentq, root

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
from .base import (
    DP,
    DUTY,
    IN,
    OUT,
    P_OUT,
    Param,
    ResultField,
    UnitOperation,
    empty_stream,
    single_inlet,
)

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

T_OUT = ResultField("T_out", "Outlet temperature", "temperature")
CONVERSION = ResultField("conversion", "Conversion of key component", "fraction")
RESIDENCE_TIME = ResultField("residence_time", "Residence time (inlet conditions)", "time")
ENERGY_GROUPS = (("T", "duty"), ("P", "dP"))


class Reactor(UnitOperation):
    """Shared energy balance and bookkeeping. Not a unit type itself (no ``type``)."""

    category = "Reactors"
    inlets = (IN,)
    outlets = (OUT,)
    exclusive_groups = ENERGY_GROUPS
    result_fields = (T_OUT, DUTY, CONVERSION)
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

    def finish(self, pkg, s: Stream, flows, P: float, hf: np.ndarray) -> Stream:
        """Outlet stream from the energy specification; records duty, T_out, conversion."""
        flows = self.clean(np.asarray(flows, dtype=float), s.total_flow)
        if self.p["T"] is not None:
            out = Stream.from_flash(pkg, flows, T=self.p["T"], P=P)
        else:
            target = self.h_flow(s, hf) + (self.p["duty"] or 0.0)
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
        T_SPEC, DUTY_SPEC, P_OUT, DP,
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
    result_fields = (T_OUT, DUTY)
    params = (
        Param("yields", "Mass yields", "component_fractions",
              help="Relative mass of each product; normalized so the outlet mass equals the "
              "reacting inlet mass."),
        Param("inerts", "Inert components", "string_list",
              help="Components that pass through unchanged."),
        T_SPEC, DUTY_SPEC, P_OUT, DP,
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
    params = (REACTIONS, PHASE, KEY, T_SPEC, DUTY_SPEC, P_OUT, DP)

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
    description = ("Minimizes the Gibbs energy of an ideal-gas mixture subject to element "
                   "balances. Needs no reactions, only which compounds may form.")
    params = (
        Param("products", "Reacting components", "string_list",
              help="Components that may form or react. Others pass through unchanged. "
              "Default: all components."),
        KEY, T_SPEC, DUTY_SPEC, P_OUT, DP,
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

        def flows_at(T: float) -> np.ndarray:
            if not len(nu):
                return n0
            g = gibbs_over_rt(pkg, T, species, why="the Gibbs reactor")
            ln_k = -nu[:, species] @ g[species]
            return n0 + equilibrium_extents(n0, nu, ln_k, ln_p) @ nu

        T = self.p["T"] if self.p["T"] is not None else (
            self.temperature_for_energy(pkg, s, P, hf, flows_at))
        return self.finish(pkg, s, flows_at(T), P, hf)

    def key_index(self, pkg):
        return super().key_index(pkg) if self.p.get("key") else None


# ---- kinetic reactors ----------------------------------------------------------------------

class KineticReactor(Reactor):
    """CSTR and PFR: rates from the reactions' kinetics."""

    result_fields = (T_OUT, DUTY, CONVERSION, RESIDENCE_TIME)

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
        PHASE, KEY, T_SPEC, DUTY_SPEC, P_OUT, DP,
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


class PFR(KineticReactor):
    type = "PFR"
    label = "Plug Flow Reactor"
    description = ("Tubular reactor without back-mixing: the balances are integrated along the "
                   "volume. Gives conversion, temperature and composition profiles.")
    params = (
        REACTIONS,
        Param("V", "Volume", "volume", min=0.0, help="Reactor volume."),
        Param("length", "Length", "length", min=0.0, help="Tube length; needs the diameter."),
        Param("diameter", "Diameter", "length", min=0.0,
              help="Tube inner diameter, used with the length."),
        PHASE, KEY, T_SPEC, DUTY_SPEC, P_OUT, DP,
    )
    exclusive_groups = (("V", "length"), *ENERGY_GROUPS)
    required_groups = (("V", "length"),)
    N_POINTS = 41

    def react(self, pkg, s, P_out):
        rxns = self.setup(pkg)
        nu, m = self._nu, len(rxns)
        if self.p["length"] is not None:
            if not self.p["diameter"]:
                raise SimulationError("the length is given but not the diameter", where=self.id,
                                      hint="set the diameter, or give the volume instead")
            V = math.pi / 4.0 * self.p["diameter"] ** 2 * self.p["length"]
        else:
            V = self.p["V"]
        if not V or V <= 0.0:
            raise SimulationError("the reactor volume must be positive", where=self.id)
        n0 = np.array(s.flows)
        scale = s.total_flow
        species = sorted({i for r in rxns for i in r.involved})
        hf = self.formation_enthalpies(pkg, [n0], species)
        fixed_T = self.p["T"]
        q_v = (self.p["duty"] or 0.0) / V
        P_in = s.P
        dPdV = (P_out - P_in) / V

        def rhs(vol, y):
            xi = y[:m] * scale
            T = fixed_T if fixed_T is not None else y[m]
            P = P_in + dPdV * vol
            flows = n0 + xi @ nu
            r = self.rates(pkg, T, P, flows)
            if fixed_T is not None:
                return r / scale
            # dH/dV = q: split into dH/dT dT/dV + sum_j dH/dxi_j r_j + dH/dP dP/dV.
            h0 = self.h_phase(pkg, T, P, flows, hf)
            dT = 1e-4 * T
            cp = (self.h_phase(pkg, T + dT, P, flows, hf) - h0) / dT
            d_xi = 1e-7 * scale
            dh = [(self.h_phase(pkg, T, P, flows + d_xi * nu[j], hf) - h0) / d_xi
                  for j in range(m)]
            dP = 1e-4 * P
            dh_dp = (self.h_phase(pkg, T, P + dP, flows, hf) - h0) / dP if dPdV else 0.0
            dTdV = (q_v - float(np.dot(dh, r)) - dh_dp * dPdV) / cp
            return np.r_[r / scale, dTdV]

        y0 = np.zeros(m) if fixed_T is not None else np.r_[np.zeros(m), s.T]
        grid = np.linspace(0.0, V, self.N_POINTS)
        sol = solve_ivp(rhs, (0.0, V), y0, method="LSODA", t_eval=grid, rtol=1e-8,
                        atol=np.r_[np.full(m, 1e-12), [1e-6]][: len(y0)])
        if not sol.success:
            raise SimulationError(f"integration along the reactor failed: {sol.message}",
                                  where=self.id,
                                  hint="check the rate constants and units, or give a "
                                  "temperature")
        xis = sol.y[:m].T * scale
        flows_profile = n0 + xis @ nu
        T_profile = np.full(len(grid), fixed_T) if fixed_T is not None else sol.y[m]

        self.results["extents"] = {r.id: float(v) for r, v in zip(rxns, xis[-1])}
        self.results["residence_time"] = self.residence_time(pkg, s, V)
        out = self.finish(pkg, s, flows_profile[-1], P_out, hf)
        self._profiles(pkg, grid, V, flows_profile, T_profile, P_in, dPdV)
        return out

    def _profiles(self, pkg, grid, V, flows, T, P_in, dPdV) -> None:
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
        if dPdV:
            series["Pressure"] = ("pressure", [float(P_in + dPdV * v) for v in grid])
        self.add_profile("conditions", "Along the reactor", x_label, x, series)
        present = [i for i in range(pkg.n) if F[:, i].max() > 0.0]
        self.add_profile(
            "composition", "Composition along the reactor", x_label, x,
            {f"{pkg.names[i]} (mole fraction)": ("fraction", [float(v) for v in z[:, i]])
             for i in present},
        )
