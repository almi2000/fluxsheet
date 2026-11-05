"""Batch reactor: a steady-state flowsheet unit that integrates one batch cycle.

The continuous feed is collected into charges, as in Aspen Plus' RBatch: each cycle takes
the charge N0 = F_in * t_cycle, fills the vessel (``charge_time``), reacts until a stop
criterion is met (at most ``batch_time``), and is discharged and cleaned (``down_time``).
The product leaves as the average stream N_end / t_cycle, so the unit fits into continuous
flowsheets and conserves mass against its feed.

The cycle time is given, follows from a given charge (t_cycle = charge / F_in), or, with
neither, is that of one vessel run back to back (t_cycle = charge + reaction + down time).
If the vessel is occupied longer than the cycle time, several vessels are needed in turn;
the unit reports how many.

Balances (enthalpies include formation enthalpies, so the heat of reaction comes out of
the balance by itself; kinetics as in the CSTR, ``Reaction.rate`` of the reacting phase):

- constant pressure: component amounts N and contents enthalpy E, dN/dt = G - F_v y_v,
  dE/dt = Q - F_v h_v. Without a vent the contents are flashed (one or two phases; vapor
  stays in the vessel) and the reaction runs in the chosen phase. With a vent the contents
  are liquid: subcooled (F_v = 0) until they reach their bubble point, then boiling, held at
  the bubble point T_b(x) with the vapor y* in equilibrium leaving as it forms. The liquid
  enthalpy H_L(N) = N h_L(T_b(x), x) then fixes the vent flow,
  F_v = (dH_L/dN . G - Q) / (dH_L/dN . y* - h_v). Events switch between the two.
- constant volume (closed rigid vessel, gas phase): N and internal energy U,
  dU/dt = Q; T and P follow from (U, V, N). The vessel volume is that of the charge at the
  feed temperature and the vessel's initial pressure.

Energy: adiabatic, fixed duty, isothermal, a temperature program (time, T table) or a
jacket with a utility (Q = UA * LMTD against the medium, see ``Utility.rate``).

The batch is intensive: it is integrated for a reference charge of 1 mol of feed
composition, and amounts and heat are scaled by S = N0 / 1 mol afterwards. Extensive
inputs (duty, UA, a fixed medium flow) are divided by S.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from ..contents_flash import ContentsFlash
from ..errors import SimulationError
from ..reactions import R
from ..stream import Stream
from .base import IN, Param, Port, ResultField, empty_stream, single_inlet
from .heat_transfer import AREA, HEAT_PARAMS
from .reactors import DUTY_SPEC, KEY, REACTIONS, T_SPEC, KineticReactor

N_POINTS = 101  # profile points
N_REF = 1.0  # mol, reference charge of the integration


@dataclass
class _Snap:
    """The contents at one instant (reference charge)."""

    T: float
    P: float
    Nt: float
    z: np.ndarray
    VF: float
    gen: np.ndarray  # mol/s formed per component
    n_vap: float = 0.0
    y_vap: np.ndarray | None = None
    h_vap: float = 0.0  # J/mol of vapor, with formation enthalpies
    H: float = 0.0  # contents enthalpy with formation enthalpies, J
    U: float = 0.0  # contents internal energy (constant volume), J
    volume: float = 0.0  # m3
    reacting: float = 0.0  # mol in the reacting phase
    fv: float = 0.0  # vent flow, mol/s


T_PROFILE = Param(
    "T_profile", "Temperature program", "records", tab="Program",
    fields=(Param("time", "Time", "time", min=0.0), Param("T", "Temperature", "temperature",
                                                          min=1.0)),
    help="Temperature over the reaction time: linear between the points, held after the "
    "last. The contents are brought to the first temperature at the start.",
)


class BatchReactor(KineticReactor):
    type = "BatchReactor"
    label = "Batch Reactor"
    description = ("Stirred batch vessel: each cycle the feed collected over the cycle time "
                   "is charged and reacts until a stop criterion; the product leaves as an "
                   "average stream. Time profiles of T, P, composition and duty.")
    inlets = (IN,)
    outlets = (Port("out", "Product (cycle average)"), Port("vent", "Vent"))
    params = (
        REACTIONS,
        Param("phase", "Reacting phase", "choice", default="liquid", options=("liquid", "gas"),
              help="Phase the reactions take place in; sets concentrations and activities."),
        KEY,
        Param("product", "Product component", "component",
              help="Component whose production rate is reported as the productivity. "
              "Defaults to the first product of the first reaction."),
        Param("batch_time", "Batch time", "time", min=0.0,
              help="Reaction time; with a stop criterion, the longest the batch may take."),
        Param("stop_conversion", "Stop at conversion", "fraction", min=0.0, max=1.0,
              tab="Conversion",
              help="Stop when the key component reaches this conversion."),
        Param("stop_mole_fraction", "Stop at mole fraction", "fraction", min=0.0, max=1.0,
              tab="Mole fraction",
              help="Stop when the mole fraction of the stop component in the contents "
              "reaches this value (rising or falling)."),
        Param("stop_component", "Stop component", "component", requires="stop_mole_fraction",
              help="Component of the mole-fraction stop criterion."),
        Param("stop_T", "Stop at temperature", "temperature", min=0.0, tab="Temperature",
              help="Stop when the contents reach this temperature (rising or falling)."),
        Param("cycle_time", "Cycle time", "time", min=0.0, tab="Cycle time",
              help="Time between two charges; the charge is the feed over one cycle. "
              "Default: one vessel back to back (charge + reaction + down time)."),
        Param("charge", "Charge per batch", "amount", min=0.0, tab="Charge",
              help="Amount charged per batch; the cycle time is the charge over the feed "
              "flow."),
        Param("charge_time", "Charge time", "time", min=0.0, default=0.0,
              help="Time to fill the vessel; no reaction is counted while charging."),
        Param("down_time", "Down time", "time", min=0.0, default=0.0,
              help="Discharge and cleaning time per cycle."),
        Param("vessel", "Vessel", "choice", default="constant_pressure",
              options=("constant_pressure", "constant_volume"),
              help="constant_pressure: the contents stay at the vessel pressure; "
              "constant_volume: a closed rigid vessel (gas phase), the pressure follows."),
        Param("P", "Vessel pressure", "pressure", min=1.0,
              help="Pressure of the contents (initial pressure of a constant-volume "
              "vessel). Defaults to the feed pressure."),
        Param("vent", "Vent vapor", "choice", default="no", options=("no", "yes"),
              help="Constant pressure, liquid phase: vapor that forms leaves into the vent "
              "outlet as it forms."),
        T_SPEC, T_PROFILE, DUTY_SPEC, AREA, *HEAT_PARAMS,
    )
    exclusive_groups = (("T", "T_profile", "duty", "area"),
                        ("stop_conversion", "stop_mole_fraction", "stop_T"),
                        ("cycle_time", "charge"))
    required_groups = (("batch_time",),)
    result_fields = (
        ResultField("T_out", "Final temperature", "temperature"),
        ResultField("P_final", "Final pressure", "pressure"),
        ResultField("conversion", "Conversion of key component", "fraction"),
        ResultField("batch_time", "Reaction time", "time"),
        ResultField("cycle_time", "Cycle time", "time"),
        ResultField("charge", "Charge per batch", "amount"),
        ResultField("vessels_needed", "Vessels needed", None),
        ResultField("duty", "Duty (cycle average)", "power"),
        ResultField("heat_per_batch", "Heat per batch", "energy"),
        ResultField("peak_duty", "Peak duty", "power"),
        ResultField("productivity", "Productivity", "mass_flow"),
        ResultField("vent_amount", "Vented per batch", "amount"),
        ResultField("max_volume", "Largest contents volume", "volume"),
        ResultField("utility_flow", "Medium flow (cycle average)", "mass_flow"),
        ResultField("UA", "UA", "thermal_conductance"),
        ResultField("heat_transfer_area", "Heat transfer area", "area"),
    )
    def __post_init__(self) -> None:
        super().__post_init__()
        self._medium = None
        self._cache: tuple | None = None
        if self.p["T_profile"] is not None:
            pts = self.p["T_profile"]
            if not pts or any(p["time"] is None or p["T"] is None for p in pts):
                raise SimulationError("every point of the temperature program needs a time "
                                      "and a temperature", where=self.id)
            times = [p["time"] for p in pts]
            if any(b <= a for a, b in itertools.pairwise(times)):
                raise SimulationError("the times of the temperature program must increase",
                                      where=self.id, hint="sort the points by time")

    # The solver's utilities.calculate() iterates the duty of a steady jacket. A batch rates
    # its jacket over time itself, so the medium utilities.attach() resolves is kept here
    # and the wrapper sees none.
    @property
    def utility(self):
        return None

    @utility.setter
    def utility(self, medium) -> None:
        self._medium = medium

    def idle_outlets(self):
        return ("vent",) if self.p["vent"] == "no" else ()

    # ---- entry -------------------------------------------------------------------------

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        key = (tuple(s.flows), s.T, s.P)
        if self._cache is not None and self._cache[0] == key:  # e.g. recycle passes
            _, outputs, results, profiles, notes = self._cache
            self.results, self.profiles, self.notes = dict(results), dict(profiles), list(notes)
            return dict(outputs)
        self.results, self.profiles = {}, {}
        outputs = self._calculate(pkg, s)
        self._cache = (key, dict(outputs), dict(self.results), dict(self.profiles),
                       list(self.notes))
        return outputs

    def _calculate(self, pkg, s: Stream):
        self.setup(pkg)  # kinetics required, as in the CSTR
        p = self.p
        self._pkg = pkg
        self._const_v = p["vessel"] == "constant_volume"
        self._vented = p["vent"] == "yes"
        if self._vented and (self._const_v or self._phase != "liquid"):
            raise SimulationError(
                "a vent needs a constant-pressure vessel with a liquid reacting phase",
                where=self.id, hint="set the vent to no, or the vessel to constant pressure "
                "and the reacting phase to liquid")
        if self._vented and (p["T"] is not None or p["T_profile"] is not None):
            raise SimulationError(
                "a vented batch boils at its bubble point, so its temperature cannot be "
                "prescribed", where=self.id,
                hint="give a duty or a jacket (or neither, adiabatic) with the vent")
        if self._const_v and self._phase != "gas":
            raise SimulationError(
                "a constant-volume vessel is modelled for gas-phase reactions only",
                where=self.id, hint="set the reacting phase to gas, or use constant pressure "
                "(the liquid then stays at the vessel pressure)")
        if self._medium is not None and p["area"] is None:
            raise SimulationError(
                "the batch reactor uses a utility through a jacket only", where=self.id,
                hint="give U and the heat transfer area of the jacket, or clear the utility")
        P = p["P"] if p["P"] is not None else s.P
        F = s.total_flow
        if F <= 0.0:
            self.results.update(duty=0.0, T_out=s.T)
            return {"out": empty_stream(pkg, s.T, P), "vent": empty_stream(pkg, s.T, P)}

        z0 = self._z0 = np.asarray(s.zs)
        self._boiling = False
        species = sorted({i for r in self._rxns for i in r.involved})
        self._hf = self.formation_enthalpies(pkg, [s.flows], species)
        self._prog = self._program()
        self._stop = self._stop_criterion(pkg, z0)
        t_max = p["batch_time"]
        if not t_max or t_max <= 0.0:
            raise SimulationError("the batch time must be positive", where=self.id)

        # Charge size before integrating: needed when the heat input is extensive.
        S = self._scale(F, t_max if self._stop is None else None)
        extensive = p["duty"] is not None or p["area"] is not None
        if extensive and S is None:
            raise SimulationError(
                "with a fixed duty or a jacket and a stop criterion, the charge must be known "
                "before the batch runs", where=self.id,
                hint="give the cycle time or the charge per batch")
        self._S = S
        self._P = P
        self._flash = ContentsFlash(pkg, s.T)
        h_feed = s.state.H + float(z0 @ self._hf)  # J/mol with formation enthalpies
        y0 = self._initial(pkg, s, z0, P, h_feed)
        sol, t_stop = self._integrate(y0, t_max)

        if S is None:
            S = self._scale(F, t_stop)
        self._S = S
        return self._finish(pkg, s, z0, sol, t_stop, S, h_feed)

    # ---- set-up ------------------------------------------------------------------------

    def _program(self):
        """(T(t), breakpoints) for prescribed temperatures, None with an energy balance."""
        if self.p["T"] is not None:
            T = self.p["T"]
            return (lambda t: T), []
        pts = self.p["T_profile"]
        if pts is None:
            return None
        times = np.array([pt["time"] for pt in pts])
        temps = np.array([pt["T"] for pt in pts])
        return (lambda t: float(np.interp(t, times, temps))), [float(t) for t in times if t > 0]

    def _stop_criterion(self, pkg, z0):
        p = self.p
        if p["stop_conversion"] is not None:
            x = p["stop_conversion"]
            if not 0.0 < x < 1.0:
                raise SimulationError(f"the stop conversion must be between 0 and 1, got {x:g}",
                                      where=self.id)
            k = self.key_index(pkg)
            if z0[k] <= 0.0:
                raise SimulationError(
                    f"the key component {pkg.names[k]!r} is not in the feed, so its conversion "
                    "cannot reach the stop value", where=self.id,
                    hint="choose a key component that is fed")
            return ("conversion", k, x)
        if p["stop_mole_fraction"] is not None:
            if not p["stop_component"]:
                raise SimulationError("the mole-fraction stop criterion needs its component",
                                      where=self.id, hint="set the stop component")
            k = self.component_index(pkg, p["stop_component"], "stop component")
            return ("mole_fraction", k, p["stop_mole_fraction"])
        if p["stop_T"] is not None:
            if p["T"] is not None:
                raise SimulationError(
                    "an isothermal batch never reaches a stop temperature", where=self.id,
                    hint="stop by time, conversion or mole fraction, or use a temperature "
                    "program or an energy balance")
            return ("temperature", None, p["stop_T"])
        return None

    def _scale(self, F: float, t_react: float | None) -> float | None:
        """S = N0 / N_REF from the cycle basis; None while the reaction time is unknown."""
        p = self.p
        if p["cycle_time"] is not None:
            if p["cycle_time"] <= 0.0:
                raise SimulationError("the cycle time must be positive", where=self.id)
            return F * p["cycle_time"] / N_REF
        if p["charge"] is not None:
            if p["charge"] <= 0.0:
                raise SimulationError("the charge must be positive", where=self.id)
            return p["charge"] / N_REF
        if t_react is None:
            return None
        return F * (p["charge_time"] + t_react + p["down_time"]) / N_REF

    def _initial(self, pkg, s: Stream, z0, P: float, h_feed: float) -> np.ndarray:
        n = pkg.n
        N0 = N_REF * z0
        self._n = n
        energy = self._prog is None
        parts = [N0]
        if self._const_v:
            st = pkg.flash(list(z0), T=s.T, P=P)
            if st.VF < 1.0 - 1e-9:
                raise SimulationError(
                    "a constant-volume vessel needs a gas charge, but the feed is not all "
                    f"vapor at {s.T - 273.15:.1f} °C and {P / 1e5:.4g} bar", where=self.id,
                    hint="heat the feed, lower the vessel pressure, or use constant pressure")
            self._V = N_REF * st.V
            self._Pg, self._Tg = P, s.T
            # Filling to P0 at the feed temperature, then closed: U0 = H(T, P0) - P0 V.
            self._U_charge = N_REF * (st.H + float(z0 @ self._hf)) - P * self._V
            if energy:
                parts.append([self._U_charge])
        elif self._vented:
            # The charge enters with the feed's enthalpy; vapor it flashes to leaves at once.
            st = pkg.flash(list(z0), P=P, H=h_feed - float(z0 @ self._hf))
            if st.VF >= 1.0 - 1e-9:
                raise SimulationError("the charge is all vapor at the vessel pressure",
                                      where=self.id, hint="raise the vessel pressure, or cool "
                                      "the feed")
            self._boiling = st.VF > 0.0
            self._T_last, self._y_last = st.T, np.asarray(st.y or z0)
            if self._boiling:
                xl, yv = np.asarray(st.x), np.asarray(st.y)
                h_v = pkg.phase_props("vapor", st.T, P, list(yv)).H + float(yv @ self._hf)
                h_l = pkg.phase_props("liquid", st.T, P, list(xl)).H + float(xl @ self._hf)
                n_v = N_REF * st.VF
                parts = [N_REF * z0 - n_v * yv, [(N_REF - n_v) * h_l], n_v * yv, [n_v * h_v]]
            else:
                parts += [[N_REF * h_feed], np.zeros(n), [0.0]]
        elif energy:
            parts.append([N_REF * h_feed])  # the charge enters with the feed's enthalpy
        if energy:
            parts.append([0.0])  # heat added so far
        self._energy = energy
        k = n
        self._iE = k if energy else None
        k += energy
        self._iNv = slice(k, k + n) if self._vented else None
        self._iHv = k + n if self._vented else None
        k += (n + 1) if self._vented else 0
        self._iQ = k if energy else None
        y0 = np.concatenate([np.asarray(v, dtype=float) for v in parts])
        RT = R * max(s.T, 200.0) * N_REF
        atol = np.full(len(y0), 1e-10 * N_REF)
        if energy:
            atol[self._iE] = atol[self._iQ] = 1e-6 * RT
        if self._vented:
            atol[self._iHv] = 1e-6 * RT
        self._atol = atol
        return y0

    # ---- contents ----------------------------------------------------------------------

    def _heat(self, T: float) -> float:
        """Heat into the reference charge, W."""
        p, S = self.p, self._S
        if p["area"] is not None:
            flow = p["utility_flow"] / S if p["utility_flow"] is not None else None
            return self._medium.rate(p["U"] * p["area"] / S, T, T, "countercurrent", flow)
        if p["duty"] is not None:
            return p["duty"] / S
        return 0.0

    def _rates(self, T, P, comp, Vm, n) -> np.ndarray:
        if n <= 0.0:
            return np.zeros(self._n)
        r = np.array([rx.rate(self._pkg, T, P, comp, Vm, self._phase) for rx in self._rxns])
        return n * Vm * (r @ self._nu)

    def _snap(self, t: float, y: np.ndarray) -> _Snap:
        N = np.maximum(y[: self._n], 0.0)
        Nt = N.sum()
        if Nt <= 1e-12 * N_REF:
            raise SimulationError("the vessel ran empty", where=self.id,
                                  hint="vent less: lower the duty or the temperature")
        z = N / Nt
        if self._const_v:
            return self._snap_closed(t, y, N, Nt, z)
        if self._vented:
            return self._snap_vented(y, N, Nt, z)
        pkg, P, hf = self._pkg, self._P, self._hf
        if self._energy:
            st = self._flash(z, P, y[self._iE] / Nt - float(z @ hf))
        else:
            st = self._flash.at_temperature(z, P, self._prog[0](t))
        T, VF = st.T, min(max(st.VF, 0.0), 1.0)
        xl = np.asarray(st.x) if st.x and VF < 1.0 else z
        yv = np.asarray(st.y) if st.y and VF > 0.0 else z
        n_vap, n_liq = Nt * VF, Nt * (1.0 - VF)
        if VF <= 0.0 or VF >= 1.0:  # one phase: the flash has its molar volume
            Vm_l = Vm_v = st.V
        else:
            Vm_l = pkg.phase_props("liquid", T, P, list(xl)).V
            Vm_v = pkg.phase_props("vapor", T, P, list(yv)).V
        if self._phase == "liquid":
            gen = self._rates(T, P, xl, Vm_l, n_liq)
            reacting = n_liq
        else:
            gen = self._rates(T, P, yv, Vm_v, n_vap)
            reacting = n_vap
        return _Snap(T=T, P=P, Nt=Nt, z=z, VF=VF, gen=gen, n_vap=n_vap, y_vap=yv,
                     H=Nt * (st.H + float(z @ hf)),
                     volume=n_liq * (Vm_l if n_liq > 0 else 0.0) + n_vap * Vm_v,
                     reacting=reacting)

    def _pv(self, T: float, Nt: float, z: np.ndarray) -> tuple[float, float]:
        """Pressure that fills the vessel volume at T, and the molar enthalpy there."""
        pkg, V = self._pkg, self._V
        P = self._Pg * (T / self._Tg)  # last pressure, carried along ideally
        for _ in range(50):
            pp = pkg.phase_props("vapor", T, P, list(z))
            f = Nt * pp.V - V
            dP = 1e-6 * P
            dfdP = (Nt * pkg.phase_props("vapor", T, P + dP, list(z)).V - Nt * pp.V) / dP
            if dfdP >= 0.0:
                raise SimulationError("no gas-phase pressure fills the vessel", where=self.id)
            step = -f / dfdP
            P_new = min(max(P + step, 0.2 * P), 5.0 * P)
            if abs(P_new - P) <= 1e-12 * P:
                P = P_new
                break
            P = P_new
        else:
            raise SimulationError("the vessel pressure did not converge", where=self.id)
        self._Pg, self._Tg = P, T
        return P, pkg.phase_props("vapor", T, P, list(z)).H

    def _internal_energy(self, T, Nt, z) -> tuple[float, float]:
        P, h = self._pv(T, Nt, z)
        return Nt * (h + float(z @ self._hf)) - P * self._V, P

    def _snap_closed(self, t, y, N, Nt, z) -> _Snap:
        if self._energy:
            U = y[self._iE]
            T = getattr(self, "_T_last", self._Tg)
            for _ in range(60):  # Newton on T for the internal energy
                u0, P = self._internal_energy(T, Nt, z)
                dT = 1e-4 * T
                cv = (self._internal_energy(T + dT, Nt, z)[0] - u0) / dT
                step = (U - u0) / cv
                step = max(min(step, 100.0), -100.0)
                T = max(T + step, 20.0)
                if abs(step) <= 1e-11 * T:
                    break
            else:
                raise SimulationError("could not find the vessel temperature", where=self.id)
            self._T_last = T
            u, P = self._internal_energy(T, Nt, z)
        else:
            T = self._prog[0](t)
            u, P = self._internal_energy(T, Nt, z)
        Vm = self._V / Nt
        gen = self._rates(T, P, z, Vm, Nt)
        return _Snap(T=T, P=P, Nt=Nt, z=z, VF=1.0, gen=gen, U=u, H=u + P * self._V,
                     volume=self._V, reacting=Nt)

    # ---- vented contents: subcooled liquid or boiling at the bubble point ---------------

    def _h_liquid(self, T: float, x: np.ndarray) -> float:
        """Molar liquid enthalpy with formation enthalpies, J/mol."""
        return self._pkg.phase_props("liquid", T, self._P, list(x)).H + float(x @ self._hf)

    def _liquid_T(self, x: np.ndarray, h: float) -> float:
        T = self._T_last
        for _ in range(60):
            h0 = self._h_liquid(T, x)
            cp = (self._h_liquid(T + 1e-3, x) - h0) / 1e-3
            if cp <= 0.0:
                break
            step = max(min((h - h0) / cp, 50.0), -50.0)
            T = max(T + step, 20.0)
            if abs(step) <= 1e-11 * T:
                return T
        raise SimulationError("could not find the liquid temperature", where=self.id)

    def _bubble(self, x: np.ndarray, T: float | None = None, y: np.ndarray | None = None):
        """Bubble temperature and vapor of liquid x at the vessel pressure, by Newton steps
        on ln(sum K x) with K from the phase models (cheap near the last bubble point)."""
        pkg, P = self._pkg, self._P
        T = self._T_last if T is None else T
        y = self._y_last if y is None else y
        xs = list(x)
        for _ in range(60):
            K = np.array(pkg.k_values(T, P, xs, list(y)))
            s = float(K @ x)
            dT = 1e-4 * T
            s2 = float(np.array(pkg.k_values(T + dT, P, xs, list(y))) @ x)
            slope = (math.log(s2) - math.log(s)) / dT
            y_new = K * x / s
            if slope <= 0.0:
                break
            step = max(min(-math.log(s) / slope, 30.0), -30.0)
            T += step
            moved = float(np.abs(y_new - y).max())
            y = y_new
            if abs(step) <= 1e-11 * T and moved <= 1e-12:
                return T, y
        raise SimulationError("could not find the bubble point of the contents",
                              where=self.id, hint="check the vessel pressure")

    def _snap_vented(self, y, N, Nt, z) -> _Snap:
        P = self._P
        if self._boiling:
            T, yv = self._bubble(z)
        else:
            T, yv = self._liquid_T(z, y[self._iE] / Nt), z
        self._T_last = T
        Vm = self._pkg.phase_props("liquid", T, P, list(z)).V
        gen = self._rates(T, P, z, Vm, Nt)
        snap = _Snap(T=T, P=P, Nt=Nt, z=z, VF=0.0, gen=gen, y_vap=yv,
                     H=Nt * self._h_liquid(T, z), volume=Nt * Vm, reacting=Nt)
        if self._boiling:
            self._y_last = yv
            snap.h_vap = (self._pkg.phase_props("vapor", T, P, list(yv)).H
                          + float(yv @ self._hf))

            def slope(v):  # dH_L/dN along v, with T on the bubble point
                scale = float(np.abs(v).max())
                if scale <= 0.0:
                    return 0.0
                eps = 1e-6 * Nt / scale
                N2 = N + eps * v
                z2 = N2 / N2.sum()
                T2, _ = self._bubble(z2, T, yv)
                return (N2.sum() * self._h_liquid(T2, z2) - snap.H) / eps

            snap.fv = (slope(gen) - self._heat(T)) / (slope(yv) - snap.h_vap)
            self._T_last = T
        return snap

    def _vent_flow(self, snap: _Snap) -> float:
        return snap.fv if self._vented and self._boiling else 0.0

    def _regime_event(self):
        """Vented batch: the subcooled liquid reaches its bubble point, or boiling stops."""
        if self._boiling:
            def g(t, y):
                return self._snap(t, y).fv

            g.direction = -1.0
        else:
            def g(t, y):
                snap = self._snap(t, y)
                return snap.T - self._bubble(snap.z, snap.T)[0]

            g.direction = 1.0
        g.terminal = True
        return g

    def _trial_rhs(self, t, y):
        """The right-hand side for the integrator. A trial point it cannot evaluate (e.g.
        amounts that a too long step drove negative) gets NaN derivatives, which makes BDF
        reject the step and retry with a shorter one."""
        try:
            return self._rhs(t, y)
        except SimulationError as exc:
            self._last_error = exc
            return np.full_like(y, np.nan)

    def _rhs(self, t, y):
        snap = self._snap(t, y)
        dy = np.zeros_like(y)
        Fv = self._vent_flow(snap)
        dy[: self._n] = snap.gen - (Fv * snap.y_vap if Fv > 0.0 else 0.0)
        if self._energy:
            Q = self._heat(snap.T)
            dy[self._iE] = Q - Fv * snap.h_vap
            dy[self._iQ] = Q
        if self._vented and Fv > 0.0:
            dy[self._iNv] = Fv * snap.y_vap
            dy[self._iHv] = Fv * snap.h_vap
        return dy

    # ---- integration -------------------------------------------------------------------

    def _events(self, y0):
        z0 = self._z0
        events = []
        if self._stop is not None:
            kind, k, target = self._stop
            if kind == "conversion":
                nk0 = self._y0k = N_REF * self._z0[k]

                def g(t, y):
                    vented = y[self._iNv][k] if self._vented else 0.0
                    return 1.0 - (max(y[k], 0.0) + vented) / nk0 - target

                g.direction = 1.0
            elif kind == "mole_fraction":
                def g(t, y):
                    N = np.maximum(y[: self._n], 0.0)
                    return N[k] / N.sum() - target

                g.direction = -1.0 if z0[k] > target else 1.0
                if z0[k] == target:
                    raise SimulationError("the charge already has the stop mole fraction",
                                          where=self.id)
            else:
                def g(t, y):
                    return self._snap(t, y).T - target

                g0 = g(0.0, y0)
                if g0 == 0.0:
                    raise SimulationError("the charge starts at the stop temperature",
                                          where=self.id)
                g.direction = -1.0 if g0 > 0.0 else 1.0
            g.terminal = True
            events.append(g)
        if self._vented:
            def empty(t, y):
                return y[: self._n].sum() - 1e-6 * N_REF

            empty.terminal, empty.direction = True, -1.0
            events.append(empty)
        return events

    def _integrate(self, y0, t_max):
        rtol = 1e-7
        events = self._events(y0)
        breaks = sorted({b for b in (self._prog[1] if self._prog else []) if 0.0 < b < t_max})
        ends = [*breaks, t_max]
        ts, ys, sols, regimes = [], [], [], []
        y, hit, t = y0, False, 0.0
        boiling = self._boiling
        for _ in range(200):  # segments: temperature-program kinks, vent regime switches
            b = next(e for e in ends if e > t)
            self._boiling = boiling
            regime = [self._regime_event()] if self._vented else []
            self._last_error = None
            sol = solve_ivp(self._trial_rhs, (t, b), y, method="BDF", rtol=rtol,
                            atol=self._atol, events=(events + regime) or None,
                            dense_output=True, max_step=t_max / 50.0)
            if not sol.success:
                why = f" ({self._last_error.message})" if self._last_error else ""
                raise SimulationError(f"the batch integration failed: {sol.message}{why}",
                                      where=self.id, hint="check the rate constants and units, "
                                      "or give a temperature")
            ts.append(sol.t)
            ys.append(sol.y)
            regimes.append(np.full(len(sol.t), boiling))
            sols.append((t, float(sol.t[-1]), sol.sol, boiling))
            y, t = sol.y[:, -1], float(sol.t[-1])
            fired = [len(e) > 0 for e in sol.t_events] if sol.status == 1 else []
            if fired and self._stop is not None and fired[0]:
                hit = True
                break
            if self._vented and fired and fired[-1]:
                boiling = not boiling  # bubble point reached, or boiling stopped
                continue
            if self._vented and fired and fired[-2]:
                raise SimulationError(
                    f"the vessel boiled dry after {t / 60:.4g} min: everything left through "
                    "the vent", where=self.id, hint="add less heat, or stop the batch earlier")
            if t >= t_max * (1.0 - 1e-12):
                break
        else:
            raise SimulationError("the vented batch switched between boiling and subcooled "
                                  "too often", where=self.id,
                                  hint="change the duty or the jacket")
        self._boiling = boiling
        self._regimes = np.concatenate(regimes)
        t_stop = float(ts[-1][-1])
        if self._stop is not None and not hit:
            kind, k, target = self._stop
            reached = self._criterion_value(t_stop, y)
            what = {"conversion": "conversion", "mole_fraction":
                    f"mole fraction of {self._names[k]}", "temperature": "temperature"}[kind]
            fmt = (lambda v: f"{v - 273.15:.2f} °C") if kind == "temperature" else (
                lambda v: f"{v:.4g}")
            raise SimulationError(
                f"the stop criterion is not reached within the batch time: the {what} is "
                f"{fmt(reached)} after {t_stop / 60:.4g} min, the target {fmt(target)}",
                where=self.id,
                hint="allow a longer batch time, or choose a target the reaction can reach "
                "(equilibrium, complete conversion)")
        self._steps = (np.concatenate(ts), np.concatenate(ys, axis=1))
        self._sols = sols
        return sols, t_stop

    def _criterion_value(self, t, y) -> float:
        kind, k, _ = self._stop
        N = np.maximum(y[: self._n], 0.0)
        if kind == "conversion":
            nk0 = self._y0k
            vented = y[self._iNv][k] if self._vented else 0.0
            return 1.0 - (N[k] + vented) / nk0
        if kind == "mole_fraction":
            return N[k] / N.sum()
        return self._snap(t, y).T

    def _dense(self, t: float) -> np.ndarray:
        """State at t from the dense output; sets the vent regime of that segment."""
        for a, b, f, boiling in self._sols:
            if t <= b + 1e-12 * max(b, 1.0):
                self._boiling = boiling
                return f(min(max(t, a), b))
        a, b, f, self._boiling = self._sols[-1]
        return f(b)

    # ---- results -----------------------------------------------------------------------

    def _heat_rate(self, t, y, snap: _Snap) -> float:
        """Instantaneous heat into the reference charge, W."""
        if self._energy:
            return self._heat(snap.T)
        # Prescribed temperature: Q = dH/dt + vent enthalpy (dU/dt in a closed vessel),
        # dH/dt by a difference along the trajectory.
        dy = self._rhs(t, y)
        eps = 1e-6 * max(self._t_stop, 1.0)
        later = self._snap(t + eps, y + eps * dy)
        if self._const_v:
            return (later.U - snap.U) / eps
        return (later.H - snap.H) / eps + self._vent_flow(snap) * snap.h_vap

    def _finish(self, pkg, s: Stream, z0, sols, t_stop, S, h_feed):
        p, n = self.p, self._n
        self._t_stop = t_stop
        F = s.total_flow
        t_cycle = S * N_REF / F
        occupied = p["charge_time"] + t_stop + p["down_time"]
        vessels = occupied / t_cycle
        if vessels > 1.0 + 1e-9:
            self.notes.append(
                f"one vessel is occupied {occupied / 3600:.3g} h per batch but the cycle time "
                f"is {t_cycle / 3600:.3g} h: {math.ceil(vessels - 1e-9)} vessels are needed "
                "in turn")
        y_end = self._steps[1][:, -1]
        end = self._snap(t_stop, y_end)

        # Profiles from the dense output, and the duty at every accepted step for its peak.
        grid = np.linspace(0.0, t_stop, N_POINTS)
        final = self._boiling
        Ys = np.column_stack([self._dense(t) for t in grid])
        regimes = []
        for t in grid:
            self._dense(t)
            regimes.append(self._boiling)
        Ys[:, 0], Ys[:, -1] = self._steps[1][:, 0], y_end
        regimes[0], regimes[-1] = bool(self._regimes[0]), final

        def at(t, y, boiling):
            self._boiling = boiling
            sn = self._snap(t, y)
            return sn, self._heat_rate(t, y, sn)

        evaluated = [at(t, Ys[:, j], regimes[j]) for j, t in enumerate(grid)]
        snaps = [e[0] for e in evaluated]
        q_grid = [e[1] for e in evaluated]
        q_all = list(q_grid)
        t_steps, y_steps = self._steps
        if len(t_steps) <= 2000:
            q_all += [at(t, y_steps[:, j], bool(self._regimes[j]))[1]
                      for j, t in enumerate(t_steps)]
        self._boiling = final
        peak = max(q_all, key=abs) if q_all else 0.0
        if self._medium is not None and any(q == 0.0 for q in q_grid[1:]):
            self.notes.append("the jacket passed no heat during part of the batch: the "
                              "contents' and the medium's temperatures crossed")
        if not self._vented and not self._const_v and any(sn.VF > 1e-9 for sn in snaps) and (
                self._phase == "liquid"):
            self.notes.append("vapor forms in the vessel and stays in it (no vent): the "
                              "vessel volume follows the contents")
        if any(sn.reacting <= 0.0 for sn in snaps):
            self.notes.append(f"the reacting phase ({self._phase}) was absent during part of "
                              "the batch, so nothing reacted then")

        # Heat per batch from the closed balance (exact, includes bringing the charge to the
        # starting temperature of a temperature program).
        Nv = y_end[self._iNv] if self._vented else np.zeros(n)
        Hv = y_end[self._iHv] if self._vented else 0.0
        if self._energy:
            heat = y_end[self._iQ]
        elif self._const_v:
            heat = end.U - self._U_charge
        else:
            heat = end.H + Hv - N_REF * h_feed

        # Outlets: the contents at the end (vented: their liquid; the vapor joins the vent).
        hf, P = self._hf, self._P
        N_end = np.maximum(y_end[:n], 0.0)
        per_s = S / t_cycle
        if self._const_v:
            out = Stream.from_flash(pkg, list(N_end * per_s), T=end.T, P=end.P)
        elif self._vented:  # the liquid as it is, also at its bubble point
            out = _liquid(pkg, N_end * per_s, end.T, P)
        else:
            h = end.H / end.Nt - float(end.z @ hf)
            out = Stream.from_flash(pkg, list(N_end * per_s), P=P, H=h)
        if Nv.sum() > 0.0:
            h_v = Hv / Nv.sum() - float(Nv @ hf) / Nv.sum()
            vent = Stream.from_flash(pkg, list(Nv * per_s), P=P, H=h_v)
        else:
            vent = empty_stream(pkg, end.T, P)

        key = self.key_index(pkg)
        conv = 1.0 - (N_end[key] + Nv[key]) / (N_REF * z0[key]) if z0[key] > 0 else None
        prod = self._product_index(pkg)
        formed = (N_end[prod] + Nv[prod] - N_REF * z0[prod]) * pkg.mws[prod] / 1000.0  # kg
        r = self.results
        r.update(
            T_out=end.T, P_final=end.P, batch_time=t_stop, cycle_time=t_cycle,
            charge=S * N_REF, vessels_needed=vessels, duty=heat * S / t_cycle,
            heat_per_batch=heat * S, peak_duty=peak * S, productivity=formed * S / t_cycle,
            vent_amount=float(Nv.sum()) * S,
            max_volume=max(sn.volume for sn in snaps) * S,
        )
        if conv is not None:
            r["conversion"] = conv
        for k, v in r.items():
            r[k] = float(v)
        if self._medium is not None:
            ut = self._medium
            UA = p["U"] * p["area"]
            q_medium = abs(y_end[self._iQ]) * S  # J per batch through the jacket
            if p["utility_flow"] is not None:
                mass = p["utility_flow"] * t_stop
            elif ut.kind == "steam":
                mass = q_medium / ut.latent
            else:
                mass = q_medium / (ut.cp * abs(ut.T_in - ut.T_out))
            r.update(UA=UA, heat_transfer_area=p["area"], utility_flow=mass / t_cycle)
        self._profiles(pkg, grid, snaps, q_grid, Ys, z0, S, key)
        return {"out": out, "vent": vent}

    def _product_index(self, pkg) -> int:
        if self.p["product"]:
            return self.component_index(pkg, self.p["product"], "product component")
        r = self._rxns[0]
        return next(i for i in (r.order or range(pkg.n)) if r.nu[i] > 0)

    def _profiles(self, pkg, grid, snaps, q_grid, Ys, z0, S, key) -> None:
        hours = grid[-1] > 2 * 3600.0
        x = [float(t / (3600.0 if hours else 60.0)) for t in grid]
        x_label = "Time [h]" if hours else "Time [min]"
        series: dict[str, tuple[str | None, list[float]]] = {}
        if key is not None and z0[key] > 0.0:
            nk0 = N_REF * z0[key]
            vented = Ys[self._iNv][key] if self._vented else np.zeros(len(grid))
            series[f"Conversion of {pkg.names[key]}"] = (
                "fraction", [float(1.0 - (max(Ys[key, j], 0.0) + vented[j]) / nk0)
                             for j in range(len(grid))])
        series["Temperature"] = ("temperature", [sn.T for sn in snaps])
        series["Pressure"] = ("pressure", [sn.P for sn in snaps])
        series["Duty"] = ("power", [float(q * S) for q in q_grid])
        if self._vented:
            series["Vented"] = ("amount", [float(Ys[self._iNv, j].sum() * S)
                                           for j in range(len(grid))])
        series["Contents volume"] = ("volume", [sn.volume * S for sn in snaps])
        self.add_profile("conditions", "Over the batch", x_label, x, series)
        present = [i for i in range(pkg.n) if max(sn.z[i] for sn in snaps) > 0.0]
        self.add_profile(
            "composition", "Contents over the batch", x_label, x,
            {f"{pkg.names[i]} (mole fraction)": ("fraction", [float(sn.z[i]) for sn in snaps])
             for i in present},
        )


def _liquid(pkg, flows: np.ndarray, T: float, P: float) -> Stream:
    """The liquid as it is in the vessel, without a re-flash across its bubble point."""
    total = float(flows.sum())
    if total <= 0.0:
        return empty_stream(pkg, T, P)
    zs = [float(f) / total for f in flows]
    return Stream([float(f) for f in flows], pkg.phase_state("liquid", T, P, zs))
