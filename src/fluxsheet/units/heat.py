"""Heating and cooling."""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from ..errors import SimulationError
from ..stream import Stream
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


class Heater(UnitOperation):
    type = "Heater"
    label = "Heater / Cooler"
    category = "Heat Exchangers"
    description = "Heats or cools a stream to a target temperature, vapor fraction, or with a "\
        "given duty."
    inlets = (IN,)
    result_fields = (DUTY,)
    outlets = (OUT,)
    params = (
        Param("T", "Outlet temperature", "temperature", min=0.0),
        Param("VF", "Outlet vapor fraction", "fraction", min=0.0, max=1.0),
        Param("duty", "Duty", "power", help="Positive = heating, negative = cooling."),
        P_OUT,
        DP,
    )
    exclusive_groups = (("T", "VF", "duty"), ("P", "dP"))
    required_groups = (("T", "VF", "duty"),)

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P = self.outlet_pressure(s.P)
        if s.total_flow <= 0.0:
            self.results["duty"] = 0.0
            return {"out": empty_stream(pkg, self.p["T"] or s.T, P)}
        if self.p["T"] is not None:
            out = Stream.from_flash(pkg, s.flows, T=self.p["T"], P=P)
        elif self.p["VF"] is not None:
            out = Stream.from_flash(pkg, s.flows, VF=self.p["VF"], P=P)
        else:
            out = Stream.from_flash(pkg, s.flows, P=P, H=s.state.H + self.p["duty"] / s.total_flow)
        self.results["duty"] = out.enthalpy_flow - s.enthalpy_flow
        return {"out": out}


# ---- two-sided exchanger -------------------------------------------------------------------

# Both temperature-duty curves are tabulated by PH flashes on a shared duty grid. Everything
# else (approach, UA, the specification solve) works on these curves, so phase changes and
# internal pinches are captured without flashing inside the root finders.
_GRID = 16  # intervals of the tables, both in duty and in temperature
_ZONES = 60  # zones of the UA integration and the approach search
_MIN_APPROACH = 1e-3  # K; closer than this needs an infinite area


def _degC(T: float) -> str:
    return f"{T - 273.15:.1f} °C"


def _kW(Q: float) -> str:
    return f"{Q / 1e3:.4g} kW"


def _log_mean(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    ratio = a / b
    with np.errstate(divide="ignore", invalid="ignore"):
        lm = (a - b) / np.log(ratio)
    return np.where(np.abs(ratio - 1.0) < 1e-6, 0.5 * (a + b), lm)


class _Curves:
    """Hot and cold temperatures as functions of the heat each side has exchanged."""

    def __init__(self, pkg, hot: Stream, cold: Stream, P_hot: float, P_cold: float,
                 countercurrent: bool):
        self.pkg, self.hot, self.cold = pkg, hot, cold
        self.P_hot, self.P_cold = P_hot, P_cold
        self.countercurrent = countercurrent
        # Neither side can pass the other's inlet temperature.
        q_hot = hot.total_flow * (hot.state.H - self._H(hot, cold.T, P_hot))
        q_cold = cold.total_flow * (self._H(cold, hot.T, P_cold) - cold.state.H)
        self.Q_cap = max(min(q_hot, q_cold), 0.0)
        self.q_hot, self.T_hot = self._table(hot, P_hot, -1.0, cold.T)
        self.q_cold, self.T_cold = self._table(cold, P_cold, 1.0, hot.T)
        self.Q_max = self._max_duty()

    def _H(self, s: Stream, T: float, P: float) -> float:
        return self.pkg.flash(s.zs, T=T, P=P).H

    def _table(self, s: Stream, P: float, sign: float, T_far: float):
        """Temperature after exchanging q (W) from 0 to Q_cap; ``sign`` -1 gives heat away.

        Points are spaced evenly in duty and in temperature, so steep stretches (a vapor
        with its low heat capacity) are resolved as well as flat ones. Bubble and dew points
        are points too: the kinks of the curve, and the flat stretch of a pure component,
        are exact rather than cut off by interpolation.
        """
        def heat(H: float) -> float:
            return sign * s.total_flow * (H - s.state.H)

        points = {}
        for dq in np.linspace(0.0, self.Q_cap, _GRID + 1):
            points[dq] = self.pkg.flash(s.zs, P=P, H=s.state.H + sign * dq / s.total_flow).T
        for T in np.linspace(s.T, T_far, _GRID + 1)[1:-1]:
            points[heat(self._H(s, T, P))] = T
        for vf in (0.0, 1.0):
            try:
                state = self.pkg.flash(s.zs, VF=vf, P=P)
            except SimulationError:  # e.g. no bubble point for a mixture rich in hydrogen
                continue
            points[heat(state.H)] = state.T
        q = np.array(sorted(dq for dq in points if 0.0 <= dq <= self.Q_cap))
        return q, np.array([points[dq] for dq in q])

    def profile(self, Q: float, n: int = _ZONES) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Positions (heat given up by the hot side since its inlet) and both temperatures.

        Besides ``n`` even steps, the positions include every grid point of both tables; the
        temperature difference is linear in between, so its minimum is found exactly.
        """
        cold_at = self.q_cold if not self.countercurrent else Q - self.q_cold
        breaks = np.concatenate([self.q_hot, cold_at])
        x = np.unique(np.concatenate([np.linspace(0.0, Q, n + 1),
                                      breaks[(breaks > 0.0) & (breaks < Q)]]))
        T_hot = np.interp(x, self.q_hot, self.T_hot)
        T_cold = np.interp(Q - x if self.countercurrent else x, self.q_cold, self.T_cold)
        return x, T_hot, T_cold

    def approach(self, Q: float) -> float:
        _, T_hot, T_cold = self.profile(Q)
        return float(np.min(T_hot - T_cold))

    def UA(self, Q: float) -> float:
        if Q <= 0.0:
            return 0.0
        x, T_hot, T_cold = self.profile(Q)
        dT = T_hot - T_cold
        if np.min(dT) <= 0.0:
            return math.inf
        return float(np.sum(np.diff(x) / _log_mean(dT[:-1], dT[1:])))

    def _max_duty(self) -> float:
        """Largest duty with a non-negative approach everywhere (the approach falls with Q)."""
        if self.Q_cap <= 0.0 or self.approach(self.Q_cap) >= 0.0:
            return self.Q_cap
        return brentq(self.approach, 0.0, self.Q_cap, xtol=1e-9 * self.Q_cap)


class HeatExchanger(UnitOperation):
    type = "HeatExchanger"
    label = "Heat Exchanger"
    category = "Heat Exchangers"
    description = "Passes heat from a hot to a cold stream, e.g. a feed/effluent exchanger. "\
        "Specify one outlet temperature, the duty, the minimum approach temperature, or UA "\
        "(rating). Phase changes are followed along the exchanger by zone analysis."
    inlets = (Port("hot_in", "Hot inlet"), Port("cold_in", "Cold inlet"))
    outlets = (Port("hot_out", "Hot outlet"), Port("cold_out", "Cold outlet"))
    params = (
        Param("T_hot_out", "Hot outlet T", "temperature", min=0.0),
        Param("T_cold_out", "Cold outlet T", "temperature", min=0.0),
        Param("duty", "Duty", "power", min=0.0, help="Heat passed from the hot to the cold side."),
        Param("min_approach", "Min. approach", "temperature_difference", min=0.0,
              help="Smallest hot-minus-cold temperature difference anywhere in the exchanger."),
        Param("UA", "UA", "thermal_conductance", min=0.0,
              help="Overall heat transfer coefficient times area (rating mode)."),
        Param("arrangement", "Flow arrangement", "choice", default="countercurrent",
              options=("countercurrent", "cocurrent")),
        Param("dP_hot", "Hot side pressure drop", "pressure_difference", min=0.0),
        Param("dP_cold", "Cold side pressure drop", "pressure_difference", min=0.0),
    )
    exclusive_groups = (("T_hot_out", "T_cold_out", "duty", "min_approach", "UA"),)
    required_groups = exclusive_groups
    result_fields = (
        DUTY,
        ResultField("LMTD", "Mean temperature difference (LMTD)", "temperature_difference"),
        ResultField("UA", "UA", "thermal_conductance"),
        ResultField("min_approach", "Minimum approach", "temperature_difference"),
    )

    def calculate(self, pkg, inlets):
        hot, cold = inlets["hot_in"][0], inlets["cold_in"][0]
        P_hot = hot.P - (self.p["dP_hot"] or 0.0)
        P_cold = cold.P - (self.p["dP_cold"] or 0.0)
        for side, P in (("hot", P_hot), ("cold", P_cold)):
            if P <= 0.0:
                raise SimulationError(f"the {side} side pressure drop exceeds its inlet pressure",
                                      hint=f"lower dP_{side}")

        if hot.total_flow <= 0.0 or cold.total_flow <= 0.0:
            self.results.update(duty=0.0, LMTD=None, UA=0.0, min_approach=hot.T - cold.T)
            return {"hot_out": Stream.from_flash(pkg, hot.flows, T=hot.T, P=P_hot),
                    "cold_out": Stream.from_flash(pkg, cold.flows, T=cold.T, P=P_cold)}
        if hot.T <= cold.T:
            raise SimulationError(
                f"the hot inlet ({_degC(hot.T)}) is not hotter than the cold inlet "
                f"({_degC(cold.T)})",
                hint="swap the streams connected to hot_in and cold_in",
            )

        countercurrent = self.p["arrangement"] == "countercurrent"
        curves = _Curves(pkg, hot, cold, P_hot, P_cold, countercurrent)
        Q = self._duty(curves)

        hot_out = Stream.from_flash(pkg, hot.flows, P=P_hot, H=hot.state.H - Q / hot.total_flow)
        cold_out = Stream.from_flash(pkg, cold.flows, P=P_cold,
                                     H=cold.state.H + Q / cold.total_flow)
        UA = curves.UA(Q)
        self.results.update(
            duty=Q, UA=UA, min_approach=curves.approach(Q),
            LMTD=Q / UA if Q > 0.0 else hot.T - cold.T,
        )
        x, T_hot, T_cold = curves.profile(Q, n=20)
        self.add_profile(
            "temperatures", "Temperature along the exchanger", "Heat transferred [kW]",
            [float(v) / 1e3 for v in x],
            {"Hot side": ("temperature", T_hot.tolist()),
             "Cold side": ("temperature", T_cold.tolist())},
        )
        return {"hot_out": hot_out, "cold_out": cold_out}

    def _duty(self, c: _Curves) -> float:
        hot, cold, p = c.hot, c.cold, self.p
        if p["duty"] is not None:
            Q = p["duty"]
        elif p["T_hot_out"] is not None:
            if p["T_hot_out"] >= hot.T:
                raise SimulationError(
                    f"hot outlet T {_degC(p['T_hot_out'])} is not below the hot inlet "
                    f"({_degC(hot.T)})",
                    hint="the hot side gives heat away; set a lower hot outlet temperature",
                )
            Q = hot.total_flow * (hot.state.H - c._H(hot, p["T_hot_out"], c.P_hot))
        elif p["T_cold_out"] is not None:
            if p["T_cold_out"] <= cold.T:
                raise SimulationError(
                    f"cold outlet T {_degC(p['T_cold_out'])} is not above the cold inlet "
                    f"({_degC(cold.T)})",
                    hint="the cold side takes heat up; set a higher cold outlet temperature",
                )
            Q = cold.total_flow * (c._H(cold, p["T_cold_out"], c.P_cold) - cold.state.H)
        elif p["min_approach"] is not None:
            inlet_gap = c.approach(0.0)
            if p["min_approach"] >= inlet_gap:
                raise SimulationError(
                    f"the inlets are only {inlet_gap:.3g} K apart, less than the minimum "
                    f"approach of {p['min_approach']:.3g} K",
                    hint=f"set a minimum approach below {inlet_gap:.3g} K, or specify the duty",
                )
            if c.approach(c.Q_max) >= p["min_approach"]:
                return c.Q_max
            return brentq(lambda Q: c.approach(Q) - p["min_approach"], 0.0, c.Q_max,
                          xtol=1e-9 * max(c.Q_max, 1.0))
        else:
            if p["UA"] == 0.0:
                return 0.0
            # UA grows without bound as the approach closes, so the root lies below Q_max.
            Q_hi = c.Q_max * (1.0 - 1e-9)
            if c.UA(Q_hi) <= p["UA"]:
                return Q_hi
            return brentq(lambda Q: c.UA(Q) - p["UA"], 0.0, Q_hi, xtol=1e-9 * max(c.Q_max, 1.0))
        self._check_feasible(c, Q)
        return Q

    def _check_feasible(self, c: _Curves, Q: float) -> None:
        if Q <= c.Q_max and (Q == 0.0 or c.approach(Q) >= _MIN_APPROACH):
            return
        Q_ok = c.Q_max
        _, T_hot, T_cold = c.profile(Q_ok)
        arrangement = self.p["arrangement"]
        limits = {
            "T_hot_out": f"the hot outlet cannot be cooled below {_degC(T_hot[-1])}",
            "T_cold_out": "the cold outlet cannot be heated above "
                          f"{_degC(T_cold[0] if c.countercurrent else T_cold[-1])}",
            "duty": f"the duty cannot exceed {_kW(Q_ok)}",
        }
        given = next(k for k in limits if self.p[k] is not None)
        hint = f"{limits[given]} in a {arrangement} exchanger (zero approach)"
        if not c.countercurrent:
            hint += "; countercurrent flow allows more heat recovery"
        raise SimulationError(
            f"temperature cross: passing {_kW(Q)} would make the hot side colder than the cold "
            f"side (at most {_kW(Q_ok)} can be transferred)",
            hint=f"{hint}; or specify a minimum approach temperature, e.g. 10 K",
        )
