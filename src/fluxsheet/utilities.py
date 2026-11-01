"""Heating and cooling utilities, and heat transfer through a jacket or coil.

Utilities are defined once per flowsheet and chosen by id on any unit with an energy
balance. Two kinds:

- ``fluid``: a heat transfer fluid (cooling water, chilled water, thermal oil) entering at
  ``T_in`` and leaving at ``T_out`` with a constant heat capacity ``cp``. The flow follows
  from the duty. A unit may fix the flow instead (``utility_flow``); then the outlet
  temperature follows.
- ``steam``: saturated steam at pressure ``P`` that condenses at constant temperature. The
  saturation temperature and the latent heat come from IAPWS-95 (``chemicals.iapws``), the
  latent heat by Clausius-Clapeyron, h_fg = T (v_g - v_l) dPsat/dT.

With a utility, a unit either reports the utility consumption for its duty, or, given U
and the area of a jacket or coil, finds the duty itself: Q = UA * LMTD, solved together with
the unit's own energy balance. On the process side a unit is either well mixed (vessels,
drums, stirred tanks: the contents are at the outlet temperature everywhere) or in plug flow
(a heater: from the inlet to the outlet temperature, counter- or cocurrent to the utility).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from chemicals.iapws import iapws95_dPsat_dT, iapws95_rhog_sat, iapws95_rhol_sat, iapws95_Tsat
from scipy.optimize import brentq

from .errors import SimulationError
from .quantities import UnitError, to_si

KINDS = ("fluid", "steam")
P_CRITICAL_WATER = 22.064e6


@dataclass(frozen=True)
class Exchange:
    """Heat exchanged with a utility. ``Q`` is the heat into the process, W."""

    Q: float
    flow: float  # utility mass flow, kg/s
    T_out: float  # utility outlet temperature, K
    LMTD: float | None  # mean temperature difference, K; None while nothing is exchanged


def _degC(T: float) -> str:
    return f"{T - 273.15:.1f} °C"


def lmtd(d1: float, d2: float) -> float | None:
    """Log mean of two end temperature differences of one sign; None if they cross or touch."""
    if d1 * d2 <= 0.0:
        return None
    if abs(d1 - d2) <= 1e-9 * abs(d1):
        return d1
    return (d1 - d2) / math.log(d1 / d2)


def _ends(T_pin, T_pout, Tu_in, Tu_out, arrangement) -> tuple[float, float]:
    """Utility minus process temperature at both ends (positive: heat flows into the process)."""
    if arrangement == "cocurrent":
        return Tu_in - T_pin, Tu_out - T_pout
    return Tu_out - T_pin, Tu_in - T_pout


@dataclass(frozen=True)
class Utility:
    id: str
    kind: str
    T_in: float  # K; for steam the saturation temperature
    T_out: float  # K
    cp: float = 0.0  # J/(kg K), fluids
    latent: float = 0.0  # J/kg, steam
    P: float | None = None  # Pa, steam

    @property
    def heats(self) -> bool:
        return self.kind == "steam" or self.T_in > self.T_out

    def describe(self) -> str:
        if self.kind == "steam":
            return f"steam {self.id} ({_degC(self.T_in)})"
        return f"{self.id} ({_degC(self.T_in)} → {_degC(self.T_out)})"

    def _outlet(self, Q: float, flow: float | None) -> tuple[float, float]:
        """Utility flow and outlet temperature for heat ``Q`` into the process."""
        if self.kind == "steam":
            return Q / self.latent, self.T_in
        if flow is None:
            return Q / (self.cp * (self.T_in - self.T_out)), self.T_out
        return flow, self.T_in - Q / (flow * self.cp) if flow > 0.0 else self.T_in

    def for_duty(self, Q: float, T_pin: float, T_pout: float, arrangement: str,
                 flow: float | None = None) -> Exchange:
        """Consumption and mean temperature difference for a known duty."""
        if abs(Q) <= 1e-9:
            return Exchange(0.0, 0.0, self.T_in, None)
        if (Q > 0.0) != self.heats and (self.kind == "steam" or flow is None):
            raise SimulationError(
                f"{self.describe()} can only {'heat' if self.heats else 'cool'}, but the unit "
                f"needs {'heating' if Q > 0 else 'cooling'} ({Q / 1e3:.4g} kW)",
                hint=f"choose a {'cooling' if Q < 0 else 'heating'} utility",
            )
        if flow is not None and flow <= 0.0:
            raise SimulationError("the utility flow must be positive",
                                  hint="give a positive medium flow, or leave it empty")
        m, T_out = self._outlet(Q, flow)
        d = _ends(T_pin, T_pout, self.T_in, T_out, arrangement)
        mean = lmtd(*d)
        if mean is None or (mean > 0.0) != (Q > 0.0):
            process = (_degC(T_pout) if T_pin == T_pout
                       else f"{_degC(T_pin)} → {_degC(T_pout)}")
            raise SimulationError(
                f"temperature cross: {self.describe()}"
                + (f", leaving at {_degC(T_out)}," if T_out != self.T_out else "")
                + f" cannot {'heat' if Q > 0 else 'cool'} the process at {process}",
                hint="choose a utility with more temperature difference, change its outlet "
                "temperature or flow, or the unit's temperature",
            )
        return Exchange(Q, m, T_out, abs(mean))

    def rate(self, UA: float, T_pin: float, T_pout: float, arrangement: str,
             flow: float | None = None) -> float:
        """Heat into the process through UA at the given process temperatures; 0 if crossed."""
        if self.kind == "steam" or flow is None:
            mean = lmtd(*_ends(T_pin, T_pout, self.T_in, self.T_out, arrangement))
            if mean is None or (mean > 0.0) != self.heats:
                return 0.0
            return UA * mean
        if flow <= 0.0 or UA <= 0.0:
            return 0.0
        # Fixed flow: the outlet temperature T balances C (T_in - T) = UA * LMTD. It lies
        # between the inlet and the process temperature it can approach at most.
        C = flow * self.cp
        limit = T_pin if arrangement != "cocurrent" else T_pout

        def balance(T: float) -> float:
            mean = lmtd(*_ends(T_pin, T_pout, self.T_in, T, arrangement)) or 0.0
            return C * (self.T_in - T) - UA * mean

        span = limit - self.T_in
        if abs(span) < 1e-9:
            return 0.0
        a, b = self.T_in + 1e-9 * span, limit - 1e-9 * span
        if balance(a) * balance(b) > 0.0:
            return 0.0
        T = brentq(balance, min(a, b), max(a, b), xtol=1e-10, rtol=1e-12)
        return C * (self.T_in - T)


# ---- flowsheet definitions ------------------------------------------------------------------

def parse_utilities(models) -> dict[str, Utility]:
    """Utilities by id from the flowsheet's ``utilities`` list."""
    out: dict[str, Utility] = {}
    for m in models or ():
        where = f"utility {m.id}"
        if m.id in out:
            raise SimulationError(f"utility {m.id!r} is defined twice", where=where)
        try:
            out[m.id] = _parse(m)
        except UnitError as exc:
            raise SimulationError(str(exc), where=where) from None
        except SimulationError as exc:
            raise exc.located(where)
    return out


def _parse(m) -> Utility:
    if m.type == "steam":
        if m.P is None:
            raise SimulationError("steam needs its pressure", hint="e.g. \"10 bar\"")
        P = to_si(m.P, "pressure")
        if not 612.0 < P < P_CRITICAL_WATER:
            raise SimulationError(f"steam pressure {P / 1e5:.4g} bar is outside the "
                                  "saturation range of water",
                                  hint="use a pressure between 0.01 and 220 bar")
        T = iapws95_Tsat(P)
        dPdT = iapws95_dPsat_dT(T)[0]
        latent = T * (1.0 / iapws95_rhog_sat(T) - 1.0 / iapws95_rhol_sat(T)) * dPdT
        return Utility(m.id, "steam", T, T, latent=latent, P=P)
    if m.type != "fluid":
        raise SimulationError(f"unknown utility type {m.type!r}",
                              hint=f"use one of: {', '.join(KINDS)}")
    if m.T_in is None or m.T_out is None:
        raise SimulationError("a utility fluid needs its inlet and outlet temperatures")
    T_in, T_out = to_si(m.T_in, "temperature"), to_si(m.T_out, "temperature")
    if abs(T_in - T_out) < 1e-6:
        raise SimulationError("inlet and outlet temperature of the utility are the same",
                              hint="a fluid utility warms up or cools down; give two "
                              "different temperatures")
    cp = to_si(m.cp, "specific_heat") if m.cp is not None else 4180.0
    if cp <= 0.0:
        raise SimulationError("the heat capacity must be positive")
    return Utility(m.id, "fluid", T_in, T_out, cp=cp)


# ---- unit side --------------------------------------------------------------------------

def attach(unit, utilities: dict[str, Utility]) -> None:
    """Resolve the unit's utility by id and check that its heat transfer spec is complete."""
    name = unit.p.get("utility")
    unit.utility = None
    if name is None:
        if unit.p.get("area") is not None:
            raise SimulationError("a jacket or coil needs a heating or cooling medium",
                                  where=unit.id, hint="choose a utility for the unit")
        return
    if name not in utilities:
        known = ", ".join(utilities) or "none defined"
        raise SimulationError(f"unknown utility {name!r}", where=unit.id,
                              hint=f"define it in the flowsheet's utilities ({known})")
    unit.utility = utilities[name]
    if unit.p.get("area") is not None and not unit.p.get("U"):
        raise SimulationError("the heat transfer area is given but not U", where=unit.id,
                              hint="set the heat transfer coefficient U")


def calculate(unit, pkg, inlets):
    """Run ``unit`` with its utility: rate a jacket or coil, or report the consumption."""
    ut = getattr(unit, "utility", None)
    if ut is None:
        return unit.calculate(pkg, inlets)
    plug = unit.process_side == "plug"
    feeds = [s for streams in inlets.values() for s in streams]
    T_in = feeds[0].T
    arrangement = unit.p.get("arrangement") or "countercurrent"
    flow = unit.p.get("utility_flow")

    def run(Q: float):
        unit.p["duty"] = Q
        outputs = unit.calculate(pkg, inlets)
        return outputs, _outlet_T(unit, outputs)

    area = unit.p.get("area")
    if area is not None:
        UA = unit.p["U"] * area
        try:
            outputs, T_out = run(0.0)
            q0 = ut.rate(UA, T_in if plug else T_out, T_out, arrangement, flow)
            if q0 != 0.0:
                def residual(Q: float) -> float:
                    T = run(Q)[1]
                    return Q - ut.rate(UA, T_in if plug else T, T, arrangement, flow)

                lo, hi = _bracket(residual, q0)
                Q = brentq(residual, min(lo, hi), max(lo, hi),
                           xtol=1e-9 * abs(q0) + 1e-6, rtol=1e-12)
                outputs, T_out = run(Q)
        finally:
            unit.p["duty"] = None
        Q = unit.results["duty"]
        if Q == 0.0 and flow is None:
            # Nothing flows through the wall: explain it if the temperatures cross.
            try:
                ut.for_duty(1.0 if ut.heats else -1.0, T_in if plug else T_out, T_out, arrangement)
            except SimulationError as exc:
                raise exc.located(unit.id)
    else:
        outputs = unit.calculate(pkg, inlets)
        Q, T_out = unit.results["duty"], _outlet_T(unit, outputs)
    try:
        ex = ut.for_duty(Q, T_in if plug else T_out, T_out, arrangement, flow)
    except SimulationError as exc:
        raise exc.located(unit.id)
    if area is not None:
        UA = unit.p["U"] * area
    else:  # the UA (and with U the area) this duty needs
        UA = abs(Q) / ex.LMTD if ex.LMTD else None
        area = UA / unit.p["U"] if UA and unit.p.get("U") else None
    unit.results.update(utility_flow=abs(ex.flow), utility_T_out=ex.T_out, LMTD=ex.LMTD,
                        UA=UA, heat_transfer_area=area)
    return outputs


def _bracket(residual, q0: float) -> tuple[float, float]:
    """Duties around the root of ``residual`` between 0 and ``q0``.

    The residual (duty minus what the jacket passes at the resulting temperature) rises
    monotonically from -q0 at zero duty and has turned at ``q0``, the heat the jacket passes
    at the adiabatic temperature. That end can be far off (a reactor cooled below its dew
    point at that duty), so where the unit cannot be calculated the step is halved.
    """
    lo, hi, bad = 0.0, q0, None
    for _ in range(60):
        try:
            f = residual(hi)
        except SimulationError:
            bad, hi = hi, lo + (hi - lo) / 2.0
            continue
        if f == 0.0 or (f > 0.0) == (q0 > 0.0):
            return lo, hi
        if bad is None:
            break
        lo, hi = hi, hi + (bad - hi) / 2.0
    raise SimulationError("the jacket or coil duty could not be found",
                          hint="check U, the area and the utility, or start from a "
                          "given temperature or duty")


def _outlet_T(unit, outputs) -> float:
    if unit.results.get("T_out") is not None:
        return unit.results["T_out"]
    streams = list(outputs.values())
    flowing = [s for s in streams if s.total_flow > 0.0]
    return (flowing or streams)[0].T
