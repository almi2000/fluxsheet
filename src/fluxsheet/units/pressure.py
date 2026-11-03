"""Pressure changers: valves, pumps, compressors.

Pumps and compressors take either the outlet pressure (or the rise) or a performance curve:
head and, optionally, efficiency against the actual inlet volumetric flow at the rated speed
(see fluxsheet.curves). At another speed N the affinity (fan) laws scale the curve with
r = N / N_rated: flow by r, head by r^2, the efficiency stays that of the corresponding point.
"""

from __future__ import annotations

import math

from scipy.optimize import brentq

from ..curves import Curve, G
from ..errors import SimulationError
from ..stream import Stream
from ..thermo_package import PropertyPackage
from .base import (
    IN,
    OUT,
    P_OUT,
    POWER,
    Param,
    ResultField,
    UnitOperation,
    empty_stream,
    single_inlet,
)
from .pressure_machines import isentropic_work

POLYTROPIC_STEPS = 12  # pressure steps of the polytropic path over the expected ratio


class Valve(UnitOperation):
    type = "Valve"
    label = "Valve"
    category = "Pressure Changers"
    description = "Isenthalpic pressure reduction (throttling). The outlet may partially "\
        "vaporize and cool."
    inlets = (IN,)
    outlets = (OUT,)
    params = (P_OUT, Param("dP", "Pressure drop", "pressure_difference", min=0.0))
    exclusive_groups = (("P", "dP"),)
    required_groups = exclusive_groups

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P = self.outlet_pressure(s.P)
        if P > s.P:
            raise SimulationError(
                f"outlet pressure ({P / 1e5:.3g} bar) is above inlet pressure "
                f"({s.P / 1e5:.3g} bar)", where=self.id,
                hint="a valve can only reduce pressure; use a pump or compressor to raise it",
            )
        if s.total_flow <= 0.0:
            return {"out": empty_stream(pkg, s.T, P)}
        return {"out": Stream.from_flash(pkg, s.flows, P=P, H=s.state.H)}


def _curve_params(head_help: str, head_unit: str) -> tuple[Param, ...]:
    return (
        Param("curve", "Performance curve", "curve", tab="Curve", options=(head_unit,),
              help=f"Head{head_help} and optionally efficiency against the actual inlet "
              "volumetric flow, at the rated speed."),
        Param("speed", "Speed", min=1e-9, requires="curve",
              help="Shaft speed, e.g. in rpm. Scales the curve by the affinity laws; needs "
              "the rated speed in the same unit."),
        Param("rated_speed", "Rated speed", min=1e-9, requires="curve",
              help="The speed the curve was measured at."),
    )


_FLOW = ResultField("volumetric_flow", "Inlet volumetric flow", "volumetric_flow")
_ETA = ResultField("efficiency", "Efficiency", "fraction")
_P = ResultField("P", "Outlet pressure", "pressure")


class _Machine(UnitOperation):
    result_fields = (POWER, _P)
    inlets = (IN,)
    outlets = (OUT,)
    exclusive_groups = (("P", "dP", "curve"),)
    required_groups = exclusive_groups
    raises_pressure = True

    def __post_init__(self) -> None:
        super().__post_init__()
        speed, rated = self.p["speed"], self.p["rated_speed"]
        if (speed is None) != (rated is None):
            raise SimulationError(
                "the speed and the rated speed go together", where=self.id,
                hint="give both (in the same unit), or neither to run at the rated speed")
        if speed is not None and self.p["curve"] is None:
            raise SimulationError("a speed needs a performance curve", where=self.id,
                                  hint="specify the machine by its curve, or remove the speed")

    def speed_ratio(self) -> float:
        return self.p["speed"] / self.p["rated_speed"] if self.p["speed"] is not None else 1.0

    def on_curve(self, s: Stream, Q: float, what: str) -> tuple[float, float]:
        """Head (J/kg) and efficiency at actual flow ``Q`` from the curve and the speed."""
        curve: Curve = self.p["curve"]
        r = self.speed_ratio()
        q = Q / r
        eta = curve.efficiency(q)
        if eta is None:
            eta = self.p["efficiency"]
        head = curve.head(q) * r * r
        if Q > 0.0 and head <= 0.0:
            self.notes.append(f"the flow is beyond the end of the {what} curve: no head left")
        return head, eta

    def range_note(self, Q: float, what: str) -> None:
        curve: Curve = self.p["curve"]
        q, r = Q / self.speed_ratio(), self.speed_ratio()
        if Q > 0.0 and not curve.inside(q):
            self.notes.append(
                f"the flow ({Q * 3600:.4g} m3/h) is outside the {what} curve "
                f"({curve.flows[0] * r * 3600:.4g} to {curve.flows[-1] * r * 3600:.4g} m3/h "
                "at this speed); the head is extrapolated")

    def calculate(self, pkg: PropertyPackage, inlets):
        s = single_inlet(inlets)
        if self.p["curve"] is not None:
            return self.calculate_curve(pkg, s)
        P = self.p["P"] if self.p["P"] is not None else s.P + self.p["dP"]
        if P < s.P:
            raise SimulationError(
                f"outlet pressure ({P / 1e5:.3g} bar) is below inlet pressure "
                f"({s.P / 1e5:.3g} bar)", where=self.id,
                hint="use a valve to reduce pressure",
            )
        self.results["P"] = P
        if s.total_flow <= 0.0:
            self.results["power"] = 0.0
            return {"out": empty_stream(pkg, s.T, P)}
        # Liquid: incompressible isentropic work, V dP (see isentropic_work).
        work = isentropic_work(pkg, s, P) / self.p["efficiency"]
        out = Stream.from_flash(pkg, s.flows, P=P, H=s.state.H + work)
        self.results["power"] = work * s.total_flow
        return {"out": out}

    def calculate_curve(self, pkg: PropertyPackage, s: Stream) -> dict[str, Stream]:
        raise NotImplementedError


def _molar_mass(pkg: PropertyPackage, s: Stream) -> float:
    """kg/mol of the mixture."""
    return sum(z * mw for z, mw in zip(s.zs, pkg.mws)) / 1000.0


class Pump(_Machine):
    type = "Pump"
    label = "Pump"
    category = "Pressure Changers"
    description = "Raises the pressure of a liquid, to a given pressure or along its head "\
        "curve (with the affinity laws at another speed). Power from an isentropic path and "\
        "efficiency. Reports the net positive suction head available."
    params = (
        P_OUT,
        Param("dP", "Pressure increase", "pressure_difference", min=0.0),
        *_curve_params(" (a liquid column, or energy per mass)", "m"),
        Param("efficiency", "Efficiency", "fraction", default=0.75, min=0.01, max=1.0,
              help="Used unless the curve has an efficiency column."),
        Param("NPSH_required", "NPSH required", "length", min=0.0,
              help="Warn when the net positive suction head available is lower."),
        Param("suction_height", "Liquid level above pump", "length",
              help="Static head on the suction side, added to the NPSH available (negative "
              "for a suction lift). Defaults to zero."),
    )
    result_fields = (
        POWER, _P,
        ResultField("head", "Head", "length"),
        _FLOW, _ETA,
        ResultField("NPSH_available", "NPSH available", "length"),
    )

    def __post_init__(self) -> None:
        super().__post_init__()
        self._vapor_pressure: dict = {}

    def calculate(self, pkg, inlets):
        for key in ("head", "volumetric_flow", "efficiency", "NPSH_available"):
            self.results.pop(key, None)  # from an earlier iteration
        out = super().calculate(pkg, inlets)
        s = single_inlet(inlets)
        if s.total_flow <= 0.0:
            return out
        try:  # reported extras: a pump specified by its pressure never fails on them
            V, M = self.liquid(pkg, s)
        except Exception:  # noqa: BLE001
            return out
        rho = M / V
        self.results["volumetric_flow"] = self.liquid_flow(s) * V
        self.results["head"] = (self.results["P"] - s.P) / (rho * G)
        self.npsh(pkg, s, rho)
        return out

    @staticmethod
    def liquid(pkg: PropertyPackage, s: Stream) -> tuple[float, float]:
        """Molar volume (m3/mol) and molar mass (kg/mol) of the liquid the pump moves: of
        the liquid part if vapor is present."""
        if 1e-9 < s.state.VF < 1.0 and s.state.x:
            V = pkg.phase_props("liquid", s.T, s.P, list(s.state.x)).V
            return V, sum(x * mw for x, mw in zip(s.state.x, pkg.mws)) / 1000.0
        return s.state.V, _molar_mass(pkg, s)

    @staticmethod
    def liquid_flow(s: Stream) -> float:
        """mol/s of liquid; all of the flow unless the inlet is two-phase."""
        return s.total_flow * (1.0 - s.state.VF) if 1e-9 < s.state.VF < 1.0 else s.total_flow

    def npsh(self, pkg: PropertyPackage, s: Stream, rho: float) -> None:
        if s.state.VF > 1e-6:
            return  # vapor at the inlet: the solver already warns about cavitation
        key = (tuple(round(z, 10) for z in s.zs), round(s.T, 6))
        P_vap = self._vapor_pressure.get(key)
        if P_vap is None:
            try:
                P_vap = pkg.saturation_point(s.zs, "bubble", T=s.T).P
            except Exception:  # noqa: BLE001 - e.g. no bubble point for a supercritical gas
                if self.p["NPSH_required"] is not None:
                    self.notes.append("no vapor pressure of the inlet liquid; NPSH available "
                                      "not calculated")
                return
            if len(self._vapor_pressure) > 32:
                self._vapor_pressure.clear()
            self._vapor_pressure[key] = P_vap
        npsh = (s.P - P_vap) / (rho * G) + (self.p["suction_height"] or 0.0)
        self.results["NPSH_available"] = npsh
        required = self.p["NPSH_required"]
        if required is not None and npsh < required:
            self.notes.append(
                f"NPSH available ({npsh:.3g} m) is below NPSH required ({required:.3g} m): the "
                "pump will cavitate (the vapor pressure is the inlet's bubble point, which "
                "dissolved light gases raise)")

    def calculate_curve(self, pkg, s):
        V, M = self.liquid(pkg, s)
        flow = self.liquid_flow(s)
        head, eta = self.on_curve(s, flow * V, "pump")
        self.range_note(flow * V, "pump")
        P = s.P + M / V * head
        self.results.update(P=P, efficiency=eta)
        if s.total_flow <= 0.0:
            self.results["power"] = 0.0
            return {"out": empty_stream(pkg, s.T, P)}
        # V dP / eta per mole of liquid, with dP = rho head; vapor (if any) just rides along.
        power = flow * M * head / eta
        self.results["power"] = power
        return {"out": Stream.from_flash(pkg, s.flows, P=P,
                                         H=s.state.H + power / s.total_flow)}


class Compressor(_Machine):
    type = "Compressor"
    label = "Compressor"
    category = "Pressure Changers"
    description = "Compresses a gas, to a given pressure or along its head curve (polytropic "\
        "or isentropic, with the fan laws at another speed). The efficiency sets power and "\
        "outlet temperature; warns outside the curve (surge, stonewall)."
    params = (
        P_OUT,
        Param("dP", "Pressure increase", "pressure_difference", min=0.0),
        *_curve_params(" (energy per mass, or as a length)", "kJ/kg"),
        Param("curve_type", "Head and efficiency", "choice", default="polytropic",
              options=("polytropic", "isentropic"), requires="curve",
              help="Whether the curve gives polytropic or isentropic head and efficiency."),
        Param("surge_flow", "Surge flow", "volumetric_flow", min=0.0, requires="curve",
              help="Lowest stable inlet flow at the rated speed. Defaults to the first curve "
              "point."),
        Param("efficiency", "Efficiency", "fraction", default=0.72, min=0.01, max=1.0,
              help="Isentropic efficiency; with a curve without an efficiency column, of the "
              "curve's kind (polytropic or isentropic)."),
    )
    result_fields = (
        POWER, _P,
        ResultField("head", "Head", "specific_energy"),
        _FLOW, _ETA,
        ResultField("pressure_ratio", "Pressure ratio"),
    )

    def calculate_curve(self, pkg, s):
        if s.total_flow <= 0.0:  # nothing to compress (e.g. a tear stream's first guess)
            self.results.update(P=s.P, power=0.0, pressure_ratio=1.0)
            return {"out": empty_stream(pkg, s.T, s.P)}
        if s.state.VF < 1.0 - 1e-6:
            raise SimulationError("the compressor inlet contains liquid", where=self.id,
                                  hint="separate the liquid first, e.g. with a flash drum")
        mw = _molar_mass(pkg, s)
        Q = s.total_flow * s.state.V
        head, eta = self.on_curve(s, Q, "compressor")
        curve: Curve = self.p["curve"]
        r = self.speed_ratio()
        q = Q / r
        surge = self.p["surge_flow"] if self.p["surge_flow"] is not None else curve.flows[0]
        if Q > 0.0 and q < surge:
            self.notes.append(
                f"the inlet flow ({Q * 3600:.4g} m3/h) is below the surge limit "
                f"({surge * r * 3600:.4g} m3/h at this speed): the compressor surges")
        elif Q > 0.0 and q > curve.flows[-1]:
            self.notes.append(
                f"the inlet flow ({Q * 3600:.4g} m3/h) is beyond the end of the curve "
                f"({curve.flows[-1] * r * 3600:.4g} m3/h at this speed): stonewall (choke)")
        P = self.discharge_pressure(pkg, s, head, eta, mw)
        work = mw * head / eta  # J/mol
        self.results.update(P=P, head=head, volumetric_flow=Q, efficiency=eta,
                            pressure_ratio=P / s.P, power=work * s.total_flow)
        return {"out": Stream.from_flash(pkg, s.flows, P=P, H=s.state.H + work)}

    def discharge_pressure(self, pkg: PropertyPackage, s: Stream, head: float, eta: float,
                           mw: float) -> float:
        """Outlet pressure at which the isentropic or polytropic head equals ``head``."""
        if head <= 0.0:
            return s.P
        # Ideal-gas estimate of the pressure ratio, to scale the search.
        k = 1.4
        exponent = (k - 1.0) / k / (eta if self.p["curve_type"] == "polytropic" else 1.0)
        RT = s.P * s.state.V  # Z R T per mole
        ratio = (1.0 + head * mw * exponent / RT) ** (1.0 / exponent)
        if self.p["curve_type"] == "isentropic":
            def excess(lnP: float) -> float:
                return isentropic_work(pkg, s, math.exp(lnP)) / mw - head
            lo, hi = math.log(s.P), math.log(s.P * ratio)
            while excess(hi) < 0.0:
                lo, hi = hi, hi + 2.0 * (hi - lo)
                if hi - math.log(s.P) > 12.0:
                    raise SimulationError("no outlet pressure gives the head of the curve",
                                          where=self.id)
            return math.exp(brentq(excess, lo, hi, xtol=1e-9))
        return self.polytropic_pressure(pkg, s, head, eta, mw, ratio)

    def polytropic_pressure(self, pkg, s, head, eta, mw, ratio) -> float:
        """Integrate the polytropic path, dH = v dP / eta, in ln P with Heun steps; the
        polytropic head is the integral of v dP. March until it reaches ``head`` and place
        the end inside the last step with a cubic Hermite interpolant (the slope d(head) /
        d(ln P) = P v is known at both ends)."""
        step = math.log(ratio) / POLYTROPIC_STEPS
        zs, x, H = s.zs, math.log(s.P), s.state.H
        pv = s.P * s.state.V / mw  # J/kg: d(head) / d(ln P)
        done = 0.0  # polytropic head so far, J/kg
        for _ in range(20 * POLYTROPIC_STEPS):
            x1 = x + step
            P1 = math.exp(x1)
            pv_guess = P1 * pkg.flash(zs, P=P1, H=H + mw * pv * step / eta).V / mw
            dh = 0.5 * (pv + pv_guess) * step
            pv1 = P1 * pkg.flash(zs, P=P1, H=H + mw * dh / eta).V / mw
            if done + dh >= head:
                def short(t: float, done=done, dh=dh, pv=pv, pv1=pv1) -> float:
                    h00, h10 = 2 * t**3 - 3 * t**2 + 1, t**3 - 2 * t**2 + t
                    h01, h11 = -2 * t**3 + 3 * t**2, t**3 - t**2
                    return (h00 * done + h10 * step * pv + h01 * (done + dh)
                            + h11 * step * pv1 - head)
                return math.exp(x + step * brentq(short, 0.0, 1.0, xtol=1e-12))
            x, H, pv, done = x1, H + mw * dh / eta, pv1, done + dh
        raise SimulationError("no outlet pressure gives the head of the curve", where=self.id)
