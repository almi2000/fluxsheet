"""Dynamic models of units: hold-up balances that start from the steady state.

A dynamic model wraps a solved steady-state unit. Its state vector is integrated over time;
``outlets`` gives the outlet streams from the state (and, for units without a hold-up of
material, from the inlets), ``derivatives`` the time derivatives of the state.

- ``Vessel`` (Flash, MixingTank): component hold-ups and enthalpy of the contents at the
  unit's pressure. The liquid leaves at a set flow (a controller usually moves it) or by
  gravity; vapor leaves as fast as it forms, with a short vapor residence time.
- ``StirredReactor`` (CSTR): full of the reacting phase; the outflow keeps the volume.
- ``Lag`` (Heater): the outlet enthalpy follows the steady-state outlet with a first-order
  time constant (the thermal mass of the exchanger), so a duty step shows up gradually.

All other units stay algebraic: they are recalculated from their inlets at every instant.
Pressure is fixed at its steady-state value (a flow-driven simulation).
"""

from __future__ import annotations

import math
from dataclasses import replace
from typing import Any, ClassVar

import numpy as np

from ..errors import SimulationError
from ..quantities import UnitError
from ..stream import Stream
from ..thermo_package import PropertyPackage
from ..units.base import Param, ResultField, UnitOperation, empty_stream

G = 9.80665  # m/s2


class DynamicModel:
    """Base class; one subclass per kind of hold-up, registered in ``DYNAMIC_MODELS``."""

    params: ClassVar[tuple[Param, ...]] = ()
    result_fields: ClassVar[tuple[ResultField, ...]] = ()
    # Outlets depend on the state alone (a hold-up of material). Such units break the
    # flowsheet into pieces that are calculated in flow order from them.
    source: ClassVar[bool] = True

    def __init__(self, unit: UnitOperation, pkg: PropertyPackage, settings: dict[str, Any]):
        self.unit, self.pkg = unit, pkg
        self.results: dict[str, float] = {}
        known = {p.name: p for p in self.params}
        for name in settings:
            if name not in known:
                raise SimulationError(
                    f"unknown dynamic setting {name!r}", where=unit.id,
                    hint=f"{unit.label} takes: {', '.join(known) or 'none'}",
                )
        self.s: dict[str, Any] = {}
        for p in self.params:
            value = settings.get(p.name, p.default)
            try:
                self.s[p.name] = None if value is None else unit._convert(p, value)
            except UnitError as exc:
                raise SimulationError(f"{p.label}: {exc}", where=unit.id) from exc
            v = self.s[p.name]
            if isinstance(v, float) and ((p.min is not None and v < p.min) or
                                         (p.max is not None and v > p.max)):
                raise SimulationError(f"{p.label} is out of range ({v:g})", where=unit.id)

    @classmethod
    def schema(cls) -> dict:
        return {"params": [p.to_dict() for p in cls.params],
                "result_fields": [{"name": f.name, "label": f.label, "quantity": f.quantity}
                                  for f in cls.result_fields]}

    # Overridden by the models:
    def initialize(self, inlets: dict[str, list[Stream]], outlets: dict[str, Stream]
                   ) -> np.ndarray:
        raise NotImplementedError

    def scale(self, x0: np.ndarray) -> np.ndarray:
        """Typical magnitude of each state, for error control."""
        return np.maximum(np.abs(x0), 1e-6 * max(np.abs(x0).max(), 1e-10))

    def outlets(self, x: np.ndarray, inlets: dict[str, list[Stream]] | None
                ) -> dict[str, Stream]:
        raise NotImplementedError

    def derivatives(self, x: np.ndarray, inlets: dict[str, list[Stream]]) -> np.ndarray:
        raise NotImplementedError

    def set(self, name: str, value: float) -> bool:
        """Change a dynamic setting (controllers, events); False if there is none called so."""
        if name in self.s:
            self.s[name] = value
            return True
        return False


def _in_totals(inlets: dict[str, list[Stream]], n: int, hf: np.ndarray | None = None
               ) -> tuple[np.ndarray, float]:
    """Component flows and enthalpy flow (with formation enthalpies if ``hf``) entering."""
    flows, h = np.zeros(n), 0.0
    for streams in inlets.values():
        for s in streams:
            f = np.asarray(s.flows)
            flows += f
            h += s.enthalpy_flow + (float(f @ hf) if hf is not None else 0.0)
    return flows, h


class _ContentsFlash:
    """PH flashes along a trajectory, which change little from one call to the next.

    A full flash (with the phase stability test) costs milliseconds and an integration needs
    thousands. While the contents stay one phase, near a temperature and composition where a
    full flash confirmed that phase, the temperature is found on that phase's enthalpy alone
    by Newton steps; the full flash runs again once they drift further.
    """

    T_WINDOW = 2.0  # K around the last confirmed state
    Z_WINDOW = 0.01  # mole fraction

    def __init__(self, pkg: PropertyPackage, T: float):
        self.pkg, self.T = pkg, T
        self.kind: str | None = None  # single phase confirmed by the last full flash
        self.checked: tuple[float, np.ndarray] | None = None

    def __call__(self, zs, P: float, H: float):
        zs = np.asarray(zs, dtype=float)
        if self.kind is not None and self._near(self.T, zs):
            T = self._single_phase_T(zs, P, H)
            if T is not None and self._near(T, zs):
                self.T = T
                return self.pkg.phase_state(self.kind, T, P, list(zs))
        st = self.pkg.flash_ph_near(list(zs), P, H, self.T)
        self.T = st.T
        single = st.phase in ("liquid", "vapor") and not st.liquid_split
        self.kind = st.phase if single else None
        self.checked = (st.T, zs.copy())
        return st

    def _near(self, T: float, zs: np.ndarray) -> bool:
        T0, z0 = self.checked
        return abs(T - T0) < self.T_WINDOW and float(np.abs(zs - z0).max()) < self.Z_WINDOW

    def _single_phase_T(self, zs: np.ndarray, P: float, H: float) -> float | None:
        pkg, kind, T = self.pkg, self.kind, self.T
        z = list(zs)
        h0 = pkg.phase_props(kind, T, P, z).H
        cp = (pkg.phase_props(kind, T + 0.01, P, z).H - h0) / 0.01
        for _ in range(8):
            if cp <= 0.0:
                return None
            step = (H - h0) / cp
            T += step
            if abs(step) < 1e-9 * T:
                return T
            h0 = pkg.phase_props(kind, T, P, z).H
        return None


# ---- vessels -----------------------------------------------------------------------------

class Vessel(DynamicModel):
    params = (
        Param("diameter", "Vessel diameter", "length", min=0.0, default=1.0,
              help="Inner diameter of the vertical cylindrical vessel."),
        Param("height", "Vessel height", "length", min=0.0, default=2.0,
              help="Height of the cylinder; the level can rise up to it."),
        Param("initial_level", "Initial level", "fraction", default=0.5, min=0.0, max=1.0,
              help="Liquid level at the start, as a fraction of the height."),
        Param("liquid_outflow", "Liquid outflow", "choice", default="fixed",
              options=("fixed", "gravity"),
              help="fixed: the liquid leaves at the set flow (a level controller usually "
              "moves it); gravity: the flow follows the square root of the level, as through "
              "an outlet valve at fixed opening."),
        Param("liquid_flow", "Liquid outflow rate", "molar_flow", min=0.0,
              help="Liquid leaving with a fixed outflow. Defaults to the steady state."),
        Param("duty", "Duty", "power",
              help="Heat added. Defaults to the steady state; a jacket with a utility "
              "follows the contents' temperature instead."),
        Param("vapor_time", "Vapor residence time", "time", default=1.0, min=1e-3,
              help="How long vapor stays in the vessel before it leaves; short."),
    )
    result_fields = (
        ResultField("level", "Liquid level", "length"),
        ResultField("level_fraction", "Liquid level (fraction of height)", "fraction"),
        ResultField("T", "Temperature", "temperature"),
        ResultField("holdup", "Liquid hold-up", "amount"),
        ResultField("liquid_volume", "Liquid volume", "volume"),
        ResultField("liquid_flow", "Liquid outflow", "molar_flow"),
        ResultField("vapor_flow", "Vapor outflow", "molar_flow"),
        ResultField("duty", "Duty", "power"),
    )
    liquid_port, vapor_port = "liquid", "vapor"

    def initialize(self, inlets, outlets):
        pkg = self.pkg
        liq = outlets[self.liquid_port]
        vap = outlets.get(self.vapor_port) or empty_stream(pkg, liq.T, liq.P)
        self.P = liq.P
        self.area = math.pi * self.s["diameter"] ** 2 / 4.0
        T = liq.T if liq.total_flow > 0 else vap.T
        self.flash = _ContentsFlash(pkg, T)
        x = np.asarray(liq.zs)
        y = np.asarray(vap.zs)
        level = self.s["initial_level"] * self.s["height"]
        v_liq = pkg.phase_props("liquid", T, self.P, list(x)).V
        n_liq = level * self.area / v_liq
        n_vap = vap.total_flow * self.s["vapor_time"]
        h_liq = liq.state.H if liq.total_flow > 0 else pkg.phase_props("liquid", T, self.P,
                                                                        list(x)).H
        self.level0 = level
        self.flow0 = liq.total_flow
        if self.s["liquid_flow"] is None:
            self.s["liquid_flow"] = liq.total_flow
        if self.s["duty"] is None:
            self.s["duty"] = self.unit.results.get("duty") or 0.0
        self._utility()
        N = n_liq * x + n_vap * y
        H = n_liq * h_liq + n_vap * (vap.state.H if vap.total_flow > 0 else 0.0)
        return np.r_[N, H]

    def _utility(self) -> None:
        """A jacket rated with a utility keeps its UA; its duty then follows the contents."""
        u = self.unit
        self.ut = getattr(u, "utility", None)
        self.UA = u.results.get("UA") if self.ut is not None else None

    def scale(self, x0):
        N = x0[:-1]
        floor = 1e-6 * max(N.sum(), 1e-6)
        h_scale = max(abs(x0[-1]), N.sum() * 8.314 * 300.0)
        return np.r_[np.maximum(N, floor), h_scale]

    def _state(self, x):
        N = np.maximum(x[:-1], 0.0)
        total = N.sum()
        if total <= 0.0:
            raise SimulationError("the vessel ran empty", where=self.unit.id,
                                  hint="check the outflow, or add a level controller")
        st = self.flash(N / total, self.P, x[-1] / total)
        return N, total, st

    def _duty(self, T: float) -> float:
        if self.UA:
            return self.ut.rate(self.UA, T, T, self.unit.p.get("arrangement") or "countercurrent",
                                self.unit.p.get("utility_flow"))
        return self.s["duty"]

    def outlets(self, x, inlets=None):
        pkg = self.pkg
        N, total, st = self._state(x)
        VF = st.VF
        n_liq, n_vap = total * (1.0 - VF), total * VF
        xs = np.asarray(st.x) if st.x else N / total
        ys = np.asarray(st.y) if st.y else N / total
        v_liq = pkg.phase_props("liquid", st.T, self.P, list(xs)).V if n_liq > 0 else 0.0
        volume = n_liq * v_liq
        level = volume / self.area
        if self.s["liquid_outflow"] == "gravity":
            F_liq = self.flow0 * math.sqrt(max(level, 0.0) / self.level0) if self.level0 else 0.0
        else:
            F_liq = self.s["liquid_flow"]
        # An empty vessel cannot deliver: the outflow fades over the last centimetre.
        F_liq *= min(1.0, max(level, 0.0) / 0.01)
        F_vap = n_vap / self.s["vapor_time"]
        h_vap = pkg.phase_props("vapor", st.T, self.P, list(ys)).H if n_vap > 0 else 0.0
        h_liq = (st.H - VF * h_vap) / (1.0 - VF) if VF < 1.0 else 0.0
        self._out = (F_liq, xs, h_liq, F_vap, ys, h_vap)
        self.results.update(level=level, level_fraction=level / self.s["height"], T=st.T,
                            holdup=n_liq, liquid_volume=volume, liquid_flow=F_liq,
                            vapor_flow=F_vap, duty=self._duty(st.T))
        if level > self.s["height"]:
            self.results["overfilled"] = 1.0
        return {
            self.liquid_port: _pinned(pkg, F_liq * xs, st.T, self.P, 0.0),
            self.vapor_port: _pinned(pkg, F_vap * ys, st.T, self.P, 1.0),
        }

    def derivatives(self, x, inlets):
        F_liq, xs, h_liq, F_vap, ys, h_vap = self._out
        f_in, h_in = _in_totals(inlets, self.pkg.n)
        dN = f_in - F_liq * xs - F_vap * ys
        dH = h_in - F_liq * h_liq - F_vap * h_vap + self.results["duty"]
        return np.r_[dN, dH]


def _pinned(pkg, flows, T: float, P: float, vf: float) -> Stream:
    """A product of one phase at the vessel's T and P (no re-flash across saturation)."""
    total = float(np.sum(flows))
    kind = "vapor" if vf == 1.0 else "liquid"
    if total <= 0.0:  # no flow: any composition will do, without a flash
        zs = [1.0 / pkg.n] * pkg.n
        return Stream([0.0] * pkg.n, pkg.phase_state(kind, T, P, zs))
    zs = [float(f) / total for f in flows]
    st = pkg.phase_state(kind, T, P, zs)
    st = replace(st, VF=vf, phase=kind, x=None if vf == 1.0 else tuple(zs),
                 y=tuple(zs) if vf == 1.0 else None)
    return Stream([float(f) for f in flows], st)


class MixingTankVessel(Vessel):
    liquid_port, vapor_port = "liquid", "vapor"


# ---- stirred reactor ---------------------------------------------------------------------

class StirredReactor(DynamicModel):
    """CSTR: the reacting phase fills the volume. States: component hold-ups and enthalpy
    including formation enthalpies, so the heat of reaction is in the balance."""

    params = (
        Param("duty", "Duty", "power",
              help="Heat added. Defaults to the steady state; a jacket with a utility "
              "follows the reactor temperature instead."),
    )
    result_fields = (
        ResultField("T", "Temperature", "temperature"),
        ResultField("holdup", "Hold-up", "amount"),
        ResultField("outflow", "Outflow", "molar_flow"),
        ResultField("duty", "Duty", "power"),
    )

    def initialize(self, inlets, outlets):
        u, pkg = self.unit, self.pkg
        out = outlets["out"]
        self.P = out.P
        self.V = u.p["V"]
        u.setup(pkg)
        feed = np.zeros(pkg.n)
        for streams in inlets.values():
            for s in streams:
                feed += np.asarray(s.flows)
        species = sorted({i for r in u._rxns for i in r.involved})
        self.hf = u.formation_enthalpies(pkg, [feed, out.flows], species)
        self.nu = u._nu
        z = np.asarray(out.zs)
        Vm = pkg.phase_props(u._kind, out.T, self.P, list(z)).V
        N = z * self.V / Vm
        self.T = out.T
        H = N.sum() * self._h(out.T, z)
        tau = self.V / max(out.total_flow * Vm, 1e-300)
        self.relax = max(1e-2 * tau, 1e-3)  # time constant that holds the volume
        if self.s["duty"] is None:
            self.s["duty"] = u.results.get("duty") or 0.0
        self.ut = getattr(u, "utility", None)
        self.UA = u.results.get("UA") if self.ut is not None else None
        return np.r_[N, H]

    def scale(self, x0):
        N = x0[:-1]
        floor = 1e-6 * max(N.sum(), 1e-6)
        return np.r_[np.maximum(N, floor), max(abs(x0[-1]), N.sum() * 8.314 * 300.0)]

    def _h(self, T: float, z: np.ndarray) -> float:
        """Molar enthalpy of the reacting phase with formation enthalpies, J/mol."""
        return self.pkg.phase_props(self.unit._kind, T, self.P, list(z)).H + float(z @ self.hf)

    def _temperature(self, z: np.ndarray, h: float) -> float:
        T = self.T
        for _ in range(50):
            f = self._h(T, z) - h
            dT = 0.01
            slope = (self._h(T + dT, z) - self._h(T, z)) / dT
            step = -f / slope
            step = max(min(step, 50.0), -50.0)
            T = max(T + step, 20.0)
            if abs(step) < 1e-7 * T:
                break
        else:
            raise SimulationError("could not find the reactor temperature", where=self.unit.id)
        self.T = T
        return T

    def outlets(self, x, inlets=None):
        u, pkg = self.unit, self.pkg
        N = np.maximum(x[:-1], 0.0)
        total = N.sum()
        if total <= 0.0:
            raise SimulationError("the reactor ran empty", where=u.id)
        z = N / total
        T = self._temperature(z, x[-1] / total)
        Vm = pkg.phase_props(u._kind, T, self.P, list(z)).V
        rates = np.array([r.rate(pkg, T, self.P, z, Vm, u._phase) for r in u._rxns])
        self._gen = self.V * rates @ self.nu  # mol/s produced per component
        self._excess = (total - self.V / Vm) / self.relax
        Q = (self.ut.rate(self.UA, T, T, u.p.get("arrangement") or "countercurrent",
                          u.p.get("utility_flow")) if self.UA else self.s["duty"])
        self._z, self._h_molar, self._Q = z, x[-1] / total, Q
        self.results.update(T=T, holdup=total, duty=Q)
        if inlets is not None:
            self._outflow(inlets)
        # The outlet is the reacting phase as it is in the reactor; no flash needed.
        return {"out": _pinned(pkg, self._F_out * z, T, self.P,
                               1.0 if u._kind == "vapor" else 0.0)}

    def _outflow(self, inlets) -> None:
        f_in, _ = _in_totals(inlets, self.pkg.n)
        self._F_out = max(f_in.sum() + self._gen.sum() + self._excess, 0.0)
        self.results["outflow"] = self._F_out

    source = False  # the outflow follows the inflow at every instant

    def derivatives(self, x, inlets):
        f_in, h_in = _in_totals(inlets, self.pkg.n, self.hf)
        dN = f_in - self._F_out * self._z + self._gen
        dH = h_in - self._F_out * self._h_molar + self._Q
        return np.r_[dN, dH]


# ---- first-order lag ---------------------------------------------------------------------

class Lag(DynamicModel):
    """Outlet enthalpy follows the steady-state outlet with a time constant."""

    params = (
        Param("time_constant", "Time constant", "time", default=60.0, min=1e-3,
              help="Thermal time constant of the exchanger: about the time the outlet "
              "temperature needs for 63 % of a step."),
    )
    result_fields = (ResultField("T_out", "Outlet temperature", "temperature"),
                     ResultField("duty", "Duty", "power"))
    source = False

    def initialize(self, inlets, outlets):
        out = outlets["out"]
        self.flash = _ContentsFlash(self.pkg, out.T)
        return np.array([out.state.H])

    def scale(self, x0):
        return np.array([max(abs(x0[0]), 8.314 * 300.0)])

    def outlets(self, x, inlets=None):
        u = self.unit
        duty = u.p.get("duty")
        f_in, h_in = _in_totals(inlets, self.pkg.n)
        F = f_in.sum()
        if duty is not None and getattr(u, "utility", None) is None and F > 0.0:
            # A set duty needs no flash for the steady outlet's enthalpy.
            P = u.outlet_pressure(min(s.P for ss in inlets.values() for s in ss))
            self._target, flows = (h_in + duty) / F, list(f_in)
        else:
            target = self.steady(inlets)["out"]
            self._target, flows, P = target.state.H, list(target.flows), target.P
            duty = u.results.get("duty", 0.0)
            if target.total_flow <= 0.0:
                return {"out": target}
        st = self.flash(np.asarray(flows) / F, P, float(x[0]))
        self.results.update(T_out=st.T, duty=duty)
        return {"out": Stream(flows, st)}

    def derivatives(self, x, inlets):
        return np.array([(self._target - x[0]) / self.s["time_constant"]])


DYNAMIC_MODELS: dict[str, type[DynamicModel]] = {
    "Flash": Vessel,
    "MixingTank": MixingTankVessel,
    "CSTR": StirredReactor,
    "Heater": Lag,
}
