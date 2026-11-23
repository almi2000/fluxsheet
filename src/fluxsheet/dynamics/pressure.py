"""Pressure-driven flows in a dynamic simulation (``"pressure_driven": true`` under dynamics).

In the default, flow-driven run the flows are those of the steady-state flowsheet (or of a
controller moving them) and every pressure stays where the steady state put it. With
``pressure_driven`` the flows between pressure nodes follow from their pressures:

- **Nodes** have a pressure: a vessel with a dynamic model (``PressureVessel``), a feed and
  a product. A vessel with gas in it at the start carries its pressure as a state: its
  contents fill the vessel, so the pressure follows from the amounts, the internal energy and
  the volume (a UV flash). A vessel full of liquid stays at its pressure, as if vented. Feeds
  and products are boundaries at their steady-state pressures.
- **Lines** run from a node's outlet through valves, pumps, compressors and heaters (one
  inlet and one outlet each) to the next node. A line with a valve that has a flow
  coefficient (``Cv`` or ``Kv``), or a pump with a performance curve, is pressure-driven: its
  flow is the one at which the pressures along it balance. Other units on it keep their
  steady-state pressure change. A line without either keeps its flow as before (a feed's
  set flow, a vessel's liquid outflow setting). The liquid leaving a vessel adds its static
  head, rho g level, to the vessel's pressure.
- **Valves** follow IEC 60534-2-1 with the flow coefficient times the inherent
  characteristic at the opening (linear, equal percentage or quick opening). Liquid:
  Q = Kv f(x) sqrt(dP / (rho / rho_water)) [m3/h, bar]. Gas: W = N6 Kv f(x) Y sqrt(x p1 rho1)
  [kg/h, bar, kg/m3], x = dP / p1, Y = 1 - x / (3 F_gamma xT), choked at x = F_gamma xT
  (F_gamma = gamma / 1.4, gamma from the steady-state gas). Flow runs downstream only; a
  reversed pressure difference gives no flow. Without an opening the valve starts at the one
  that passes the steady-state flow, so an undisturbed run stays at the steady state.
- **Pumps** with a curve deliver the flow at which their head (affinity laws at the current
  speed) meets the line; a pump without a curve keeps its steady-state pressure rise.

The square roots of the sizing equations are rounded off for pressure differences below
1e-4 of the source's steady pressure (laminar-like), which keeps the Jacobian finite when the
pressures across a valve equalize.
"""

from __future__ import annotations

import math
from collections import OrderedDict
from dataclasses import dataclass, field, replace
from types import SimpleNamespace

import numpy as np
from scipy.optimize import brentq

from ..errors import SimulationError
from ..stream import Stream
from ..units.base import Param
from .units import G, Vessel, _pinned

R_GAS = 8.314462618  # J/mol/K
KV_PER_CV = 1.0 / 1.156  # Kv (m3/h, bar) = 0.865 Cv (US gpm, psi)
N6_KV = 31.6  # IEC 60534-2-1 gas equation with Kv: W in kg/h, p in bar, rho in kg/m3
# (27.3 with Cv)
RHO_WATER = 1000.0  # kg/m3, the reference density in the definition of Kv
ROUNDING = 1e-4  # of the source pressure: where the square root of dP is rounded off

# Units a pressure-driven line may run through: one inlet, one outlet, the flow unchanged.
LINE_UNITS = ("Valve", "Pump", "Compressor", "Heater")


# ---- valve equations -----------------------------------------------------------------------

def characteristic(kind: str, x: float, rangeability: float) -> float:
    """Inherent characteristic f(x): fraction of the full flow coefficient at opening x.

    Equal percentage is R^(x-1) shifted by (1-x)/R so that it closes tightly at x = 0;
    above about x = 0.2 it is within 1 % of R^(x-1)."""
    x = min(max(x, 0.0), 1.0)
    if kind == "linear":
        return x
    if kind == "quick_opening":
        return math.sqrt(x)
    if kind == "equal_percentage":
        return rangeability ** (x - 1.0) - (1.0 - x) / rangeability
    raise ValueError(f"unknown valve characteristic {kind!r}")


def opening_for(kind: str, f: float, rangeability: float) -> float:
    """The opening x at which characteristic(x) = f (0 <= f <= 1)."""
    if f <= 0.0:
        return 0.0
    if f >= 1.0:
        return 1.0
    if kind == "linear":
        return f
    if kind == "quick_opening":
        return f * f
    return brentq(lambda x: characteristic(kind, x, rangeability) - f, 0.0, 1.0,
                  xtol=1e-15, rtol=1e-15, maxiter=200)


def _root(d: float, delta: float) -> float:
    """sqrt(d) for d >> delta, rounded off to a finite slope at d = 0; zero for d <= 0."""
    if d <= 0.0:
        return 0.0
    return d / (d * d + delta * delta) ** 0.25


def liquid_valve_flow(kv: float, dp: float, rho: float, delta: float = 0.0) -> float:
    """Volumetric flow, m3/s, of a liquid of density ``rho`` through a valve whose effective
    flow coefficient is ``kv`` (Kv times the characteristic) at pressure drop ``dp`` (Pa):
    Q [m3/h] = Kv sqrt(dP [bar] / (rho / 1000))."""
    return kv / 3600.0 * _root(dp, delta) * math.sqrt(RHO_WATER / (rho * 1e5))


def gas_valve_flow(kv: float, p1: float, p2: float, rho1: float, gamma: float, xT: float,
                   delta: float = 0.0) -> float:
    """Mass flow, kg/s, of a gas through a valve (IEC 60534-2-1, turbulent, with choking).

    ``p1``, ``p2`` in Pa, ``rho1`` the inlet density in kg/m3, ``gamma`` the heat capacity
    ratio, ``xT`` the valve's pressure differential ratio factor at choked flow."""
    if p1 <= 0.0:
        return 0.0
    choked = gamma / 1.4 * xT
    x = min((p1 - p2) / p1, choked)
    if x <= 0.0:
        return 0.0
    Y = 1.0 - x / (3.0 * choked)
    return N6_KV / 3600.0 * kv * Y * _root(x, delta / p1) * math.sqrt(p1 / 1e5 * rho1)


# ---- vessels with a pressure ---------------------------------------------------------------

class _UVFlash:
    """Temperature and pressure of contents that fill a fixed volume.

    Given the composition z, the molar internal energy u and the molar volume v, solves
    h(T, P) - P v = u and v(T, P) = v by Newton steps in (T, ln P) on T-P flashes, which are
    smooth to some 1e-14. While the contents are one phase, the steps run on that phase's
    properties alone and one T-P flash at the answer confirms the phase (as in
    ``ContentsFlash``). The answer is converged to 1e-12 and better, so it does not depend on
    where the steps started: the integrator's finite-difference Jacobian would see that.
    """

    T_WINDOW = 5.0  # K, P_WINDOW: in ln P; Z_WINDOW: mole fraction, around a confirmed phase
    P_WINDOW = 0.1
    Z_WINDOW = 0.01

    def __init__(self, pkg, T: float, P: float, z: np.ndarray, kind: str | None):
        self.pkg, self.T, self.P = pkg, T, P
        self.kind = kind
        self.checked = (T, P, np.asarray(z, dtype=float).copy())
        self.J: dict[str, np.ndarray | None] = {"single": None, "flash": None}

    def __call__(self, z: np.ndarray, u: float, v: float):
        pkg, zl = self.pkg, list(z)
        if self.kind is not None and self._near(self.T, self.P, z):
            sol = self._newton(zl, u, v, "single")
            if sol is not None and self._near(*sol, z):
                st = pkg.flash(zl, T=sol[0], P=sol[1])
                if st.phase == self.kind and not st.liquid_split:
                    self.T, self.P = sol
                    return pkg.phase_state(self.kind, sol[0], sol[1], zl)
        sol = self._newton(zl, u, v, "flash")
        if sol is None:
            raise SimulationError(
                "could not find the pressure of the vessel contents",
                hint="the vessel may be overfilled with liquid or almost empty; check the "
                "valves and the level")
        self.T, self.P = sol
        st = pkg.flash(zl, T=sol[0], P=sol[1])
        single = st.phase in ("liquid", "vapor") and not st.liquid_split
        self.kind = st.phase if single else None
        self.checked = (sol[0], sol[1], np.asarray(z, dtype=float).copy())
        if single:  # the same answer the shortcut gives next time, to the last bits
            again = self._newton(zl, u, v, "single")
            if again is not None:
                self.T, self.P = again
                return pkg.phase_state(self.kind, again[0], again[1], zl)
        return st

    def _near(self, T: float, P: float, z: np.ndarray) -> bool:
        T0, P0, z0 = self.checked
        return (abs(T - T0) < self.T_WINDOW and abs(math.log(P / P0)) < self.P_WINDOW
                and float(np.abs(np.asarray(z) - z0).max()) < self.Z_WINDOW)

    def _props(self, how: str, zl: list[float], T: float, P: float) -> tuple[float, float]:
        if how == "single":
            ph = self.pkg.phase_props(self.kind, T, P, zl)
            return ph.H, ph.V
        st = self.pkg.flash(zl, T=T, P=P)
        return st.H, st.V

    def _newton(self, zl, u: float, v: float, how: str) -> tuple[float, float] | None:
        def resid(T: float, lnP: float) -> np.ndarray:
            P = math.exp(lnP)
            h, vm = self._props(how, zl, T, P)
            return np.array([(h - P * v - u) / (R_GAS * T), vm / v - 1.0])

        def jacobian(T: float, lnP: float, f: np.ndarray) -> np.ndarray:
            dT, dl = 1e-6 * T, 1e-6
            return np.column_stack([(resid(T + dT, lnP) - f) / dT,
                                    (resid(T, lnP + dl) - f) / dl])

        T, lnP = self.T, math.log(self.P)
        try:
            f = resid(T, lnP)
            J = self.J[how]
            fresh = J is None
            if fresh:
                J = jacobian(T, lnP, f)
            last = math.inf
            for _ in range(60):
                step = np.linalg.solve(J, -f)
                scale = min(1.0, 20.0 / max(abs(step[0]), 1e-300),
                            0.3 / max(abs(step[1]), 1e-300))
                Tn, ln_n = T + scale * step[0], lnP + scale * step[1]
                if Tn <= 1.0:
                    Tn = 0.5 * T
                fn = resid(Tn, ln_n)
                small = abs(Tn - T) < 1e-12 * T and abs(ln_n - lnP) < 1e-12
                T, lnP, norm, f = Tn, ln_n, float(np.abs(fn).max()), fn
                if small:
                    self.J[how] = J
                    return T, math.exp(lnP)
                if norm > 0.5 * last and not fresh:  # slow chord steps: renew the Jacobian
                    J, fresh = jacobian(T, lnP, f), True
                elif norm > 0.5 * last:
                    fresh = False
                last = norm
        except (SimulationError, np.linalg.LinAlgError, ValueError, OverflowError):
            pass
        self.J[how] = None
        return None


class PressureVessel(Vessel):
    """A ``Vessel`` in a pressure-driven run.

    With gas in it at the start, its states are the component amounts and the internal
    energy U = H - P V of the contents, which fill the vessel's volume; the pressure follows
    from them. Without gas it stays at its pressure. Lines with a valve or a pump curve set
    the outflows (``forced``); the other outlets follow the vessel's own settings, the vapor
    of a vessel with gas at ``vapor_flow``.
    """

    params = (
        *Vessel.params,
        Param("vapor_flow", "Vapor outflow rate", "molar_flow", min=0.0,
              help="With pressure-driven flows: vapor leaving a vessel with gas in it when "
              "the vapor line has no valve with a flow coefficient. Defaults to the steady "
              "state."),
    )

    def __init__(self, unit, pkg, settings):
        super().__init__(unit, pkg, settings)
        self.given = set(settings)
        self.forced: dict[str, float] = {}
        self._memo: OrderedDict = OrderedDict()

    def initialize(self, inlets, outlets):
        liq = outlets[self.liquid_port]
        vap = outlets.get(self.vapor_port)
        self.gas = vap is not None and vap.total_flow > 0.0
        if self.gas and liq.total_flow <= 0.0 and "initial_level" not in self.given:
            self.s["initial_level"] = 0.0  # a gas holder: no liquid unless asked for
        if not self.gas:
            x0 = super().initialize(inlets, outlets)
            self.volume = self.area * self.s["height"]
            if self.s["vapor_flow"] is None:
                self.s["vapor_flow"] = 0.0
            return x0
        pkg = self.pkg
        self.P = vap.P
        self.area = math.pi * self.s["diameter"] ** 2 / 4.0
        T = vap.T
        level = self.s["initial_level"] * self.s["height"]
        gas_volume = self.area * (self.s["height"] - level)
        if gas_volume <= 1e-3 * self.area * self.s["height"]:
            raise SimulationError(
                "the vessel is full of liquid at the start, so it has no gas to hold its "
                "pressure", where=self.unit.id, hint="lower initial_level")
        x, y = np.asarray(liq.zs), np.asarray(vap.zs)
        n_liq = 0.0
        h_liq = 0.0
        if level > 0.0:
            lp = pkg.phase_props("liquid", T, self.P, list(x))
            n_liq = level * self.area / lp.V
            h_liq = liq.state.H if liq.total_flow > 0 else lp.H
        vp = pkg.phase_props("vapor", T, self.P, list(y))
        n_vap = gas_volume / vp.V
        self.level0, self.flow0 = level, liq.total_flow
        self.n_vap0 = n_vap
        for name, value in (("liquid_flow", liq.total_flow), ("vapor_flow", vap.total_flow),
                            ("duty", self.unit.results.get("duty") or 0.0)):
            if self.s[name] is None:
                self.s[name] = value
        self._utility()
        N = n_liq * x + n_vap * y
        H = n_liq * h_liq + n_vap * vp.H
        total = N.sum()
        z = N / total
        # The volume the model's own mixture takes at the start, so the steady state is its
        # exact answer (the phases' volumes add up to it within the volume translation).
        if n_liq > 0.0:
            st = pkg.flash(list(z), T=T, P=self.P)
            kind = st.phase if st.phase in ("liquid", "vapor") and not st.liquid_split else None
            self.volume = total * st.V
        else:
            kind = "vapor"
            self.volume = gas_volume
        self.uv = _UVFlash(pkg, T, self.P, z, kind)
        return np.r_[N, H - self.P * self.volume]

    # -- contents --

    def _contents(self, x: np.ndarray) -> SimpleNamespace:
        key = x.tobytes()
        hit = self._memo.get(key)
        if hit is not None:
            self._memo.move_to_end(key)
            return hit
        pkg = self.pkg
        if self.gas:
            N = np.maximum(x[:-1], 0.0)
            total = N.sum()
            if total <= 0.0:
                raise SimulationError("the vessel ran empty", where=self.unit.id)
            st = self.uv(N / total, x[-1] / total, self.volume / total)
            P = st.P
        else:
            N, total, st = self._state(x)
            P = self.P
        VF = st.VF
        n_liq, n_vap = total * (1.0 - VF), total * VF
        xs = np.asarray(st.x) if st.x else N / total
        ys = np.asarray(st.y) if st.y else N / total
        v_liq = pkg.phase_props("liquid", st.T, P, list(xs)).V if n_liq > 0 else 0.0
        h_vap = pkg.phase_props("vapor", st.T, P, list(ys)).H if n_vap > 0 else 0.0
        h_liq = (st.H - VF * h_vap) / (1.0 - VF) if VF < 1.0 else 0.0
        c = SimpleNamespace(st=st, T=st.T, P=P, total=total, n_liq=n_liq, n_vap=n_vap, xs=xs,
                            ys=ys, v_liq=v_liq, h_liq=h_liq, h_vap=h_vap,
                            level=n_liq * v_liq / self.area)
        self._memo[key] = c
        if len(self._memo) > 8:
            self._memo.popitem(last=False)
        self.now = c
        return c

    def prepare(self, x: np.ndarray) -> None:
        """Contents at state ``x``, for the lines to read before the outlets are formed."""
        self.now = self._contents(x)

    def inlet_pressure(self) -> float:
        """Pressure at the inlets (above the liquid)."""
        return self.now.P

    def port_conditions(self, port: str) -> tuple[float, float, float]:
        """Pressure, density (kg/m3) and molar mass (kg/mol) of what leaves at ``port``."""
        c, pkg = self.now, self.pkg
        if port == self.liquid_port:
            v = c.v_liq or pkg.phase_props("liquid", c.T, c.P, list(c.xs)).V
            M = float(c.xs @ np.asarray(pkg.mws)) / 1000.0
            return c.P + M / v * G * max(c.level, 0.0), M / v, M
        v = pkg.phase_props("vapor", c.T, c.P, list(c.ys)).V
        M = float(c.ys @ np.asarray(pkg.mws)) / 1000.0
        return c.P, M / v, M

    def outlets(self, x, inlets=None):
        pkg = self.pkg
        c = self._contents(x)
        self.now = c
        level = c.level
        if self.liquid_port in self.forced:
            F_liq = self.forced[self.liquid_port]
        elif self.s["liquid_outflow"] == "gravity":
            F_liq = self.flow0 * math.sqrt(max(level, 0.0) / self.level0) if self.level0 else 0.0
        else:
            F_liq = self.s["liquid_flow"]
        F_liq *= min(1.0, max(level, 0.0) / 0.01)  # an empty vessel cannot deliver
        if self.vapor_port in self.forced:
            F_vap = self.forced[self.vapor_port]
        elif self.gas:
            F_vap = self.s["vapor_flow"]
        else:
            F_vap = c.n_vap / self.s["vapor_time"]
        if self.gas:  # nor can one whose gas has condensed
            F_vap *= min(1.0, c.n_vap / (1e-4 * self.n_vap0))
        self._out = (F_liq, c.xs, c.h_liq, F_vap, c.ys, c.h_vap)
        self.results.update(level=level, level_fraction=level / self.s["height"], T=c.T,
                            holdup=c.n_liq, liquid_volume=c.n_liq * c.v_liq,
                            liquid_flow=F_liq, vapor_flow=F_vap, duty=self._duty(c.T), P=c.P)
        if level > self.s["height"]:
            self.results["overfilled"] = 1.0
        return {
            self.liquid_port: _pinned(pkg, F_liq * c.xs, c.T, c.P, 0.0),
            self.vapor_port: _pinned(pkg, F_vap * c.ys, c.T, c.P, 1.0),
        }


PRESSURE_MODELS = {"Flash": PressureVessel, "MixingTank": PressureVessel}


# ---- lines between nodes -------------------------------------------------------------------

@dataclass
class _Line:
    source: str  # a PressureVessel or a Feed
    port: str
    target: str  # a node with inlet_pressure(), or a Product
    units: list[str]  # in flow order
    phase: str  # "liquid" or "vapor"
    F0: float  # steady-state molar flow
    valve: str | None = None
    pump: str | None = None  # a pump with a curve
    fixed_before: float = 0.0  # steady pressure rise of the other units, before the valve
    fixed_after: float = 0.0  # ... and after it
    pump_before: bool = True
    P_target: float | None = None  # a product's pressure
    k_line: float = 0.0  # Pa / (mol/s)^2: line resistance of a line with a pump only
    kv: float = 0.0
    kind: str = "linear"
    rangeability: float = 50.0
    xT: float = 0.7
    gamma: float = 1.4
    delta: float = 1.0  # Pa, see ROUNDING
    now: dict = field(default_factory=dict)  # the last evaluation


class Network:
    """The pressure-driven lines of a dynamic run, set up from the steady state."""

    def __init__(self, sim):
        self.sim = sim
        self.lines: list[_Line] = []
        self.vessels = {uid: m for uid, m in sim.dynamic.items()
                        if isinstance(m, PressureVessel)}
        self.feed_streams: dict[str, Stream] = {}
        self.valve_p2: dict[str, float] = {}
        self._find_lines()

    # -- set-up --

    def _find_lines(self) -> None:
        sim = self.sim
        conn = {c.id: c for c in sim.connections}
        into: dict[str, int] = {}
        for c in sim.connections:
            if c.dst_unit is not None:
                into[c.dst_unit] = into.get(c.dst_unit, 0) + 1
        starts = [(uid, m.liquid_port) for uid, m in self.vessels.items()]
        starts += [(uid, m.vapor_port) for uid, m in self.vessels.items()]
        starts += [(uid, "out") for uid, u in sim.units.items() if u.type == "Feed"]
        on_line: set[str] = set()
        for uid, port in starts:
            units, streams, target = [], [], None
            sid = sim.by_source.get((uid, port))
            while sid is not None and sid in conn:
                streams.append(sid)
                dst = conn[sid].dst_unit
                if dst is None:
                    break
                if dst in sim.dynamic:
                    target = dst if hasattr(sim.dynamic[dst], "inlet_pressure") else None
                    break
                unit = sim.units[dst]
                if unit.type == "Product":
                    target = dst
                    break
                if (unit.type not in LINE_UNITS or len(unit.inlets) != 1
                        or len(unit.outlets) != 1 or into.get(dst, 0) != 1 or dst in units):
                    break
                units.append(dst)
                sid = sim.by_source.get((dst, unit.outlets[0].name))
            valves = [u for u in units if sim.units[u].type == "Valve" and _kv(sim.units[u])]
            pumps = [u for u in units
                     if sim.units[u].type == "Pump" and sim.units[u].p.get("curve") is not None]
            if not valves and not pumps:
                continue
            if target is None:
                continue  # the valves are reported below
            if len(valves) > 1 or len(pumps) > 1:
                raise SimulationError(
                    f"the line from {uid} has more than one "
                    f"{'valve with a flow coefficient' if len(valves) > 1 else 'pump curve'}",
                    where=", ".join(valves if len(valves) > 1 else pumps),
                    hint="give the flow coefficient to one valve only; the others keep their "
                    "steady-state pressure drop")
            self.lines.append(self._line(uid, port, target, units, streams,
                                         valves[0] if valves else None,
                                         pumps[0] if pumps else None))
            on_line.update(units)
        for uid, unit in sim.units.items():
            if unit.type == "Valve" and _kv(unit) and uid not in on_line:
                sim.warnings.append(
                    f"{uid}: the valve's flow does not follow from pressures: it does not sit "
                    "on a line from a vessel with a dynamic model or a feed to such a vessel "
                    "or a product, through valves, pumps, compressors and heaters only")

    def _line(self, source, port, target, units, streams, valve, pump) -> _Line:
        sim = self.sim
        snap = sim.snapshot
        if source in self.vessels:
            vessel = self.vessels[source]
            phase = "liquid" if port == vessel.liquid_port else "vapor"
            if phase == "vapor" and not vessel.gas:
                raise SimulationError(
                    f"{source} holds no gas at the start, so its vapor line cannot be "
                    "pressure-driven", where=valve or pump,
                    hint="remove the flow coefficient, or start the vessel with vapor in it")
        else:
            vf = snap[streams[0]].state.VF
            if 1e-9 < vf < 1.0 - 1e-9:
                raise SimulationError(
                    f"the feed {source} is two-phase; a valve on it needs one phase",
                    where=valve or pump)
            phase = "vapor" if vf >= 0.5 else "liquid"
        line = _Line(source, port, target, units, phase, snap[streams[0]].total_flow,
                     valve=valve, pump=pump)
        if pump is not None and phase == "vapor":
            raise SimulationError("a pump curve on a gas line", where=pump,
                                  hint="pumps move liquid; use a compressor for gas")
        # Steady pressure change of each unit; the others keep theirs (as a pressure drop or
        # rise, so they follow their inlet pressure).
        after = False
        for u, s_in, s_out in zip(units, streams, streams[1:]):
            unit = sim.units[u]
            rise = snap[s_out].P - snap[s_in].P
            if u == valve:
                after = True
                continue
            if u == pump:
                line.pump_before = not after
                continue
            if after:
                line.fixed_after += rise
            else:
                line.fixed_before += rise
            if "dP" in unit.p and "P" in unit.p:
                if unit.type in ("Pump", "Compressor") and rise >= 0.0:
                    unit.p["P"], unit.p["dP"] = None, rise
                elif unit.type in ("Valve", "Heater") and rise <= 0.0:
                    unit.p["P"], unit.p["dP"] = None, -rise
        if target not in sim.dynamic:  # a product: it stays at its steady-state pressure
            line.P_target = snap[streams[-1]].P
        line.delta = ROUNDING * snap[streams[0]].P
        if valve is not None:
            v = sim.units[valve]
            line.kv = _kv(v)
            line.kind = v.p.get("characteristic") or "linear"
            line.rangeability = v.p.get("rangeability") or 50.0
            line.xT = v.p.get("xT") or 0.7
            if phase == "vapor":
                line.gamma = _gamma(sim.sections.package(source), snap[streams[0]])
        return line

    def start(self, x0: dict[str, np.ndarray]) -> None:
        """Valve openings (unless given) that pass the steady-state flows, and the line
        resistance of lines with a pump curve only, so an undisturbed run stays put.
        ``x0``: the initial state of each vessel."""
        sim = self.sim
        for uid, m in self.vessels.items():
            m.prepare(x0[uid])
        for line in self.lines:
            if line.valve is None:
                up = self._source_conditions(line)
                excess = self._pressures(line, line.F0, up)
                drive = excess[0] - excess[1]
                if line.F0 <= 0.0 or drive < -1e-9 * up[0]:
                    raise SimulationError(
                        f"pump {line.pump} cannot deliver the steady-state flow from "
                        f"{line.source} to {line.target}: it falls "
                        f"{-drive / 1e5:.4g} bar short", where=line.pump,
                        hint="check the pump curve and the pressures, or add a valve with a "
                        "flow coefficient to the line")
                line.k_line = max(drive, 0.0) / line.F0 ** 2
                continue
            valve = sim.units[line.valve]
            given = valve.p.get("opening")
            if given is None:
                up = self._source_conditions(line)
                p1, p2 = self._pressures(line, line.F0, up)
                full = self._valve_flow(line, 1.0, p1, p2, up)
                if line.F0 > 0.0 and full <= 0.0:
                    raise SimulationError(
                        f"no pressure difference drives the steady-state flow through "
                        f"{line.valve}: {p1 / 1e5:.5g} bar before it, {p2 / 1e5:.5g} bar after",
                        where=line.valve,
                        hint="the vessel or product downstream is at a higher pressure than "
                        "the one upstream; check their pressures")
                need = line.F0 / full if full > 0.0 else 0.0
                if need > 1.0 + 1e-9:
                    raise SimulationError(
                        f"the valve is too small for the steady-state flow: fully open it "
                        f"passes {full / line.F0:.3g} of it", where=line.valve,
                        hint=f"raise the flow coefficient to at least {need:.3g} times its value")
                valve.p["opening"] = opening_for(line.kind, need, line.rangeability)
            else:
                F = self._solve(line)
                if abs(F - line.F0) > 1e-6 * max(line.F0, 1e-12):
                    sim.warnings.append(
                        f"{line.valve}: at the given opening it passes {F:.6g} mol/s, the "
                        f"steady state {line.F0:.6g} mol/s, so the run does not start at rest")

    def check_controllers(self) -> None:
        """A vessel outflow that a line sets cannot also be a controller's output."""
        for cid, ctrl in self.sim.controllers.items():
            m = ctrl.model.manipulate
            for line in self.lines:
                if line.source == m.unit and (
                        (m.param == "liquid_flow" and line.phase == "liquid")
                        or (m.param == "vapor_flow" and line.phase == "vapor")):
                    raise SimulationError(
                        f"{m.param} of {m.unit} follows from the pressures in a pressure-driven "
                        f"run (line through {line.valve or line.pump})", where=cid,
                        hint=f"manipulate the opening of {line.valve} instead"
                        if line.valve else "manipulate the pump speed instead")

    def dependencies(self) -> dict[str, set[str]]:
        """State blocks that the outflow of each line's source depends on, besides its own."""
        deps: dict[str, set[str]] = {}
        manipulators: dict[str, set[str]] = {}
        for cid, ctrl in self.sim.controllers.items():
            manipulators.setdefault(ctrl.model.manipulate.unit, set()).add(cid)
        for line in self.lines:
            d = deps.setdefault(line.source, set())
            if line.target in self.sim.dynamic:
                d.add(line.target)
            for u in line.units:
                d |= manipulators.get(u, set())
        return deps

    # -- at each instant --

    def update(self, x: np.ndarray) -> None:
        """Pressures from the states, then the flow of every line."""
        sim = self.sim
        for uid, m in self.vessels.items():
            m.prepare(x[sim.slices[uid]])
        for line in self.lines:
            up = self._source_conditions(line)
            F = self._solve(line, up)
            p1, p2 = self._pressures(line, F, up)
            if line.source in self.vessels:
                self.vessels[line.source].forced[line.port] = F
            else:
                s = line.now["feed"]
                self.feed_streams[line.source] = Stream([z * F for z in s.zs], s.state)
            if line.valve is not None:
                self.valve_p2[line.valve] = p2
            line.now.update(F=F, p1=p1, p2=p2)

    def set_valve(self, uid: str, inlets: dict[str, list[Stream]]) -> None:
        """Before the sweep calculates a line's valve: its outlet at the line's pressure."""
        unit = self.sim.units[uid]
        s = inlets["in"][0]
        unit.p["P"], unit.p["dP"] = min(self.valve_p2[uid], s.P), None

    @staticmethod
    def isenthalpic(inlet: Stream, outputs: dict[str, Stream]) -> dict[str, Stream]:
        """The valve's outlet with exactly its inlet's enthalpy: the flash stops within its
        tolerance, and the vessel downstream would integrate that."""
        out = outputs["out"]
        if out.total_flow <= 0.0:
            return outputs
        return {"out": Stream(out.flows, replace(out.state, H=inlet.state.H,
                                                 H_reactive=inlet.state.H_reactive))}

    def valve_results(self) -> dict[str, dict[str, float]]:
        out = {}
        for line in self.lines:
            if line.valve is not None and line.now:
                out[line.valve] = {
                    "opening": float(self.sim.units[line.valve].p["opening"]),
                    "dP": line.now["p1"] - line.now["p2"],
                }
        return out

    # -- one line --

    def _source_conditions(self, line: _Line) -> tuple[float, float, float]:
        """Pressure, density and molar mass at the start of the line."""
        if line.source in self.vessels:
            return self.vessels[line.source].port_conditions(line.port)
        s = self.sim._calc(line.source)["out"]
        line.now["feed"] = s
        M = sum(z * mw for z, mw in zip(s.zs, self.sim.pkg.mws)) / 1000.0
        return s.P, M / s.state.V, M

    def _target_pressure(self, line: _Line) -> float:
        if line.P_target is not None:
            return line.P_target
        return self.sim.dynamic[line.target].inlet_pressure()

    def _pump_rise(self, line: _Line, F: float, rho: float, M: float) -> float:
        pump = self.sim.units[line.pump]
        r = pump.speed_ratio()
        Q = F * M / rho
        return rho * pump.p["curve"].head(Q / r) * r * r

    def _pressures(self, line: _Line, F: float, up) -> tuple[float, float]:
        """Pressure before and after the valve at flow F (with a pump only: at the end of
        the line, and what the target holds against it)."""
        P_up, rho, M = up
        p1 = P_up + line.fixed_before
        p2 = self._target_pressure(line) - line.fixed_after + line.k_line * F * F
        if line.pump is not None:
            rise = self._pump_rise(line, F, rho, M)
            if line.pump_before or line.valve is None:
                p1 += rise
            else:
                p2 -= rise
        return p1, p2

    def _valve_flow(self, line: _Line, f: float, p1: float, p2: float, up) -> float:
        """mol/s through the valve at characteristic f (fraction of Kv)."""
        P_up, rho, M = up
        if line.phase == "liquid":
            return liquid_valve_flow(f * line.kv, p1 - p2, rho, line.delta) * rho / M
        rho1 = rho * p1 / P_up  # the gas expands a little on its way to the valve
        return gas_valve_flow(f * line.kv, p1, p2, rho1, line.gamma, line.xT, line.delta) / M

    def _solve(self, line: _Line, up=None) -> float:
        """The line's flow, at which the pressures balance."""
        if up is None:
            up = self._source_conditions(line)
        if line.valve is not None:
            f = characteristic(line.kind, self.sim.units[line.valve].p["opening"],
                               line.rangeability)

            def excess(F: float) -> float:
                return self._valve_flow(line, f, *self._pressures(line, F, up), up) - F

            if line.pump is None and line.k_line == 0.0:
                return self._valve_flow(line, f, *self._pressures(line, 0.0, up), up)
        else:
            def excess(F: float) -> float:
                p1, p2 = self._pressures(line, F, up)
                return p1 - p2
        e0 = excess(0.0)
        if e0 <= 0.0:
            return 0.0
        hi = max(line.F0, 1e-9)
        for _ in range(80):
            if excess(hi) < 0.0:
                break
            hi *= 2.0
        else:
            raise SimulationError(f"no flow balances the pressures on the line from "
                                  f"{line.source}", where=line.valve or line.pump)
        return brentq(excess, 0.0, hi, xtol=1e-14 * hi, rtol=1e-15, maxiter=200)


def _kv(unit) -> float:
    """The valve's Kv (m3/h at 1 bar of water), from Kv or Cv; 0 if it has neither."""
    if unit.p.get("Kv") is not None:
        return float(unit.p["Kv"])
    if unit.p.get("Cv") is not None:
        return float(unit.p["Cv"]) * KV_PER_CV
    return 0.0


def _gamma(pkg, s: Stream) -> float:
    """Heat capacity ratio cp / (cp - R) of a gas stream at its state."""
    zs = list(s.zs)
    dT = 1e-3 * s.T
    cp = (pkg.phase_props("vapor", s.T + dT, s.P, zs).H
          - pkg.phase_props("vapor", s.T - dT, s.P, zs).H) / (2.0 * dT)
    return cp / (cp - R_GAS)
