"""Pipe: pressure drop from friction, height change and acceleration, with optional heat loss.

One-dimensional homogeneous flow at constant mass flux G (kg/m2/s), with v the specific
volume, marched in segments:

    momentum   dP = -f dL/D G^2 v / 2 - g dz / v - G^2 dv
    energy     h + (G v)^2 / 2 + g z changes by the heat taken up

The last momentum term is the pressure the fluid needs to speed up as it expands. It is
negligible for liquids and for gases well below the speed of sound, and it is what chokes a
gas or flashing flow: no pressure at the outlet can pass more than the flow at Mach 1.

A stream carries no velocity, so its enthalpy counts as the total (stagnation) enthalpy:
inside the pipe the static enthalpy is lower by the kinetic energy, and the outlet stream
gets the kinetic energy back as heat at the outlet pressure, as when the pipe ends in a
vessel or a wider line. The flowsheet's energy balance therefore closes on enthalpies.
"""

from __future__ import annotations

import math

from scipy.optimize import brentq

from ..errors import SimulationError
from ..stream import Stream
from .base import DUTY, IN, OUT, Param, ResultField, UnitOperation, empty_stream, single_inlet

G = 9.80665  # m/s2
_ACCELERATION_TOL = 1e-7  # of the pressure: the acceleration pressure drop has converged
_KE_TOL = 1e-6  # J/mol: the kinetic energy has converged
_MAX_ACCELERATION_STEPS = 200  # successive substitution converges like Mach^2


def friction_factor(Re: float, rel_roughness: float) -> float:
    """Darcy friction factor: 64/Re in laminar flow, Swamee-Jain in turbulent flow."""
    if Re < 2300.0:
        return 64.0 / max(Re, 1e-9)
    return 0.25 / math.log10(rel_roughness / 3.7 + 5.74 / Re ** 0.9) ** 2


class Pipe(UnitOperation):
    type = "Pipe"
    label = "Pipe"
    category = "Pressure Changers"
    description = "Pressure drop along a pipe: friction (Darcy-Weisbach; laminar below "\
        "Re 2300, else the Swamee-Jain friction factor), the height change and the "\
        "acceleration of an expanding gas, with optional heat loss to the surroundings. "\
        "Reports the outlet Mach number and stops when the flow chokes. Two-phase flow is "\
        "treated as homogeneous."
    inlets = (IN,)
    outlets = (OUT,)
    params = (
        Param("length", "Length", "length", min=0.0),
        Param("diameter", "Inner diameter", "length", min=0.0),
        Param("roughness", "Wall roughness", "length", default="0.046 mm", min=0.0,
              help="Absolute roughness; 0.046 mm is commercial steel."),
        Param("elevation", "Height change", "length", default=0.0,
              help="Outlet minus inlet height; positive when the pipe rises."),
        Param("U", "Heat loss coefficient U", "heat_transfer_coefficient", min=0.0,
              help="Through the wall, per inner wall area. Leave empty for an adiabatic pipe."),
        Param("T_ambient", "Surrounding temperature", "temperature", min=0.0, requires="U"),
        Param("segments", "Segments", "integer", default=20, min=1, max=500),
    )
    exclusive_groups = (("length",), ("diameter",))
    required_groups = exclusive_groups
    result_fields = (
        ResultField("dP", "Pressure drop", "pressure_difference"),
        DUTY,
        ResultField("velocity_in", "Inlet velocity", "velocity"),
        ResultField("velocity_out", "Outlet velocity", "velocity"),
        ResultField("reynolds", "Reynolds number (inlet)"),
        ResultField("friction_factor", "Darcy friction factor (inlet)"),
        ResultField("mach", "Mach number (outlet)"),
    )

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        L, D = self.p["length"], self.p["diameter"]
        if D <= 0.0 or L <= 0.0:
            raise SimulationError("length and diameter must be positive", where=self.id)
        if self.p["U"] and self.p["T_ambient"] is None:
            raise SimulationError("the heat loss coefficient is given but not the "
                                  "surrounding temperature", where=self.id)
        if s.total_flow <= 0.0:
            self.results.update(dP=0.0, duty=0.0)
            return {"out": empty_stream(pkg, s.T, s.P)}
        n = self.p["segments"]
        area = math.pi * D * D / 4.0
        mw = sum(f * m for f, m in zip(s.flows, pkg.mws)) / 1000.0 / s.total_flow  # kg/mol
        m_dot = s.total_flow * mw  # kg/s
        flux = m_dot / area  # kg/m2/s, constant along the pipe
        dL, dz = L / n, (self.p["elevation"] or 0.0) / n
        eps = (self.p["roughness"] or 0.0) / D

        def kinetic(v: float) -> float:
            """Kinetic energy, J/mol, at specific volume v (m3/kg)."""
            return 0.5 * mw * (flux * v) ** 2

        def choked(at: float, u: float) -> SimulationError:
            return SimulationError(
                f"the flow chokes after {at:.4g} m of the {L:.4g} m pipe: it reaches the "
                f"speed of sound (velocity {u:.3g} m/s), so the pipe cannot pass this flow",
                where=self.id,
                hint="use a larger diameter, a shorter pipe, a lower flow or a higher inlet "
                "pressure")

        # H is the total enthalpy (J/mol), st the static state in the pipe.
        H, Q = s.state.H, 0.0
        st = s.state
        for _ in range(60):  # the static inlet state, lower by the kinetic energy it has
            ke = kinetic(st.V / mw)
            st = pkg.flash_ph_near(s.zs, s.P, H - ke, st.T)
            if abs(kinetic(st.V / mw) - ke) <= _KE_TOL:
                break
        else:
            raise choked(0.0, flux * st.V / mw)
        x, P_prof, T_prof, u_prof = [0.0], [st.P], [st.T], []
        expansion = 0.0  # change of v over the last segment: the first guess for the next
        for k in range(n):
            v1 = st.V / mw
            u = flux * v1
            Re = flux * D / self._viscosity(pkg, st, s.zs)
            f = friction_factor(Re, eps)
            if k == 0:
                self.results.update(velocity_in=u, reynolds=Re, friction_factor=f)
            u_prof.append(u)
            P_wall = st.P - f * dL / D * flux * u / 2.0 - G * dz / v1  # friction and height
            H_adiabatic = H - mw * G * dz
            if not P_wall > 0.01 * s.P:
                raise SimulationError(
                    f"the pressure is used up after {(k + 1) * dL:.4g} m of the {L:.4g} m "
                    f"pipe (velocity {u:.3g} m/s)", where=self.id,
                    hint="use a larger diameter, a shorter pipe or a lower flow")
            # Successive substitution for the specific volume at the end of the segment. It
            # starts from the expansion of the segment before, so that a pass accepted
            # within the tolerances has the acceleration in it already.
            v2 = max(v1 + expansion, 0.5 * v1)
            for _ in range(_MAX_ACCELERATION_STEPS):
                P = P_wall - flux * flux * (v2 - v1)
                if not P > 0.01 * s.P:
                    raise choked((k + 1) * dL, flux * v2)
                # Static enthalpy after the pressure, elevation and velocity change, before
                # heat transfer.
                base_H = H_adiabatic - kinetic(v2)
                if self.p["U"]:
                    ambient = self.p["T_ambient"]
                    conductance = self.p["U"] * math.pi * D * dL
                    ambient_H = pkg.flash(s.zs, T=ambient, P=P).H

                    def residual(h, pressure=P, guess=st.T, base=base_H,
                                 ua=conductance, ambient_T=ambient):
                        state = pkg.flash_ph_near(s.zs, pressure, h, guess)
                        return s.total_flow * (h - base) - ua * (ambient_T - state.T)

                    low, high = sorted((base_H, ambient_H))
                    h_static = brentq(residual, low, high) if high > low else low
                else:
                    h_static = base_H
                new = pkg.flash_ph_near(s.zs, P, h_static, st.T)
                v_new = new.V / mw
                if (flux * flux * abs(v_new - v2) <= _ACCELERATION_TOL * P
                        and abs(kinetic(v_new) - kinetic(v2)) <= _KE_TOL):
                    break
                v2 = v_new
            else:
                raise choked((k + 1) * dL, flux * v2)
            st, expansion = new, v_new - v1
            q = s.total_flow * (st.H - base_H) if self.p["U"] else 0.0
            Q += q
            H = H_adiabatic + q / s.total_flow
            x.append((k + 1) * dL)
            P_prof.append(st.P)
            T_prof.append(st.T)
        u_out = flux * st.V / mw
        u_prof.append(u_out)
        mach = self._mach(pkg, s.zs, st, mw, flux)
        if mach is not None and mach > 0.9:
            self.notes.append(f"the flow leaves at Mach {mach:.2f}, close to choking: a little "
                              "more flow or length and the pipe cannot pass it")
        self.results.update(dP=s.P - st.P, duty=Q, velocity_out=u_out, mach=mach)
        self.add_profile("pipe", "Along the pipe", "Length (m)", x,
                         {"Pressure": ("pressure", P_prof),
                          "Temperature": ("temperature", T_prof),
                          "Velocity": ("velocity", u_prof)})
        # The outlet stream at the total enthalpy: its kinetic energy turns into heat.
        try:
            out = pkg.flash(s.zs, P=st.P, H=H)
        except SimulationError:
            out = pkg.flash_ph_near(s.zs, st.P, H, st.T)
        return {"out": Stream(list(s.flows), out)}

    @staticmethod
    def _mach(pkg, zs, st, mw: float, flux: float) -> float | None:
        """Mach number of the homogeneous flow: G sqrt(-dv/dP) at constant entropy, with
        (dV/dP)_S = (dV/dP)_H + V (dV/dH)_P from two enthalpy flashes. None for a liquid
        the property method treats as incompressible."""
        dP, dH = 1e-3 * st.P, 1e-3 * 8.314462618 * st.T
        try:
            lower = pkg.flash_ph_near(zs, st.P - dP, st.H, st.T)
            hotter = pkg.flash_ph_near(zs, st.P, st.H + dH, st.T)
        except SimulationError:
            return None
        dV_dP = (st.V - lower.V) / dP + st.V * (hotter.V - st.V) / (hotter.H - st.H)
        return flux * math.sqrt(-dV_dP / mw) if dV_dP < 0.0 else None

    @staticmethod
    def _viscosity(pkg, st, zs) -> float:
        """Mixture viscosity, Pa s; two-phase flow by the McAdams mass-quality average."""
        mws = pkg.mws

        def mass_fractions(xs):
            m = [x * mw for x, mw in zip(xs, mws)]
            return [v / sum(m) for v in m]

        props = pkg.props
        mu_l = mu_g = None
        if st.VF < 1.0:
            xs = list(st.x or zs)
            mu_l = props.ViscosityLiquidMixture(T=st.T, P=st.P, zs=xs, ws=mass_fractions(xs))
        if st.VF > 0.0:
            ys = list(st.y or zs)
            mu_g = props.ViscosityGasMixture(T=st.T, P=st.P, zs=ys, ws=mass_fractions(ys))
        if mu_l is not None and mu_g is not None:
            m_vap = st.VF * sum(y * mw for y, mw in zip(st.y, mws))
            m_liq = (1.0 - st.VF) * sum(x * mw for x, mw in zip(st.x, mws))
            quality = m_vap / (m_vap + m_liq)
            mu = 1.0 / (quality / mu_g + (1.0 - quality) / mu_l)
        else:
            mu = mu_l if mu_g is None else mu_g
        if not mu or mu <= 0.0:
            raise SimulationError("no viscosity data for this mixture",
                                  hint="the pipe needs the liquid and gas viscosities of "
                                  "the components")
        return float(mu)
