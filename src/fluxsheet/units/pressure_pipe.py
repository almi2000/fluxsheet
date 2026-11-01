"""Pipe: pressure drop from friction and height change, with optional heat loss."""

from __future__ import annotations

import math

from ..errors import SimulationError
from ..stream import Stream
from .base import DUTY, IN, OUT, Param, ResultField, UnitOperation, empty_stream, single_inlet

G = 9.80665  # m/s2


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
        "Re 2300, else the Swamee-Jain friction factor) and the height change, with optional "\
        "heat loss to the surroundings. Two-phase flow is treated as homogeneous."
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
        dL, dz = L / n, (self.p["elevation"] or 0.0) / n
        eps = (self.p["roughness"] or 0.0) / D
        st, H, Q = s.state, s.state.H, 0.0
        x, P_prof, T_prof, u_prof = [0.0], [st.P], [st.T], []
        for k in range(n):
            rho = mw / st.V
            u = m_dot / (rho * area)
            Re = rho * u * D / self._viscosity(pkg, st, s.zs)
            f = friction_factor(Re, eps)
            if k == 0:
                self.results.update(velocity_in=u, reynolds=Re, friction_factor=f)
            u_prof.append(u)
            P = st.P - f * dL / D * rho * u * u / 2.0 - rho * G * dz
            if P <= 0.01 * s.P:
                raise SimulationError(
                    f"the pressure is used up after {(k + 1) * dL:.4g} m of the {L:.4g} m "
                    f"pipe (velocity {u:.3g} m/s)", where=self.id,
                    hint="use a larger diameter, a shorter pipe or a lower flow")
            q = 0.0
            if self.p["U"]:
                q = self.p["U"] * math.pi * D * dL * (self.p["T_ambient"] - st.T)
            Q += q
            H += (q - m_dot * G * dz) / s.total_flow  # heat in, lifting work out
            st = pkg.flash_ph_near(s.zs, P, H, st.T)
            x.append((k + 1) * dL)
            P_prof.append(st.P)
            T_prof.append(st.T)
        u_out = m_dot * st.V / (mw * area)
        u_prof.append(u_out)
        self.results.update(dP=s.P - st.P, duty=Q, velocity_out=u_out)
        self.add_profile("pipe", "Along the pipe", "Length (m)", x,
                         {"Pressure": ("pressure", P_prof),
                          "Temperature": ("temperature", T_prof),
                          "Velocity": ("velocity", u_prof)})
        return {"out": Stream(list(s.flows), st)}

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
