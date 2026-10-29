"""Liquid-liquid separators."""

from __future__ import annotations

from ..errors import SimulationError
from ..stream import Stream
from .base import (
    DP,
    DUTY,
    IN,
    P_OUT,
    Param,
    Port,
    ResultField,
    UnitOperation,
    empty_stream,
    single_inlet,
)


class Decanter(UnitOperation):
    type = "Decanter"
    label = "Decanter"
    category = "Separators"
    description = "Settles a liquid into its two liquid phases in equilibrium, e.g. an organic "\
        "layer over water, and vents any vapor. Adiabatic unless a temperature or duty is "\
        "given. Needs a property method that predicts liquid-liquid splits (NRTL, UNIQUAC or "\
        "UNIFAC). A feed that stays one liquid leaves entirely through the heavy outlet."
    inlets = (IN,)
    outlets = (Port("vapor", "Vapor"), Port("light", "Light"), Port("heavy", "Heavy"))
    params = (
        Param("T", "Temperature", "temperature", min=0.0),
        Param("duty", "Duty", "power", help="Heat added. Defaults to zero (adiabatic)."),
        P_OUT,
        DP,
    )
    exclusive_groups = (("T", "duty"), ("P", "dP"))
    result_fields = (
        DUTY,
        ResultField("heavy_fraction", "Heavy fraction of liquid"),
        ResultField("vapor_fraction", "Vapor fraction"),
    )

    def calculate(self, pkg, inlets):
        if not pkg.supports_lle:
            raise SimulationError(
                f"the property method {pkg.method!r} cannot predict two liquid phases, so the "
                "decanter would never separate anything",
                where=self.id,
                hint="set the flowsheet's property method to NRTL, UNIQUAC or UNIFAC",
            )
        s = single_inlet(inlets)
        P = self.outlet_pressure(s.P)
        if s.total_flow <= 0.0:
            T = self.p["T"] or s.T
            self.results.update(duty=0.0, heavy_fraction=0.0, vapor_fraction=0.0)
            return {port: empty_stream(pkg, T, P) for port in ("light", "heavy", "vapor")}
        if self.p["T"] is not None:
            mix = Stream.from_flash(pkg, s.flows, T=self.p["T"], P=P)
        else:
            duty = self.p["duty"] or 0.0
            mix = Stream.from_flash(pkg, s.flows, P=P, H=s.state.H + duty / s.total_flow)
        st, F = mix.state, mix.total_flow
        vapor = [F * st.VF * y for y in st.y] if st.y else [0.0] * pkg.n
        light = [F * st.light_fraction * x for x in st.x_light] if st.liquid_split \
            else [0.0] * pkg.n
        # The heavy layer takes the rest, so the component balance closes exactly.
        heavy = [max(f - v - li, 0.0) for f, v, li in zip(mix.flows, vapor, light)]

        out = {}
        for port, flows, kind in (("light", light, "liquid"), ("heavy", heavy, "liquid"),
                                  ("vapor", vapor, "vapor")):
            total = sum(flows)
            if total <= 1e-12 * F:
                out[port] = empty_stream(pkg, st.T, P)
                continue
            # Each layer is one phase at the decanter's T and P; no second flash, which could
            # split a saturated layer again.
            zs = [f / total for f in flows]
            out[port] = Stream(flows, pkg.phase_state(kind, st.T, P, zs))
        liquid = sum(light) + sum(heavy)
        self.results.update(
            duty=sum(o.enthalpy_flow for o in out.values()) - s.enthalpy_flow,
            heavy_fraction=sum(heavy) / liquid if liquid > 0 else 0.0,
            vapor_fraction=st.VF,
        )
        return out
