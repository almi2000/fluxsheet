"""Phase and component separators."""

from __future__ import annotations

from dataclasses import replace

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


class Flash(UnitOperation):
    type = "Flash"
    label = "Flash drum"
    category = "Separators"
    description = "Separates a stream into vapor and liquid in equilibrium. Adiabatic unless a "\
        "temperature or duty is given."
    inlets = (IN,)
    result_fields = (DUTY, ResultField("vapor_fraction", "Vapor fraction"))
    outlets = (Port("vapor", "Vapor"), Port("liquid", "Liquid"))
    params = (
        Param("T", "Temperature", "temperature", min=0.0),
        Param("duty", "Duty", "power", help="Heat added. Defaults to zero (adiabatic)."),
        P_OUT,
        DP,
    )
    exclusive_groups = (("T", "duty"), ("P", "dP"))

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P = self.outlet_pressure(s.P)
        if s.total_flow <= 0.0:
            T = self.p["T"] or s.T
            self.results.update(duty=0.0, vapor_fraction=0.0)
            return {"vapor": empty_stream(pkg, T, P), "liquid": empty_stream(pkg, T, P)}
        if self.p["T"] is not None:
            mix = Stream.from_flash(pkg, s.flows, T=self.p["T"], P=P)
            duty = mix.enthalpy_flow - s.enthalpy_flow
        else:
            duty = self.p["duty"] or 0.0
            mix = Stream.from_flash(pkg, s.flows, P=P, H=s.state.H + duty / s.total_flow)
        st, F = mix.state, mix.total_flow
        vapor = [F * st.VF * y for y in st.y] if st.y else [0.0] * pkg.n
        liquid = [F * (1.0 - st.VF) * x for x in st.x] if st.x else [0.0] * pkg.n
        out = {}
        for port, flows, vf in (("vapor", vapor, 1.0), ("liquid", liquid, 0.0)):
            out[port] = (
                Stream.from_flash(pkg, flows, T=st.T, P=P) if sum(flows) > 1e-12 * F
                else empty_stream(pkg, st.T, P)
            )
            # Pin each product to its phase: at saturation, a T-P flash may land on either side.
            out[port] = Stream(out[port].flows, _with_phase(out[port].state, vf, flows))
        self.results.update(duty=duty, vapor_fraction=st.VF)
        return out


def _with_phase(state, vf, flows):
    total = sum(flows)
    zs = tuple(f / total for f in flows) if total > 0 else state.x or state.y
    return replace(
        state, VF=vf, phase="vapor" if vf == 1.0 else "liquid",
        x=None if vf == 1.0 else zs, y=zs if vf == 1.0 else None,
    )


class ComponentSplitter(UnitOperation):
    type = "ComponentSplitter"
    label = "Component splitter"
    category = "Separators"
    description = "Idealized separation: sends a chosen fraction of each component to the top. "\
        "Useful for conceptual design before a rigorous column is set up."
    inlets = (IN,)
    result_fields = (DUTY,)
    outlets = (Port("top", "Top"), Port("bottom", "Bottom"))
    params = (
        Param("fractions", "Fraction of each component to top", "component_fractions",
              help="Components not listed go entirely to the bottom."),
        Param("T_top", "Top temperature", "temperature", help="Defaults to the inlet temperature."),
        Param("T_bottom", "Bottom temperature", "temperature",
              help="Defaults to the inlet temperature."),
    )
    exclusive_groups = (("fractions",),)
    required_groups = exclusive_groups

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        fr = [self.p["fractions"].get(name, 0.0) for name in pkg.names]
        top = [x * f for x, f in zip(s.flows, fr)]
        bottom = [x - t for x, t in zip(s.flows, top)]
        out = {}
        for port, flows, T in (("top", top, self.p["T_top"]), ("bottom", bottom, self.p["T_bottom"])):
            T = T if T is not None else s.T
            out[port] = Stream.from_flash(pkg, flows, T=T, P=s.P) if sum(flows) > 0 else \
                empty_stream(pkg, T, s.P)
        self.results["duty"] = sum(o.enthalpy_flow for o in out.values()) - s.enthalpy_flow
        return out
