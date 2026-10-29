"""Heating and cooling."""

from __future__ import annotations

from ..stream import Stream
from .base import DP, DUTY, IN, OUT, P_OUT, Param, UnitOperation, empty_stream, single_inlet


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
