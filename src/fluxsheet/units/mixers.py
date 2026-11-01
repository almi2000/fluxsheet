"""Mixing and splitting streams."""

from __future__ import annotations

from ..stream import Stream
from .base import IN, OUT, Param, Port, ResultField, UnitOperation, mix_streams, single_inlet
from .separators import Flash


class Mixer(UnitOperation):
    type = "Mixer"
    label = "Mixer"
    category = "Mixers & Splitters"
    description = "Combines any number of streams adiabatically. Outlet pressure defaults to the "\
        "lowest inlet pressure."
    inlets = (Port("in", "Inlets", multiple=True, separate=True),)
    outlets = (OUT,)
    params = (Param("P", "Outlet pressure", "pressure",
                    help="Defaults to the lowest inlet pressure.", min=1.0),)

    def calculate(self, pkg, inlets):
        return {"out": mix_streams(pkg, inlets["in"], self.p["P"])}


class Splitter(UnitOperation):
    type = "Splitter"
    label = "Splitter"
    category = "Mixers & Splitters"
    description = "Divides one stream into two of identical composition and state."
    inlets = (IN,)
    outlets = (Port("out1", "Outlet 1"), Port("out2", "Outlet 2"))
    params = (Param("split", "Fraction to outlet 1", "fraction", default=0.5, min=0.0, max=1.0),)

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        f = self.p["split"]
        return {
            "out1": Stream([x * f for x in s.flows], s.state),
            "out2": Stream([x * (1.0 - f) for x in s.flows], s.state),
        }


class MixingTank(Flash):
    type = "MixingTank"
    label = "Mixing tank"
    category = "Mixers & Splitters"
    description = "Stirred tank that blends any number of feeds into one liquid, with a vent "        "for any vapor. Adiabatic unless a temperature, a duty or a heated or cooled jacket "        "is given. With the liquid volume it reports the mean residence time."
    outlets = (Port("liquid", "Outlet"), Port("vapor", "Vent"))
    params = (
        *Flash.params,
        Param("V", "Liquid volume", "volume", min=0.0,
              help="Liquid hold-up of the tank; gives the mean residence time."),
    )
    result_fields = (*Flash.result_fields, ResultField("residence_time", "Residence time", "time"))

    def idle_outlets(self) -> tuple[str, ...]:
        """The vent may stay unconnected while no vapor forms."""
        return ("vapor",)

    def calculate(self, pkg, inlets):
        out = super().calculate(pkg, inlets)
        liquid = out["liquid"]
        volume_flow = liquid.total_flow * liquid.state.V
        if self.p["V"] is not None and volume_flow > 0.0:
            self.results["residence_time"] = self.p["V"] / volume_flow
        return out
