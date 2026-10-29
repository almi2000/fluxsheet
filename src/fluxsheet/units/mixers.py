"""Mixing and splitting streams."""

from __future__ import annotations

from ..stream import Stream
from .base import IN, OUT, Param, Port, UnitOperation, empty_stream, single_inlet


class Mixer(UnitOperation):
    type = "Mixer"
    label = "Mixer"
    category = "Mixers & Splitters"
    description = "Combines any number of streams adiabatically. Outlet pressure defaults to the "\
        "lowest inlet pressure."
    inlets = (Port("in", "Inlets", multiple=True),)
    outlets = (OUT,)
    params = (Param("P", "Outlet pressure", "pressure",
                    help="Defaults to the lowest inlet pressure.", min=1.0),)

    def calculate(self, pkg, inlets):
        streams = inlets["in"]
        P = self.p["P"] if self.p["P"] is not None else min(s.P for s in streams)
        flows = [sum(s.flows[i] for s in streams) for i in range(pkg.n)]
        total = sum(flows)
        if total <= 0.0:
            return {"out": empty_stream(pkg, streams[0].T, P)}
        H = sum(s.enthalpy_flow for s in streams) / total
        return {"out": Stream.from_flash(pkg, flows, P=P, H=H)}


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
