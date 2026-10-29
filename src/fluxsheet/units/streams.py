"""Feeds and products: where material enters and leaves the process."""

from __future__ import annotations

from ..errors import SimulationError
from ..stream import Stream
from .base import IN, OUT, Param, UnitOperation


class Feed(UnitOperation):
    type = "Feed"
    label = "Feed"
    category = "Streams"
    description = "A material stream entering the process, at a given state and flow rate."
    outlets = (OUT,)
    params = (
        Param("T", "Temperature", "temperature", help="Feed temperature.", min=0.0),
        Param("VF", "Vapor fraction", "fraction",
              help="Instead of a temperature: 0 = saturated liquid, 1 = saturated vapor.",
              min=0.0, max=1.0),
        Param("P", "Pressure", "pressure", help="Feed pressure.", min=1.0),
        Param("flows", "Component flows", "composition",
              help="Molar flow of each component, e.g. {\"water\": \"10 kmol/h\"}."),
    )
    exclusive_groups = (("T", "VF"), ("P",), ("flows",))
    required_groups = exclusive_groups

    def calculate(self, pkg, inlets):
        flows_in = self.p["flows"]
        unknown = [name for name in flows_in if name not in pkg.names]
        if unknown:
            raise SimulationError(
                f"{', '.join(unknown)} {'is' if len(unknown) == 1 else 'are'} not in the "
                "flowsheet's component list", where=self.id,
                hint="add the compound to the component list, or remove it from this feed",
            )
        flows = [flows_in.get(name, 0.0) for name in pkg.names]
        if sum(flows) <= 0.0:
            raise SimulationError("feed has zero total flow", where=self.id,
                                  hint="give at least one component a positive flow")
        spec = {"P": self.p["P"]}
        spec.update({"T": self.p["T"]} if self.p["T"] is not None else {"VF": self.p["VF"]})
        return {"out": Stream.from_flash(pkg, flows, **spec)}


class Product(UnitOperation):
    type = "Product"
    label = "Product"
    category = "Streams"
    description = "A material stream leaving the process."
    inlets = (IN,)

    def calculate(self, pkg, inlets):
        return {}
