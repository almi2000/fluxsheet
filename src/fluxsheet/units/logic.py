"""Logic blocks that manipulate streams without any physics: multiply, duplicate, select.

They do not conserve mass by design (``conserves_mass = False``), so the solver's mass
balance treats them as boundaries, like a Product followed by a Feed.
"""

from __future__ import annotations

from ..errors import SimulationError
from ..stream import Stream
from .base import IN, OUT, Param, Port, UnitOperation, single_inlet

CATEGORY = "Manipulators"


def _copy(s: Stream, factor: float = 1.0) -> Stream:
    return Stream([f * factor for f in s.flows], s.state)


class Mult(UnitOperation):
    type = "Mult"
    label = "Multiplier"
    category = CATEGORY
    description = "Multiplies the flows of a stream by a factor, at the same temperature, "\
        "pressure and composition. For scaling a part of the process, not physical."
    inlets = (IN,)
    outlets = (OUT,)
    params = (Param("factor", "Factor", default=1.0, min=0.0,
                    help="Every component flow is multiplied by this."),)
    conserves_mass = False

    def calculate(self, pkg, inlets):
        return {"out": _copy(single_inlet(inlets), self.p["factor"])}


class Dupl(UnitOperation):
    type = "Dupl"
    label = "Duplicator"
    category = CATEGORY
    description = "Copies a stream to two to four outlets, each identical to the inlet. For "\
        "studying one stream in several ways (e.g. two separation options side by side)."
    inlets = (IN,)
    outlets = tuple(Port(f"out{k}", f"Copy {k}") for k in range(1, 5))
    conserves_mass = False
    drop_unconnected = True  # unconnected copies are not results and not in the balance

    def idle_outlets(self) -> tuple[str, ...]:
        return ("out3", "out4")

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        return {port.name: _copy(s) for port in self.outlets}


class Selector(UnitOperation):
    type = "Selector"
    label = "Selector"
    category = CATEGORY
    description = "Passes one of its inlet streams through unchanged; the others end here. "\
        "For switching between alternative feeds or routes."
    inlets = (Port("in", "Inlets", multiple=True, separate=True),)
    outlets = (OUT,)
    params = (Param("selected", "Selected stream", "stream",
                    help="Id of the inlet stream to pass through."),)
    exclusive_groups = (("selected",),)
    required_groups = exclusive_groups
    conserves_mass = False

    def calculate(self, pkg, inlets):
        ids = self.inlet_ids.get("in", [])
        sid = self.p["selected"]
        if sid not in ids:
            raise SimulationError(
                f"stream {sid!r} does not enter the selector", where=self.id,
                hint=f"select one of: {', '.join(ids) or 'no inlet streams'}")
        return {"out": _copy(inlets["in"][ids.index(sid)])}
