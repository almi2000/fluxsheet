"""Pressure changers: valves, pumps, compressors."""

from __future__ import annotations

from ..errors import SimulationError
from ..stream import Stream
from ..thermo_package import PropertyPackage
from .base import IN, OUT, P_OUT, POWER, Param, UnitOperation, empty_stream, single_inlet


class Valve(UnitOperation):
    type = "Valve"
    label = "Valve"
    category = "Pressure Changers"
    description = "Isenthalpic pressure reduction (throttling). The outlet may partially "\
        "vaporize and cool."
    inlets = (IN,)
    outlets = (OUT,)
    params = (P_OUT, Param("dP", "Pressure drop", "pressure_difference", min=0.0))
    exclusive_groups = (("P", "dP"),)
    required_groups = exclusive_groups

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P = self.outlet_pressure(s.P)
        if P > s.P:
            raise SimulationError(
                f"outlet pressure ({P / 1e5:.3g} bar) is above inlet pressure "
                f"({s.P / 1e5:.3g} bar)", where=self.id,
                hint="a valve can only reduce pressure; use a pump or compressor to raise it",
            )
        if s.total_flow <= 0.0:
            return {"out": empty_stream(pkg, s.T, P)}
        return {"out": Stream.from_flash(pkg, s.flows, P=P, H=s.state.H)}


class _Machine(UnitOperation):
    result_fields = (POWER,)
    inlets = (IN,)
    outlets = (OUT,)
    exclusive_groups = (("P", "dP"),)
    required_groups = exclusive_groups
    raises_pressure = True

    def calculate(self, pkg: PropertyPackage, inlets):
        s = single_inlet(inlets)
        P = self.p["P"] if self.p["P"] is not None else s.P + self.p["dP"]
        if P < s.P:
            raise SimulationError(
                f"outlet pressure ({P / 1e5:.3g} bar) is below inlet pressure "
                f"({s.P / 1e5:.3g} bar)", where=self.id,
                hint="use a valve to reduce pressure",
            )
        if s.total_flow <= 0.0:
            self.results["power"] = 0.0
            return {"out": empty_stream(pkg, s.T, P)}
        if s.state.VF < 1e-9:
            # Liquid: incompressible isentropic work, V dP. Exact enough and avoids an
            # isentropic flash on a liquid, which is ill-conditioned.
            ideal_work = s.state.V * (P - s.P)
        else:
            ideal_work = pkg.flash(s.zs, P=P, S=s.state.S).H - s.state.H
        work = ideal_work / self.p["efficiency"]
        out = Stream.from_flash(pkg, s.flows, P=P, H=s.state.H + work)
        self.results["power"] = work * s.total_flow
        return {"out": out}


class Pump(_Machine):
    type = "Pump"
    label = "Pump"
    category = "Pressure Changers"
    description = "Raises the pressure of a liquid. Power from an isentropic path and efficiency."
    params = (
        P_OUT,
        Param("dP", "Pressure increase", "pressure_difference", min=0.0),
        Param("efficiency", "Efficiency", "fraction", default=0.75, min=0.01, max=1.0),
    )


class Compressor(_Machine):
    type = "Compressor"
    label = "Compressor"
    category = "Pressure Changers"
    description = "Compresses a gas. Isentropic efficiency sets power and outlet temperature."
    params = (
        P_OUT,
        Param("dP", "Pressure increase", "pressure_difference", min=0.0),
        Param("efficiency", "Isentropic efficiency", "fraction", default=0.72, min=0.01, max=1.0),
    )
