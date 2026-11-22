"""A unit operation: cool (or heat) a vapor to its dew point plus a margin."""

from fluxsheet.stream import Stream
from fluxsheet.units.base import (
    DUTY,
    IN,
    OUT,
    P_OUT,
    Param,
    UnitOperation,
    empty_stream,
    single_inlet,
)


class DewPointCooler(UnitOperation):
    type = "DewPointCooler"
    label = "Dew-point cooler"
    category = "Heat Exchangers"
    description = ("Cools a vapor to its dew point at the outlet pressure, plus a superheat "
                   "margin, so nothing condenses yet. Example plug-in unit.")
    inlets = (IN,)
    outlets = (OUT,)
    params = (
        Param("superheat", "Superheat above the dew point", "temperature_difference",
              default=0.0, min=0.0),
        P_OUT,
    )
    result_fields = (DUTY,)
    energy_results = ("duty",)  # the duty is the only energy crossing the boundary
    process_side = "plug"

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P = self.p["P"] or s.P
        if s.total_flow <= 0.0:
            self.results["duty"] = 0.0
            return {"out": empty_stream(pkg, s.T, P)}
        zs = [f / s.total_flow for f in s.flows]
        T = pkg.saturation_point(zs, "dew", P=P).T + self.p["superheat"]
        out = Stream.from_flash(pkg, s.flows, T=T, P=P)
        self.results["duty"] = out.enthalpy_flow - s.enthalpy_flow
        return {"out": out}
