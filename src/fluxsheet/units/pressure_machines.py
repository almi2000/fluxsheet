"""Turbines and multistage compressors (pumps, single compressors and valves: pressure.py)."""

from __future__ import annotations

import numpy as np

from ..errors import SimulationError
from ..stream import Stream
from ..thermo_package import PropertyPackage
from .base import (
    IN,
    OUT,
    P_OUT,
    POWER,
    Param,
    Port,
    ResultField,
    UnitOperation,
    empty_stream,
    single_inlet,
)


def isentropic_work(pkg: PropertyPackage, s: Stream, P: float) -> float:
    """Ideal (isentropic) enthalpy change per mole from the inlet state to pressure P."""
    if s.state.VF < 1e-9:
        # Liquid: incompressible, V dP. Exact enough and avoids an ill-conditioned flash.
        return s.state.V * (P - s.P)
    return pkg.flash(s.zs, P=P, S=s.state.S).H - s.state.H


class Turbine(UnitOperation):
    type = "Turbine"
    label = "Turbine / expander"
    category = "Pressure Changers"
    description = "Expands a gas or steam (or a liquid, as a hydraulic turbine) to a lower "\
        "pressure and recovers work. Isentropic efficiency sets the power and the outlet "\
        "state, which may be wet. The power is negative: it is produced."
    inlets = (IN,)
    outlets = (OUT,)
    params = (
        P_OUT,
        Param("dP", "Pressure drop", "pressure_difference", min=0.0),
        Param("efficiency", "Isentropic efficiency", "fraction", default=0.8, min=0.01, max=1.0),
    )
    exclusive_groups = (("P", "dP"),)
    required_groups = exclusive_groups
    result_fields = (POWER, ResultField("vapor_fraction", "Outlet vapor fraction"))

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P = self.outlet_pressure(s.P)
        if P >= s.P:
            raise SimulationError(
                f"outlet pressure ({P / 1e5:.3g} bar) is not below the inlet pressure "
                f"({s.P / 1e5:.3g} bar)", where=self.id,
                hint="a turbine expands; use a compressor or pump to raise the pressure")
        if s.total_flow <= 0.0:
            self.results.update(power=0.0)
            return {"out": empty_stream(pkg, s.T, P)}
        work = isentropic_work(pkg, s, P) * self.p["efficiency"]
        out = Stream.from_flash(pkg, s.flows, P=P, H=s.state.H + work)
        self.results.update(power=work * s.total_flow, vapor_fraction=out.state.VF)
        return {"out": out}


class MultistageCompressor(UnitOperation):
    type = "MultistageCompressor"
    label = "Multistage compressor"
    category = "Pressure Changers"
    description = "Compresses a gas in stages of equal pressure ratio with intercoolers "\
        "between them; liquid that condenses in a cooler drains to the condensate outlet. "\
        "Needs less power than one stage and keeps the discharge temperature down."
    inlets = (IN,)
    outlets = (OUT, Port("condensate", "Condensate"))
    params = (
        P_OUT,
        Param("stages", "Stages", "integer", default=3, min=1, max=12),
        Param("efficiency", "Isentropic efficiency", "fraction", default=0.75, min=0.01,
              max=1.0, help="Of every stage."),
        Param("T_cooler", "Intercooler outlet temperature", "temperature", min=0.0,
              help="The gas is cooled to this temperature between stages. Leave empty for "
              "no intercooling."),
        Param("aftercooler", "Aftercooler", "choice", default="no", options=("no", "yes"),
              help="Cool the discharge to the intercooler temperature as well."),
        Param("dP_cooler", "Pressure drop per cooler", "pressure_difference", default=0.0,
              min=0.0),
    )
    exclusive_groups = (("P",),)
    required_groups = exclusive_groups
    result_fields = (
        POWER,
        ResultField("cooler_duty", "Cooler duty", "power"),
        ResultField("pressure_ratio", "Pressure ratio per stage"),
        ResultField("T_max", "Highest discharge temperature", "temperature"),
    )

    def idle_outlets(self) -> tuple[str, ...]:
        return ("condensate",)

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P_out, n = self.p["P"], self.p["stages"]
        if P_out <= s.P:
            raise SimulationError(
                f"outlet pressure ({P_out / 1e5:.3g} bar) is not above the inlet pressure "
                f"({s.P / 1e5:.3g} bar)", where=self.id, hint="use a valve or a turbine")
        if s.total_flow <= 0.0:
            self.results.update(power=0.0, cooler_duty=0.0)
            return {"out": empty_stream(pkg, s.T, P_out),
                    "condensate": empty_stream(pkg, s.T, s.P)}
        if s.state.VF < 1.0 - 1e-9:
            raise SimulationError("the compressor inlet contains liquid", where=self.id,
                                  hint="separate the liquid first, e.g. with a flash drum")
        cooled = [k for k in range(1, n + 1)
                  if self.p["T_cooler"] and (k < n or self.p["aftercooler"] == "yes")]
        dP_c = self.p["dP_cooler"] or 0.0
        # Equal pressure ratios; the stages also make up the coolers' pressure drops.
        ratio = ((P_out + dP_c * len(cooled)) / s.P) ** (1.0 / n)
        gas, P = s, s.P
        drained, H_drained = np.zeros(pkg.n), 0.0
        power = cooler = 0.0
        T_prof, P_prof, w_prof = [], [], []
        for k in range(1, n + 1):
            P *= ratio
            work = isentropic_work(pkg, gas, P) / self.p["efficiency"]
            gas = Stream.from_flash(pkg, gas.flows, P=P, H=gas.state.H + work)
            power += work * gas.total_flow
            T_prof.append(gas.T)
            P_prof.append(P)
            w_prof.append(work * gas.total_flow)
            if k not in cooled:
                continue
            P -= dP_c
            mix = Stream.from_flash(pkg, gas.flows, T=self.p["T_cooler"], P=P)
            cooler += mix.enthalpy_flow - gas.enthalpy_flow
            st = mix.state
            if st.VF >= 1.0 or not st.x:
                gas = mix
                continue
            if st.VF <= 0.0:
                raise SimulationError(
                    f"the gas condenses completely in the cooler after stage {k}",
                    where=self.id, hint="raise the intercooler temperature")
            vapor = mix.total_flow * st.VF * np.array(st.y)
            liquid = np.maximum(np.array(mix.flows) - vapor, 0.0)
            gas = Stream.from_flash(pkg, list(vapor), T=st.T, P=P)
            drained += liquid
            H_drained += mix.enthalpy_flow - gas.enthalpy_flow
        if drained.sum() > 1e-12 * s.total_flow:
            # The drains are collected at the first cooler's pressure.
            P_drain = s.P * ratio - dP_c
            condensate = Stream.from_flash(pkg, list(drained), P=P_drain,
                                           H=H_drained / drained.sum())
        else:
            condensate = empty_stream(pkg, s.T, s.P)
        self.results.update(power=power, cooler_duty=cooler, pressure_ratio=ratio,
                            T_max=max(T_prof))
        self.add_profile("stages", "Stages", "Stage", list(range(1, n + 1)),
                         {"Discharge temperature": ("temperature", T_prof),
                          "Discharge pressure": ("pressure", P_prof),
                          "Power": ("power", w_prof)})
        return {"out": gas, "condensate": condensate}
