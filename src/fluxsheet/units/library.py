"""The built-in unit operations."""

from __future__ import annotations

from dataclasses import replace

from ..errors import SimulationError
from ..stream import Stream
from ..thermo_package import PropertyPackage
from .base import Param, Port, UnitOperation, empty_stream

_IN = Port("in", "Inlet")
_OUT = Port("out", "Outlet")

_P_OUT = Param("P", "Outlet pressure", "pressure", help="Absolute outlet pressure.", min=1.0)
_DP = Param("dP", "Pressure drop", "pressure_difference",
            help="Pressure loss across the unit. Defaults to zero.", min=0.0)


def _single(inlets: dict[str, list[Stream]], port: str = "in") -> Stream:
    return inlets[port][0]


class Feed(UnitOperation):
    type = "Feed"
    label = "Feed"
    category = "Streams"
    description = "A material stream entering the process, at a given state and flow rate."
    outlets = (_OUT,)
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
    inlets = (_IN,)

    def calculate(self, pkg, inlets):
        return {}


class Mixer(UnitOperation):
    type = "Mixer"
    label = "Mixer"
    category = "Mixers & Splitters"
    description = "Combines any number of streams adiabatically. Outlet pressure defaults to the "\
        "lowest inlet pressure."
    inlets = (Port("in", "Inlets", multiple=True),)
    outlets = (_OUT,)
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
    inlets = (_IN,)
    outlets = (Port("out1", "Outlet 1"), Port("out2", "Outlet 2"))
    params = (Param("split", "Fraction to outlet 1", "fraction", default=0.5, min=0.0, max=1.0),)

    def calculate(self, pkg, inlets):
        s = _single(inlets)
        f = self.p["split"]
        return {
            "out1": Stream([x * f for x in s.flows], s.state),
            "out2": Stream([x * (1.0 - f) for x in s.flows], s.state),
        }


class ComponentSplitter(UnitOperation):
    type = "ComponentSplitter"
    label = "Component splitter"
    category = "Separators"
    description = "Idealized separation: sends a chosen fraction of each component to the top. "\
        "Useful for conceptual design before a rigorous column is set up."
    inlets = (_IN,)
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
        s = _single(inlets)
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


class Heater(UnitOperation):
    type = "Heater"
    label = "Heater / Cooler"
    category = "Heat Exchangers"
    description = "Heats or cools a stream to a target temperature, vapor fraction, or with a "\
        "given duty."
    inlets = (_IN,)
    outlets = (_OUT,)
    params = (
        Param("T", "Outlet temperature", "temperature", min=0.0),
        Param("VF", "Outlet vapor fraction", "fraction", min=0.0, max=1.0),
        Param("duty", "Duty", "power", help="Positive = heating, negative = cooling."),
        _P_OUT,
        _DP,
    )
    exclusive_groups = (("T", "VF", "duty"), ("P", "dP"))
    required_groups = (("T", "VF", "duty"),)

    def calculate(self, pkg, inlets):
        s = _single(inlets)
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


class Flash(UnitOperation):
    type = "Flash"
    label = "Flash drum"
    category = "Separators"
    description = "Separates a stream into vapor and liquid in equilibrium. Adiabatic unless a "\
        "temperature or duty is given."
    inlets = (_IN,)
    outlets = (Port("vapor", "Vapor"), Port("liquid", "Liquid"))
    params = (
        Param("T", "Temperature", "temperature", min=0.0),
        Param("duty", "Duty", "power", help="Heat added. Defaults to zero (adiabatic)."),
        _P_OUT,
        _DP,
    )
    exclusive_groups = (("T", "duty"), ("P", "dP"))

    def calculate(self, pkg, inlets):
        s = _single(inlets)
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


class Valve(UnitOperation):
    type = "Valve"
    label = "Valve"
    category = "Pressure Changers"
    description = "Isenthalpic pressure reduction (throttling). The outlet may partially "\
        "vaporize and cool."
    inlets = (_IN,)
    outlets = (_OUT,)
    params = (_P_OUT, Param("dP", "Pressure drop", "pressure_difference", min=0.0))
    exclusive_groups = (("P", "dP"),)
    required_groups = exclusive_groups

    def calculate(self, pkg, inlets):
        s = _single(inlets)
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
    inlets = (_IN,)
    outlets = (_OUT,)
    exclusive_groups = (("P", "dP"),)
    required_groups = exclusive_groups
    raises_pressure = True

    def calculate(self, pkg: PropertyPackage, inlets):
        s = _single(inlets)
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
        _P_OUT,
        Param("dP", "Pressure increase", "pressure_difference", min=0.0),
        Param("efficiency", "Efficiency", "fraction", default=0.75, min=0.01, max=1.0),
    )


class Compressor(_Machine):
    type = "Compressor"
    label = "Compressor"
    category = "Pressure Changers"
    description = "Compresses a gas. Isentropic efficiency sets power and outlet temperature."
    params = (
        _P_OUT,
        Param("dP", "Pressure increase", "pressure_difference", min=0.0),
        Param("efficiency", "Isentropic efficiency", "fraction", default=0.72, min=0.01, max=1.0),
    )


UNIT_TYPES: dict[str, type[UnitOperation]] = {
    cls.type: cls
    for cls in (Feed, Product, Mixer, Splitter, Heater, Flash, ComponentSplitter, Valve, Pump,
                Compressor)
}
