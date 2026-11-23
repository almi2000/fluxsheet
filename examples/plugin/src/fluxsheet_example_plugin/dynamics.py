"""A dynamic model for the dew-point cooler: its outlet lags behind the steady state.

Registered through the ``fluxsheet.dynamic_models`` entry point as {unit type: model}. In a
dynamic run the cooler's outlet enthalpy h follows the steady-state outlet h_ss (the unit's
own ``calculate`` at the inlets and parameters of the moment) with a first-order time
constant tau, the thermal mass of the exchanger:

    dh/dt = (h_ss - h) / tau

so after a step of the superheat the duty approaches its new value as 1 - exp(-t / tau).
"""

import numpy as np

from fluxsheet.contents_flash import ContentsFlash
from fluxsheet.dynamics.units import DynamicModel
from fluxsheet.stream import Stream
from fluxsheet.units.base import Param, ResultField


class CoolerLag(DynamicModel):
    params = (
        Param("time_constant", "Time constant", "time", default=30.0, min=1e-3,
              help="Thermal time constant of the cooler: the outlet covers 63 % of a step "
              "in this time."),
    )
    result_fields = (ResultField("T_out", "Outlet temperature", "temperature"),
                     ResultField("duty", "Duty", "power"))
    # The outlet depends on the inlet at every instant (no hold-up of material), so the
    # simulator gives the model ``self.steady(inlets)``: the unit solved at those inlets.
    source = False

    def initialize(self, inlets, outlets):
        out = outlets["out"]
        self.flash = ContentsFlash(self.pkg, out.T)  # PH flashes along the trajectory
        return np.array([out.state.H])

    def scale(self, x0):
        return np.array([max(abs(x0[0]), 8.314 * 300.0)])

    def outlets(self, x, inlets=None):
        target = self.steady(inlets)["out"]
        self._target = target.state.H
        h_in = sum(s.enthalpy_flow for streams in inlets.values() for s in streams)
        if target.total_flow <= 0.0:
            self.results.update(T_out=target.T, duty=0.0)
            return {"out": target}
        st = self.flash(np.asarray(target.zs), target.P, float(x[0]))
        out = Stream(list(target.flows), st)
        self.results.update(T_out=st.T, duty=out.enthalpy_flow - h_in)
        return {"out": out}

    def derivatives(self, x, inlets):
        return np.array([(self._target - x[0]) / self.s["time_constant"]])


MODELS = {"DewPointCooler": CoolerLag}
