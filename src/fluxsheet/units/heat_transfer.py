"""Parameters and results for heating and cooling with a utility; the model is in
``fluxsheet.utilities``. Units with an energy balance add these to their own."""

from .base import Param, ResultField

UTILITY = Param("utility", "Heating / cooling medium", "utility",
                help="A utility from the flowsheet's list. The unit reports how much of it "
                "it uses; with U and an area, the heat exchanged through a jacket or coil "
                "sets the duty.")
UTILITY_FLOW = Param("utility_flow", "Medium flow", "mass_flow", min=0.0, requires="utility",
                     help="Fluid utilities: fixes the medium flow, so its outlet temperature "
                     "follows. Defaults to the flow that brings it to the utility's outlet "
                     "temperature.")
U = Param("U", "Heat transfer coefficient U", "heat_transfer_coefficient", min=0.0,
          requires="utility",
          help="Overall coefficient through the wall, e.g. a few hundred W/m²K for a "
          "water-cooled jacket on a stirred liquid, more for condensing steam.")
AREA = Param("area", "Heat transfer area", "area", min=0.0, tab="Jacket / coil",
             help="Area of the jacket or coil. The duty then follows from U, the area, the "
             "utility and the unit's own temperature.")
ARRANGEMENT = Param("arrangement", "Flow arrangement", "choice", default="countercurrent",
                    options=("countercurrent", "cocurrent"), requires="utility",
                    help="Direction of the utility relative to the process stream.")
#: Parameters of a well-mixed unit; the area goes into its energy specification group.
HEAT_PARAMS = (UTILITY, UTILITY_FLOW, U)

HEAT_RESULTS = (
    ResultField("utility_flow", "Medium flow", "mass_flow"),
    ResultField("utility_T_out", "Medium outlet temperature", "temperature"),
    ResultField("LMTD", "Mean temperature difference", "temperature_difference"),
    ResultField("UA", "UA", "thermal_conductance"),
    ResultField("heat_transfer_area", "Heat transfer area", "area"),
)
