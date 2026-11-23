"""Parameters and results for heating and cooling with a utility; the model is in
``fluxsheet.utilities``. Units with an energy balance add these to their own.

U is either given or, opt-in with ``U_from``, found from the geometry and film coefficients
(``fluxsheet.film_coefficients``). The geometry parameters are shown once ``U_from`` is set.
"""

from ..film_coefficients import TUBE_CORRELATIONS
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
          "water-cooled jacket on a stirred liquid, more for condensing steam. Or leave it "
          "empty and find it from the geometry (U from film coefficients).")
AREA = Param("area", "Heat transfer area", "area", min=0.0, tab="Jacket / coil",
             help="Area of the jacket or coil. The duty then follows from U, the area, the "
             "utility and the unit's own temperature.")
ARRANGEMENT = Param("arrangement", "Flow arrangement", "choice", default="countercurrent",
                    options=("countercurrent", "cocurrent"), requires="utility",
                    help="Direction of the utility relative to the process stream.")


# ---- U from film coefficients -------------------------------------------------------------

def u_from(options: tuple[str, ...], help: str) -> Param:
    return Param("U_from", "U from film coefficients", "choice", options=options,
                 help=help + " Instead of a given U: 1/U = 1/h_i + R_f,i + wall + R_f,o + "
                 "1/h_o, with the film coefficients from correlations and the fluid "
                 "properties from the property package.")


def _film(name: str, label: str, quantity: str | None, help: str, **kw) -> Param:
    return Param(name, label, quantity, help=help, requires="U_from", **kw)


VESSEL_U_FROM = u_from(
    ("jacket", "coil"),
    "Stirred vessel with a jacket or a helical coil (Chilton-Drew-Jebens type on the "
    "process side). U is referred to the vessel wall (jacket) or the outer coil surface.")
VESSEL_DIAMETER = _film("vessel_diameter", "Vessel diameter", "length",
                        "Inner diameter of the vessel, the length of its Nusselt number.",
                        min=0.0)
IMPELLER_DIAMETER = _film("impeller_diameter", "Impeller diameter", "length",
                          "Impeller diameter d in Re = ρ N d² / μ.", min=0.0)
IMPELLER_SPEED = _film("impeller_speed", "Impeller speed", "rotational_speed",
                       "Revolutions per time.", min=0.0)
IMPELLER_CONSTANT = _film(
    "impeller_constant", "Impeller constant C", None,
    "C in Nu = C Re^a Pr^(1/3) (μ/μw)^0.14. Defaults to the Chilton-Drew-Jebens paddle "
    "values: 0.36 at the jacket, 0.87 on a coil.", min=0.0)
IMPELLER_EXPONENT = _film(
    "impeller_exponent", "Impeller exponent a", None,
    "Exponent of Re. Defaults to 2/3 at the jacket, 0.62 on a coil.", min=0.0)
PROCESS_FILM = _film(
    "process_film", "Process side", "choice",
    "stirred: the agitated contents (Chilton-Drew-Jebens). boiling: the contents boil on "
    "the heated wall, nucleate pool boiling (Mostinski; simplified, no critical heat flux "
    "check).", default="stirred", options=("stirred", "boiling"))
JACKET_GAP = _film(
    "jacket_gap", "Jacket gap", "length",
    "Width of an annular jacket a liquid utility flows through (as an annulus). Not needed "
    "for steam, or with a given medium-side film coefficient.", min=0.0)
COIL_DIAMETER = _film(
    "coil_diameter", "Coil diameter", "length",
    "Diameter of the coil's helix, for the enhancement 1 + 3.5 d_i/D inside it. Without "
    "it the coil is rated as a straight tube.", min=0.0)
TUBE_OD = _film("tube_outer_diameter", "Tube outer diameter", "length",
                "Outer diameter of the tubes (of the coil).", min=0.0)
WALL_THICKNESS = _film("wall_thickness", "Wall thickness", "length",
                       "Tube or vessel wall. Tubes: the inner diameter is the outer one "
                       "minus twice this. A vessel without it has no wall resistance.",
                       min=0.0)
WALL_CONDUCTIVITY = _film("wall_conductivity", "Wall thermal conductivity",
                          "thermal_conductivity",
                          "Defaults to 16 W/m·K (stainless steel); carbon steel is about "
                          "50.", default=16.0, min=0.0)
FOULING_PROCESS = _film("fouling_process", "Fouling, process side", "fouling_resistance",
                        "Fouling resistance on the process side's surface, e.g. 0.0002 "
                        "m²K/W. A glass lining can be entered here as its s/λ.", min=0.0)
FOULING_MEDIUM = _film("fouling_medium", "Fouling, medium side", "fouling_resistance",
                       "Fouling resistance on the medium side's surface, e.g. 0.0002 m²K/W "
                       "for treated cooling water.", min=0.0)
H_MEDIUM = _film(
    "h_medium", "Film coefficient, medium side", "heat_transfer_coefficient",
    "Given film coefficient of the utility side. Without it: condensing steam by Nusselt "
    "film condensation; a liquid utility is rated as water flowing through the jacket gap, "
    "the coil, the annulus or the shell.", min=0.0)
TUBE_CORRELATION = _film("tube_correlation", "Tube correlation", "choice",
                         "Turbulent flow in tubes and annuli. Laminar flow (Re < 2300) has "
                         "Nu = 3.66, the transition up to 10⁴ is interpolated.",
                         default="gnielinski", options=TUBE_CORRELATIONS)

#: Parameters of a stirred vessel's film coefficients.
VESSEL_FILM_PARAMS = (VESSEL_U_FROM, VESSEL_DIAMETER, IMPELLER_DIAMETER, IMPELLER_SPEED,
                      IMPELLER_CONSTANT, IMPELLER_EXPONENT, PROCESS_FILM, JACKET_GAP,
                      COIL_DIAMETER, TUBE_OD, WALL_THICKNESS, WALL_CONDUCTIVITY,
                      FOULING_PROCESS, FOULING_MEDIUM, H_MEDIUM, TUBE_CORRELATION)

# Tube bundles and double pipes (heaters, two-sided exchangers).
TUBE_U_FROM_HELP = ("double_pipe: one or more tubes inside an outer pipe, the outer fluid "
                    "in the annulus. shell_and_tube: a tube bundle in a baffled shell, "
                    "the shell side by Kern's method. U is referred to the outer tube "
                    "surface, so the area is tubes × π d_o × length.")
TUBES = _film("tubes", "Number of tubes", "integer", "Tubes in parallel (all passes).",
              default=1, min=1)
TUBE_PASSES = _film("tube_passes", "Tube passes", "integer",
                    "Passes of a shell-and-tube bundle; the tubes of one pass carry the "
                    "whole tube-side flow. The LMTD is not corrected for multipass flow.",
                    default=1, min=1)
ANNULUS_DIAMETER = _film("annulus_diameter", "Outer pipe inner diameter", "length",
                         "Double pipe: inner diameter of the outer pipe.", min=0.0)
SHELL_DIAMETER = _film("shell_diameter", "Shell inner diameter", "length",
                       "Shell and tube.", min=0.0)
TUBE_PITCH = _film("tube_pitch", "Tube pitch", "length",
                   "Centre-to-centre distance of the tubes.", min=0.0)
PITCH_LAYOUT = _film("pitch_layout", "Tube layout", "choice", "Tube pitch layout.",
                     default="triangular", options=("triangular", "square"))
BAFFLE_SPACING = _film("baffle_spacing", "Baffle spacing", "length",
                       "Distance between the shell's baffles.", min=0.0)
TUBE_GEOMETRY = (TUBES, TUBE_PASSES, TUBE_OD, WALL_THICKNESS, WALL_CONDUCTIVITY,
                 ANNULUS_DIAMETER, SHELL_DIAMETER, TUBE_PITCH, PITCH_LAYOUT, BAFFLE_SPACING,
                 TUBE_CORRELATION)
#: A heater: the process stream in the tubes, the utility outside.
HEATER_FILM_PARAMS = (u_from(("double_pipe", "shell_and_tube"), "The process stream flows "
                             "in the tubes, the utility outside. " + TUBE_U_FROM_HELP),
                      *TUBE_GEOMETRY, FOULING_PROCESS, FOULING_MEDIUM, H_MEDIUM)

#: Parameters of a well-mixed unit; the area goes into its energy specification group.
HEAT_PARAMS = (UTILITY, UTILITY_FLOW, U, *VESSEL_FILM_PARAMS)

FILM_RESULTS = (
    ResultField("U", "U from film coefficients", "heat_transfer_coefficient"),
    ResultField("h_process", "Film coefficient, process side", "heat_transfer_coefficient"),
    ResultField("h_medium", "Film coefficient, medium side", "heat_transfer_coefficient"),
    ResultField("Re_process", "Reynolds number, process side"),
    ResultField("Nu_process", "Nusselt number, process side"),
    ResultField("Re_medium", "Reynolds number, medium side"),
    ResultField("Nu_medium", "Nusselt number, medium side"),
    ResultField("T_wall", "Wall temperature, process side", "temperature"),
)

HEAT_RESULTS = (
    ResultField("utility_flow", "Medium flow", "mass_flow"),
    ResultField("utility_T_out", "Medium outlet temperature", "temperature"),
    ResultField("LMTD", "Mean temperature difference", "temperature_difference"),
    ResultField("UA", "UA", "thermal_conductance"),
    ResultField("heat_transfer_area", "Heat transfer area", "area"),
    *FILM_RESULTS,
)
