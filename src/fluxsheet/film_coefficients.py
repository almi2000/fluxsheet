"""Film coefficients from correlations, and the overall heat transfer coefficient from them.

The overall coefficient U is a series of resistances, all referred to one reference area
A_ref (which one is up to the unit, see below)::

    1/U = (A_ref/A_i) (1/h_i + R_f,i) + R_wall + (A_ref/A_o) (R_f,o + 1/h_o)

with h the film coefficients of the two sides, R_f their fouling resistances (m²K/W on
their own surface) and R_wall the wall's conduction resistance referred to A_ref: for a
tube wall d_ref ln(d_o/d_i) / (2 k_w), for a plane wall (a vessel jacket) s / k_w. The
units refer U to these areas:

- a heater, a two-sided exchanger, a coil: the outer tube surface, π d_o L per tube;
- a vessel jacket: the vessel wall (plane, inner and outer area taken as equal);
- a tubular reactor: the inner tube surface π D L, the area it already uses.

Correlations (all SI, dimensionless numbers on the stated length):

- **Inside tubes**: fully developed laminar flow Nu = 3.66 (constant wall temperature) or
  4.36 (constant heat flux) up to Re = 2300; Gnielinski (1976),
  Nu = (f/8)(Re - 1000) Pr / (1 + 12.7 sqrt(f/8)(Pr^(2/3) - 1)), f = (0.790 ln Re -
  1.64)^-2, from Re = 10^4; in between linear in Re from the laminar value to Gnielinski
  at 10^4, as the VDI Heat Atlas (G1) recommends, so Nu is continuous. Dittus-Boelter,
  Nu = 0.023 Re^0.8 Pr^n (n = 0.4 heating the fluid, 0.3 cooling it), replaces Gnielinski
  as an option, with the same laminar and transition treatment. No entrance effects and no
  property-variation correction: a long tube, conservative.
- **Annulus** (double pipe, heat through the inner tube): on the hydraulic diameter
  D_o - d_i; laminar Nu = 3.66 + 1.2 (d_i/D_o)^-0.8 (VDI G2, outer wall insulated);
  turbulent the tube correlation times 0.86 (d_i/D_o)^-0.16 (Petukhov and Roizen).
- **Helical coil, inside**: the tube correlation times 1 + 3.5 d_i/D_coil (Jeschke).
- **Shell side, Kern (1950)**: Nu = 0.36 Re^0.55 Pr^(1/3) (mu/mu_w)^0.14 on Kern's
  equivalent diameter, Re = G_s D_e / mu with the cross-flow area between two baffles
  A_s = (p_t - d_o) D_s B / p_t; meant for 2000 < Re < 10^6, a segmental-baffle shell
  without leakage or bypass corrections (Bell-Delaware would be lower).
- **Stirred vessel** (Chilton, Drew and Jebens, 1944): Nu = C Re^a Pr^(1/3)
  (mu/mu_w)^0.14 with Nu = h D_T / k on the vessel diameter, Re = rho N d² / mu on the
  impeller diameter d and speed N (1/s): C = 0.36, a = 2/3 at the jacket wall and C =
  0.87, a = 0.62 on a coil (paddle stirrer, the original constants); other impellers have
  other constants.
- **Film condensation** (Nusselt, laminar film): h = C [g rho_l (rho_l - rho_v) k_l³ h'_fg
  / (mu_l dT L)]^(1/4), C = 0.943 on a vertical wall of height L and 0.729 outside a
  horizontal tube of diameter L, h'_fg = h_fg + 0.68 cp_l dT (Rohsenow), liquid
  properties at the film temperature. A wavy or turbulent film transfers more; this is the
  conservative laminar value. *Simplified.*
- **Nucleate boiling** (Mostinski, 1963): h = 0.104 P_c^0.69 q^0.7 (1.8 p_r^0.17 + 4
  p_r^1.2 + 10 p_r^10), P_c in bar, q in W/m², reduced pressure p_r = P / P_c; a mixture
  takes its mole-averaged critical pressure. Fully developed nucleate pool boiling only, no
  critical heat flux check. *Simplified.*

Film coefficients that depend on the wall (the viscosity at the wall, the temperature drop
of a condensing film, the heat flux of a boiling one) are found together with the wall
temperatures by a damped fixed-point iteration (:func:`overall`).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

from .errors import SimulationError

G = 9.80665  # m/s²
RE_LAMINAR = 2300.0
RE_TURBULENT = 1e4
NU_LAMINAR = {"T": 3.66, "q": 4.36}  # fully developed laminar flow in a round tube
TUBE_CORRELATIONS = ("gnielinski", "dittus-boelter")
#: Chilton-Drew-Jebens constants (C, exponent of Re) for a paddle stirrer.
STIRRED = {"jacket": (0.36, 2.0 / 3.0), "coil": (0.87, 0.62)}
KERN_RE_RANGE = (2e3, 1e6)


# ---- inside tubes ---------------------------------------------------------------------------

def friction_factor(Re: float) -> float:
    """Darcy friction factor of a smooth tube, Filonenko / Petukhov: (0.790 ln Re - 1.64)^-2."""
    return (0.790 * math.log(Re) - 1.64) ** -2


def gnielinski(Re: float, Pr: float) -> float:
    """Gnielinski's turbulent Nusselt number (fully developed, no wall correction)."""
    f8 = friction_factor(Re) / 8.0
    return f8 * (Re - 1000.0) * Pr / (1.0 + 12.7 * math.sqrt(f8) * (Pr ** (2.0 / 3.0) - 1.0))


def dittus_boelter(Re: float, Pr: float, heating: bool = True) -> float:
    """Dittus-Boelter: 0.023 Re^0.8 Pr^n, n = 0.4 if the fluid is heated, 0.3 if cooled."""
    return 0.023 * Re ** 0.8 * Pr ** (0.4 if heating else 0.3)


def nu_tube(Re: float, Pr: float, correlation: str = "gnielinski", heating: bool = True,
            boundary: str = "T", laminar: float | None = None) -> tuple[float, str]:
    """Mean Nusselt number in a tube (or on a hydraulic diameter) and the flow regime.

    ``boundary`` picks the laminar value: "T" constant wall temperature (3.66), "q"
    constant heat flux (4.36); ``laminar`` overrides it (e.g. for an annulus).
    """
    if correlation not in TUBE_CORRELATIONS:
        raise ValueError(f"unknown tube correlation {correlation!r}")
    lam = NU_LAMINAR[boundary] if laminar is None else laminar

    def turbulent(re: float) -> float:
        return gnielinski(re, Pr) if correlation == "gnielinski" else dittus_boelter(
            re, Pr, heating)

    if Re <= RE_LAMINAR:
        return lam, "laminar"
    if Re >= RE_TURBULENT:
        return turbulent(Re), "turbulent"
    g = (Re - RE_LAMINAR) / (RE_TURBULENT - RE_LAMINAR)
    return (1.0 - g) * lam + g * turbulent(RE_TURBULENT), "transition"


def nu_annulus(Re: float, Pr: float, a: float, correlation: str = "gnielinski",
               heating: bool = True) -> tuple[float, str]:
    """Annulus, heat through the inner wall, ``a`` = d_i / D_o; Re on D_o - d_i."""
    if not 0.0 < a < 1.0:
        raise ValueError("the annulus needs 0 < d_i/D_o < 1")
    laminar = 3.66 + 1.2 * a ** -0.8
    factor = 0.86 * a ** -0.16
    if Re <= RE_LAMINAR:
        return laminar, "laminar"

    def turbulent(re: float) -> float:
        return factor * nu_tube(re, Pr, correlation, heating)[0]

    if Re >= RE_TURBULENT:
        return turbulent(Re), "turbulent"
    g = (Re - RE_LAMINAR) / (RE_TURBULENT - RE_LAMINAR)
    return (1.0 - g) * laminar + g * turbulent(RE_TURBULENT), "transition"


def coil_factor(d_i: float, D_coil: float) -> float:
    """Enhancement inside a helical coil, 1 + 3.5 d_i / D_coil (Jeschke)."""
    return 1.0 + 3.5 * d_i / D_coil


# ---- shell side (Kern) ----------------------------------------------------------------------

def kern_equivalent_diameter(d_o: float, pitch: float, layout: str = "triangular") -> float:
    """Kern's equivalent diameter: four times the free area around a tube over its wetted
    perimeter (for triangular pitch per half tube, with Kern's 0.86 p_t row height)."""
    if pitch <= d_o:
        raise ValueError("the tube pitch must be larger than the tube outer diameter")
    if layout == "square":
        return 4.0 * (pitch ** 2 - math.pi * d_o ** 2 / 4.0) / (math.pi * d_o)
    if layout == "triangular":
        return 4.0 * (0.5 * pitch * 0.86 * pitch - 0.5 * math.pi * d_o ** 2 / 4.0) / (
            0.5 * math.pi * d_o)
    raise ValueError(f"unknown tube layout {layout!r}")


def kern_crossflow_area(D_s: float, d_o: float, pitch: float, baffle_spacing: float) -> float:
    """Shell-side cross-flow area at the shell centre line: (p_t - d_o) D_s B / p_t."""
    return (pitch - d_o) * D_s * baffle_spacing / pitch


def nu_kern(Re: float, Pr: float, mu_ratio: float = 1.0) -> float:
    """Kern's shell-side Nusselt number on the equivalent diameter."""
    return 0.36 * Re ** 0.55 * Pr ** (1.0 / 3.0) * mu_ratio ** 0.14


# ---- stirred vessels ------------------------------------------------------------------------

def stirred_reynolds(rho: float, N: float, d: float, mu: float) -> float:
    """Impeller Reynolds number rho N d² / mu (N in revolutions per second)."""
    return rho * N * d ** 2 / mu


def nu_stirred(Re: float, Pr: float, C: float, a: float, mu_ratio: float = 1.0) -> float:
    """Chilton-Drew-Jebens type: C Re^a Pr^(1/3) (mu/mu_w)^0.14, Nu on the vessel diameter."""
    return C * Re ** a * Pr ** (1.0 / 3.0) * mu_ratio ** 0.14


# ---- condensation and boiling (simplified) -------------------------------------------------

def h_condensation(rho_l: float, rho_v: float, mu_l: float, k_l: float, cp_l: float,
                   h_fg: float, dT: float, length: float, geometry: str = "vertical") -> float:
    """Mean Nusselt film condensation coefficient, W/(m² K), for a wall ``dT`` below the
    saturation temperature: a vertical wall of height ``length`` or the outside of a
    horizontal tube of diameter ``length``."""
    C = {"vertical": 0.943, "horizontal_tube": 0.729}[geometry]
    dT = max(dT, 1e-6)
    h_fg_mod = h_fg + 0.68 * cp_l * dT
    return C * (G * rho_l * (rho_l - rho_v) * k_l ** 3 * h_fg_mod
                / (mu_l * dT * length)) ** 0.25


def h_nucleate_boiling(q: float, P: float, Pc: float) -> float:
    """Mostinski's nucleate pool boiling coefficient, W/(m² K), at heat flux ``q`` (W/m²),
    pressure ``P`` and critical pressure ``Pc`` (Pa)."""
    pr = P / Pc
    F = 1.8 * pr ** 0.17 + 4.0 * pr ** 1.2 + 10.0 * pr ** 10
    return 0.104 * (Pc / 1e5) ** 0.69 * max(q, 1e-9) ** 0.7 * F


# ---- walls ----------------------------------------------------------------------------------

def tube_wall_resistance(d_i: float, d_o: float, k_w: float, d_ref: float) -> float:
    """Conduction resistance of a tube wall, m²K/W on the surface of diameter ``d_ref``."""
    return d_ref * math.log(d_o / d_i) / (2.0 * k_w)


# ---- overall coefficient --------------------------------------------------------------------

@dataclass(frozen=True)
class Film:
    """One side's film coefficient and how it was found."""

    h: float  # W/(m² K), on this side's own surface
    correlation: str
    Re: float | None = None
    Nu: float | None = None
    regime: str | None = None
    given: bool = False  # entered, not from a correlation


@dataclass(frozen=True)
class Side:
    """One side of the wall.

    ``film(T_wall, q)`` gives the film at a wall temperature on this side and a heat flux
    (W/m², magnitude) through this side's own surface; ``wall_dependent`` says whether it
    uses them. ``area_ratio`` is A_ref over this side's area, ``fouling`` its fouling
    resistance on its own area.
    """

    T: float  # bulk temperature, K
    film: Callable[[float, float], Film]
    area_ratio: float = 1.0
    fouling: float = 0.0
    wall_dependent: bool = False


@dataclass(frozen=True)
class Overall:
    U: float  # W/(m² K) on the reference area
    inner: Film
    outer: Film
    T_wall_inner: float
    T_wall_outer: float


def given_film(h: float) -> Callable[[float, float], Film]:
    film = Film(h, "given", given=True)
    return lambda T_wall, q: film


def overall(inner: Side, outer: Side, wall: float = 0.0, max_iter: int = 200) -> Overall:
    """U from the two sides and the wall resistance ``wall`` (m²K/W on the reference area).

    Heat flows from the hotter bulk to the colder one. Where a film depends on its wall
    temperature or on the heat flux, those are iterated (damped by half) until U changes
    by less than 1e-10.
    """
    dT = outer.T - inner.T
    sign = 1.0 if dT >= 0.0 else -1.0
    Tw_i, Tw_o = inner.T + 0.5 * dT, outer.T - 0.5 * dT
    q = None  # heat flux on the reference area
    U_old = None
    iterate = inner.wall_dependent or outer.wall_dependent
    for _ in range(max_iter if iterate else 1):
        q_ref = q if q is not None else max(abs(dT), 1e-3) * 1e3  # a start for a boiling film
        fi = inner.film(Tw_i, q_ref * inner.area_ratio)
        fo = outer.film(Tw_o, q_ref * outer.area_ratio)
        if not (fi.h > 0.0 and fo.h > 0.0):
            raise SimulationError("a film coefficient is not positive",
                                  hint="check the geometry and the flows")
        r_i = inner.area_ratio * (1.0 / fi.h + inner.fouling)
        r_o = outer.area_ratio * (1.0 / fo.h + outer.fouling)
        U = 1.0 / (r_i + wall + r_o)
        q_new = U * abs(dT)
        Tw_i_new = inner.T + sign * q_new * inner.area_ratio / fi.h
        Tw_o_new = outer.T - sign * q_new * outer.area_ratio / fo.h
        if not iterate:
            return Overall(U, fi, fo, Tw_i_new, Tw_o_new)
        if U_old is not None and abs(U - U_old) <= 1e-10 * U:
            return Overall(U, fi, fo, Tw_i_new, Tw_o_new)
        U_old = U
        q = q_new if q is None else 0.5 * (q + q_new)
        Tw_i, Tw_o = 0.5 * (Tw_i + Tw_i_new), 0.5 * (Tw_o + Tw_o_new)
    raise SimulationError("the wall temperatures of the film coefficients did not converge",
                          hint="give the heat transfer coefficient U instead")
