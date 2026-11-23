"""Film coefficients, the overall coefficient U from them, and units rated from geometry.

Sources of the reference values:

- Kraume, Transportvorgänge in der Verfahrenstechnik: Aufgaben und Lösungen (Springer
  Vieweg), Example 5.12 (double pipe) and Example 18.2 (glass-lined stirred vessel).
- Kern, Process Heat Transfer (1950), Fig. 28: equivalent diameters of square pitch.
- Incropera et al., Fundamentals of Heat and Mass Transfer: laminar Nu = 3.66 and 4.36
  (Eqs. 8.55, 8.53), Nusselt film condensation constants 0.943 and 0.729 (Eqs. 10.30,
  10.41), Dittus-Boelter (Eq. 8.60).
- Petukhov (1970) and Cooper (1984), independent correlations Gnielinski and Mostinski
  are compared with.
- Chilton, Drew and Jebens, Ind. Eng. Chem. 36 (1944) 510: stirred-vessel constants.

Everything else is checked against the same equations worked by hand here, or against a
closed-form solution (effectiveness-NTU, the exponential approach to condensing steam).
"""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from fluxsheet import film_coefficients as fc
from fluxsheet import solve
from fluxsheet.film_rating import Fluid, annulus_film, tube_film
from fluxsheet.flowsheet import UtilityModel
from fluxsheet.thermo_package import TransportProps, get_package
from fluxsheet.utilities import lmtd, parse_utilities

MW_WATER = 18.01528e-3
STEAM = {"id": "LPS", "type": "steam", "P": "3 bar"}
CW = {"id": "CW", "type": "fluid", "T_in": "25 degC", "T_out": "35 degC", "cp": "4.18 kJ/kgK"}
OIL = {"id": "OIL", "type": "fluid", "T_in": "250 degC", "T_out": "200 degC",
       "cp": "2.3 kJ/kgK"}


# ---- tubes ----------------------------------------------------------------------------------

def test_laminar_values_are_the_textbook_ones():
    # Incropera, Eqs. 8.55 and 8.53: fully developed, constant wall temperature or flux.
    assert fc.nu_tube(1500.0, 5.0)[0] == 3.66
    assert fc.nu_tube(1500.0, 5.0, boundary="q") == (4.36, "laminar")
    assert fc.nu_tube(1500.0, 5.0, "dittus-boelter")[0] == 3.66


def test_gnielinski_by_hand():
    # Re = 1e4, Pr = 0.7: f = (0.790 ln 1e4 - 1.64)^-2 = 5.63614^-2 = 0.0314798,
    # Nu = (f/8)(Re - 1000) Pr / (1 + 12.7 sqrt(f/8) (Pr^(2/3) - 1))
    #    = 0.00393497 * 9000 * 0.7 / (1 + 12.7 * 0.0627293 * (0.788374 - 1)) = 29.817
    assert fc.friction_factor(1e4) == pytest.approx(0.0314798, rel=1e-5)
    assert fc.gnielinski(1e4, 0.7) == pytest.approx(29.817, rel=1e-4)
    # Re = 1e5, Pr = 7: f = (0.790 * 11.5129 - 1.64)^-2 = 0.0179870,
    # Nu = 0.00224837 * 99000 * 7 / (1 + 12.7 * 0.0474170 * 2.65913) = 599.07
    assert fc.gnielinski(1e5, 7.0) == pytest.approx(599.07, rel=1e-4)
    assert fc.nu_tube(1e5, 7.0) == (fc.gnielinski(1e5, 7.0), "turbulent")


def test_transition_runs_linearly_from_laminar_to_gnielinski():
    Pr = 3.0
    turbulent = fc.gnielinski(1e4, Pr)
    assert fc.nu_tube(2300.0, Pr)[0] == pytest.approx(3.66)
    assert fc.nu_tube(1e4, Pr)[0] == pytest.approx(turbulent)
    mid = 0.5 * (2300.0 + 1e4)
    assert fc.nu_tube(mid, Pr) == (pytest.approx(0.5 * (3.66 + turbulent)), "transition")


def petukhov(Re, Pr):
    """Petukhov (1970), the equation Gnielinski extended to lower Re."""
    f8 = (1.82 * math.log10(Re) - 1.64) ** -2 / 8.0
    return f8 * Re * Pr / (1.07 + 12.7 * math.sqrt(f8) * (Pr ** (2 / 3) - 1.0))


@pytest.mark.parametrize("Re", [1e5, 3e5, 1e6])
@pytest.mark.parametrize("Pr", [0.7, 3.0, 20.0, 100.0])
def test_gnielinski_agrees_with_petukhov_at_high_reynolds_numbers(Re, Pr):
    # Gnielinski replaced Re by Re - 1000 and 1.07 by 1 to reach the transition region; in
    # Petukhov's range the two stay within about 8 % (the largest at Pr 0.7).
    assert fc.gnielinski(Re, Pr) == pytest.approx(petukhov(Re, Pr), rel=0.08)


@pytest.mark.parametrize("Re", [1e4, 3e4, 1e5])
@pytest.mark.parametrize("Pr", [0.7, 3.0, 10.0])
def test_dittus_boelter_is_within_its_known_scatter_of_gnielinski(Re, Pr):
    # Incropera: Dittus-Boelter errors reach 25 %.
    assert fc.dittus_boelter(Re, Pr) == pytest.approx(fc.gnielinski(Re, Pr), rel=0.25)


def test_dittus_boelter_by_hand():
    # 0.023 * 1e4^0.8 * 0.7^0.4 = 0.023 * 1584.89 * 0.867080 = 31.606; cooling: Pr^0.3.
    assert fc.dittus_boelter(1e4, 0.7) == pytest.approx(31.606, rel=1e-4)
    assert fc.dittus_boelter(1e4, 0.7, heating=False) == pytest.approx(
        0.023 * 1e4 ** 0.8 * 0.7 ** 0.3, rel=1e-12)


def test_annulus_by_hand():
    a = 0.5
    assert fc.nu_annulus(1000.0, 7.0, a) == (pytest.approx(3.66 + 1.2 * 2 ** 0.8), "laminar")
    assert fc.nu_annulus(5e4, 7.0, a)[0] == pytest.approx(
        0.86 * 2 ** 0.16 * fc.gnielinski(5e4, 7.0))


def test_coil_factor():
    assert fc.coil_factor(0.025, 0.5) == pytest.approx(1.175)


# ---- Kraume, Example 5.12: water in a double pipe ------------------------------------------

class ConstantProperties:
    """A package stand-in with Kraume's water: rho 1000, cp 4180, nu 1e-6, lambda 0.62."""

    def transport_properties(self, kind, T, P, zs):
        return TransportProps(rho=1000.0, cp=4180.0, mu=1e-3, k=0.62)


KRAUME_WATER = Fluid(ConstantProperties(), "liquid", (1.0,), 1e5)


def test_double_pipe_reynolds_numbers_match_kraume():
    # 1 m3/h in the inner tube (d_i 50 mm), and 30/28 of it in the annulus to the 100 mm
    # outer pipe (the inner tube's outer diameter 56 mm): Re 7074 and 2429 (Tab. 5.1, 5.2),
    # on the hydraulic diameter 0.044 m; Pr = 6.74.
    m_i = 1000.0 / 3600.0
    m_a = m_i * 30.0 / 28.0
    inner = tube_film(KRAUME_WATER, 300.0, m_i / (math.pi / 4 * 0.05 ** 2), 0.05,
                      "gnielinski", True)(300.0, 0.0)
    outer = annulus_film(KRAUME_WATER, 300.0, m_a / (math.pi / 4 * (0.1 ** 2 - 0.056 ** 2)),
                         0.056, 0.1, "gnielinski", False)(300.0, 0.0)
    assert inner.Re == pytest.approx(7074, abs=1)
    assert outer.Re == pytest.approx(2429, abs=1)
    assert KRAUME_WATER.at(300.0).Pr == pytest.approx(6.74, abs=0.005)


@pytest.mark.parametrize("h_i, h_o, k", [(875.0, 939.0, 415.0), (4554.0, 1313.0, 847.0)])
def test_overall_coefficient_of_a_plane_wall_matches_kraume(h_i, h_o, k):
    # Tab. 5.3: the film coefficients of both pipe pairs and k from Eq. 5.170 with the
    # 3 mm steel wall (15 W/mK) as a plane wall.
    o = fc.overall(fc.Side(300.0, fc.given_film(h_i)), fc.Side(320.0, fc.given_film(h_o)),
                   wall=0.003 / 15.0)
    assert o.U == pytest.approx(k, abs=0.5)


def test_overall_coefficient_of_a_tube_wall_by_hand():
    d_i, d_o, k_w, h_i, h_o, R_i, R_o = 0.02, 0.025, 16.0, 2000.0, 5000.0, 2e-4, 1e-4
    inner = fc.Side(300.0, fc.given_film(h_i), d_o / d_i, R_i)
    outer = fc.Side(400.0, fc.given_film(h_o), 1.0, R_o)
    o = fc.overall(inner, outer, fc.tube_wall_resistance(d_i, d_o, k_w, d_o))
    # 1/U = 1.25 (1/2000 + 2e-4) + 0.025 ln 1.25 / 32 + 1e-4 + 1/5000
    #     = 8.75e-4 + 1.743309e-4 + 3e-4 = 1.3493309e-3
    assert o.U == pytest.approx(1.0 / 1.3493309e-3, rel=1e-6)
    # The heat flux 100 K * U crosses each film: T_wall = T + q / h on each side's area.
    q = o.U * 100.0
    assert o.T_wall_inner == pytest.approx(300.0 + q * 1.25 / h_i, rel=1e-12)
    assert o.T_wall_outer == pytest.approx(400.0 - q / h_o, rel=1e-12)


# ---- shell side, Kern -------------------------------------------------------------------------

def test_kern_equivalent_diameters():
    inch = 0.0254
    # Kern, Fig. 28: 1 in tubes on 1 1/4 in and 3/4 in tubes on 1 in square pitch.
    assert fc.kern_equivalent_diameter(inch, 1.25 * inch, "square") / inch == pytest.approx(
        0.99, abs=0.005)
    assert fc.kern_equivalent_diameter(0.75 * inch, inch, "square") / inch == pytest.approx(
        0.95, abs=0.005)
    # Triangular pitch: the free area of half a tube in a triangle of height 0.86 p, over
    # half its perimeter; for 1 in on 1 1/4 in: 4 (0.671875 - 0.392699) / 1.570796 in.
    assert fc.kern_equivalent_diameter(inch, 1.25 * inch) / inch == pytest.approx(
        4 * (0.5 * 1.25 * 0.86 * 1.25 - math.pi / 8) / (math.pi / 2), rel=1e-12)


def test_kern_by_hand():
    # Cross-flow area (25.4 - 19.05)/25.4 * 0.5 * 0.2 m = 0.025 m²; Nu = 0.36 Re^0.55 Pr^1/3.
    assert fc.kern_crossflow_area(0.5, 0.01905, 0.0254, 0.2) == pytest.approx(0.025)
    assert fc.nu_kern(1e4, 5.0) == pytest.approx(0.36 * 158.489 * 1.709976, rel=1e-5)
    assert fc.nu_kern(1e4, 5.0, mu_ratio=2.0) == pytest.approx(
        fc.nu_kern(1e4, 5.0) * 2 ** 0.14)


# ---- stirred vessels --------------------------------------------------------------------------

def test_stirred_vessel_matches_kraume_example_18_2():
    # D = 1.76 m, impeller d = 1.17 m at 90 1/min, rho 1496, eta 1.85 mPa s, cp 2303,
    # lambda 0.488, C_turb 0.33: Re = 1.66e6, Pr = 8.73, Nu = 9.53e3, alpha_i = 2640.
    Re = fc.stirred_reynolds(1496.0, 1.5, 1.17, 1.85e-3)
    Pr = 2303.0 * 1.85e-3 / 0.488
    Nu = fc.nu_stirred(Re, Pr, 0.33, 2.0 / 3.0)
    assert Re == pytest.approx(1.66e6, rel=0.005)
    assert Pr == pytest.approx(8.73, abs=0.005)
    assert Nu == pytest.approx(9.53e3, rel=0.002)
    h_i = Nu * 0.488 / 1.76
    assert h_i == pytest.approx(2640, abs=5)
    # The enamel (1.3 mm, 1.16 W/mK) and steel (20 mm, 52 W/mK) layers and 1000 W/m2K in
    # the jacket: k = 347 W/m2K.
    wall = 0.0013 / 1.16 + 0.020 / 52.0
    o = fc.overall(fc.Side(350.0, fc.given_film(h_i)), fc.Side(290.0, fc.given_film(1000.0)),
                   wall)
    assert o.U == pytest.approx(347, abs=0.5)


def test_chilton_drew_jebens_constants():
    assert fc.STIRRED == {"jacket": (0.36, 2 / 3), "coil": (0.87, 0.62)}
    assert fc.nu_stirred(1e5, 4.0, 0.36, 2 / 3, mu_ratio=0.5) == pytest.approx(
        0.36 * 1e5 ** (2 / 3) * 4 ** (1 / 3) * 0.5 ** 0.14)


# ---- condensation and boiling -------------------------------------------------------------

# Saturated water at 1 atm, condensate near 95 °C (Incropera, Table A.6).
RHO_L, RHO_V, MU_L, K_L, CP_L, H_FG = 961.9, 0.596, 2.97e-4, 0.677, 4212.0, 2.257e6


def test_film_condensation_on_a_vertical_wall_solves_nusselts_film_equations():
    # Nusselt's film: mass flow per width Gamma = g rho (rho - rho_v) delta³ / (3 mu) grows
    # by the heat conducted through it, dGamma/dx = k dT / (delta h'). With delta from Gamma,
    # Gamma^(1/3) dGamma = A dx, so Gamma(L) = (4/3 A L)^(3/4), and the mean coefficient is
    # h' Gamma(L) / (L dT).
    dT, L = 10.0, 1.0
    h_mod = H_FG + 0.68 * CP_L * dT
    A = K_L * dT / h_mod * (9.80665 * RHO_L * (RHO_L - RHO_V) / (3 * MU_L)) ** (1 / 3)
    gamma = (4 / 3 * A * L) ** 0.75
    h = fc.h_condensation(RHO_L, RHO_V, MU_L, K_L, CP_L, H_FG, dT, L, "vertical")
    assert h == pytest.approx(h_mod * gamma / (L * dT), rel=1e-3)  # 0.943 is rounded


def test_film_condensation_on_a_horizontal_tube_solves_nusselts_film_equations():
    # Around a tube of radius R, the film drains with g sin(phi): Gamma^(1/3) dGamma =
    # A R sin(phi)^(1/3) dphi on each side, so Gamma(pi) = (4/3 A R I)^(3/4) with I the
    # integral of sin^(1/3) from 0 to pi; both sides over the circumference 2 pi R.
    dT, D = 10.0, 0.025
    R = D / 2
    h_mod = H_FG + 0.68 * CP_L * dT
    A = K_L * dT / h_mod * (9.80665 * RHO_L * (RHO_L - RHO_V) / (3 * MU_L)) ** (1 / 3)
    integral = quad(lambda phi: math.sin(phi) ** (1 / 3), 0.0, math.pi)[0]
    gamma = (4 / 3 * A * R * integral) ** 0.75
    h_ref = 2 * gamma * h_mod / (2 * math.pi * R * dT)
    h = fc.h_condensation(RHO_L, RHO_V, MU_L, K_L, CP_L, H_FG, dT, D, "horizontal_tube")
    assert h == pytest.approx(h_ref, rel=2e-3)  # Nusselt's 0.728, Incropera's 0.729


def cooper(q, pr, M):
    """Cooper (1984) nucleate boiling, smooth surface (R_p = 1 micron), W/m2K."""
    return 55.0 * pr ** 0.12 * (-math.log10(pr)) ** -0.55 * M ** -0.5 * q ** 0.67


def test_mostinski_by_hand_and_against_cooper_for_water():
    # Water at 1 atm, q = 1e5 W/m2: P_c = 220.64 bar, p_r = 0.0045923;
    # F = 1.8 p_r^0.17 + 4 p_r^1.2 + 10 p_r^10 = 0.720805 + 0.006259 = 0.727064,
    # h = 0.104 * 220.64^0.69 * 1e5^0.7 * F = 0.104 * 41.4135 * 3162.28 * 0.727064.
    h = fc.h_nucleate_boiling(1e5, 101325.0, 220.64e5)
    assert h == pytest.approx(0.104 * 41.4135 * 3162.28 * 0.727064, rel=1e-5)
    # Cooper's correlation, fitted to other data, gives 9530 W/m2K: within 5 %.
    assert h == pytest.approx(cooper(1e5, 101325.0 / 220.64e5, 18.015), rel=0.05)


def test_wall_dependent_films_iterate_to_a_consistent_wall():
    # Steam condensing outside a tube, water inside: the condensate film's dT is the one its
    # coefficient was evaluated at, and the heat flux is the same through both films.
    w = Fluid(get_package(("water",), "ideal"), "liquid", (1.0,), 5e5)

    def condensing(T_wall, q):
        dT = 400.0 - T_wall
        return fc.Film(fc.h_condensation(RHO_L, RHO_V, MU_L, K_L, CP_L, H_FG, dT, 0.025,
                                          "horizontal_tube"), "test")

    inner = fc.Side(330.0, tube_film(w, 330.0, 800.0, 0.02, "gnielinski", True), 1.25)
    outer = fc.Side(400.0, condensing, 1.0, wall_dependent=True)
    o = fc.overall(inner, outer, fc.tube_wall_resistance(0.02, 0.025, 16.0, 0.025))
    q = o.U * 70.0
    assert o.outer.h == pytest.approx(condensing(o.T_wall_outer, q).h, rel=1e-8)
    assert o.outer.h * (400.0 - o.T_wall_outer) == pytest.approx(q, rel=1e-8)
    assert o.inner.h * (o.T_wall_inner - 330.0) == pytest.approx(q * 1.25, rel=1e-8)


def test_a_boiling_film_iterates_to_its_own_heat_flux():
    film = (lambda T_wall, q: fc.Film(fc.h_nucleate_boiling(q, 101325.0, 220.64e5), "b"))
    o = fc.overall(fc.Side(373.15, film, wall_dependent=True),
                   fc.Side(410.0, fc.given_film(8000.0)), wall=5e-4)
    q = o.U * (410.0 - 373.15)
    assert o.inner.h == pytest.approx(fc.h_nucleate_boiling(q, 101325.0, 220.64e5), rel=1e-7)


# ---- units rated from their geometry ----------------------------------------------------------

def utility(model):
    return next(iter(parse_utilities([UtilityModel(**model)]).values()))


def one_unit(unit, *, flows=None, T="25 degC", P="3 bar", components=("water",),
             outlets=("out",), utilities=(STEAM, CW, OIL), method="ideal"):
    r = solve({
        "components": list(components), "property_method": method,
        "utilities": list(utilities),
        "units": [{"id": "F", "type": "Feed",
                   "params": {"T": T, "P": P, "flows": flows or {"water": "20 kmol/h"}}},
                  {"id": "U", **unit}],
        "streams": [{"id": "IN", "from": "F.out", "to": "U.in"}]
        + [{"id": port.upper(), "from": f"U.{port}"} for port in outlets],
    })
    assert r.ok, r.error
    return r


def water_at(T, P=3e5):
    return get_package(("water",), "ideal").transport_properties("liquid", T, P, [1.0])


def gnielinski_tube(Re, Pr, laminar=3.66, factor=1.0):
    """The tube correlation worked by hand: laminar up to 2300, Gnielinski (times
    ``factor``) from 1e4, linear in between."""
    if Re <= 2300:
        return laminar
    f8 = (0.790 * math.log(max(Re, 1e4)) - 1.64) ** -2 / 8
    turb = factor * f8 * (max(Re, 1e4) - 1000) * Pr / (
        1 + 12.7 * math.sqrt(f8) * (Pr ** (2 / 3) - 1))
    if Re >= 1e4:
        return turb
    g = (Re - 2300) / 7700
    return (1 - g) * laminar + g * turb


def annulus_by_hand(Re, Pr, a):
    return gnielinski_tube(Re, Pr, 3.66 + 1.2 * a ** -0.8, 0.86 * a ** -0.16)


DOUBLE_PIPE = {"U_from": "double_pipe", "tube_outer_diameter": "56 mm",
               "wall_thickness": "3 mm", "annulus_diameter": "100 mm", "wall_conductivity": 15}


def kraume_exchanger(**params):
    """Kraume's Example 5.12: 1 m3/h of water from 20 °C in the inner tube, hot water from
    60 °C in the annulus at 30/28 of that flow, 45.9 m long."""
    m_c = 1000.0 / 3600.0
    m_h = m_c * 30.0 / 28.0
    area = math.pi * 0.056 * 45.9
    r = solve({
        "components": ["water"], "property_method": "ideal",
        "units": [
            {"id": "H", "type": "Feed", "params": {"T": "60 degC", "P": "3 bar",
                                                   "flows": {"water": m_h / MW_WATER}}},
            {"id": "C", "type": "Feed", "params": {"T": "20 degC", "P": "3 bar",
                                                   "flows": {"water": m_c / MW_WATER}}},
            {"id": "X", "type": "HeatExchanger",
             "params": {"area": area, **DOUBLE_PIPE, **params}}],
        "streams": [{"id": "HI", "from": "H.out", "to": "X.hot_in"},
                    {"id": "CI", "from": "C.out", "to": "X.cold_in"},
                    {"id": "HO", "from": "X.hot_out"}, {"id": "CO", "from": "X.cold_out"}]})
    assert r.ok, r.error
    return r, m_h, m_c, area


def test_double_pipe_exchanger_rated_from_geometry_matches_the_closed_form():
    r, m_h, m_c, area = kraume_exchanger()
    u = r.units["X"]
    T_hi, T_ci = 333.15, 293.15
    T_ho, T_co = r.streams["HO"]["T"], r.streams["CO"]["T"]
    Q = u["duty"]
    # The zone analysis needs the UA the film coefficients give.
    assert u["UA"] == pytest.approx(u["U"] * area, rel=1e-6)
    assert u["heat_transfer_area"] == pytest.approx(area)
    # Effectiveness-NTU of a countercurrent double pipe at that U, with each side's mean
    # heat capacity (liquid water: within 0.1 % constant over 20-60 °C).
    C_h, C_c = Q / (T_hi - T_ho), Q / (T_co - T_ci)
    C_min, C_max = min(C_h, C_c), max(C_h, C_c)
    NTU, Cr = u["U"] * area / C_min, C_min / C_max
    eps = (1 - math.exp(-NTU * (1 - Cr))) / (1 - Cr * math.exp(-NTU * (1 - Cr)))
    assert Q == pytest.approx(eps * C_min * (T_hi - T_ci), rel=1e-3)

    # U by hand from the water's properties at each side's mean temperature: Gnielinski in
    # the inner tube (d_i 50 mm), the annulus 56/100 mm on its hydraulic diameter 44 mm,
    # the steel wall on the outer tube surface.
    cold, hot = water_at(0.5 * (T_ci + T_co)), water_at(0.5 * (T_hi + T_ho))
    Re_i = m_c / (math.pi / 4 * 0.05 ** 2) * 0.05 / cold.mu
    h_i = gnielinski_tube(Re_i, cold.Pr) * cold.k / 0.05
    Re_o = m_h / (math.pi / 4 * (0.1 ** 2 - 0.056 ** 2)) * 0.044 / hot.mu
    h_o = annulus_by_hand(Re_o, hot.Pr, 0.56) * hot.k / 0.044
    U = 1 / (0.056 / 0.05 / h_i + 0.056 * math.log(0.056 / 0.05) / 30 + 1 / h_o)
    assert u["Re_tube"] == pytest.approx(Re_i, rel=1e-9)
    assert u["Re_shell"] == pytest.approx(Re_o, rel=1e-9)
    assert u["h_tube"] == pytest.approx(h_i, rel=1e-9)
    assert u["h_shell"] == pytest.approx(h_o, rel=1e-9)
    assert u["U"] == pytest.approx(U, rel=1e-9)
    # Both flows are in the transition region (Kraume, who counts Re > 2040 as turbulent,
    # gets 415 W/m2K with his correlations).
    assert 2300 < Re_o < Re_i < 1e4


def test_swapping_the_sides_swaps_the_films():
    r = kraume_exchanger(tube_side="hot")[0]
    u = r.units["X"]
    hot = r.streams["HI"]["T"]
    assert hot > r.streams["CI"]["T"]
    m_h = 1000.0 / 3600.0 * 30 / 28
    T_hot = 0.5 * (hot + r.streams["HO"]["T"])
    assert u["Re_tube"] == pytest.approx(
        m_h / (math.pi / 4 * 0.05 ** 2) * 0.05 / water_at(T_hot).mu, rel=1e-9)


def test_exchanger_with_a_temperature_reports_the_area_its_films_need():
    r, _, _, area = kraume_exchanger()
    T_ho = r.streams["HO"]["T"]
    # The same exchanger given its hot outlet temperature needs the same area.
    spec = solve({
        "components": ["water"], "property_method": "ideal",
        "units": [
            {"id": "H", "type": "Feed", "params": {"T": "60 degC", "P": "3 bar",
                                                   "flows": {"water": 1 / 3.6 * 30 / 28
                                                             / MW_WATER}}},
            {"id": "C", "type": "Feed", "params": {"T": "20 degC", "P": "3 bar",
                                                   "flows": {"water": 1 / 3.6 / MW_WATER}}},
            {"id": "X", "type": "HeatExchanger", "params": {"T_hot_out": T_ho, **DOUBLE_PIPE}}],
        "streams": [{"id": "HI", "from": "H.out", "to": "X.hot_in"},
                    {"id": "CI", "from": "C.out", "to": "X.cold_in"},
                    {"id": "HO", "from": "X.hot_out"}, {"id": "CO", "from": "X.cold_out"}]})
    assert spec.ok, spec.error
    assert spec.units["X"]["heat_transfer_area"] == pytest.approx(area, rel=1e-5)


HEATER_TUBE = {"utility": "LPS", "area": "2 m2", "U_from": "double_pipe",
               "tube_outer_diameter": "25 mm", "wall_thickness": "2 mm",
               "annulus_diameter": "50 mm"}


def test_steam_heater_rated_from_geometry_follows_the_exponential_approach():
    r = one_unit({"type": "Heater", "params": HEATER_TUBE})
    u = r.units["U"]
    T_in, T_out, Ts = r.streams["IN"]["T"], r.streams["OUT"]["T"], utility(STEAM).T_in
    UA = u["U"] * 2.0
    assert u["UA"] == pytest.approx(UA, rel=1e-12)
    assert u["duty"] == pytest.approx(UA * lmtd(Ts - T_in, Ts - T_out), rel=1e-6)
    C = u["duty"] / (T_out - T_in)
    assert T_out == pytest.approx(Ts - (Ts - T_in) * math.exp(-UA / C), abs=1e-5)
    # Inside: Gnielinski at the mean temperature, 20 kmol/h of water in a 21 mm tube.
    m = 20 / 3.6 * MW_WATER
    w = water_at(0.5 * (T_in + T_out))
    Re = m / (math.pi / 4 * 0.021 ** 2) * 0.021 / w.mu
    assert u["Re_process"] == pytest.approx(Re, rel=1e-9) and Re > 1e4
    assert u["h_process"] == pytest.approx(gnielinski_tube(Re, w.Pr) * w.k / 0.021, rel=1e-9)
    # Outside: Nusselt condensation on the 25 mm tube, at the wall temperature it reports.
    assert u["h_medium"] > 5 * u["h_process"]
    assert T_in < u["T_wall"] < Ts


def test_a_given_U_reproduces_the_rated_one():
    rated = one_unit({"type": "Heater", "params": HEATER_TUBE})
    given = one_unit({"type": "Heater", "params": {"utility": "LPS", "area": "2 m2",
                                                    "U": rated.units["U"]["U"]}})
    assert given.units["U"]["duty"] == pytest.approx(rated.units["U"]["duty"], rel=1e-8)
    assert given.streams["OUT"]["T"] == pytest.approx(rated.streams["OUT"]["T"], abs=1e-7)
    assert "h_process" not in given.units["U"] or given.units["U"]["h_process"] is None


def test_dittus_boelter_option():
    r = one_unit({"type": "Heater", "params": {**HEATER_TUBE,
                                                "tube_correlation": "dittus-boelter"}})
    u = r.units["U"]
    w = water_at(0.5 * (r.streams["IN"]["T"] + r.streams["OUT"]["T"]))
    assert u["Nu_process"] == pytest.approx(
        0.023 * u["Re_process"] ** 0.8 * w.Pr ** 0.4, rel=1e-9)


def test_shell_and_tube_cooler_with_cooling_water_on_the_shell_side():
    r = one_unit({"type": "Heater", "params": {
        "utility": "CW", "area": "6 m2", "U_from": "shell_and_tube", "tubes": 40,
        "tube_passes": 2, "tube_outer_diameter": "19.05 mm", "wall_thickness": "1.65 mm",
        "shell_diameter": "0.3 m", "tube_pitch": "23.8 mm", "pitch_layout": "square",
        "baffle_spacing": "0.15 m", "fouling_medium": "0.0002 m2K/W"}},
        flows={"water": "100 kmol/h"}, T="80 degC")
    u = r.units["U"]
    assert u["duty"] < 0.0
    # Shell side by hand: Kern with the cooling water at its mean temperature; the flow
    # follows from the duty and the utility's 10 K rise.
    m_cw = -u["duty"] / (4180 * 10)
    assert u["utility_flow"] == pytest.approx(m_cw, rel=1e-9)
    cw = water_at(303.15, 5e5)
    D_e = 4 * (0.0238 ** 2 - math.pi * 0.01905 ** 2 / 4) / (math.pi * 0.01905)
    A_s = (0.0238 - 0.01905) * 0.3 * 0.15 / 0.0238
    Re = m_cw / A_s * D_e / cw.mu
    assert u["Re_medium"] == pytest.approx(Re, rel=1e-9)
    mu_ratio = cw.mu / water_at(u["T_wall"] + 0.0, 5e5).mu  # wall correction
    assert 0.5 < mu_ratio < 1.5
    Nu_no_wall = 0.36 * Re ** 0.55 * cw.Pr ** (1 / 3)
    assert u["Nu_medium"] == pytest.approx(Nu_no_wall, rel=0.08)
    # Tube side: 40 tubes in 2 passes, 20 carry the flow.
    m = 100 / 3.6 * MW_WATER
    d_i = 0.01905 - 0.0033
    w = water_at(0.5 * (r.streams["IN"]["T"] + r.streams["OUT"]["T"]))
    assert u["Re_process"] == pytest.approx(m / 20 / (math.pi / 4 * d_i ** 2) * d_i / w.mu,
                                            rel=1e-9)
    # U on the outer surface includes the fouling resistance.
    R = (0.01905 / d_i / u["h_process"] + 0.01905 * math.log(0.01905 / d_i) / 32
         + 0.0002 + 1 / u["h_medium"])
    assert u["U"] == pytest.approx(1 / R, rel=1e-9)


JACKET = {"U_from": "jacket", "vessel_diameter": "1.2 m", "impeller_diameter": "0.4 m",
          "impeller_speed": "120 rpm", "wall_thickness": "8 mm", "h_medium": "1500 W/m2K"}


def test_jacketed_drum_by_chilton_drew_jebens():
    r = one_unit({"type": "Flash", "params": {"utility": "CW", "area": "3 m2", **JACKET}},
                 T="80 degC", outlets=("vapor", "liquid"))
    u = r.units["U"]
    T = r.streams["LIQUID"]["T"]
    w = water_at(T)
    Re = w.rho * 2.0 * 0.4 ** 2 / w.mu
    mu_ratio = w.mu / water_at(u["T_wall"]).mu
    h = 0.36 * Re ** (2 / 3) * w.Pr ** (1 / 3) * mu_ratio ** 0.14 * w.k / 1.2
    assert u["Re_process"] == pytest.approx(Re, rel=1e-9)
    assert u["h_process"] == pytest.approx(h, rel=1e-6)
    assert u["U"] == pytest.approx(1 / (1 / h + 0.008 / 16 + 1 / 1500), rel=1e-6)
    # The jacket passes U A LMTD at the drum temperature.
    assert -u["duty"] == pytest.approx(u["U"] * 3.0 * u["LMTD"], rel=1e-6)
    # The wall is the temperature that passes the heat flux through the process film.
    assert u["T_wall"] == pytest.approx(T + u["duty"] / 3.0 / h, rel=1e-6)


def test_jacketed_drum_at_a_given_temperature_reports_the_area():
    r = one_unit({"type": "Flash", "params": {"utility": "CW", "T": "50 degC", **JACKET}},
                 T="80 degC", outlets=("vapor", "liquid"))
    u = r.units["U"]
    assert u["heat_transfer_area"] == pytest.approx(-u["duty"] / (u["U"] * u["LMTD"]),
                                                    rel=1e-9)


def test_cooling_water_in_the_jacket_gap_and_in_a_coil():
    params = {k: v for k, v in JACKET.items() if k != "h_medium"}
    jacket = one_unit({"type": "Flash", "params": {"utility": "CW", "area": "3 m2",
                                                   "jacket_gap": "30 mm", **params}},
                      T="80 degC", outlets=("vapor", "liquid")).units["U"]
    # The annular jacket between 1.216 and 1.276 m, on its hydraulic diameter 60 mm.
    m = jacket["utility_flow"]
    cw = water_at(303.15, 5e5)
    Re = m / (math.pi / 4 * (1.276 ** 2 - 1.216 ** 2)) * 0.06 / cw.mu
    assert jacket["Re_medium"] == pytest.approx(Re, rel=1e-9)
    r = one_unit({"type": "Flash", "params": {
        "utility": "CW", "area": "2 m2", **params, "U_from": "coil",
        "tube_outer_diameter": "33.7 mm", "wall_thickness": "2.6 mm",
        "coil_diameter": "0.9 m", "utility_flow": "1 kg/s"}},
        T="80 degC", outlets=("vapor", "liquid"))
    coil = r.units["U"]
    d_i = 0.0285
    T_cw = 0.5 * (298.15 + coil["utility_T_out"])
    cw = water_at(T_cw, 5e5)
    Re = 1.0 / (math.pi / 4 * d_i ** 2) * d_i / cw.mu
    assert coil["Re_medium"] == pytest.approx(Re, rel=1e-6)
    assert coil["Nu_medium"] == pytest.approx(
        gnielinski_tube(Re, cw.Pr) * (1 + 3.5 * d_i / 0.9), rel=1e-6)
    # The coil's process side: Chilton-Drew-Jebens with the coil constants.
    w = water_at(r.streams["LIQUID"]["T"])
    Re = w.rho * 2.0 * 0.4 ** 2 / w.mu
    mu_ratio = w.mu / water_at(coil["T_wall"]).mu
    assert coil["Nu_process"] == pytest.approx(
        0.87 * Re ** 0.62 * w.Pr ** (1 / 3) * mu_ratio ** 0.14, rel=1e-6)


def test_boiling_drum_on_a_steam_jacket():
    r = one_unit({"type": "Flash", "params": {
        "utility": "LPS", "area": "2 m2", "U_from": "jacket", "vessel_diameter": "1.2 m",
        "process_film": "boiling", "wall_thickness": "8 mm"}},
        components=("water", "ethylene glycol"),
        flows={"water": "20 kmol/h", "ethylene glycol": "5 kmol/h"}, T="95 degC",
        P="1 atm", outlets=("vapor", "liquid"))
    u = r.units["U"]
    assert r.streams["VAPOR"]["molar_flow"] > 0.0
    # Both temperatures are uniform: the heat flux through the boiling film is Q / A.
    q = u["duty"] / 2.0
    x = r.streams["LIQUID"]["liquid_composition"]
    pkg = get_package(("water", "ethylene glycol"), "ideal")
    Pc = sum(xi * pc for xi, pc in zip([x["water"], x["ethylene glycol"]], pkg.consts.Pcs))
    assert u["h_process"] == pytest.approx(fc.h_nucleate_boiling(q, 101325.0, Pc), rel=1e-6)
    # The steam condenses on the jacket wall as on a vertical wall 1.2 m high.
    assert 3000 < u["h_medium"] < 20000


@pytest.mark.parametrize("unit, message", [
    ({"type": "Heater", "params": {"utility": "LPS", "area": "1 m2", "U": 300,
                                   **{k: v for k, v in HEATER_TUBE.items()
                                      if k not in ("utility", "area")}}}, "over-specified"),
    ({"type": "Heater", "params": {"utility": "LPS", "area": "1 m2", "U_from": "double_pipe",
                                   "tube_outer_diameter": "25 mm", "wall_thickness": "2 mm"}},
     "outer pipe inner diameter"),
    ({"type": "Heater", "params": {"utility": "OIL", "area": "1 m2",
                                   **{k: v for k, v in HEATER_TUBE.items()
                                      if k not in ("utility", "area")}}}, "computed for water"),
    ({"type": "Heater", "params": {**HEATER_TUBE, "annulus_diameter": "20 mm"}},
     "wider than the tube"),
    ({"type": "Flash", "params": {"utility": "CW", "area": "1 m2", "U_from": "jacket",
                                  "vessel_diameter": "1 m"}}, "impeller diameter"),
    ({"type": "Flash", "params": {"utility": "CW", "area": "1 m2", "U_from": "jacket",
                                  "vessel_diameter": "1 m", "impeller_diameter": "0.3 m",
                                  "impeller_speed": "100 rpm"}}, "jacket gap"),
])
def test_film_mistakes_are_explained(unit, message):
    r = solve({
        "components": ["water"], "property_method": "ideal", "utilities": [STEAM, CW, OIL],
        "units": [{"id": "F", "type": "Feed", "params": {"T": "60 degC", "P": "3 bar",
                                                         "flows": {"water": "20 kmol/h"}}},
                  {"id": "U", **unit}],
        "streams": [{"id": "IN", "from": "F.out", "to": "U.in"}]
        + [{"id": p, "from": f"U.{p}"} for p in (("vapor", "liquid") if unit["type"] == "Flash"
                                                 else ("out",))]})
    assert not r.ok
    assert message in r.error["message"] + (r.error.get("hint") or "")


# ---- tubular reactor --------------------------------------------------------------------------

INERT = ["argon", "n-butane", "isobutane"]
NO_RATE = {"id": "ISO", "stoichiometry": {"n-butane": -1, "isobutane": 1},
           "kinetics": {"k0": 1.0}}
PFR_CW = {"id": "CW", "type": "fluid", "T_in": "300 K", "T_out": "320 K", "cp": "4180 J/kgK"}


def argon_pfr(params):
    r = solve({
        "components": INERT, "property_method": "ideal", "reactions": [NO_RATE],
        "utilities": [PFR_CW],
        "units": [{"id": "F", "type": "Feed",
                   "params": {"T": 600.0, "P": "2 bar", "flows": {"argon": 10.0}}},
                  {"id": "R", "type": "PFR",
                   "params": {"length": 10.0, "diameter": 0.05, "phase": "gas",
                              "reactions": ["ISO"], "utility": "CW",
                              "utility_flow": "0.5 kg/s", **params}}],
        "streams": [{"id": "IN", "from": "F.out", "to": "R.in"},
                    {"id": "OUT", "from": "R.out"}]})
    assert r.ok, r.error
    return r.units["R"]


def profile(unit, label):
    c = unit["profiles"]["conditions"]
    return np.array(c["x"]), np.array(next(s["values"] for s in c["series"]
                                           if s["label"] == label))


@pytest.mark.parametrize("arrangement", ["cocurrent", "countercurrent"])
def test_tubular_reactor_evaluates_U_along_the_tube(arrangement):
    u = argon_pfr({"U_from": "double_pipe", "annulus_diameter": "80 mm",
                   "wall_thickness": "3 mm", "arrangement": arrangement})
    x, U = profile(u, "U")
    _, T = profile(u, "Temperature")
    _, T_m = profile(u, "Medium temperature")
    pkg = get_package(tuple(INERT), "ideal")
    G = 10.0 * 39.948e-3 / (math.pi / 4 * 0.05 ** 2)
    m_cw = 0.5
    for i in (0, 20, 40):
        gas = pkg.transport_properties("vapor", T[i], 2e5, [1.0, 0.0, 0.0])
        Re = G * 0.05 / gas.mu
        h_i = gnielinski_tube(Re, gas.Pr) * gas.k / 0.05
        cw = water_at(T_m[i], 5e5)
        Re_o = m_cw / (math.pi / 4 * (0.08 ** 2 - 0.056 ** 2)) * 0.024 / cw.mu
        h_o = annulus_by_hand(Re_o, cw.Pr, 0.056 / 0.08) * cw.k / 0.024
        # On the inner surface: 1/U = 1/h_i + D ln(d_o/D) / (2 k_w) + (D/d_o) / h_o.
        U_ref = 1 / (1 / h_i + 0.05 * math.log(0.056 / 0.05) / 32 + 0.05 / 0.056 / h_o)
        assert U[i] == pytest.approx(U_ref, rel=1e-6)
    # The heat through the wall is the integral of U pi D (T_m - T) along the tube.
    Q = np.trapezoid(U * math.pi * 0.05 * (T_m - T), x)
    assert u["wall_duty"] == pytest.approx(Q, rel=2e-3)
    assert u["U"] == pytest.approx(np.trapezoid(U, x) / 10.0, rel=1e-9)
    assert u["UA"] == pytest.approx(u["U"] * math.pi * 0.05 * 10.0, rel=1e-9)


def test_tubular_reactor_with_a_given_outside_coefficient():
    u = argon_pfr({"U_from": "tube", "h_medium": "800 W/m2K"})
    _, U = profile(u, "U")
    _, T = profile(u, "Temperature")
    gas = get_package(tuple(INERT), "ideal").transport_properties("vapor", T[0], 2e5,
                                                                  [1.0, 0.0, 0.0])
    Re = 10.0 * 39.948e-3 / (math.pi / 4 * 0.05 ** 2) * 0.05 / gas.mu
    h_i = gnielinski_tube(Re, gas.Pr) * gas.k / 0.05
    # No wall thickness: no wall resistance, and the outer surface is the inner one.
    assert U[0] == pytest.approx(1 / (1 / h_i + 1 / 800.0), rel=1e-6)
    assert u["h_medium"] == pytest.approx(800.0)


@pytest.mark.parametrize("params, message", [
    ({"U_from": "tube"}, "annulus it flows in"),
    ({"U_from": "tube", "U": 100, "h_medium": 500}, "over-specified"),
    ({"U_from": "tube", "h_medium": 500, "area_per_length": 0.2}, "inner tube surface"),
    ({"U_from": "double_pipe", "annulus_diameter": "50 mm", "h_medium": 500},
     "wider than the tube"),
])
def test_tubular_reactor_film_mistakes_are_explained(params, message):
    with pytest.raises(AssertionError, match=message):
        argon_pfr(params)
