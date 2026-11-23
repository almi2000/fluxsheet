"""U of a unit from its geometry: builds both sides of the wall for the correlations of
``fluxsheet.film_coefficients`` from the unit's parameters, its streams and its utility.

A unit with ``U_from`` set finds U this way instead of taking it as a number. The process
side's properties come from the unit's property package
(``PropertyPackage.transport_properties``), the utility side's from water: a liquid
utility is taken to be water (its heat capacity must be within 10 % of water's), steam
condenses to water. Other utilities need their film coefficient given (``h_medium``).
"""

from __future__ import annotations

import functools
import math
from collections.abc import Callable
from dataclasses import dataclass

from chemicals.iapws import iapws95_rhog_sat

from . import film_coefficients as fc
from .errors import SimulationError
from .utilities import _ends, lmtd

WATER_CP = 4186.0  # J/(kg K), to recognize a liquid utility as water
WATER_P = 5e5  # Pa, liquid water for its properties (they barely depend on it)
VESSEL_GEOMETRIES = ("jacket", "coil")
TUBE_GEOMETRIES = ("double_pipe", "shell_and_tube")

FilmFn = Callable[[float, float], fc.Film]


def _degC(T: float) -> str:
    return f"{T - 273.15:.1f} °C"


@functools.lru_cache(maxsize=1)
def _water_package():
    from .thermo_package import PropertyPackage  # built once, outside the flowsheet cache

    return PropertyPackage(["water"], "ideal")


@dataclass(frozen=True)
class Fluid:
    """A fluid of fixed composition and phase, its properties evaluated by temperature."""

    pkg: object
    kind: str  # "liquid" or "vapor"
    zs: tuple[float, ...]
    P: float

    def at(self, T: float):
        return self.pkg.transport_properties(self.kind, T, self.P, list(self.zs))

    def mu_ratio(self, mu: float, T_wall: float) -> float:
        """mu / mu_wall; 1 where the wall viscosity cannot be evaluated."""
        try:
            return mu / self.at(T_wall).mu
        except SimulationError:
            return 1.0


def water() -> Fluid:
    return Fluid(_water_package(), "liquid", (1.0,), WATER_P)


def stream_fluid(pkg, zs, T: float, P: float, notes: list[str], label: str) -> Fluid:
    """The phase a stream of composition ``zs`` has at T and P. A two-phase stream is rated
    as its larger phase: the single-phase correlations do not cover boiling or condensing
    flow."""
    st = pkg.flash(list(zs), T=T, P=P)
    if st.phase == "vapor" or st.x is None:
        return Fluid(pkg, "vapor", tuple(st.y or zs), P)
    if st.phase == "liquid" or st.y is None:
        return Fluid(pkg, "liquid", tuple(st.x), P)
    kind = "vapor" if st.VF >= 0.5 else "liquid"
    note = (f"{label} is two-phase at {_degC(T)} (vapor fraction {st.VF:.2f}); its film "
            f"coefficient is for single-phase {kind} flow")
    if note not in notes:
        notes.append(note)
    return Fluid(pkg, kind, tuple(st.y if kind == "vapor" else st.x), P)


def mass_flow(pkg, stream) -> float:
    """kg/s"""
    return sum(f * mw for f, mw in zip(stream.flows, pkg.mws)) / 1e3


# ---- films ------------------------------------------------------------------------------------
# Each returns ``film(T_wall, q)`` for fluxsheet.film_coefficients.Side; ``T`` is the bulk
# temperature the fluid's properties are taken at.

def tube_film(fluid: Fluid, T: float, mass_flux: float, d: float, correlation: str,
              heated: bool, factor: float = 1.0, label: str = "tube") -> FilmFn:
    """Inside a tube of diameter ``d``; ``factor`` multiplies Nu (a coil)."""
    pr = fluid.at(T)
    Re = mass_flux * d / pr.mu
    Nu, regime = fc.nu_tube(Re, pr.Pr, correlation, heated)
    Nu *= factor
    name = "Gnielinski" if correlation == "gnielinski" else "Dittus-Boelter"
    result = fc.Film(Nu * pr.k / d, f"{name}, {label}", Re, Nu, regime)
    return lambda T_wall, q: result


def annulus_film(fluid: Fluid, T: float, mass_flux: float, d_i: float, D_o: float,
                 correlation: str, heated: bool) -> FilmFn:
    """In the annulus between a tube of outer diameter ``d_i`` and a pipe of inner diameter
    ``D_o``, heat through the inner wall."""
    D_h = D_o - d_i
    pr = fluid.at(T)
    Re = mass_flux * D_h / pr.mu
    Nu, regime = fc.nu_annulus(Re, pr.Pr, d_i / D_o, correlation, heated)
    result = fc.Film(Nu * pr.k / D_h, "annulus", Re, Nu, regime)
    return lambda T_wall, q: result


def kern_film(fluid: Fluid, T: float, mass_flow: float, d_o: float, pitch: float,
              layout: str, D_s: float, baffle: float, notes: list[str]) -> FilmFn:
    """Shell side of a baffled bundle, Kern's method with the wall viscosity correction."""
    D_e = fc.kern_equivalent_diameter(d_o, pitch, layout)
    G = mass_flow / fc.kern_crossflow_area(D_s, d_o, pitch, baffle)
    pr = fluid.at(T)
    Re = G * D_e / pr.mu
    lo, hi = fc.KERN_RE_RANGE
    if not lo <= Re <= hi:
        note = (f"the shell-side Reynolds number {Re:.3g} is outside the range of Kern's "
                f"correlation ({lo:.0f} to {hi:.0e})")
        if note not in notes:
            notes.append(note)

    def film(T_wall: float, q: float) -> fc.Film:
        Nu = fc.nu_kern(max(Re, 1e-12), pr.Pr, fluid.mu_ratio(pr.mu, T_wall))
        return fc.Film(Nu * pr.k / D_e, "Kern shell side", Re, Nu,
                       "turbulent" if lo <= Re else "below range")

    return film


def stirred_film(fluid: Fluid, T: float, D_T: float, d: float, N: float, C: float,
                 a: float) -> FilmFn:
    pr = fluid.at(T)
    Re = fc.stirred_reynolds(pr.rho, N, d, pr.mu)

    def film(T_wall: float, q: float) -> fc.Film:
        Nu = fc.nu_stirred(Re, pr.Pr, C, a, fluid.mu_ratio(pr.mu, T_wall))
        return fc.Film(Nu * pr.k / D_T, "Chilton-Drew-Jebens stirred vessel", Re, Nu,
                       "turbulent" if Re > 1e4 else "transition" if Re > 10.0 else "laminar")

    return film


def condensing_film(T_sat: float, latent: float, length: float, geometry: str) -> FilmFn:
    """Steam condensing on the wall (Nusselt film), the condensate's properties at the film
    temperature."""
    rho_v = iapws95_rhog_sat(T_sat)
    w = water()

    def film(T_wall: float, q: float) -> fc.Film:
        dT = max(T_sat - T_wall, 1e-3)
        pr = w.at(T_sat - dT / 2.0)
        h = fc.h_condensation(pr.rho, rho_v, pr.mu, pr.k, pr.cp, latent, dT, length,
                              geometry)
        return fc.Film(h, "Nusselt film condensation (laminar film; simplified)", None,
                       h * length / pr.k, "condensing")

    return film


def boiling_film(pkg, zs, P: float) -> FilmFn:
    Pc = sum(z * pc for z, pc in zip(zs, pkg.consts.Pcs))

    def film(T_wall: float, q: float) -> fc.Film:
        return fc.Film(fc.h_nucleate_boiling(q, P, Pc),
                       "Mostinski nucleate pool boiling (simplified)", regime="boiling")

    return film


# ---- geometry -----------------------------------------------------------------------------------

def _param(unit, name: str) -> float:
    value = unit.p.get(name)
    if value is None or value <= 0.0:
        label = next(p.label for p in unit.params if p.name == name)
        raise SimulationError(f"U from film coefficients needs the {label[0].lower()}"
                              f"{label[1:]}", where=unit.id,
                              hint=f"set {name}, or give U instead")
    return value


@dataclass
class Tubes:
    """Tubes of a double pipe or a shell-and-tube bundle, from the unit's parameters."""

    kind: str  # "double_pipe" or "shell_and_tube"
    n: int  # tubes
    passes: int
    d_o: float
    d_i: float
    wall: float  # m²K/W on the outer tube surface
    correlation: str
    D_a: float | None = None  # inner diameter of the outer pipe
    shell: tuple[float, str, float, float] | None = None  # pitch, layout, D_s, baffles

    @classmethod
    def of(cls, unit) -> Tubes:
        p = unit.p
        kind = p["U_from"]
        n = int(p.get("tubes") or 1)
        passes = int(p.get("tube_passes") or 1) if kind == "shell_and_tube" else 1
        if n % passes:
            raise SimulationError(f"{n} tubes do not divide evenly into {passes} passes",
                                  where=unit.id, hint="give a multiple of the passes")
        d_o = _param(unit, "tube_outer_diameter")
        d_i = d_o - 2.0 * _param(unit, "wall_thickness")
        if d_i <= 0.0:
            raise SimulationError("the tube wall is thicker than the tube radius",
                                  where=unit.id, hint="check the wall thickness")
        wall = fc.tube_wall_resistance(d_i, d_o, p.get("wall_conductivity") or 16.0, d_o)
        tubes = cls(kind, n, passes, d_o, d_i, wall, p.get("tube_correlation") or "gnielinski")
        if kind == "double_pipe":
            tubes.D_a = _param(unit, "annulus_diameter")
            if tubes.D_a <= d_o:
                raise SimulationError("the outer pipe must be wider than the tube",
                                      where=unit.id, hint="the annulus diameter is the outer "
                                      "pipe's inner diameter")
        else:
            pitch = _param(unit, "tube_pitch")
            if pitch <= d_o:
                raise SimulationError("the tube pitch must exceed the tube outer diameter",
                                      where=unit.id)
            tubes.shell = (pitch, p.get("pitch_layout") or "triangular",
                           _param(unit, "shell_diameter"), _param(unit, "baffle_spacing"))
        return tubes

    @property
    def ratio_inside(self) -> float:
        """Outer (reference) over inner tube surface."""
        return self.d_o / self.d_i

    def inside(self, fluid: Fluid, T: float, m: float, heated: bool) -> FilmFn:
        """Film in the tubes for the total tube-side mass flow ``m`` (kg/s)."""
        per_tube = m / (self.n / self.passes)
        G = per_tube / (math.pi / 4.0 * self.d_i ** 2)
        return tube_film(fluid, T, G, self.d_i, self.correlation, heated)

    def outside(self, fluid: Fluid, T: float, m: float, heated: bool,
                notes: list[str]) -> tuple[FilmFn, bool]:
        """Film in the annuli (one per tube) or on the shell side, and whether it depends
        on the wall."""
        if self.kind == "double_pipe":
            G = m / self.n / (math.pi / 4.0 * (self.D_a ** 2 - self.d_o ** 2))
            return annulus_film(fluid, T, G, self.d_o, self.D_a, self.correlation,
                                heated), False
        pitch, layout, D_s, B = self.shell
        return kern_film(fluid, T, m, self.d_o, pitch, layout, D_s, B, notes), True


def check_water(unit, ut) -> None:
    if ut.kind == "fluid" and unit.p.get("h_medium") is None and (
            abs(ut.cp - WATER_CP) > 0.1 * WATER_CP):
        raise SimulationError(
            f"the medium side's film coefficient is computed for water, but {ut.describe()} "
            f"has a heat capacity of {ut.cp / 1e3:.3g} kJ/kgK", where=unit.id,
            hint="give the medium side's film coefficient h_medium")


def utility_side(ut, flow: float | None, Q: float) -> tuple[float, float]:
    """Mean temperature and mass flow (kg/s; 0 for steam) of a utility that passes ``Q``
    into the process, at its fixed flow or at the flow that brings it to its outlet
    temperature."""
    if ut.kind == "steam":
        return ut.T_in, 0.0
    if flow is not None:
        T_out = ut.T_in - Q / (flow * ut.cp) if flow > 0.0 else ut.T_in
        return 0.5 * (ut.T_in + T_out), flow
    m = abs(Q) / (ut.cp * abs(ut.T_in - ut.T_out))
    return 0.5 * (ut.T_in + (ut.T_out if Q != 0.0 else ut.T_in)), m


class FilmRating:
    """U of a unit with a utility (a vessel with a jacket or coil, or a heater) from its
    geometry; see ``U_from``.

    ``overall(pkg, outputs, feed, T_pin, T_pout, Q)`` rates the wall at the process
    temperatures (inlet and outlet; the same for a well-mixed unit) and the duty ``Q`` into
    the process, which sets a liquid utility's flow when it is not fixed. U is referred to
    the area the unit uses (fluxsheet.film_coefficients).
    """

    def __init__(self, unit, ut, flow: float | None):
        self.unit, self.ut, self.flow = unit, ut, flow
        p = unit.p
        self.geometry = p["U_from"]
        self.notes: list[str] = []
        self.last: fc.Overall | None = None
        self.h_given = p.get("h_medium")
        check_water(unit, ut)
        k_w = p.get("wall_conductivity") or 16.0
        t = p.get("wall_thickness") or 0.0
        self.correlation = p.get("tube_correlation") or "gnielinski"
        self.arrangement = p.get("arrangement") or "countercurrent"
        self.tubes: Tubes | None = None
        if self.geometry in VESSEL_GEOMETRIES:
            self.D_T = _param(unit, "vessel_diameter")
            if p.get("process_film") != "boiling":
                self.d_imp = _param(unit, "impeller_diameter")
                self.N = _param(unit, "impeller_speed")
                C, a = fc.STIRRED[self.geometry]
                self.C = p.get("impeller_constant") or C
                self.a = p.get("impeller_exponent") or a
            if self.geometry == "jacket":
                self.wall, self.ratio_medium = t / k_w, 1.0
                if self.h_given is None and ut.kind == "fluid":
                    D_i = self.D_T + 2.0 * t
                    self.jacket = (D_i, D_i + 2.0 * _param(unit, "jacket_gap"))
            else:
                self.d_o = _param(unit, "tube_outer_diameter")
                self.d_i = self.d_o - 2.0 * _param(unit, "wall_thickness")
                if self.d_i <= 0.0:
                    raise SimulationError("the coil's wall is thicker than its radius",
                                          where=unit.id, hint="check the wall thickness")
                self.wall = fc.tube_wall_resistance(self.d_i, self.d_o, k_w, self.d_o)
                self.ratio_medium = self.d_o / self.d_i
            self.ratio_process = 1.0
        elif self.geometry in TUBE_GEOMETRIES:
            self.tubes = Tubes.of(unit)
            self.wall = self.tubes.wall
            self.ratio_process, self.ratio_medium = self.tubes.ratio_inside, 1.0
        else:
            raise SimulationError(f"U from {self.geometry!r} is not available for this unit",
                                  where=unit.id)

    def _note(self, text: str) -> None:
        if text not in self.notes:
            self.notes.append(text)

    def _process(self, pkg, outputs, feed, T_pin: float, T_pout: float,
                 heated: bool) -> tuple[float, FilmFn, bool]:
        unit = self.unit
        if self.tubes is None:
            T, P, kind, zs = contents(outputs, T_pout)
            if unit.p.get("process_film") == "boiling":
                if not heated:
                    raise SimulationError("a boiling process side must be heated",
                                          where=unit.id, hint="set the process side to "
                                          "stirred, or heat the vessel")
                return T, boiling_film(pkg, zs, P), True
            if kind == "vapor":
                self._note("the vessel holds no liquid; the stirred-vessel correlation is "
                           "meant for agitated liquids")
            fluid = Fluid(pkg, kind, tuple(zs), P)
            return T, stirred_film(fluid, T, self.D_T, self.d_imp, self.N, self.C,
                                   self.a), True
        T = 0.5 * (T_pin + T_pout)
        fluid = stream_fluid(pkg, feed.zs, T, feed.P, self.notes, "the process stream")
        return T, self.tubes.inside(fluid, T, mass_flow(pkg, feed), heated), False

    def _medium(self, Q: float) -> tuple[float, FilmFn, bool]:
        ut = self.ut
        T, m = utility_side(ut, self.flow, Q)
        if self.h_given is not None:
            return T, fc.given_film(self.h_given), False
        if ut.kind == "steam":
            if self.geometry == "jacket":  # a vertical wall as high as the vessel is wide
                return T, condensing_film(ut.T_in, ut.latent, self.D_T, "vertical"), True
            if self.geometry == "coil":
                self._note("steam condensing inside the coil is rated with Nusselt's formula "
                           "for the outside of a horizontal tube (an estimate)")
                return T, condensing_film(ut.T_in, ut.latent, self.d_i, "horizontal_tube"), True
            return T, condensing_film(ut.T_in, ut.latent, self.tubes.d_o,
                                      "horizontal_tube"), True
        heated = not ut.heats  # the utility warms up when it cools the process
        w = water()
        if self.geometry == "jacket":
            d_i, D_o = self.jacket
            G = m / (math.pi / 4.0 * (D_o ** 2 - d_i ** 2))
            return T, annulus_film(w, T, G, d_i, D_o, self.correlation, heated), False
        if self.geometry == "coil":
            D_c = self.unit.p.get("coil_diameter")
            factor = fc.coil_factor(self.d_i, D_c) if D_c else 1.0
            G = m / (math.pi / 4.0 * self.d_i ** 2)
            return T, tube_film(w, T, G, self.d_i, self.correlation, heated, factor,
                                "coil"), False
        film, dependent = self.tubes.outside(w, T, m, heated, self.notes)
        return T, film, dependent

    def overall(self, pkg, outputs, feed, T_pin: float, T_pout: float, Q: float) -> fc.Overall:
        heated = self.ut.heats if Q == 0.0 else Q > 0.0
        T_p, f_p, dep_p = self._process(pkg, outputs, feed, T_pin, T_pout, heated)
        T_m, f_m, dep_m = self._medium(Q)
        # The properties are taken at each side's mean temperature, but the mean heat flux
        # is U times the log mean temperature difference: the walls are placed for that.
        ut = self.ut
        if ut.kind == "steam":
            T_u_out = ut.T_in
        elif self.flow is not None:
            T_u_out = ut.T_in - Q / (self.flow * ut.cp) if self.flow > 0.0 else ut.T_in
        else:
            T_u_out = ut.T_out if Q != 0.0 else ut.T_in
        mean = lmtd(*_ends(T_pin, T_pout, ut.T_in, T_u_out, self.arrangement))
        if mean is not None:
            if ut.kind == "steam":
                T_p = T_m - mean
            else:
                T_m = T_p + mean
        p = self.unit.p
        inner = fc.Side(T_p, f_p, self.ratio_process, p.get("fouling_process") or 0.0, dep_p)
        outer = fc.Side(T_m, f_m, self.ratio_medium, p.get("fouling_medium") or 0.0, dep_m)
        try:
            self.last = fc.overall(inner, outer, self.wall)
        except SimulationError as exc:
            raise exc.located(self.unit.id)
        return self.last

    def results(self) -> dict:
        return overall_results(self.last, "process", "medium")


def overall_results(o: fc.Overall | None, inner: str, outer: str) -> dict:
    if o is None:
        return {}
    return {"U": o.U, f"h_{inner}": o.inner.h, f"h_{outer}": o.outer.h,
            f"Re_{inner}": o.inner.Re, f"Nu_{inner}": o.inner.Nu,
            f"Re_{outer}": o.outer.Re, f"Nu_{outer}": o.outer.Nu,
            "T_wall": o.T_wall_inner}


def contents(outputs, T: float) -> tuple[float, float, str, list[float]]:
    """Temperature, pressure, phase and composition of a vessel's agitated contents: its
    liquid if it has one, else its vapor."""
    streams = list(outputs.values())
    n = len(streams[0].flows)
    liquid, vapor = [0.0] * n, [0.0] * n
    for s in streams:
        F, st = s.total_flow, s.state
        if F <= 0.0:
            continue
        x = st.x if st.x is not None else s.zs
        y = st.y if st.y is not None else s.zs
        for i in range(n):
            liquid[i] += F * (1.0 - st.VF) * x[i]
            vapor[i] += F * st.VF * y[i]
    P = min(s.P for s in streams)
    for kind, flows in (("liquid", liquid), ("vapor", vapor)):
        total = sum(flows)
        if total > 0.0:
            return T, P, kind, [f / total for f in flows]
    return T, P, "liquid", list(streams[0].zs)


class TubeWall:
    """U along the wall of a tubular reactor at the local conditions of both sides,
    referred to the inner tube surface pi D (the area the reactor already uses).

    Inside: the reacting stream in an empty tube of diameter D. Outside: a given film
    coefficient, steam condensing on the tube, or (``double_pipe``) the medium flowing in
    the annulus: a liquid utility as water, or the medium stream itself.
    """

    def __init__(self, unit, pkg, feed, ut=None, coolant=None):
        p = unit.p
        self.unit, self.pkg, self.ut = unit, pkg, ut
        self.kind = "vapor" if p["phase"] == "gas" else "liquid"
        self.D = p["diameter"]
        t = p.get("wall_thickness") or 0.0
        self.d_o = self.D + 2.0 * t
        k_w = p.get("wall_conductivity") or 16.0
        self.wall = fc.tube_wall_resistance(self.D, self.d_o, k_w, self.D) if t else 0.0
        self.ratio_medium = self.D / self.d_o
        self.correlation = p.get("tube_correlation") or "gnielinski"
        self.G = mass_flow(pkg, feed) / (math.pi / 4.0 * self.D ** 2)
        self.h_given = p.get("h_medium")
        self.notes: list[str] = []
        self.annulus = p["U_from"] == "double_pipe"
        self.medium_fluid = None
        if self.annulus:
            self.D_a = _param(unit, "annulus_diameter")
            if self.D_a <= self.d_o:
                raise SimulationError("the outer pipe must be wider than the tube",
                                      where=unit.id, hint="the annulus diameter is the outer "
                                      "pipe's inner diameter; the tube's outer diameter is "
                                      "D plus twice the wall thickness")
        if self.h_given is not None:
            return
        if coolant is not None:
            if not self.annulus:
                raise SimulationError(
                    "the medium stream's film coefficient needs the annulus it flows in",
                    where=unit.id, hint="choose double_pipe with the annulus diameter, or "
                    "give the medium side's film coefficient h_medium")
            st = coolant.state
            kind = "vapor" if st.VF >= 0.5 else "liquid"
            if 0.0 < st.VF < 1.0:
                self.notes.append(f"the medium stream enters two-phase; its film coefficient "
                                  f"is for single-phase {kind} flow")
            zs = (st.y if kind == "vapor" else st.x) or coolant.zs
            self.medium_fluid = Fluid(pkg, kind, tuple(zs), coolant.P)
            self.medium_flow = mass_flow(pkg, coolant)
        elif ut.kind == "fluid":
            check_water(unit, ut)
            if not self.annulus:
                raise SimulationError(
                    "the utility's film coefficient needs the annulus it flows in",
                    where=unit.id, hint="choose double_pipe with the annulus diameter, or "
                    "give the medium side's film coefficient h_medium")
            self.medium_fluid = water()

    def overall(self, T: float, P: float, flows, T_m: float, medium_flow: float | None = None
                ) -> fc.Overall:
        """``medium_flow``: kg/s of a liquid utility (a stream medium has its own)."""
        total = sum(flows)
        zs = [max(f, 0.0) / total for f in flows]
        process = Fluid(self.pkg, self.kind, tuple(zs), P)
        inner = fc.Side(T, tube_film(process, T, self.G, self.D, self.correlation, T_m > T),
                        1.0, self.unit.p.get("fouling_process") or 0.0)
        heated = T_m < T  # the medium warms up
        dependent = False
        if self.h_given is not None:
            film = fc.given_film(self.h_given)
        elif self.medium_fluid is None:  # steam
            film, dependent = condensing_film(self.ut.T_in, self.ut.latent, self.d_o,
                                              "horizontal_tube"), True
        else:
            m = self.medium_flow if medium_flow is None else medium_flow
            G = m / (math.pi / 4.0 * (self.D_a ** 2 - self.d_o ** 2))
            film = annulus_film(self.medium_fluid, T_m, G, self.d_o, self.D_a,
                                self.correlation, heated)
        outer = fc.Side(T_m, film, self.ratio_medium,
                        self.unit.p.get("fouling_medium") or 0.0, dependent)
        try:
            return fc.overall(inner, outer, self.wall)
        except SimulationError as exc:
            raise exc.located(self.unit.id)
