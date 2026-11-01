"""Evaporators.

The thin-film evaporator (wiped or falling film) runs a liquid film down a heated wall; the
vapor is drawn off as it forms, with little contact with the liquid further down. It is
modeled in plug flow along the heated length: the length is cut into segments, each adds
its share of the duty and flashes off what boils, and the vapor leaves the film. The film
first heats up to its bubble point, then evaporates while it gets richer in the heavy
components and hotter. This is a differential (Rayleigh) evaporation in space, so it
separates better than a single flash at the same fraction evaporated.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from ..errors import SimulationError
from ..stream import Stream
from .base import DUTY, IN, Param, Port, ResultField, UnitOperation, empty_stream, single_inlet
from .heat_transfer import AREA, HEAT_PARAMS, HEAT_RESULTS


class ThinFilmEvaporator(UnitOperation):
    type = "ThinFilmEvaporator"
    label = "Thin-film evaporator"
    category = "Separators"
    description = "Wiped or falling film evaporator: the liquid runs down a heated wall as a "\
        "thin film and the vapor is drawn off as it forms. Calculated in plug flow along the "\
        "heated length, so it separates better than a single flash. Usual under vacuum for "\
        "heat-sensitive or high-boiling liquids."
    inlets = (IN,)
    outlets = (Port("vapor", "Vapor"), Port("residue", "Residue"))
    params = (
        Param("VF", "Fraction evaporated", "fraction", min=0.0, max=0.999,
              help="Moles of vapor per mole of feed; the duty follows."),
        Param("duty", "Duty", "power", min=0.0, help="Heat put into the film."),
        AREA,
        *HEAT_PARAMS,
        Param("P", "Pressure", "pressure", min=1.0,
              help="Operating pressure, often a vacuum. Defaults to the feed pressure."),
        Param("segments", "Segments", "integer", default=20, min=2, max=400,
              help="Steps along the heated length; more give a finer profile."),
    )
    exclusive_groups = (("VF", "duty", "area"),)
    required_groups = exclusive_groups
    process_side = "plug"
    result_fields = (
        DUTY,
        ResultField("vapor_fraction", "Fraction evaporated", "fraction"),
        ResultField("T_residue", "Residue temperature", "temperature"),
        *HEAT_RESULTS,
    )

    def calculate(self, pkg, inlets):
        s = single_inlet(inlets)
        P = self.p["P"] or s.P
        if s.total_flow <= 0.0:
            self.results.update(duty=0.0, vapor_fraction=0.0, T_residue=s.T)
            return {"vapor": empty_stream(pkg, s.T, P), "residue": empty_stream(pkg, s.T, P)}
        if self.p["VF"] is not None:
            Q = self._duty_for(pkg, s, P, self.p["VF"])
        else:
            Q = self.p["duty"] or 0.0
        return self._run(pkg, s, P, Q, record=True)

    def _evaporate(self, pkg, s: Stream, P: float, Q: float):
        """Film along the wall: liquid and vapor flows, enthalpies and the profile."""
        n = self.p["segments"]
        liquid = np.array(s.flows, dtype=float)
        H_liquid = s.enthalpy_flow  # W
        vapor = np.zeros(pkg.n)
        H_vapor = 0.0
        F = s.total_flow
        z, T, evaporated, xs = [], [], [], []
        T_film = s.T
        # Segment 0 is the inlet: a feed above its bubble point flashes as it enters.
        for k in range(n + 1):
            add = Q / n if k else 0.0
            total = liquid.sum()
            if total <= 1e-12 * F:  # all evaporated: the rest superheats the vapor
                H_vapor += add
                continue
            st = pkg.flash_ph_near(list(liquid / total), P, (H_liquid + add) / total, T_film)
            T_film = st.T
            if st.VF > 0.0 and st.y:
                v = total * st.VF * np.array(st.y)
                hv = pkg.phase_props("vapor", st.T, P, list(st.y)).H if st.VF < 1.0 else st.H
                vapor += v
                H_vapor += v.sum() * hv
                liquid = np.maximum(liquid - v, 0.0)
                H_liquid = H_liquid + add - v.sum() * hv
            else:
                H_liquid += add
            z.append(k / n)
            T.append(st.T)
            evaporated.append(vapor.sum() / F)
            total = liquid.sum()
            xs.append(liquid / total if total > 0 else np.zeros(pkg.n))
        return liquid, H_liquid, vapor, H_vapor, (z, T, evaporated, xs)

    def _duty_for(self, pkg, s: Stream, P: float, VF: float) -> float:
        """Duty that evaporates the fraction ``VF``: secant steps from the duty of a single
        flash (a film needs a little more, its residue ends hotter), then Brent's method."""
        def short(Q: float) -> float:
            return self._evaporate(pkg, s, P, Q)[2].sum() / s.total_flow - VF

        flashed = short(0.0)
        if flashed > 0.0:
            raise SimulationError(
                f"the feed already flashes {flashed + VF:.3g} of itself at "
                f"{P / 1e5:.4g} bar, more than the {VF:.3g} asked for", where=self.id,
                hint="raise the pressure, cool the feed, or ask for more evaporation")
        Q0 = s.total_flow * (pkg.flash(s.zs, P=P, VF=VF).H - s.state.H)
        Q1 = 1.05 * Q0
        f0, f1 = short(Q0), short(Q1)
        for _ in range(12):
            if abs(f1) < 1e-7:
                return Q1
            if f1 == f0:
                break
            Q0, f0, Q1 = Q1, f1, Q1 - f1 * (Q1 - Q0) / (f1 - f0)
            if Q1 <= 0.0:
                break
            f1 = short(Q1)
        # Enough heat to reach the dew point of the whole feed, then more if the residue
        # still has to get hotter than that.
        Q_hi = max(s.total_flow * (pkg.flash(s.zs, P=P, VF=1.0).H - s.state.H), 1.0)
        while short(Q_hi) < 0.0:
            Q_hi *= 2.0
            if Q_hi > 1e12:
                raise SimulationError("could not reach the fraction evaporated", where=self.id,
                                      hint="ask for less evaporation")
        return brentq(short, 0.0, Q_hi, xtol=1e-9 * Q_hi, rtol=1e-12)

    def _run(self, pkg, s: Stream, P: float, Q: float, record: bool):
        liquid, H_liquid, vapor, H_vapor, (z, T, evaporated, xs) = self._evaporate(pkg, s, P, Q)
        F = s.total_flow
        T_end = T[-1] if T else s.T
        residue = (Stream.from_flash(pkg, list(liquid), P=P, H=H_liquid / liquid.sum())
                   if liquid.sum() > 1e-12 * F else empty_stream(pkg, T_end, P))
        out_vapor = (Stream.from_flash(pkg, list(vapor), P=P, H=H_vapor / vapor.sum())
                     if vapor.sum() > 1e-12 * F else empty_stream(pkg, T_end, P))
        self.results.update(duty=float(Q), vapor_fraction=float(vapor.sum() / F),
                            T_residue=residue.T)
        if record:
            present = [i for i in range(pkg.n) if s.flows[i] > 0.0]
            self.add_profile(
                "film", "Along the heated wall", "Heated length (fraction)", z,
                {"Film temperature": ("temperature", T),
                 "Evaporated": ("fraction", evaporated),
                 **{f"x {pkg.names[i]}": ("fraction", [x[i] for x in xs]) for i in present}},
            )
        return {"vapor": out_vapor, "residue": residue}
