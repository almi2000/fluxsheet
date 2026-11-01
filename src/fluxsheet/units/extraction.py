"""Liquid-liquid extraction column.

Countercurrent equilibrium stages at constant temperature: the feed enters stage 1, the
solvent the last stage. On every stage the two incoming liquids settle into two liquid
phases in equilibrium (vapor-liquid-liquid flash); the phase richer in the solvent is the
extract and flows towards stage 1, the other is the raffinate and flows on towards the last
stage. The stages are solved by sweeps down and up the column (symmetric Gauss-Seidel),
each stage being a decanter, until the flows stop changing.
"""

from __future__ import annotations

import numpy as np

from ..errors import SimulationError
from ..stream import Stream
from .base import DUTY, Param, Port, ResultField, UnitOperation

MAX_SWEEPS = 300


class ExtractionColumn(UnitOperation):
    type = "ExtractionColumn"
    label = "Extraction column"
    category = "Columns"
    description = "Countercurrent liquid-liquid extraction on equilibrium stages at constant "\
        "temperature: the feed enters stage 1, the solvent the last stage; the extract "\
        "leaves at the feed end, the raffinate at the solvent end. Needs a property method "\
        "that predicts two liquid phases (NRTL, UNIQUAC or UNIFAC)."
    inlets = (Port("feed", "Feed (stage 1)", multiple=True),
              Port("solvent", "Solvent (last stage)", multiple=True))
    outlets = (Port("extract", "Extract"), Port("raffinate", "Raffinate"))
    params = (
        Param("stages", "Number of stages", "integer", default=5, min=1, max=60),
        Param("T", "Temperature", "temperature", min=0.0,
              help="Of every stage. Defaults to the feed temperature."),
        Param("P", "Pressure", "pressure", min=1.0,
              help="Defaults to the lower of the feed and solvent pressures."),
        Param("solvent_component", "Solvent component", "component",
              help="The extract is the liquid richer in it. Defaults to the main component of "
              "the solvent stream."),
    )
    result_fields = (DUTY, ResultField("sweeps", "Sweeps to converge"))

    def calculate(self, pkg, inlets):
        if not pkg.supports_lle:
            raise SimulationError(
                f"the property method {pkg.method!r} cannot predict two liquid phases",
                where=self.id, hint="set the property method to NRTL, UNIQUAC or UNIFAC")
        feed, solvent = inlets["feed"][0], inlets["solvent"][0]
        T = self.p["T"] or feed.T
        P = self.p["P"] or min(feed.P, solvent.P)
        n = self.p["stages"]
        if feed.total_flow <= 0.0 or solvent.total_flow <= 0.0:
            raise SimulationError("the extraction column needs both a feed and a solvent flow",
                                  where=self.id)
        key = (self.component_index(pkg, self.p["solvent_component"], "solvent component")
               if self.p["solvent_component"] else int(np.argmax(solvent.flows)))
        f, s = np.array(feed.flows), np.array(solvent.flows)
        scale = f.sum() + s.sum()
        R = [f.copy() for _ in range(n)]  # raffinate leaving each stage (towards stage n)
        E = [s.copy() for _ in range(n)]  # extract leaving each stage (towards stage 1)
        split = [True] * n

        def stage(j: int) -> None:
            inflow = (R[j - 1] if j else f) + (E[j + 1] if j + 1 < n else s)
            R[j], E[j], split[j] = self._settle(pkg, inflow, T, P, key, j + 1)

        sweeps = 0
        for sweeps in range(1, MAX_SWEEPS + 1):
            before = np.concatenate(E + R)
            for j in range(n):
                stage(j)
            for j in reversed(range(n)):
                stage(j)
            if np.abs(np.concatenate(E + R) - before).max() < 1e-10 * scale:
                break
        else:
            raise SimulationError(f"the stages did not converge in {MAX_SWEEPS} sweeps",
                                  where=self.id, hint="try fewer stages or more solvent")
        if not all(split):
            bad = ", ".join(str(j + 1) for j in range(n) if not split[j])
            raise SimulationError(
                f"no second liquid phase on stage {bad}: feed and solvent mix completely",
                where=self.id, hint="use more solvent, a less soluble solvent or another "
                "temperature; or the property method does not predict the liquid-liquid "
                "split of these components (compare with UNIFAC)")
        extract = Stream.from_flash(pkg, list(E[0]), T=T, P=P)
        raffinate = Stream.from_flash(pkg, list(R[-1]), T=T, P=P)
        h_in = sum(mix.enthalpy_flow for mix in (feed, solvent))
        self.results.update(duty=extract.enthalpy_flow + raffinate.enthalpy_flow - h_in,
                            sweeps=sweeps)
        present = [i for i in range(pkg.n) if f[i] + s[i] > 0.0]
        stages = list(range(1, n + 1))
        self.add_profile(
            "stages", "Stage profile", "Stage", stages,
            {"Raffinate flow": ("molar_flow", [float(r.sum()) for r in R]),
             "Extract flow": ("molar_flow", [float(e.sum()) for e in E]),
             **{f"Raffinate x {pkg.names[i]}": ("fraction", [float(r[i] / r.sum()) for r in R])
                for i in present}},
        )
        return {"extract": extract, "raffinate": raffinate}

    def _settle(self, pkg, flows: np.ndarray, T: float, P: float, key: int, stage: int):
        """Raffinate and extract flows of one stage, and whether two liquids formed."""
        total = flows.sum()
        st = pkg.flash(list(flows / total), T=T, P=P)
        if st.VF > 1e-9:
            raise SimulationError(f"liquid boils on stage {stage}", where=self.id,
                                  hint="lower the temperature or raise the pressure")
        if not st.liquid_split:
            # One liquid: it goes with the solvent if it is mostly solvent.
            mostly_solvent = flows[key] / total > 0.5
            empty = np.zeros_like(flows)
            return (empty, flows.copy(), False) if mostly_solvent else (flows.copy(), empty, False)
        light = total * st.light_fraction * np.array(st.x_light)
        heavy = np.maximum(flows - light, 0.0)
        light_is_extract = st.x_light[key] >= st.x_heavy[key]
        return (heavy, light, True) if light_is_extract else (light, heavy, True)

