"""Flashes of a vessel's contents along a time trajectory.

Shared by the dynamic hold-up models (dynamics/units.py) and the batch reactor
(units/batch.py), which both integrate the contents of a vessel over time.
"""

from __future__ import annotations

import numpy as np

from .thermo_package import PropertyPackage


class ContentsFlash:
    """PH flashes along a trajectory, which change little from one call to the next.

    A full flash (with the phase stability test) costs milliseconds and an integration needs
    thousands. While the contents stay one phase, near a temperature and composition where a
    full flash confirmed that phase, the temperature is found on that phase's enthalpy alone
    by Newton steps; the full flash runs again once they drift further.
    """

    T_WINDOW = 2.0  # K around the last confirmed state
    Z_WINDOW = 0.01  # mole fraction

    def __init__(self, pkg: PropertyPackage, T: float):
        self.pkg, self.T = pkg, T
        self.kind: str | None = None  # single phase confirmed by the last full flash
        self.checked: tuple[float, np.ndarray] | None = None

    def __call__(self, zs, P: float, H: float):
        zs = np.asarray(zs, dtype=float)
        if self.kind is not None and self._near(self.T, zs):
            T = self._single_phase_T(zs, P, H)
            if T is not None and self._near(T, zs):
                self.T = T
                return self.pkg.phase_state(self.kind, T, P, list(zs))
        st = self.pkg.flash_ph_near(list(zs), P, H, self.T)
        self.T = st.T
        single = st.phase in ("liquid", "vapor") and not st.liquid_split
        self.kind = st.phase if single else None
        self.checked = (st.T, zs.copy())
        return st

    def _near(self, T: float, zs: np.ndarray) -> bool:
        T0, z0 = self.checked
        return abs(T - T0) < self.T_WINDOW and float(np.abs(zs - z0).max()) < self.Z_WINDOW

    def _single_phase_T(self, zs: np.ndarray, P: float, H: float) -> float | None:
        pkg, kind, T = self.pkg, self.kind, self.T
        z = list(zs)
        h0 = pkg.phase_props(kind, T, P, z).H
        cp = (pkg.phase_props(kind, T + 0.01, P, z).H - h0) / 0.01
        for _ in range(8):
            if cp <= 0.0:
                return None
            step = (H - h0) / cp
            T += step
            if abs(step) < 1e-9 * T:
                return T
            h0 = pkg.phase_props(kind, T, P, z).H
        return None

    def at_temperature(self, zs, P: float, T: float):
        """T-P flash of the contents, with the same single-phase shortcut: near a state where
        a full flash confirmed one phase, that phase is evaluated directly."""
        zs = np.asarray(zs, dtype=float)
        if self.kind is not None and self._near(T, zs):
            return self.pkg.phase_state(self.kind, T, P, list(zs))
        st = self.pkg.flash(list(zs), T=T, P=P)
        self.T = st.T
        single = st.phase in ("liquid", "vapor") and not st.liquid_split
        self.kind = st.phase if single else None
        self.checked = (st.T, zs.copy())
        return st
