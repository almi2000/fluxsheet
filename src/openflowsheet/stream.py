"""Material streams: component flows plus the thermodynamic state they are in."""

from __future__ import annotations

from dataclasses import dataclass

from .thermo_package import PhaseState, PropertyPackage


@dataclass
class Stream:
    flows: list[float]  # mol/s per component
    state: PhaseState

    @property
    def total_flow(self) -> float:
        return sum(self.flows)

    @property
    def zs(self) -> list[float]:
        total = self.total_flow
        if total <= 0.0:
            return [1.0 / len(self.flows)] * len(self.flows)
        return [f / total for f in self.flows]

    @property
    def T(self) -> float:
        return self.state.T

    @property
    def P(self) -> float:
        return self.state.P

    @property
    def enthalpy_flow(self) -> float:
        """Total enthalpy flow, W."""
        return self.state.H * self.total_flow

    @classmethod
    def from_flash(cls, pkg: PropertyPackage, flows: list[float], **spec: float) -> Stream:
        flows = [max(f, 0.0) for f in flows]
        total = sum(flows)
        zs = [f / total for f in flows] if total > 0 else [1.0 / len(flows)] * len(flows)
        return cls(flows=flows, state=pkg.flash(zs, **spec))

    def to_dict(self, pkg: PropertyPackage) -> dict:
        mws = pkg.mws
        st = self.state
        return {
            "T": st.T,
            "P": st.P,
            "vapor_fraction": st.VF,
            "phase": st.phase,
            "molar_flow": self.total_flow,
            "mass_flow": sum(f * mw for f, mw in zip(self.flows, mws)) / 1000.0,  # kg/s
            "enthalpy": st.H,
            "enthalpy_flow": self.enthalpy_flow,
            "component_flows": dict(zip(pkg.names, self.flows)),
            "mole_fractions": dict(zip(pkg.names, self.zs)),
            "liquid_composition": dict(zip(pkg.names, st.x)) if st.x else None,
            "vapor_composition": dict(zip(pkg.names, st.y)) if st.y else None,
        }
