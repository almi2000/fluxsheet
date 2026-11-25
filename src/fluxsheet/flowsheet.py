"""The flowsheet file format.

A flowsheet is plain JSON, meant to be read, diffed and versioned like source code::

    {
      "format": "fluxsheet/1",
      "name": "Ethanol flash",
      "components": ["water", "ethanol"],
      "property_method": "auto",
      "units": [
        {"id": "FEED", "type": "Feed",
         "params": {"T": "25 degC", "P": "1 atm", "flows": {"water": "60 kmol/h", ...}}},
        {"id": "H1", "type": "Heater", "params": {"T": "85 degC"}},
        {"id": "C1", "type": "DistillationColumn", "property_method": "nrtl", ...},
        ...
      ],
      "streams": [
        {"id": "S1", "from": "FEED.out", "to": "H1.in"},
        ...
      ]
    }

Values with units are optional: plain numbers are read as SI (K, Pa, mol/s, W).
On a recycle line, ``circulation`` sets the steady mass flow and, when needed,
its mole fractions. ``recycle_guess`` supplies starting values for the tear iteration.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    FiniteFloat,
    PositiveFloat,
    PositiveInt,
    field_validator,
    model_validator,
)

from .dynamics.models import DynamicsModel
from .study_models import (
    CalculatorModel,
    CostingModel,
    DesignSpecModel,
    OptimizationModel,
    SensitivityModel,
)

FORMAT = "fluxsheet/1"


class UnitModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str
    type: str
    params: dict[str, Any] = Field(default_factory=dict)
    # Overrides the flowsheet's property_method for this unit; see property_sections.py.
    property_method: str | None = None
    ui: dict[str, Any] | None = None  # editor layout (position etc.); ignored by the solver

    @field_validator("id")
    @classmethod
    def _valid_id(cls, v: str) -> str:
        if not v or "." in v or v != v.strip():
            raise ValueError(f"unit id {v!r} must be non-empty, without dots or edge spaces")
        return v


class CirculationModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    mass_flow: float | str
    mole_fractions: dict[str, float] | None = None


class RecycleGuessModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    T: float | str | None = None
    P: float | str | None = None
    flows: dict[str, float | str] | None = None


class StreamModel(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    id: str
    source: str = Field(alias="from")
    target: str | None = Field(default=None, alias="to")
    circulation: CirculationModel | None = None
    recycle_guess: RecycleGuessModel | None = None

    @field_validator("source", "target")
    @classmethod
    def _endpoint(cls, v: str | None) -> str | None:
        if v is not None and v.count(".") != 1:
            raise ValueError(f"stream endpoint {v!r} must look like 'UNIT.port'")
        return v


class SolverOptions(BaseModel):
    model_config = ConfigDict(extra="forbid")

    max_iterations: int = 100
    tolerance: float = 1e-6
    acceleration: Literal["wegstein", "broyden", "none"] = "wegstein"


class RateLawModel(BaseModel):
    """Power law ``r = k0 exp(-Ea/RT) prod(c_i^order_i)``; see reactions.py for units."""

    model_config = ConfigDict(extra="forbid")

    k0: float
    Ea: float | str = 0.0
    orders: dict[str, float] | None = None  # default: the stoichiometric coefficients
    T_ref: float | str | None = None  # if set, k0 is the rate constant at T_ref


class KineticsModel(RateLawModel):
    basis: str = "concentration"  # or "partial_pressure"
    reverse: RateLawModel | None = None
    equilibrium_limited: bool = False


class KeqModel(BaseModel):
    """User equilibrium constant: ln K = A + B/T + C ln(T) + D T (T in K)."""

    model_config = ConfigDict(extra="forbid")

    A: float = 0.0
    B: float = 0.0
    C: float = 0.0
    D: float = 0.0


class ReactionModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str
    stoichiometry: dict[str, float]  # negative for reactants, positive for products
    kinetics: KineticsModel | None = None
    keq: KeqModel | None = None


class UtilityModel(BaseModel):
    """A heating or cooling medium; see utilities.py."""

    model_config = ConfigDict(extra="forbid")

    id: str
    type: str = "fluid"  # or "steam"
    T_in: float | str | None = None  # fluid
    T_out: float | str | None = None  # fluid
    cp: float | str | None = None  # fluid, default that of water
    P: float | str | None = None  # steam, saturated


class CustomCompoundModel(BaseModel):
    """A user-defined compound; missing constants are estimated (custom_compounds.py)."""

    model_config = ConfigDict(extra="forbid")

    name: str
    formula: str
    smiles: str | None = None
    joback_groups: dict[str, int] | None = None  # group name or number -> count
    unifac_groups: dict[str, int] | None = None  # UNIFAC (Dortmund) subgroups
    constants: dict[str, float | str] | None = None  # Tc, Pc, omega, Tb, MW, Hf, Gf, ...
    vapor_pressure: str | None = None  # "lee-kesler" (default) or "ambrose-walton"
    antoine: dict[str, float | str] | None = None  # A, B, C: log10(P/Pa) = A - B/(T/K + C)


class BinaryParameterModel(BaseModel):
    """Directional coefficients in component order, using the thermo model equations."""

    model_config = ConfigDict(extra="forbid", frozen=True, str_strip_whitespace=True,
                              allow_inf_nan=False)

    method: Literal["nrtl", "uniquac", "wilson"]
    components: tuple[str, str]
    convention: Literal["thermo"]
    source: str = Field(min_length=1)
    T_range_K: tuple[PositiveFloat, PositiveFloat]
    P_range_Pa: tuple[PositiveFloat, PositiveFloat] | None = None
    a: tuple[FiniteFloat, FiniteFloat] = (0.0, 0.0)
    b: tuple[FiniteFloat, FiniteFloat] = (0.0, 0.0)
    alpha: FiniteFloat | None = Field(default=None, ge=0)

    @model_validator(mode="after")
    def _check_parameters(self) -> BinaryParameterModel:
        if not all(self.components) or self.components[0] == self.components[1]:
            raise ValueError("binary parameters require two different component names")
        if (self.method == "nrtl") != (self.alpha is not None):
            raise ValueError("alpha is required for NRTL and only applies to NRTL")
        if not self.model_fields_set.intersection({"a", "b"}):
            raise ValueError("supply at least one of a or b for the binary pair")
        for limits in (self.T_range_K, self.P_range_Pa):
            if limits is not None and limits[0] > limits[1]:
                raise ValueError("validity range must be ordered as [minimum, maximum]")
        return self


class GroupAssignmentModel(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    component: str = Field(min_length=1)
    convention: Literal["dortmund"]
    source: str = Field(min_length=1)
    groups: dict[PositiveInt, PositiveInt] = Field(min_length=1)


class PropertyDataModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    binary_parameters: list[BinaryParameterModel] = Field(default_factory=list)
    unifac_groups: list[GroupAssignmentModel] = Field(default_factory=list)
    strict: bool = False


class FlowsheetModel(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    format: str = FORMAT
    name: str = "Untitled flowsheet"
    description: str = ""
    components: list[str]
    property_method: str = "auto"
    # Dissolved gases of activity methods; None = automatic, see henry.py.
    henry_components: list[str] | None = None
    custom_compounds: list[CustomCompoundModel] | None = None  # see custom_compounds.py
    property_data: PropertyDataModel | None = None
    units: list[UnitModel]
    streams: list[StreamModel] = Field(default_factory=list)
    reactions: list[ReactionModel] | None = None  # referenced by id from reactor units
    utilities: list[UtilityModel] | None = None  # referenced by id from units; utilities.py
    solver: SolverOptions = Field(default_factory=SolverOptions)
    design_specs: list[DesignSpecModel] | None = None  # solved around the flowsheet; specs.py
    sensitivity: list[SensitivityModel] | None = None  # run on request; sensitivity.py
    calculators: list[CalculatorModel] | None = None  # set parameters from formulas; calculator.py
    dynamics: DynamicsModel | None = None  # run on request; dynamics/
    costing: CostingModel | None = None  # equipment and utility costs; costing.py
    optimizations: list[OptimizationModel] | None = None  # run on request; optimization.py

    @model_validator(mode="after")
    def _check_property_data(self) -> FlowsheetModel:
        if self.property_data is None:
            return self
        components = set(self.components)
        assignments = set()
        for entry in self.property_data.unifac_groups:
            if entry.component not in components or entry.component in assignments:
                raise ValueError(f"UNIFAC assignment for unknown or repeated {entry.component!r}")
            assignments.add(entry.component)
        pairs = set()
        for entry in self.property_data.binary_parameters:
            pair = (entry.method, frozenset(entry.components))
            if not pair[1] <= components or pair in pairs:
                raise ValueError(f"{entry.method} parameters for unknown or repeated pair "
                                 f"{' / '.join(entry.components)}")
            pairs.add(pair)
        return self

    @classmethod
    def load(cls, path: str | Path) -> FlowsheetModel:
        return cls.model_validate(json.loads(Path(path).read_text(encoding="utf-8")))

    def dump(self, path: str | Path) -> None:
        Path(path).write_text(self.to_json() + "\n", encoding="utf-8")

    def to_json(self) -> str:
        return self.model_dump_json(by_alias=True, exclude_none=True, indent=2)
