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
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

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


class StreamModel(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    id: str
    source: str = Field(alias="from")
    target: str | None = Field(default=None, alias="to")

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


class FlowsheetModel(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    format: str = FORMAT
    name: str = "Untitled flowsheet"
    description: str = ""
    components: list[str]
    property_method: str = "auto"
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

    @classmethod
    def load(cls, path: str | Path) -> FlowsheetModel:
        return cls.model_validate(json.loads(Path(path).read_text(encoding="utf-8")))

    def dump(self, path: str | Path) -> None:
        Path(path).write_text(self.to_json() + "\n", encoding="utf-8")

    def to_json(self) -> str:
        return self.model_dump_json(by_alias=True, exclude_none=True, indent=2)
