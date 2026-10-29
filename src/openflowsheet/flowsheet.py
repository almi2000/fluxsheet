"""The flowsheet file format.

A flowsheet is plain JSON, meant to be read, diffed and versioned like source code::

    {
      "format": "openflowsheet/1",
      "name": "Ethanol flash",
      "components": ["water", "ethanol"],
      "property_method": "auto",
      "units": [
        {"id": "FEED", "type": "Feed",
         "params": {"T": "25 degC", "P": "1 atm", "flows": {"water": "60 kmol/h", ...}}},
        {"id": "H1", "type": "Heater", "params": {"T": "85 degC"}},
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
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator

FORMAT = "openflowsheet/1"


class UnitModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str
    type: str
    params: dict[str, Any] = Field(default_factory=dict)
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
    acceleration: str = "wegstein"  # or "none"


class FlowsheetModel(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    format: str = FORMAT
    name: str = "Untitled flowsheet"
    description: str = ""
    components: list[str]
    property_method: str = "auto"
    units: list[UnitModel]
    streams: list[StreamModel] = Field(default_factory=list)
    solver: SolverOptions = Field(default_factory=SolverOptions)

    @classmethod
    def load(cls, path: str | Path) -> FlowsheetModel:
        return cls.model_validate(json.loads(Path(path).read_text(encoding="utf-8")))

    def dump(self, path: str | Path) -> None:
        Path(path).write_text(self.to_json() + "\n", encoding="utf-8")

    def to_json(self) -> str:
        return self.model_dump_json(by_alias=True, exclude_none=True, indent=2)
