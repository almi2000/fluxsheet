"""File format of design specifications and sensitivity studies.

Both refer to flowsheet values in the same way::

    "design_specs": [
      {"id": "PURITY",
       "target": {"stream": "DIST", "property": "mole_fraction", "component": "benzene"},
       "value": 0.995,
       "vary": {"unit": "COLUMN", "param": "reflux_ratio", "lower": 1.2, "upper": 5}}
    ],
    "sensitivity": [
      {"id": "REFLUX",
       "vary": [{"unit": "COLUMN", "param": "reflux_ratio", "lower": 1.2, "upper": 4,
                 "points": 9}],
       "outputs": [{"stream": "DIST", "property": "mole_fraction", "component": "benzene"},
                   {"unit": "COLUMN", "property": "reboiler_duty"}]}
    ]

A target or output is a stream property (see ``specs.STREAM_PROPERTIES``) or a scalar unit
result such as a duty. Bounds and values take units like any parameter ("10 kmol/h").
See specs.py and sensitivity.py for how they are solved.
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field, model_validator


class VariableRef(BaseModel):
    """A computed value: ``property`` of a stream, or the result ``property`` of a unit."""

    model_config = ConfigDict(extra="forbid")

    stream: str | None = None
    unit: str | None = None
    property: str
    component: str | None = None  # for mole_fraction and component_flow

    @model_validator(mode="after")
    def _one_owner(self) -> VariableRef:
        if (self.stream is None) == (self.unit is None):
            raise ValueError("a variable belongs to exactly one of 'stream' or 'unit'")
        return self


class VariedParam(BaseModel):
    """A unit parameter that a design spec or sensitivity study changes between bounds."""

    model_config = ConfigDict(extra="forbid")

    unit: str
    param: str
    component: str | None = None  # one entry of a per-component parameter, e.g. Feed flows
    lower: float | str
    upper: float | str


class DesignSpecModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str
    target: VariableRef
    value: float | str
    vary: VariedParam
    tolerance: float | str | None = None  # absolute, in the target's units; default relative
    max_iterations: int = Field(default=30, ge=2, le=200)
    active: bool = True


class SensitivityVary(VariedParam):
    points: int = Field(default=11, ge=2, le=101)


class SensitivityModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str
    vary: list[SensitivityVary] = Field(min_length=1, max_length=2)
    outputs: list[VariableRef] = Field(min_length=1)
