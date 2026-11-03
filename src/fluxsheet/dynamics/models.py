"""File format of dynamic simulations, under ``"dynamics"`` in a flowsheet::

    "dynamics": {
      "t_end": "2 h", "points": 241,
      "units": {"TANK": {"diameter": "1.5 m", "height": "3 m", "initial_level": 0.5}},
      "controllers": [
        {"id": "LC1", "measure": {"unit": "TANK", "property": "level"},
         "manipulate": {"unit": "TANK", "param": "liquid_flow"},
         "action": "direct", "gain": 2, "integral_time": "5 min",
         "op_range": ["0 kmol/h", "200 kmol/h"]}
      ],
      "events": [
        {"time": "10 min", "unit": "FEED", "param": "flows", "component": "water",
         "value": "120 kmol/h"},
        {"time": "1 h", "controller": "LC1", "value": "1.8 m", "ramp": "5 min"}
      ]
    }

``units`` holds the dynamic settings of units with a hold-up (sizes, initial level, time
constants; see ``dynamics.units``). Controllers and events address parameters like design
specs do. Values take units; plain numbers are SI.
"""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from ..study_models import VariableRef


class ManipulatedParam(BaseModel):
    """A unit parameter or dynamic setting that a controller or event changes."""

    model_config = ConfigDict(extra="forbid")

    unit: str
    param: str
    component: str | None = None


class ControllerModel(BaseModel):
    """A PID controller. Gain is in % of output range per % of measurement range."""

    model_config = ConfigDict(extra="forbid")

    id: str
    measure: VariableRef
    manipulate: ManipulatedParam
    setpoint: float | str | None = None  # default: the measured steady-state value
    # "reverse": output rises when the measurement falls below the setpoint (a heater on a
    # temperature); "direct": it rises with the measurement (an outlet valve on a level).
    action: Literal["reverse", "direct"] = "reverse"
    gain: float = Field(default=1.0, ge=0.0)
    integral_time: float | str | None = None  # None: no integral action
    derivative_time: float | str | None = None
    measurement_lag: float | str = 1.0  # first-order sensor lag, s
    op_range: tuple[float | str, float | str]
    pv_range: tuple[float | str, float | str] | None = None  # default: setpoint ± 50 %
    active: bool = True


class EventModel(BaseModel):
    """At ``time``, set a unit parameter (or a controller's setpoint) to ``value``, at once or
    as a linear ramp over ``ramp``."""

    model_config = ConfigDict(extra="forbid")

    time: float | str
    unit: str | None = None
    param: str | None = None
    component: str | None = None
    controller: str | None = None
    value: float | str
    ramp: float | str | None = None

    @model_validator(mode="after")
    def _one_target(self) -> EventModel:
        if (self.controller is None) == (self.unit is None):
            raise ValueError("an event changes either a unit parameter or a controller setpoint")
        if self.unit is not None and self.param is None:
            raise ValueError("an event on a unit needs 'param'")
        return self


class DynamicsModel(BaseModel):
    model_config = ConfigDict(extra="forbid")

    t_end: float | str
    points: int = Field(default=201, ge=2, le=2001)
    units: dict[str, dict[str, Any]] = Field(default_factory=dict)
    controllers: list[ControllerModel] = Field(default_factory=list)
    events: list[EventModel] = Field(default_factory=list)
    rtol: float = Field(default=1e-5, gt=0.0, lt=0.1)
