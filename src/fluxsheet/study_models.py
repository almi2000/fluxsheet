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

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

# Flowsheet totals of the costing (costing.py) that a VariableRef without owner can name.
COSTING_TOTALS = ("tac", "capital_cost", "operating_cost", "bare_module_cost", "purchased_cost")


class VariableRef(BaseModel):
    """A computed value: ``property`` of a stream, the result ``property`` of a unit, or,
    with neither, a flowsheet total of the costing such as ``{"property": "tac"}``."""

    model_config = ConfigDict(extra="forbid")

    stream: str | None = None
    unit: str | None = None
    property: str
    component: str | None = None  # for mole_fraction and component_flow

    @model_validator(mode="after")
    def _one_owner(self) -> VariableRef:
        if self.stream is not None and self.unit is not None:
            raise ValueError("a variable belongs to one of 'stream' or 'unit', not both")
        if self.stream is None and self.unit is None and self.property not in COSTING_TOTALS:
            raise ValueError("a variable belongs to a 'stream' or a 'unit', unless it is a "
                             f"costing total: {', '.join(COSTING_TOTALS)}")
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


class ParamRef(BaseModel):
    """The current input value of a unit parameter, e.g. a feed flow."""

    model_config = ConfigDict(extra="forbid")

    unit: str
    param: str
    component: str | None = None


class Assignment(BaseModel):
    """Set a unit parameter to the value of an expression (SI units)."""

    model_config = ConfigDict(extra="forbid")

    unit: str
    param: str
    component: str | None = None
    expression: str


class CalculatorModel(BaseModel):
    """A calculator block: named variables read from the flowsheet, and parameters set from
    expressions over them, e.g. a make-up flow that follows the purge::

        {"id": "MAKEUP",
         "variables": {"purge": {"stream": "PURGE", "property": "component_flow",
                                 "component": "methanol"}},
         "assign": [{"unit": "FEED2", "param": "flows", "component": "methanol",
                     "expression": "purge"}]}

    See calculator.py for the expression language and how it is solved.
    """

    model_config = ConfigDict(extra="forbid")

    id: str
    variables: dict[str, VariableRef | ParamRef] = Field(default_factory=dict)
    assign: list[Assignment] = Field(min_length=1)
    active: bool = True


class OptimizationVary(VariedParam):
    """A decision variable of an optimization: a unit parameter between bounds."""

    start: float | str | None = None  # default: the parameter's value, else mid-range


class ConstraintModel(BaseModel):
    """``expression`` (SI units) must be ``<=``, ``>=`` or ``=`` to ``value``."""

    model_config = ConfigDict(extra="forbid")

    expression: str
    relation: Literal["<=", ">=", "="] = "<="
    value: float = 0.0


class OptimizationModel(BaseModel):
    """Minimize or maximize an expression over flowsheet values by varying unit parameters,
    like Aspen's Optimization block::

        {"id": "MIN_TAC",
         "variables": {"tac": {"property": "tac"},
                       "x": {"stream": "DIST", "property": "mole_fraction",
                             "component": "benzene"}},
         "objective": "tac", "sense": "minimize",
         "vary": [{"unit": "COLUMN", "param": "reflux_ratio", "lower": 1.2, "upper": 5}],
         "constraints": [{"expression": "x", "relation": ">=", "value": 0.99}]}

    Variables and expressions are those of calculator blocks; ``{"property": "tac"}`` reads
    the costing summary (costing.py). See optimization.py for the method.
    """

    model_config = ConfigDict(extra="forbid")

    id: str
    variables: dict[str, VariableRef | ParamRef] = Field(default_factory=dict)
    objective: str
    sense: Literal["minimize", "maximize"] = "minimize"
    vary: list[OptimizationVary] = Field(min_length=1, max_length=10)
    constraints: list[ConstraintModel] = Field(default_factory=list)
    method: Literal["auto", "slsqp", "cobyla"] = "auto"
    max_iterations: int = Field(default=50, ge=1, le=500)
    tolerance: float = Field(default=1e-5, gt=0.0, lt=1.0)  # on the scaled objective
    step: float = Field(default=1e-3, gt=0.0, le=0.1)  # finite differences, fraction of range


class UnitCostingModel(BaseModel):
    """Per-unit costing choices."""

    model_config = ConfigDict(extra="forbid")

    material: Literal["CS", "SS"] | None = None  # default: the flowsheet's
    U: float | str | None = None  # heat transfer coefficient for exchanger areas
    exclude: bool = False  # leave the unit out (capital and utilities)


class CostingModel(BaseModel):
    """Equipment and utility costing; see costing.py. Prices in US$ (``$/GJ``, ``$/kWh``)."""

    model_config = ConfigDict(extra="forbid")

    cepci: float = Field(default=800.0, gt=0.0)  # cost index of the estimate's year
    hours_per_year: float = Field(default=8000.0, gt=0.0, le=8760.0)
    interest_rate: float = Field(default=0.10, ge=0.0, le=1.0)
    plant_life: float = Field(default=10.0, gt=0.0, le=100.0)  # years
    annualization_factor: float | None = Field(default=None, gt=0.0, le=10.0)  # 1/yr
    electricity_price: float = Field(default=0.07, ge=0.0)  # $/kWh
    heating_price: float | None = Field(default=None, ge=0.0)  # $/GJ; default by temperature
    cooling_price: float | None = Field(default=None, ge=0.0)  # $/GJ; default by temperature
    utility_prices: dict[str, float] = Field(default_factory=dict)  # $/GJ per utility id
    material: Literal["CS", "SS"] = "CS"
    U: float | str = 500.0  # W/(m2 K), exchangers without a given area
    tray_efficiency: float = Field(default=0.7, gt=0.0, le=1.0)
    units: dict[str, UnitCostingModel] = Field(default_factory=dict)
