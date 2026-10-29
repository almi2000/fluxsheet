"""Fluxsheet: an open, friendly chemical process simulator."""

from .errors import SimulationError
from .flowsheet import FlowsheetModel
from .solver import SolveResult, solve

__version__ = "0.1.0"
__all__ = ["FlowsheetModel", "SimulationError", "SolveResult", "__version__", "solve"]
