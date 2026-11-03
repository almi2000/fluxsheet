"""Dynamic simulation: hold-ups, controllers and disturbances over time.

Starts from the steady state and integrates the hold-up balances of the units that have a
dynamic model (see units.py) with a stiff integrator; see simulator.py.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..flowsheet import FlowsheetModel


def simulate(model: FlowsheetModel | dict) -> dict:
    """Run the flowsheet's dynamic simulation; see simulator.simulate."""
    from .simulator import simulate as run

    return run(model)
