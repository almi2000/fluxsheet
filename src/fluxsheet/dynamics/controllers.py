"""PID controllers in percent of range, as on a control system faceplate.

The measurement passes a sensor lag (two first-order stages of half the lag each); the
filtered value ``m`` is a state, so the controller output depends on states only and never
forms an algebraic loop with the process, and so is its rate of change for the derivative.
The integral of the error is the third state. With the output at a limit, the integral is
pulled back by the amount the output is clipped (back-calculation anti-windup, tracking time
Ti), so it cannot wind up; unlike switching the integration off, this keeps the equations
smooth, which the stiff integrator's Jacobian needs. The derivative acts on the measurement, not
the error, so setpoint steps do not kick the output.

    e   = ±(SP - m) / span(PV) * 100            (+ for reverse, - for direct action)
    OP% = bias + K (e + I / Ti ∓ Td d(m%)/dt)   clipped to 0 … 100
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .models import ControllerModel


@dataclass
class Controller:
    model: ControllerModel
    pv_lo: float
    pv_hi: float
    op_lo: float
    op_hi: float
    setpoint: float
    bias: float  # % of the output range at the start
    Ti: float | None
    Td: float | None
    lag: float

    n_states = 3  # measurement after the first and second lag stage, integral of the error

    def initial(self, measured: float) -> np.ndarray:
        return np.array([measured, measured, 0.0])

    def scale(self) -> np.ndarray:
        span = self.pv_hi - self.pv_lo
        size = max(abs(self.setpoint), span)
        return np.array([size, size, 100.0 * max(self.Ti or 1.0, 1.0)])

    def _error(self, m: float) -> float:
        e = (self.setpoint - m) / (self.pv_hi - self.pv_lo) * 100.0
        return e if self.model.action == "reverse" else -e

    def output(self, x: np.ndarray) -> tuple[float, float]:
        """Output in engineering units (SI), and in % of its range (before clipping)."""
        m1, m, integral = x
        op = self.bias + self.model.gain * self._error(m)
        if self.Ti:
            op += self.model.gain * integral / self.Ti
        if self.Td:
            dm = (m1 - m) / (0.5 * self.lag) / (self.pv_hi - self.pv_lo) * 100.0
            op -= self.model.gain * self.Td * (dm if self.model.action == "reverse" else -dm)
        clipped = min(max(op, 0.0), 100.0)
        return self.op_lo + clipped / 100.0 * (self.op_hi - self.op_lo), op

    def derivatives(self, x: np.ndarray, measured: float) -> np.ndarray:
        m1, m, _ = x
        e = self._error(m)
        _, op = self.output(x)
        clipped = min(max(op, 0.0), 100.0) - op  # % the limits take off the output
        half = 0.5 * self.lag
        integral = e + clipped / self.model.gain if self.Ti and self.model.gain else 0.0
        return np.array([(measured - m1) / half, (m1 - m) / half, integral])
