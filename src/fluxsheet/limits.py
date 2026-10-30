"""Resource limits for a shared, public server.

Off by default. With ``FLUXSHEET_DEMO=1`` (the online demo) every calculation request
- waits for its turn: one calculation at a time, so a small server stays responsive,
- has a time budget (``FLUXSHEET_MAX_SECONDS``, default 60 s), checked before every unit
  calculation, so recycles, design specs and sensitivity studies all stop in time,
- and the flowsheet itself must stay within a modest size.
Locally (``fluxsheet serve``, the CLI, Python) nothing is limited.
"""

from __future__ import annotations

import contextlib
import contextvars
import os
import threading
import time
from collections.abc import Iterator

from .errors import SimulationError

DEMO = os.environ.get("FLUXSHEET_DEMO") == "1"
MAX_SECONDS = float(os.environ.get("FLUXSHEET_MAX_SECONDS", "60"))
MAX_BODY_BYTES = 256 * 1024
MAX_UNITS = 60
MAX_COMPONENTS = 20
MAX_STUDY_POINTS = 60
QUEUE_WAIT_SECONDS = 20.0

_LOCAL_HINT = "install Fluxsheet locally (pip install \"fluxsheet[server]\"), where there is no limit"

_deadline: contextvars.ContextVar[float | None] = contextvars.ContextVar("deadline", default=None)
_turn = threading.Semaphore(1)


class Busy(Exception):
    """Another calculation is running and did not finish within the queue wait."""


def check_deadline(where: str | None = None) -> None:
    """Raise if the current request has used up its time budget (no-op outside the demo)."""
    deadline = _deadline.get()
    if deadline is not None and time.monotonic() > deadline:
        raise SimulationError(
            f"the calculation took longer than {MAX_SECONDS:g} s, the limit on this demo server",
            where=where, hint=f"simplify the flowsheet or loosen its specs, or {_LOCAL_HINT}",
        )


@contextlib.contextmanager
def calculation(enabled: bool | None = None) -> Iterator[None]:
    """Run one request's calculation under the demo limits (one at a time, time budget)."""
    if not (DEMO if enabled is None else enabled):
        yield
        return
    if not _turn.acquire(timeout=QUEUE_WAIT_SECONDS):
        raise Busy
    token = _deadline.set(time.monotonic() + MAX_SECONDS)
    try:
        yield
    finally:
        _deadline.reset(token)
        _turn.release()


def check_model_size(model, enabled: bool | None = None) -> None:
    """Reject flowsheets beyond the demo's size limits with a clear message."""
    if not (DEMO if enabled is None else enabled):
        return
    checks = [
        (len(model.units), MAX_UNITS, "units"),
        (len(model.components), MAX_COMPONENTS, "components"),
    ]
    for study in getattr(model, "sensitivity", None) or []:
        points = 1
        for v in study.vary:
            points *= v.points
        checks.append((points, MAX_STUDY_POINTS, f"points in sensitivity study {study.id}"))
    for count, limit, what in checks:
        if count > limit:
            raise SimulationError(
                f"this flowsheet has {count} {what}; the demo server allows at most {limit}",
                hint=_LOCAL_HINT,
            )
