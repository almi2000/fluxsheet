"""Enforce conservation diagnostics for every successful stationary solve."""

import copy
import math

import pytest

from fluxsheet.solver import Solver


@pytest.fixture
def expected_balance_exceptions():
    """Map exact diagnostic keys to a reason for deliberately exceptional cases."""
    return {}


@pytest.fixture(autouse=True)
def enforce_balance_diagnostics(monkeypatch, expected_balance_exceptions):
    original = Solver.solve
    snapshots = []

    def recording_solve(self):
        result = original(self)
        if result.ok:
            snapshots.append(copy.deepcopy(result.balance_checks))
        return result

    monkeypatch.setattr(Solver, "solve", recording_solve)
    yield
    observed = set()
    errors = []
    for checks in snapshots:
        if not checks:
            errors.append("successful solve has no balance diagnostics")
            continue
        for quantity in ("mass", "energy"):
            if not any(c["scope"] == "flowsheet" and c["quantity"] == quantity for c in checks):
                errors.append(f"missing flowsheet {quantity} diagnostic")
        for c in checks:
            status = c["status"]
            key = (c["scope"], c["quantity"], c["component"], status)
            if status == "ok":
                if not all(isinstance(c[field], (float, int)) and math.isfinite(c[field])
                           for field in ("residual", "tolerance")) or c["tolerance"] < 0:
                    errors.append(f"invalid numerical diagnostic: {c}")
            elif status in ("warning", "not_checkable"):
                observed.add(key)
                if not expected_balance_exceptions.get(key):
                    errors.append(f"unexpected balance diagnostic: {c}")
            else:
                errors.append(f"unknown balance status: {c}")
    missing = set(expected_balance_exceptions) - observed
    if missing:
        errors.append(f"declared balance exceptions not observed: {sorted(missing, key=str)}")
    assert not errors, "\n".join(errors)
