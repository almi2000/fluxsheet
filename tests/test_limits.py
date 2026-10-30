"""Demo-server limits (FLUXSHEET_DEMO=1): time budget, one calculation at a time, sizes."""

import threading
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from fluxsheet import FlowsheetModel, limits, solve
from fluxsheet.errors import SimulationError
from fluxsheet.server import app

EXAMPLES = Path(__file__).parent.parent / "examples"


def test_limits_are_off_by_default():
    assert not limits.DEMO
    with limits.calculation():
        limits.check_deadline()  # no budget set, never raises
    r = solve(FlowsheetModel.load(EXAMPLES / "gas_recycle.json"))
    assert r.ok


def test_time_budget_stops_a_long_calculation(monkeypatch):
    monkeypatch.setattr(limits, "MAX_SECONDS", 0.0)
    with limits.calculation(enabled=True):
        r = solve(FlowsheetModel.load(EXAMPLES / "gas_recycle.json"))
    assert not r.ok
    assert "longer than 0 s" in r.error["message"]
    assert "locally" in r.error["hint"]


def test_one_calculation_at_a_time(monkeypatch):
    monkeypatch.setattr(limits, "QUEUE_WAIT_SECONDS", 0.05)
    inside, release = threading.Event(), threading.Event()

    def hold():
        with limits.calculation(enabled=True):
            inside.set()
            release.wait(5)

    t = threading.Thread(target=hold)
    t.start()
    inside.wait(5)
    try:
        with pytest.raises(limits.Busy), limits.calculation(enabled=True):
            pass
    finally:
        release.set()
        t.join()
    with limits.calculation(enabled=True):  # free again
        pass


def test_model_size_limits():
    model = FlowsheetModel.load(EXAMPLES / "methanol_plant.json")
    limits.check_model_size(model, enabled=True)  # the examples fit
    model.components = [f"c{i}" for i in range(limits.MAX_COMPONENTS + 1)]
    with pytest.raises(SimulationError, match="at most"):
        limits.check_model_size(model, enabled=True)


def test_demo_server_rejects_large_bodies_and_reports_limits(monkeypatch):
    monkeypatch.setattr(limits, "DEMO", True)
    client = TestClient(app)
    big = "x" * (limits.MAX_BODY_BYTES + 1)
    r = client.post("/api/solve", content=big, headers={"Content-Type": "application/json"})
    assert r.status_code == 413

    model = FlowsheetModel.load(EXAMPLES / "ethanol_flash.json").model_dump(by_alias=True)
    model["components"] = [f"c{i}" for i in range(limits.MAX_COMPONENTS + 1)]
    body = client.post("/api/solve", json=model).json()
    assert not body["ok"] and "demo server allows" in body["error"]["message"]
