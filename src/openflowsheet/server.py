"""HTTP API used by the web editor. Start with ``openflowsheet serve``."""

from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import __version__
from .compounds import UnknownCompoundError, lookup, search
from .errors import SimulationError
from .flowsheet import FlowsheetModel
from .solver import solve
from .thermo_package import METHODS, suggest_method
from .units import UNIT_TYPES

app = FastAPI(title="OpenFlowsheet", version=__version__)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_methods=["*"],
    allow_headers=["*"],
)


class ComponentList(BaseModel):
    components: list[str]


@app.get("/api/health")
def health() -> dict:
    return {"status": "ok", "version": __version__}


@app.get("/api/compounds")
def compounds(q: str, limit: int = 15) -> list[dict]:
    return [c.to_dict() for c in search(q, limit=min(limit, 50))]


@app.get("/api/compounds/{identifier}")
def compound(identifier: str) -> dict:
    try:
        return lookup(identifier).to_dict()
    except UnknownCompoundError as exc:
        raise HTTPException(404, str(exc)) from exc


@app.get("/api/methods")
def methods() -> list[dict]:
    return [m.__dict__ for m in METHODS.values()]


@app.post("/api/methods/suggest")
def suggest(body: ComponentList) -> dict:
    if not body.components:
        raise HTTPException(400, "add at least one component")
    try:
        return suggest_method(body.components).__dict__
    except (UnknownCompoundError, SimulationError) as exc:
        raise HTTPException(400, str(exc)) from exc


@app.get("/api/unit-types")
def unit_types() -> list[dict]:
    return [cls.schema() for cls in UNIT_TYPES.values()]


_examples = Path(__file__).resolve().parents[2] / "examples"


@app.get("/api/examples")
def examples() -> list[dict]:
    out = []
    for path in sorted(_examples.glob("*.json")):
        model = FlowsheetModel.load(path)
        out.append({"id": path.stem, "name": model.name, "description": model.description})
    return out


@app.get("/api/examples/{name}")
def example(name: str) -> dict:
    path = _examples / f"{name}.json"
    if not name.replace("_", "").isalnum() or not path.is_file():
        raise HTTPException(404, f"no example called {name!r}")
    return FlowsheetModel.load(path).model_dump(by_alias=True, exclude_none=True)


@app.post("/api/solve")
def solve_flowsheet(model: FlowsheetModel) -> dict:
    try:
        return solve(model).to_dict()
    except UnknownCompoundError as exc:
        return {"ok": False, "error": {"message": str(exc), "where": None,
                                       "hint": "search for the compound and pick it from the list"}}


# Serve the built web editor (web/dist) from the same origin when it exists.
_dist = Path(__file__).resolve().parents[2] / "web" / "dist"
if _dist.is_dir():
    app.mount("/", StaticFiles(directory=_dist, html=True), name="web")
