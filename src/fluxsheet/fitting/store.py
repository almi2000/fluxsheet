"""Versioned files of fitted binary parameters, shipped in ``fluxsheet/data/binary_parameters``.

One JSON file per model (``nrtl.json``, ``uniquac.json``). Each pair records the parameters
in the component order it lists, the temperature and pressure range of the data, the
number of points, the fit statistics, a comparison with other parameter sets on the same
points, and the DOI of every source article. The files are written by
``fluxsheet fit-binary --write`` and meant to be reviewed like code.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import cache
from pathlib import Path

FORMAT_VERSION = 1
DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "binary_parameters"

_HEADERS = {
    "nrtl": {
        "equation": "tau_ij = a_ij + b_ij / T (T in K), G_ij = exp(-alpha tau_ij); "
                    "the thermo.NRTL convention (tau_as, tau_bs, alpha_cs)",
    },
    "uniquac": {
        "equation": "tau_ij = exp(a_ij + b_ij / T) (T in K), r and q from UNIFAC groups; "
                    "the thermo.UNIQUAC convention (tau_as, tau_bs)",
    },
}
_COMMON = {
    "description": "Binary parameters regressed to vapor-liquid equilibrium data from the "
                   "NIST TRC ThermoML Archive by `fluxsheet fit-binary`. Phase model: "
                   "modified Raoult's law with ideal vapor and the thermo vapor-pressure "
                   "correlations, as in Fluxsheet's activity-coefficient property packages.",
    "data_source": "NIST TRC ThermoML Archive, https://trc.nist.gov/ThermoML/ "
                   "(doi:10.18434/mds2-2422). Only derived parameters and references are "
                   "stored here; see docs/data-sources.md for the terms.",
}


@dataclass(frozen=True)
class FittedPair:
    """A stored parameter set, oriented to the order ``(cas_i, cas_j)`` it was asked for."""

    a_ij: float
    a_ji: float
    b_ij: float
    b_ji: float
    alpha: float | None
    record: dict  # the stored entry, for provenance

    @property
    def dois(self) -> list[str]:
        return [s["doi"] for s in self.record.get("sources", [])]

    def provenance(self) -> dict:
        r = self.record
        return {"references": self.dois, "T_range_K": r.get("T_range_K"),
                "n_points": r.get("n_points"), "fit": r.get("fit")}


def path_for(model: str, directory: Path | None = None) -> Path:
    return (directory or DATA_DIR) / f"{model}.json"


def load(model: str, directory: Path | None = None) -> dict:
    path = path_for(model, directory)
    if not path.exists():
        return {"format_version": FORMAT_VERSION, "model": model, **_HEADERS.get(model, {}),
                **_COMMON, "pairs": []}
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("format_version") != FORMAT_VERSION:
        raise ValueError(f"{path}: unsupported format_version {data.get('format_version')}")
    return data


@cache
def _index(model: str, directory: Path | None) -> dict[tuple[str, str], dict]:
    return {tuple(c["cas"] for c in e["components"]): e
            for e in load(model, directory).get("pairs", [])}


def lookup(model: str, cas_i: str, cas_j: str, directory: Path | None = None) -> FittedPair | None:
    """The fitted parameters for a pair in either order, or None."""
    index = _index(model, directory)
    if (e := index.get((cas_i, cas_j))) is not None:
        return FittedPair(e["a12"], e["a21"], e["b12"], e["b21"], e.get("alpha"), e)
    if (e := index.get((cas_j, cas_i))) is not None:
        return FittedPair(e["a21"], e["a12"], e["b21"], e["b12"], e.get("alpha"), e)
    return None


def save_entry(model: str, entry: dict, directory: Path | None = None) -> Path:
    """Insert or replace a pair's entry (either component order), keeping the file sorted."""
    data = load(model, directory)
    data.update({"format_version": FORMAT_VERSION, "model": model, **_HEADERS.get(model, {}),
                 **_COMMON})
    key = {c["cas"] for c in entry["components"]}
    pairs = [e for e in data.get("pairs", []) if {c["cas"] for c in e["components"]} != key]
    pairs.append(entry)
    pairs.sort(key=lambda e: tuple(c["name"] for c in e["components"]))
    data["pairs"] = pairs
    path = path_for(model, directory)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8",
                    newline="\n")
    _index.cache_clear()
    return path
