"""Compound lookup on top of the open ``chemicals`` database (~76,000 compounds).

Whether a compound can actually be simulated depends on whether critical data and
correlations exist for it; :func:`search` reports that as ``simulation_ready``.
"""

from __future__ import annotations

import bisect
from dataclasses import asdict, dataclass
from functools import cache, lru_cache

from chemicals.critical import Pc_methods, Tc_methods
from chemicals.identifiers import pubchem_db, search_chemical


@dataclass(frozen=True)
class Compound:
    cas: str
    name: str
    formula: str
    mw: float
    smiles: str
    simulation_ready: bool

    def to_dict(self) -> dict:
        return asdict(self)


class UnknownCompoundError(LookupError):
    pass


@cache
def _name_index() -> tuple[list[str], dict[str, object]]:
    pubchem_db.finish_loading()
    names = sorted(pubchem_db.name_index)
    return names, pubchem_db.name_index


def _is_simulation_ready(cas: str) -> bool:
    return bool(Tc_methods(cas)) and bool(Pc_methods(cas))


def _from_metadata(meta) -> Compound:
    return Compound(
        cas=meta.CASs,
        name=meta.common_name or meta.iupac_name,
        formula=meta.formula,
        mw=meta.MW,
        smiles=meta.smiles,
        simulation_ready=_is_simulation_ready(meta.CASs),
    )


@lru_cache(maxsize=4096)
def lookup(identifier: str) -> Compound:
    """Resolve a name, CAS number, formula or SMILES to a single compound."""
    try:
        meta = search_chemical(identifier.strip())
    except Exception as exc:  # chemicals raises a bare Exception when nothing matches
        raise UnknownCompoundError(
            f"no compound called {identifier!r} in the database; try its CAS number"
        ) from exc
    return _from_metadata(meta)


def search(query: str, limit: int = 20) -> list[Compound]:
    """Prefix search over names and synonyms, best matches first."""
    query = query.strip().lower()
    if not query:
        return []
    names, index = _name_index()
    hits: dict[str, tuple[tuple, object]] = {}

    # Exact CAS / formula hits come first.
    for finder in (pubchem_db.search_CAS, pubchem_db.search_formula):
        try:
            meta = finder(query.upper() if finder is pubchem_db.search_formula else query)
        except ValueError:  # not a CAS number
            meta = None
        if meta:
            hits[meta.CASs] = ((0, 0, 0), meta)

    start = bisect.bisect_left(names, query)
    for name in names[start : start + 2000]:
        if not name.startswith(query):
            break
        meta = index[name]
        is_common = (meta.common_name or "").lower() == name
        rank = (0 if name == query else 1, 0 if is_common else 1, len(name))
        if meta.CASs not in hits or rank < hits[meta.CASs][0]:
            hits[meta.CASs] = (rank, meta)

    ranked = sorted(hits.values(), key=lambda item: item[0])[:limit]
    return [_from_metadata(meta) for _, meta in ranked]
