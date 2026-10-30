"""From two compound names to a fitted, documented parameter set (``fluxsheet fit-binary``)."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from thermo.interaction_parameters import IPDB

from ..compounds import UnknownCompoundError, lookup
from . import archive
from .models import Pair
from .regression import FitOptions, FitResult, compare, fit_binary
from .thermoml import ThermoMLCompound, ThermoMLError, VLEDataset, parse_thermoml

_CHEMSEP = {"nrtl": "ChemSep NRTL", "uniquac": "ChemSep UNIQUAC"}


@dataclass
class Collected:
    pair: Pair
    datasets: list[tuple[VLEDataset, str]]  # oriented to the pair, with citations
    files: int  # ThermoML files read
    skipped: list[str]  # files that could not be downloaded or parsed


@lru_cache(maxsize=1024)
def _cas_of(compound: ThermoMLCompound) -> str | None:
    """CAS number of a ThermoML compound: by InChIKey, then InChI, then its names."""
    candidates = []
    if compound.inchikey:
        candidates.append(f"InChIKey={compound.inchikey}")
    if compound.inchi:
        candidates.append(compound.inchi)
    candidates.extend(compound.names)
    for identifier in candidates:
        try:
            found = lookup(identifier)
        except (UnknownCompoundError, ValueError):
            continue
        if compound.formula and found.formula and _formula(found.formula) != _formula(
                compound.formula):
            continue  # a name that resolves to a different compound
        return found.cas
    return None


def _formula(formula: str) -> str:
    return formula.replace(" ", "")


def _inchikey(cas: str) -> str:
    from chemicals.identifiers import search_chemical

    return search_chemical(cas).InChI_key


def datasets_for_pair(pair: Pair, paths, skipped: list[str]) -> list[tuple[VLEDataset, str]]:
    """Binary VLE data sets of exactly this pair in the given ThermoML files, oriented."""
    out = []
    for path in paths:
        try:
            report = parse_thermoml(Path(path))
        except (ThermoMLError, OSError) as exc:
            skipped.append(f"{path}: {exc}")
            continue
        for ds in report.datasets:
            cas = tuple(_cas_of(c) for c in ds.compounds)
            if cas == pair.cas:
                out.append((ds, report.citation))
            elif cas == pair.cas[::-1]:
                out.append((ds.swapped(), report.citation))
    return out


def collect(component_a: str, component_b: str, cache: Path | None = None,
            max_files: int = 0, offline: bool = False, progress=None) -> Collected:
    a, b = lookup(component_a), lookup(component_b)
    pair = Pair.from_cas(a.cas, b.cas)
    cache = cache or archive.cache_dir()
    hits = archive.search((_inchikey(a.cas), _inchikey(b.cas)), cache, offline=offline)
    if max_files:
        hits = hits[:max_files]
    paths, skipped = [], []
    for i, (doi, _) in enumerate(hits):
        if progress:
            progress(f"[{i + 1}/{len(hits)}] {doi}")
        try:
            paths.append(archive.fetch(doi, cache, offline=offline))
        except (archive.ArchiveError, OSError) as exc:
            skipped.append(f"{doi}: {exc}")
    return Collected(pair, datasets_for_pair(pair, paths, skipped), len(paths), skipped)


def chemsep_parameters(model: str, pair: Pair):
    """ChemSep's (a12, a21, b12, b21) and alpha for the pair as Fluxsheet uses them, or None."""
    table = _CHEMSEP[model]
    if not IPDB.has_ip_specific(table, list(pair.cas), "bij"):
        return None
    b = IPDB.get_ip_asymmetric_matrix(table, list(pair.cas), "bij")
    alpha = (IPDB.get_ip_asymmetric_matrix(table, list(pair.cas), "alphaij")[0][1]
             if model == "nrtl" else None)
    return (0.0, 0.0, b[0][1], b[1][0]), alpha


def fit(collected: Collected, options: FitOptions) -> FitResult:
    """Fit and add the ChemSep set and the UNIFAC estimate as comparisons, on the same points."""
    from ..thermo_package import _estimate_binary

    chemsep = chemsep_parameters(options.model, collected.pair)
    result = fit_binary(collected.pair, collected.datasets, options,
                        extra_starts=[chemsep[0]] if chemsep else [])
    if chemsep:
        compare(result, "ChemSep", chemsep[0], chemsep[1])
    est = _estimate_binary(options.model, *collected.pair.cas)
    if est:
        compare(result, "UNIFAC estimate", est[:4], est[4])
    return result
