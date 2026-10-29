"""Find and download ThermoML files from the NIST TRC ThermoML Archive.

Only the archive's published endpoints are used:

* search: ``https://trc.nist.gov/ThermoML-API/objects?query=...`` (the API behind
  https://trc.nist.gov/ThermoML/Search), queried by the compounds' standard InChIKeys;
* files: ``https://trc.nist.gov/ThermoML/<DOI>.xml``, one per journal article.

Downloads are cached outside the repository, by default in the user cache directory
(``%LOCALAPPDATA%\\fluxsheet\\thermoml`` on Windows, ``~/.cache/fluxsheet/thermoml``
elsewhere; ``FLUXSHEET_THERMOML_CACHE`` overrides it), so the raw data are never
committed. Files are treated as untrusted: DOIs are validated before they become paths,
responses are size-limited, and parsing happens in :mod:`.thermoml` with ``defusedxml``.
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

SEARCH_URL = "https://trc.nist.gov/ThermoML-API/objects"
FILE_URL = "https://trc.nist.gov/ThermoML/{doi}.xml"
MAX_DOWNLOAD = 30_000_000  # bytes
_DOI = re.compile(r"^10\.\d{4,9}/[A-Za-z0-9._()\-]+$")
_USER_AGENT = "fluxsheet-fit-binary (+https://github.com/almi2000/fluxsheet)"


class ArchiveError(RuntimeError):
    pass


def cache_dir() -> Path:
    env = os.environ.get("FLUXSHEET_THERMOML_CACHE")
    if env:
        return Path(env)
    if sys.platform == "win32" and os.environ.get("LOCALAPPDATA"):
        return Path(os.environ["LOCALAPPDATA"]) / "fluxsheet" / "thermoml"
    base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
    return Path(base) / "fluxsheet" / "thermoml"


def _check_doi(doi: str) -> str:
    if not _DOI.match(doi) or ".." in doi:
        raise ArchiveError(f"refusing unexpected DOI {doi!r}")
    return doi


def file_path(doi: str, cache: Path) -> Path:
    prefix, suffix = _check_doi(doi).split("/", 1)
    return cache / prefix / f"{suffix}.xml"


def _get(url: str, timeout: float = 120.0) -> bytes:
    if not url.startswith(("https://trc.nist.gov/ThermoML/", SEARCH_URL)):
        raise ArchiveError(f"refusing to download from {url}")
    request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        data = response.read(MAX_DOWNLOAD + 1)
    if len(data) > MAX_DOWNLOAD:
        raise ArchiveError(f"{url} is larger than {MAX_DOWNLOAD} bytes")
    return data


def search(inchikeys: tuple[str, str], cache: Path, offline: bool = False,
           page_size: int = 200) -> list[tuple[str, int]]:
    """Articles that report binary phase-equilibrium compositions and mention both compounds.

    Returns ``(doi, binary composition points)`` pairs, most data first. The archive's
    summary does not say *which* binary was measured, so files still need parsing.
    """
    key = "_".join(sorted(inchikeys))
    path = cache / "search" / f"{key}.json"
    if path.exists():
        return [tuple(item) for item in json.loads(path.read_text(encoding="utf-8"))]
    if offline:
        raise ArchiveError(f"no cached search for {inchikeys} and offline mode is on")
    query = "type:TRCTml4 AND (" + " AND ".join(f'"{k}"' for k in inchikeys) + ")"
    hits: dict[str, int] = {}
    page = 0
    while True:
        params = urllib.parse.urlencode({
            "query": query, "pageNum": page, "pageSize": page_size,
            "requestContext": json.dumps({"Citation": True, "data_summary": True}),
        })
        data = json.loads(_get(f"{SEARCH_URL}?{params}"))
        results = data.get("results") or []
        for r in results:
            content = r.get("content") or {}
            doi = (content.get("Citation") or {}).get("sDOI")
            binary = (content.get("data_summary") or {}).get("binary") or {}
            comp = binary.get("CompositionAtPhaseEquilibrium")
            if doi and comp and _DOI.match(doi):
                hits[doi] = int(comp.get("data_points", 0))
        page += 1
        if not results or page * page_size >= int(data.get("size", 0)):
            break
    ranked = sorted(hits.items(), key=lambda item: (-item[1], item[0]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(ranked), encoding="utf-8")
    return ranked


def fetch(doi: str, cache: Path, offline: bool = False, pause: float = 0.3) -> Path:
    """The cached ThermoML file of an article, downloading it first if needed."""
    path = file_path(doi, cache)
    if path.exists():
        return path
    if offline:
        raise ArchiveError(f"{doi} is not cached and offline mode is on")
    data = _get(FILE_URL.format(doi=doi))
    if not data.lstrip().startswith(b"<?xml"):
        raise ArchiveError(f"{doi}: the archive did not return an XML document")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".part")
    tmp.write_bytes(data)
    tmp.replace(path)
    time.sleep(pause)  # be gentle with a public service
    return path
