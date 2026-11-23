"""Download the DWSIM sample simulations registered in tests/data/dwsim/real-samples.json.

    python scripts/fetch_dwsim_samples.py --cache-dir PATH

The files come from DWSIM's own repository at a fixed commit and are checked against their
SHA-256. They are GPL-licensed and stay outside this repository; point the tests at them with
FLUXSHEET_DWSIM_SAMPLES=PATH. Exit status 0 when every file is present and verified, 2
otherwise.
"""

import argparse
import hashlib
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "tests" / "data" / "dwsim" / "real-samples.json"
MAX_BYTES = 5_000_000


def url(manifest: dict, name: str) -> str:
    repo = manifest["source"].removeprefix("https://github.com/")
    path = urllib.parse.quote(f"{manifest['directory']}/{name}")
    return f"https://raw.githubusercontent.com/{repo}/{manifest['commit']}/{path}"


def verified(path: Path, sample: dict) -> bool:
    return path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == sample["sha256"]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cache-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    cache = args.cache_dir.resolve()
    if cache.is_relative_to(ROOT):
        parser.exit(2, "the sample cache must be outside the repository (GPL files)\n")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    cache.mkdir(parents=True, exist_ok=True)
    failed = 0
    for sample in manifest["samples"]:
        target = cache / sample["file"]
        if not verified(target, sample):
            with urllib.request.urlopen(url(manifest, sample["file"]), timeout=60) as r:
                data = r.read(MAX_BYTES + 1)
            if len(data) > MAX_BYTES or hashlib.sha256(data).hexdigest() != sample["sha256"]:
                print(f"{sample['file']}: size or SHA-256 mismatch, not saved", file=sys.stderr)
                failed += 1
                continue
            target.write_bytes(data)
        print(f"{sample['file']}: ok")
    return 2 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
