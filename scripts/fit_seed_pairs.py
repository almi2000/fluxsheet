"""Refit the seed set of binary parameters from NIST ThermoML data and print a comparison.

    python scripts/fit_seed_pairs.py            # fit, print a Markdown table
    python scripts/fit_seed_pairs.py --write    # also update src/fluxsheet/data/binary_parameters

Downloads go to the user cache directory (see fluxsheet.fitting.archive), never into the
repository. The first run reads about a thousand articles and takes a while; later runs use
the cache (add --offline to make sure nothing is downloaded).
"""

from __future__ import annotations

import argparse
import sys

from fluxsheet.fitting import pipeline, store
from fluxsheet.fitting.regression import FitOptions

SEED_PAIRS = [
    ("water", "ethanol"),
    ("water", "methanol"),
    ("methanol", "ethanol"),
    ("acetone", "water"),
    ("water", "1-butanol"),
    ("ethanol", "ethyl acetate"),
    ("benzene", "toluene"),
    ("toluene", "methylcyclohexane"),
]


def _fmt(stats: dict, key: str) -> str:
    value = stats.get(key)
    return "-" if value is None else f"{value:.3g}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--models", default="nrtl,uniquac")
    args = parser.parse_args()
    print("| Pair | Model | Points (sets) | Fit AAD T [K] / P [%] / y | ChemSep AAD T / P / y "
          "| UNIFAC est. AAD T / P / y |")
    print("|---|---|---|---|---|---|")
    for a, b in SEED_PAIRS:
        collected = pipeline.collect(a, b, offline=args.offline,
                                     progress=lambda m: print(m, file=sys.stderr))
        for model in args.models.split(","):
            try:
                result = pipeline.fit(collected, FitOptions(model=model))
            except (ValueError, RuntimeError) as exc:
                print(f"| {a}/{b} | {model} | {exc} | | | |")
                continue
            entry = result.entry()
            cells = [" / ".join(_fmt(s, k) for k in ("AAD_T_K", "AAD_P_pct", "AAD_y"))
                     for s in (result.stats, result.comparison.get("ChemSep", {}),
                               result.comparison.get("UNIFAC estimate", {}))]
            print(f"| {a}/{b} | {model.upper()} | {entry['n_points']} ({entry['n_datasets']}) | "
                  + " | ".join(cells) + " |", flush=True)
            if args.write:
                store.save_entry(model, entry)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
