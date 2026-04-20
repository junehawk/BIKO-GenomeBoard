"""NARD2 allele frequency query module.

NARD2 (nard.macrogen.com) is the National Archiving of Reference Database
version 2, providing Korean population allele frequencies from large-scale
whole-genome sequencing by Macrogen.  This module loads a pre-extracted TSV
file into an in-memory dict with thread-safe singleton caching.

TSV format: chrom  pos  ref  alt  frequency
"""

import threading
from typing import Optional

from scripts.common.config import get
from scripts.common.models import Variant

_NARD2_CACHE: dict = {}
_NARD2_LOCK = threading.Lock()


def _load_nard2(path: str) -> dict:
    with _NARD2_LOCK:
        if path in _NARD2_CACHE:
            return _NARD2_CACHE[path]
    data = {}
    try:
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    key = f"{parts[0]}:{parts[1]}:{parts[2]}>{parts[3]}"
                    data[key] = float(parts[4])
    except FileNotFoundError:
        import logging

        logging.getLogger(__name__).warning(f"NARD2 file not found: {path}")
    with _NARD2_LOCK:
        _NARD2_CACHE[path] = data
    return data


def query_nard2(variant: Variant, nard2_path: str = None) -> Optional[float]:
    if nard2_path is None:
        nard2_path = get("paths.nard2_freq", "data/nard2_freq.tsv")
    data = _load_nard2(nard2_path)
    key = f"{variant.chrom}:{variant.pos}:{variant.ref}>{variant.alt}"
    return data.get(key)


if __name__ == "__main__":
    import json
    import sys

    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: python -m scripts.population.query_nard2 'chr17:7675088 C>A' [nard2_path]"}))
        sys.exit(1)
    v = Variant.from_string(sys.argv[1])
    path = sys.argv[2] if len(sys.argv) > 2 else "data/nard2_freq.tsv"
    result = query_nard2(v, path)
    print(json.dumps({"variant": v.variant_id, "nard2_freq": result}, indent=2))
