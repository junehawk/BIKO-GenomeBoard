"""Korea4K allele frequency query module.

Korea4K (korea4k.10000genomes.org) is a Korean population whole-genome
sequencing database containing allele frequencies derived from 4,000+
Korean individuals.  This module loads a pre-extracted TSV file into an
in-memory dict with thread-safe singleton caching.

TSV format: chrom  pos  ref  alt  frequency
"""

import threading
from typing import Optional
from scripts.common.models import Variant
from scripts.common.config import get

_KOREA4K_CACHE: dict = {}
_KOREA4K_LOCK = threading.Lock()


def _load_korea4k(path: str) -> dict:
    with _KOREA4K_LOCK:
        if path in _KOREA4K_CACHE:
            return _KOREA4K_CACHE[path]
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

        logging.getLogger(__name__).warning(f"Korea4K file not found: {path}")
    with _KOREA4K_LOCK:
        _KOREA4K_CACHE[path] = data
    return data


def query_korea4k(variant: Variant, korea4k_path: str = None) -> Optional[float]:
    if korea4k_path is None:
        korea4k_path = get("paths.korea4k_freq", "data/korea4k_freq.tsv")
    data = _load_korea4k(korea4k_path)
    key = f"{variant.chrom}:{variant.pos}:{variant.ref}>{variant.alt}"
    return data.get(key)


if __name__ == "__main__":
    import sys
    import json

    if len(sys.argv) < 2:
        print(
            json.dumps(
                {"error": "Usage: python -m scripts.korean_pop.query_korea4k 'chr17:7675088 C>A' [korea4k_path]"}
            )
        )
        sys.exit(1)
    v = Variant.from_string(sys.argv[1])
    path = sys.argv[2] if len(sys.argv) > 2 else "data/korea4k_freq.tsv"
    result = query_korea4k(v, path)
    print(json.dumps({"variant": v.variant_id, "korea4k_freq": result}, indent=2))
