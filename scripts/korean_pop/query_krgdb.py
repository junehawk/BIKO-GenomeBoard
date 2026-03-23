import threading
from typing import Optional
from scripts.common.models import Variant
from scripts.common.config import get

_KRGDB_CACHE: dict = {}
_KRGDB_LOCK = threading.Lock()


def _load_krgdb(path: str) -> dict:
    with _KRGDB_LOCK:
        if path in _KRGDB_CACHE:
            return _KRGDB_CACHE[path]
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

        logging.getLogger(__name__).warning(f"KRGDB file not found: {path}")
    with _KRGDB_LOCK:
        _KRGDB_CACHE[path] = data
    return data


def query_krgdb(variant: Variant, krgdb_path: str = None) -> Optional[float]:
    if krgdb_path is None:
        krgdb_path = get("paths.krgdb", "data/krgdb_freq.tsv")
    data = _load_krgdb(krgdb_path)
    key = f"{variant.chrom}:{variant.pos}:{variant.ref}>{variant.alt}"
    return data.get(key)


if __name__ == "__main__":
    import sys
    import json

    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: python -m scripts.korean_pop.query_krgdb 'chr17:7577120 G>A' [krgdb_path]"}))
        sys.exit(1)
    v = Variant.from_string(sys.argv[1])
    path = sys.argv[2] if len(sys.argv) > 2 else "data/krgdb_freq.tsv"
    result = query_krgdb(v, path)
    print(json.dumps({"variant": v.variant_id, "krgdb_freq": result}, indent=2))
