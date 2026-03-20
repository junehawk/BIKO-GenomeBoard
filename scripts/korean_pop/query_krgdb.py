from typing import Optional
from scripts.common.models import Variant

_KRGDB_CACHE: dict = {}

def _load_krgdb(path: str) -> dict:
    if path in _KRGDB_CACHE:
        return _KRGDB_CACHE[path]
    data = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                key = f"{parts[0]}:{parts[1]}:{parts[2]}>{parts[3]}"
                data[key] = float(parts[4])
    _KRGDB_CACHE[path] = data
    return data

def query_krgdb(variant: Variant, krgdb_path: str = "data/krgdb_freq.tsv") -> Optional[float]:
    data = _load_krgdb(krgdb_path)
    key = f"{variant.chrom}:{variant.pos}:{variant.ref}>{variant.alt}"
    return data.get(key)
