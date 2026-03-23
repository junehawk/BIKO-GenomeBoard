# scripts/pharma/korean_pgx.py
import json
import threading
from pathlib import Path
from typing import Optional
from scripts.common.models import Variant, PgxResult
from scripts.common.config import get

_PGX_DATA = None
_PGX_LOCK = threading.Lock()


def _load_pgx_data() -> list:
    global _PGX_DATA
    with _PGX_LOCK:
        if _PGX_DATA is None:
            path = get("paths.pgx_table") or str(Path(__file__).parent.parent.parent / "data" / "korean_pgx_table.json")
            try:
                with open(path) as f:
                    _PGX_DATA = json.load(f)["genes"]
            except (FileNotFoundError, json.JSONDecodeError, KeyError):
                import logging

                logging.getLogger(__name__).warning(f"PGx data file not found or invalid: {path}")
                _PGX_DATA = []
    return _PGX_DATA


PGX_GENES = set(get("pgx.genes", ["CYP2D6", "CYP2C19", "CYP2C9", "HLA-B", "HLA-A", "NUDT15", "TPMT", "DPYD"]))


def check_korean_pgx(variant: Variant) -> Optional[PgxResult]:
    if variant.gene not in PGX_GENES:
        return None

    pgx_data = _load_pgx_data()
    for entry in pgx_data:
        if entry["gene"] == variant.gene:
            korean_prev = entry.get("korean_freq", 0)
            western_prev = entry.get("western_freq", 0)
            phenotype = ""
            if variant.gene == "CYP2C19":
                phenotype = "Intermediate Metabolizer (*2 carrier)"
            elif variant.gene == "HLA-B":
                phenotype = "HLA-B*5701 carrier — abacavir hypersensitivity risk"
            elif variant.gene == "NUDT15":
                phenotype = "NUDT15 intermediate metabolizer (p.R139C carrier)"
            return PgxResult(
                gene=variant.gene,
                star_allele=entry.get("variant", ""),
                phenotype=phenotype,
                cpic_level=entry["cpic_level"],
                korean_prevalence=korean_prev,
                western_prevalence=western_prev,
                clinical_impact=entry["clinical_impact"],
                cpic_recommendation="",
            )
    return None


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: python -m scripts.pharma.korean_pgx 'chr10:96541616 G>A' [gene]"}))
        sys.exit(1)
    v = Variant.from_string(sys.argv[1])
    if len(sys.argv) > 2:
        v.gene = sys.argv[2]
    result = check_korean_pgx(v)
    if result:
        print(
            json.dumps(
                {
                    "gene": result.gene,
                    "korean_flag": result.korean_flag,
                    "cpic_level": result.cpic_level,
                    "clinical_impact": result.clinical_impact,
                },
                indent=2,
                ensure_ascii=False,
            )
        )
    else:
        print(json.dumps({"result": None}))
