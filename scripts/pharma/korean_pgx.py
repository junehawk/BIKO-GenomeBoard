# scripts/pharma/korean_pgx.py
import json
from pathlib import Path
from typing import Optional
from scripts.common.models import Variant, PgxResult

_PGX_DATA = None

def _load_pgx_data() -> list:
    global _PGX_DATA
    if _PGX_DATA is None:
        path = Path(__file__).parent.parent.parent / "data" / "korean_pgx_table.json"
        with open(path) as f:
            _PGX_DATA = json.load(f)["genes"]
    return _PGX_DATA

PGX_GENES = {"CYP2D6", "CYP2C19", "CYP2C9", "HLA-B", "NUDT15"}

def check_korean_pgx(variant: Variant) -> Optional[PgxResult]:
    if variant.gene not in PGX_GENES:
        return None

    pgx_data = _load_pgx_data()
    for entry in pgx_data:
        if entry["gene"] == variant.gene:
            korean_prev = entry.get("korean_freq", 0)
            western_prev = entry.get("western_freq", 0)
            return PgxResult(
                gene=variant.gene,
                star_allele=entry.get("variant", ""),
                phenotype="",  # to be filled by LLM agent
                cpic_level=entry["cpic_level"],
                korean_prevalence=korean_prev,
                western_prevalence=western_prev,
                clinical_impact=entry["clinical_impact"],
                cpic_recommendation="",
            )
    return None
