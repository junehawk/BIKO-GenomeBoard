# scripts/pharma/korean_pgx.py
"""Korean-population pharmacogenomics (PGx) analysis.

Provides two levels of API:

* ``check_korean_pgx(variant)`` -- per-variant SNV matching against the
  local ``korean_pgx_table.json`` (builtin engine, always available).
* ``get_pgx_results(variants, germline_vcf, config)`` -- unified entry
  point that prefers PharmCAT when a germline VCF is supplied and
  PharmCAT is installed, falling back to the builtin engine otherwise.
"""

import json
import logging
import threading
from pathlib import Path
from typing import List, Optional

from scripts.common.config import get
from scripts.common.models import PgxResult, Variant

logger = logging.getLogger(__name__)

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
                logger.warning("PGx data file not found or invalid: %s", path)
                _PGX_DATA = []
    return _PGX_DATA


PGX_GENES = set(get("pgx.genes", ["CYP2D6", "CYP2C19", "CYP2C9", "HLA-B", "HLA-A", "NUDT15", "TPMT", "DPYD"]))


def check_korean_pgx(variant: Variant) -> Optional[PgxResult]:
    """Return a :class:`PgxResult` if *variant* matches a known PGx gene."""
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
            elif variant.gene == "CYP3A5":
                phenotype = "CYP3A5 expressor (*1 carrier) — higher tacrolimus dose needed"
            elif variant.gene == "UGT1A1":
                phenotype = "UGT1A1 poor metabolizer (*6 carrier) — irinotecan toxicity risk"
            elif variant.gene == "SLCO1B1":
                phenotype = "SLCO1B1 decreased function (*15 carrier) — statin myopathy risk"
            elif variant.gene == "VKORC1":
                phenotype = "VKORC1 low-dose warfarin phenotype (-1639A carrier)"
            elif variant.gene == "G6PD":
                phenotype = "G6PD deficient — rasburicase contraindicated"
            elif variant.gene == "IFNL3":
                phenotype = "IFNL3 favorable responder (CC genotype)"
            elif variant.gene == "CYP1A2":
                phenotype = "CYP1A2 poor metabolizer (*1C carrier)"
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


# ---------------------------------------------------------------------------
# PharmCAT → PgxResult conversion
# ---------------------------------------------------------------------------


def _convert_pharmcat_to_pgx_hits(pharmcat_result) -> list[dict]:
    """Convert a :class:`PharmCATResult` into the legacy pgx_hits list format.

    Each entry mirrors the dict produced by ``PgxResult`` serialisation in
    ``orchestrate.py`` so the report templates work unchanged.
    """
    hits: list[dict] = []
    for gene, diplotype in pharmcat_result.diplotypes.items():
        phenotype = pharmcat_result.phenotypes.get(gene, "")
        # Collect drug recommendations for this gene
        gene_drugs = [r for r in pharmcat_result.drug_recommendations if r.get("gene") == gene]
        cpic_rec = ""
        cpic_level = ""
        classification = ""
        if gene_drugs:
            cpic_rec = gene_drugs[0].get("guideline", "")
            classification = gene_drugs[0].get("classification", "")
            cpic_level = "A" if classification == "Actionable PGx" else "B"

        hits.append(
            {
                "gene": gene,
                "star_allele": diplotype,
                "phenotype": phenotype,
                "cpic_level": cpic_level,
                "korean_prevalence": 0.0,
                "western_prevalence": 0.0,
                "clinical_impact": classification,
                "cpic_recommendation": cpic_rec,
                "korean_flag": False,
            }
        )
    return hits


# ---------------------------------------------------------------------------
# Unified PGx entry point
# ---------------------------------------------------------------------------


def get_pgx_results(
    variants: List[Variant],
    germline_vcf: str | None = None,
    config: dict | None = None,
) -> dict:
    """Unified PGx entry point.

    Returns a dict with keys compatible with ``report_data``:

    * ``pgx_hits`` -- list of dicts (legacy format for templates)
    * ``pgx_source`` -- ``"pharmcat"`` | ``"builtin"`` | ``"builtin_limited"``
    * ``pharmcat_version`` -- version string (empty when builtin)
    * ``germline_provided`` -- whether a germline VCF was supplied
    * ``warnings`` -- informational warnings
    """
    warnings: list[str] = []

    # --- PharmCAT path (when germline VCF provided) -----------------------
    if germline_vcf:
        try:
            from scripts.pharma.pharmcat_runner import is_pharmcat_available, run_pharmcat

            if is_pharmcat_available(config):
                pharmcat_result = run_pharmcat(germline_vcf, config=config)
                if pharmcat_result is not None:
                    hits = _convert_pharmcat_to_pgx_hits(pharmcat_result)
                    return {
                        "pgx_hits": hits,
                        "pgx_source": "pharmcat",
                        "pharmcat_version": pharmcat_result.version,
                        "germline_provided": True,
                        "warnings": pharmcat_result.warnings,
                    }
                else:
                    warnings.append("PharmCAT execution failed; falling back to builtin PGx")
            else:
                warnings.append("PharmCAT not available (Java 17+ or JAR missing); falling back to builtin PGx")
        except ImportError:
            warnings.append("pharmcat_runner module not found; falling back to builtin PGx")

        logger.warning("PharmCAT not available, falling back to builtin PGx")

    # --- Builtin fallback -------------------------------------------------
    hits: list[dict] = []
    for v in variants:
        result = check_korean_pgx(v)
        if result:
            hits.append(
                {
                    "gene": result.gene,
                    "star_allele": result.star_allele,
                    "phenotype": result.phenotype,
                    "cpic_level": result.cpic_level,
                    "korean_prevalence": result.korean_prevalence,
                    "western_prevalence": result.western_prevalence,
                    "clinical_impact": result.clinical_impact,
                    "cpic_recommendation": result.cpic_recommendation,
                    "korean_flag": result.korean_flag,
                }
            )

    source = "builtin"
    if germline_vcf:
        # Germline was provided but PharmCAT was not available
        source = "builtin_limited"

    return {
        "pgx_hits": hits,
        "pgx_source": source,
        "pharmcat_version": "",
        "germline_provided": bool(germline_vcf),
        "warnings": warnings,
    }


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
