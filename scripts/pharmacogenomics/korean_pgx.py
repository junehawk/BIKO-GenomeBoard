# scripts/pharma/korean_pgx.py
"""Korean-population pharmacogenomics (PGx) analysis.

Provides two levels of API:

* ``check_korean_pgx(variant)`` -- per-variant SNV matching against the
  local ``korean_pgx_table.json`` (builtin engine, always available).
* ``get_pgx_results(variants, germline_vcf, config)`` -- unified entry
  point that prefers PharmCAT when a germline VCF is supplied and
  PharmCAT is installed, falling back to the builtin engine otherwise.

Phenotype strings are fully JSON-driven via the ``default_phenotype``
field on each gene entry in ``data/korean_pgx_table.json``. The builtin
engine no longer carries any hardcoded gene → phenotype elif chain;
editing the JSON is sufficient to extend or correct the curated set.
"""

from __future__ import annotations

import json
import logging
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional

from scripts.common.config import get
from scripts.common.models import PgxResult, Variant
from scripts.common.types import PgxHit, PgxResultsDict

logger = logging.getLogger(__name__)

_PGX_DATA = None
_PGX_LOCK = threading.Lock()

# Fallback PGx gene set used only when the JSON table cannot be loaded
# (file missing, malformed JSON, etc.). Keeps check_korean_pgx() callable
# during degraded operation instead of raising.
_FALLBACK_PGX_GENES = [
    "CYP2D6",
    "CYP2C19",
    "CYP2C9",
    "HLA-B",
    "HLA-A",
    "NUDT15",
    "TPMT",
    "DPYD",
    "CYP3A5",
    "UGT1A1",
    "SLCO1B1",
    "VKORC1",
    "CYP1A2",
    "G6PD",
    "IFNL3",
    "CYP2B6",
    "CYP4F2",
    "ABCG2",
    "NAT2",
    "CACNA1S",
    "CFTR",
    "CYP3A4",
    "MT-RNR1",
    "RYR1",
]


def _load_pgx_data() -> List[Dict[str, Any]]:
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


def _get_pgx_genes() -> set[str]:
    """Derive the active PGx gene set.

    Priority:
      1. Genes present in ``data/korean_pgx_table.json`` (authoritative).
      2. ``pgx.genes`` config list (legacy override).
      3. Hardcoded 24-gene fallback list (only when both above fail).
    """
    pgx_data = _load_pgx_data()
    if pgx_data:
        return {entry["gene"] for entry in pgx_data if entry.get("gene")}
    configured = get("pgx.genes", _FALLBACK_PGX_GENES)
    return set(configured) if configured else set(_FALLBACK_PGX_GENES)


PGX_GENES = _get_pgx_genes()


def check_korean_pgx(variant: Variant) -> Optional[PgxResult]:
    """Return a :class:`PgxResult` if *variant* matches a known PGx gene.

    Phenotype is sourced from the JSON ``default_phenotype`` field; no
    hardcoded gene → phenotype mapping remains in this module.
    """
    if variant.gene not in PGX_GENES:
        return None

    pgx_data = _load_pgx_data()
    for entry in pgx_data:
        if entry["gene"] == variant.gene:
            return PgxResult(
                gene=variant.gene,
                star_allele=entry.get("variant", ""),
                phenotype=entry.get("default_phenotype", ""),
                cpic_level=entry.get("cpic_level", ""),
                korean_prevalence=entry.get("korean_freq", 0),
                western_prevalence=entry.get("western_freq", 0),
                clinical_impact=entry.get("clinical_impact", ""),
                cpic_recommendation="",
            )
    return None


# ---------------------------------------------------------------------------
# PharmCAT → PgxResult conversion
# ---------------------------------------------------------------------------


def _convert_pharmcat_to_pgx_hits(pharmcat_result: Any) -> List[PgxHit]:
    """Convert a :class:`PharmCATResult` into the legacy pgx_hits list format.

    PharmCAT provides accurate diplotype + phenotype from the actual germline
    genotype. The builtin PGx table (``korean_pgx_table.json``) provides
    CPIC level, Korean/Western population prevalence, and clinical impact
    that PharmCAT does not carry. We merge both sources: PharmCAT for the
    genotype call, builtin for the population/guideline metadata.

    Phenotype priority: PharmCAT diplotype-derived phenotype >
    builtin ``default_phenotype`` > empty string.
    """
    # Build a gene→entry lookup from the builtin PGx table so we can
    # enrich PharmCAT results with CPIC level + Korean prevalence.
    builtin_by_gene: dict[str, dict] = {}
    for entry in _load_pgx_data():
        builtin_by_gene[entry["gene"]] = entry

    hits: list[dict] = []
    for gene, diplotype in pharmcat_result.diplotypes.items():
        phenotype = pharmcat_result.phenotypes.get(gene, "")

        # Skip genes with no meaningful result (Unknown/Unknown, No Result)
        if phenotype in ("No Result", "") and "Unknown" in diplotype:
            continue

        # Merge with builtin metadata
        builtin = builtin_by_gene.get(gene, {})
        cpic_level = builtin.get("cpic_level", "")
        korean_prev = builtin.get("korean_freq", 0.0)
        western_prev = builtin.get("western_freq", 0.0)
        clinical_impact = builtin.get("clinical_impact", "")

        # Phenotype precedence: PharmCAT diplotype-derived > builtin default > ""
        if not phenotype:
            phenotype = builtin.get("default_phenotype", "")

        # Drug recommendations from PharmCAT (if available)
        gene_drugs = [r for r in pharmcat_result.drug_recommendations if gene in r.get("gene", "")]
        cpic_rec = ""
        if gene_drugs:
            cpic_rec = gene_drugs[0].get("guideline", "")
            if not clinical_impact:
                clinical_impact = gene_drugs[0].get("classification", "")

        # Korean flag: set if this gene has notable Korean-specific prevalence
        korean_flag = korean_prev > 0 and abs(korean_prev - western_prev) / max(western_prev, 0.01) > 0.5

        hits.append(
            {
                "gene": gene,
                "star_allele": diplotype,
                "phenotype": phenotype,
                "cpic_level": cpic_level,
                "korean_prevalence": korean_prev,
                "western_prevalence": western_prev,
                "clinical_impact": clinical_impact,
                "cpic_recommendation": cpic_rec,
                "korean_flag": korean_flag,
            }
        )
    return hits


# ---------------------------------------------------------------------------
# Unified PGx entry point
# ---------------------------------------------------------------------------


def get_pgx_results(
    variants: List[Variant],
    germline_vcf: Optional[str] = None,
    config: Optional[Dict[str, Any]] = None,
) -> PgxResultsDict:
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
            from scripts.pharmacogenomics.pharmcat_runner import is_pharmcat_available, run_pharmcat

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
    hits = []  # type: List[PgxHit]
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
        print(json.dumps({"error": "Usage: python -m scripts.pharmacogenomics.korean_pgx 'chr10:96541616 G>A' [gene]"}))
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
