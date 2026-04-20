"""Shared TypedDict / type alias definitions for BIKO GenomeBoard.

This module is the single home for the dict-shaped contracts that flow
between pipeline stages — primarily the per-variant record dict written
by :func:`scripts.pipeline.classify.build_variant_records` and consumed
by the report templates, the Clinical Board, and the curators.

These TypedDicts are deliberately ``total=False`` because the pipeline
populates fields lazily depending on mode (cancer vs rare-disease) and
data availability (e.g. CIViC enrichment is cancer-only). Access through
``dict.get(...)`` remains the safest pattern at call sites.

Keep dataclass-shaped models (``Variant``, ``StructuralVariant``,
``AcmgEvidence``, ``ClassificationResult``, ``PgxResult``,
``CuratedTreatment``, ``AgentOpinion``, ``BoardOpinion``,
``CancerBoardOpinion``, ``PharmCATResult``) where they live today —
this module does not re-declare them as TypedDicts.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, TypedDict

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

# A domain sheet is the per-agent, per-mode formatted text block produced
# by :func:`scripts.clinical_board.domain_sheets.build_domain_sheet`. It is
# fed verbatim into the LLM prompt, so ``str`` is the contract.
DomainSheet = str

# HPO term hits keyed by gene symbol (output of HPO matcher).
HpoResults = Dict[str, List[str]]


# ---------------------------------------------------------------------------
# Curator + CIViC dict contracts
# ---------------------------------------------------------------------------


class CivicEvidence(TypedDict, total=False):
    """A single CIViC predictive-evidence row.

    Mirrors the per-row dict assembled by
    :func:`scripts.db.query_civic.get_predictive_evidence_for_tier`.
    """

    gene: str
    variant: str
    disease: str
    therapies: str
    evidence_type: str
    evidence_level: str
    evidence_direction: str
    significance: str
    pmid: str
    rating: float


class CivicEvidenceBundle(TypedDict, total=False):
    """Wrapper returned by :func:`get_predictive_evidence_for_tier`.

    ``match_level`` is one of ``"variant"``, ``"gene"``, ``"none"``.
    """

    match_level: str
    evidence: List[CivicEvidence]


# ---------------------------------------------------------------------------
# PGx dict contracts (output of get_pgx_results)
# ---------------------------------------------------------------------------


class PgxHit(TypedDict, total=False):
    """A single PGx finding row, used in ``report_data["pgx_hits"]``."""

    gene: str
    star_allele: str
    phenotype: str
    cpic_level: str
    korean_prevalence: float
    western_prevalence: float
    clinical_impact: str
    cpic_recommendation: str
    korean_flag: bool


class PgxResultsDict(TypedDict, total=False):
    """Return shape of :func:`scripts.pharmacogenomics.korean_pgx.get_pgx_results`."""

    pgx_hits: List[PgxHit]
    pgx_source: str  # "pharmcat" | "builtin" | "builtin_limited"
    pharmcat_version: str
    germline_provided: bool
    warnings: List[str]


# ---------------------------------------------------------------------------
# Board selection audit
# ---------------------------------------------------------------------------


class BoardSelectionMetadata(TypedDict, total=False):
    """Audit trail emitted by ``select_board_variants``."""

    total_input: int
    selected: int
    criteria_summary: str
    rejected_reasons: Dict[str, int]


# ---------------------------------------------------------------------------
# Per-variant record (the central pipeline ↔ report contract)
# ---------------------------------------------------------------------------


class VariantRecord(TypedDict, total=False):
    """Shape of a single dict in ``report_data["variants"]``.

    Built by :func:`scripts.pipeline.classify.build_variant_records` and
    consumed downstream by:

    * the report templates (``templates/cancer/report.html``,
      ``templates/rare-disease/report.html``),
    * the Clinical Board variant selector + curators,
    * the narrative scrubber, and
    * the knowledge-base sink.

    Fields are conditionally present:

    * Tier I/II/III/IV fields are populated only in cancer mode.
    * ``hpo_score``, ``omim_*``, ``inheritance``, ``clingen_validity`` are
      populated only in rare-disease mode.
    * ``selection_reason_list`` is populated when the variant selector
      admits the variant via a carve-out (de novo, HPO, hotspot, ...).
    """

    # ── Identity ──────────────────────────────────────────────────────
    variant: str  # canonical "chr:pos:ref>alt" id
    gene: Optional[str]
    chrom: str
    pos: int
    ref: str
    alt: str

    # ── Classification ────────────────────────────────────────────────
    classification: str  # "Pathogenic" | "Likely Pathogenic" | "VUS" | ...
    acmg_codes: List[str]
    conflict: bool
    clinvar_override: bool
    clinvar_override_reason: str
    original_engine_classification: Optional[str]
    clinvar_significance: str
    clinvar_id: Optional[str]
    review_status: Optional[str]

    # ── Population frequency ──────────────────────────────────────────
    gnomad_all: Optional[float]
    gnomad_eas: Optional[float]
    krgdb_freq: Optional[float]
    korea4k_freq: Optional[float]
    nard2_freq: Optional[float]
    korean_flag: str

    # ── Annotation ────────────────────────────────────────────────────
    hgvsc: str
    hgvsp: str
    consequence: str
    transcript: str
    impact: str
    sift: str
    polyphen: str
    in_silico: Dict[str, Any]

    # ── Trio / inheritance ────────────────────────────────────────────
    variant_inheritance: Optional[str]
    confirmed_denovo: bool
    selection_reason_list: List[str]
    source: str  # "primary" | "germline_inherited" | "germline_pgx"

    # ── Cancer-mode fields (AMP tiering + OncoKB) ─────────────────────
    tier: int
    tier_label: str
    tier_evidence_source: str
    civic_match_level: str
    civic_evidence: List[CivicEvidence]
    cancer_gene_type: str
    oncokb_level: str

    # ── Rare-disease enrichment ───────────────────────────────────────
    hpo_score: float
    matching_hpo: List[str]
    omim_mim: str
    omim_phenotypes: List[str]
    inheritance: str
    clingen_validity: str

    # ── Per-domain agent payload (clinical, korean_pop) ───────────────
    agents: Dict[str, Dict[str, Any]]


# ---------------------------------------------------------------------------
# Public re-exports
# ---------------------------------------------------------------------------

__all__ = [
    "BoardSelectionMetadata",
    "CivicEvidence",
    "CivicEvidenceBundle",
    "DomainSheet",
    "HpoResults",
    "PgxHit",
    "PgxResultsDict",
    "VariantRecord",
]
