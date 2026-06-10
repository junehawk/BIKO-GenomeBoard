"""ACMG classification, variant record assembly, and summary generation."""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Tuple

from scripts.classification.acmg_engine import (
    apply_clinvar_override,
    apply_hotspot_conflict_reconciliation,
    check_clinvar_conflict,
    classify_variant,
)
from scripts.enrichment.hpo_matcher import calculate_hpo_score, get_matching_hpo_terms
from scripts.enrichment.oncokb import get_cancer_gene_info
from scripts.enrichment.query_clingen import get_gene_validity
from scripts.enrichment.query_omim import query_omim
from scripts.common.config import get
from scripts.common.gene_knowledge import get_gene_info
from scripts.common.models import AcmgEvidence, StructuralVariant, Variant
from scripts.common.types import HpoResults, VariantRecord

logger = logging.getLogger(__name__)


def _resolve_ps1_pm5_exclusivity(evidences: List[AcmgEvidence]) -> List[AcmgEvidence]:
    """Enforce ACMG PS1 / PM5 mutual exclusivity on a variant's evidence list.

    PS1 (same amino-acid change as an established pathogenic variant) and PM5
    (a *different* change at a residue with a known pathogenic variant) describe
    incompatible facts — a variant cannot be both. When both co-occur (e.g. an
    InterVar-supplied PS1 alongside the engine's self-computed PM5), drop PS1 and
    keep PM5: PM5 preserves the A4 ClinVar-conflict reconciliation gate
    (PM1+PM5) and is the conservative (Moderate, not Strong) contribution, so a
    single contradictory pair can never manufacture a Pathogenic call
    (v2.7 review CLAS-04).
    """
    has_ps1 = any(e.code.upper().startswith("PS1") for e in evidences)
    has_pm5 = any(e.code.upper().startswith("PM5") for e in evidences)
    if has_ps1 and has_pm5:
        return [e for e in evidences if not e.code.upper().startswith("PS1")]
    return evidences


# ClinGen gene-validity classifications that count as an "established" gene-disease
# link for the PS2/PM6 de-novo gate (T1-08). Limited / Disputed / Refuted / No Known
# Disease Relationship are NOT established.
_ESTABLISHED_VALIDITY = {"definitive", "strong", "moderate"}


def _gene_has_disease_association(gene: Optional[str], cache: Dict[str, bool]) -> bool:
    """True when ``gene`` has an established gene-disease association (ClinGen validity
    Definitive/Strong/Moderate, or an OMIM disease phenotype). Memoised in ``cache``
    for the current batch. Used to gate PS2/PM6 (T1-08)."""
    if not gene:
        return False
    if gene in cache:
        return cache[gene]
    result = False
    validity = get_gene_validity(gene)
    if validity and validity.strip().lower() in _ESTABLISHED_VALIDITY:
        result = True
    else:
        omim = query_omim(gene)
        if isinstance(omim, dict) and omim.get("phenotypes"):
            result = True
    cache[gene] = result
    return result


def classify_variants(
    variants: List[Variant],
    db_results: Dict[str, Dict[str, Any]],
    freq_results: Dict[str, Dict[str, Any]],
    intervar_data: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Run ACMG classification for a list of variants.

    Returns classification_results dict keyed by variant_id.
    """
    # Lazy imports for optional modules
    try:
        from scripts.classification.in_silico import (
            generate_pp3_bp4,
            get_configured_thresholds,
            parse_in_silico_from_csq,
        )

        _has_in_silico = True
        _in_silico_thresholds = get_configured_thresholds()  # config-tunable PP3/BP4 (T5-02)
    except ImportError:
        _has_in_silico = False
        _in_silico_thresholds = None
    try:
        from scripts.classification.evidence_collector import (
            collect_additional_evidence,
            collect_denovo_evidence,
        )

        _has_evidence_collector = True
    except ImportError:
        _has_evidence_collector = False
    try:
        from scripts.annotation.parse_intervar import get_intervar_evidence

        _has_intervar = True
    except ImportError:
        _has_intervar = False

    # v2.3-T6: Batch-prefetch ClinVar protein positions for self-computed PM5.
    # The ClinVar local DB now carries an hgvsp column (build_clinvar_db.py
    # schema_version=2). For every distinct gene in this batch we ask
    # query_local_clinvar.get_clinvar_pathogenic_positions for the set of
    # residues with at least one ClinVar P/LP entry, and pass that set into
    # collect_additional_evidence per variant. evidence_collector then fires
    # PM5 when the variant's protein position is in the set, which in turn
    # unblocks the full A4 ClinVar-conflict override gate (PM1+PM5) on real
    # data. Falls back gracefully to an empty set if the DB lacks the column
    # (legacy build) or is absent — the InterVar upstream PM5 path still
    # works in that case.
    try:
        from scripts.storage.query_local_clinvar import (
            get_clinvar_pathogenic_changes,
            get_clinvar_pathogenic_positions,
        )

        _has_clinvar_pm5 = True
    except ImportError:
        _has_clinvar_pm5 = False

    if _has_clinvar_pm5:
        _genes = {v.gene for v in variants if getattr(v, "gene", None)}
        clinvar_pos_by_gene: dict = {g: get_clinvar_pathogenic_positions(g) for g in _genes}
        # T1-07: per-residue pathogenic AA changes enable the refined PM5 gate
        # (different-AA-only). Falls back to position-only when unavailable.
        clinvar_changes_by_gene: dict = {g: get_clinvar_pathogenic_changes(g) for g in _genes}
    else:
        clinvar_pos_by_gene = {}
        clinvar_changes_by_gene = {}

    classification_results = {}
    _denovo_assoc_cache: Dict[str, bool] = {}  # per-batch memo for the PS2/PM6 gene gate (T1-08)
    for variant in variants:
        db = db_results[variant.variant_id]
        freq = freq_results[variant.variant_id]

        evidences = []
        for code in db["clinvar"].get("acmg_codes", []):
            evidences.append(AcmgEvidence(code=code, source="clinvar", description=""))
        for code in freq.get("acmg_codes", []):
            evidences.append(AcmgEvidence(code=code, source="freq_comparison", description=""))

        # In silico PP3/BP4 evidence
        if _has_in_silico and variant.in_silico:
            pp3_bp4_codes = generate_pp3_bp4(
                parse_in_silico_from_csq(variant.in_silico), thresholds=_in_silico_thresholds
            )
            for code in pp3_bp4_codes:
                evidences.append(AcmgEvidence(code=code, source="in_silico", description=""))

        # Additional evidence from variant annotations (PVS1, PM1, PM4, PP2, BP7, etc.)
        # gene_info (pLI / missense_z) is required for PVS1/PVS1_Strong/PP2/BP1 to fire —
        # without it those codes were silently dead on real data.
        # clinvar_pathogenic_positions enables the self-computed PM5 path: the
        # ClinVar local DB now carries an hgvsp column (v2.3-T6), so we look up
        # the set of residues in this gene with a known ClinVar P/LP entry and
        # let evidence_collector fire PM5 when the current variant's position
        # is in the set. Without this set PM5 only reaches A4 via the InterVar
        # upstream path when intervar_data is supplied; with it both paths are
        # live. clinvar_pos_by_gene falls back to an empty set on legacy DBs.
        if _has_evidence_collector:
            gene_info = get_gene_info(variant.gene) if variant.gene else None
            clinvar_pos = clinvar_pos_by_gene.get(getattr(variant, "gene", None)) or set()
            clinvar_changes = clinvar_changes_by_gene.get(getattr(variant, "gene", None)) or {}
            extra_codes = collect_additional_evidence(
                variant,
                gene_info=gene_info,
                clinvar_pathogenic_positions=clinvar_pos,
                clinvar_pathogenic_changes=clinvar_changes,
            )
            for code in extra_codes:
                evidences.append(AcmgEvidence(code=code, source="evidence_collector", description=""))

            # v1 de novo support — PS2/PM6 from trio INFO flags.
            # Kept as a separate call per spec Q4 so PS2/PM6 can never leak
            # into the ClinVar-conflict override path in acmg_engine.
            # T1-08: PS2/PM6 require an established gene-disease link — compute it
            # only for de-novo-flagged variants (rare) to avoid per-variant lookups.
            _denovo_flag = getattr(variant, "inheritance", None) in (
                "de_novo",
                "confirmed_de_novo",
            ) or bool(getattr(variant, "confirmed_denovo", False))
            _gene_assoc = _gene_has_disease_association(variant.gene, _denovo_assoc_cache) if _denovo_flag else True
            denovo_codes = collect_denovo_evidence(variant, gene_has_disease_association=_gene_assoc)
            for code in denovo_codes:
                evidences.append(AcmgEvidence(code=code, source="denovo", description=""))

        # InterVar evidence (if provided)
        if _has_intervar and intervar_data:
            intervar_codes = get_intervar_evidence(variant, intervar_data)
            for code in intervar_codes:
                # Avoid duplicates — InterVar may overlap with self-collected evidence
                existing_codes = {e.code for e in evidences}
                if code not in existing_codes:
                    evidences.append(AcmgEvidence(code=code, source="intervar", description=""))

        evidences = _resolve_ps1_pm5_exclusivity(evidences)
        classification = classify_variant(evidences, gene=variant.gene)
        clinvar_sig = db["clinvar"].get("clinvar_significance", "Not Found")
        review_status = db["clinvar"].get("review_status", "")
        conflict = check_clinvar_conflict(classification.classification, clinvar_sig)
        classification.conflict = conflict

        # A4: narrow hotspot+PM5 reconciliation when ClinVar is Conflicting.
        # Fires only when all four gate conditions hold (engine ≥LP, ClinVar
        # Conflicting, PM1, PM5). Populates classification.clinvar_override_reason
        # when it fires; leaves the engine classification untouched. When the
        # reason is populated we skip the legacy apply_clinvar_override so the
        # conflict cannot silently downgrade the already-earned LP/P.
        apply_hotspot_conflict_reconciliation(
            classification,
            clinvar_significance=clinvar_sig,
            gene=getattr(variant, "gene", "") or "",
            hgvsp=getattr(variant, "hgvsp", "") or "",
        )

        original_classification = classification.classification
        if classification.clinvar_override_reason:
            final_classification = original_classification
        else:
            final_classification = apply_clinvar_override(original_classification, clinvar_sig, review_status)
        if final_classification != original_classification:
            classification.classification = final_classification
            classification.clinvar_override = True
            classification.original_engine_classification = original_classification
        else:
            classification.clinvar_override = False
            classification.original_engine_classification = None

        classification_results[variant.variant_id] = classification

    return classification_results


def build_variant_records(
    variants: List[Variant],
    db_results: Dict[str, Dict[str, Any]],
    freq_results: Dict[str, Dict[str, Any]],
    classification_results: Dict[str, Any],
    mode: str,
    hpo_results: HpoResults,
) -> List[VariantRecord]:
    """Build the variant_records list from per-variant results."""
    from scripts.storage.query_civic import get_predictive_evidence_for_tier
    from scripts.somatic.amp_tiering import amp_assign_tier

    variant_records = []
    for variant in variants:
        db = db_results[variant.variant_id]
        freq = freq_results[variant.variant_id]
        classification = classification_results[variant.variant_id]
        clinvar = db["clinvar"]

        variant_records.append(
            {
                "variant": variant.variant_id,
                "gene": variant.gene,
                "chrom": variant.chrom,
                "pos": variant.pos,
                "ref": variant.ref,
                "alt": variant.alt,
                "classification": classification.classification,
                "acmg_codes": classification.evidence_codes,
                "conflict": classification.conflict,
                "clinvar_override": getattr(classification, "clinvar_override", False),
                "clinvar_override_reason": getattr(classification, "clinvar_override_reason", "") or "",
                "original_engine_classification": getattr(classification, "original_engine_classification", None),
                "clinvar_significance": clinvar.get("clinvar_significance", "Not Found"),
                "clinvar_id": clinvar.get("clinvar_id"),
                "review_status": clinvar.get("review_status"),
                "gnomad_all": db["gnomad"].get("gnomad_all"),
                "gnomad_eas": db["gnomad"].get("gnomad_eas"),
                "kova_freq": db.get("kova_freq"),
                "kova_homozygote": db.get("kova_homozygote"),
                "korean_flag": freq.get("korean_flag", ""),
                # Annotation fields
                "hgvsc": variant.hgvsc or "",
                "hgvsp": variant.hgvsp or "",
                "consequence": variant.consequence or "",
                "transcript": variant.transcript or "",
                "impact": variant.impact or "",
                "sift": variant.sift or "",
                "polyphen": variant.polyphen or "",
                # In silico prediction scores
                "in_silico": variant.in_silico or {},
                # Trio-aware de novo context (populated from VCF INFO flags by
                # parse_vcf). Kept under a ``variant_inheritance`` key because
                # the rare-disease enrichment path below writes a gene-level
                # ``inheritance`` field (from OMIM AD/AR/XL) — the two are
                # semantically different and must not collide.
                "variant_inheritance": getattr(variant, "inheritance", None),
                "confirmed_denovo": bool(getattr(variant, "confirmed_denovo", False)),
                # v2.4 QW-A: selection reason audit channel. The
                # clinical_board.variant_selector populates this list when a
                # variant hits the de novo / HPO / hotspot carve-outs
                # (VUS_denovo_neurodev, VUS_denovo_splice, VUS_HPO_match, ...).
                # Defaults to [] so legacy Variant objects without the attribute
                # render an empty badge row in the template instead of raising
                # UndefinedError.
                "selection_reason_list": list(getattr(variant, "selection_reason_list", []) or []),
                # Variant origin badge: "primary" (main VCF), "germline_inherited"
                # (from --germline target extraction), or "germline_pgx".
                "source": getattr(variant, "source", None) or "primary",
                "agents": {
                    "clinical": clinvar,
                    "korean_pop": {
                        "gnomad_all": db["gnomad"].get("gnomad_all"),
                        "gnomad_eas": db["gnomad"].get("gnomad_eas"),
                        "kova_freq": db.get("kova_freq"),
                        "kova_homozygote": db.get("kova_homozygote"),
                        "korean_flag": freq.get("korean_flag", ""),
                        "api_available": db["gnomad"].get("api_available", False),
                    },
                },
            }
        )

    # Assign tiers
    if mode == "cancer":
        tiering_strategy = get("somatic.tiering_strategy", "B")
    else:
        tiering_strategy = "C"

    for v_result in variant_records:
        gene = v_result.get("gene", "")
        cls = v_result.get("classification", "VUS")
        hgvsp = v_result.get("hgvsp", "")

        if tiering_strategy in ("A", "B"):
            civic_evidence = get_predictive_evidence_for_tier(gene, hgvsp)
        else:
            civic_evidence = None

        tier_result = amp_assign_tier(
            classification=cls,
            gene=gene,
            hgvsp=hgvsp,
            strategy=tiering_strategy,
            civic_evidence=civic_evidence,
        )
        v_result["tier"] = tier_result.tier
        v_result["tier_label"] = tier_result.tier_label
        v_result["tier_evidence_source"] = tier_result.evidence_source
        v_result["civic_match_level"] = tier_result.civic_match_level
        v_result["civic_evidence"] = tier_result.civic_evidence

        cancer_info = get_cancer_gene_info(gene)
        if cancer_info:
            v_result["cancer_gene_type"] = cancer_info.get("type", "")
            v_result["oncokb_level"] = cancer_info.get("level", "")
        else:
            v_result["cancer_gene_type"] = ""
            v_result["oncokb_level"] = ""

    # Rare disease enrichment
    if mode == "rare-disease":
        for v_result in variant_records:
            gene = v_result.get("gene", "")
            v_result["hpo_score"] = calculate_hpo_score(gene, hpo_results)
            v_result["matching_hpo"] = get_matching_hpo_terms(gene, hpo_results)
            omim = query_omim(gene)
            if omim:
                v_result["omim_mim"] = omim.get("mim", "")
                v_result["omim_phenotypes"] = omim.get("phenotypes", [])
                v_result["inheritance"] = omim.get("inheritance", "")
            else:
                v_result["omim_mim"] = ""
                v_result["omim_phenotypes"] = []
                v_result["inheritance"] = ""
            v_result["clingen_validity"] = get_gene_validity(gene) or ""

    return variant_records


def build_summary(variant_records: List[VariantRecord]) -> Dict[str, int]:
    """Build summary statistics from variant records."""
    summary = {
        "total": len(variant_records),
        "pathogenic": 0,
        "likely_pathogenic": 0,
        "vus": 0,
        "likely_benign": 0,
        "benign": 0,
        "drug_response": 0,
        "risk_factor": 0,
    }
    for v in variant_records:
        cls = v.get("classification", "VUS").lower()
        if cls == "pathogenic":
            summary["pathogenic"] += 1
        elif cls == "likely pathogenic":
            summary["likely_pathogenic"] += 1
        elif cls == "drug response":
            summary["drug_response"] += 1
        elif cls == "risk factor":
            summary["risk_factor"] += 1
        elif cls == "likely benign":
            summary["likely_benign"] += 1
        elif cls == "benign":
            summary["benign"] += 1
        else:
            summary["vus"] += 1
    return summary


def sv_to_dict(sv: StructuralVariant) -> Dict[str, Any]:
    """Convert StructuralVariant to template-friendly dict."""
    return {
        "annotsv_id": sv.annotsv_id,
        "chrom": sv.chrom,
        "start": sv.start,
        "end": sv.end,
        "length": sv.length,
        "sv_type": sv.sv_type,
        "sample_id": sv.sample_id,
        "acmg_class": sv.acmg_class,
        "acmg_label": sv.acmg_label,
        "acmg_scored": sv.acmg_scored,
        "ranking_score": sv.ranking_score,
        "cytoband": sv.cytoband,
        "gene_name": sv.gene_name,
        "gene_count": sv.gene_count,
        "gene_details": sv.gene_details,
        "size_display": sv.size_display,
        "phenotypes": sv.phenotypes,
        "evidence_source": sv.evidence_source,
        "p_gain_phen": sv.p_gain_phen,
        "p_loss_phen": sv.p_loss_phen,
        "p_gain_hpo": sv.p_gain_hpo,
        "p_loss_hpo": sv.p_loss_hpo,
        "b_gain_af_max": sv.b_gain_af_max,
        "b_loss_af_max": sv.b_loss_af_max,
        "omim_morbid": sv.omim_morbid,
        "is_pathogenic": sv.is_pathogenic,
    }


def split_variants_for_display(
    variant_records: List[VariantRecord],
    hide_vus: bool,
    mode: str = "cancer",
) -> Tuple[
    List[VariantRecord],
    List[VariantRecord],
    List[VariantRecord],
    int,
    List[VariantRecord],
    List[VariantRecord],
]:
    """Split variant records into detailed and omitted lists for template rendering.

    Returns (tier1, tier2, tier3, tier4_count, detailed, omitted).
    """
    tier1_variants = [v for v in variant_records if v.get("tier") == 1]
    tier2_variants = [v for v in variant_records if v.get("tier") == 2]
    tier3_variants = [v for v in variant_records if v.get("tier") == 3]
    tier4_count = sum(1 for v in variant_records if v.get("tier") == 4)

    if hide_vus:
        _significant = {"pathogenic", "likely pathogenic", "drug response", "risk factor"}
        detailed_variants = [
            v for v in variant_records if v.get("tier") in (1, 2) or v.get("classification", "").lower() in _significant
        ]
        omitted_variants = [v for v in variant_records if v not in detailed_variants]
    else:
        detailed_variants = list(variant_records)
        omitted_variants = []

    # Sort detail pages by tier then classification rank
    _tier_sort = {1: 0, 2: 1, 3: 2, 4: 3}
    _cls_sort = {
        "pathogenic": 0,
        "likely pathogenic": 1,
        "drug response": 2,
        "risk factor": 3,
        "vus": 4,
        "likely benign": 5,
        "benign": 6,
    }
    if mode == "rare-disease":
        # Rare-disease reports are organised by ACMG classification, not AMP tier
        # (tier is a borrowed cancer concept on this path), so a Pathogenic disease
        # variant must precede a Drug-Response / Risk-Factor PGx finding even though
        # the latter is assigned AMP Tier I. Sort classification-first.
        detailed_variants.sort(
            key=lambda v: (
                _cls_sort.get(v.get("classification", "VUS").lower(), 4),
                _tier_sort.get(v.get("tier", 4), 4),
            )
        )
    else:
        detailed_variants.sort(
            key=lambda v: (
                _tier_sort.get(v.get("tier", 4), 4),
                _cls_sort.get(v.get("classification", "VUS").lower(), 4),
            )
        )

    return tier1_variants, tier2_variants, tier3_variants, tier4_count, detailed_variants, omitted_variants
