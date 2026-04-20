"""Variant selector for Clinical Board — tiered filter per AMP 2017 / ACMG 2015.

Applies deterministic selection criteria to cut the full variant list down to
a clinically-defensible subset for the Clinical Board agents. No LLM in the
selection path. All criteria implement _workspace/variant-selector/00_clinical_review.md
corrections (hotspot MUST, no "protein-impacting alone", soft caps, no fallback,
SF-gene exclusion in rare-disease).

Proxies and known gaps
----------------------
- **Variant-level OncoKB oncogenicity is not stored in BIKO yet.** We use
  ``is_cancer_gene`` + ``oncokb_level in {"1","2"}`` (gene-level actionable) as
  a pragmatic proxy for "oncogenic in CGC Tier 1 gene". Upgrade to true
  variant-level OncoKB curation is a future enhancement.
- **HPO scoring is an integer count, not Resnik-normalized.** The clinical
  review's target threshold is ``hpo_score_min = 0.3`` on a Resnik-normalized
  scale; we use ``hpo_score >= 1`` from the current HPO matcher until the
  matcher is upgraded.
- **ACMG SF v3.2 opt-in is not modeled in BIKO yet.** VUS in SF genes are
  silently excluded in rare-disease mode to avoid unintended secondary-findings
  reporting (Miller DT et al. 2023, PMID 37347242). When opt-in is added, move
  SF admission to a separate ``_select_sf_candidates`` path.

Reference: _workspace/variant-selector/00_clinical_review.md
"""

from __future__ import annotations

from typing import Optional

from scripts.enrichment.oncokb import is_cancer_gene
from scripts.common.config import get
from scripts.common.ddg2p_panel import is_admitted_neurodev_gene
from scripts.db.query_civic import extract_protein_position, is_hotspot

# ACMG SF v3.2 (Miller DT et al. 2023, Genet Med, PMID 37347242) —
# conservative partial list covering highly actionable cardiac, cancer, and
# metabolic conditions. VUS in these genes are silently excluded from the
# rare-disease MAY list until an explicit SF opt-in flow is modeled.
_ACMG_SF_V32 = frozenset(
    {
        "BRCA1",
        "BRCA2",
        "PALB2",
        "TP53",
        "STK11",
        "MLH1",
        "MSH2",
        "MSH6",
        "PMS2",
        "EPCAM",
        "APC",
        "MUTYH",
        "BMPR1A",
        "SMAD4",
        "VHL",
        "MEN1",
        "RET",
        "NF2",
        "TSC1",
        "TSC2",
        "PTEN",
        "DICER1",
        "RB1",
        "SDHAF2",
        "SDHB",
        "SDHC",
        "SDHD",
        "MAX",
        "TMEM127",
        "WT1",
        "KCNQ1",
        "KCNH2",
        "SCN5A",
        "RYR2",
        "CASQ2",
        "TRDN",
        "CALM1",
        "CALM2",
        "CALM3",
        "MYH7",
        "MYBPC3",
        "TNNI3",
        "TNNT2",
        "TPM1",
        "MYL2",
        "MYL3",
        "ACTC1",
        "PLN",
        "LMNA",
        "DSP",
        "PKP2",
        "DSG2",
        "DSC2",
        "TMEM43",
        "FBN1",
        "TGFBR1",
        "TGFBR2",
        "SMAD3",
        "ACTA2",
        "MYH11",
        "COL3A1",
        "LDLR",
        "APOB",
        "PCSK9",
        "RYR1",
        "CACNA1S",
        "ATP7B",
        "BTD",
        "OTC",
        "GAA",
        "HNF1A",
        "GLA",
        "TTR",
    }
)

_TSG_LOF_CONSEQUENCES = frozenset(
    {
        "stop_gained",
        "frameshift_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "start_lost",
    }
)

_INFRAME_INDEL_CONSEQUENCES = frozenset(
    {
        "inframe_insertion",
        "inframe_deletion",
    }
)

_PROTEIN_IMPACTING_CONSEQUENCES = frozenset(
    {
        "missense_variant",
        "stop_gained",
        "stop_lost",
        "frameshift_variant",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "start_lost",
        "inframe_insertion",
        "inframe_deletion",
    }
)

# Inverse of scripts/intake/parse_annotation.py::format_consequence. The real
# pipeline stores the BIKO-formatted short label (e.g. "Missense") on the
# variant dict, but the selector's gate is expressed in raw VEP SO terms. Both
# forms must resolve to the same membership check or Tier III MAY-arm variants
# silently disappear from the cancer board.
_FORMATTED_TO_SO = {
    "missense": "missense_variant",
    "nonsense": "stop_gained",
    "nonsense / stop gain": "stop_gained",
    "frameshift": "frameshift_variant",
    "splice donor": "splice_donor_variant",
    "splice acceptor": "splice_acceptor_variant",
    "in-frame deletion": "inframe_deletion",
    "in-frame insertion": "inframe_insertion",
    "synonymous": "synonymous_variant",
    "start loss": "start_lost",
    "stop loss": "stop_lost",
    "intronic": "intron_variant",
    "splice region": "splice_region_variant",
}


def _canonical_consequence(raw: str) -> str:
    """Return canonical VEP SO term for either a raw term or a BIKO-formatted
    label. Unknown terms pass through unchanged."""
    if not raw:
        return ""
    key = raw.strip().lower()
    if "&" in key:
        key = key.split("&", 1)[0].strip()
    return _FORMATTED_TO_SO.get(key, key)


# v2.2 B1: splice-region / synonymous variants are rescued (admitted) only when
# SpliceAI delta_max >= 0.2, per Tavtigian et al. 2023 PP3-moderate threshold.
_SPLICE_RESCUE_CONSEQUENCES = frozenset(
    {
        "synonymous_variant",
        "splice_region_variant",
    }
)
_SPLICEAI_RESCUE_THRESHOLD = 0.2

# v2.2 B2: MMR/Lynch panel carve-out — protein-impacting VUS in any of these
# genes are admitted to the cancer board regardless of hotspot / TSG / HPO
# status. Clinical rationale: Lynch syndrome screening + pembrolizumab MSI-H
# tumour-agnostic eligibility require MMR-VUS visibility even when no other
# driver signal is present. The B1 consequence gate still applies — there is
# no Lynch exception for intronic / UTR variants.
_MMR_LYNCH_GENES = frozenset({"MLH1", "MSH2", "MSH6", "PMS2", "EPCAM"})

_P_LP = frozenset({"Pathogenic", "Likely Pathogenic"})
_BENIGN = frozenset({"Benign", "Likely Benign"})
_MUST_TIERS = frozenset({"Tier I", "Tier II"})

# AMP 2017 Table 2 ordering (lower = earlier in output list). v2.2 B2 inserts
# VUS_MMR_Lynch between VUS_hotspot and VUS_TSG_LoF at priority 6. v1 de novo
# support (spec Q3.1) inserts ``VUS_denovo_neurodev`` and ``VUS_denovo_splice``
# at priority 5 — between VUS_hotspot (cancer MAY) and VUS_MMR_Lynch (6).
# Rationale: de novo in a DDG2P-admitted NDD gene is a stronger diagnostic
# signal than an MMR-panel VUS flagged for Lynch screening, but weaker than a
# cancer-specific hotspot VUS.
_REASON_PRIORITY = {
    "P_LP": 0,
    "Tier_I": 1,
    "Tier_II": 2,
    "Tier_III_hotspot": 3,
    "Tier_III_oncokb_gene": 4,
    "VUS_hotspot": 5,
    "VUS_denovo_neurodev": 5,
    "VUS_denovo_splice": 5,
    "VUS_MMR_Lynch": 6,
    "VUS_TSG_LoF": 7,
    "VUS_indel_hotspot": 8,
    "VUS_HPO_match": 9,
}

_MAY_REASONS_CANCER = (
    "VUS_hotspot",
    "VUS_MMR_Lynch",
    "VUS_TSG_LoF",
    "VUS_indel_hotspot",
)
_MAY_REASONS_RARE = (
    "VUS_HPO_match",
    "VUS_denovo_neurodev",
    "VUS_denovo_splice",
)

# v1 de novo support — pLI / missense-Z thresholds for the "constrained gene
# not yet in DDG2P" branch (spec Q3.1 OR-clause). Both required: pLI alone
# admits too many housekeeping genes per spec rationale.
_DENOVO_CONSTRAINT_PLI = 0.9
_DENOVO_CONSTRAINT_MISSENSE_Z = 3.09

# Inheritance labels that trigger the de novo rare-disease carve-out.
_DENOVO_INHERITANCE_LABELS = frozenset({"de_novo", "confirmed_de_novo"})


def _passes_consequence_gate(v: dict) -> bool:
    """v2.2 B1 protein-impacting consequence gate.

    Returns True if the variant's primary VEP consequence is a coding-effect
    class, or if it is a splice-region / synonymous variant rescued by
    SpliceAI delta_max >= 0.2 (read from ``v["in_silico"]["spliceai_max"]``).

    The ``P_LP`` must-reason branch bypasses this gate unconditionally
    (deep-intronic ClinVar-Pathogenic splice variants must still pass).
    """
    consequence = _canonical_consequence(v.get("consequence") or "")
    if consequence in _PROTEIN_IMPACTING_CONSEQUENCES:
        return True
    if consequence in _SPLICE_RESCUE_CONSEQUENCES:
        in_silico = v.get("in_silico") or {}
        raw = in_silico.get("spliceai_max") if isinstance(in_silico, dict) else None
        try:
            if raw is not None and float(raw) >= _SPLICEAI_RESCUE_THRESHOLD:
                return True
        except (TypeError, ValueError):
            return False
    return False


def select_board_variants(
    variants: list[dict],
    mode: str,
    report_data: Optional[dict] = None,
    config_overrides: Optional[dict] = None,
) -> tuple[list[dict], dict]:
    """Return (selected_variants_in_priority_order, selection_metadata).

    Parameters
    ----------
    variants : list[dict]
        Variant records as produced by ``orchestrate._build_variant_records``.
        Expected fields: ``gene``, ``classification``, ``tier``, ``hgvsp``,
        ``consequence``, ``cancer_gene_type``, ``oncokb_level``, ``hpo_score``,
        ``gnomad_af``, ``variant_type``.
    mode : str
        ``"cancer"`` or ``"rare-disease"``.
    report_data : dict, optional
        Full pipeline output (used for TMB footnote detection).
    config_overrides : dict, optional
        Per-call overrides, keys: ``max_cancer_board_variants``,
        ``max_rare_disease_board_variants``, ``max_cancer_may_include``,
        ``max_rare_disease_may_include``, ``rare_disease_vus_max_gnomad_af``.

    Returns
    -------
    selected, metadata : tuple[list[dict], dict]
        ``selected`` is a list of *copies* of the input variants, each with a
        ``selection_reason`` field. Order follows AMP 2017 Table 2 priority.
        ``metadata`` contains audit fields (see module docstring).
    """
    overrides = config_overrides or {}
    total_input = len(variants)

    if mode == "cancer":
        selected, must_n, may_n, excluded_n, truncated, n_dropped = _select_cancer(variants, overrides)
        empty_reason = "no Tier I/II variants, no ClinVar P/LP, no qualifying hotspot/driver VUS"
        criteria_summary = (
            "P/LP + AMP Tier I/II + Tier III hotspot/OncoKB-gene "
            "+ VUS hotspot/TSG-LoF/indel-hotspot (soft caps, no fallback)"
        )
    elif mode == "rare-disease":
        selected, must_n, may_n, excluded_n, truncated, n_dropped = _select_rare_disease(variants, overrides)
        empty_reason = "no P/LP variants, no HPO-matched VUS above threshold"
        criteria_summary = "ACMG P/LP + HPO-matched rare VUS (SF genes excluded, soft caps, no fallback)"
    else:
        raise ValueError(f"Unknown mode: {mode!r} (expected 'cancer' or 'rare-disease')")

    tmb_high = False
    if mode == "cancer" and report_data:
        tmb = report_data.get("tmb") or {}
        try:
            if float(tmb.get("score", 0) or 0) >= 10.0:
                tmb_high = True
        except (TypeError, ValueError):
            tmb_high = False

    by_reason: dict[str, int] = {}
    for v in selected:
        by_reason[v["selection_reason"]] = by_reason.get(v["selection_reason"], 0) + 1

    metadata = {
        "mode": mode,
        "total_input": total_input,
        "selected": len(selected),
        "must_included": must_n,
        "may_included": may_n,
        "excluded": excluded_n,
        "truncated": truncated,
        "n_dropped": n_dropped,
        # Hard cap is never applied (MUST items always pass); exposed for audit.
        "hard_cap_applied": False,
        "empty": len(selected) == 0,
        "empty_reason": empty_reason if len(selected) == 0 else "",
        "tmb_high_footnote": tmb_high,
        "criteria_summary": criteria_summary,
        "by_selection_reason": by_reason,
    }
    return selected, metadata


# ---------------------------------------------------------------------------
# Cancer mode
# ---------------------------------------------------------------------------


def _select_cancer(variants: list[dict], overrides: dict) -> tuple[list[dict], int, int, int, bool, int]:
    max_total = int(
        overrides.get(
            "max_cancer_board_variants",
            get("clinical_board.variant_selection.max_cancer_board_variants", 30),
        )
    )
    max_may = int(
        overrides.get(
            "max_cancer_may_include",
            get("clinical_board.variant_selection.max_cancer_may_include", 10),
        )
    )

    must: list[dict] = []
    may: list[dict] = []
    excluded = 0

    for v in variants:
        classification = v.get("classification", "")

        # Benign/LB always excluded from the cancer board regardless of tier.
        if classification in _BENIGN:
            excluded += 1
            continue

        reason = _cancer_must_reason(v)
        if reason:
            must.append(_tag(v, reason))
            continue

        may_reason = _cancer_may_reason(v)
        if may_reason:
            may.append(_tag(v, may_reason))
        else:
            excluded += 1

    # Soft cap: never truncate MUST. MAY gets truncated first to max_may,
    # then further if total exceeds max_total.
    n_may_total = len(may)
    may_sorted = sorted(may, key=lambda v: _REASON_PRIORITY[v["selection_reason"]])
    may_kept = may_sorted[:max_may]
    dropped_by_may_cap = n_may_total - len(may_kept)

    dropped_by_total_cap = 0
    if len(must) + len(may_kept) > max_total:
        allowed_may = max(0, max_total - len(must))
        dropped_by_total_cap = len(may_kept) - allowed_may
        may_kept = may_kept[:allowed_may]

    truncated = (dropped_by_may_cap + dropped_by_total_cap) > 0
    n_dropped = dropped_by_may_cap + dropped_by_total_cap

    combined = must + may_kept
    combined.sort(key=lambda v: (_REASON_PRIORITY[v["selection_reason"]],))
    return combined, len(must), len(may_kept), excluded, truncated, n_dropped


def _cancer_must_reason(v: dict) -> Optional[str]:
    classification = v.get("classification", "")
    tier = v.get("tier", "")

    # P/LP bypasses the v2.2 B1 consequence gate unconditionally —
    # a deep-intronic ClinVar-Pathogenic splice variant must still pass.
    if classification in _P_LP:
        return "P_LP"

    if not _passes_consequence_gate(v):
        return None

    if tier == "Tier I":
        return "Tier_I"
    if tier == "Tier II":
        return "Tier_II"
    if tier == "Tier III":
        gene = v.get("gene", "")
        pos = extract_protein_position(v.get("hgvsp", "") or "")
        if gene and pos is not None and is_hotspot(gene, pos):
            return "Tier_III_hotspot"
        if gene and is_cancer_gene(gene):
            level = str(v.get("oncokb_level", "") or "")
            if level in {"1", "2"}:
                return "Tier_III_oncokb_gene"
    return None


def _cancer_may_reason(v: dict) -> Optional[str]:
    if v.get("classification", "") != "VUS":
        return None

    if not _passes_consequence_gate(v):
        return None

    gene = v.get("gene", "")
    consequence = _canonical_consequence(v.get("consequence", "") or "")
    pos = extract_protein_position(v.get("hgvsp", "") or "")
    hotspot_hit = bool(gene) and pos is not None and is_hotspot(gene, pos)

    if hotspot_hit and consequence not in _INFRAME_INDEL_CONSEQUENCES:
        return "VUS_hotspot"

    # v2.2 B2: MMR/Lynch panel carve-out — any protein-impacting VUS in an
    # MMR gene is admitted regardless of hotspot / TSG / HPO status. Placed
    # after VUS_hotspot so a hotspot-positive MMR variant still takes the
    # higher-priority hotspot reason.
    if gene in _MMR_LYNCH_GENES:
        return "VUS_MMR_Lynch"

    if v.get("cancer_gene_type", "") == "TSG" and consequence in _TSG_LOF_CONSEQUENCES:
        return "VUS_TSG_LoF"

    if consequence in _INFRAME_INDEL_CONSEQUENCES and gene and is_cancer_gene(gene) and hotspot_hit:
        return "VUS_indel_hotspot"

    # NOTE: generic "protein-impacting consequence alone" is explicitly REJECTED
    # per clinical review — it floods the board with passenger missense.
    return None


# ---------------------------------------------------------------------------
# Rare-disease mode
# ---------------------------------------------------------------------------


def _select_rare_disease(variants: list[dict], overrides: dict) -> tuple[list[dict], int, int, int, bool, int]:
    max_total = int(
        overrides.get(
            "max_rare_disease_board_variants",
            get("clinical_board.variant_selection.max_rare_disease_board_variants", 20),
        )
    )
    max_may = int(
        overrides.get(
            "max_rare_disease_may_include",
            get("clinical_board.variant_selection.max_rare_disease_may_include", 10),
        )
    )
    max_af = float(
        overrides.get(
            "rare_disease_vus_max_gnomad_af",
            get("clinical_board.variant_selection.rare_disease_vus_max_gnomad_af", 0.01),
        )
    )

    must: list[dict] = []
    may: list[dict] = []
    excluded = 0

    for v in variants:
        classification = v.get("classification", "")
        if classification in _P_LP:
            must.append(_tag(v, "P_LP"))
            continue

        if classification != "VUS":
            excluded += 1
            continue

        gene = v.get("gene", "")
        if gene in _ACMG_SF_V32:
            # SF genes excluded from diagnostic VUS admission path.
            excluded += 1
            continue

        # v1 de novo carve-out (spec Q3.1 / Q3.2). Checked **before** the
        # HPO-score gate so a DDG2P-admitted de novo missense VUS reaches the
        # board even when no HPO terms were supplied. Collision point 3 from
        # the spec: if a later arm also matches, we still want the de novo
        # reason code to be recorded — the reason-priority ordering picks the
        # strongest admission reason, but the board sees both in the
        # ``selection_reason_list`` audit field. See _collect_rare_reasons.
        denovo_reason = _rare_denovo_reason(v)

        try:
            hpo_score = int(v.get("hpo_score", 0) or 0)
        except (TypeError, ValueError):
            hpo_score = 0

        af = v.get("gnomad_af")
        af_ok = True
        if af is not None:
            try:
                if float(af) >= max_af:
                    af_ok = False
            except (TypeError, ValueError):
                pass  # unparseable AF → treat as unknown (pass)

        if not af_ok:
            excluded += 1
            continue

        # De novo carve-out admission path. Bypasses the HPO-score gate but
        # still respects the AF gate and the B1 consequence gate (checked
        # inside _rare_denovo_reason).
        if denovo_reason is not None:
            tagged = _tag(v, denovo_reason)
            # Collision point 3: record the auxiliary HPO match so the board
            # renderer can surface both reasons without the selector silently
            # dropping one.
            if hpo_score >= 1:
                tagged["selection_reason_list"] = [denovo_reason, "VUS_HPO_match"]
            else:
                tagged["selection_reason_list"] = [denovo_reason]
            may.append(tagged)
            continue

        # Existing VUS + HPO admission path — unchanged.
        if hpo_score < 1:
            excluded += 1
            continue

        consequence = _canonical_consequence(v.get("consequence", "") or "")
        if consequence not in _PROTEIN_IMPACTING_CONSEQUENCES:
            excluded += 1
            continue

        may.append(_tag(v, "VUS_HPO_match"))

    n_may_total = len(may)
    may_kept = may[:max_may]
    dropped_by_may_cap = n_may_total - len(may_kept)

    dropped_by_total_cap = 0
    if len(must) + len(may_kept) > max_total:
        allowed_may = max(0, max_total - len(must))
        dropped_by_total_cap = len(may_kept) - allowed_may
        may_kept = may_kept[:allowed_may]

    truncated = (dropped_by_may_cap + dropped_by_total_cap) > 0
    n_dropped = dropped_by_may_cap + dropped_by_total_cap

    combined = must + may_kept
    combined.sort(key=lambda v: _REASON_PRIORITY[v["selection_reason"]])
    return combined, len(must), len(may_kept), excluded, truncated, n_dropped


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _tag(v: dict, reason: str) -> dict:
    """Return a shallow copy of ``v`` with ``selection_reason`` set."""
    out = dict(v)
    out["selection_reason"] = reason
    return out


# ---------------------------------------------------------------------------
# v1 de novo carve-out (spec Q3.1 / Q3.2)
# ---------------------------------------------------------------------------


def _is_denovo(v: dict) -> bool:
    """Return True if the variant record carries a de novo INFO flag.

    Accepts either a ``variant_inheritance`` key (canonical — populated by
    classify.build_variant_records from the Variant dataclass) or a
    ``confirmed_denovo`` truthy flag. The legacy ``inheritance`` key is NOT
    consulted because the rare-disease enrichment path uses that key for
    OMIM gene-level AD/AR/XL, and mixing the two semantics would cause
    false positives.
    """
    inh = v.get("variant_inheritance")
    if inh in _DENOVO_INHERITANCE_LABELS:
        return True
    if bool(v.get("confirmed_denovo")):
        return True
    return False


def _gene_is_constrained(gene: str) -> bool:
    """Return True if ``gene`` has pLI >= 0.9 AND missense Z >= 3.09.

    Sourced from the gnomAD v4.1 constraint metrics SQLite DB built by
    ``scripts/db/build_gnomad_constraint_db.py``. Both conditions are
    required per spec Q3.1 — pLI alone admits too many housekeeping genes,
    and the missense-Z gate matches Karczewski 2020 (PMID 32461654) top-1%
    constraint.

    Prior to v2.3-T8 this read from ``gene_knowledge.get_gene_info``,
    which never carried numeric pLI / missense_z fields → the OR-branch
    of the de novo carve-out (DDG2P-admitted OR constrained) was dead
    code. The constraint DB now lives in its own SQLite to keep the
    LLM-narrative gene_knowledge.json schema clean.

    Falls through to False (rather than raising) when the gnomAD
    constraint DB is unavailable — query_gnomad_constraint logs a single
    warning per run via its own availability cache.
    """
    if not gene:
        return False
    # Lazy import: keeps the variant_selector test fixture independent of the
    # db module so unit tests that only exercise the cancer arm need not
    # build a constraint DB stub.
    from scripts.db.query_gnomad_constraint import is_constrained

    try:
        return is_constrained(gene)
    except Exception:  # noqa: BLE001 — graceful degrade on DB errors
        return False


def _rare_denovo_reason(v: dict) -> Optional[str]:
    """Return the de novo carve-out reason code or ``None``.

    Gate order (spec Q3.1 / Q3.2):

    1. Variant has a de novo INFO flag (``_is_denovo``).
    2. Gene is DDG2P-admitted OR (pLI >= 0.9 AND missense Z >= 3.09).
    3. Consequence passes the B1 gate OR is rescued by SpliceAI >= 0.2.

    Protein-impacting de novo → ``"VUS_denovo_neurodev"``.
    Non-coding de novo rescued by SpliceAI → ``"VUS_denovo_splice"``.
    """
    if not _is_denovo(v):
        return None

    gene = v.get("gene", "") or ""
    if not is_admitted_neurodev_gene(gene) and not _gene_is_constrained(gene):
        return None

    consequence = _canonical_consequence(v.get("consequence", "") or "")

    # Protein-impacting path — straight admission under the B1 gate.
    if consequence in _PROTEIN_IMPACTING_CONSEQUENCES:
        if _passes_consequence_gate(v):
            return "VUS_denovo_neurodev"
        return None

    # Splice-region / synonymous rescue — B1 gate handles the SpliceAI check.
    if consequence in _SPLICE_RESCUE_CONSEQUENCES:
        if _passes_consequence_gate(v):
            return "VUS_denovo_splice"
        return None

    # Deep-intronic rescue — spec Q3.2 opens a narrow SpliceAI >= 0.2 window
    # for ``intron_variant`` that the generic B1 gate does not cover.
    if _passes_intron_splice_rescue(v):
        return "VUS_denovo_splice"

    return None


def _passes_intron_splice_rescue(v: dict) -> bool:
    """Return True if a de novo intronic variant has SpliceAI delta_max >= 0.2.

    Spec Q3.2 rescue path for deep-intronic de novo variants — the generic B1
    consequence gate does not admit ``intron_variant`` even with SpliceAI
    signal, so the de novo arm opens a narrow window here.
    """
    consequence = _canonical_consequence(v.get("consequence") or "")
    if consequence != "intron_variant":
        return False
    in_silico = v.get("in_silico") or {}
    raw = in_silico.get("spliceai_max") if isinstance(in_silico, dict) else None
    try:
        return raw is not None and float(raw) >= _SPLICEAI_RESCUE_THRESHOLD
    except (TypeError, ValueError):
        return False
