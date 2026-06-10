"""Grounding scrubber (v2.8) — annotate factual entities not in the briefing.

The narrative_scrubber enforces the curate-then-narrate contract for *drugs*
(every treatment row must trace to a curated_id). Empirically (v2.8 实测) the
LLM also fabricates *genes* the case does not contain — e.g. it headlined a
"KRAS G12D driver" on a KRAS-absent case. Prompt-level grounding clauses help
but do not fully stop this under generation pressure, so this module adds the
same deterministic enforcement BIKO already uses for drugs, extended to genes.

v2.8 is **annotate-only**: it does NOT strip or rewrite the narrative. It scans
the board opinion's prose for known gene symbols that are NOT in the case variant
set and appends a non-fatal note to ``opinion.grounding_flags`` so the reviewing
researcher sees "these genes were named but are not in the case — verify". This
is the conservative first step; a stricter strip/down-confidence mode can layer
on later.
"""

from __future__ import annotations

import logging
import re
from typing import Any, Iterable, List, Set

logger = logging.getLogger(__name__)

# High-precision vocabulary of gene symbols an LLM is likely to name in a cancer
# or rare-disease board. Kept curated (not a full HGNC dump) so the intersection
# with free text stays high-signal — a token only flags if it is BOTH a known
# gene here AND absent from the case. Unioned at runtime with gene_knowledge
# keys when that table is available (best-effort).
_DRIVER_GENES: Set[str] = {
    # Cancer drivers / oncogenes & TSGs commonly hallucinated
    "KRAS",
    "NRAS",
    "HRAS",
    "BRAF",
    "EGFR",
    "ALK",
    "ROS1",
    "MET",
    "RET",
    "ERBB2",
    "ERBB3",
    "PIK3CA",
    "PIK3R1",
    "AKT1",
    "PTEN",
    "MTOR",
    "TP53",
    "MDM2",
    "RB1",
    "CDKN2A",
    "CCND1",
    "CDK4",
    "MYC",
    "MYCN",
    "KIT",
    "PDGFRA",
    "FLT3",
    "KMT2A",
    "NPM1",
    "JAK2",
    "IDH1",
    "IDH2",
    "SMAD4",
    "APC",
    "CTNNB1",
    "VHL",
    "NF1",
    "NF2",
    "FGFR1",
    "FGFR2",
    "FGFR3",
    "FGFR4",
    "KEAP1",
    "STK11",
    "ARID1A",
    "ARID1B",
    "SMARCA4",
    "BAP1",
    "MEN1",
    "GNAS",
    "GNAQ",
    "GNA11",
    "DDR2",
    "MAP2K1",
    "NTRK1",
    "NTRK2",
    "NTRK3",
    "ESR1",
    "AR",
    "POLE",
    "POLD1",
    # DNA-repair / HBOC / Lynch
    "BRCA1",
    "BRCA2",
    "PALB2",
    "ATM",
    "CHEK2",
    "BARD1",
    "BRIP1",
    "RAD51C",
    "RAD51D",
    "MLH1",
    "MSH2",
    "MSH6",
    "PMS2",
    "EPCAM",
    "MUTYH",
    # Cardiomyopathy / arrhythmia / connective tissue (rare disease)
    "MYH7",
    "MYBPC3",
    "TNNT2",
    "TNNI3",
    "LMNA",
    "KCNQ1",
    "KCNH2",
    "SCN5A",
    "RYR2",
    "FBN1",
    "TGFBR1",
    "TGFBR2",
    "COL3A1",
    # Neurodevelopmental / metabolic / other common rare-disease genes
    "CFTR",
    "DMD",
    "SMN1",
    "HTT",
    "MECP2",
    "FMR1",
    "SCN1A",
    "CHD7",
    "CHD8",
    "PTPN11",
    "SOS1",
    "RAF1",
    "CDH1",
    "RET",
}

_GENE_TOKEN_RE = re.compile(r"\b[A-Z][A-Z0-9]{1,9}\b")

# Uppercase tokens that look gene-like but are clinical abbreviations / never a
# gene we want to flag. Defensive — the known-gene gate already removes most.
_NON_GENE_TOKENS: Set[str] = {
    "TMB",
    "MSI",
    "HRD",
    "DNA",
    "RNA",
    "VUS",
    "ACMG",
    "AMP",
    "CNV",
    "SV",
    "SNV",
    "PARP",
    "TKI",
    "NSCLC",
    "PDAC",
    "CNS",
    "FDA",
    "NCCN",
    "ESMO",
    "CPIC",
    "HPO",
    "OMIM",
    "PMID",
    "NCT",
    "WES",
    "WGS",
    "VAF",
    "LOH",
    "CR",
    "PR",
    "OS",
    "PFS",
    "HER2",
    "PD",
    "PDL1",
    "IHC",
    "FISH",
    "AA",
    "ID",
}

_KNOWN_GENES_CACHE: Set[str] | None = None


def _known_genes() -> Set[str]:
    """Curated driver set ∪ gene_knowledge keys (best-effort), all upper-cased."""
    global _KNOWN_GENES_CACHE
    if _KNOWN_GENES_CACHE is not None:
        return _KNOWN_GENES_CACHE
    genes = set(_DRIVER_GENES)
    try:
        from scripts.common.gene_knowledge import _load_knowledge

        genes |= {str(g).upper() for g in (_load_knowledge() or {}).keys() if g}
    except Exception:  # noqa: BLE001 — gene_knowledge is optional; curated set still works
        pass
    _KNOWN_GENES_CACHE = genes
    return genes


def _briefing_genes(report_data: dict) -> Set[str]:
    """Upper-cased set of gene symbols actually present in the case."""
    out: Set[str] = set()
    for key in ("_board_variants", "variants"):
        for v in report_data.get(key, []) or []:
            gene = (v.get("gene") if isinstance(v, dict) else getattr(v, "gene", None)) or ""
            if gene:
                out.add(str(gene).upper())
    return out


_STR_FIELDS = (
    "therapeutic_headline",
    "therapeutic_implications",
    "therapeutic_evidence",
    "immunotherapy_eligibility",
    "primary_diagnosis",
    "primary_diagnosis_evidence",
)
_LIST_FIELDS = (
    "actionable_findings",
    "clinical_actions",
    "monitoring_plan",
    "key_findings",
    "recommendations",
    "follow_up",
    "dissenting_opinions",
)
_DICT_LIST_FIELDS = (("differential_diagnoses", ("diagnosis", "evidence")),)


def _iter_text(opinion: Any) -> Iterable[str]:
    """Yield every prose string from a Cancer or rare-disease board opinion.

    Treatment-row drug names are intentionally skipped — those are governed by
    the narrative_scrubber + curate-then-narrate contract, not this module.
    """
    for f in _STR_FIELDS:
        val = getattr(opinion, f, None)
        if isinstance(val, str) and val:
            yield val
    for f in _LIST_FIELDS:
        for item in getattr(opinion, f, None) or []:
            if isinstance(item, str) and item:
                yield item
    for f, keys in _DICT_LIST_FIELDS:
        for item in getattr(opinion, f, None) or []:
            if isinstance(item, dict):
                for k in keys:
                    if isinstance(item.get(k), str) and item[k]:
                        yield item[k]
    # Domain-agent opinions are rendered verbatim (render._render_agent_opinions_section),
    # so an off-briefing gene a domain agent fabricates in its finding/recommendation/
    # concern must be scanned too — not just the Chair synthesis (BOARD-03 / T5-04).
    for agent_op in getattr(opinion, "agent_opinions", None) or []:
        for finding in getattr(agent_op, "findings", None) or []:
            if isinstance(finding, dict):
                for k in ("finding", "evidence"):
                    if isinstance(finding.get(k), str) and finding[k]:
                        yield finding[k]
            elif isinstance(finding, str) and finding:
                yield finding
        for f in ("recommendations", "concerns"):
            for item in getattr(agent_op, f, None) or []:
                if isinstance(item, str) and item:
                    yield item


def find_offbriefing_genes(texts: Iterable[str], briefing_genes: Iterable[str]) -> List[str]:
    """Return known gene symbols appearing in ``texts`` but not in the case.

    The shared detection primitive used by both :func:`annotate_grounding`
    (opinion objects) and the fabrication probe (serialised report JSON). A
    token flags only if it is BOTH a known gene and absent from the briefing.
    """
    known = _known_genes()
    briefing = {str(g).upper() for g in (briefing_genes or [])}
    mentioned: Set[str] = set()
    for text in texts:
        if text:
            mentioned |= set(_GENE_TOKEN_RE.findall(str(text)))
    return sorted((mentioned & known) - briefing - _NON_GENE_TOKENS)


# Tumour type / histology / stage terms. The Chair empirically invents a tumour
# type ("Stage III NSCLC") on cases that carry NO clinical note stating one — it
# infers it from the genes. When the briefing has no clinical note, any such term
# in the narrative is ungrounded and gets flagged (v2.8 Phase 4).
_TUMOR_TYPE_RE = re.compile(
    r"\bstage\s+(?:iv|iii|ii|i|0|[0-4])\b"
    r"|\b(?:NSCLC|SCLC|PDAC|HCC|CRC|TNBC|AML|CML|CLL|ALL|DLBCL|GBM)\b"
    r"|\b[a-z]*(?:adenocarcinoma|carcinoma|sarcoma|melanoma|lymphoma|leukaemia|leukemia|blastoma|glioma|mesothelioma|myeloma)\b",
    re.IGNORECASE,
)


def find_offbriefing_tumor_terms(texts: Iterable[str], has_clinical_note: bool) -> List[str]:
    """Return tumour-type / histology / stage terms in ``texts`` that are NOT
    grounded by a clinical note. Empty when a clinical note is present (the note
    is the grounding source for the tumour type)."""
    if has_clinical_note:
        return []
    found: List[str] = []
    seen: Set[str] = set()
    for text in texts:
        for m in _TUMOR_TYPE_RE.finditer(str(text or "")):
            term = m.group(0).strip()
            key = term.lower()
            if key not in seen:
                seen.add(key)
                found.append(term)
    return found


def annotate_grounding(opinion: Any, report_data: dict) -> dict:
    """Append grounding notes for entities named in ``opinion`` but not grounded
    in the case briefing — (a) known gene symbols absent from the case variant
    set, and (b) tumour type / histology / stage when the briefing carries no
    clinical note stating one. Mutates ``opinion.grounding_flags`` in place.

    Returns ``{"offbriefing_genes": [...], "offbriefing_tumor_terms": [...]}`` for
    caller logging. Annotate-only — never strips or rewrites the narrative (v2.8).
    """
    if not hasattr(opinion, "grounding_flags"):
        return {"offbriefing_genes": [], "offbriefing_tumor_terms": []}

    briefing = _briefing_genes(report_data)
    texts = list(_iter_text(opinion))
    offbriefing: List[str] = find_offbriefing_genes(texts, briefing)
    has_note = bool((report_data.get("clinical_note") or "").strip())
    tumor_terms: List[str] = find_offbriefing_tumor_terms(texts, has_note)

    if offbriefing:
        opinion.grounding_flags.append(
            "⚠ Grounding: gene symbol(s) named in the board narrative but NOT in the "
            "case variant set — verify against the source data before use: " + ", ".join(offbriefing)
        )
        logger.warning(
            "[grounding_scrubber] off-briefing gene(s) in board narrative: %s (case genes: %s)",
            ", ".join(offbriefing),
            ", ".join(sorted(briefing)) or "(none)",
        )
    if tumor_terms:
        opinion.grounding_flags.append(
            "⚠ Grounding: tumour type / stage named with no clinical note stating it — "
            "the briefing carries no tumour-type source, so verify or remove: " + ", ".join(tumor_terms)
        )
        logger.warning(
            "[grounding_scrubber] ungrounded tumour-type/stage term(s) (no clinical note): %s",
            ", ".join(tumor_terms),
        )
    return {"offbriefing_genes": offbriefing, "offbriefing_tumor_terms": tumor_terms}


__all__ = ["annotate_grounding", "find_offbriefing_genes", "find_offbriefing_tumor_terms"]
