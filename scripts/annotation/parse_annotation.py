"""Parse VEP CSQ or SnpEff ANN fields from pre-annotated VCFs."""

import logging
import re
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def parse_csq_header(header_line: str) -> List[str]:
    """Extract CSQ field names from VCF ##INFO=<ID=CSQ...> header.
    Returns list of field names like ['Allele', 'Consequence', 'IMPACT', 'SYMBOL', ...]
    """
    if "Format:" in header_line:
        fmt = header_line.split("Format:")[1].strip().rstrip('">')
        return [f.strip() for f in fmt.split("|")]
    return []


def parse_ann_header(header_line: str) -> List[str]:
    """Extract ANN field names from VCF ##INFO=<ID=ANN...> header."""
    # SnpEff ANN has a standard field order
    return [
        "Allele",
        "Annotation",
        "Annotation_Impact",
        "Gene_Name",
        "Gene_ID",
        "Feature_Type",
        "Feature_ID",
        "Transcript_BioType",
        "Rank",
        "HGVS.c",
        "HGVS.p",
        "cDNA.pos/cDNA.length",
        "CDS.pos/CDS.length",
        "AA.pos/AA.length",
        "Distance",
        "ERRORS/WARNINGS/INFO",
    ]


def _pick_best_transcript(entries: List[Dict[str, str]], gene: Optional[str] = None) -> Optional[Dict[str, str]]:
    """Pick the best annotation entry — prefer MANE/canonical, matching gene, highest impact."""
    if not entries:
        return None

    # Filter by gene if specified
    if gene:
        gene_matches = [e for e in entries if e.get("gene", "").upper() == gene.upper()]
        if gene_matches:
            entries = gene_matches

    # Prefer entries with MANE_SELECT or CANONICAL=YES
    for e in entries:
        if e.get("mane_select") or e.get("canonical") == "YES":
            return e

    # Sort by impact: HIGH > MODERATE > LOW > MODIFIER
    impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
    entries = sorted(entries, key=lambda e: impact_order.get(e.get("impact", "MODIFIER"), 3))

    return entries[0]


def parse_csq_value(
    csq_string: str, fields: List[str], gene: Optional[str] = None, alt: Optional[str] = None
) -> Optional[Dict[str, str]]:
    """Parse a CSQ INFO value into an annotation dict. Picks the best transcript.

    ``alt`` is the called ALT allele(s) from the VCF (comma-separated for
    multi-allelic lines). When supplied, entries are restricted to those whose
    VEP ``Allele`` field matches a called allele, so a multi-allelic locus is
    never annotated with a *different* allele's consequence (v2.7 review
    CROS-03). The filter is conservative: if no entry matches (e.g. an indel
    minimal-representation mismatch) all entries are kept, so annotation is
    never silently lost.
    """
    entries = []
    for entry_str in csq_string.split(","):
        values = entry_str.split("|")
        entry = {}
        for i, val in enumerate(values):
            if i < len(fields):
                entry[fields[i].lower()] = val
        # Map to standard keys
        # Extract rsID from Existing_variation (e.g., "rs12345" or "rs12345&COSM67890")
        existing_var = entry.get("existing_variation", "")
        rsid = ""
        if existing_var:
            for var_id in existing_var.split("&"):
                if var_id.startswith("rs"):
                    rsid = var_id
                    break

        mapped = {
            "allele": entry.get("allele", ""),
            "gene": entry.get("symbol", ""),
            "consequence": entry.get("consequence", ""),
            "impact": entry.get("impact", ""),
            "transcript": entry.get("feature", ""),
            "hgvsc": entry.get("hgvsc", ""),
            "hgvsp": entry.get("hgvsp", ""),
            "sift": entry.get("sift", ""),
            "polyphen": entry.get("polyphen", ""),
            "canonical": entry.get("canonical", ""),
            "mane_select": entry.get("mane_select", ""),
            "rsid": rsid,
            "clin_sig": entry.get("clin_sig", ""),
            "gmaf": entry.get("gmaf", ""),
            # In silico prediction scores (from VEP dbNSFP plugin)
            "revel_score": entry.get("revel_score", ""),
            "cadd_phred": entry.get("cadd_phred", ""),
            "am_class": entry.get("am_class", ""),
            "am_pathogenicity": entry.get("am_pathogenicity", ""),
            "spliceai_pred_ds_ag": entry.get("spliceai_pred_ds_ag", ""),
            "spliceai_pred_ds_al": entry.get("spliceai_pred_ds_al", ""),
            "spliceai_pred_ds_dg": entry.get("spliceai_pred_ds_dg", ""),
            "spliceai_pred_ds_dl": entry.get("spliceai_pred_ds_dl", ""),
            "domains": entry.get("domains", ""),
        }
        entries.append(mapped)

    # Multi-allelic safety: scope to the called ALT allele(s) when known so a
    # multi-allelic locus is never annotated with another allele's consequence.
    if alt:
        called = {a.strip() for a in str(alt).split(",") if a.strip()}
        matched = [e for e in entries if e.get("allele", "") in called]
        if matched:
            entries = matched

    return _pick_best_transcript(entries, gene)


def parse_ann_value(
    ann_string: str, fields: List[str], gene: Optional[str] = None, alt: Optional[str] = None
) -> Optional[Dict[str, str]]:
    """Parse a SnpEff ANN INFO value into annotation dict.

    When ``alt`` is the single called ALT allele, ANN entries are scoped to it so a
    multi-allelic locus is never annotated with another allele's (higher-impact)
    consequence — mirrors ``parse_csq_value`` (INTK-05).
    """
    entries = []
    for entry_str in ann_string.split(","):
        values = entry_str.split("|")
        entry = {}
        for i, val in enumerate(values):
            if i < len(fields):
                entry[fields[i].lower()] = val.strip()
        mapped = {
            "allele": entry.get("allele", ""),
            "gene": entry.get("gene_name", ""),
            "consequence": entry.get("annotation", ""),
            "impact": entry.get("annotation_impact", ""),
            "transcript": entry.get("feature_id", ""),
            "hgvsc": entry.get("hgvs.c", ""),
            "hgvsp": entry.get("hgvs.p", ""),
            "sift": "",
            "polyphen": "",
        }
        entries.append(mapped)

    if alt:
        called = {a.strip() for a in str(alt).split(",") if a.strip()}
        matched = [e for e in entries if e.get("allele", "") in called]
        if matched:
            entries = matched

    return _pick_best_transcript(entries, gene)


def format_consequence(consequence: str) -> str:
    """Convert VEP consequence term to human-readable format.
    e.g., 'missense_variant' -> 'Missense', 'frameshift_variant' -> 'Frameshift'
    """
    mapping = {
        "missense_variant": "Missense",
        "nonsense": "Nonsense",
        "stop_gained": "Nonsense / Stop gain",
        "frameshift_variant": "Frameshift",
        "splice_donor_variant": "Splice donor",
        "splice_acceptor_variant": "Splice acceptor",
        "splice_region_variant": "Splice region",
        "inframe_deletion": "In-frame deletion",
        "inframe_insertion": "In-frame insertion",
        "synonymous_variant": "Synonymous",
        "start_lost": "Start loss",
        "stop_lost": "Stop loss",
        "intron_variant": "Intronic",
        "5_prime_UTR_variant": "5' UTR",
        "3_prime_UTR_variant": "3' UTR",
        "intergenic_variant": "Intergenic",
    }
    if not consequence:
        return ""
    # Handle multiple consequences (e.g., "missense_variant&splice_region_variant")
    parts = consequence.split("&")
    readable = []
    for part in parts:
        readable.append(mapping.get(part, part.replace("_", " ").title()))
    return " / ".join(readable[:2])  # Show max 2


# ---------------------------------------------------------------------------
# Canonical consequence resolution — single source of truth
# ---------------------------------------------------------------------------
#
# The real pipeline stores ``variant.consequence`` as a BIKO short label
# produced by :func:`format_consequence` (e.g. ``"Missense"``,
# ``"Splice region"``, or the compound ``"Missense / Splice region"``), while
# tests and some upstream callers pass raw VEP SO terms (``"missense_variant"``,
# ``"missense_variant&splice_region_variant"``). Both forms must resolve to the
# same canonical VEP SO term or downstream membership checks become silent dead
# code on production data (see v2.7 review — BOAR-01 / CROS-03 / CROS-11).
#
# This map and the two helpers below are the ONE place that knowledge lives.
# ``variant_selector``, ``evidence_collector``, and ``tmb`` all delegate here
# instead of carrying their own (previously divergent) copies.

# Inverse of ``format_consequence`` — lowercased BIKO label → canonical SO term.
# Includes both the canonical 2-word forms emitted by ``format_consequence`` and
# the title-cased / legacy variants seen in older data so the round-trip is
# bullet-proof regardless of which label form a caller supplies.
_FORMATTED_TO_SO: Dict[str, str] = {
    "missense": "missense_variant",
    "nonsense": "stop_gained",
    "nonsense / stop gain": "stop_gained",
    "stop gain": "stop_gained",
    "frameshift": "frameshift_variant",
    "splice donor": "splice_donor_variant",
    "splice acceptor": "splice_acceptor_variant",
    "splice region": "splice_region_variant",
    "splice region variant": "splice_region_variant",
    "in-frame deletion": "inframe_deletion",
    "inframe deletion": "inframe_deletion",
    "in-frame insertion": "inframe_insertion",
    "inframe insertion": "inframe_insertion",
    "synonymous": "synonymous_variant",
    "start loss": "start_lost",
    "stop loss": "stop_lost",
    "intronic": "intron_variant",
    "5' utr": "5_prime_utr_variant",
    "3' utr": "3_prime_utr_variant",
    "intergenic": "intergenic_variant",
}

# A consequence label may be compound — joined by " / " (``format_consequence``)
# or "&" (raw VEP). Split on either delimiter.
_CONSEQUENCE_SPLIT = re.compile(r"\s*/\s*|&")


def canonical_consequence_parts(label: str) -> List[str]:
    """Return every consequence component of ``label`` as canonical SO terms.

    Handles raw VEP SO terms, BIKO short labels, and compound labels joined by
    " / " or "&". Components that are not in the inverse map pass through
    unchanged (they are already SO terms, or genuinely unknown). Order is
    preserved, so element 0 is the primary (most-severe) consequence.
    """
    if not label:
        return []
    parts = []
    for raw in _CONSEQUENCE_SPLIT.split(label.strip().lower()):
        part = raw.strip()
        if part:
            parts.append(_FORMATTED_TO_SO.get(part, part))
    return parts


def canonical_consequence(label: str) -> str:
    """Return the primary (most-severe) consequence of ``label`` as an SO term.

    Empty input → ``""``. See :func:`canonical_consequence_parts`.
    """
    parts = canonical_consequence_parts(label)
    return parts[0] if parts else ""
