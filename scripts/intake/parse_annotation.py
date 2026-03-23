"""Parse VEP CSQ or SnpEff ANN fields from pre-annotated VCFs."""
import logging
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
        "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID",
        "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank",
        "HGVS.c", "HGVS.p", "cDNA.pos/cDNA.length", "CDS.pos/CDS.length",
        "AA.pos/AA.length", "Distance", "ERRORS/WARNINGS/INFO"
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


def parse_csq_value(csq_string: str, fields: List[str], gene: Optional[str] = None) -> Optional[Dict[str, str]]:
    """Parse a CSQ INFO value into annotation dict. Picks best transcript."""
    entries = []
    for entry_str in csq_string.split(","):
        values = entry_str.split("|")
        entry = {}
        for i, val in enumerate(values):
            if i < len(fields):
                entry[fields[i].lower()] = val
        # Map to standard keys
        mapped = {
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
        }
        entries.append(mapped)

    return _pick_best_transcript(entries, gene)


def parse_ann_value(ann_string: str, fields: List[str], gene: Optional[str] = None) -> Optional[Dict[str, str]]:
    """Parse a SnpEff ANN INFO value into annotation dict."""
    entries = []
    for entry_str in ann_string.split(","):
        values = entry_str.split("|")
        entry = {}
        for i, val in enumerate(values):
            if i < len(fields):
                entry[fields[i].lower()] = val.strip()
        mapped = {
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
