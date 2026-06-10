"""Parse AnnotSV TSV output into StructuralVariant objects."""

import csv
import logging
from typing import Dict, List, Optional

from scripts.common.models import StructuralVariant

logger = logging.getLogger(__name__)


def _safe_int(val: str, default: int = 0) -> int:
    try:
        return int(val) if val else default
    except (ValueError, TypeError):
        return default


def _safe_float(val: str, default: Optional[float] = None) -> Optional[float]:
    try:
        return float(val) if val else default
    except (ValueError, TypeError):
        return default


def parse_annotsv(tsv_path: str) -> List[StructuralVariant]:
    """Parse AnnotSV TSV output into StructuralVariant objects.
    Returns one StructuralVariant per SV (from full rows),
    enriched with gene details (from split rows).
    """
    try:
        # Try utf-8 first, fall back to latin-1
        try:
            f = open(tsv_path, encoding="utf-8")
            content = f.read()
            f.close()
        except UnicodeDecodeError:
            f = open(tsv_path, encoding="latin-1")
            content = f.read()
            f.close()

        if not content.strip():
            return []

        # Detect delimiter (tab or pipe)
        first_line = content.split("\n")[0]
        delimiter = "|" if "|" in first_line and "\t" not in first_line else "\t"

        lines = content.strip().split("\n")
        reader = csv.DictReader(lines, delimiter=delimiter)

        # Header validation: a renamed/reordered or wrong-format file would otherwise
        # parse to garbage/empty silently. Warn loudly when core AnnotSV columns are
        # absent so the reviewer can tell "no SVs" from "parser could not read this file" (T4-06).
        _required = {"AnnotSV_ID", "SV_chrom", "SV_type"}
        _fieldnames = set(reader.fieldnames or [])
        _missing = _required - _fieldnames
        if _missing:
            logger.warning(
                "AnnotSV TSV %s is missing expected column(s): %s — the header may be a "
                "different AnnotSV version/format; SV parsing may be incomplete.",
                tsv_path,
                ", ".join(sorted(_missing)),
            )

        # First pass: collect full rows
        full_rows: Dict[str, Dict] = {}
        split_rows: Dict[str, List[Dict]] = {}

        for row in reader:
            annotsv_id = row.get("AnnotSV_ID", "")
            mode = row.get("Annotation_mode", "full")

            if mode == "full":
                full_rows[annotsv_id] = row
                if annotsv_id not in split_rows:
                    split_rows[annotsv_id] = []
            elif mode == "split":
                if annotsv_id not in split_rows:
                    split_rows[annotsv_id] = []
                split_rows[annotsv_id].append(row)

        # Build StructuralVariant objects
        results = []
        for annotsv_id, row in full_rows.items():
            gene_details = []
            for sr in split_rows.get(annotsv_id, []):
                gene_details.append(
                    {
                        "gene": sr.get("Gene_name", ""),
                        "transcript": sr.get("Tx", ""),
                        "cds_percent": _safe_float(sr.get("Overlapped_CDS_percent"), 0.0),
                        "frameshift": sr.get("Frameshift", ""),
                        "location": sr.get("Location", ""),
                        "exon_count": _safe_int(sr.get("Exon_count")),
                        "hi": _safe_int(sr.get("HI")),
                        "ts": _safe_int(sr.get("TS")),
                        "pli": _safe_float(sr.get("GnomAD_pLI"), 0.0),
                        "omim_morbid": (sr.get("OMIM_morbid") or "").lower() == "yes",
                    }
                )

            # Distinguish "AnnotSV scored this SV" from "AnnotSV could not score it"
            # (blank / NA / out-of-range ACMG_class — BND, complex, annotation gap).
            # Unscored SVs keep the display default 3 but are flagged so they reach the
            # reviewer as 'unscored — manual review' rather than silently as a VUS (T2-06).
            raw_acmg = (row.get("ACMG_class") or "").strip()
            try:
                acmg_class = int(raw_acmg)
                acmg_scored = 1 <= acmg_class <= 5
                if not acmg_scored:
                    acmg_class = 3
            except (ValueError, TypeError):
                acmg_class = 3
                acmg_scored = False

            sv = StructuralVariant(
                annotsv_id=annotsv_id,
                chrom=row.get("SV_chrom", ""),
                start=_safe_int(row.get("SV_start")),
                end=_safe_int(row.get("SV_end")),
                length=_safe_int(row.get("SV_length")),
                sv_type=row.get("SV_type", ""),
                sample_id=row.get("Samples_ID", ""),
                acmg_class=acmg_class,
                acmg_scored=acmg_scored,
                # AnnotSV v3.x renamed AnnotSV_ranking -> AnnotSV_ranking_score; accept either (T4-06).
                ranking_score=_safe_float(row.get("AnnotSV_ranking") or row.get("AnnotSV_ranking_score"), 0.0) or 0.0,
                cytoband=row.get("CytoBand", ""),
                gene_name=row.get("Gene_name", ""),
                gene_count=_safe_int(row.get("Gene_count")),
                gene_details=gene_details,
                p_gain_phen=row.get("P_gain_phen", ""),
                p_gain_hpo=row.get("P_gain_hpo", ""),
                p_gain_source=row.get("P_gain_source", ""),
                p_loss_phen=row.get("P_loss_phen", ""),
                p_loss_hpo=row.get("P_loss_hpo", ""),
                p_loss_source=row.get("P_loss_source", ""),
                b_gain_af_max=_safe_float(row.get("B_gain_AFmax")),
                b_loss_af_max=_safe_float(row.get("B_loss_AFmax")),
                omim_morbid=(row.get("OMIM_morbid") or "").lower() == "yes",
            )
            results.append(sv)

        # Sort: pathogenic first (class 5,4,3,2,1)
        results.sort(key=lambda s: (-s.acmg_class, s.gene_name))
        logger.info(f"Parsed {len(results)} structural variants from {tsv_path}")
        return results

    except FileNotFoundError:
        logger.error(f"AnnotSV file not found: {tsv_path}")
        return []
    except Exception as e:
        logger.error(f"Error parsing AnnotSV file {tsv_path}: {e}")
        return []
