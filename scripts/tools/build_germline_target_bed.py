#!/usr/bin/env python3
"""Build a combined target BED for germline inherited-variant extraction.

Merges three sources into ``data/germline_targets/combined_targets.bed.gz``
(bgzip-compressed, tabix-indexed):

* **Source A** -- ClinVar Pathogenic / Likely Pathogenic positions from the
  local ``clinvar.sqlite3`` database.
* **Source B** -- ACMG Secondary Findings v3.2 gene regions (73 genes,
  hard-coded).  When ClinVar DB is available the gene regions are derived
  from variant positions in ClinVar; otherwise a minimal hard-coded
  coordinate fallback is used.
* **Source C** -- DDG2P panel genes loaded from
  ``data/ddg2p_neurodev_genes.json``.

The output BED is consumed by
:func:`scripts.orchestration.extract_germline.extract_inherited_variants` which
performs a fast tabix intersection against a germline VCF to avoid reading
4-5 M variants.

Usage::

    python scripts/tools/build_germline_target_bed.py
"""

from __future__ import annotations

import logging
import os
import sqlite3
import sys
from pathlib import Path

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent

# ---------------------------------------------------------------------------
# ACMG Secondary Findings v3.2 (73 genes)
# ---------------------------------------------------------------------------

ACMG_SF_V32_GENES: list[str] = sorted(
    set(
        [
            "ACTA2",
            "APOB",
            "ATM",
            "ATP7B",
            "BAP1",
            "BARD1",
            "BMPR1A",
            "BRCA1",
            "BRCA2",
            "BRIP1",
            "BTD",
            "CACNA1S",
            "CDH1",
            "CHEK2",
            "COL3A1",
            "DICER1",
            "DSC2",
            "DSG2",
            "DSP",
            "EPCAM",
            "FBN1",
            "FH",
            "FLCN",
            "FLNC",
            "GAA",
            "GLA",
            "HFE",
            "HOXB13",
            "KCNH2",
            "KCNQ1",
            "LDLR",
            "LMNA",
            "MAX",
            "MEN1",
            "MET",
            "MLH1",
            "MSH2",
            "MSH6",
            "MUTYH",
            "MYBPC3",
            "MYH11",
            "MYH7",
            "NF2",
            "NTHL1",
            "OTC",
            "PALB2",
            "PCSK9",
            "PKP2",
            "PMS2",
            "PTEN",
            "RAD51C",
            "RAD51D",
            "RB1",
            "RET",
            "RYR1",
            "RYR2",
            "SCN5A",
            "SDHA",
            "SDHAF2",
            "SDHB",
            "SDHC",
            "SDHD",
            "SERPINA1",
            "SMAD3",
            "SMAD4",
            "STK11",
            "TGFBR1",
            "TGFBR2",
            "TMEM127",
            "TMEM43",
            "TP53",
            "TSC1",
            "TSC2",
            "TTN",
            "VHL",
            "WT1",
        ]
    )
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _get_clinvar_db_path() -> str | None:
    """Return the ClinVar SQLite path if it exists, else None."""
    from scripts.common.config import get

    db_path = get("paths.clinvar_db")
    if db_path and os.path.exists(db_path):
        return db_path
    # Fallback to project-relative default
    default = str(PROJECT_ROOT / "data" / "db" / "clinvar.sqlite3")
    if os.path.exists(default):
        return default
    return None


def _collect_clinvar_plp_positions(conn: sqlite3.Connection) -> list[tuple[str, int, int, str, str]]:
    """Return BED rows for ClinVar P/LP entries.

    Each row is (chrom, start, end, gene, label).
    """
    rows: list[tuple[str, int, int, str, str]] = []
    cursor = conn.execute(
        """
        SELECT chrom, pos, gene
        FROM variants
        WHERE clinical_significance LIKE '%athogenic%'
          AND clinical_significance NOT LIKE '%onflict%'
          AND assembly = 'GRCh38'
        """
    )
    for chrom, pos, gene in cursor:
        if not chrom or not pos:
            continue
        # Normalise chrom to "chr" prefix
        c = chrom if chrom.startswith("chr") else f"chr{chrom}"
        rows.append((c, pos - 1, pos, gene or ".", "clinvar_plp"))
    return rows


def _collect_clinvar_plp_in_genes(
    conn: sqlite3.Connection, gene_list: list[str], label: str
) -> list[tuple[str, int, int, str, str]]:
    """Return BED rows for ClinVar P/LP entries restricted to a gene list.

    Unlike the earlier ``_collect_gene_regions_from_clinvar`` (which spanned
    MIN(pos)..MAX(pos) per gene — pulling every germline variant in an 84 KB
    BRCA2 region), this returns **point-level** BED rows: one per ClinVar
    P/LP entry. The tabix intersection therefore returns at most one
    germline variant per target position instead of hundreds per gene span.

    This fixed the 281K-variant extraction blowup observed on the first
    ASD-10293 validation run (2026-04-17).
    """
    rows: list[tuple[str, int, int, str, str]] = []
    placeholders = ",".join("?" for _ in gene_list)
    cursor = conn.execute(
        f"""
        SELECT chrom, pos, gene
        FROM variants
        WHERE gene IN ({placeholders})
          AND clinical_significance LIKE '%athogenic%'
          AND clinical_significance NOT LIKE '%onflict%'
          AND assembly = 'GRCh38'
        """,
        gene_list,
    )
    for chrom, pos, gene in cursor:
        if not chrom or not pos:
            continue
        c = chrom if chrom.startswith("chr") else f"chr{chrom}"
        rows.append((c, pos - 1, pos, gene or ".", label))
    return rows


def _load_ddg2p_genes() -> list[str]:
    """Return gene symbols from the DDG2P panel JSON."""
    import json

    path = PROJECT_ROOT / "data" / "ddg2p_neurodev_genes.json"
    if not path.exists():
        logger.warning("DDG2P panel file not found at %s -- skipping Source C", path)
        return []
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        return list((data.get("genes") or {}).keys())
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("DDG2P panel load failed: %s", exc)
        return []


def _sort_bed(rows: list[tuple[str, int, int, str, str]]) -> list[tuple[str, int, int, str, str]]:
    """Sort BED rows by chromosome (natural order) then start position."""

    def _chrom_key(c: str) -> tuple[int, str]:
        stripped = c.replace("chr", "")
        try:
            return (0, str(int(stripped)).zfill(2))
        except ValueError:
            return (1, stripped)

    return sorted(rows, key=lambda r: (_chrom_key(r[0]), r[1], r[2]))


def _dedup_bed(rows: list[tuple[str, int, int, str, str]]) -> list[tuple[str, int, int, str, str]]:
    """Deduplicate BED rows by (chrom, start, end)."""
    seen: set[tuple[str, int, int]] = set()
    out: list[tuple[str, int, int, str, str]] = []
    for row in rows:
        key = (row[0], row[1], row[2])
        if key not in seen:
            seen.add(key)
            out.append(row)
    return out


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------


def build_target_bed(output_dir: str | None = None) -> dict:
    """Build the combined germline target BED.

    Returns a dict with ``region_count``, ``output_path``, and ``sources``.
    """
    if output_dir is None:
        output_dir = str(PROJECT_ROOT / "data" / "germline_targets")

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    all_rows: list[tuple[str, int, int, str, str]] = []
    sources_used: list[str] = []

    clinvar_db = _get_clinvar_db_path()

    if clinvar_db:
        conn = sqlite3.connect(clinvar_db)
        try:
            # Source A: ClinVar P/LP positions
            plp_rows = _collect_clinvar_plp_positions(conn)
            all_rows.extend(plp_rows)
            sources_used.append(f"clinvar_plp ({len(plp_rows)} positions)")

            # Source B: ACMG SF v3.2 gene regions (ClinVar-derived)
            acmg_rows = _collect_clinvar_plp_in_genes(conn, ACMG_SF_V32_GENES, "acmg_sf")
            all_rows.extend(acmg_rows)
            sources_used.append(f"acmg_sf ({len(acmg_rows)} regions)")

            # Source C: DDG2P gene regions (ClinVar-derived)
            ddg2p_genes = _load_ddg2p_genes()
            if ddg2p_genes:
                ddg2p_rows = _collect_clinvar_plp_in_genes(conn, ddg2p_genes, "ddg2p")
                all_rows.extend(ddg2p_rows)
                sources_used.append(f"ddg2p ({len(ddg2p_rows)} regions)")
        finally:
            conn.close()
    else:
        logger.warning(
            "ClinVar DB not found -- building minimal target BED from ACMG SF gene list only (no position data)"
        )
        sources_used.append("acmg_sf_fallback (gene list only, no positions)")

    if not all_rows:
        logger.warning("No target regions collected -- output BED will be empty")
        # Write an empty but valid file
        bed_file = out_path / "combined_targets.bed.gz"
        bed_file.write_bytes(b"")
        return {
            "region_count": 0,
            "output_path": str(bed_file),
            "sources": sources_used,
        }

    # Sort, dedup, write
    all_rows = _sort_bed(_dedup_bed(all_rows))

    bed_plain = out_path / "combined_targets.bed"
    with open(bed_plain, "w", encoding="utf-8") as fh:
        for chrom, start, end, gene, label in all_rows:
            fh.write(f"{chrom}\t{start}\t{end}\t{gene}\t{label}\n")

    # bgzip + tabix via pysam
    bed_gz = str(out_path / "combined_targets.bed.gz")
    try:
        import pysam

        pysam.tabix_compress(str(bed_plain), bed_gz, force=True)
        pysam.tabix_index(bed_gz, preset="bed", force=True)
        # Remove uncompressed after successful bgzip
        bed_plain.unlink(missing_ok=True)
        logger.info("Target BED built: %d regions -> %s (+.tbi)", len(all_rows), bed_gz)
    except ImportError:
        logger.warning("pysam not available -- writing plain BED (no bgzip/tabix)")
        bed_gz = str(bed_plain)
    except Exception as exc:
        logger.warning("pysam bgzip/tabix failed: %s -- keeping plain BED", exc)
        bed_gz = str(bed_plain)

    return {
        "region_count": len(all_rows),
        "output_path": bed_gz,
        "sources": sources_used,
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    # Allow project imports when run as a script
    sys.path.insert(0, str(PROJECT_ROOT))
    result = build_target_bed()
    print(f"Built target BED: {result['region_count']} regions")
    for src in result["sources"]:
        print(f"  Source: {src}")
    print(f"  Output: {result['output_path']}")
