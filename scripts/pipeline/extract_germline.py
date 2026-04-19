"""Extract inherited pathogenic variants from a germline VCF.

Uses tabix intersection with a pre-built target BED to avoid reading the
full 4-5 M variant germline file.  Each returned :class:`Variant` carries
``source = "germline_inherited"`` so downstream code (report templates,
``build_variant_records``) can distinguish inherited variants from the
primary (somatic / de novo) input.

Public API
----------
.. autofunction:: extract_inherited_variants
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import List, Optional, Set

from scripts.common.models import Variant

logger = logging.getLogger(__name__)

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
_DEFAULT_TARGET_BED = str(_PROJECT_ROOT / "data" / "germline_targets" / "combined_targets.bed.gz")

# ---------------------------------------------------------------------------
# VCF line parser (lightweight -- avoids cyvcf2 dependency)
# ---------------------------------------------------------------------------

# Minimal VEP CSQ / SnpEff ANN parsing for gene + consequence + HGVS.
_CSQ_GENE_IDX: Optional[int] = None
_CSQ_CONSEQUENCE_IDX: Optional[int] = None
_CSQ_HGVSC_IDX: Optional[int] = None
_CSQ_HGVSP_IDX: Optional[int] = None


def _parse_vcf_line(line: str, fallback_gene: str = "") -> Optional[Variant]:
    """Parse a single VCF data line into a :class:`Variant`.

    Returns None if the line cannot be parsed (header, malformed, etc.).
    """
    if line.startswith("#"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 5:
        return None

    chrom = parts[0]
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    try:
        pos = int(parts[1])
    except (ValueError, IndexError):
        return None
    rsid_raw = parts[2] if len(parts) > 2 else "."
    rsid = rsid_raw if rsid_raw != "." else None
    ref = parts[3]
    alt = parts[4].split(",")[0]  # take first ALT allele

    gene = None
    consequence = None
    hgvsc = None
    hgvsp = None

    # Try to extract gene from INFO field
    info = parts[7] if len(parts) > 7 else ""

    # VEP CSQ
    csq_match = re.search(r"CSQ=([^;\t]+)", info)
    if csq_match:
        csq_val = csq_match.group(1).split(",")[0]  # first transcript
        csq_fields = csq_val.split("|")
        # Standard VEP CSQ field order (may vary, but common defaults)
        if len(csq_fields) > 3:
            gene = csq_fields[3] if csq_fields[3] else gene
        if len(csq_fields) > 1:
            consequence = csq_fields[1] if csq_fields[1] else consequence
        if len(csq_fields) > 10:
            hgvsc = csq_fields[10] if csq_fields[10] else hgvsc
        if len(csq_fields) > 11:
            hgvsp = csq_fields[11] if csq_fields[11] else hgvsp

    # SnpEff ANN fallback
    if not gene:
        ann_match = re.search(r"ANN=([^;\t]+)", info)
        if ann_match:
            ann_val = ann_match.group(1).split(",")[0]
            ann_fields = ann_val.split("|")
            if len(ann_fields) > 3:
                gene = ann_fields[3] if ann_fields[3] else gene
            if len(ann_fields) > 1:
                consequence = ann_fields[1] if ann_fields[1] else consequence

    # GENEINFO fallback (ClinVar-style VCF)
    if not gene:
        gi_match = re.search(r"GENEINFO=([^;:\t]+)", info)
        if gi_match:
            gene = gi_match.group(1)

    # Target BED fallback — when the germline VCF is unannotated (no CSQ /
    # ANN / GENEINFO) the only source of gene identity is the target BED
    # row that selected this position. Without this, inherited ClinVar
    # P/LP variants render as "NONE | Pathogenic" in the report.
    if not gene and fallback_gene:
        gene = fallback_gene

    return Variant(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        gene=gene or None,
        rsid=rsid,
        consequence=consequence,
        hgvsc=hgvsc,
        hgvsp=hgvsp,
        source="germline_inherited",
    )


# ---------------------------------------------------------------------------
# Target BED loader
# ---------------------------------------------------------------------------


def _load_target_regions(target_bed: str) -> list[tuple[str, int, int, str]]:
    """Load target regions from a BED file (plain or bgzipped).

    Returns a list of (chrom, start, end, gene) tuples. The gene column
    (4th BED column) is used as a fallback when the germline VCF's INFO
    field lacks CSQ / ANN / GENEINFO annotation — without this, inherited
    variants showed up in the report as "NONE | Pathogenic" with no way
    to tell which gene they belonged to.
    """
    regions: list[tuple[str, int, int, str]] = []
    path = Path(target_bed)
    if not path.exists():
        return regions

    import gzip

    opener = gzip.open if target_bed.endswith(".gz") else open
    mode = "rt" if target_bed.endswith(".gz") else "r"

    try:
        with opener(target_bed, mode) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                gene = parts[3] if len(parts) >= 4 else ""
                if gene == ".":
                    gene = ""
                regions.append((chrom, start, end, gene))
    except Exception as exc:
        logger.warning("Failed to load target BED %s: %s", target_bed, exc)

    return regions


# ---------------------------------------------------------------------------
# gnomAD AF lookup
# ---------------------------------------------------------------------------


def _get_gnomad_af(variant: Variant) -> Optional[float]:
    """Return gnomAD global AF for *variant*, or None if unavailable."""
    try:
        from scripts.db.query_tabix_gnomad import query_tabix_gnomad

        result = query_tabix_gnomad(variant)
        if result and result.get("gnomad_all") is not None:
            return result["gnomad_all"]
    except Exception:
        pass

    try:
        # ClinVar local may have allele frequency in some builds -- skip
        pass
    except Exception:
        pass

    return None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def extract_inherited_variants(
    germline_vcf: str,
    target_bed: str = _DEFAULT_TARGET_BED,
    max_gnomad_af: float = 0.01,
    primary_variant_ids: Optional[Set[str]] = None,
) -> List[Variant]:
    """Extract known pathogenic / clinically significant variants from a germline VCF.

    Uses tabix intersection with the pre-built target BED to avoid reading
    the full 4-5 M variant germline file.  Filters by gnomAD AF <
    *max_gnomad_af* to exclude common benign variants that happen to overlap
    target positions.

    Parameters
    ----------
    germline_vcf:
        Path to the germline VCF (bgzip + tabix indexed recommended).
    target_bed:
        Path to the combined target BED (bgzip + tabix indexed).
    max_gnomad_af:
        Maximum gnomAD global allele frequency.  Variants at or above this
        threshold are excluded.
    primary_variant_ids:
        Set of variant IDs already present in the primary VCF.  Inherited
        variants that duplicate a primary variant are skipped (primary takes
        precedence).

    Returns
    -------
    List of :class:`Variant` objects, each with ``source="germline_inherited"``.
    Returns an empty list when the target BED or germline VCF is unavailable
    (with a WARNING log).
    """
    if primary_variant_ids is None:
        primary_variant_ids = set()

    # --- Validate inputs ---------------------------------------------------
    if not Path(target_bed).exists():
        logger.warning(
            "Germline target BED not found at %s -- skipping inherited variant extraction. "
            "Run `python scripts/tools/build_germline_target_bed.py` to build it.",
            target_bed,
        )
        return []

    if not Path(germline_vcf).exists():
        logger.warning("Germline VCF not found: %s", germline_vcf)
        return []

    # --- Try tabix-based extraction ----------------------------------------
    variants = _extract_via_tabix(germline_vcf, target_bed, primary_variant_ids)
    if variants is None:
        # tabix not available or failed -- try linear scan with region filter
        variants = _extract_via_linear_scan(germline_vcf, target_bed, primary_variant_ids)

    # --- AF filter ---------------------------------------------------------
    filtered: List[Variant] = []
    for v in variants:
        af = _get_gnomad_af(v)
        if af is not None and af >= max_gnomad_af:
            logger.debug("Excluding common variant %s (gnomAD AF=%.4f)", v.variant_id, af)
            continue
        filtered.append(v)

    # Safety cap: if extraction still returns an unreasonable number of
    # variants after AF filtering, something is wrong with the target BED
    # (e.g. gene-span regions instead of point-level). Cap and warn rather
    # than feeding 100K+ variants into the classification pipeline.
    _MAX_INHERITED = 5000
    if len(filtered) > _MAX_INHERITED:
        logger.warning(
            "Germline extraction returned %d variants (> safety cap %d). "
            "This usually means the target BED contains gene-span regions "
            "instead of point-level ClinVar P/LP positions. Truncating to "
            "the first %d. Rebuild targets with: "
            "python scripts/tools/build_germline_target_bed.py",
            len(filtered),
            _MAX_INHERITED,
            _MAX_INHERITED,
        )
        filtered = filtered[:_MAX_INHERITED]

    logger.info(
        "Germline extraction: %d variants from tabix/scan, %d after AF filter",
        len(variants),
        len(filtered),
    )
    return filtered


def _extract_via_tabix(
    germline_vcf: str,
    target_bed: str,
    primary_ids: Set[str],
) -> Optional[List[Variant]]:
    """Try pysam.TabixFile-based extraction.  Returns None if pysam unavailable."""
    try:
        import pysam
    except ImportError:
        logger.warning("pysam not available -- falling back to linear scan")
        return None

    # Check for tabix index
    tbi = germline_vcf + ".tbi"
    csi = germline_vcf + ".csi"
    if not (Path(tbi).exists() or Path(csi).exists()):
        logger.warning(
            "Germline VCF %s has no tabix index (.tbi/.csi) -- "
            "falling back to linear scan. Index with: tabix -p vcf %s",
            germline_vcf,
            germline_vcf,
        )
        return None

    regions = _load_target_regions(target_bed)
    if not regions:
        logger.warning("Target BED is empty -- no regions to query")
        return []

    variants: List[Variant] = []
    seen_ids: set[str] = set()

    try:
        tbx = pysam.TabixFile(germline_vcf)
    except Exception as exc:
        logger.warning("Cannot open germline VCF via tabix: %s", exc)
        return None

    try:
        contigs = set(tbx.contigs)
        for chrom, start, end, gene in regions:
            # Try both chr-prefixed and bare chromosome names
            query_chrom = None
            if chrom in contigs:
                query_chrom = chrom
            elif chrom.replace("chr", "") in contigs:
                query_chrom = chrom.replace("chr", "")
            elif f"chr{chrom}" in contigs:
                query_chrom = f"chr{chrom}"

            if query_chrom is None:
                continue

            try:
                for line in tbx.fetch(query_chrom, start, end):
                    v = _parse_vcf_line(line, fallback_gene=gene)
                    if v is None:
                        continue
                    if v.variant_id in primary_ids or v.variant_id in seen_ids:
                        continue
                    seen_ids.add(v.variant_id)
                    variants.append(v)
            except ValueError:
                # Region not found in index
                continue
    finally:
        tbx.close()

    return variants


def _extract_via_linear_scan(
    germline_vcf: str,
    target_bed: str,
    primary_ids: Set[str],
) -> List[Variant]:
    """Fallback: linear scan of germline VCF filtered against target regions.

    This is O(n) in the VCF size and should only be used when tabix is
    unavailable.  For a 4-5 M variant VCF this may take tens of seconds.
    """
    import gzip

    regions = _load_target_regions(target_bed)
    if not regions:
        return []

    # Build a lookup: chrom -> sorted list of (start, end, gene)
    region_map: dict[str, list[tuple[int, int, str]]] = {}
    for chrom, start, end, gene in regions:
        bare = chrom.replace("chr", "")
        for c in (chrom, bare, f"chr{bare}"):
            region_map.setdefault(c, []).append((start, end, gene))
    for k in region_map:
        region_map[k].sort()

    def _overlapping_gene(chrom: str, pos: int) -> Optional[str]:
        """Return the target-BED gene for the first region overlapping pos,
        or None if pos falls outside all target regions for chrom."""
        intervals = region_map.get(chrom)
        if not intervals:
            return None
        # BED is 0-based half-open, VCF pos is 1-based
        vcf_0 = pos - 1
        for s, e, g in intervals:
            if s <= vcf_0 < e:
                return g or ""
            if s > vcf_0:
                break
        return None

    logger.info("Linear scan of germline VCF %s (this may be slow for large files)", germline_vcf)
    opener = gzip.open if germline_vcf.endswith(".gz") or germline_vcf.endswith(".bgz") else open
    mode = "rt" if germline_vcf.endswith((".gz", ".bgz")) else "r"

    variants: List[Variant] = []
    seen_ids: set[str] = set()

    try:
        with opener(germline_vcf, mode) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t", 3)
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                try:
                    pos = int(parts[1])
                except ValueError:
                    continue

                target_gene = _overlapping_gene(chrom, pos)
                if target_gene is None:
                    continue

                v = _parse_vcf_line(line, fallback_gene=target_gene)
                if v is None:
                    continue
                if v.variant_id in primary_ids or v.variant_id in seen_ids:
                    continue
                seen_ids.add(v.variant_id)
                variants.append(v)
    except Exception as exc:
        logger.warning("Linear scan of germline VCF failed: %s", exc)

    return variants
