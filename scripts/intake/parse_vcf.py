# scripts/intake/parse_vcf.py
import gzip
import logging
import os
import re
from typing import List, Optional, Tuple

from scripts.common.models import Variant
from scripts.intake.parse_annotation import (
    format_consequence,
    parse_ann_header,
    parse_ann_value,
    parse_csq_header,
    parse_csq_value,
)

logger = logging.getLogger(__name__)

# Default DGQ threshold from triodenovo (see the header comment of ASD trio
# VCFs: `--minDQ 7`). Variants with lower DGQ are treated as noise and not
# labelled de novo even if the parental genotypes look homozygous reference.
_DEFAULT_DGQ_THRESHOLD = 7.0


def _open_vcf(vcf_path: str):
    """Open a VCF file as text, transparently handling gzip compression.

    Detects compressed input by file-extension (`.gz`, `.bgz`) and falls back
    to plain-text open otherwise. Returns a text-mode file handle so callers
    can iterate with `for line in handle:` unchanged.

    Real-pipeline VCFs (cohort / trio / batch) are almost always shipped as
    bgzipped `.vcf.gz`; without this helper the parser raises
    `UnicodeDecodeError: 0x8b` on the first gzip magic byte.
    """
    lowered = vcf_path.lower()
    if lowered.endswith(".gz") or lowered.endswith(".bgz"):
        return gzip.open(vcf_path, "rt", encoding="utf-8")
    return open(vcf_path, encoding="utf-8")


def _detect_trio_proband(vcf_path: str, sample_ids: List[str]) -> Tuple[Optional[int], List[int]]:
    """Return (proband_index, parent_indices) for a trio VCF, or (None, []).

    Heuristic proband detection (v1 — PED file support is v2):
      1. Trio detection requires exactly 3 samples in the #CHROM header.
      2. For each sample ID, extract any ≥3-digit numeric suffix and check
         whether that digit string appears in the VCF filename. If one and
         only one sample matches, it's the proband — the user convention
         is that trio VCFs are named after the proband sample
         (e.g. IBS-ASD-9703-...triodenovo.vcf includes the proband 9703
         plus parents 9701 and 9702).
      3. If zero or multiple IDs match, fall back to "last sample is
         proband" with a visible warning so the reviewer can decide to
         rename the file or supply a PED later.

    Only applies to exactly-3-sample VCFs. Singleton / quartet / cohort
    VCFs return (None, []) so the caller skips trio de novo detection
    entirely.
    """
    if len(sample_ids) != 3:
        return None, []

    fname = os.path.basename(vcf_path).lower()
    matches: List[int] = []
    for idx, sid in enumerate(sample_ids):
        sid_lower = sid.lower()
        if sid_lower in fname:
            matches.append(idx)
            continue
        for num in re.findall(r"\d{3,}", sid_lower):
            if num in fname:
                matches.append(idx)
                break

    if len(matches) == 1:
        proband_idx = matches[0]
        parent_idxs = [i for i in range(3) if i != proband_idx]
        logger.info(
            "Trio VCF detected: proband=%s (matched to filename), parents=%s",
            sample_ids[proband_idx],
            [sample_ids[i] for i in parent_idxs],
        )
        return proband_idx, parent_idxs

    logger.warning(
        "Trio VCF with %d samples but no unique filename-to-sample match. "
        "Falling back to 'last sample is proband' convention "
        "(proband=%s). Rename the VCF to include the proband ID, or add "
        "PED-file support (future ticket), for an explicit mapping.",
        len(sample_ids),
        sample_ids[-1],
    )
    return len(sample_ids) - 1, [0, 1]


def _gt_value(sample_col: str, gt_idx: int) -> str:
    """Extract GT sub-field from a VCF sample column by format index."""
    parts = sample_col.split(":")
    if 0 <= gt_idx < len(parts):
        return parts[gt_idx]
    return ""


def _normalize_gt(gt: str, ref: str, alt: str) -> str:
    """Normalize a GT value to the VCF-standard integer form (0/0, 0/1, ...).

    Most variant callers emit GT as integer indices into REF/ALT. However
    ``polymutt`` rewrites GT to the actual nucleotide base (e.g. ``G/G``,
    ``A/G``) — this is present in our real ASD trio fixtures' VCF header:
    ``##Note=VCF file modified by polymutt. Updated fields include: ... GT ...``.

    Without normalisation the trio de novo detector silently fails on every
    polymutt-produced row (``G/G`` never matches ``"0/0"``). We normalise
    here by mapping each allele token back to its REF/ALT index, so the
    downstream helpers can assume standard integer GT.

    Unknown tokens (not ``.``, not REF, not in the ALT comma-list) are
    returned unchanged so we fail safe rather than silently reclassifying
    a mystery allele as reference.
    """
    if not gt or gt in (".", "./.", ".|."):
        return gt
    for delim in ("/", "|"):
        if delim in gt:
            alleles = gt.split(delim)
            # Already in integer form → no work to do.
            if all(a.isdigit() or a == "." for a in alleles):
                return gt
            alt_list = alt.split(",") if alt else []
            norm: List[str] = []
            for a in alleles:
                if a == ".":
                    norm.append(".")
                elif a == ref:
                    norm.append("0")
                elif a in alt_list:
                    norm.append(str(alt_list.index(a) + 1))
                else:
                    norm.append(a)
            return delim.join(norm)
    # Haploid / single-allele form.
    if gt.isdigit() or gt == ".":
        return gt
    if gt == ref:
        return "0"
    if alt and gt in alt.split(","):
        return str(alt.split(",").index(gt) + 1)
    return gt


def _is_ref_homozygous(gt: str) -> bool:
    """Return True if GT looks like 0/0, 0|0, or the missing-data case.

    Missing-data ("./.", ".|.") is conservatively treated as "not evidence
    of a de novo event" — we want at most a PM6 call, not a false positive
    driven by a missing parental call.
    """
    if not gt or gt in (".", "./.", ".|."):
        return False
    for delim in ("/", "|"):
        if delim in gt:
            alleles = gt.split(delim)
            return all(a == "0" for a in alleles)
    return gt == "0"


def _has_alt_allele(gt: str) -> bool:
    """Return True if GT carries at least one non-reference allele."""
    if not gt or gt in (".", "./.", ".|."):
        return False
    for delim in ("/", "|"):
        if delim in gt:
            alleles = gt.split(delim)
            return any(a not in ("0", ".", "") for a in alleles)
    return gt not in ("0", ".", "")


def _detect_trio_denovo_from_genotype(
    format_keys: List[str],
    proband_sample: str,
    p1_sample: str,
    p2_sample: str,
    ref: str = "",
    alt: str = "",
    dgq_threshold: float = _DEFAULT_DGQ_THRESHOLD,
) -> Optional[str]:
    """Return ``"de_novo"`` when the trio GT pattern implies assumed de novo.

    Per ACMG/AMP 2015 + ClinGen SVI 2018 (PMID 30311383), trio-callset
    genotype concordance alone is **not** confirmed de novo — sample
    swap / non-paternity / germline mosaicism cannot be ruled out without
    an explicit parental identity test. This helper therefore only
    returns the "assumed de novo" label (PM6 downstream), never
    "confirmed". The PS2 upgrade path remains gated behind the
    `CONFIRMED_DN` / `PS2` INFO flag or the future
    `rare_disease.denovo.paternity_confirmed` config opt-in.

    DGQ quality gate: triodenovo's default threshold is 7 (see the
    `--minDQ 7` in its command line in our fixture VCFs). If the proband
    sample carries a DGQ sub-field below threshold, the variant is
    treated as noise and the label is NOT set — better to miss a low-
    quality candidate than flood the rare-disease selector with trio
    caller false positives.
    """
    if "GT" not in format_keys:
        return None
    gt_idx = format_keys.index("GT")
    # Normalise to standard integer form so polymutt-modified VCFs
    # (which emit GT as nucleotides: G/G, A/G, ...) go through the
    # same index-based homozygous / alt checks as vanilla GATK output.
    p1_gt = _normalize_gt(_gt_value(p1_sample, gt_idx), ref, alt)
    p2_gt = _normalize_gt(_gt_value(p2_sample, gt_idx), ref, alt)
    proband_gt = _normalize_gt(_gt_value(proband_sample, gt_idx), ref, alt)

    if not (_is_ref_homozygous(p1_gt) and _is_ref_homozygous(p2_gt)):
        return None
    if not _has_alt_allele(proband_gt):
        return None

    # Optional DGQ quality gate. Absence of DGQ is OK — fall through.
    if "DGQ" in format_keys:
        dgq_idx = format_keys.index("DGQ")
        proband_parts = proband_sample.split(":")
        if 0 <= dgq_idx < len(proband_parts):
            try:
                dgq = float(proband_parts[dgq_idx])
            except ValueError:
                dgq = None
            if dgq is not None and dgq < dgq_threshold:
                return None

    return "de_novo"


def parse_vcf(vcf_path: str, ped_path: Optional[str] = None) -> List[Variant]:
    """Parse a VCF file into a list of Variant objects.
    Uses simple text parsing to avoid cyvcf2 dependency issues in tests.
    Supports VEP CSQ and SnpEff ANN annotation fields if present.
    Accepts plain-text `.vcf` as well as gzip/bgzip `.vcf.gz` / `.vcf.bgz`.

    Trio-aware de novo detection (v2.3-T9): if the #CHROM header exposes
    sample columns the parser identifies proband + parents and, per
    variant, checks whether the parent genotypes are homozygous reference
    and the proband carries an alt allele. Such variants are labelled
    ``inheritance="de_novo"`` so downstream collect_denovo_evidence
    fires PM6 even on raw triodenovo output that lacks INFO-level DN /
    DENOVO flags. PS2 is never inferred from genotype concordance alone
    — see _detect_trio_denovo_from_genotype for the ACMG rationale.

    Proband resolution precedence (v2.4 Quick Win C):
      1. If ``ped_path`` is provided → **strict** PED resolution via
         :func:`scripts.intake.parse_ped.resolve_trio`. Supports N≥3 sample
         VCFs (quartet, cohort). Fails loud with ``ValueError`` when the
         PED cannot resolve a trio against the VCF samples — no silent
         fallback, because mis-assigning the proband silently propagates
         into de novo misclassification.
      2. Otherwise → legacy filename heuristic
         (:func:`_detect_trio_proband`), 3-sample VCFs only.
    """
    variants = []
    csq_fields: List[str] = []
    ann_fields: List[str] = []
    sample_ids: List[str] = []
    proband_idx: Optional[int] = None
    parent_idxs: List[int] = []

    try:
        f_handle = _open_vcf(vcf_path)
    except FileNotFoundError:
        logger.error(f"VCF file not found: {vcf_path}")
        return variants
    except OSError as e:
        logger.error(f"VCF file could not be opened ({vcf_path}): {e}")
        return variants

    with f_handle as f:
        for line in f:
            # Parse annotation format headers before data lines
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = parse_csq_header(line)
                logger.debug(f"Detected VEP CSQ header with {len(csq_fields)} fields")
            elif line.startswith("##INFO=<ID=ANN"):
                ann_fields = parse_ann_header(line)
                logger.debug(f"Detected SnpEff ANN header with {len(ann_fields)} fields")

            # Capture #CHROM header to discover sample columns and map
            # proband + parents for trio-aware de novo detection on the
            # data lines that follow. When ``ped_path`` is supplied we
            # use PED-based resolution in STRICT mode (fails loud if the
            # PED can't be matched against the VCF samples). Otherwise
            # we fall back to the legacy filename heuristic, which only
            # applies to exactly-3-sample VCFs.
            if line.startswith("#CHROM"):
                header_fields = line.rstrip("\n").split("\t")
                if len(header_fields) > 9:
                    sample_ids = header_fields[9:]
                    if ped_path:
                        # Import locally so repos without PED helpers in
                        # an older checkout don't break on import-time.
                        from scripts.intake.parse_ped import parse_ped as _parse_ped
                        from scripts.intake.parse_ped import resolve_trio as _resolve_trio

                        ped_entries = _parse_ped(ped_path)
                        proband_idx, parent_idxs = _resolve_trio(ped_entries, sample_ids)
                        if proband_idx is None:
                            raise ValueError(
                                f"PED file {ped_path} did not resolve a trio against "
                                f"VCF samples {sample_ids}. Strict mode — no fallback. "
                                f"Remove --ped to use the filename heuristic, or fix "
                                f"the PED so one affected individual has both parents "
                                f"listed with IDs present in the VCF."
                            )
                        logger.info(
                            "PED-based trio: proband=%s, parents=%s (from %s)",
                            sample_ids[proband_idx],
                            [sample_ids[i] for i in parent_idxs],
                            ped_path,
                        )
                    else:
                        proband_idx, parent_idxs = _detect_trio_proband(vcf_path, sample_ids)
                continue

            if line.startswith("#"):
                continue

            try:
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                chrom = fields[0] if fields[0].startswith("chr") else f"chr{fields[0]}"
                pos = int(fields[1])
                rsid = fields[2] if len(fields) > 2 and fields[2] != "." else None
                ref = fields[3]
                alt = fields[4]
                gene = None
                # Trio-aware de novo INFO flags. ``DN`` / ``DENOVO`` mark an
                # assumed de novo (PM6 at classification time); ``CONFIRMED_DN``
                # and ``PS2`` mark a parental-identity-verified de novo
                # (PS2_Moderate). Absence of any flag leaves ``inheritance`` at
                # ``None`` — ACMG 2015 says absence is not evidence of
                # inheritance, only absence of information.
                inheritance_label: Optional[str] = None
                confirmed_denovo = False
                if len(fields) > 7:
                    info = fields[7]
                    for item in info.split(";"):
                        if item.startswith("Gene="):
                            gene = item.split("=")[1]
                        # Boolean flags appear bare or as KEY=1 (VCF spec).
                        key = item.split("=", 1)[0]
                        val = item.split("=", 1)[1] if "=" in item else "1"
                        if key in ("DN", "DENOVO") and val not in ("0", "false", "False", ""):
                            if inheritance_label is None:
                                inheritance_label = "de_novo"
                        elif key in ("CONFIRMED_DN", "PS2") and val not in ("0", "false", "False", ""):
                            confirmed_denovo = True
                            inheritance_label = "confirmed_de_novo"

                # Trio-aware de novo fallback (v2.3-T9). INFO flags take
                # priority — we only consult FORMAT/GT when an upstream trio
                # caller (triodenovo, DeepTrio, etc.) did not tag the INFO
                # field but did emit trio genotypes. Genotype concordance
                # alone is at most assumed de novo (PM6); the CONFIRMED_DN
                # / PS2 INFO path is still the only way to reach
                # PS2_Moderate downstream.
                if proband_idx is not None and inheritance_label is None and len(fields) >= 10:
                    format_col = fields[8] if len(fields) > 8 else ""
                    format_keys = format_col.split(":") if format_col else []
                    sample_cols = fields[9:]
                    if proband_idx < len(sample_cols) and all(p < len(sample_cols) for p in parent_idxs):
                        trio_label = _detect_trio_denovo_from_genotype(
                            format_keys,
                            sample_cols[proband_idx],
                            sample_cols[parent_idxs[0]],
                            sample_cols[parent_idxs[1]],
                            ref=ref,
                            alt=alt,
                        )
                        if trio_label:
                            inheritance_label = trio_label

                # Parse annotations from INFO field
                annotation = None
                if len(fields) > 7:
                    info = fields[7]
                    for item in info.split(";"):
                        if item.startswith("CSQ=") and csq_fields:
                            annotation = parse_csq_value(item[4:], csq_fields, gene)
                            break
                        elif item.startswith("ANN=") and ann_fields:
                            annotation = parse_ann_value(item[4:], ann_fields, gene)
                            break

                variant = Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, gene=gene, rsid=rsid)
                if inheritance_label is not None:
                    variant.inheritance = inheritance_label
                if confirmed_denovo:
                    variant.confirmed_denovo = True
                if annotation:
                    variant.hgvsc = annotation.get("hgvsc") or None
                    variant.hgvsp = annotation.get("hgvsp") or None
                    raw_consequence = annotation.get("consequence", "")
                    variant.consequence = format_consequence(raw_consequence) or None
                    variant.transcript = annotation.get("transcript") or None
                    variant.impact = annotation.get("impact") or None
                    variant.sift = annotation.get("sift") or None
                    variant.polyphen = annotation.get("polyphen") or None
                    # In silico prediction scores (from VEP dbNSFP plugin)
                    _in_silico = {}
                    for _isf in (
                        "revel_score",
                        "cadd_phred",
                        "am_class",
                        "am_pathogenicity",
                        "spliceai_pred_ds_ag",
                        "spliceai_pred_ds_al",
                        "spliceai_pred_ds_dg",
                        "spliceai_pred_ds_dl",
                        "domains",
                    ):
                        _val = annotation.get(_isf)
                        if _val:
                            _in_silico[_isf] = _val
                    if _in_silico:
                        variant.in_silico = _in_silico
                    # Override gene from annotation if not in INFO
                    if not gene and annotation.get("gene"):
                        variant.gene = annotation["gene"]
                    # Fill rsID from CSQ Existing_variation if VCF ID column was "."
                    if not variant.rsid and annotation.get("rsid"):
                        variant.rsid = annotation["rsid"]

                variants.append(variant)
            except (ValueError, IndexError) as e:
                logger.warning(f"Skipping malformed VCF line: {line.strip()!r} — {e}")

    if len(variants) > 1000:
        logger.warning(f"VCF contains {len(variants)} variants (>1000). Consider filtering first.")

    # v2.3-T9: trio de novo detection telemetry so the reviewer can see
    # whether the FORMAT/GT fallback caught anything on a raw triodenovo
    # VCF without INFO flags.
    if proband_idx is not None:
        dn_count = sum(1 for v in variants if getattr(v, "inheritance", None) in ("de_novo", "confirmed_de_novo"))
        confirmed_count = sum(1 for v in variants if getattr(v, "confirmed_denovo", False))
        if dn_count:
            logger.info(
                "Trio de novo detection: %d variants flagged as de novo (%d confirmed). Proband sample was %s.",
                dn_count,
                confirmed_count,
                sample_ids[proband_idx] if sample_ids else "<unknown>",
            )

    # v2.3-T10: fail-loud warning when the VCF carries no annotation. BIKO's
    # classification, variant selector, and rare-disease / cancer tiering
    # arms all depend on gene + consequence + HGVSp being populated. On a
    # raw call-set (no VEP / SnpEff / ANNOVAR upstream) every variant falls
    # through to VUS, producing a technically-valid but clinically useless
    # report. Previously this failure mode was silent — the reviewer got a
    # green-ish report with 0 Tier I variants and no hint that the pipeline
    # was starved for input. This warning makes the gap visible in the log
    # stream so orchestrate.py, CI, and the human reviewer all see the
    # root cause immediately. The one feature that still works on
    # unannotated input is the trio FORMAT/GT de novo detection above; the
    # ACMG classification arms and the rare-disease selector require real
    # annotation to do anything meaningful.
    if variants:
        gene_hits = sum(1 for v in variants if v.gene)
        consequence_hits = sum(1 for v in variants if getattr(v, "consequence", None))
        if gene_hits == 0 and consequence_hits == 0:
            logger.warning(
                "⚠️  BIKO: %d variants parsed but NONE carry annotation "
                "(no Gene=, CSQ=, or ANN= INFO field; 0 gene + 0 consequence). "
                "Classification will fall through to 100%% VUS and the "
                "variant selector will not admit any variant. Run VEP / "
                "SnpEff / ANNOVAR on this VCF first. Trio FORMAT de novo "
                "detection is the only v2.3 feature that still works on "
                "raw callsets.",
                len(variants),
            )
        elif gene_hits / max(len(variants), 1) < 0.1:
            logger.warning(
                "⚠️  BIKO: only %d of %d variants (%.1f%%) carry a gene "
                "annotation. Classification coverage will be very limited. "
                "Verify the VCF was run through a full annotator before "
                "relying on these results.",
                gene_hits,
                len(variants),
                100.0 * gene_hits / len(variants),
            )

    return variants


if __name__ == "__main__":
    import json
    import sys

    vcf_path = sys.argv[1] if len(sys.argv) > 1 else "data/sample_vcf/demo_variants.vcf"
    variants = parse_vcf(vcf_path)
    print(json.dumps([{"variant_id": v.variant_id, "gene": v.gene} for v in variants], indent=2))
