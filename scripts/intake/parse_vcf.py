# scripts/intake/parse_vcf.py
import logging
from typing import List
from scripts.common.models import Variant
from scripts.intake.parse_annotation import (
    parse_csq_header, parse_ann_header,
    parse_csq_value, parse_ann_value,
    format_consequence,
)

logger = logging.getLogger(__name__)

def parse_vcf(vcf_path: str) -> List[Variant]:
    """Parse a VCF file into a list of Variant objects.
    Uses simple text parsing to avoid cyvcf2 dependency issues in tests.
    Supports VEP CSQ and SnpEff ANN annotation fields if present.
    """
    variants = []
    csq_fields: List[str] = []
    ann_fields: List[str] = []

    try:
        f_handle = open(vcf_path)
    except FileNotFoundError:
        logger.error(f"VCF file not found: {vcf_path}")
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
                if len(fields) > 7:
                    info = fields[7]
                    for item in info.split(";"):
                        if item.startswith("Gene="):
                            gene = item.split("=")[1]

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
                if annotation:
                    variant.hgvsc = annotation.get("hgvsc") or None
                    variant.hgvsp = annotation.get("hgvsp") or None
                    raw_consequence = annotation.get("consequence", "")
                    variant.consequence = format_consequence(raw_consequence) or None
                    variant.transcript = annotation.get("transcript") or None
                    variant.impact = annotation.get("impact") or None
                    variant.sift = annotation.get("sift") or None
                    variant.polyphen = annotation.get("polyphen") or None
                    # Override gene from annotation if not in INFO
                    if not gene and annotation.get("gene"):
                        variant.gene = annotation["gene"]

                variants.append(variant)
            except (ValueError, IndexError) as e:
                logger.warning(f"Skipping malformed VCF line: {line.strip()!r} — {e}")

    if len(variants) > 1000:
        logger.warning(f"VCF contains {len(variants)} variants (>1000). Consider filtering first.")
    return variants


if __name__ == "__main__":
    import sys, json
    vcf_path = sys.argv[1] if len(sys.argv) > 1 else "data/sample_vcf/demo_variants.vcf"
    variants = parse_vcf(vcf_path)
    print(json.dumps([{"variant_id": v.variant_id, "gene": v.gene} for v in variants], indent=2))
