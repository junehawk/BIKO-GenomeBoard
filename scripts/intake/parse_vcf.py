# scripts/intake/parse_vcf.py
import logging
from typing import List
from scripts.common.models import Variant

logger = logging.getLogger(__name__)

def parse_vcf(vcf_path: str) -> List[Variant]:
    """Parse a VCF file into a list of Variant objects.
    Uses simple text parsing to avoid cyvcf2 dependency issues in tests.
    """
    variants = []
    try:
        f_handle = open(vcf_path)
    except FileNotFoundError:
        logger.error(f"VCF file not found: {vcf_path}")
        return variants
    with f_handle as f:
        for line in f:
            if line.startswith("#"):
                continue
            try:
                fields = line.strip().split("\t")
                if len(fields) < 5:
                    continue
                chrom = fields[0] if fields[0].startswith("chr") else f"chr{fields[0]}"
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                gene = None
                if len(fields) > 7:
                    info = fields[7]
                    for item in info.split(";"):
                        if item.startswith("Gene="):
                            gene = item.split("=")[1]
                variants.append(Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, gene=gene))
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
