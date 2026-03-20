# tests/test_intake.py
import pytest
from pathlib import Path
from scripts.intake.parse_vcf import parse_vcf
from scripts.intake.parse_text import parse_text

DEMO_VCF = Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf"

def test_parse_vcf_returns_variants():
    variants = parse_vcf(str(DEMO_VCF))
    assert len(variants) == 10

def test_parse_vcf_first_variant():
    variants = parse_vcf(str(DEMO_VCF))
    v = variants[0]
    assert v.chrom == "chr17"
    assert v.pos == 7577120
    assert v.ref == "G"
    assert v.alt == "A"
    assert v.gene == "TP53"

def test_parse_vcf_warns_over_1000(tmp_path, caplog):
    # 1001 variant VCF to trigger warning
    vcf_file = tmp_path / "big.vcf"
    lines = ["##fileformat=VCFv4.1\n",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"]
    for i in range(1001):
        lines.append(f"chr1\t{i+1}\t.\tA\tT\t.\tPASS\tGene=TEST\n")
    vcf_file.write_text("".join(lines))
    variants = parse_vcf(str(vcf_file))
    assert len(variants) == 1001
    assert "1000" in caplog.text  # warning message

def test_parse_text_single():
    variants = parse_text("chr17:7577120 G>A")
    assert len(variants) == 1
    assert variants[0].variant_id == "chr17:7577120:G>A"

def test_parse_text_multiple():
    text = "chr17:7577120 G>A\nchr13:32337326 C>T"
    variants = parse_text(text)
    assert len(variants) == 2

def test_parse_text_comma_separated():
    text = "chr17:7577120 G>A, chr13:32337326 C>T"
    variants = parse_text(text)
    assert len(variants) == 2
