# tests/test_intake.py
from pathlib import Path

from scripts.intake.parse_text import parse_text
from scripts.intake.parse_vcf import parse_vcf

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
    assert v.rsid == "rs28934578"


def test_parse_vcf_rsid_captured():
    """Verify rsID is captured from VCF column 3 for all variants."""
    variants = parse_vcf(str(DEMO_VCF))
    assert all(v.rsid is not None for v in variants)
    assert variants[1].rsid == "rs80358981"  # BRCA2


def test_parse_vcf_rsid_dot_becomes_none(tmp_path):
    """VCF lines with '.' in ID column should produce rsid=None."""
    vcf_file = tmp_path / "no_rsid.vcf"
    vcf_file.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t12345\t.\tA\tT\t.\tPASS\tGene=TEST\n"
    )
    variants = parse_vcf(str(vcf_file))
    assert len(variants) == 1
    assert variants[0].rsid is None


def test_parse_vcf_warns_over_1000(tmp_path, caplog):
    # 1001 variant VCF to trigger warning
    vcf_file = tmp_path / "big.vcf"
    lines = ["##fileformat=VCFv4.1\n", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"]
    for i in range(1001):
        lines.append(f"chr1\t{i + 1}\t.\tA\tT\t.\tPASS\tGene=TEST\n")
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


# ── gzip / bgzip VCF support ──────────────────────────────────────────────


def test_parse_vcf_gzipped_produces_same_variants(tmp_path):
    """Regression: production rare-disease / trio VCFs are shipped as
    .vcf.gz. Without transparent gzip open, parse_vcf would raise
    `UnicodeDecodeError: 0x8b in position 1` on the gzip magic byte.
    """
    import gzip

    plain = parse_vcf(str(DEMO_VCF))

    gz_path = tmp_path / "demo_variants.vcf.gz"
    with open(DEMO_VCF, "rb") as src, gzip.open(gz_path, "wb") as dst:
        dst.write(src.read())

    gzipped = parse_vcf(str(gz_path))

    assert len(gzipped) == len(plain) == 10
    assert [v.variant_id for v in gzipped] == [v.variant_id for v in plain]
    assert gzipped[0].chrom == "chr17"
    assert gzipped[0].pos == 7577120


def test_parse_vcf_bgzipped_extension_accepted(tmp_path):
    """The `.bgz` extension (bgzip output from tabix) must also be
    recognized — the bgzip format is gzip-compatible so `gzip.open` handles
    it transparently.
    """
    import gzip

    bgz_path = tmp_path / "demo_variants.vcf.bgz"
    with open(DEMO_VCF, "rb") as src, gzip.open(bgz_path, "wb") as dst:
        dst.write(src.read())

    variants = parse_vcf(str(bgz_path))
    assert len(variants) == 10


def test_parse_vcf_missing_file_returns_empty():
    """parse_vcf must not raise on missing files — it logs and returns []."""
    variants = parse_vcf("/nonexistent/path/to/nothing.vcf.gz")
    assert variants == []
