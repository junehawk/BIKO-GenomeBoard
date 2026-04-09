# tests/test_intervar.py
"""Tests for InterVar TSV parser."""
import pytest
from pathlib import Path

from scripts.common.models import Variant
from scripts.intake.parse_intervar import parse_intervar, get_intervar_evidence

SAMPLE_FILE = str(Path(__file__).parent.parent / "data" / "sample_intervar" / "sample_intervar_output.tsv")


# ── parse_intervar ────────────────────────────────────────────────────────────

def test_parse_intervar_basic():
    """Parse sample file — verify variant count and specific evidence codes."""
    data = parse_intervar(SAMPLE_FILE)
    assert len(data) == 15

    # TP53 R248L: PM1, PM2, PM5, PP2, PP3
    tp53 = data["chr17:7675088:C>A"]
    assert "PM1" in tp53
    assert "PM2" in tp53
    assert "PM5" in tp53
    assert "PP2" in tp53
    assert "PP3" in tp53
    assert "PVS1" not in tp53

    # BRCA1 pathogenic: PVS1, PS1, PM2
    brca1 = data["chr17:43093449:C>T"]
    assert "PVS1" in brca1
    assert "PS1" in brca1
    assert "PM2" in brca1

    # BRCA2 VUS: PM2, PP3
    brca2 = data["chr13:32340300:G>A"]
    assert brca2 == ["PM2", "PP3"]


def test_parse_intervar_benign_variant():
    """Benign variant should have BA1 only."""
    data = parse_intervar(SAMPLE_FILE)
    bard1 = data["chr2:215645464:C>G"]
    assert bard1 == ["BA1"]


def test_parse_intervar_likely_benign():
    """Likely benign should have BS1 + BP4 or BS1 + BP7."""
    data = parse_intervar(SAMPLE_FILE)
    # ATM D1853N: BS1 + BP4
    atm = data["chr11:108236086:G>A"]
    assert "BS1" in atm
    assert "BP4" in atm
    # CDH1 synonymous: BS1 + BP7
    cdh1 = data["chr16:68842399:G>A"]
    assert "BS1" in cdh1
    assert "BP7" in cdh1


def test_parse_intervar_empty(tmp_path):
    """Empty file returns empty dict."""
    empty_file = tmp_path / "empty.tsv"
    empty_file.write_text("")
    data = parse_intervar(str(empty_file))
    assert data == {}


def test_parse_intervar_no_file():
    """Nonexistent file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        parse_intervar("/nonexistent/path/intervar.tsv")


def test_parse_intervar_all_zero(tmp_path):
    """Variant with all evidence codes at 0 returns empty list."""
    tsv = tmp_path / "all_zero.tsv"
    header = "#Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tExonicFunc.refGene\tAAChange.refGene\tInterVar: InterVar and target\tPVS1\tPS1\tPS2\tPS3\tPS4\tPM1\tPM2\tPM3\tPM4\tPM5\tPM6\tPP1\tPP2\tPP3\tPP4\tPP5\tBA1\tBS1\tBS2\tBS3\tBS4\tBP1\tBP2\tBP3\tBP4\tBP5\tBP6\tBP7\n"
    row = "chr1\t100\t100\tA\tG\texonic\tFAKE\tnonsynonymous SNV\tFAKE:NM_001:exon1:c.1A>G:p.K1E\tUncertain significance\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
    tsv.write_text(header + row)
    data = parse_intervar(str(tsv))
    assert "chr1:100:A>G" in data
    assert data["chr1:100:A>G"] == []


def test_parse_intervar_header_only(tmp_path):
    """File with only header returns empty dict."""
    tsv = tmp_path / "header_only.tsv"
    header = "#Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tExonicFunc.refGene\tAAChange.refGene\tInterVar: InterVar and target\tPVS1\tPS1\n"
    tsv.write_text(header)
    data = parse_intervar(str(tsv))
    assert data == {}


def test_parse_intervar_bom(tmp_path):
    """File with UTF-8 BOM should parse correctly."""
    tsv = tmp_path / "bom.tsv"
    header = "\ufeff#Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tExonicFunc.refGene\tAAChange.refGene\tInterVar: InterVar and target\tPVS1\tPS1\tPS2\tPS3\tPS4\tPM1\tPM2\tPM3\tPM4\tPM5\tPM6\tPP1\tPP2\tPP3\tPP4\tPP5\tBA1\tBS1\tBS2\tBS3\tBS4\tBP1\tBP2\tBP3\tBP4\tBP5\tBP6\tBP7\n"
    row = "chr1\t500\t500\tG\tT\texonic\tTP53\tnonsynonymous SNV\tTP53:NM_000546:exon5:c.1A>G:p.K1E\tPathogenic\t1\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
    tsv.write_text(header + row, encoding="utf-8")
    data = parse_intervar(str(tsv))
    assert "chr1:500:G>T" in data
    assert data["chr1:500:G>T"] == ["PVS1", "PS1"]


def test_parse_intervar_no_chr_prefix(tmp_path):
    """Chromosome without 'chr' prefix should be normalised."""
    tsv = tmp_path / "no_chr.tsv"
    header = "#Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tExonicFunc.refGene\tAAChange.refGene\tInterVar: InterVar and target\tPVS1\tPS1\tPS2\tPS3\tPS4\tPM1\tPM2\tPM3\tPM4\tPM5\tPM6\tPP1\tPP2\tPP3\tPP4\tPP5\tBA1\tBS1\tBS2\tBS3\tBS4\tBP1\tBP2\tBP3\tBP4\tBP5\tBP6\tBP7\n"
    row = "17\t7675088\t7675088\tC\tA\texonic\tTP53\tnonsynonymous SNV\tTP53:NM_000546:exon7:c.743G>T:p.R248L\tLikely pathogenic\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
    tsv.write_text(header + row)
    data = parse_intervar(str(tsv))
    assert "chr17:7675088:C>A" in data


# ── get_intervar_evidence ─────────────────────────────────────────────────────

def test_get_intervar_evidence_found():
    """Variant present in parsed data returns its codes."""
    data = parse_intervar(SAMPLE_FILE)
    v = Variant(chrom="chr17", pos=7675088, ref="C", alt="A", gene="TP53")
    codes = get_intervar_evidence(v, data)
    assert "PM1" in codes
    assert "PM2" in codes


def test_get_intervar_evidence_not_found():
    """Variant absent from parsed data returns empty list."""
    data = parse_intervar(SAMPLE_FILE)
    v = Variant(chrom="chr1", pos=999999, ref="A", alt="T")
    codes = get_intervar_evidence(v, data)
    assert codes == []


def test_get_intervar_evidence_empty_data():
    """Empty intervar_data dict returns empty list."""
    v = Variant(chrom="chr17", pos=7675088, ref="C", alt="A")
    assert get_intervar_evidence(v, {}) == []


def test_get_intervar_evidence_none_data():
    """None intervar_data returns empty list."""
    v = Variant(chrom="chr17", pos=7675088, ref="C", alt="A")
    assert get_intervar_evidence(v, None) == []
