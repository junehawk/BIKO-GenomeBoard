# tests/test_expanded_pgx.py
import json
from pathlib import Path

import pytest

from scripts.common.models import Variant
from scripts.pharma.korean_pgx import check_korean_pgx


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _pgx_table():
    path = Path(__file__).parent.parent / "data" / "korean_pgx_table.json"
    with open(path) as f:
        return json.load(f)


def _config_genes():
    import yaml
    path = Path(__file__).parent.parent / "config.yaml"
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return cfg["pgx"]["genes"]


# ---------------------------------------------------------------------------
# 1. Table completeness
# ---------------------------------------------------------------------------

def test_pgx_table_has_12_genes():
    data = _pgx_table()
    assert len(data["genes"]) == 12


# ---------------------------------------------------------------------------
# 2. New gene phenotypes
# ---------------------------------------------------------------------------

def test_check_korean_pgx_cyp3a5():
    variant = Variant(chrom="chr7", pos=99245146, ref="A", alt="G", gene="CYP3A5")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.gene == "CYP3A5"
    assert result.cpic_level == "A"
    assert "tacrolimus" in result.phenotype.lower()


def test_check_korean_pgx_ugt1a1():
    variant = Variant(chrom="chr2", pos=234668879, ref="G", alt="A", gene="UGT1A1")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.gene == "UGT1A1"
    assert result.cpic_level == "A"
    assert "irinotecan" in result.phenotype.lower()


def test_check_korean_pgx_vkorc1():
    """VKORC1 -1639A allele freq 0.90 vs 0.37 → ratio ~2.43 → korean_flag=True."""
    variant = Variant(chrom="chr16", pos=31096368, ref="G", alt="A", gene="VKORC1")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.gene == "VKORC1"
    assert result.korean_flag is True
    assert "warfarin" in result.phenotype.lower()


def test_check_korean_pgx_slco1b1():
    """SLCO1B1 *15 freq 0.16 vs 0.15 → ratio ~1.07 → korean_flag=False (no significant enrichment)."""
    variant = Variant(chrom="chr12", pos=21178615, ref="T", alt="C", gene="SLCO1B1")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.gene == "SLCO1B1"
    assert result.korean_flag is False
    assert "statin" in result.phenotype.lower()


def test_check_korean_pgx_g6pd():
    variant = Variant(chrom="chrX", pos=154531391, ref="G", alt="T", gene="G6PD")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.gene == "G6PD"
    assert result.cpic_level == "A"
    assert "rasburicase" in result.phenotype.lower()


def test_check_korean_pgx_ifnl3():
    variant = Variant(chrom="chr19", pos=39738525, ref="C", alt="T", gene="IFNL3")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.gene == "IFNL3"
    assert result.cpic_level == "A"
    assert "ifnl3" in result.phenotype.lower()


def test_check_korean_pgx_cyp1a2():
    variant = Variant(chrom="chr15", pos=74749576, ref="G", alt="A", gene="CYP1A2")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.gene == "CYP1A2"
    assert "cyp1a2" in result.phenotype.lower()


# ---------------------------------------------------------------------------
# 3. Config completeness
# ---------------------------------------------------------------------------

def test_config_pgx_genes_includes_new():
    genes = _config_genes()
    expected = [
        "CYP2D6", "CYP2C19", "CYP2C9", "HLA-B", "HLA-A",
        "NUDT15", "TPMT", "DPYD",
        "CYP3A5", "UGT1A1", "SLCO1B1", "VKORC1", "CYP1A2", "G6PD", "IFNL3",
    ]
    for gene in expected:
        assert gene in genes, f"{gene} missing from config pgx.genes"
    assert len(genes) == 15


# ---------------------------------------------------------------------------
# 4. Existing genes still work
# ---------------------------------------------------------------------------

def test_existing_pgx_genes_still_work_cyp2c19():
    variant = Variant(chrom="chr10", pos=96541616, ref="G", alt="A", gene="CYP2C19")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.gene == "CYP2C19"
    assert result.phenotype == "Intermediate Metabolizer (*2 carrier)"
    assert result.korean_flag is True
    assert result.cpic_level == "A"


def test_existing_pgx_genes_still_work_hla_b():
    variant = Variant(chrom="chr6", pos=31353875, ref="C", alt="A", gene="HLA-B")
    result = check_korean_pgx(variant)
    assert result is not None
    assert result.phenotype == "HLA-B*5701 carrier — abacavir hypersensitivity risk"
    assert result.cpic_level == "A"
