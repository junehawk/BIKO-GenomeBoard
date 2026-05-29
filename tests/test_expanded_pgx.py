# tests/test_expanded_pgx.py
import json
from pathlib import Path

import pytest

from scripts.common.models import Variant
from scripts.pharmacogenomics.korean_pgx import PGX_GENES, check_korean_pgx

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


def _default_phenotype(gene: str) -> str:
    """The curated default_phenotype for *gene* from the JSON table.

    v2.7 PGX-02: the gene-only builtin path no longer surfaces this as the
    patient's phenotype (it would falsely assert a genotype) — it is the
    curated reference the PharmCAT path uses once a diplotype is actually
    called. These tests assert the curated data is intact here, separately
    from check_korean_pgx's (now non-attributive) behavior.
    """
    for e in _pgx_table()["genes"]:
        if e.get("gene") == gene:
            return e.get("default_phenotype", "")
    return ""


def _assert_builtin_gene_level(result, gene: str) -> None:
    """The builtin (gene-symbol-only) path must be non-attributive (PGX-02)."""
    assert result is not None, f"{gene} check_korean_pgx returned None"
    assert result.gene == gene
    assert result.star_allele == "", f"{gene} builtin path must not assert a star allele"
    assert gene.lower() in result.phenotype.lower()
    assert "pharmcat" in result.phenotype.lower()
    assert "carrier" not in result.phenotype.lower(), f"{gene} must not claim a carrier from a gene match"
    assert "metabolizer" not in result.phenotype.lower(), f"{gene} must not assert a metabolizer phenotype"


# ---------------------------------------------------------------------------
# 1. Table completeness
# ---------------------------------------------------------------------------


def test_pgx_table_has_24_genes():
    data = _pgx_table()
    assert len(data["genes"]) == 24


# ---------------------------------------------------------------------------
# 2. New gene phenotypes
# ---------------------------------------------------------------------------


# Each gene's builtin result must be gene-level / non-attributive (PGX-02),
# while gene-level facts (cpic_level, korean_flag) are preserved, and the drug
# keyword still lives in the curated JSON default_phenotype (PharmCAT path).
def test_check_korean_pgx_cyp3a5():
    variant = Variant(chrom="chr7", pos=99245146, ref="A", alt="G", gene="CYP3A5")
    result = check_korean_pgx(variant)
    _assert_builtin_gene_level(result, "CYP3A5")
    assert result.cpic_level == "A"
    assert "tacrolimus" in _default_phenotype("CYP3A5").lower()


def test_check_korean_pgx_ugt1a1():
    variant = Variant(chrom="chr2", pos=234668879, ref="G", alt="A", gene="UGT1A1")
    result = check_korean_pgx(variant)
    _assert_builtin_gene_level(result, "UGT1A1")
    assert result.cpic_level == "A"
    assert "irinotecan" in _default_phenotype("UGT1A1").lower()


def test_check_korean_pgx_vkorc1():
    """VKORC1 -1639A allele freq 0.90 vs 0.37 → ratio ~2.43 → korean_flag=True."""
    variant = Variant(chrom="chr16", pos=31096368, ref="G", alt="A", gene="VKORC1")
    result = check_korean_pgx(variant)
    _assert_builtin_gene_level(result, "VKORC1")
    assert result.korean_flag is True
    assert "warfarin" in _default_phenotype("VKORC1").lower()


def test_check_korean_pgx_slco1b1():
    """SLCO1B1 *15 freq 0.16 vs 0.15 → ratio ~1.07 → korean_flag=False (no significant enrichment)."""
    variant = Variant(chrom="chr12", pos=21178615, ref="T", alt="C", gene="SLCO1B1")
    result = check_korean_pgx(variant)
    _assert_builtin_gene_level(result, "SLCO1B1")
    assert result.korean_flag is False
    assert "statin" in _default_phenotype("SLCO1B1").lower()


def test_check_korean_pgx_g6pd():
    variant = Variant(chrom="chrX", pos=154531391, ref="G", alt="T", gene="G6PD")
    result = check_korean_pgx(variant)
    _assert_builtin_gene_level(result, "G6PD")
    assert result.cpic_level == "A"
    assert "rasburicase" in _default_phenotype("G6PD").lower()


def test_check_korean_pgx_ifnl3():
    variant = Variant(chrom="chr19", pos=39738525, ref="C", alt="T", gene="IFNL3")
    result = check_korean_pgx(variant)
    _assert_builtin_gene_level(result, "IFNL3")
    assert result.cpic_level == "A"
    assert "ifnl3" in _default_phenotype("IFNL3").lower()


def test_check_korean_pgx_cyp1a2():
    variant = Variant(chrom="chr15", pos=74749576, ref="G", alt="A", gene="CYP1A2")
    result = check_korean_pgx(variant)
    _assert_builtin_gene_level(result, "CYP1A2")
    assert "cyp1a2" in _default_phenotype("CYP1A2").lower()


# ---------------------------------------------------------------------------
# 3. Config completeness
# ---------------------------------------------------------------------------


def test_config_pgx_genes_includes_new():
    genes = _config_genes()
    expected = [
        "CYP2D6",
        "CYP2C19",
        "CYP2C9",
        "HLA-B",
        "HLA-A",
        "NUDT15",
        "TPMT",
        "DPYD",
        "CYP3A5",
        "UGT1A1",
        "SLCO1B1",
        "VKORC1",
        "CYP1A2",
        "G6PD",
        "IFNL3",
        "CYP2B6",
        "CYP4F2",
        "ABCG2",
        "NAT2",
        "CACNA1S",
        "CFTR",
        "CYP3A4",
        "MT-RNR1",
        "RYR1",
    ]
    for gene in expected:
        assert gene in genes, f"{gene} missing from config pgx.genes"
    assert len(genes) == 24


# ---------------------------------------------------------------------------
# 4. Existing genes still work
# ---------------------------------------------------------------------------


def test_existing_pgx_genes_still_work_cyp2c19():
    variant = Variant(chrom="chr10", pos=96541616, ref="G", alt="A", gene="CYP2C19")
    result = check_korean_pgx(variant)
    _assert_builtin_gene_level(result, "CYP2C19")
    assert result.korean_flag is True
    assert result.cpic_level == "A"
    # Curated diplotype phenotype remains in the JSON for the PharmCAT path.
    assert _default_phenotype("CYP2C19") == "Intermediate Metabolizer (*2 carrier)"


def test_existing_pgx_genes_still_work_hla_b():
    variant = Variant(chrom="chr6", pos=31353875, ref="C", alt="A", gene="HLA-B")
    result = check_korean_pgx(variant)
    _assert_builtin_gene_level(result, "HLA-B")
    assert result.cpic_level == "A"
    assert _default_phenotype("HLA-B") == "HLA-B*5701 carrier — abacavir hypersensitivity risk"


# ---------------------------------------------------------------------------
# 5. JSON-driven phenotype — regression + new coverage
# ---------------------------------------------------------------------------

# Snapshot of the legacy 12-gene phenotype strings. These MUST remain
# byte-identical to the hardcoded elif chain that existed before the
# v2.4 Quick Win B refactor — any accidental change is a regression.
_LEGACY_12_GENE_PHENOTYPE_SNAPSHOTS: dict[str, str] = {
    "CYP2C19": "Intermediate Metabolizer (*2 carrier)",
    "HLA-B": "HLA-B*5701 carrier — abacavir hypersensitivity risk",
    "NUDT15": "NUDT15 intermediate metabolizer (p.R139C carrier)",
    "CYP3A5": "CYP3A5 expressor (*1 carrier) — higher tacrolimus dose needed",
    "UGT1A1": "UGT1A1 poor metabolizer (*6 carrier) — irinotecan toxicity risk",
    "SLCO1B1": "SLCO1B1 decreased function (*15 carrier) — statin myopathy risk",
    "VKORC1": "VKORC1 low-dose warfarin phenotype (-1639A carrier)",
    "G6PD": "G6PD deficient — rasburicase contraindicated",
    "IFNL3": "IFNL3 favorable responder (CC genotype)",
    "CYP1A2": "CYP1A2 poor metabolizer (*1C carrier)",
}


@pytest.mark.parametrize("gene,expected_phenotype", sorted(_LEGACY_12_GENE_PHENOTYPE_SNAPSHOTS.items()))
def test_legacy_12_genes_curated_phenotype_unchanged(gene: str, expected_phenotype: str):
    """Regression: the curated default_phenotype strings in the JSON must stay
    byte-identical (the PharmCAT path surfaces them). v2.7 PGX-02: the gene-only
    builtin path no longer ASSERTS these — it returns a non-attributive
    gene-level result — but the curated data itself must not drift."""
    assert _default_phenotype(gene) == expected_phenotype, (
        f"{gene} curated phenotype drifted from snapshot: {_default_phenotype(gene)!r} != {expected_phenotype!r}"
    )
    # And the builtin path must NOT assert this phenotype on a gene match.
    variant = Variant(chrom="chr1", pos=1, ref="A", alt="T", gene=gene)
    _assert_builtin_gene_level(check_korean_pgx(variant), gene)


# Keyword probes for the 12 newly-curated genes. Each entry asserts
# that a non-empty phenotype is returned AND that it contains a
# clinically meaningful keyword — this catches both dead-code regressions
# (phenotype missing entirely) and empty-string drift.
_NEW_12_GENE_KEYWORDS: list[tuple[str, str]] = [
    ("DPYD", "fluoropyrimidine"),
    ("TPMT", "thiopurine"),
    ("HLA-A", "carbamazepine"),
    ("CYP2B6", "efavirenz"),
    ("CYP4F2", "warfarin"),
    ("ABCG2", "rosuvastatin"),
    ("NAT2", "isoniazid"),
    ("CACNA1S", "volatile anesthetics"),
    ("CFTR", "ivacaftor"),
    ("CYP3A4", "tacrolimus"),
    ("MT-RNR1", "aminoglycosides"),
    ("RYR1", "volatile anesthetics"),
]


@pytest.mark.parametrize("gene,expected_keyword", _NEW_12_GENE_KEYWORDS)
def test_new_12_genes_curated_phenotype_populated(gene: str, expected_keyword: str):
    """Every newly-curated gene must have a non-empty curated default_phenotype
    in the JSON containing the expected drug/mechanism keyword, and the builtin
    path must return a (non-attributive) gene-level result for it."""
    curated = _default_phenotype(gene)
    assert curated != "", f"{gene} has empty curated default_phenotype in the JSON"
    assert expected_keyword.lower() in curated.lower(), (
        f"{gene} curated phenotype missing expected keyword {expected_keyword!r}: {curated!r}"
    )
    variant = Variant(chrom="chr1", pos=1, ref="A", alt="T", gene=gene)
    _assert_builtin_gene_level(check_korean_pgx(variant), gene)


def test_pgx_genes_set_covers_all_24():
    """PGX_GENES must include every gene in the curated 24-gene panel."""
    expected = {
        "CYP2D6",
        "CYP2C19",
        "CYP2C9",
        "HLA-B",
        "HLA-A",
        "NUDT15",
        "TPMT",
        "DPYD",
        "CYP3A5",
        "UGT1A1",
        "SLCO1B1",
        "VKORC1",
        "CYP1A2",
        "G6PD",
        "IFNL3",
        "CYP2B6",
        "CYP4F2",
        "ABCG2",
        "NAT2",
        "CACNA1S",
        "CFTR",
        "CYP3A4",
        "MT-RNR1",
        "RYR1",
    }
    missing = expected - PGX_GENES
    assert not missing, f"PGX_GENES missing expected genes: {sorted(missing)}"


def test_check_korean_pgx_unknown_gene_returns_none():
    """Genes outside the PGx panel return None rather than a stub PgxResult."""
    variant = Variant(chrom="chr1", pos=1, ref="A", alt="T", gene="FOOBAR")
    assert check_korean_pgx(variant) is None


def test_default_phenotype_field_present_all_entries():
    """Every JSON entry must have a non-empty ``default_phenotype``.
    Without this, the JSON-driven phenotype path silently degrades
    to an empty string, which is indistinguishable from a missing gene
    in the report template."""
    path = Path(__file__).parent.parent / "data" / "korean_pgx_table.json"
    with open(path) as f:
        data = json.load(f)
    for entry in data["genes"]:
        gene = entry.get("gene", "<unknown>")
        dp = entry.get("default_phenotype")
        assert dp is not None, f"{gene} missing default_phenotype field"
        assert isinstance(dp, str), f"{gene} default_phenotype is not a string: {type(dp)}"
        assert dp.strip(), f"{gene} default_phenotype is empty/whitespace"
