# tests/conftest.py
"""Shared pytest fixtures for BIKO GenomeBoard tests.

Naming convention for new tests (adopt incrementally — do not rename
existing tests):

    test_<module>_<behavior>_<condition>

Examples:
    test_pgx_table_includes_all_pharmcat_genes
    test_clinvar_override_skips_when_pm5_absent

Fixture conventions:
    sample_*       — plain dataclass instances (no I/O)
    mock_*         — dict payloads that mimic external API shapes
    *_vcf_path     — on-disk VCFs written to tmp_path (function-scoped)
    *_db           — on-disk SQLite DBs written to tmp_path
"""

from pathlib import Path

import pytest

from scripts.common.models import FrequencyData, Variant

_REPO_ROOT = Path(__file__).resolve().parent.parent
_DEMO_ANNOTATED_VCF = _REPO_ROOT / "data" / "sample_vcf" / "demo_variants_grch38_annotated.vcf"


@pytest.fixture
def sample_tp53_variant():
    return Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")


@pytest.fixture
def sample_brca2_variant():
    return Variant(chrom="chr13", pos=32337326, ref="C", alt="T", gene="BRCA2")


@pytest.fixture
def sample_cyp2c19_variant():
    return Variant(chrom="chr10", pos=96541616, ref="G", alt="A", gene="CYP2C19")


@pytest.fixture
def mock_clinvar_pathogenic():
    return {
        "clinical_significance": {"description": "Pathogenic"},
        "gene": {"symbol": "TP53"},
        "trait_set": [{"trait_name": "Li-Fraumeni syndrome"}],
        "review_status": "criteria provided, multiple submitters, no conflicts",
        "variation_id": "12375",
        "molecular_consequence": "nonsense",
    }


@pytest.fixture
def mock_clinvar_not_found():
    return None


@pytest.fixture
def mock_gnomad_response():
    return {"data": {"variant": {"genome": {"af": 0.0002, "populations": [{"id": "eas", "af": 0.0003}]}}}}


@pytest.fixture
def mock_gnomad_not_found():
    return {"data": {"variant": None}}


@pytest.fixture
def sample_frequency_data():
    return FrequencyData(kova=0.0001, gnomad_eas=0.0003, gnomad_all=0.0002)


@pytest.fixture
def demo_annotated_vcf_path() -> str:
    """Absolute path to the canonical pre-annotated demo VCF.

    This is the fixture to use instead of hard-coding
    `"data/sample_vcf/demo_variants_grch38_annotated.vcf"`; currently
    4 test files repeat that literal path.
    """
    return str(_DEMO_ANNOTATED_VCF)


@pytest.fixture
def sample_vcf_path(tmp_path) -> Path:
    """Write a minimal 3-variant VCF (TP53 / BRCA2 / CYP2C19) to tmp_path.

    Function-scoped so each test gets a fresh file. Use this when a
    test only needs a parseable VCF and does not care about annotations.
    """
    vcf = tmp_path / "sample.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr10>\n"
        "##contig=<ID=chr13>\n"
        "##contig=<ID=chr17>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr17\t7577120\t.\tG\tA\t100\tPASS\t.\n"
        "chr13\t32337326\t.\tC\tT\t100\tPASS\t.\n"
        "chr10\t96541616\t.\tG\tA\t100\tPASS\t.\n"
    )
    return vcf


@pytest.fixture
def trio_vcf_path(tmp_path) -> Path:
    """Write a 3-sample trio VCF (PROBAND / FATHER / MOTHER) to tmp_path.

    Proband is 0/1 heterozygous, both parents are 0/0 — i.e. a de novo
    genotype pattern suitable for exercising the PS2/PM6 path.
    """
    vcf = tmp_path / "trio.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr17>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND\tFATHER\tMOTHER\n"
        "chr17\t7577120\t.\tG\tA\t100\tPASS\t.\tGT\t0/1\t0/0\t0/0\n"
    )
    return vcf


@pytest.fixture
def mock_clinvar_db(tmp_path) -> Path:
    """Write a minimal 2-row ClinVar SQLite DB matching the local schema.

    Covers the TP53 R175H Pathogenic fixture plus one BRCA2 Benign row.
    Scope is function so tests can mutate the DB without cross-talk.
    """
    import sqlite3

    db = tmp_path / "clinvar.sqlite3"
    conn = sqlite3.connect(db)
    conn.execute(
        "CREATE TABLE clinvar ("
        "chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, gene TEXT, "
        "clinical_significance TEXT, review_status TEXT, "
        "variation_id TEXT, molecular_consequence TEXT)"
    )
    conn.executemany(
        "INSERT INTO clinvar VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        [
            ("chr17", 7577120, "G", "A", "TP53", "Pathogenic", "criteria provided", "12375", "missense"),
            ("chr13", 32337326, "C", "T", "BRCA2", "Benign", "criteria provided", "99999", "synonymous"),
        ],
    )
    conn.commit()
    conn.close()
    return db
