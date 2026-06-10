# tests/test_models.py
import pytest

from scripts.common.models import AcmgEvidence, FrequencyData, PgxResult, Variant


def test_variant_from_string():
    v = Variant.from_string("chr17:7577120 G>A")
    assert v.chrom == "chr17"
    assert v.pos == 7577120
    assert v.ref == "G"
    assert v.alt == "A"


def test_variant_from_string_with_chr_prefix():
    v = Variant.from_string("17:7577120 G>A")
    assert v.chrom == "chr17"  # auto chr prefix


def test_variant_id():
    v = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    assert v.variant_id == "chr17:7577120:G>A"


def test_acmg_evidence_creation():
    e = AcmgEvidence(code="PVS1", source="clinical_geneticist", description="null variant")
    assert e.strength == "very_strong"  # PVS1 → very_strong auto mapping


def test_acmg_evidence_strength_mapping():
    assert AcmgEvidence(code="PS1", source="test", description="").strength == "strong"
    assert AcmgEvidence(code="PM2_Supporting", source="test", description="").strength == "supporting"
    assert AcmgEvidence(code="BA1", source="test", description="").strength == "stand_alone"
    assert AcmgEvidence(code="BS1", source="test", description="").strength == "strong"


def test_frequency_data():
    f = FrequencyData(kova=0.001, gnomad_eas=0.002, gnomad_all=0.003)
    # Enrichment is now KOVA / gnomAD EAS (Korean vs broader East Asian)
    assert f.korean_vs_global_ratio() == pytest.approx(0.5, rel=0.01)
    assert f.korean_max == pytest.approx(0.001)


def test_pgx_result():
    p = PgxResult(
        gene="CYP2C19",
        star_allele="*2",
        phenotype="Poor Metabolizer",
        cpic_level="A",
        korean_prevalence=0.15,
        western_prevalence=0.03,
    )
    assert p.korean_flag is True  # >= 5x difference


def test_structural_variant_creation():
    from scripts.common.models import StructuralVariant

    sv = StructuralVariant(
        annotsv_id="CNV_ERBB2_AMP",
        chrom="chr17",
        start=37844393,
        end=37884925,
        length=40532,
        sv_type="DUP",
        sample_id="Sample1",
        acmg_class=5,
        ranking_score=0.99,
        cytoband="17q12",
        gene_name="ERBB2",
        gene_count=1,
    )
    assert sv.sv_type == "DUP"
    assert sv.acmg_class == 5
    assert sv.is_pathogenic
    assert not sv.is_benign
    assert sv.size_display == "40.5 kb"


def test_structural_variant_dosage_sensitive_del():
    from scripts.common.models import StructuralVariant

    sv = StructuralVariant(
        annotsv_id="test",
        chrom="chr17",
        start=1,
        end=100000,
        length=100000,
        sv_type="DEL",
        sample_id="S1",
        acmg_class=3,
        ranking_score=0.2,
        cytoband="17q21",
        gene_name="BRCA1",
        gene_count=1,
    )
    sv.gene_details = [{"gene": "BRCA1", "hi": 3, "ts": 0, "pli": 0.98}]
    assert sv.is_dosage_sensitive(mode="cancer")


def test_structural_variant_not_dosage_sensitive():
    from scripts.common.models import StructuralVariant

    sv = StructuralVariant(
        annotsv_id="test",
        chrom="chr1",
        start=1,
        end=100,
        length=100,
        sv_type="DUP",
        sample_id="S1",
        acmg_class=3,
        ranking_score=0.1,
        cytoband="1p36",
        gene_name="GENE1",
        gene_count=1,
    )
    sv.gene_details = [{"gene": "GENE1", "hi": 0, "ts": 0, "pli": 0.05}]
    assert not sv.is_dosage_sensitive(mode="cancer")


def _class3_sv(sv_type: str, gene_details: list) -> "object":
    """Build a Class-3 SV with the given split-row gene_details for dosage tests."""
    from scripts.common.models import StructuralVariant

    sv = StructuralVariant(
        annotsv_id="test",
        chrom="chr1",
        start=1,
        end=100,
        length=100,
        sv_type=sv_type,
        sample_id="S1",
        acmg_class=3,
        ranking_score=0.1,
        cytoband="1p36",
        gene_name="GENE1",
        gene_count=1,
    )
    sv.gene_details = gene_details
    return sv


def test_dosage_sensitive_ignores_clingen_sentinel_codes():
    """ClinGen HI/TS codes 30 (AR gene) and 40 (dosage sensitivity unlikely) are
    NOT real dosage-evidence scores — they must not pass the hi/ts >= threshold gate."""
    # DEL with HI=40 ("dosage sensitivity unlikely") — ClinGen says NOT haploinsufficient.
    assert not _class3_sv("DEL", [{"gene": "G", "hi": 40, "ts": 0, "pli": 0.0}]).is_dosage_sensitive(mode="cancer")
    # DEL with HI=30 ("gene associated with autosomal-recessive phenotype") — not HI.
    assert not _class3_sv("DEL", [{"gene": "G", "hi": 30, "ts": 0, "pli": 0.0}]).is_dosage_sensitive(mode="cancer")
    # DUP with TS=40 — not triplosensitive.
    assert not _class3_sv("DUP", [{"gene": "G", "hi": 0, "ts": 40, "pli": 0.0}]).is_dosage_sensitive(mode="cancer")
    # Same sentinels must also be rejected in the more permissive rare-disease mode (thresh=1).
    assert not _class3_sv("DEL", [{"gene": "G", "hi": 40, "ts": 0, "pli": 0.0}]).is_dosage_sensitive(
        mode="rare-disease"
    )


def test_dosage_sensitive_accepts_emerging_evidence_code():
    """Real ClinGen evidence scores in range (cancer >=2, rare-disease >=1) still flag DS."""
    # HI=2 ("emerging evidence") DEL — supporting at the cancer threshold.
    assert _class3_sv("DEL", [{"gene": "G", "hi": 2, "ts": 0, "pli": 0.0}]).is_dosage_sensitive(mode="cancer")
    # HI=1 ("little evidence") DEL — below cancer thresh (2) but at rare-disease thresh (1).
    assert not _class3_sv("DEL", [{"gene": "G", "hi": 1, "ts": 0, "pli": 0.0}]).is_dosage_sensitive(mode="cancer")
    assert _class3_sv("DEL", [{"gene": "G", "hi": 1, "ts": 0, "pli": 0.0}]).is_dosage_sensitive(mode="rare-disease")
