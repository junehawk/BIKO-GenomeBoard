# tests/test_integration.py
from pathlib import Path

import pytest

from scripts.classification.acmg_engine import classify_variant
from scripts.clinical.query_clinvar import query_clinvar
from scripts.common.models import AcmgEvidence, FrequencyData
from scripts.counselor.generate_pdf import generate_report_html
from scripts.intake.parse_vcf import parse_vcf
from scripts.korean_pop.compare_freq import compare_frequencies
from scripts.korean_pop.query_gnomad import query_gnomad
from scripts.korean_pop.query_krgdb import _KRGDB_CACHE, query_krgdb
from scripts.pharma.korean_pgx import check_korean_pgx

pytestmark = pytest.mark.integration

DEMO_VCF_PATH = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf")


def test_full_pipeline_with_demo_vcf(mocker, tmp_path):
    """End-to-end: demo VCF -> parse -> query (mocked) -> classify -> report"""
    # Clear KRGDB cache to avoid stale data between test runs
    _KRGDB_CACHE.clear()

    # Mock ClinVar
    mocker.patch(
        "scripts.clinical.query_clinvar._search_clinvar_variant",
        return_value={
            "clinical_significance": {"description": "Pathogenic"},
            "gene": {"symbol": "TP53"},
            "review_status": "criteria provided, multiple submitters, no conflicts",
            "variation_id": "12375",
        },
    )
    # Mock gnomAD
    mocker.patch(
        "scripts.korean_pop.query_gnomad._graphql_query",
        return_value={"data": {"variant": {"genome": {"af": 0.0002, "populations": [{"id": "eas", "af": 0.0003}]}}}},
    )
    # Create temp KRGDB data
    krgdb_file = tmp_path / "krgdb_freq.tsv"
    krgdb_file.write_text("chr17\t7577120\tG\tA\t0.0001\n")

    # 1. Parse
    variants = parse_vcf(DEMO_VCF_PATH)
    assert len(variants) > 0

    # 2. For each variant: run all agents
    variant_results = []
    pgx_results = []
    for v in variants[:3]:
        # Clinical Geneticist
        clinical = query_clinvar(v)
        # Korean Pop Geneticist
        gnomad = query_gnomad(v)
        krgdb_freq = query_krgdb(v, str(krgdb_file))
        freq_data = FrequencyData(
            krgdb=krgdb_freq,
            gnomad_eas=gnomad["gnomad_eas"],
            gnomad_all=gnomad["gnomad_all"],
        )
        freq_result = compare_frequencies(freq_data)
        # Pharmacogenomicist
        pgx = check_korean_pgx(v)
        if pgx:
            pgx_results.append(pgx)

        # 3. Collect ACMG codes and classify
        all_codes = []
        for c in clinical["acmg_codes"]:
            all_codes.append(AcmgEvidence(code=c, source="clinical", description=""))
        for c in freq_result["acmg_codes"]:
            all_codes.append(AcmgEvidence(code=c, source="korean_pop", description=""))
        classification = classify_variant(all_codes, gene=v.gene)

        variant_results.append(
            {
                "variant": v.variant_id,
                "gene": v.gene,
                "classification": classification.classification,
                "acmg_codes": classification.evidence_codes,
                "agents": {
                    "clinical": clinical,
                    "korean_pop": {"korean_flag": freq_result["korean_flag"]},
                },
                "conflict": classification.conflict,
            }
        )

    # 4. Generate report
    report_data = {
        "sample_id": "TEST_001",
        "date": "2026-03-20",
        "variants": variant_results,
        "pgx_results": [{"gene": p.gene, "phenotype": p.phenotype, "cpic_level": p.cpic_level} for p in pgx_results],
        "summary": {
            "total": len(variant_results),
            "pathogenic": sum(1 for r in variant_results if r["classification"] == "Pathogenic"),
            "vus": sum(1 for r in variant_results if r["classification"] == "VUS"),
            "benign": sum(1 for r in variant_results if r["classification"] in ("Benign", "Likely Benign")),
        },
        "db_versions": {"clinvar": "2026-03-15", "gnomad": "4.0", "krgdb": "2026-03-01"},
    }
    html = generate_report_html(report_data)
    assert "BIKO GenomeBoard" in html
    assert "TP53" in html
    assert "Research Use Only" in html
