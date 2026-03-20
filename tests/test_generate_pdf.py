# tests/test_generate_pdf.py
from scripts.counselor.generate_pdf import generate_report_html, generate_pdf
from scripts.common.models import Variant


def test_generate_html():
    report_data = {
        "sample_id": "TEST_001",
        "date": "2026-03-20",
        "variants": [
            {
                "variant": "chr17:7577120:G>A",
                "gene": "TP53",
                "classification": "Pathogenic",
                "acmg_codes": ["PVS1", "PS1"],
                "agents": {
                    "clinical": {"clinvar_significance": "Pathogenic"},
                    "korean_pop": {"korean_flag": "한국인 희귀 변이"},
                },
                "conflict": False,
            }
        ],
        "pgx_results": [],
        "summary": {"total": 1, "pathogenic": 1, "vus": 0, "benign": 0},
        "db_versions": {"clinvar": "2026-03-15", "gnomad": "4.0"},
    }
    html = generate_report_html(report_data)
    assert "GenomeBoard" in html
    assert "TP53" in html
    assert "Pathogenic" in html
    assert "Research Use Only" in html
