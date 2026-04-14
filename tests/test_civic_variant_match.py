import pytest
from scripts.db.build_civic_db import build_db
from scripts.db.query_civic import reset_civic_connection


@pytest.fixture
def civic_db(tmp_path, monkeypatch):
    """Build a temp CIViC DB with known test data."""
    # Create minimal TSV files
    civic_dir = tmp_path / "civic"
    civic_dir.mkdir()

    # Gene summaries
    (civic_dir / "GeneSummaries.tsv").write_text(
        "feature_id\tname\tdescription\tfeature_aliases\tentrez_id\tfeature_type\n"
        "1\tBRAF\tBRAF proto-oncogene\tBRAF1\t673\tGene\n"
        "2\tKRAS\tKRAS proto-oncogene\t-\t3845\tGene\n"
    )

    # Variant summaries
    (civic_dir / "VariantSummaries.tsv").write_text(
        "variant_id\tfeature_name\tvariant\tvariant_types\tchromosome\tstart\tstop\treference_bases\tvariant_bases\tentrez_id\n"
        "1\tBRAF\tV600E\tmissense\t7\t140753336\t140753336\tA\tT\t673\n"
        "2\tKRAS\tG12D\tmissense\t12\t25398284\t25398284\tG\tA\t3845\n"
    )

    # Evidence
    (civic_dir / "ClinicalEvidenceSummaries.tsv").write_text(
        "evidence_id\tmolecular_profile\tdisease\ttherapies\tevidence_type\tevidence_direction\tevidence_level\tsignificance\tevidence_statement\tcitation_id\tcitation\tnct_ids\trating\tvariant_origin\n"
        "1\tBRAF V600E\tMelanoma\tVemurafenib\tPredictive\tSupports\tA\tSensitivity/Response\tBRAF V600E predicts response to vemurafenib\t20979469\tChapman 2011\t\t5\tSomatic\n"
        "2\tBRAF V600E\tMelanoma\tDabrafenib\tPredictive\tSupports\tA\tSensitivity/Response\tBRAF V600E predicts response to dabrafenib\t22608338\tHauschild 2012\t\t4\tSomatic\n"
        "3\tBRAF V600E\tMelanoma\t\tPrognostic\tSupports\tB\tPoor Outcome\tBRAF V600E associated with poor prognosis\t26287849\tLong 2015\t\t3\tSomatic\n"
        "4\tKRAS G12D\tPancreatic cancer\tErlotinib\tPredictive\tSupports\tB\tResistance\tKRAS G12D predicts resistance\t18316791\tLievre 2008\t\t4\tSomatic\n"
        "5\tBRAF AMPLIFICATION\tColorectal cancer\tCetuximab\tPredictive\tSupports\tC\tResistance\tBRAF amp predicts resistance\t\tSmith 2019\t\t2\tSomatic\n"
    )

    db_path = str(tmp_path / "civic.sqlite3")
    build_db(str(civic_dir), db_path)

    monkeypatch.setattr("scripts.db.query_civic.get", lambda k, d=None: db_path if k == "paths.civic_db" else d)
    reset_civic_connection()

    yield db_path

    reset_civic_connection()


def test_get_predictive_evidence_variant_match(civic_db):
    """Variant-specific 매치 시 match_level='variant' 반환."""
    from scripts.db.query_civic import get_predictive_evidence_for_tier

    result = get_predictive_evidence_for_tier("BRAF", "p.Val600Glu")
    assert result["match_level"] == "variant"
    assert len(result["evidence"]) >= 1
    assert all(e["evidence_type"] == "Predictive" for e in result["evidence"])
    assert result["evidence"][0]["evidence_level"] == "A"


def test_get_predictive_evidence_gene_fallback(civic_db):
    """Variant 못 찾으면 gene-level fallback, match_level='gene'."""
    from scripts.db.query_civic import get_predictive_evidence_for_tier

    result = get_predictive_evidence_for_tier("BRAF", "p.Lys601Glu")
    assert result["match_level"] == "gene"


def test_get_predictive_evidence_no_match(civic_db):
    """CIViC에 없는 유전자는 match_level='none'."""
    from scripts.db.query_civic import get_predictive_evidence_for_tier

    result = get_predictive_evidence_for_tier("FAKEGENE", "p.Ala1Val")
    assert result["match_level"] == "none"
    assert result["evidence"] == []


def test_predictive_only_no_prognostic(civic_db):
    """Predictive evidence만 반환, Prognostic 제외."""
    from scripts.db.query_civic import get_predictive_evidence_for_tier

    result = get_predictive_evidence_for_tier("BRAF", "p.Val600Glu")
    for e in result["evidence"]:
        assert e["evidence_type"] == "Predictive"
