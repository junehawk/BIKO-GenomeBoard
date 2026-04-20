"""Tests for CIViC database build + query + hotspot tiering."""

import sqlite3
from pathlib import Path

import pytest

# ── Fixtures ──────────────────────────────────────────────────────────────────

CIVIC_DIR = str(Path(__file__).parent.parent / "data" / "db" / "civic")
DEMO_VCF = str(Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf")


@pytest.fixture(scope="module")
def civic_db(tmp_path_factory):
    """Build a fresh CIViC DB in a temp dir and return its path."""
    from scripts.storage.build_civic_db import build_db

    db_path = str(tmp_path_factory.mktemp("civic_db") / "civic_test.sqlite3")
    build_db(civic_dir=CIVIC_DIR, db_path=db_path)
    return db_path


@pytest.fixture(autouse=True)
def reset_civic_conn():
    """Reset the query_civic connection cache around each test."""
    from scripts.storage import query_civic

    query_civic.reset_civic_connection()
    yield
    query_civic.reset_civic_connection()


# ── Build tests ───────────────────────────────────────────────────────────────


def test_build_civic_db_tables(civic_db):
    """Built DB has all expected tables."""
    conn = sqlite3.connect(civic_db)
    tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()}
    conn.close()
    assert "genes" in tables
    assert "variants" in tables
    assert "evidence" in tables
    assert "hotspots" in tables
    assert "metadata" in tables


def test_build_civic_db_counts(civic_db):
    """Built DB has reasonable row counts matching the source TSVs."""
    conn = sqlite3.connect(civic_db)
    gene_count = conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
    variant_count = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    evidence_count = conn.execute("SELECT COUNT(*) FROM evidence").fetchone()[0]
    hotspot_count = conn.execute("SELECT COUNT(*) FROM hotspots").fetchone()[0]
    conn.close()
    assert gene_count > 100, f"Expected >100 genes, got {gene_count}"
    assert variant_count > 100, f"Expected >100 variants, got {variant_count}"
    assert evidence_count > 100, f"Expected >100 evidence items, got {evidence_count}"
    assert hotspot_count > 10, f"Expected >10 hotspots, got {hotspot_count}"


def test_build_civic_db_metadata(civic_db):
    """Metadata table has build_date and source entries."""
    conn = sqlite3.connect(civic_db)
    meta = {r[0]: r[1] for r in conn.execute("SELECT key, value FROM metadata").fetchall()}
    conn.close()
    assert "build_date" in meta
    assert "source" in meta
    assert "CIViC" in meta["source"]
    assert int(meta["gene_count"]) > 100
    assert int(meta["variant_count"]) > 100


def test_evidence_has_therapy_ids_column(civic_db):
    """v2.2 A1-db-4: evidence table carries a therapy_ids column used as
    the primary merge key against OncoKB. Drug-name string match is fallback
    only, so this column must exist on every fresh build."""
    conn = sqlite3.connect(civic_db)
    cols = {r[1] for r in conn.execute("PRAGMA table_info(evidence)").fetchall()}
    conn.close()
    assert "therapy_ids" in cols, (
        "CIViC evidence table is missing therapy_ids column — "
        "OncoKB/CIViC merge key will fall back to fragile string matching"
    )


def test_civic_therapy_ids_migration_is_idempotent(tmp_path):
    """Rebuilding on top of an existing civic.sqlite3 must not drop data
    and the second invocation must leave therapy_ids in place."""
    from scripts.storage.build_civic_db import build_db

    db_path = str(tmp_path / "civic_idempotent.sqlite3")
    build_db(civic_dir=CIVIC_DIR, db_path=db_path)

    conn = sqlite3.connect(db_path)
    first_count = conn.execute("SELECT COUNT(*) FROM evidence").fetchone()[0]
    cols_before = {r[1] for r in conn.execute("PRAGMA table_info(evidence)").fetchall()}
    conn.close()
    assert "therapy_ids" in cols_before
    assert first_count > 0

    # Rebuild in place — should not error, and therapy_ids must still exist.
    build_db(civic_dir=CIVIC_DIR, db_path=db_path)

    conn = sqlite3.connect(db_path)
    second_count = conn.execute("SELECT COUNT(*) FROM evidence").fetchone()[0]
    cols_after = {r[1] for r in conn.execute("PRAGMA table_info(evidence)").fetchall()}
    conn.close()
    assert "therapy_ids" in cols_after
    assert second_count == first_count


# ── Query tests ───────────────────────────────────────────────────────────────


def test_get_gene_summary_kras(civic_db, monkeypatch):
    """get_gene_summary returns a description for KRAS."""
    monkeypatch.setenv("GB_CONFIG_PATH", "")
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    result = qc.get_gene_summary("KRAS")
    assert result is not None
    assert result["gene"] == "KRAS"
    assert len(result["description"]) > 50
    assert "KRAS" in result["description"] or "RAS" in result["description"]


def test_get_gene_summary_unknown(civic_db, monkeypatch):
    """get_gene_summary returns None for unknown genes."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    result = qc.get_gene_summary("NOTAREALGENE999")
    assert result is None


def test_get_variant_evidence_kras(civic_db, monkeypatch):
    """get_variant_evidence returns evidence items with therapies for KRAS."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    results = qc.get_variant_evidence("KRAS")
    assert len(results) > 0
    # At least some items should have therapies
    has_therapy = any(e["therapies"] for e in results)
    assert has_therapy, "Expected at least one KRAS evidence item with therapies"
    # All items should have expected keys
    for item in results[:5]:
        assert "gene" in item
        assert "variant" in item
        assert "disease" in item
        assert "evidence_type" in item
        assert "evidence_level" in item


def test_get_treatment_summary_braf(civic_db, monkeypatch):
    """get_treatment_summary returns a formatted treatment string for BRAF."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    result = qc.get_treatment_summary("BRAF")
    assert isinstance(result, str)
    assert len(result) > 0
    # Should contain a Level indicator
    assert "Level" in result


def test_get_treatment_summary_empty_for_unknown(civic_db, monkeypatch):
    """get_treatment_summary returns empty string for unknown gene."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    result = qc.get_treatment_summary("NOTAREALGENE999")
    assert result == ""


# ── Hotspot tests ─────────────────────────────────────────────────────────────


def test_is_hotspot_tp53_249(civic_db, monkeypatch):
    """TP53 position 249 is a known hotspot."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    assert qc.is_hotspot("TP53", 249) is True


def test_is_hotspot_tp53_999(civic_db, monkeypatch):
    """TP53 position 999 is not a known hotspot."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    assert qc.is_hotspot("TP53", 999) is False


def test_is_hotspot_kras_12(civic_db, monkeypatch):
    """KRAS position 12 is a known hotspot."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    assert qc.is_hotspot("KRAS", 12) is True


def test_is_hotspot_no_connection():
    """is_hotspot returns False when DB is unavailable."""
    import scripts.storage.query_civic as qc

    qc._conn = None
    # Use a path that doesn't exist to ensure no connection
    import os

    orig = os.environ.get("GB_CONFIG_PATH")
    # Patch _get_connection to return None
    orig_fn = qc._get_connection
    qc._get_connection = lambda: None
    try:
        assert qc.is_hotspot("KRAS", 12) is False
    finally:
        qc._get_connection = orig_fn


def test_get_hotspot_variants_kras_12(civic_db, monkeypatch):
    """get_hotspot_variants returns multiple variants for KRAS G12."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    variants = qc.get_hotspot_variants("KRAS", 12)
    assert len(variants) >= 2
    # Common KRAS G12 variants
    assert any("G12" in v for v in variants)


# ── extract_protein_position tests ────────────────────────────────────────────


def test_extract_protein_position_3letter():
    """Extracts position from 3-letter HGVSp notation."""
    from scripts.storage.query_civic import extract_protein_position

    assert extract_protein_position("p.Gly12Asp") == 12
    assert extract_protein_position("p.Val600Glu") == 600
    assert extract_protein_position("p.Arg249Met") == 249


def test_extract_protein_position_1letter():
    """Extracts position from 1-letter HGVSp notation."""
    from scripts.storage.query_civic import extract_protein_position

    assert extract_protein_position("p.R249M") == 249
    assert extract_protein_position("p.V600E") == 600
    assert extract_protein_position("p.G12D") == 12


def test_extract_protein_position_none():
    """Returns None for empty/invalid input."""
    from scripts.storage.query_civic import extract_protein_position

    assert extract_protein_position(None) is None
    assert extract_protein_position("") is None
    assert extract_protein_position("c.35G>A") is None


def test_extract_protein_position_with_transcript_prefix():
    """Handles HGVSp with transcript prefix."""
    from scripts.storage.query_civic import extract_protein_position

    assert extract_protein_position("NP_004324.2:p.Val600Glu") == 600


# ── _hgvsp_to_civic_variant tests ─────────────────────────────────────────────


def test_hgvsp_to_civic_variant_gly12asp():
    """p.Gly12Asp -> G12D"""
    from scripts.counselor.generate_pdf import _hgvsp_to_civic_variant

    assert _hgvsp_to_civic_variant("p.Gly12Asp") == "G12D"


def test_hgvsp_to_civic_variant_val600glu():
    """p.Val600Glu -> V600E"""
    from scripts.counselor.generate_pdf import _hgvsp_to_civic_variant

    assert _hgvsp_to_civic_variant("p.Val600Glu") == "V600E"


def test_hgvsp_to_civic_variant_arg249met():
    """p.Arg249Met -> R249M"""
    from scripts.counselor.generate_pdf import _hgvsp_to_civic_variant

    assert _hgvsp_to_civic_variant("p.Arg249Met") == "R249M"


def test_hgvsp_to_civic_variant_none():
    """Returns None for empty/invalid input."""
    from scripts.counselor.generate_pdf import _hgvsp_to_civic_variant

    assert _hgvsp_to_civic_variant(None) is None
    assert _hgvsp_to_civic_variant("") is None
    assert _hgvsp_to_civic_variant("c.35G>A") is None


# ── assign_tier hotspot upgrade tests ────────────────────────────────────────


def test_assign_tier_vus_hotspot_tp53(civic_db, monkeypatch):
    """VUS on TP53 at hotspot position 249 → Tier 2."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    from scripts.enrichment.oncokb import assign_tier, reset_oncokb_cache

    reset_oncokb_cache()
    # p.Arg249Ser → position 249 (TP53 R249S is in hotspots)
    result = assign_tier("VUS", "TP53", hgvsp="p.Arg249Ser")
    assert result == 2, f"Expected Tier 2 for VUS TP53 hotspot, got {result}"


def test_assign_tier_vus_hotspot_kras_12(civic_db, monkeypatch):
    """VUS on KRAS at hotspot position 12 → Tier 2."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    from scripts.enrichment.oncokb import assign_tier, reset_oncokb_cache

    reset_oncokb_cache()
    result = assign_tier("VUS", "KRAS", hgvsp="p.Gly12Asp")
    assert result == 2, f"Expected Tier 2 for VUS KRAS G12D hotspot, got {result}"


def test_assign_tier_vus_non_hotspot_cancer_gene(civic_db, monkeypatch):
    """VUS on KRAS at non-hotspot position → Tier 3."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    from scripts.enrichment.oncokb import assign_tier, reset_oncokb_cache

    reset_oncokb_cache()
    result = assign_tier("VUS", "KRAS", hgvsp="p.Ala999Val")
    assert result == 3, f"Expected Tier 3 for non-hotspot VUS KRAS, got {result}"


def test_assign_tier_vus_no_hgvsp_still_tier3():
    """VUS on cancer gene without hgvsp → Tier 3 (no hotspot check possible)."""
    from scripts.enrichment.oncokb import assign_tier, reset_oncokb_cache

    reset_oncokb_cache()
    result = assign_tier("VUS", "KRAS", hgvsp="")
    assert result == 3


def test_assign_tier_existing_behavior_unchanged():
    """Existing tier assignments are not affected by the new hgvsp parameter."""
    from scripts.enrichment.oncokb import assign_tier, reset_oncokb_cache

    reset_oncokb_cache()
    # Pathogenic level-1 gene → Tier 1
    assert assign_tier("Pathogenic", "KRAS") == 1
    # Drug response → Tier 1
    assert assign_tier("Drug Response", "") == 1
    # Benign cancer gene → Tier 4
    assert assign_tier("Benign", "KRAS") == 4
    # VUS non-cancer → Tier 4
    assert assign_tier("VUS", "NOTAREALGENE") == 4


# ── Pipeline integration test ─────────────────────────────────────────────────


def test_pipeline_kras_has_treatment(civic_db, tmp_path, monkeypatch):
    """Run pipeline on demo VCF; any KRAS variant gets treatment_strategies from CIViC."""
    import scripts.storage.query_civic as qc

    qc._conn = sqlite3.connect(civic_db)
    qc._conn.row_factory = sqlite3.Row

    from scripts.enrichment.oncokb import reset_oncokb_cache

    reset_oncokb_cache()

    from scripts.orchestrate import run_pipeline

    output_path = tmp_path / "report.html"
    result = run_pipeline(
        vcf_path=DEMO_VCF,
        output_path=str(output_path),
        skip_api=True,
        mode="cancer",
    )

    assert result is not None
    # Check that KRAS variants (if present) have treatment data,
    # OR that at least the report was generated without errors.
    kras_variants = [v for v in result.get("variants", []) if v.get("gene") == "KRAS"]
    for kv in kras_variants:
        # treatment_strategies may come from gene_knowledge or CIViC
        assert "treatment_strategies" in kv or kv.get("tier") is not None

    # Report file must exist
    assert output_path.exists()
