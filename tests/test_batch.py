# tests/test_batch.py
"""Tests for batch pipeline: discovery, deduplication, and parallel report generation."""

import csv
from pathlib import Path

import pytest

BATCH_DIR = Path(__file__).parent.parent / "data" / "sample_vcf" / "batch_test"


# ── Discovery tests ────────────────────────────────────────────────────────────


def test_discover_samples_from_directory():
    """_discover_samples finds all VCFs in batch_test directory."""
    from scripts.orchestrate import _discover_samples

    samples = _discover_samples(str(BATCH_DIR))

    assert len(samples) == 3
    sample_ids = {s["sample_id"] for s in samples}
    assert sample_ids == {"sample_001", "sample_002", "sample_003"}
    for s in samples:
        assert Path(s["vcf_path"]).exists(), f"VCF not found: {s['vcf_path']}"


def test_discover_samples_from_csv(tmp_path):
    """_discover_samples reads a manifest CSV with sample_id and vcf_path columns."""
    from scripts.orchestrate import _discover_samples

    manifest = tmp_path / "manifest.csv"
    rows = [
        {"sample_id": "S001", "vcf_path": str(BATCH_DIR / "sample_001.vcf")},
        {"sample_id": "S002", "vcf_path": str(BATCH_DIR / "sample_002.vcf")},
    ]
    with open(manifest, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_id", "vcf_path"])
        writer.writeheader()
        writer.writerows(rows)

    samples = _discover_samples(str(manifest))

    assert len(samples) == 2
    assert samples[0]["sample_id"] == "S001"
    assert samples[1]["sample_id"] == "S002"
    assert samples[0]["vcf_path"] == str(BATCH_DIR / "sample_001.vcf")


def test_discover_samples_invalid_path():
    """_discover_samples raises ValueError for a non-directory, non-CSV path."""
    from scripts.orchestrate import _discover_samples

    with pytest.raises(ValueError, match="Batch path must be"):
        _discover_samples("/tmp/nonexistent_file.txt")


# ── Deduplication tests ────────────────────────────────────────────────────────


def test_collect_unique_variants():
    """_collect_unique_variants correctly deduplicates shared variants.

    sample_001 has TP53, BRCA2, CFTR
    sample_002 has CFTR, CYP2C19, HLA-B      (CFTR shared with 001)
    sample_003 has TP53, PALB2, PTPN11       (TP53 shared with 001)

    Total = 9 variant references, Unique = 7 (TP53 and CFTR appear twice).
    """
    from scripts.orchestrate import _collect_unique_variants, _discover_samples

    samples = _discover_samples(str(BATCH_DIR))
    unique_variants, sample_map = _collect_unique_variants(samples)

    total = sum(len(keys) for keys in sample_map.values())

    # 3 samples × 3 variants each = 9 total references
    assert total == 9
    # TP53 shared between 001 and 003, CFTR shared between 001 and 002 → 7 unique
    assert len(unique_variants) == 7
    assert total > len(unique_variants), "Expected deduplication (total > unique)"


def test_batch_pipeline_dedup_count():
    """run_batch_pipeline reports unique_variants < total_variants due to shared TP53/CFTR."""
    from scripts.orchestrate import run_batch_pipeline

    result = run_batch_pipeline(
        batch_path=str(BATCH_DIR),
        output_dir=None,  # skip report generation
        mode="cancer",
        workers=1,
        skip_api=True,
        _skip_reports=True,
    )

    assert result["total_variants"] == 9
    assert result["unique_variants"] == 7
    assert result["cache_hits"] == 2  # TP53 + CFTR each seen a second time


# ── End-to-end tests ───────────────────────────────────────────────────────────


def test_batch_pipeline_end_to_end(tmp_path):
    """run_batch_pipeline processes 3 VCFs and generates 3 HTML reports."""
    from scripts.orchestrate import run_batch_pipeline

    out_dir = str(tmp_path / "batch_out")

    result = run_batch_pipeline(
        batch_path=str(BATCH_DIR),
        output_dir=out_dir,
        mode="cancer",
        workers=2,
        skip_api=True,
    )

    assert result["samples_processed"] == 3
    assert len(result["reports_generated"]) == 3
    assert result["errors"] == []

    # All reports must actually exist and contain BIKO GenomeBoard HTML
    for report_path in result["reports_generated"]:
        p = Path(report_path)
        assert p.exists(), f"Report not found: {report_path}"
        assert "BIKO GenomeBoard" in p.read_text(encoding="utf-8")


def test_batch_output_dir_created(tmp_path):
    """run_batch_pipeline creates the output directory if it does not exist."""
    from scripts.orchestrate import run_batch_pipeline

    out_dir = tmp_path / "deep" / "nested" / "batch_out"
    assert not out_dir.exists()

    run_batch_pipeline(
        batch_path=str(BATCH_DIR),
        output_dir=str(out_dir),
        mode="cancer",
        workers=1,
        skip_api=True,
    )

    assert out_dir.exists()


def test_batch_pipeline_csv_manifest(tmp_path):
    """run_batch_pipeline works with a CSV manifest instead of a directory."""
    from scripts.orchestrate import run_batch_pipeline

    manifest = tmp_path / "manifest.csv"
    with open(manifest, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_id", "vcf_path"])
        writer.writeheader()
        writer.writerow({"sample_id": "S001", "vcf_path": str(BATCH_DIR / "sample_001.vcf")})
        writer.writerow({"sample_id": "S002", "vcf_path": str(BATCH_DIR / "sample_002.vcf")})

    out_dir = str(tmp_path / "out")
    result = run_batch_pipeline(
        batch_path=str(manifest),
        output_dir=out_dir,
        mode="cancer",
        workers=1,
        skip_api=True,
    )

    assert result["samples_processed"] == 2
    assert len(result["reports_generated"]) == 2
    assert result["errors"] == []


def test_batch_result_schema():
    """run_batch_pipeline return dict contains all required keys."""
    from scripts.orchestrate import run_batch_pipeline

    result = run_batch_pipeline(
        batch_path=str(BATCH_DIR),
        output_dir=None,
        mode="cancer",
        workers=1,
        skip_api=True,
        _skip_reports=True,
    )

    for key in (
        "samples_processed",
        "total_variants",
        "unique_variants",
        "cache_hits",
        "reports_generated",
        "errors",
        "elapsed_seconds",
    ):
        assert key in result, f"Missing key in batch result: {key}"
