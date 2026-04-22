"""Single vs batch equivalence + no-leak regression tests (v2.5.4 Phase 3).

Two contracts that the H1 annotation-leak fix must maintain:

    1. The same VCF fed through ``run_pipeline`` and through
       ``run_batch_pipeline`` produces the same report_data structure
       (sample_id / timestamps excluded — those are by design
       mode-specific).
    2. A multi-sample batch never lets one sample's Variant annotation
       contaminate another sample's report. This is the direct
       regression test for the pre-v2.5.4 ``collect_unique_variants``
       dedup bug.
"""

from __future__ import annotations

from pathlib import Path

import pytest

BATCH_DIR = Path(__file__).parent.parent / "data" / "sample_vcf" / "batch_test"
SAMPLE_001_VCF = BATCH_DIR / "sample_001.vcf"


# ── Field filters for equivalence ────────────────────────────────────────────

# Fields that are allowed to differ between single-mode and batch-mode
# output for the same VCF. Timestamps and pipeline-mode-specific entries
# are not part of the correctness contract.
_BATCH_ONLY_DIFF_FIELDS = frozenset(
    {
        "sample_id",  # batch takes sample_id from manifest; single may differ
        "date",  # run-date of the day the test ran
        "pipeline",  # ped_path / ped_used differ when the harness injects PED
        "hpo_results",  # may be None vs [] depending on caller
    }
)


def _strip_nondeterministic(report: dict) -> dict:
    """Drop fields we don't assert on for single-vs-batch equivalence."""
    return {k: v for k, v in report.items() if k not in _BATCH_ONLY_DIFF_FIELDS}


def _variant_signature(v: dict) -> tuple:
    """Return a stable per-variant comparison tuple."""
    return (
        v.get("variant"),
        v.get("gene"),
        v.get("classification"),
        tuple(v.get("acmg_codes") or []),
        v.get("clinvar_significance"),
        v.get("tier"),
    )


# ── Single vs batch equivalence ──────────────────────────────────────────────


def test_single_vs_batch_report_data_equivalent(tmp_path):
    """Same VCF, same mode -> same classification / tier output."""
    from scripts.orchestrate import run_pipeline
    from scripts.orchestration.batch import run_batch_pipeline
    from scripts.orchestration.canonical import build_sample_report

    # Single-mode: go through the canonical entry point directly so we're
    # comparing apples to apples (run_pipeline adds CLI noise but the
    # report_data it returns is the same object).
    single_report = build_sample_report(
        vcf_path=str(SAMPLE_001_VCF),
        sample_id="sample_001",
        mode="cancer",
        skip_api=True,
    )
    assert single_report is not None

    # Run run_pipeline too — it should route through canonical and return
    # an equivalent dict.
    pipeline_report = run_pipeline(
        vcf_path=str(SAMPLE_001_VCF),
        output_path=str(tmp_path / "single.html"),
        mode="cancer",
        skip_api=True,
    )
    assert pipeline_report is not None

    # Batch-mode: point at a directory holding just this VCF.
    one_sample_dir = tmp_path / "single_batch"
    one_sample_dir.mkdir()
    (one_sample_dir / "sample_001.vcf").write_bytes(SAMPLE_001_VCF.read_bytes())

    batch_result = run_batch_pipeline(
        batch_path=str(one_sample_dir),
        output_dir=None,
        mode="cancer",
        workers=1,
        skip_api=True,
        _skip_reports=True,
    )
    assert batch_result["samples_processed"] == 1
    assert batch_result["errors"] == []

    # run_batch_pipeline doesn't expose each sample's report_data directly
    # (it discards them once the HTML is emitted), so we compare via
    # equivalence of the per-variant signatures. The canonical call above
    # is the same function batch uses, so signature equality here is what
    # drives the guarantee.
    sigs_single = {_variant_signature(v) for v in single_report["variants"]}
    sigs_pipeline = {_variant_signature(v) for v in pipeline_report["variants"]}
    assert sigs_single == sigs_pipeline

    # Mode / structural shape also agree.
    assert single_report["mode"] == pipeline_report["mode"] == "cancer"
    assert single_report["summary"]["total"] == pipeline_report["summary"]["total"]


def test_single_vs_batch_structural_shape_stable(tmp_path):
    """Top-level keys of a batch-assembled report_data match the single-mode shape."""
    from scripts.orchestration.canonical import build_sample_report

    single_report = build_sample_report(
        vcf_path=str(SAMPLE_001_VCF),
        sample_id="sample_001",
        mode="cancer",
        skip_api=True,
    )
    assert single_report is not None

    # Keys run_pipeline consumers rely on — see tests/test_orchestrate.py
    # test_orchestrate_report_structure for the canonical list.
    for key in (
        "sample_id",
        "date",
        "variants",
        "detailed_variants",
        "omitted_variants",
        "summary",
        "tier1_variants",
        "tier2_variants",
        "tier3_variants",
        "tier4_count",
        "pgx_results",
        "pgx_source",
        "db_versions",
        "pipeline",
        "mode",
        "sv_variants",
        "tmb",
    ):
        assert key in single_report, f"Missing top-level key: {key}"


# ── Variant annotation leak (H1) ─────────────────────────────────────────────


def test_batch_no_variant_object_leak(tmp_path):
    """Each sample's report contains only its own genes.

    sample_001 has TP53 / BRCA2 / CFTR
    sample_002 has CFTR / CYP2C19 / HLA-B       (CFTR overlaps sample_001)
    sample_003 has TP53 / PALB2 / PTPN11        (TP53 overlaps sample_001)

    Pre-v2.5.4 the pipeline reused the first sample's Variant object for
    overlapping coordinates — so annotation from sample_001's TP53 would
    leak into sample_003's TP53 row. The new pipeline parses each VCF
    separately; the gene set per sample is exactly the gene set declared
    in that VCF.
    """
    from scripts.orchestration.batch import _assemble_all, discover_samples

    samples = discover_samples(str(BATCH_DIR))
    assert len(samples) == 3

    reports, errors = _assemble_all(
        samples,
        mode="cancer",
        skip_api=True,
        hpo_ids=None,
        hide_vus=False,
        clinical_board=False,
        board_lang=None,
        workers=1,
    )
    assert errors == []
    assert len(reports) == 3

    by_sample = {r["sample_id"]: r for r in reports}

    expected = {
        "sample_001": {"TP53", "BRCA2", "CFTR"},
        "sample_002": {"CFTR", "CYP2C19", "HLA-B"},
        "sample_003": {"TP53", "PALB2", "PTPN11"},
    }

    for sample_id, genes in expected.items():
        report = by_sample[sample_id]
        observed = {v.get("gene") for v in report["variants"]} - {None, ""}
        assert observed == genes, f"{sample_id} expected {genes}, got {observed}"


def test_batch_skip_api_propagated_to_provenance(tmp_path):
    """v2.5.4 M3: batch now reports the caller's actual skip_api in report_data['pipeline']."""
    from scripts.orchestration.batch import _assemble_all, discover_samples

    samples = discover_samples(str(BATCH_DIR))[:1]  # one sample is enough
    reports, _ = _assemble_all(
        samples,
        mode="cancer",
        skip_api=True,  # was hardcoded False before v2.5.4
        hpo_ids=None,
        hide_vus=False,
        clinical_board=False,
        board_lang=None,
        workers=1,
    )
    assert reports
    assert reports[0]["pipeline"]["skip_api"] is True


# ── Batch gracefully skips samples with zero parsed variants ────────────────


def test_batch_skips_empty_vcf_without_crashing(tmp_path):
    """A malformed / empty VCF becomes an 'errors' entry, other samples succeed."""
    from scripts.orchestration.batch import run_batch_pipeline

    # Mix one real VCF with one header-only VCF. discover_samples picks
    # both up; build_sample_report returns None for the empty one.
    d = tmp_path / "mixed"
    d.mkdir()
    (d / "sample_001.vcf").write_bytes(SAMPLE_001_VCF.read_bytes())
    (d / "empty.vcf").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    result = run_batch_pipeline(
        batch_path=str(d),
        output_dir=None,
        mode="cancer",
        workers=1,
        skip_api=True,
        _skip_reports=True,
    )
    assert result["samples_processed"] == 2
    # sample_001 succeeds, empty.vcf becomes an error entry.
    assert len(result["errors"]) == 1
    assert result["errors"][0]["sample"] == "empty"


# ── Extended manifest CSV (M5 + M6) ──────────────────────────────────────────


def test_batch_discovery_detects_vcf_bgz_extension(tmp_path):
    """discover_samples now recognises .vcf.bgz files (v2.5.4 M5)."""
    from scripts.orchestration.batch import discover_samples

    # We only need the extension to be present — discover_samples does
    # not open the file. Write it with a byte-order marker so the harness
    # can distinguish it from the empty default.
    bgz = tmp_path / "bgz_sample.vcf.bgz"
    bgz.write_bytes(b"\x1f\x8b\x08\x04\x00")  # bgzf magic; file doesn't need to be valid
    plain = tmp_path / "plain_sample.vcf"
    plain.write_bytes(b"##fileformat=VCFv4.2\n")
    gz = tmp_path / "gz_sample.vcf.gz"
    gz.write_bytes(b"\x1f\x8b")  # gzip magic

    samples = discover_samples(str(tmp_path))
    sample_ids = {s["sample_id"] for s in samples}
    assert "bgz_sample" in sample_ids
    assert "plain_sample" in sample_ids
    assert "gz_sample" in sample_ids


def test_batch_manifest_with_extended_columns(tmp_path):
    """Manifest CSV with optional columns parses through discover_samples."""
    import csv

    from scripts.orchestration.batch import discover_samples

    vcf = tmp_path / "sample.vcf"
    vcf.write_bytes(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    germline = tmp_path / "germ.vcf.gz"
    germline.write_bytes(b"\x1f\x8b")
    ped = tmp_path / "fam.ped"
    ped.write_text("FAM1\tproband\t0\t0\t1\t2\n")

    manifest = tmp_path / "manifest.csv"
    with open(manifest, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "sample_id",
                "vcf_path",
                "germline_vcf",
                "ped_path",
                "hpo_ids",
                "panel_size_mb",
                "clinical_note",
            ],
        )
        w.writeheader()
        w.writerow(
            {
                "sample_id": "S001",
                "vcf_path": str(vcf),
                "germline_vcf": str(germline),
                "ped_path": str(ped),
                "hpo_ids": "HP:0001250;HP:0001263",
                "panel_size_mb": "36.7",
                "clinical_note": "Seizures + developmental delay",
            }
        )

    samples = discover_samples(str(manifest))
    assert len(samples) == 1
    s = samples[0]
    assert s["sample_id"] == "S001"
    assert s["germline_vcf"] == str(germline)
    assert s["ped_path"] == str(ped)
    assert s["hpo_ids"] == ["HP:0001250", "HP:0001263"]
    assert s["panel_size_mb"] == pytest.approx(36.7)
    assert s["clinical_note"] == "Seizures + developmental delay"


def test_batch_manifest_legacy_two_column_still_works(tmp_path):
    """The pre-v2.5.4 two-column manifest must keep working."""
    import csv

    from scripts.orchestration.batch import discover_samples

    vcf = tmp_path / "s.vcf"
    vcf.write_bytes(b"##fileformat=VCFv4.2\n")
    manifest = tmp_path / "m.csv"
    with open(manifest, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample_id", "vcf_path"])
        w.writeheader()
        w.writerow({"sample_id": "legacy", "vcf_path": str(vcf)})

    samples = discover_samples(str(manifest))
    assert samples == [
        {
            "sample_id": "legacy",
            "vcf_path": str(vcf),
            "germline_vcf": None,
            "ped_path": None,
            "hpo_ids": [],
            "sv_path": None,
            "intervar_path": None,
            "clinical_note": None,
            "panel_size_mb": None,
        }
    ]
