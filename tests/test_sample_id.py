"""sample_id normalization — single mode and batch mode must agree.

v2.5.4 M4 fix coverage. Before the fix:
    * single mode derived sample_id via Path(vcf_path).stem.upper()
    * batch mode stripped .vcf / .vcf.gz manually and kept the original
      case.

These disagreed for ``patient_001.vcf`` (single: ``PATIENT_001``;
batch: ``patient_001``) — which broke regression suites that compared
single-mode output to batch-mode output. The new helper
``scripts.orchestration.canonical.normalize_sample_id`` is the single
source of truth; both modes import it.
"""

from __future__ import annotations

import pytest


# ── normalize_sample_id unit tests ───────────────────────────────────────────


def test_normalize_strips_plain_vcf():
    from scripts.orchestration.canonical import normalize_sample_id

    assert normalize_sample_id("/tmp/patient.vcf") == "patient"


def test_normalize_strips_bgzip_vcf():
    from scripts.orchestration.canonical import normalize_sample_id

    assert normalize_sample_id("/tmp/patient.vcf.gz") == "patient"


def test_normalize_strips_bgz_vcf():
    from scripts.orchestration.canonical import normalize_sample_id

    assert normalize_sample_id("/tmp/patient.vcf.bgz") == "patient"


def test_normalize_preserves_case():
    from scripts.orchestration.canonical import normalize_sample_id

    assert normalize_sample_id("/tmp/Patient_ABC_123.vcf") == "Patient_ABC_123"


def test_normalize_preserves_case_bgzip():
    from scripts.orchestration.canonical import normalize_sample_id

    assert normalize_sample_id("/tmp/MixedCASE.vcf.gz") == "MixedCASE"


def test_normalize_leaves_unknown_extension_alone():
    """A filename without a recognised VCF extension is returned as-is."""
    from scripts.orchestration.canonical import normalize_sample_id

    assert normalize_sample_id("/tmp/sample.bam") == "sample.bam"


def test_normalize_handles_bare_filename():
    from scripts.orchestration.canonical import normalize_sample_id

    assert normalize_sample_id("patient_001.vcf") == "patient_001"


def test_normalize_prefers_longer_extension():
    """``.vcf.bgz`` must take precedence over ``.vcf``/``.vcf.gz``."""
    from scripts.orchestration.canonical import normalize_sample_id

    # All three endings appear but only .vcf.bgz should be stripped.
    assert normalize_sample_id("my.vcf.bgz") == "my"


# ── Cross-mode consistency (single vs batch) ─────────────────────────────────


@pytest.mark.parametrize(
    "filename, expected",
    [
        ("patient_001.vcf", "patient_001"),
        ("Patient_XYZ.vcf.gz", "Patient_XYZ"),
        ("sample.vcf.bgz", "sample"),
        ("DEMO-KR-001.vcf", "DEMO-KR-001"),
    ],
)
def test_single_and_batch_agree_on_sample_id(filename: str, expected: str, tmp_path):
    """run_pipeline (default sample_id) and batch discovery produce the same id."""
    from scripts.orchestration.batch import discover_samples
    from scripts.orchestration.canonical import normalize_sample_id

    # Write a minimal but valid-looking VCF stub so discover_samples picks
    # it up. discover_samples does not parse the file — it just returns
    # the derived sample_id — so content does not matter here.
    vcf = tmp_path / filename
    vcf.write_bytes(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # Single-mode path — the default code inside run_pipeline is exactly
    # the normalize_sample_id helper.
    single_id = normalize_sample_id(str(vcf))
    assert single_id == expected

    # Batch-mode path — discover_samples writes the derived id.
    samples = discover_samples(str(tmp_path))
    matching = [s for s in samples if s["vcf_path"].endswith(filename)]
    assert len(matching) == 1
    assert matching[0]["sample_id"] == expected


def test_sample_id_no_longer_uppercased():
    """v2.5.4 breaking change: .stem.upper() is gone."""
    from scripts.orchestration.canonical import normalize_sample_id

    # Pre-v2.5.4 this returned "SAMPLE.VCF" (.stem.upper() on "sample.vcf.gz"
    # left the inner ".vcf" and uppercased everything). The new helper
    # returns "sample" — lowercase preserved, extension fully stripped.
    assert normalize_sample_id("sample.vcf.gz") == "sample"
    assert normalize_sample_id("sample.vcf.gz") != "SAMPLE.VCF"
