"""Tests for scripts/db/build_pm1_hotspots.py — v2.2 A3.

Exercises the three-layer merge (cancerhotspots + ClinGen baseline + overrides),
JSON schema invariants, idempotency, and exclusions handling.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from scripts.db import build_pm1_hotspots as bph


PMID_RE = re.compile(r"PMID\s+\d+")
ALLOWED_STRENGTHS = {"moderate", "supporting"}


# ── Schema / invariant tests ──────────────────────────────────────────────────


def test_default_build_against_repo_overrides(tmp_path):
    """Building with the canonical data/pm1_hotspot_overrides.yaml yields a
    schema-conformant JSON payload that covers the v2.2 minimum gene set."""
    output = tmp_path / "pm1.json"
    payload = bph.build_pm1_hotspots(output_path=output)

    assert output.exists()
    on_disk = json.loads(output.read_text())
    assert on_disk == payload

    # Top-level schema fields
    assert payload["$schema"] == "https://biko/schema/pm1_hotspot_domains.v1.json"
    assert payload["version"]
    assert payload["build_date"]
    assert payload["source_refs"]
    assert payload["source_hash"].startswith("sha256:")
    assert payload["record_count"] == sum(len(v) for v in payload["genes"].values())

    # v2.2 minimum gene set
    required_genes = {"TP53", "KRAS", "NRAS", "HRAS", "BRAF", "PIK3CA", "EGFR", "IDH1", "IDH2"}
    assert required_genes.issubset(payload["genes"].keys())


def test_every_entry_has_pmid_source(tmp_path):
    """A3-db-4 invariant: every hotspot entry must carry a PMID-qualified
    source. Entries without a PMID reference are rejected by _validate."""
    output = tmp_path / "pm1.json"
    payload = bph.build_pm1_hotspots(output_path=output)

    violations: list[str] = []
    for gene, entries in payload["genes"].items():
        for entry in entries:
            source = entry.get("source", "")
            if not PMID_RE.search(source):
                violations.append(f"{gene} {entry.get('range')}: {source!r}")
    assert not violations, f"Entries without PMID source: {violations}"


def test_every_entry_valid_strength_and_range(tmp_path):
    """Strength must be moderate|supporting; range must be [start, end] with
    start <= end. Both invariants are enforced by _validate but we also
    assert them here so a direct YAML edit cannot sneak malformed rows past."""
    output = tmp_path / "pm1.json"
    payload = bph.build_pm1_hotspots(output_path=output)

    for gene, entries in payload["genes"].items():
        for entry in entries:
            assert entry["strength"] in ALLOWED_STRENGTHS, (
                f"{gene} {entry['range']}: strength={entry['strength']!r}"
            )
            rng = entry["range"]
            assert isinstance(rng, list) and len(rng) == 2
            assert isinstance(rng[0], int) and isinstance(rng[1], int)
            assert rng[0] <= rng[1]


def test_tp53_block_moderate_and_locked(tmp_path):
    """Clinical lock: TP53 245-249 must be moderate and locked=true.
    This is the load-bearing override from clinical-advisor §1.3."""
    output = tmp_path / "pm1.json"
    payload = bph.build_pm1_hotspots(output_path=output)

    tp53_entries = payload["genes"]["TP53"]
    block = next(
        (e for e in tp53_entries if e["range"] == [245, 249]),
        None,
    )
    assert block is not None, "TP53 245-249 block missing"
    assert block["strength"] == "moderate"
    assert block["locked"] is True


def test_deterministic_source_hash_over_runs(tmp_path):
    """Two builds with the same inputs must yield the same source_hash.
    Ensures the hash is stable for version-manager drift detection."""
    out_a = tmp_path / "a.json"
    out_b = tmp_path / "b.json"
    payload_a = bph.build_pm1_hotspots(output_path=out_a)
    payload_b = bph.build_pm1_hotspots(output_path=out_b)
    assert payload_a["source_hash"] == payload_b["source_hash"]
    assert payload_a["genes"] == payload_b["genes"]


def test_idempotent_rebuild_in_place(tmp_path):
    """Rebuilding over an existing output path preserves record_count."""
    output = tmp_path / "pm1.json"
    first = bph.build_pm1_hotspots(output_path=output)
    second = bph.build_pm1_hotspots(output_path=output)
    assert first["record_count"] == second["record_count"] > 0


# ── Layer-level merge tests ───────────────────────────────────────────────────


def test_cancerhotspots_stub_is_skipped(tmp_path):
    """A 404 stub for cancerhotspots_v2_single.tsv must not raise; the
    build proceeds using baseline + overrides."""
    stub = tmp_path / "stub.tsv"
    stub.write_text('{"status":404,"error":"Not Found"}')
    output = tmp_path / "pm1.json"
    payload = bph.build_pm1_hotspots(
        cancerhotspots_tsv=stub,
        output_path=output,
    )
    assert payload["record_count"] > 0
    assert "TP53" in payload["genes"]


def test_missing_cancerhotspots_file_is_skipped(tmp_path):
    """An entirely absent cancerhotspots path is treated as an empty layer."""
    missing = tmp_path / "does_not_exist.tsv"
    output = tmp_path / "pm1.json"
    payload = bph.build_pm1_hotspots(
        cancerhotspots_tsv=missing,
        output_path=output,
    )
    assert payload["record_count"] > 0


def test_overrides_exclusions_remove_entries(tmp_path):
    """When the overrides YAML lists an exclusion, the merged output must
    drop the matching entry even if a higher layer introduced it.

    Uses an add-then-exclude dance on a non-required residue so the build
    still satisfies REQUIRED_COVERAGE.
    """
    custom_overrides = tmp_path / "overrides.yaml"
    custom_overrides.write_text(
        """
version: "1.0-test"
additions:
  TP53:
    - range: [999, 999]
      domain: "throwaway test entry"
      strength: supporting
      source: "PMID 1"
exclusions:
  TP53:
    - range: [999, 999]
      reason: "test exclusion"
"""
    )
    output = tmp_path / "pm1.json"
    payload = bph.build_pm1_hotspots(
        overrides_path=custom_overrides,
        output_path=output,
    )
    tp53 = payload["genes"].get("TP53", [])
    assert not any(e["range"] == [999, 999] for e in tp53), (
        "TP53 999 should be excluded by the overrides YAML"
    )
    # Baseline TP53 entries (175, 245-249, 273, 282) must remain untouched.
    assert any(e["range"] == [175, 175] for e in tp53)
    assert any(e["range"] == [245, 249] for e in tp53)


def test_override_addition_overwrites_baseline(tmp_path):
    """A higher-layer addition with the same `range` as a baseline entry
    must replace the baseline entry (locked flag and rationale win)."""
    custom_overrides = tmp_path / "overrides.yaml"
    custom_overrides.write_text(
        """
version: "1.0-test"
additions:
  KRAS:
    - range: [12, 13]
      domain: "TEST OVERRIDE DOMAIN"
      strength: moderate
      source: "PMID 99999999"
      rationale: "test-only override"
      locked: true
exclusions: {}
"""
    )
    output = tmp_path / "pm1.json"
    payload = bph.build_pm1_hotspots(
        overrides_path=custom_overrides,
        output_path=output,
    )
    kras = payload["genes"]["KRAS"]
    entry = next((e for e in kras if e["range"] == [12, 13]), None)
    assert entry is not None
    assert entry["domain"] == "TEST OVERRIDE DOMAIN"
    assert entry["source"] == "PMID 99999999"
    assert entry["locked"] is True


# ── Negative validation ──────────────────────────────────────────────────────


def test_validate_rejects_missing_pmid():
    with pytest.raises(ValueError, match="does not contain a PMID"):
        bph._validate({"FOO": [{"range": [1, 1], "strength": "moderate", "source": ""}]})


def test_validate_rejects_bad_strength():
    with pytest.raises(ValueError, match="strength must be one of"):
        bph._validate({
            "FOO": [{"range": [1, 1], "strength": "strong", "source": "PMID 1"}]
        })


def test_validate_rejects_reversed_range():
    with pytest.raises(ValueError, match="start > end"):
        bph._validate({
            "FOO": [{"range": [10, 5], "strength": "moderate", "source": "PMID 1"}]
        })


# ── Coverage assertion (clinical-advisor 2026-04-14) ─────────────────────────


def test_assert_coverage_passes_on_default_build(tmp_path):
    """The default build must satisfy every REQUIRED_COVERAGE pair."""
    output = tmp_path / "pm1.json"
    payload = bph.build_pm1_hotspots(output_path=output)
    # If build_pm1_hotspots returned without raising, coverage held. Double-
    # check directly so regressions here surface as a coverage failure, not a
    # downstream build failure.
    bph._assert_coverage(payload["genes"])


def test_assert_coverage_rejects_missing_residue():
    """Dropping a required residue must fail _assert_coverage with a clear error."""
    genes = {
        "TP53": [
            {"range": [245, 249], "strength": "moderate", "source": "PMID 30224644"},
            # Intentionally missing 175, 273, 282
        ],
    }
    with pytest.raises(ValueError, match="missing required coverage"):
        bph._assert_coverage(genes)


def test_assert_coverage_rejects_missing_gene():
    """A gene with zero entries must fail the coverage check."""
    with pytest.raises(ValueError, match="KRAS:12"):
        bph._assert_coverage({
            "TP53": [
                {"range": [175, 175], "strength": "moderate", "source": "PMID 30224644"},
                {"range": [245, 249], "strength": "moderate", "source": "PMID 30224644"},
                {"range": [273, 273], "strength": "moderate", "source": "PMID 30224644"},
                {"range": [282, 282], "strength": "moderate", "source": "PMID 30224644"},
            ],
        })


def test_source_refs_excludes_cancerhotspots_when_stub(tmp_path):
    """When the cancerhotspots layer is empty (stub/missing), source_refs must
    NOT list cancerhotspots_v2 — clinical-advisor requires source_refs to
    reflect only layers that actually contributed."""
    stub = tmp_path / "stub.tsv"
    stub.write_text('{"status":404,"error":"Not Found"}')
    output = tmp_path / "pm1.json"
    payload = bph.build_pm1_hotspots(
        cancerhotspots_tsv=stub,
        output_path=output,
    )
    assert not any("cancerhotspots" in ref for ref in payload["source_refs"])
    # And the YAML-based anchors must be present.
    assert any("pm1_hotspot_overrides.yaml" in ref for ref in payload["source_refs"])
    assert any("ClinGen_TP53_VCEP" in ref for ref in payload["source_refs"])


def test_stub_fallback_logs_warning(tmp_path, caplog):
    """The 404-stub fallback must emit a WARNING-level log with the
    clinical-advisor-specified message fragment."""
    stub = tmp_path / "stub.tsv"
    stub.write_text('{"status":404,"error":"Not Found"}')
    output = tmp_path / "pm1.json"
    with caplog.at_level("WARNING", logger="scripts.db.build_pm1_hotspots"):
        bph.build_pm1_hotspots(
            cancerhotspots_tsv=stub,
            output_path=output,
        )
    messages = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
    assert any("cancerhotspots_v2 TSV unavailable" in m for m in messages), messages
    assert any("clinical-advisor pre-approved" in m for m in messages), messages
