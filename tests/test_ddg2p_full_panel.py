"""Tests for the full EBI Gene2Phenotype DDG2P ingest (v2.3-T7).

These tests only run when ``data/ddg2p_neurodev_genes.json`` has been
populated with the full DDG2P snapshot — i.e. the
``scripts/tools/build_ddg2p_table.py`` build has succeeded with network
access. When the v1 30-gene fallback is in place (or no panel exists at
all), every test in this file is skipped so CI without the EBI download
stays green.
"""

from __future__ import annotations

import json
import os

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PANEL_PATH = os.path.join(REPO_ROOT, "data", "ddg2p_neurodev_genes.json")


def _load_panel() -> dict | None:
    if not os.path.exists(PANEL_PATH):
        return None
    try:
        with open(PANEL_PATH, "r", encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


_PANEL = _load_panel()
_GENE_COUNT = len((_PANEL or {}).get("genes") or {})

pytestmark = pytest.mark.skipif(
    _GENE_COUNT < 100,
    reason=(
        f"DDG2P full panel not present (gene_count={_GENE_COUNT}); "
        "run scripts/tools/build_ddg2p_table.py to enable these tests"
    ),
)


def test_panel_gene_count_reasonable():
    """The full DDG2P admitted set should land somewhere in the 500-3000 range.

    EBI's DDG2P sits around ~2000-2300 admitted (definitive/strong/moderate)
    genes as of 2025-2026. The bounds are deliberately wide so monthly drops
    do not flake the test, while still catching obvious truncation bugs.
    """
    assert 500 <= _GENE_COUNT <= 3000, f"unexpected DDG2P gene count: {_GENE_COUNT}"


def test_panel_contains_well_known_neurodev_genes():
    """All v1 hand-curated genes plus a handful of additional well-known
    neurodev genes must be present in the full ingest."""
    expected = [
        # Original v1 30-gene set
        "SCN1A",
        "SCN2A",
        "SCN8A",
        "KCNQ2",
        "KCNQ3",
        "STXBP1",
        "CDKL5",
        "FOXG1",
        "MECP2",
        "GRIN1",
        "GRIN2A",
        "GRIN2B",
        "ARID1B",
        "SYNGAP1",
        "CHD2",
        "CHD8",
        "TBR1",
        "TCF4",
        "GABRB3",
        "DYRK1A",
        "ANKRD11",
        "ADNP",
        "SETD5",
        "POGZ",
        "ASH1L",
        "SHANK3",
        "NRXN1",
        "PTEN",
        "TSC1",
        "TSC2",
        # Additional well-known neurodev genes that should be in DDG2P
        "GNAO1",
        "SCN3A",
        "STX1B",
        "CACNA1A",
        "NF1",
        "FMR1",
        "MEF2C",
        "UBE3A",
        "HNRNPU",
        "GABRA1",
    ]
    genes = _PANEL["genes"]
    missing = [g for g in expected if g not in genes]
    assert not missing, f"missing expected genes: {missing}"


def test_panel_rejects_limited_confidence():
    """The admitted set must not contain ``limited``, ``disputed``, or
    ``refuted`` records — those are filtered out at build time."""
    rejected = {"limited", "disputed", "refuted"}
    bad = [gene for gene, rec in _PANEL["genes"].items() if (rec.get("confidence") or "").lower() in rejected]
    assert not bad, f"rejected-confidence genes leaked into admission set: {bad[:10]}"


def test_panel_admission_confidences_field():
    """The ``admission_confidences`` metadata array must exactly match the
    spec Q2 list: definitive / strong / moderate."""
    assert _PANEL["admission_confidences"] == ["definitive", "strong", "moderate"]


def test_build_note_has_source_url():
    """build_note + source_url metadata must be populated so reports can cite
    the upstream EBI snapshot."""
    assert _PANEL.get("build_note"), "build_note is empty"
    src = _PANEL.get("source_url", "")
    assert "ebi.ac.uk" in src, f"source_url does not look like an EBI URL: {src!r}"
    assert _PANEL.get("source") == "Gene2Phenotype DDG2P"
    assert _PANEL.get("license") == "CC0"


def test_per_gene_record_shape():
    """Each gene record exposes the v1 schema fields the loader reads."""
    sample = _PANEL["genes"]["SCN1A"]
    for key in ("confidence", "disease", "allelic_requirement"):
        assert key in sample, f"SCN1A record missing {key}"
    assert sample["confidence"] in {"definitive", "strong", "moderate"}


def test_panel_loader_consumes_full_set():
    """``scripts.common.ddg2p_panel`` must still operate on the full ingest
    after the v1 30-gene file is replaced — graceful degradation should not
    trigger and admitted lookups should succeed."""
    from scripts.common import ddg2p_panel

    ddg2p_panel._reset_for_tests()
    assert ddg2p_panel.is_admitted_neurodev_gene("SCN1A") is True
    assert ddg2p_panel.is_admitted_neurodev_gene("GNAO1") is True
    assert ddg2p_panel.is_admitted_neurodev_gene("XYZNOTAREALGENE") is False
    assert ddg2p_panel.is_admitted_neurodev_gene(None) is False
    info = ddg2p_panel.get_neurodev_info("SCN1A")
    assert info is not None
    assert info.get("confidence") in {"definitive", "strong", "moderate"}


def test_version_manager_registration():
    """``version_manager.get_version('DDG2P')`` must return the new entry
    populated from the panel JSON."""
    from scripts.db.version_manager import get_version

    meta = get_version("DDG2P", skip_api=True)
    assert meta is not None, "version_manager did not return a DDG2P entry"
    assert meta["source"] == "local_json"
    assert meta["license"] == "CC0"
    assert meta["upstream"] == "Gene2Phenotype DDG2P"
    assert meta["gene_count"] == _GENE_COUNT
    assert "ebi.ac.uk" in meta["source_url"]
    assert meta["admission_confidences"] == ["definitive", "strong", "moderate"]
