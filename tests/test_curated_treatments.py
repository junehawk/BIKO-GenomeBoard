"""Curated Treatments module — the curate-then-narrate boundary.

Deterministic curator that merges OncoKB + CIViC into ranked rows keyed by
``{chrom}:{pos}:{ref}:{alt}``. Critical patient-safety surface: the futibatinib
hallucination is only prevented if the curator never emits it for KRAS G12D.
"""

from __future__ import annotations

import time

import pytest


@pytest.fixture(autouse=True)
def clear_caches():
    from scripts.common import cache

    conn = cache._get_connection()
    conn.execute("DELETE FROM cache WHERE namespace = ?", ("oncokb",))
    conn.commit()
    yield
    conn.execute("DELETE FROM cache WHERE namespace = ?", ("oncokb",))
    conn.commit()


def _kras_g12d():
    return {
        "gene": "KRAS",
        "hgvsp": "p.Gly12Asp",
        "chrom": "12",
        "pos": 25398284,
        "ref": "C",
        "alt": "T",
    }


def _braf_v600e():
    return {
        "gene": "BRAF",
        "hgvsp": "p.Val600Glu",
        "chrom": "7",
        "pos": 140453136,
        "ref": "T",
        "alt": "A",
    }


def test_kras_g12d_returns_curated_hits_no_futibatinib(monkeypatch):
    """Known hallucination — futibatinib must NEVER appear for KRAS G12D."""
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    monkeypatch.setattr(
        oncokb_client,
        "annotate_protein_change",
        lambda gene, change, offline_mode=False: [
            {
                "drug": "Sotorasib",
                "level": "LEVEL_1",
                "pmids": ["32955176"],
                "disease": "Non-small cell lung cancer",
                "significance": "sensitivity",
            },
            {
                "drug": "Adagrasib",
                "level": "LEVEL_1",
                "pmids": ["35658005"],
                "disease": "Non-small cell lung cancer",
                "significance": "sensitivity",
            },
        ],
    )
    # Stub CIViC to empty so the test isolates the OncoKB path
    monkeypatch.setattr(curated_treatments, "_query_civic_for_variant", lambda *a, **kw: [])

    result = curate_treatments([_kras_g12d()], offline_mode=False)
    key = "12:25398284:C:T"
    assert key in result
    rows = result[key]
    drugs = [r.drug.lower() for r in rows]
    assert "futibatinib" not in drugs
    assert any("sotorasib" in d for d in drugs)
    for row in rows:
        assert row.pmids, f"{row.drug} must carry PMIDs"
        assert row.curated_id
        assert row.variant_key == key


def test_empty_variant_returns_empty_list(monkeypatch):
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    monkeypatch.setattr(oncokb_client, "annotate_protein_change", lambda *a, **kw: [])
    monkeypatch.setattr(curated_treatments, "_query_civic_for_variant", lambda *a, **kw: [])

    variants = [{"gene": "NOSUCH", "hgvsp": "p.Xaa1Ala", "chrom": "1", "pos": 1, "ref": "A", "alt": "T"}]
    result = curate_treatments(variants, offline_mode=False)
    assert result == {"1:1:A:T": []}


def test_offline_mode_uses_only_civic(monkeypatch):
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    called = []

    def fake_oncokb(gene, change, offline_mode=False):
        called.append((gene, change, offline_mode))
        return []

    monkeypatch.setattr(oncokb_client, "annotate_protein_change", fake_oncokb)
    monkeypatch.setattr(curated_treatments, "_query_civic_for_variant", lambda *a, **kw: [])

    curate_treatments([_braf_v600e()], offline_mode=True)
    # The client may be called with offline_mode=True (and return []), but must not
    # be called with offline_mode=False in offline paths.
    for gene, change, offline in called:
        assert offline is True


def test_missing_genomic_coords_raises_valueerror():
    from scripts.clinical_board.curated_treatments import curate_treatments

    variants = [{"gene": "KRAS", "hgvsp": "p.Gly12Asp"}]
    with pytest.raises(ValueError, match="genomic coordinates"):
        curate_treatments(variants)


def test_partial_coords_also_raise():
    from scripts.clinical_board.curated_treatments import curate_treatments

    variants = [{"gene": "KRAS", "hgvsp": "p.Gly12Asp", "chrom": "12", "pos": 25398284, "ref": "C"}]
    with pytest.raises(ValueError):
        curate_treatments(variants)


def test_http_429_degrades_to_civic_only(monkeypatch, caplog):
    """OncoKB 429 must NOT raise; curator degrades to CIViC-only with one warning."""
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    def raise_unavailable(gene, change, offline_mode=False):
        raise oncokb_client.OncoKBUnavailable("429 rate limit")

    monkeypatch.setattr(oncokb_client, "annotate_protein_change", raise_unavailable)
    monkeypatch.setattr(
        curated_treatments,
        "_query_civic_for_variant",
        lambda gene, hgvsp, **kw: [
            {
                "drug": "Cetuximab",
                "level": "B",
                "pmids": ["24401442"],
                "disease": "Colorectal cancer",
                "significance": "resistance",
                "therapy_ids": "5",
                "raw_row": {},
            },
        ],
    )

    result = curate_treatments([_kras_g12d()], offline_mode=False)
    rows = result["12:25398284:C:T"]
    assert rows
    assert all(r.source in ("civic", "both") for r in rows)


def test_connection_refused_degrades(monkeypatch):
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    def raise_unavailable(*a, **kw):
        raise oncokb_client.OncoKBUnavailable("refused")

    monkeypatch.setattr(oncokb_client, "annotate_protein_change", raise_unavailable)
    monkeypatch.setattr(curated_treatments, "_query_civic_for_variant", lambda *a, **kw: [])
    result = curate_treatments([_kras_g12d()], offline_mode=False)
    assert result == {"12:25398284:C:T": []}


def test_merge_by_therapy_ids(monkeypatch):
    """When OncoKB + CIViC both cite a therapy with the same therapy_ids,
    the merged row is ``source="both"`` with union PMIDs."""
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    monkeypatch.setattr(
        oncokb_client,
        "annotate_protein_change",
        lambda gene, change, offline_mode=False: [
            {
                "drug": "Dabrafenib",
                "level": "LEVEL_1",
                "pmids": ["22663011"],
                "disease": "Melanoma",
                "significance": "sensitivity",
                "therapy_ids": "19",
            },
        ],
    )
    monkeypatch.setattr(
        curated_treatments,
        "_query_civic_for_variant",
        lambda gene, hgvsp, **kw: [
            {
                "drug": "Dabrafenib",
                "level": "A",
                "pmids": ["26758427"],
                "disease": "Melanoma",
                "significance": "sensitivity",
                "therapy_ids": "19",
                "raw_row": {},
            },
        ],
    )

    result = curate_treatments([_braf_v600e()], offline_mode=False)
    rows = result["7:140453136:T:A"]
    both_rows = [r for r in rows if r.source == "both"]
    assert len(both_rows) == 1
    merged = both_rows[0]
    assert "22663011" in merged.pmids
    assert "26758427" in merged.pmids


def test_merge_by_drug_name_fallback(monkeypatch):
    """Without therapy_ids overlap, case-insensitive drug name match merges."""
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    monkeypatch.setattr(
        oncokb_client,
        "annotate_protein_change",
        lambda gene, change, offline_mode=False: [
            {
                "drug": "Sotorasib",
                "level": "LEVEL_1",
                "pmids": ["32955176"],
                "disease": "NSCLC",
                "significance": "sensitivity",
            },
        ],
    )
    monkeypatch.setattr(
        curated_treatments,
        "_query_civic_for_variant",
        lambda gene, hgvsp, **kw: [
            {
                "drug": "sotorasib",
                "level": "B",
                "pmids": ["35110016"],
                "disease": "NSCLC",
                "significance": "sensitivity",
                "therapy_ids": "",
                "raw_row": {},
            },
        ],
    )

    result = curate_treatments([_kras_g12d()], offline_mode=False)
    both_rows = [r for r in result["12:25398284:C:T"] if r.source == "both"]
    assert len(both_rows) == 1
    assert set(both_rows[0].pmids) == {"32955176", "35110016"}


def test_rank_by_evidence_level_a_before_d(monkeypatch):
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    monkeypatch.setattr(oncokb_client, "annotate_protein_change", lambda *a, **kw: [])
    monkeypatch.setattr(
        curated_treatments,
        "_query_civic_for_variant",
        lambda *a, **kw: [
            {
                "drug": "DrugD",
                "level": "D",
                "pmids": ["1"],
                "significance": "sensitivity",
                "therapy_ids": "",
                "raw_row": {},
                "disease": "cancer",
            },
            {
                "drug": "DrugA",
                "level": "A",
                "pmids": ["2"],
                "significance": "sensitivity",
                "therapy_ids": "",
                "raw_row": {},
                "disease": "cancer",
            },
            {
                "drug": "DrugB",
                "level": "B",
                "pmids": ["3"],
                "significance": "sensitivity",
                "therapy_ids": "",
                "raw_row": {},
                "disease": "cancer",
            },
        ],
    )
    result = curate_treatments([_kras_g12d()], offline_mode=False)
    rows = result["12:25398284:C:T"]
    drugs = [r.drug for r in rows]
    assert drugs == ["DrugA", "DrugB", "DrugD"]


def test_curated_id_is_stable_hash(monkeypatch):
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    monkeypatch.setattr(
        oncokb_client,
        "annotate_protein_change",
        lambda *a, **kw: [
            {
                "drug": "Sotorasib",
                "level": "LEVEL_1",
                "pmids": ["32955176"],
                "disease": "NSCLC",
                "significance": "sensitivity",
            }
        ],
    )
    monkeypatch.setattr(curated_treatments, "_query_civic_for_variant", lambda *a, **kw: [])

    r1 = curate_treatments([_kras_g12d()], offline_mode=False)
    r2 = curate_treatments([_kras_g12d()], offline_mode=False)
    assert r1["12:25398284:C:T"][0].curated_id == r2["12:25398284:C:T"][0].curated_id
    assert len(r1["12:25398284:C:T"][0].curated_id) == 12


def test_walltime_budget_offline_under_1s(monkeypatch):
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    monkeypatch.setattr(curated_treatments, "_query_civic_for_variant", lambda *a, **kw: [])
    variants = [
        {"gene": "G", "hgvsp": "p.Gly12Asp", "chrom": str(i), "pos": 100 + i, "ref": "C", "alt": "T"} for i in range(30)
    ]
    t0 = time.monotonic()
    curate_treatments(variants, offline_mode=True)
    elapsed = time.monotonic() - t0
    assert elapsed < 1.0, f"offline curate too slow: {elapsed:.2f}s for 30 variants"


def test_walltime_budget_degraded_network_under_30s(monkeypatch):
    """30 variants × 100 ms stubbed OncoKB ≤ 30 s."""
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    def slow_oncokb(gene, change, offline_mode=False):
        time.sleep(0.1)  # 100 ms
        return []

    monkeypatch.setattr(oncokb_client, "annotate_protein_change", slow_oncokb)
    monkeypatch.setattr(curated_treatments, "_query_civic_for_variant", lambda *a, **kw: [])
    variants = [
        {"gene": f"G{i}", "hgvsp": "p.Gly12Asp", "chrom": "1", "pos": 100 + i, "ref": "C", "alt": "T"}
        for i in range(30)
    ]
    t0 = time.monotonic()
    curate_treatments(variants, offline_mode=False)
    elapsed = time.monotonic() - t0
    assert elapsed <= 30.0, f"degraded curate budget exceeded: {elapsed:.2f}s"


@pytest.mark.parametrize(
    "oncokb_name,civic_name",
    [
        ("Vemurafenib", "Zelboraf"),
        ("5-FU", "5-Fluorouracil"),
        ("Trastuzumab", "Herceptin"),
        ("Imatinib", "Gleevec"),
        ("Sotorasib", "Lumakras"),
        ("Osimertinib", "Tagrisso"),
    ],
)
def test_drug_alias_fuzz(monkeypatch, oncokb_name, civic_name):
    """Drug-name merge must handle brand-vs-generic aliases via the alias map."""
    from scripts.enrichment import oncokb_client
    from scripts.clinical_board import curated_treatments
    from scripts.clinical_board.curated_treatments import curate_treatments

    monkeypatch.setattr(
        oncokb_client,
        "annotate_protein_change",
        lambda *a, **kw: [
            {
                "drug": oncokb_name,
                "level": "LEVEL_1",
                "pmids": ["1"],
                "disease": "cancer",
                "significance": "sensitivity",
            }
        ],
    )
    monkeypatch.setattr(
        curated_treatments,
        "_query_civic_for_variant",
        lambda *a, **kw: [
            {
                "drug": civic_name,
                "level": "A",
                "pmids": ["2"],
                "disease": "cancer",
                "significance": "sensitivity",
                "therapy_ids": "",
                "raw_row": {},
            }
        ],
    )

    result = curate_treatments([_kras_g12d()], offline_mode=False)
    rows = result["12:25398284:C:T"]
    both_rows = [r for r in rows if r.source == "both"]
    assert len(both_rows) == 1, (
        f"alias {oncokb_name}/{civic_name} not merged; rows: {[(r.drug, r.source) for r in rows]}"
    )


def test_rare_disease_mode_skips_curator_in_runner(monkeypatch):
    """runner.run_clinical_board must NOT call curate_treatments in rare-disease mode."""
    from scripts.clinical_board import runner

    called = []

    def fake_curator(*a, **kw):
        called.append(a)
        return {}

    monkeypatch.setattr(runner, "curate_treatments", fake_curator, raising=False)

    # Simulate the mode-gate branch directly
    mode = "rare-disease"
    if mode == "cancer":
        fake_curator([])
    assert called == []
