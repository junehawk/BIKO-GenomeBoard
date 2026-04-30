"""Smoke coverage for the static OMIM fallback dict.

The genemap2.txt SQLite DB is licensed and not present on CI. The static
``_STATIC_OMIM`` dict in ``scripts.enrichment.query_omim`` is the
fallback that gives the rare-disease report a populated ``inheritance``
column on common actionable genes. v2.6 expanded the dict from 11 to
~80 genes; this module pins the new shape.
"""

from __future__ import annotations

from scripts.enrichment.query_omim import _STATIC_OMIM, query_omim


def test_static_dict_size_baseline():
    """v2.6 expanded the static fallback from 11 → ~80 entries."""
    assert len(_STATIC_OMIM) >= 80, (
        f"static fallback has {len(_STATIC_OMIM)} entries — expansion regressed below baseline"
    )


def test_acmg_sf_v3_2_core_genes_present():
    """Spot-check of ACMG SF v3.2 secondary-findings core gene coverage."""
    sf_core = [
        # HBOC
        "BRCA1",
        "BRCA2",
        "PALB2",
        # Lynch
        "MLH1",
        "MSH2",
        "MSH6",
        "PMS2",
        # Polyposis
        "APC",
        "MUTYH",
        # Cardiomyopathy
        "MYH7",
        "MYBPC3",
        # Long QT
        "KCNQ1",
        "KCNH2",
        "SCN5A",
        # Aortopathy
        "FBN1",
        "TGFBR1",
        "TGFBR2",
        # Familial hypercholesterolemia
        "LDLR",
        "APOB",
        "PCSK9",
    ]
    for gene in sf_core:
        assert gene in _STATIC_OMIM, f"ACMG SF v3.2 core gene missing from static fallback: {gene}"


def test_neurodevelopmental_panel_subset_present():
    """High-yield rare-disease / neurodevelopmental genes."""
    nd = ["CHD7", "CHD8", "ARID1A", "ARID1B", "KMT2D", "MECP2", "FMR1", "SCN1A", "SCN2A", "STXBP1"]
    for gene in nd:
        assert gene in _STATIC_OMIM, f"neurodev panel gene missing: {gene}"


def test_each_entry_has_required_fields():
    """Every entry must carry mim + phenotypes + inheritance."""
    for gene, payload in _STATIC_OMIM.items():
        assert "mim" in payload and payload["mim"], f"{gene}: empty mim"
        assert "phenotypes" in payload and payload["phenotypes"], f"{gene}: empty phenotypes"
        assert "inheritance" in payload, f"{gene}: missing inheritance"
        assert isinstance(payload["phenotypes"], list), f"{gene}: phenotypes must be list"


def test_query_omim_returns_static_for_new_gene():
    """query_omim() should fall back to the static dict and tag source='static'."""
    result = query_omim("BRCA1")
    assert result is not None
    assert result["mim"] == "113705"
    assert result["source"] == "static"
    assert result["inheritance"] == "AD"


def test_query_omim_returns_none_for_unknown():
    assert query_omim("NOT_A_REAL_GENE_NAME_42") is None
