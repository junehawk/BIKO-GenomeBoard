"""Tests for data/korean_pgx_table.json completeness and validity.

Validates that the builtin PGx table covers all PharmCAT 3.2.0
pharmacogenes and that every entry has valid fields.
"""

import json
from pathlib import Path

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

_TABLE_PATH = Path(__file__).parent.parent / "data" / "korean_pgx_table.json"


def _load_table() -> dict:
    with open(_TABLE_PATH) as f:
        return json.load(f)


# PharmCAT 3.2.0 reports 23 pharmacogenes.  The builtin PGx table must
# contain an entry for every one so that the PharmCAT-to-builtin merge
# path (_convert_pharmcat_to_pgx_hits) can provide CPIC level and
# Korean prevalence for all genes.
# PharmCAT 3.2.0 standard 23-gene panel.
_PHARMCAT_GENES = frozenset(
    [
        "ABCG2",
        "CACNA1S",
        "CFTR",
        "CYP2B6",
        "CYP2C9",
        "CYP2C19",
        "CYP2D6",
        "CYP3A4",
        "CYP3A5",
        "CYP4F2",
        "DPYD",
        "G6PD",
        "HLA-A",
        "HLA-B",
        "IFNL3",
        "MT-RNR1",
        "NAT2",
        "NUDT15",
        "RYR1",
        "SLCO1B1",
        "TPMT",
        "UGT1A1",
        "VKORC1",
    ]
)


# ---------------------------------------------------------------------------
# 1. PharmCAT coverage — every PharmCAT gene must be in the builtin table
# ---------------------------------------------------------------------------


def test_all_pharmcat_genes_in_builtin():
    """All 23 PharmCAT pharmacogenes should have entries in the builtin table."""
    data = _load_table()
    builtin_genes = {entry["gene"] for entry in data["genes"]}

    missing = _PHARMCAT_GENES - builtin_genes
    assert not missing, f"PharmCAT genes missing from builtin table: {sorted(missing)}"


# ---------------------------------------------------------------------------
# 2. CPIC level — every entry must have a non-empty cpic_level
# ---------------------------------------------------------------------------


def test_cpic_level_not_empty():
    """Every PGx table entry must have a non-empty cpic_level."""
    data = _load_table()
    for entry in data["genes"]:
        gene = entry["gene"]
        cpic = entry.get("cpic_level", "")
        assert cpic, f"{gene} has empty cpic_level"
        assert cpic in ("A", "B", "C", "D"), f"{gene} has invalid cpic_level '{cpic}'"


# ---------------------------------------------------------------------------
# 3. Frequency validation — korean_freq and western_freq in [0, 1]
# ---------------------------------------------------------------------------


def test_korean_freq_reasonable():
    """korean_freq must be a float in [0, 1]."""
    data = _load_table()
    for entry in data["genes"]:
        gene = entry["gene"]
        freq = entry.get("korean_freq")
        assert freq is not None, f"{gene} missing korean_freq"
        assert isinstance(freq, (int, float)), f"{gene} korean_freq is not numeric: {freq}"
        assert 0 <= freq <= 1, f"{gene} korean_freq out of range: {freq}"


def test_western_freq_reasonable():
    """western_freq must be a float in [0, 1]."""
    data = _load_table()
    for entry in data["genes"]:
        gene = entry["gene"]
        freq = entry.get("western_freq")
        assert freq is not None, f"{gene} missing western_freq"
        assert isinstance(freq, (int, float)), f"{gene} western_freq is not numeric: {freq}"
        assert 0 <= freq <= 1, f"{gene} western_freq out of range: {freq}"


# ---------------------------------------------------------------------------
# 4. Required fields — every entry must have the schema keys
# ---------------------------------------------------------------------------


def test_required_fields_present():
    """Every entry must have gene, korean_freq, western_freq, clinical_impact, cpic_level, source."""
    data = _load_table()
    required = {"gene", "korean_freq", "western_freq", "clinical_impact", "cpic_level", "source"}
    for entry in data["genes"]:
        gene = entry.get("gene", "<unknown>")
        missing = required - set(entry.keys())
        assert not missing, f"{gene} missing required fields: {missing}"


# ---------------------------------------------------------------------------
# 5. No duplicate genes
# ---------------------------------------------------------------------------


def test_no_duplicate_genes():
    """Each gene should appear exactly once in the table."""
    data = _load_table()
    genes = [entry["gene"] for entry in data["genes"]]
    duplicates = [g for g in genes if genes.count(g) > 1]
    assert not duplicates, f"Duplicate genes in PGx table: {set(duplicates)}"
