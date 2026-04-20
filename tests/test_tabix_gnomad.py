# tests/test_tabix_gnomad.py
"""Tests for the gnomAD tabix query module (query_tabix_gnomad.py)."""

from unittest.mock import MagicMock

import pytest

pytest.importorskip("pysam")

from scripts.common.models import Variant

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def reset_tabix_state():
    """Reset module-level tabix file cache before each test."""
    import scripts.storage.query_tabix_gnomad as qmod

    qmod.close()
    yield
    qmod.close()


# ---------------------------------------------------------------------------
# 1. test_find_vcf_for_chrom
# ---------------------------------------------------------------------------


def test_find_vcf_for_chrom(tmp_path, monkeypatch):
    """_find_vcf_for_chrom locates files matching the standard gnomAD naming pattern."""
    import scripts.storage.query_tabix_gnomad as qmod

    # Create a fake VCF bgz + tbi pair
    vcf_file = tmp_path / "gnomad.exomes.v4.1.sites.chr17.vcf.bgz"
    tbi_file = tmp_path / "gnomad.exomes.v4.1.sites.chr17.vcf.bgz.tbi"
    vcf_file.touch()
    tbi_file.touch()

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))

    result = qmod._find_vcf_for_chrom("chr17")
    assert result == str(vcf_file)


def test_find_vcf_for_chrom_strips_chr_prefix(tmp_path, monkeypatch):
    """_find_vcf_for_chrom accepts 'chr17' and '17' as equivalent input."""
    import scripts.storage.query_tabix_gnomad as qmod

    vcf_file = tmp_path / "gnomad.exomes.v4.1.sites.chr7.vcf.bgz"
    tbi_file = tmp_path / "gnomad.exomes.v4.1.sites.chr7.vcf.bgz.tbi"
    vcf_file.touch()
    tbi_file.touch()

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))

    result = qmod._find_vcf_for_chrom("chr7")
    assert result is not None
    assert "chr7" in result


def test_find_vcf_for_chrom_not_found(tmp_path, monkeypatch):
    """_find_vcf_for_chrom returns None when no matching file exists."""
    import scripts.storage.query_tabix_gnomad as qmod

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))

    result = qmod._find_vcf_for_chrom("chr99")
    assert result is None


def test_find_vcf_for_chrom_glob_fallback(tmp_path, monkeypatch):
    """_find_vcf_for_chrom uses glob fallback for non-standard file names."""
    import scripts.storage.query_tabix_gnomad as qmod

    # Non-standard name — won't match direct patterns but will match glob
    vcf_file = tmp_path / "gnomad.exomes.v4.1.custom.chr22.vcf.bgz"
    tbi_file = tmp_path / "gnomad.exomes.v4.1.custom.chr22.vcf.bgz.tbi"
    vcf_file.touch()
    tbi_file.touch()

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))

    result = qmod._find_vcf_for_chrom("chr22")
    assert result == str(vcf_file)


# ---------------------------------------------------------------------------
# 2. test_parse_info_field
# ---------------------------------------------------------------------------


def test_parse_info_field_basic():
    """_parse_info_field correctly splits key=value pairs."""
    from scripts.storage.query_tabix_gnomad import _parse_info_field

    info = "AF=0.001;AF_eas=0.01;AF_afr=0.0005;DP=1000"
    result = _parse_info_field(info)

    assert result["AF"] == "0.001"
    assert result["AF_eas"] == "0.01"
    assert result["AF_afr"] == "0.0005"
    assert result["DP"] == "1000"


def test_parse_info_field_flag():
    """_parse_info_field handles flag fields (no = sign) as 'true'."""
    from scripts.storage.query_tabix_gnomad import _parse_info_field

    info = "PASS;AF=0.05;COMMON"
    result = _parse_info_field(info)

    assert result["PASS"] == "true"
    assert result["COMMON"] == "true"
    assert result["AF"] == "0.05"


def test_parse_info_field_empty():
    """_parse_info_field handles empty or dot INFO field."""
    from scripts.storage.query_tabix_gnomad import _parse_info_field

    result = _parse_info_field(".")
    assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# 3. test_extract_af_single_alt
# ---------------------------------------------------------------------------


def test_extract_af_single_alt():
    """_extract_af returns correct AFs for a single-alt variant."""
    from scripts.storage.query_tabix_gnomad import _extract_af

    info_dict = {
        "AF": "0.001234",
        "AF_eas": "0.01",
        "AF_afr": "0.0005",
        "AF_amr": "0.002",
        "AF_nfe": "0.0008",
        "AF_sas": "0.0003",
    }
    result = _extract_af(info_dict, alt_idx=0)

    assert result["af_global"] == pytest.approx(0.001234)
    assert result["af_eas"] == pytest.approx(0.01)
    assert result["af_afr"] == pytest.approx(0.0005)
    assert result["af_amr"] == pytest.approx(0.002)
    assert result["af_nfe"] == pytest.approx(0.0008)
    assert result["af_sas"] == pytest.approx(0.0003)


def test_extract_af_missing_keys():
    """_extract_af returns None for absent INFO keys."""
    from scripts.storage.query_tabix_gnomad import _extract_af

    info_dict = {"AF": "0.005"}
    result = _extract_af(info_dict, alt_idx=0)

    assert result["af_global"] == pytest.approx(0.005)
    assert result["af_eas"] is None
    assert result["af_afr"] is None


# ---------------------------------------------------------------------------
# 4. test_extract_af_multi_alt
# ---------------------------------------------------------------------------


def test_extract_af_multi_alt_first():
    """_extract_af picks the first allele frequency for alt_idx=0."""
    from scripts.storage.query_tabix_gnomad import _extract_af

    info_dict = {
        "AF": "0.001,0.05",
        "AF_eas": "0.01,0.08",
    }
    result = _extract_af(info_dict, alt_idx=0)

    assert result["af_global"] == pytest.approx(0.001)
    assert result["af_eas"] == pytest.approx(0.01)


def test_extract_af_multi_alt_second():
    """_extract_af picks the second allele frequency for alt_idx=1."""
    from scripts.storage.query_tabix_gnomad import _extract_af

    info_dict = {
        "AF": "0.001,0.05",
        "AF_eas": "0.01,0.08",
    }
    result = _extract_af(info_dict, alt_idx=1)

    assert result["af_global"] == pytest.approx(0.05)
    assert result["af_eas"] == pytest.approx(0.08)


def test_extract_af_multi_alt_out_of_bounds():
    """_extract_af falls back to last element when alt_idx exceeds list length."""
    from scripts.storage.query_tabix_gnomad import _extract_af

    info_dict = {"AF": "0.001,0.05"}
    result = _extract_af(info_dict, alt_idx=5)

    # alt_idx=5 >= len(["0.001","0.05"]) → falls back to parts[0]
    assert result["af_global"] == pytest.approx(0.001)


# ---------------------------------------------------------------------------
# 5. test_get_available_chromosomes
# ---------------------------------------------------------------------------


def test_get_available_chromosomes(tmp_path, monkeypatch):
    """get_available_chromosomes returns chromosomes with VCF files present."""
    import scripts.storage.query_tabix_gnomad as qmod

    # Create fake VCF bgz files for chr1, chr2, chrX
    for chrom in ["chr1", "chr2", "chrX"]:
        (tmp_path / f"gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz").touch()

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))

    chroms = qmod.get_available_chromosomes()

    assert "chr1" in chroms
    assert "chr2" in chroms
    assert "chrX" in chroms
    assert len(chroms) == 3


def test_get_available_chromosomes_empty_dir(tmp_path, monkeypatch):
    """get_available_chromosomes returns empty list when no VCF files present."""
    import scripts.storage.query_tabix_gnomad as qmod

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))

    chroms = qmod.get_available_chromosomes()
    assert chroms == []


def test_get_available_chromosomes_no_dir(tmp_path, monkeypatch):
    """get_available_chromosomes returns empty list when directory doesn't exist."""
    import scripts.storage.query_tabix_gnomad as qmod

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path / "nonexistent"))

    chroms = qmod.get_available_chromosomes()
    assert chroms == []


# ---------------------------------------------------------------------------
# 6. test_get_db_version_not_available
# ---------------------------------------------------------------------------


def test_get_db_version_not_available(tmp_path, monkeypatch):
    """get_db_version returns not_available when VCF dir is missing."""
    import scripts.storage.query_tabix_gnomad as qmod

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path / "nonexistent"))

    meta = qmod.get_db_version()
    assert meta["source"] == "not_available"


def test_get_db_version_empty_dir(tmp_path, monkeypatch):
    """get_db_version returns not_available when VCF dir exists but has no files."""
    import scripts.storage.query_tabix_gnomad as qmod

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))

    meta = qmod.get_db_version()
    assert meta["source"] == "not_available"


def test_get_db_version_with_files(tmp_path, monkeypatch):
    """get_db_version infers version and metadata from VCF filenames."""
    import scripts.storage.query_tabix_gnomad as qmod

    vcf_file = tmp_path / "gnomad.exomes.v4.1.sites.chr1.vcf.bgz"
    vcf_file.touch()

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))

    meta = qmod.get_db_version()

    assert meta["source"] == "tabix"
    assert meta["gnomad_version"] == "4.1"
    assert meta["assembly"] == "GRCh38"
    assert meta["chromosomes"] == 1
    assert "build_date" in meta
    assert "vcf_dir" in meta


# ---------------------------------------------------------------------------
# 7. test_query_tabix_gnomad_no_files
# ---------------------------------------------------------------------------


def test_query_tabix_gnomad_no_files(tmp_path, monkeypatch):
    """query_tabix_gnomad returns api_available=False when no VCF files found."""
    import scripts.storage.query_tabix_gnomad as qmod

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path / "nonexistent"))

    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
    result = qmod.query_tabix_gnomad(variant)

    assert result["gnomad_all"] is None
    assert result["gnomad_eas"] is None
    assert result["api_available"] is False


def test_query_tabix_gnomad_no_tbi(tmp_path, monkeypatch):
    """query_tabix_gnomad returns api_available=False when .tbi index is missing."""
    import scripts.storage.query_tabix_gnomad as qmod

    # Create VCF but no .tbi
    vcf_file = tmp_path / "gnomad.exomes.v4.1.sites.chr17.vcf.bgz"
    vcf_file.touch()

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))

    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
    result = qmod.query_tabix_gnomad(variant)

    assert result["gnomad_all"] is None
    assert result["gnomad_eas"] is None
    assert result["api_available"] is False


def test_query_tabix_gnomad_exact_match(tmp_path, monkeypatch):
    """query_tabix_gnomad returns correct AFs when tabix returns a matching record."""
    import scripts.storage.query_tabix_gnomad as qmod

    # Create VCF and .tbi stubs
    vcf_file = tmp_path / "gnomad.exomes.v4.1.sites.chr17.vcf.bgz"
    tbi_file = tmp_path / "gnomad.exomes.v4.1.sites.chr17.vcf.bgz.tbi"
    vcf_file.touch()
    tbi_file.touch()

    # Build a fake VCF record line (tab-separated)
    fake_record = "\t".join(
        [
            "chr17",
            "7577120",
            "rs28934578",
            "G",
            "A",
            ".",
            "PASS",
            "AF=0.000123;AF_eas=0.001;AF_afr=0.00005;AF_amr=0.0002;AF_nfe=0.0001;AF_sas=0.00003",
        ]
    )

    mock_tbx = MagicMock()
    mock_tbx.fetch.return_value = [fake_record]

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))
    monkeypatch.setattr(qmod, "_get_tabix", lambda chrom: mock_tbx)

    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
    result = qmod.query_tabix_gnomad(variant)

    assert result["api_available"] is True
    assert result["gnomad_all"] == pytest.approx(0.000123)
    assert result["gnomad_eas"] == pytest.approx(0.001)
    assert result["gnomad_afr"] == pytest.approx(0.00005)


def test_query_tabix_gnomad_no_match(tmp_path, monkeypatch):
    """query_tabix_gnomad returns None AFs with api_available=True when no records found."""
    import scripts.storage.query_tabix_gnomad as qmod

    mock_tbx = MagicMock()
    mock_tbx.fetch.return_value = []

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))
    monkeypatch.setattr(qmod, "_get_tabix", lambda chrom: mock_tbx)

    variant = Variant(chrom="chr17", pos=9999999, ref="A", alt="T")
    result = qmod.query_tabix_gnomad(variant)

    assert result["gnomad_all"] is None
    assert result["gnomad_eas"] is None
    assert result["api_available"] is True


def test_query_tabix_gnomad_ref_mismatch(tmp_path, monkeypatch):
    """query_tabix_gnomad skips records where ref allele doesn't match."""
    import scripts.storage.query_tabix_gnomad as qmod

    # Record has ref=C but we query ref=G
    fake_record = "\t".join(["chr17", "7577120", ".", "C", "A", ".", "PASS", "AF=0.001;AF_eas=0.01"])

    mock_tbx = MagicMock()
    mock_tbx.fetch.return_value = [fake_record]

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))
    monkeypatch.setattr(qmod, "_get_tabix", lambda chrom: mock_tbx)

    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    result = qmod.query_tabix_gnomad(variant)

    assert result["gnomad_all"] is None
    assert result["api_available"] is True


def test_query_tabix_gnomad_multi_alt_correct_index(tmp_path, monkeypatch):
    """query_tabix_gnomad picks the correct AF for a multi-allelic ALT match."""
    import scripts.storage.query_tabix_gnomad as qmod

    # Two ALTs: T and C; we want the C allele (alt_idx=1)
    fake_record = "\t".join(["chr17", "7577120", ".", "G", "T,C", ".", "PASS", "AF=0.001,0.05;AF_eas=0.01,0.09"])

    mock_tbx = MagicMock()
    mock_tbx.fetch.return_value = [fake_record]

    monkeypatch.setattr(qmod, "_get_vcf_dir", lambda: str(tmp_path))
    monkeypatch.setattr(qmod, "_get_tabix", lambda chrom: mock_tbx)

    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="C")
    result = qmod.query_tabix_gnomad(variant)

    assert result["api_available"] is True
    assert result["gnomad_all"] == pytest.approx(0.05)
    assert result["gnomad_eas"] == pytest.approx(0.09)


# ---------------------------------------------------------------------------
# 8. test_close
# ---------------------------------------------------------------------------


def test_close_clears_cache(monkeypatch):
    """close() clears the _tabix_files cache without raising."""
    import scripts.storage.query_tabix_gnomad as qmod

    mock_tbx = MagicMock()
    qmod._tabix_files["chr1"] = mock_tbx

    qmod.close()

    assert qmod._tabix_files == {}
    mock_tbx.close.assert_called_once()
