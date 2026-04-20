from scripts.common.models import Variant
from scripts.population.query_nard2 import query_nard2


def test_query_nard2_found(tmp_path):
    nard2_file = tmp_path / "nard2_freq.tsv"
    nard2_file.write_text("chr17\t7675088\tC\tA\t0.00006\n")
    variant = Variant(chrom="chr17", pos=7675088, ref="C", alt="A")
    result = query_nard2(variant, str(nard2_file))
    assert result == 0.00006


def test_query_nard2_not_found(tmp_path):
    nard2_file = tmp_path / "nard2_freq.tsv"
    nard2_file.write_text("chr17\t7675088\tC\tA\t0.00006\n")
    variant = Variant(chrom="chr1", pos=99999, ref="A", alt="T")
    result = query_nard2(variant, str(nard2_file))
    assert result is None


def test_query_nard2_no_file(tmp_path):
    """Missing NARD2 file should log a warning and return None (not raise)."""
    from scripts.population.query_nard2 import _NARD2_CACHE

    missing_path = str(tmp_path / "nonexistent.tsv")
    # Remove from cache if previously cached (test isolation)
    _NARD2_CACHE.pop(missing_path, None)
    variant = Variant(chrom="chr17", pos=7675088, ref="C", alt="A")
    result = query_nard2(variant, missing_path)
    assert result is None


def test_nard2_cache_singleton(tmp_path):
    """Verify that loading the same file twice reuses the cached dict."""
    from scripts.population.query_nard2 import _NARD2_CACHE, _load_nard2

    nard2_file = tmp_path / "nard2_cache_test.tsv"
    nard2_file.write_text("chr10\t94781859\tG\tA\t0.29\n")
    path = str(nard2_file)
    _NARD2_CACHE.pop(path, None)

    data1 = _load_nard2(path)
    data2 = _load_nard2(path)
    assert data1 is data2  # same object, not a reload
