from scripts.common.models import Variant
from scripts.korean_pop.query_korea4k import query_korea4k


def test_query_korea4k_found(tmp_path):
    korea4k_file = tmp_path / "korea4k_freq.tsv"
    korea4k_file.write_text("chr17\t7675088\tC\tA\t0.00008\n")
    variant = Variant(chrom="chr17", pos=7675088, ref="C", alt="A")
    result = query_korea4k(variant, str(korea4k_file))
    assert result == 0.00008


def test_query_korea4k_not_found(tmp_path):
    korea4k_file = tmp_path / "korea4k_freq.tsv"
    korea4k_file.write_text("chr17\t7675088\tC\tA\t0.00008\n")
    variant = Variant(chrom="chr1", pos=99999, ref="A", alt="T")
    result = query_korea4k(variant, str(korea4k_file))
    assert result is None


def test_query_korea4k_no_file(tmp_path):
    """Missing Korea4K file should log a warning and return None (not raise)."""
    from scripts.korean_pop.query_korea4k import _KOREA4K_CACHE

    missing_path = str(tmp_path / "nonexistent.tsv")
    # Remove from cache if previously cached (test isolation)
    _KOREA4K_CACHE.pop(missing_path, None)
    variant = Variant(chrom="chr17", pos=7675088, ref="C", alt="A")
    result = query_korea4k(variant, missing_path)
    assert result is None


def test_korea4k_cache_singleton(tmp_path):
    """Verify that loading the same file twice reuses the cached dict."""
    from scripts.korean_pop.query_korea4k import _KOREA4K_CACHE, _load_korea4k

    korea4k_file = tmp_path / "korea4k_cache_test.tsv"
    korea4k_file.write_text("chr10\t94781859\tG\tA\t0.31\n")
    path = str(korea4k_file)
    _KOREA4K_CACHE.pop(path, None)

    data1 = _load_korea4k(path)
    data2 = _load_korea4k(path)
    assert data1 is data2  # same object, not a reload
