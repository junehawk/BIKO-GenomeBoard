from scripts.korean_pop.query_krgdb import query_krgdb
from scripts.common.models import Variant

def test_query_krgdb_found(tmp_path):
    krgdb_file = tmp_path / "krgdb_freq.tsv"
    krgdb_file.write_text("chr17\t7577120\tG\tA\t0.0001\n")
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    result = query_krgdb(variant, str(krgdb_file))
    assert result == 0.0001

def test_query_krgdb_not_found(tmp_path):
    krgdb_file = tmp_path / "krgdb_freq.tsv"
    krgdb_file.write_text("chr17\t7577120\tG\tA\t0.0001\n")
    variant = Variant(chrom="chr1", pos=99999, ref="A", alt="T")
    result = query_krgdb(variant, str(krgdb_file))
    assert result is None

def test_query_krgdb_missing_file(tmp_path):
    """Missing KRGDB file should log a warning and return None (not raise)."""
    from scripts.korean_pop.query_krgdb import _KRGDB_CACHE
    missing_path = str(tmp_path / "nonexistent.tsv")
    # Remove from cache if previously cached (test isolation)
    _KRGDB_CACHE.pop(missing_path, None)
    variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
    result = query_krgdb(variant, missing_path)
    assert result is None
