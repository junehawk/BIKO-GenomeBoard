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
