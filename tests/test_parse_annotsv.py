def test_parse_cancer_annotsv():
    """Cancer AnnotSV TSV 파싱 — 8 full SVs."""
    from scripts.intake.parse_annotsv import parse_annotsv

    svs = parse_annotsv("data/sample_sv/cancer_somatic_annotsv.tsv")
    assert len(svs) == 8  # 8 unique SVs (full rows)
    # ERBB2 amplification
    erbb2 = next(sv for sv in svs if "ERBB2" in sv.gene_name)
    assert erbb2.sv_type == "DUP"
    assert erbb2.acmg_class == 5
    assert erbb2.cytoband == "17q12"
    assert len(erbb2.gene_details) >= 1
    assert erbb2.gene_details[0]["gene"] == "ERBB2"


def test_parse_rare_disease_annotsv():
    """Rare disease AnnotSV TSV 파싱 — 7 full SVs."""
    from scripts.intake.parse_annotsv import parse_annotsv

    svs = parse_annotsv("data/sample_sv/rare_disease_annotsv.tsv")
    assert len(svs) == 7
    # 22q11 deletion spans 3 genes
    del22q = next(sv for sv in svs if "22q11" in sv.cytoband)
    assert del22q.gene_count == 3
    assert "TBX1" in del22q.gene_name
    assert len(del22q.gene_details) == 3  # TBX1, COMT, HIRA


def test_parse_split_rows_populate_gene_details():
    """Split rows가 gene_details에 채워짐."""
    from scripts.intake.parse_annotsv import parse_annotsv

    svs = parse_annotsv("data/sample_sv/cancer_somatic_annotsv.tsv")
    erbb2 = next(sv for sv in svs if "ERBB2" in sv.gene_name)
    gd = erbb2.gene_details[0]
    assert gd["gene"] == "ERBB2"
    assert gd["transcript"] == "NM_004448.4"
    assert gd["cds_percent"] == 100.0
    assert gd["hi"] == 2
    assert gd["ts"] == 3


def test_parse_sv_types():
    """DEL, DUP, INV 타입이 올바르게 파싱됨."""
    from scripts.intake.parse_annotsv import parse_annotsv

    svs = parse_annotsv("data/sample_sv/cancer_somatic_annotsv.tsv")
    types = {sv.sv_type for sv in svs}
    assert "DEL" in types
    assert "DUP" in types
    assert "INV" in types


def test_parse_empty_file(tmp_path):
    """빈 파일은 빈 리스트 반환."""
    from scripts.intake.parse_annotsv import parse_annotsv

    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    assert parse_annotsv(str(empty)) == []


def test_parse_benign_has_af():
    """Benign SV는 AF 정보를 가짐."""
    from scripts.intake.parse_annotsv import parse_annotsv

    svs = parse_annotsv("data/sample_sv/cancer_somatic_annotsv.tsv")
    benign = next(sv for sv in svs if sv.acmg_class == 1)
    assert benign.b_loss_af_max is not None or benign.b_gain_af_max is not None
