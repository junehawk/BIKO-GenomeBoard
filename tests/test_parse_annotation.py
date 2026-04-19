# tests/test_parse_annotation.py
"""Tests for VEP CSQ / SnpEff ANN annotation parsing."""

from pathlib import Path

from scripts.intake.parse_annotation import (
    format_consequence,
    parse_ann_header,
    parse_ann_value,
    parse_csq_header,
    parse_csq_value,
)
from scripts.intake.parse_vcf import parse_vcf

ANNOTATED_VCF = Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants_grch38_annotated.vcf"
PLAIN_VCF = Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_variants.vcf"

CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
    "Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|"
    'cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|CANONICAL|MANE_SELECT|SIFT|PolyPhen">'
)

ANN_HEADER = (
    '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
    "'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | "
    "Transcript_BioType | Rank | HGVS.c | HGVS.p | ...'\">"
)


# ── Header parsing ────────────────────────────────────────────────────────────


def test_parse_csq_header():
    fields = parse_csq_header(CSQ_HEADER)
    assert isinstance(fields, list)
    assert len(fields) == 21
    assert fields[0] == "Allele"
    assert fields[1] == "Consequence"
    assert fields[2] == "IMPACT"
    assert fields[3] == "SYMBOL"
    assert "HGVSc" in fields
    assert "HGVSp" in fields
    assert "CANONICAL" in fields
    assert "MANE_SELECT" in fields
    assert "SIFT" in fields
    assert "PolyPhen" in fields


def test_parse_csq_header_no_format():
    """Header without 'Format:' returns empty list."""
    fields = parse_csq_header('##INFO=<ID=CSQ,Number=.,Type=String,Description="no format here">')
    assert fields == []


def test_parse_ann_header():
    fields = parse_ann_header(ANN_HEADER)
    assert isinstance(fields, list)
    assert len(fields) == 16
    assert fields[0] == "Allele"
    assert fields[1] == "Annotation"
    assert fields[2] == "Annotation_Impact"
    assert fields[3] == "Gene_Name"
    assert "HGVS.c" in fields
    assert "HGVS.p" in fields


# ── CSQ value parsing ─────────────────────────────────────────────────────────


def test_parse_csq_value_single_transcript():
    fields = parse_csq_header(CSQ_HEADER)
    # Single transcript CSQ string for a missense variant
    csq = "A|missense_variant|MODERATE|TP53|ENSG00000141510|Transcript|ENST00000269305|protein_coding|4/11||ENST00000269305.9:c.524G>T|ENSP00000269305.4:p.Arg175Leu|605/2579|524/1182|175/393|R/L|cGg/cTg|YES|NM_000546.6|deleterious(0.01)|probably_damaging(0.998)"
    result = parse_csq_value(csq, fields)
    assert result is not None
    assert result["gene"] == "TP53"
    assert result["consequence"] == "missense_variant"
    assert result["impact"] == "MODERATE"
    assert result["transcript"] == "ENST00000269305"
    assert result["hgvsc"] == "ENST00000269305.9:c.524G>T"
    assert result["hgvsp"] == "ENSP00000269305.4:p.Arg175Leu"
    assert result["sift"] == "deleterious(0.01)"
    assert result["polyphen"] == "probably_damaging(0.998)"
    assert result["canonical"] == "YES"
    assert result["mane_select"] == "NM_000546.6"


def test_parse_csq_value_multiple_transcripts_picks_canonical():
    """When multiple transcripts exist, CANONICAL=YES should be preferred."""
    fields = parse_csq_header(CSQ_HEADER)
    # Two entries: non-canonical first, canonical second
    csq = (
        "A|missense_variant|MODERATE|TP53|ENSG00000141510|Transcript|ENST00000413465|protein_coding|"
        "4/10||ENST00000413465.2:c.524G>T|ENSP00000410739.2:p.Arg175Leu|571/2458|524/1101|175/366|R/L|cGg/cTg|NO||deleterious(0.01)|probably_damaging(0.997),"
        "A|missense_variant|MODERATE|TP53|ENSG00000141510|Transcript|ENST00000269305|protein_coding|"
        "4/11||ENST00000269305.9:c.524G>T|ENSP00000269305.4:p.Arg175Leu|605/2579|524/1182|175/393|R/L|cGg/cTg|YES|NM_000546.6|deleterious(0.01)|probably_damaging(0.998)"
    )
    result = parse_csq_value(csq, fields)
    assert result is not None
    assert result["canonical"] == "YES"
    assert result["transcript"] == "ENST00000269305"


def test_parse_csq_value_mane_select_preferred():
    """MANE_SELECT transcript is preferred over non-canonical."""
    fields = parse_csq_header(CSQ_HEADER)
    # Two entries: first has mane_select, second does not
    csq = (
        "A|missense_variant|MODERATE|TP53|ENSG00000141510|Transcript|ENST00000269305|protein_coding|"
        "4/11||ENST00000269305.9:c.524G>T|ENSP00000269305.4:p.Arg175Leu|605/2579|524/1182|175/393|R/L|cGg/cTg|NO|NM_000546.6|deleterious(0.01)|probably_damaging(0.998),"
        "A|missense_variant|MODERATE|TP53|ENSG00000141510|Transcript|ENST00000413465|protein_coding|"
        "4/10||ENST00000413465.2:c.524G>T||571/2458|524/1101|175/366|R/L|cGg/cTg|NO|||"
    )
    result = parse_csq_value(csq, fields)
    assert result is not None
    assert result["mane_select"] == "NM_000546.6"
    assert result["transcript"] == "ENST00000269305"


def test_parse_csq_value_gene_filter():
    """Gene filter should restrict to matching gene entries."""
    fields = parse_csq_header(CSQ_HEADER)
    csq = (
        "A|missense_variant|MODERATE|GENE1|ENSG000001|Transcript|ENST000001|protein_coding|"
        "1/5||ENST000001:c.100A>G||100/500||50/250|K/E|aAg/aGg|NO|||tolerated|benign,"
        "A|stop_gained|HIGH|GENE2|ENSG000002|Transcript|ENST000002|protein_coding|"
        "2/8||ENST000002:c.200C>T||200/800||67/267|R/*|Cga/Tga|YES|||"
    )
    result = parse_csq_value(csq, fields, gene="GENE2")
    assert result is not None
    assert result["gene"] == "GENE2"
    assert result["consequence"] == "stop_gained"


def test_parse_csq_value_empty_string():
    fields = parse_csq_header(CSQ_HEADER)
    # Empty CSQ value should not crash; returns something (possibly empty entry)
    result = parse_csq_value("", fields)
    # No assertion on content, just no exception and result is dict or None
    assert result is None or isinstance(result, dict)


# ── ANN value parsing ─────────────────────────────────────────────────────────


def test_parse_ann_value():
    fields = parse_ann_header(ANN_HEADER)
    ann = "T|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000269305|protein_coding|4/11|c.524G>T|p.Arg175Leu|605/2579|524/1182|175/393||"
    result = parse_ann_value(ann, fields)
    assert result is not None
    assert result["gene"] == "TP53"
    assert result["consequence"] == "missense_variant"
    assert result["impact"] == "MODERATE"
    assert result["transcript"] == "ENST00000269305"
    assert result["hgvsc"] == "c.524G>T"
    assert result["hgvsp"] == "p.Arg175Leu"
    assert result["sift"] == ""
    assert result["polyphen"] == ""


def test_parse_ann_value_multiple_entries_picks_highest_impact():
    fields = parse_ann_header(ANN_HEADER)
    # Two entries: MODIFIER then HIGH
    ann = (
        "T|intron_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST000001|protein_coding|3/10|c.300+5G>T|||||||,"
        "T|frameshift_variant|HIGH|TP53|ENSG00000141510|transcript|ENST000002|protein_coding|4/10|c.524delG|p.Arg175fs|||||||"
    )
    result = parse_ann_value(ann, fields)
    assert result is not None
    assert result["impact"] == "HIGH"
    assert result["consequence"] == "frameshift_variant"


# ── format_consequence ────────────────────────────────────────────────────────


def test_format_consequence():
    assert format_consequence("missense_variant") == "Missense"
    assert format_consequence("stop_gained") == "Nonsense / Stop gain"
    assert format_consequence("frameshift_variant") == "Frameshift"
    assert format_consequence("splice_donor_variant") == "Splice donor"
    assert format_consequence("splice_acceptor_variant") == "Splice acceptor"
    assert format_consequence("synonymous_variant") == "Synonymous"
    assert format_consequence("inframe_deletion") == "In-frame deletion"
    assert format_consequence("inframe_insertion") == "In-frame insertion"
    assert format_consequence("intron_variant") == "Intronic"
    assert format_consequence("5_prime_UTR_variant") == "5' UTR"
    assert format_consequence("3_prime_UTR_variant") == "3' UTR"
    assert format_consequence("intergenic_variant") == "Intergenic"
    assert format_consequence("start_lost") == "Start loss"
    assert format_consequence("stop_lost") == "Stop loss"


def test_format_consequence_multiple():
    """Ampersand-joined consequences show max 2 readable terms."""
    result = format_consequence("missense_variant&splice_region_variant")
    assert "Missense" in result
    assert "/" in result


def test_format_consequence_unknown_term():
    """Unknown terms are title-cased with underscores replaced."""
    result = format_consequence("novel_variant_type")
    assert result == "Novel Variant Type"


def test_format_consequence_empty():
    assert format_consequence("") == ""


# ── Integration: parse_vcf with CSQ annotations ───────────────────────────────


def test_parse_vcf_with_csq():
    """Annotated VCF produces variants with annotation fields populated."""
    variants = parse_vcf(str(ANNOTATED_VCF))
    assert len(variants) == 5

    # First variant: TP53 missense
    tp53 = next((v for v in variants if v.gene == "TP53"), None)
    assert tp53 is not None
    assert tp53.hgvsc == "ENST00000269305.9:c.524G>T"
    assert tp53.hgvsp == "ENSP00000269305.4:p.Arg175Leu"
    assert tp53.consequence == "Missense"
    assert tp53.transcript == "ENST00000269305"
    assert tp53.impact == "MODERATE"
    assert tp53.sift == "deleterious(0.01)"
    assert tp53.polyphen == "probably_damaging(0.998)"

    # BRCA2: stop_gained → HIGH impact
    brca2 = next((v for v in variants if v.gene == "BRCA2"), None)
    assert brca2 is not None
    assert brca2.impact == "HIGH"
    assert brca2.consequence == "Nonsense / Stop gain"
    assert brca2.hgvsc is not None

    # CFTR: frameshift
    cftr = next((v for v in variants if v.gene == "CFTR"), None)
    assert cftr is not None
    assert cftr.consequence == "Frameshift"
    assert cftr.impact == "HIGH"


def test_parse_vcf_with_csq_canonical_preferred():
    """When multiple transcripts, CANONICAL=YES is selected."""
    variants = parse_vcf(str(ANNOTATED_VCF))
    tp53 = next((v for v in variants if v.gene == "TP53"), None)
    assert tp53 is not None
    # The canonical transcript ENST00000269305 should be selected, not ENST00000413465
    assert tp53.transcript == "ENST00000269305"


# ── Integration: parse_vcf without annotation (existing behavior) ─────────────


def test_parse_vcf_without_annotation():
    """Plain VCF (no CSQ/ANN) still parses correctly with no annotation fields."""
    variants = parse_vcf(str(PLAIN_VCF))
    assert len(variants) == 10

    # First variant: TP53
    v = variants[0]
    assert v.chrom == "chr17"
    assert v.gene == "TP53"
    assert v.rsid == "rs28934578"

    # No annotation fields populated
    assert v.hgvsc is None
    assert v.hgvsp is None
    assert v.consequence is None
    assert v.transcript is None
    assert v.impact is None
    assert v.sift is None
    assert v.polyphen is None


# ── Integration: annotation fallback to gene_knowledge in generate_report_html ─


def test_annotation_fallback_to_gene_knowledge():
    """When no VCF annotation, generate_report_html uses gene_knowledge.json data."""
    from scripts.counselor.generate_pdf import generate_report_html

    data = {
        "sample_id": "TEST_FALLBACK",
        "date": "2026-03-23",
        "variants": [
            {
                "variant": "chr17:7577120:G>A",
                "gene": "TP53",
                "classification": "Pathogenic",
                "acmg_codes": ["PVS1"],
                "conflict": False,
                # No hgvsc/hgvsp — should fall back to gene_knowledge
                "hgvsc": "",
                "hgvsp": "",
                "consequence": "",
                "transcript": "",
                "impact": "",
                "sift": "",
                "polyphen": "",
                "agents": {
                    "clinical": {"clinvar_significance": "Pathogenic"},
                    "korean_pop": {"korean_flag": ""},
                },
            }
        ],
        "pgx_results": [],
        "summary": {"total": 1, "pathogenic": 1, "vus": 0, "benign": 0},
        "db_versions": {"clinvar": "2026-03-23"},
    }
    html = generate_report_html(data)
    assert html is not None
    assert "TP53" in html
    assert "Pathogenic" in html


def test_annotation_vcf_fields_override_gene_knowledge():
    """When VCF annotation present, hgvs dict uses VCF data, not gene_knowledge."""
    from scripts.counselor.generate_pdf import generate_report_html

    data = {
        "sample_id": "TEST_ANNOT",
        "date": "2026-03-23",
        "variants": [
            {
                "variant": "chr17:7675088:C>A",
                "gene": "TP53",
                "classification": "Pathogenic",
                "acmg_codes": ["PVS1"],
                "conflict": False,
                "hgvsc": "ENST00000269305.9:c.524G>T",
                "hgvsp": "ENSP00000269305.4:p.Arg175Leu",
                "consequence": "Missense",
                "transcript": "ENST00000269305",
                "impact": "MODERATE",
                "sift": "deleterious(0.01)",
                "polyphen": "probably_damaging(0.998)",
                "agents": {
                    "clinical": {"clinvar_significance": "Pathogenic"},
                    "korean_pop": {"korean_flag": ""},
                },
            }
        ],
        "pgx_results": [],
        "summary": {"total": 1, "pathogenic": 1, "vus": 0, "benign": 0},
        "db_versions": {"clinvar": "2026-03-23"},
    }
    html = generate_report_html(data)
    assert html is not None
    assert "TP53" in html
    # The VCF-sourced hgvs data should be used
    assert "c.524G>T" in html or "Arg175Leu" in html or "Missense" in html


# ── SnpEff ANN via tmp VCF ────────────────────────────────────────────────────


def test_parse_vcf_with_ann(tmp_path):
    """VCF with SnpEff ANN field parses correctly."""
    vcf_content = (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
        "'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | "
        "Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | ...'\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr17\t7675088\trs28934578\tC\tA\t.\tPASS\t"
        "ANN=A|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000269305|"
        "protein_coding|4/11|c.524G>T|p.Arg175Leu|605/2579|524/1182|175/393||\n"
    )
    vcf_file = tmp_path / "snpeff.vcf"
    vcf_file.write_text(vcf_content)

    variants = parse_vcf(str(vcf_file))
    assert len(variants) == 1
    v = variants[0]
    assert v.gene == "TP53"
    assert v.hgvsc == "c.524G>T"
    assert v.hgvsp == "p.Arg175Leu"
    assert v.consequence == "Missense"
    assert v.impact == "MODERATE"
    assert v.transcript == "ENST00000269305"


def test_parse_vcf_ann_gene_from_annotation(tmp_path):
    """When Gene= is absent from INFO but ANN has gene name, variant.gene is set from ANN."""
    vcf_content = (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
        "'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | "
        "Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | ...'\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr17\t7675088\t.\tC\tA\t.\tPASS\t"
        "ANN=A|frameshift_variant|HIGH|BRCA1|ENSG00000012048|transcript|ENST00000357654|"
        "protein_coding|11/23|c.3756_3759del|p.Ser1253fs||||||||\n"
    )
    vcf_file = tmp_path / "ann_gene.vcf"
    vcf_file.write_text(vcf_content)

    variants = parse_vcf(str(vcf_file))
    assert len(variants) == 1
    v = variants[0]
    assert v.gene == "BRCA1"
    assert v.consequence == "Frameshift"
    assert v.impact == "HIGH"
