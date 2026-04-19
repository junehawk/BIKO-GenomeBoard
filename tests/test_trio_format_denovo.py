"""Regression tests for v2.3-T9 (trio FORMAT/GT de novo detection) and
v2.3-T10 (annotation-missing fail-loud warning).

Ticket 1 shipped INFO-level DN / DENOVO / CONFIRMED_DN / PS2 parsing
under the assumption that an upstream trio caller would tag the INFO
field. End-to-end validation against the ASD triodenovo fixtures on
2026-04-16 exposed that the real-world caller (triodenovo 0.05) emits
the de novo signal in FORMAT DQ / DGQ plus 3-sample GT, NOT in INFO.
This module exercises the FORMAT fallback and the new fail-loud path.
"""

from __future__ import annotations

import logging

from scripts.intake.parse_vcf import parse_vcf

_HEADER = (
    "##fileformat=VCFv4.1\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=DQ,Number=1,Type=Float,Description="Denovo Quality">\n'
    '##FORMAT=<ID=DGQ,Number=1,Type=Integer,Description="Denovo Genotype Quality">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n'
)


def _write_trio_vcf(
    tmp_path, filename, records, father_id="SAMPLE-9701", mother_id="SAMPLE-9702", child_id="SAMPLE-9703"
):
    path = tmp_path / filename
    lines = [_HEADER]
    lines.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{father_id}\t{mother_id}\t{child_id}\n")
    for rec in records:
        lines.append(rec)
    path.write_text("".join(lines))
    return path


# ─── T9 core — trio GT de novo detection ──────────────────────────────────


def test_trio_format_denovo_parents_refhom_child_het_is_detected(tmp_path, caplog):
    """Father 0/0 + Mother 0/0 + Child 0/1 (with passing DGQ) → de novo."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-denovo.vcf",
        [
            "chr1\t100\t.\tA\tT\t100\tPASS\tDP=30\tGT:DQ:DGQ:DP\t0/0:12.5:100:40\t0/0:12.5:100:42\t0/1:12.5:100:35\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert len(variants) == 1
    assert variants[0].inheritance == "de_novo"
    # Trio GT alone never promotes to confirmed (spec Q1.2).
    assert variants[0].confirmed_denovo is False


def test_trio_format_denovo_low_dgq_is_rejected(tmp_path):
    """DGQ below triodenovo default threshold (7) → not flagged as de novo."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-noise.vcf",
        [
            "chr1\t200\t.\tA\tT\t100\tPASS\tDP=30\tGT:DQ:DGQ:DP\t0/0:5:5:40\t0/0:5:5:42\t0/1:5:5:35\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance is None


def test_trio_format_denovo_inherited_not_flagged(tmp_path):
    """Father 0/1 + Mother 0/0 + Child 0/1 → inherited from father, not de novo."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-inh.vcf",
        [
            "chr1\t300\t.\tA\tT\t100\tPASS\tDP=30\tGT:DQ:DGQ:DP\t0/1:12:100:40\t0/0:12:100:42\t0/1:12:100:35\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance is None


def test_trio_format_denovo_parent_missing_gt_not_flagged(tmp_path):
    """Missing parental GT (./.) is conservatively treated as 'not evidence'."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-miss.vcf",
        [
            "chr1\t400\t.\tA\tT\t100\tPASS\tDP=30\tGT:DQ:DGQ:DP\t./.:12:100:40\t0/0:12:100:42\t0/1:12:100:35\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance is None


def test_trio_format_denovo_no_dgq_still_works(tmp_path):
    """When FORMAT has no DGQ the quality gate is skipped — label still set."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-nodgq.vcf",
        [
            "chr1\t500\t.\tA\tT\t100\tPASS\tDP=30\tGT:DP\t0/0:40\t0/0:42\t0/1:35\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance == "de_novo"


def test_trio_format_denovo_phased_bar_separator(tmp_path):
    """Phased genotypes (0|0, 0|1) are also recognised."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-phased.vcf",
        [
            "chr1\t600\t.\tA\tT\t100\tPASS\tDP=30\tGT:DGQ\t0|0:100\t0|0:100\t0|1:100\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance == "de_novo"


def test_trio_format_denovo_polymutt_nucleotide_gt(tmp_path):
    """polymutt-modified VCFs emit GT as actual nucleotides (G/G, A/G)
    instead of the VCF-standard integer form (0/0, 0/1). Without GT
    normalisation the trio detector silently fails on every polymutt
    row. Regression for the ASD-9703 validation finding on 2026-04-16.
    """
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-polymutt.vcf",
        [
            # REF=G ALT=A, parents G/G G/G, child A/G — nucleotide form
            "chr1\t101268\t.\tG\tA\t100\tPASS\tDP=30\tGT:DGQ\tG/G:100\tG/G:100\tA/G:100\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance == "de_novo", (
        "polymutt nucleotide GT (G/G + G/G + A/G) should be normalised and recognised as a de novo call"
    )


def test_trio_format_denovo_polymutt_phased_nucleotide(tmp_path):
    """Phased nucleotide GT (G|G, A|G) must also normalise."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-polymutt-phased.vcf",
        [
            "chr1\t101500\t.\tG\tA\t100\tPASS\tDP=30\tGT:DGQ\tG|G:100\tG|G:100\tA|G:100\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance == "de_novo"


def test_trio_format_denovo_polymutt_inherited_not_flagged(tmp_path):
    """Father nucleotide A/G (heterozygous) → child inherited, not de novo."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-polymutt-inh.vcf",
        [
            "chr1\t101700\t.\tG\tA\t100\tPASS\tDP=30\tGT:DGQ\tA/G:100\tG/G:100\tA/G:100\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance is None


def test_trio_format_denovo_homozygous_child(tmp_path):
    """Father 0/0 + Mother 0/0 + Child 1/1 is still a de novo call (rare but
    valid — two independent de novo events on both alleles or uniparental
    disomy). We accept it; downstream reviewer decides the interpretation."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-hom.vcf",
        [
            "chr1\t700\t.\tA\tT\t100\tPASS\tDP=30\tGT:DGQ\t0/0:100\t0/0:100\t1/1:100\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance == "de_novo"


# ─── Proband-index heuristic ─────────────────────────────────────────────


def test_trio_proband_detected_by_filename_match(tmp_path):
    """Filename contains 9703 → sample-9703 is the proband even if listed last.

    The heuristic must not rely on column order — verified here by placing
    the proband in the MIDDLE column and checking that a father-at-col0 /
    mother-at-col2 arrangement is still correctly identified.
    """
    path = tmp_path / "IBS-ASD-9703-blood-wgs-ILLUMINA.fmarked_denovo.vcf"
    path.write_text(
        _HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        "\tIBS-ASD-9701-blood\tIBS-ASD-9703-blood\tIBS-ASD-9702-blood\n"
        + "chr1\t900\t.\tA\tT\t100\tPASS\tDP=30\tGT:DGQ\t0/0:100\t0/1:100\t0/0:100\n"
    )
    variants = parse_vcf(str(path))
    assert variants[0].inheritance == "de_novo"


def test_trio_proband_fallback_to_last_with_warning(tmp_path, caplog):
    """No sample ID matches the filename → fallback to last column + warning."""
    path = tmp_path / "generic.vcf"
    path.write_text(
        _HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tfoo\tbar\tbaz\n"
        + "chr1\t1000\t.\tA\tT\t100\tPASS\tDP=30\tGT:DGQ\t0/0:100\t0/0:100\t0/1:100\n"
    )
    with caplog.at_level(logging.WARNING, logger="scripts.intake.parse_vcf"):
        variants = parse_vcf(str(path))
    assert variants[0].inheritance == "de_novo"
    assert any("falling back" in r.getMessage().lower() for r in caplog.records)


def test_non_trio_vcf_skips_trio_detection(tmp_path):
    """Singleton VCF (1 sample) must not trigger trio detection."""
    path = tmp_path / "singleton.vcf"
    path.write_text(
        _HEADER
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband\n"
        + "chr1\t1100\t.\tA\tT\t100\tPASS\tDP=30\tGT:DGQ\t0/1:100\n"
    )
    variants = parse_vcf(str(path))
    assert variants[0].inheritance is None


# ─── INFO-flag precedence ────────────────────────────────────────────────


def test_info_confirmed_dn_flag_wins_over_trio_gt(tmp_path):
    """A CONFIRMED_DN INFO flag promotes to confirmed_de_novo even if the
    trio GT pattern alone would have fired PM6-only."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-confirmed.vcf",
        [
            "chr1\t1200\t.\tA\tT\t100\tPASS\tCONFIRMED_DN=1\tGT:DGQ\t0/0:100\t0/0:100\t0/1:100\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance == "confirmed_de_novo"
    assert variants[0].confirmed_denovo is True


def test_trio_gt_does_not_downgrade_info_confirmed(tmp_path):
    """Trio FORMAT path must never overwrite an upstream 'confirmed' label."""
    vcf = _write_trio_vcf(
        tmp_path,
        "SAMPLE-9703-confirmed-gt-inherited.vcf",
        [
            # INFO says confirmed; GT disagrees (father 0/1 looks inherited).
            # Spec: INFO wins — this is a case where the PED-level identity
            # check overrode the callset genotype. We trust INFO.
            "chr1\t1300\t.\tA\tT\t100\tPASS\tCONFIRMED_DN=1\tGT:DGQ\t0/1:100\t0/0:100\t0/1:100\n",
        ],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].confirmed_denovo is True


# ─── T10 annotation fail-loud ────────────────────────────────────────────


def test_t10_warns_when_vcf_has_no_annotation(tmp_path, caplog):
    """Raw call-set VCF with no Gene=, CSQ=, or ANN= → visible warning."""
    path = tmp_path / "unannotated.vcf"
    path.write_text(
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t2000\t.\tA\tT\t100\tPASS\tDP=30\n"
        "chr2\t3000\t.\tG\tC\t100\tPASS\tDP=40\n"
    )
    with caplog.at_level(logging.WARNING, logger="scripts.intake.parse_vcf"):
        parse_vcf(str(path))
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "NONE carry annotation" in msgs
    assert "Run VEP" in msgs or "SnpEff" in msgs


def test_t10_silent_when_vcf_has_gene_info(tmp_path, caplog):
    """Annotated VCF must not emit the fail-loud warning."""
    path = tmp_path / "annotated.vcf"
    path.write_text(
        "##fileformat=VCFv4.1\n"
        '##INFO=<ID=Gene,Number=1,Type=String,Description="Gene symbol">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr17\t7577120\trs28934578\tG\tA\t.\tPASS\tGene=TP53\n"
    )
    with caplog.at_level(logging.WARNING, logger="scripts.intake.parse_vcf"):
        parse_vcf(str(path))
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "NONE carry annotation" not in msgs


def test_t10_warns_on_low_gene_coverage(tmp_path, caplog):
    """If <10% of variants carry a gene annotation, surface a partial-cover
    warning — mirrors the full-miss path but tells the reviewer the input
    is partially annotated."""
    path = tmp_path / "partial.vcf"
    lines = [
        "##fileformat=VCFv4.1\n",
        '##INFO=<ID=Gene,Number=1,Type=String,Description="Gene symbol">\n',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        "chr1\t100\t.\tA\tT\t.\tPASS\tGene=TP53\n",
    ]
    # Pad with 20 un-gene'd variants so 1/21 ≈ 4.8% is annotated.
    for i in range(20):
        lines.append(f"chr1\t{200 + i}\t.\tA\tT\t.\tPASS\tDP=30\n")
    path.write_text("".join(lines))
    with caplog.at_level(logging.WARNING, logger="scripts.intake.parse_vcf"):
        parse_vcf(str(path))
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "only" in msgs.lower() and "gene annotation" in msgs.lower()
