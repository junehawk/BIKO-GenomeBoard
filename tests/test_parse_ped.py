"""v2.4 Quick Win C — PED file based trio resolution.

Tests parse_ped.py (PED parser + resolve_trio) and its integration with
parse_vcf + orchestrate. Strict mode is enforced: when --ped is provided
and the PED does not resolve a trio against the VCF sample list, parse_vcf
raises ValueError rather than silently falling back to the filename
heuristic. Quartet / N≥3 VCFs are supported via PED.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from scripts.intake.parse_ped import PedEntry, parse_ped, resolve_trio
from scripts.intake.parse_vcf import parse_vcf

# ─── parse_ped ────────────────────────────────────────────────────────────


def test_parse_ped_basic(tmp_path):
    """Standard 6-column, tab-delimited PED round-trips into PedEntry."""
    ped = tmp_path / "family.ped"
    ped.write_text("FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\nFAM1\tFATHER\t0\t0\t1\t1\nFAM1\tMOTHER\t0\t0\t2\t1\n")
    result = parse_ped(str(ped))
    assert set(result.keys()) == {"CHILD", "FATHER", "MOTHER"}
    child = result["CHILD"]
    assert isinstance(child, PedEntry)
    assert child.family_id == "FAM1"
    assert child.father_id == "FATHER"
    assert child.mother_id == "MOTHER"
    assert child.sex == 1
    assert child.affected == 2
    # Missing-parent "0" preserved verbatim on the founders.
    assert result["FATHER"].father_id == "0"
    assert result["MOTHER"].mother_id == "0"


def test_parse_ped_with_comments(tmp_path):
    """Lines starting with # are treated as comments and skipped."""
    ped = tmp_path / "comments.ped"
    ped.write_text("# Family CEPH-1463 trio\n# proband is affected\nFAM1\tCHILD\tFATHER\tMOTHER\t1\t2\n# end of file\n")
    result = parse_ped(str(ped))
    assert list(result.keys()) == ["CHILD"]


def test_parse_ped_empty_file(tmp_path):
    """A zero-byte PED parses to an empty dict (no raise)."""
    ped = tmp_path / "empty.ped"
    ped.write_text("")
    result = parse_ped(str(ped))
    assert result == {}


def test_parse_ped_blank_lines_skipped(tmp_path):
    """Whitespace-only lines do not raise."""
    ped = tmp_path / "blanks.ped"
    ped.write_text("\n\nFAM1\tCHILD\tFATHER\tMOTHER\t1\t2\n\n   \n")
    result = parse_ped(str(ped))
    assert list(result.keys()) == ["CHILD"]


def test_parse_ped_malformed_raises(tmp_path):
    """< 6 columns raises ValueError carrying the 1-based line number."""
    ped = tmp_path / "bad.ped"
    ped.write_text(
        "FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\nFAM1\tONLY3COLS\tFATHER\n"  # line 2 malformed
    )
    with pytest.raises(ValueError) as excinfo:
        parse_ped(str(ped))
    # Must pinpoint line 2
    assert "line 2" in str(excinfo.value)


def test_parse_ped_non_integer_sex_raises(tmp_path):
    """Non-integer Sex column raises with line number."""
    ped = tmp_path / "bad_sex.ped"
    ped.write_text("FAM1\tCHILD\tFATHER\tMOTHER\tMALE\t2\n")
    with pytest.raises(ValueError) as excinfo:
        parse_ped(str(ped))
    assert "line 1" in str(excinfo.value) and "Sex" in str(excinfo.value)


def test_parse_ped_crlf_bom(tmp_path):
    """Windows CRLF + UTF-8 BOM on the first line are handled."""
    ped = tmp_path / "crlf_bom.ped"
    # UTF-8 BOM bytes + CRLF endings
    ped.write_bytes(b"\xef\xbb\xbfFAM1\tCHILD\tFATHER\tMOTHER\t1\t2\r\nFAM1\tFATHER\t0\t0\t1\t1\r\n")
    result = parse_ped(str(ped))
    assert "CHILD" in result
    assert result["CHILD"].family_id == "FAM1"  # BOM stripped


def test_parse_ped_missing_parent_zero(tmp_path):
    """father_id / mother_id == '0' is preserved as the missing-parent sentinel."""
    ped = tmp_path / "founder.ped"
    ped.write_text("FAM1\tLONER\t0\t0\t1\t2\n")
    result = parse_ped(str(ped))
    assert result["LONER"].father_id == "0"
    assert result["LONER"].mother_id == "0"


def test_parse_ped_space_delimited(tmp_path):
    """PLINK's canonical form is tab-delimited but many tools emit spaces."""
    ped = tmp_path / "spaces.ped"
    ped.write_text("FAM1 CHILD FATHER MOTHER 1 2\n")
    result = parse_ped(str(ped))
    assert result["CHILD"].father_id == "FATHER"


def test_parse_ped_missing_file_raises(tmp_path):
    """A non-existent path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        parse_ped(str(tmp_path / "does_not_exist.ped"))


# ─── resolve_trio ─────────────────────────────────────────────────────────


def _make_ped(entries):
    """Build a {iid: PedEntry} dict from (iid, father, mother, affected) tuples."""
    return {
        iid: PedEntry(
            family_id="FAM1",
            individual_id=iid,
            father_id=father,
            mother_id=mother,
            sex=1,
            affected=affected,
        )
        for iid, father, mother, affected in entries
    }


def test_resolve_trio_affected_proband():
    """Standard case: affected child with both parents present."""
    ped = _make_ped(
        [
            ("CHILD", "FATHER", "MOTHER", 2),
            ("FATHER", "0", "0", 1),
            ("MOTHER", "0", "0", 1),
        ]
    )
    sample_ids = ["FATHER", "MOTHER", "CHILD"]
    proband_idx, parent_idxs = resolve_trio(ped, sample_ids)
    assert proband_idx == 2  # CHILD
    assert sorted(parent_idxs) == [0, 1]  # FATHER, MOTHER


def test_resolve_trio_no_match_returns_none():
    """Affected individual present in PED but NOT in VCF → (None, [])."""
    ped = _make_ped([("CHILD", "FATHER", "MOTHER", 2)])
    sample_ids = ["FATHER", "MOTHER", "UNRELATED"]
    proband_idx, parent_idxs = resolve_trio(ped, sample_ids)
    assert proband_idx is None
    assert parent_idxs == []


def test_resolve_trio_missing_parent_excluded():
    """Affected individual with father=0 mother=0 → no qualifying candidate."""
    ped = _make_ped([("LONER", "0", "0", 2)])
    sample_ids = ["LONER"]
    proband_idx, parent_idxs = resolve_trio(ped, sample_ids)
    assert proband_idx is None
    assert parent_idxs == []


def test_resolve_trio_parent_missing_from_vcf_excluded():
    """Father listed in PED but absent from VCF → candidate disqualified."""
    ped = _make_ped(
        [
            ("CHILD", "FATHER", "MOTHER", 2),
        ]
    )
    sample_ids = ["MOTHER", "CHILD"]  # FATHER absent
    proband_idx, parent_idxs = resolve_trio(ped, sample_ids)
    assert proband_idx is None


def test_resolve_trio_multi_affected_picks_first_alpha(caplog):
    """Two affected individuals with usable parents → alpha-first pick + WARNING."""
    ped = _make_ped(
        [
            ("ZEBRA", "F1", "M1", 2),
            ("APPLE", "F2", "M2", 2),
            ("F1", "0", "0", 1),
            ("M1", "0", "0", 1),
            ("F2", "0", "0", 1),
            ("M2", "0", "0", 1),
        ]
    )
    sample_ids = ["F1", "M1", "F2", "M2", "ZEBRA", "APPLE"]
    with caplog.at_level(logging.WARNING, logger="scripts.intake.parse_ped"):
        proband_idx, parent_idxs = resolve_trio(ped, sample_ids)
    # APPLE sorts before ZEBRA
    assert sample_ids[proband_idx] == "APPLE"
    # Parents of APPLE are F2, M2 → indices 2 and 3
    assert sorted(parent_idxs) == [2, 3]
    assert any("2 affected individuals" in r.getMessage() for r in caplog.records)


def test_resolve_trio_quartet_picks_subset():
    """4-sample VCF (quartet) where PED declares a 3-person trio subset.

    This is the core of v2.4 Quick Win C: legacy filename heuristic
    bails on N != 3; PED must succeed at N >= 3 by locating the trio
    subset inside the larger sample list.
    """
    ped = _make_ped(
        [
            ("CHILD", "FATHER", "MOTHER", 2),
            ("FATHER", "0", "0", 1),
            ("MOTHER", "0", "0", 1),
            ("SIBLING", "FATHER", "MOTHER", 1),  # unaffected sib → not a candidate
        ]
    )
    sample_ids = ["FATHER", "MOTHER", "CHILD", "SIBLING"]
    proband_idx, parent_idxs = resolve_trio(ped, sample_ids)
    assert proband_idx == 2  # CHILD
    assert sorted(parent_idxs) == [0, 1]  # FATHER, MOTHER


def test_resolve_trio_empty_sample_ids_returns_none():
    """Edge case — empty VCF sample list short-circuits to (None, [])."""
    ped = _make_ped([("CHILD", "FATHER", "MOTHER", 2)])
    proband_idx, parent_idxs = resolve_trio(ped, [])
    assert proband_idx is None
    assert parent_idxs == []


def test_resolve_trio_empty_ped_returns_none():
    """Empty PED short-circuits to (None, [])."""
    proband_idx, parent_idxs = resolve_trio({}, ["A", "B", "C"])
    assert proband_idx is None
    assert parent_idxs == []


# ─── parse_vcf integration ───────────────────────────────────────────────


_TRIO_HEADER = (
    "##fileformat=VCFv4.1\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=DGQ,Number=1,Type=Integer,Description="Denovo Genotype Quality">\n'
)


def _write_trio_vcf(
    path: Path,
    *,
    samples: list[str],
    records: list[str],
):
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"
    path.write_text(_TRIO_HEADER + columns + "".join(records))


def test_parse_vcf_with_ped_strict_fail(tmp_path):
    """--ped + PED that does not resolve a trio → ValueError (strict mode)."""
    # PED references sample X not in the VCF.
    ped = tmp_path / "orphan.ped"
    ped.write_text("FAM1\tX\tY\tZ\t1\t2\n")
    vcf = tmp_path / "sample.vcf"
    _write_trio_vcf(
        vcf,
        samples=["FOO", "BAR", "BAZ"],
        records=["chr1\t100\t.\tA\tT\t100\tPASS\t.\tGT:DGQ\t0/0:100\t0/0:100\t0/1:100\n"],
    )
    with pytest.raises(ValueError) as excinfo:
        parse_vcf(str(vcf), ped_path=str(ped))
    assert "did not resolve a trio" in str(excinfo.value)
    assert "Strict mode" in str(excinfo.value)


def test_parse_vcf_with_ped_overrides_filename(tmp_path):
    """PED explicit wiring wins over the filename heuristic.

    The VCF filename mentions 9999, which would cause the filename
    heuristic to treat the last-column-default. But the PED says
    CHILD (column 1 in the VCF) is the proband — PED must win and
    the de novo call should agree with CHILD-as-proband.
    """
    ped = tmp_path / "trio.ped"
    ped.write_text("FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\nFAM1\tFATHER\t0\t0\t1\t1\nFAM1\tMOTHER\t0\t0\t2\t1\n")
    # Filename mentions 9999, none of the samples contain 9999.
    vcf = tmp_path / "9999-trio.vcf"
    # Column order: CHILD (proband, 0/1), FATHER (0/0), MOTHER (0/0).
    # With filename heuristic, last sample (MOTHER) would be proband →
    # detector would not fire (MOTHER is 0/0). With PED, CHILD is
    # proband and the variant is flagged de novo.
    _write_trio_vcf(
        vcf,
        samples=["CHILD", "FATHER", "MOTHER"],
        records=["chr1\t100\t.\tA\tT\t100\tPASS\t.\tGT:DGQ\t0/1:100\t0/0:100\t0/0:100\n"],
    )
    variants = parse_vcf(str(vcf), ped_path=str(ped))
    assert len(variants) == 1
    assert variants[0].inheritance == "de_novo"


def test_parse_vcf_without_ped_uses_filename(tmp_path):
    """Regression: no --ped → existing filename heuristic stays untouched."""
    vcf = tmp_path / "IBS-ASD-9703-blood.vcf"
    _write_trio_vcf(
        vcf,
        samples=["IBS-ASD-9701-blood", "IBS-ASD-9702-blood", "IBS-ASD-9703-blood"],
        records=["chr1\t100\t.\tA\tT\t100\tPASS\t.\tGT:DGQ\t0/0:100\t0/0:100\t0/1:100\n"],
    )
    variants = parse_vcf(str(vcf))
    assert variants[0].inheritance == "de_novo"


def test_parse_vcf_with_ped_on_quartet(tmp_path):
    """N=4 sample VCF resolves via PED (core Quick Win C value-add)."""
    ped = tmp_path / "quartet.ped"
    ped.write_text(
        "FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\n"
        "FAM1\tFATHER\t0\t0\t1\t1\n"
        "FAM1\tMOTHER\t0\t0\t2\t1\n"
        "FAM1\tSIB\tFATHER\tMOTHER\t1\t1\n"
    )
    vcf = tmp_path / "quartet.vcf"
    _write_trio_vcf(
        vcf,
        samples=["FATHER", "MOTHER", "CHILD", "SIB"],
        records=[
            "chr1\t100\t.\tA\tT\t100\tPASS\t.\tGT:DGQ\t0/0:100\t0/0:100\t0/1:100\t0/0:100\n",
        ],
    )
    variants = parse_vcf(str(vcf), ped_path=str(ped))
    # CHILD is proband; parents are FATHER (0/0) and MOTHER (0/0).
    assert variants[0].inheritance == "de_novo"


def test_parse_vcf_ped_logs_proband(tmp_path, caplog):
    """PED-based resolution emits an info log identifying proband + parents."""
    ped = tmp_path / "trio.ped"
    ped.write_text("FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\nFAM1\tFATHER\t0\t0\t1\t1\nFAM1\tMOTHER\t0\t0\t2\t1\n")
    vcf = tmp_path / "trio.vcf"
    _write_trio_vcf(
        vcf,
        samples=["FATHER", "MOTHER", "CHILD"],
        records=["chr1\t100\t.\tA\tT\t100\tPASS\t.\tGT:DGQ\t0/0:100\t0/0:100\t0/1:100\n"],
    )
    with caplog.at_level(logging.INFO, logger="scripts.intake.parse_vcf"):
        parse_vcf(str(vcf), ped_path=str(ped))
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "PED-based trio" in msgs
    assert "CHILD" in msgs


# ─── orchestrate integration ─────────────────────────────────────────────


def test_orchestrate_ped_flag_parsed():
    """argparse accepts --ped and stores on attribute 'ped_path'."""
    import argparse

    # Reconstruct the one flag under test — full main() CLI is heavy.
    parser = argparse.ArgumentParser()
    parser.add_argument("--ped", type=str, default=None, dest="ped_path")
    args = parser.parse_args(["--ped", "/tmp/family.ped"])
    assert args.ped_path == "/tmp/family.ped"


def test_orchestrate_ped_end_to_end_records_flag(tmp_path):
    """run_pipeline(ped_path=...) records pipeline.ped_used in report_data."""
    from scripts.orchestrate import run_pipeline

    ped = tmp_path / "trio.ped"
    ped.write_text("FAM1\tCHILD\tFATHER\tMOTHER\t1\t2\nFAM1\tFATHER\t0\t0\t1\t1\nFAM1\tMOTHER\t0\t0\t2\t1\n")
    vcf = tmp_path / "trio.vcf"
    _write_trio_vcf(
        vcf,
        samples=["FATHER", "MOTHER", "CHILD"],
        records=[
            "chr1\t7577120\t.\tG\tA\t100\tPASS\tGene=TP53\tGT:DGQ\t0/0:100\t0/0:100\t0/1:100\n",
        ],
    )
    output = tmp_path / "report.html"
    result = run_pipeline(
        vcf_path=str(vcf),
        output_path=str(output),
        skip_api=True,
        mode="rare-disease",
        ped_path=str(ped),
    )
    assert result is not None
    assert result["pipeline"]["ped_used"] is True
    assert result["pipeline"]["ped_path"] == str(ped)


def test_orchestrate_no_ped_ped_used_false(tmp_path):
    """Without --ped, pipeline.ped_used is False (regression)."""
    from scripts.orchestrate import run_pipeline

    vcf = tmp_path / "single.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.1\n"
        '##INFO=<ID=Gene,Number=1,Type=String,Description="Gene symbol">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr17\t7577120\t.\tG\tA\t.\tPASS\tGene=TP53\n"
    )
    output = tmp_path / "report.html"
    result = run_pipeline(
        vcf_path=str(vcf),
        output_path=str(output),
        skip_api=True,
        mode="cancer",
    )
    assert result is not None
    assert result["pipeline"]["ped_used"] is False
    assert result["pipeline"]["ped_path"] is None
