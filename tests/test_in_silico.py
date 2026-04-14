"""Tests for in silico prediction scoring and PP3/BP4 evidence generation."""

from pathlib import Path

import pytest

from scripts.classification.in_silico import (
    InSilicoScores,
    format_scores_for_display,
    generate_pp3_bp4,
    parse_in_silico_from_csq,
)

# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

DEMO_VCF = Path(__file__).parent.parent / "data" / "sample_vcf" / "demo_with_in_silico.vcf"


def _csq_dict(**kwargs: str) -> dict:
    """Build a lowercased CSQ field dict with only the given keys populated."""
    return {k: v for k, v in kwargs.items()}


# ═══════════════════════════════════════════════════════════════════════════
# InSilicoScores — parsing from CSQ
# ═══════════════════════════════════════════════════════════════════════════


class TestParseFromCsq:
    """Tests for parse_in_silico_from_csq."""

    def test_parse_from_csq_full(self):
        """All in silico fields present → all attributes populated."""
        csq = _csq_dict(
            revel_score="0.95",
            cadd_phred="35.0",
            am_class="likely_pathogenic",
            am_pathogenicity="0.97",
            spliceai_pred_ds_ag="0.01",
            spliceai_pred_ds_al="0.80",
            spliceai_pred_ds_dg="0.02",
            spliceai_pred_ds_dl="0.00",
            sift="deleterious(0.01)",
            polyphen="probably_damaging(0.998)",
        )
        scores = parse_in_silico_from_csq(csq)

        assert scores.revel == pytest.approx(0.95)
        assert scores.cadd_phred == pytest.approx(35.0)
        assert scores.alphamissense_score == pytest.approx(0.97)
        assert scores.alphamissense_class == "likely_pathogenic"
        assert scores.spliceai_ds_ag == pytest.approx(0.01)
        assert scores.spliceai_ds_al == pytest.approx(0.80)
        assert scores.spliceai_ds_dg == pytest.approx(0.02)
        assert scores.spliceai_ds_dl == pytest.approx(0.00)
        assert scores.spliceai_max == pytest.approx(0.80)
        assert scores.sift == "deleterious(0.01)"
        assert scores.polyphen == "probably_damaging(0.998)"

    def test_parse_from_csq_partial(self):
        """Only some fields present → others are None."""
        csq = _csq_dict(revel_score="0.75", cadd_phred="28.0")
        scores = parse_in_silico_from_csq(csq)

        assert scores.revel == pytest.approx(0.75)
        assert scores.cadd_phred == pytest.approx(28.0)
        assert scores.alphamissense_score is None
        assert scores.alphamissense_class is None
        assert scores.spliceai_max is None
        assert scores.spliceai_ds_ag is None
        assert scores.sift is None
        assert scores.polyphen is None

    def test_parse_from_csq_empty(self):
        """No in silico fields at all → all None."""
        scores = parse_in_silico_from_csq({})

        assert scores.revel is None
        assert scores.cadd_phred is None
        assert scores.alphamissense_score is None
        assert scores.alphamissense_class is None
        assert scores.spliceai_max is None
        assert scores.sift is None
        assert scores.polyphen is None

    def test_parse_missing_vep_sentinels(self):
        """VEP missing-value sentinels (., -, empty) → None."""
        csq = _csq_dict(
            revel_score=".",
            cadd_phred="-",
            am_class="",
            am_pathogenicity="NA",
            spliceai_pred_ds_ag="None",
        )
        scores = parse_in_silico_from_csq(csq)

        assert scores.revel is None
        assert scores.cadd_phred is None
        assert scores.alphamissense_class is None
        assert scores.alphamissense_score is None
        assert scores.spliceai_ds_ag is None

    def test_parse_invalid_float(self):
        """Non-numeric strings for float fields → None."""
        csq = _csq_dict(revel_score="not_a_number", cadd_phred="abc")
        scores = parse_in_silico_from_csq(csq)

        assert scores.revel is None
        assert scores.cadd_phred is None

    def test_spliceai_max_calculation(self):
        """spliceai_max = max of four delta scores."""
        csq = _csq_dict(
            spliceai_pred_ds_ag="0.10",
            spliceai_pred_ds_al="0.30",
            spliceai_pred_ds_dg="0.55",
            spliceai_pred_ds_dl="0.20",
        )
        scores = parse_in_silico_from_csq(csq)

        assert scores.spliceai_max == pytest.approx(0.55)

    def test_spliceai_max_partial(self):
        """spliceai_max calculated from available scores only."""
        csq = _csq_dict(
            spliceai_pred_ds_ag="0.10",
            spliceai_pred_ds_al=".",
            spliceai_pred_ds_dg="0.05",
        )
        scores = parse_in_silico_from_csq(csq)

        assert scores.spliceai_max == pytest.approx(0.10)
        assert scores.spliceai_ds_al is None

    def test_parse_alternative_revel_key(self):
        """Accept 'revel' as well as 'revel_score' CSQ field name."""
        csq = _csq_dict(revel="0.88")
        scores = parse_in_silico_from_csq(csq)
        assert scores.revel == pytest.approx(0.88)

    def test_parse_negative_cadd(self):
        """Negative CADD phred scores are valid (raw CADD can be negative)."""
        csq = _csq_dict(cadd_phred="-2.5")
        scores = parse_in_silico_from_csq(csq)
        assert scores.cadd_phred == pytest.approx(-2.5)


# ═══════════════════════════════════════════════════════════════════════════
# PP3 / BP4 evidence generation
# ═══════════════════════════════════════════════════════════════════════════


class TestPP3BP4Generation:
    """Tests for generate_pp3_bp4."""

    # --- REVEL-based evidence ---

    def test_pp3_strong_revel(self):
        """REVEL >= 0.932 → PP3_Strong."""
        scores = InSilicoScores(revel=0.95)
        result = generate_pp3_bp4(scores)
        assert result == ["PP3_Strong"]

    def test_pp3_strong_revel_boundary(self):
        """REVEL == 0.932 (exactly at threshold) → PP3_Strong."""
        scores = InSilicoScores(revel=0.932)
        result = generate_pp3_bp4(scores)
        assert result == ["PP3_Strong"]

    def test_pp3_moderate_revel(self):
        """REVEL >= 0.644 (below Strong) → PP3_Moderate."""
        scores = InSilicoScores(revel=0.70)
        result = generate_pp3_bp4(scores)
        assert result == ["PP3_Moderate"]

    def test_pp3_moderate_revel_boundary(self):
        """REVEL == 0.644 (exactly at threshold) → PP3_Moderate."""
        scores = InSilicoScores(revel=0.644)
        result = generate_pp3_bp4(scores)
        assert result == ["PP3_Moderate"]

    def test_bp4_strong_revel(self):
        """REVEL <= 0.016 → BP4_Strong."""
        scores = InSilicoScores(revel=0.01)
        result = generate_pp3_bp4(scores)
        assert result == ["BP4_Strong"]

    def test_bp4_strong_revel_boundary(self):
        """REVEL == 0.016 (exactly at threshold) → BP4_Strong."""
        scores = InSilicoScores(revel=0.016)
        result = generate_pp3_bp4(scores)
        assert result == ["BP4_Strong"]

    def test_bp4_moderate_revel(self):
        """REVEL <= 0.183 (above Strong BP4) → BP4_Moderate."""
        scores = InSilicoScores(revel=0.10)
        result = generate_pp3_bp4(scores)
        assert result == ["BP4_Moderate"]

    def test_bp4_moderate_revel_boundary(self):
        """REVEL == 0.183 (exactly at threshold) → BP4_Moderate."""
        scores = InSilicoScores(revel=0.183)
        result = generate_pp3_bp4(scores)
        assert result == ["BP4_Moderate"]

    def test_revel_indeterminate(self):
        """REVEL in the middle (0.3) — no evidence."""
        scores = InSilicoScores(revel=0.30)
        result = generate_pp3_bp4(scores)
        assert result == []

    # --- SpliceAI-based evidence ---

    def test_pp3_splice_strong(self):
        """SpliceAI max >= 0.5 → PP3_Strong."""
        scores = InSilicoScores(
            spliceai_ds_al=0.85,
            spliceai_ds_ag=0.10,
            spliceai_max=0.85,
        )
        result = generate_pp3_bp4(scores)
        assert result == ["PP3_Strong"]

    def test_pp3_splice_strong_boundary(self):
        """SpliceAI max == 0.5 (exactly at threshold) → PP3_Strong."""
        scores = InSilicoScores(
            spliceai_ds_dg=0.50,
            spliceai_max=0.50,
        )
        result = generate_pp3_bp4(scores)
        assert result == ["PP3_Strong"]

    def test_pp3_splice_moderate(self):
        """SpliceAI max >= 0.2 (below Strong) → PP3_Moderate."""
        scores = InSilicoScores(
            spliceai_ds_ag=0.30,
            spliceai_max=0.30,
        )
        result = generate_pp3_bp4(scores)
        assert result == ["PP3_Moderate"]

    def test_bp4_splice(self):
        """SpliceAI max < 0.1 → BP4_Supporting."""
        scores = InSilicoScores(
            spliceai_ds_ag=0.02,
            spliceai_ds_al=0.01,
            spliceai_ds_dg=0.03,
            spliceai_ds_dl=0.01,
            spliceai_max=0.03,
        )
        result = generate_pp3_bp4(scores)
        assert result == ["BP4_Supporting"]

    def test_spliceai_indeterminate(self):
        """SpliceAI in indeterminate range (0.1-0.2) — no evidence from splice."""
        scores = InSilicoScores(
            spliceai_ds_ag=0.15,
            spliceai_max=0.15,
        )
        result = generate_pp3_bp4(scores)
        assert result == []

    # --- Fallback: CADD + AlphaMissense ---

    def test_pp3_fallback_cadd_am(self):
        """No REVEL, CADD >= 25 + AM pathogenic → PP3_Supporting."""
        scores = InSilicoScores(
            cadd_phred=32.5,
            alphamissense_class="likely_pathogenic",
            alphamissense_score=0.92,
        )
        result = generate_pp3_bp4(scores)
        assert result == ["PP3_Supporting"]

    def test_fallback_cadd_only_no_am(self):
        """No REVEL, high CADD but no AM class → no evidence."""
        scores = InSilicoScores(cadd_phred=32.5)
        result = generate_pp3_bp4(scores)
        assert result == []

    def test_fallback_am_only_no_cadd(self):
        """No REVEL, AM pathogenic but no CADD → no evidence."""
        scores = InSilicoScores(
            alphamissense_class="likely_pathogenic",
            alphamissense_score=0.92,
        )
        result = generate_pp3_bp4(scores)
        assert result == []

    def test_fallback_cadd_low_am_pathogenic(self):
        """No REVEL, CADD < 25 + AM pathogenic → no PP3 (CADD not high enough)."""
        scores = InSilicoScores(
            cadd_phred=20.0,
            alphamissense_class="likely_pathogenic",
        )
        result = generate_pp3_bp4(scores)
        assert result == []

    def test_fallback_bp4_benign(self):
        """No REVEL, low CADD + AM benign → BP4_Supporting."""
        scores = InSilicoScores(
            cadd_phred=8.0,
            alphamissense_class="likely_benign",
            alphamissense_score=0.05,
        )
        result = generate_pp3_bp4(scores)
        assert result == ["BP4_Supporting"]

    # --- Missing / edge cases ---

    def test_missing_scores(self):
        """All None scores → empty evidence list."""
        scores = InSilicoScores()
        result = generate_pp3_bp4(scores)
        assert result == []

    def test_no_evidence_indeterminate_zone(self):
        """Scores in indeterminate zone → no evidence."""
        scores = InSilicoScores(
            revel=0.40,
            spliceai_max=0.15,
        )
        result = generate_pp3_bp4(scores)
        assert result == []

    def test_mutual_exclusivity(self):
        """PP3 and BP4 cannot coexist — pathogenic wins when in conflict."""
        # REVEL moderate pathogenic + SpliceAI benign → PP3 wins
        scores = InSilicoScores(
            revel=0.70,
            spliceai_ds_ag=0.02,
            spliceai_ds_al=0.01,
            spliceai_ds_dg=0.03,
            spliceai_ds_dl=0.01,
            spliceai_max=0.03,
        )
        result = generate_pp3_bp4(scores)
        # Should return PP3, not BP4
        assert len(result) == 1
        assert result[0].startswith("PP3")

    def test_mutual_exclusivity_splice_pp3_revel_bp4(self):
        """SpliceAI strong PP3 + REVEL benign → PP3 wins (pathogenic takes precedence)."""
        scores = InSilicoScores(
            revel=0.01,
            spliceai_ds_al=0.80,
            spliceai_max=0.80,
        )
        result = generate_pp3_bp4(scores)
        assert len(result) == 1
        assert result[0] == "PP3_Strong"

    def test_both_pp3_picks_strongest(self):
        """REVEL moderate + SpliceAI strong → strongest PP3 wins."""
        scores = InSilicoScores(
            revel=0.70,
            spliceai_ds_dg=0.60,
            spliceai_max=0.60,
        )
        result = generate_pp3_bp4(scores)
        assert result == ["PP3_Strong"]

    def test_both_bp4_picks_strongest(self):
        """REVEL strong BP4 + SpliceAI supporting BP4 → strongest BP4 wins."""
        scores = InSilicoScores(
            revel=0.01,
            spliceai_ds_ag=0.02,
            spliceai_max=0.02,
        )
        result = generate_pp3_bp4(scores)
        assert result == ["BP4_Strong"]

    # --- Custom thresholds ---

    def test_custom_thresholds(self):
        """Override REVEL thresholds via thresholds dict."""
        scores = InSilicoScores(revel=0.80)

        # With default thresholds, 0.80 → PP3_Moderate
        default_result = generate_pp3_bp4(scores)
        assert default_result == ["PP3_Moderate"]

        # With custom threshold lowering PP3_Strong to 0.75
        custom = {"revel_pp3_strong": 0.75}
        custom_result = generate_pp3_bp4(scores, thresholds=custom)
        assert custom_result == ["PP3_Strong"]

    def test_custom_thresholds_partial_override(self):
        """Partial threshold override — non-overridden keys keep defaults."""
        scores = InSilicoScores(revel=0.01)

        custom = {"revel_pp3_strong": 0.99}
        result = generate_pp3_bp4(scores, thresholds=custom)
        # BP4 thresholds unchanged, REVEL 0.01 ≤ 0.016 → BP4_Strong
        assert result == ["BP4_Strong"]

    def test_custom_spliceai_thresholds(self):
        """Override SpliceAI thresholds."""
        scores = InSilicoScores(spliceai_ds_ag=0.35, spliceai_max=0.35)

        # Default: 0.35 is between moderate (0.2) and strong (0.5) → PP3_Moderate
        default_result = generate_pp3_bp4(scores)
        assert default_result == ["PP3_Moderate"]

        # Custom: lower strong to 0.3
        custom = {"spliceai_pp3_strong": 0.3}
        custom_result = generate_pp3_bp4(scores, thresholds=custom)
        assert custom_result == ["PP3_Strong"]


# ═══════════════════════════════════════════════════════════════════════════
# Display formatting
# ═══════════════════════════════════════════════════════════════════════════


class TestFormatScores:
    """Tests for format_scores_for_display."""

    def test_format_scores_full(self):
        """All scores present → formatted dict with all keys."""
        scores = InSilicoScores(
            revel=0.89,
            cadd_phred=32.0,
            alphamissense_class="likely_pathogenic",
            alphamissense_score=0.95,
            spliceai_max=0.85,
            spliceai_ds_al=0.85,
            spliceai_ds_ag=0.10,
            spliceai_ds_dg=0.02,
            spliceai_ds_dl=0.01,
            sift="deleterious(0.01)",
            polyphen="probably_damaging(0.998)",
        )
        result = format_scores_for_display(scores)

        assert result["REVEL"] == "0.890"
        assert result["CADD"] == "32.0"
        assert "likely_pathogenic" in result["AlphaMissense"]
        assert "(0.950)" in result["AlphaMissense"]
        assert "0.85" in result["SpliceAI"]
        assert "AL" in result["SpliceAI"]
        assert result["SIFT"] == "deleterious(0.01)"
        assert result["PolyPhen"] == "probably_damaging(0.998)"

    def test_format_scores_partial(self):
        """Only REVEL present → only REVEL in output."""
        scores = InSilicoScores(revel=0.65)
        result = format_scores_for_display(scores)

        assert "REVEL" in result
        assert result["REVEL"] == "0.650"
        assert "CADD" not in result
        assert "AlphaMissense" not in result
        assert "SpliceAI" not in result
        assert "SIFT" not in result
        assert "PolyPhen" not in result

    def test_format_scores_empty(self):
        """No scores → empty dict."""
        scores = InSilicoScores()
        result = format_scores_for_display(scores)
        assert result == {}

    def test_format_scores_am_class_only(self):
        """AlphaMissense class but no score → shows class without score."""
        scores = InSilicoScores(alphamissense_class="ambiguous")
        result = format_scores_for_display(scores)
        assert result["AlphaMissense"] == "ambiguous"

    def test_format_scores_am_score_only(self):
        """AlphaMissense score but no class → shows score only."""
        scores = InSilicoScores(alphamissense_score=0.45)
        result = format_scores_for_display(scores)
        assert result["AlphaMissense"] == "(0.450)"

    def test_format_spliceai_max_label_dg(self):
        """SpliceAI label reflects the delta score that is the max (DG)."""
        scores = InSilicoScores(
            spliceai_max=0.60,
            spliceai_ds_ag=0.10,
            spliceai_ds_al=0.20,
            spliceai_ds_dg=0.60,
            spliceai_ds_dl=0.05,
        )
        result = format_scores_for_display(scores)
        assert "DG" in result["SpliceAI"]


# ═══════════════════════════════════════════════════════════════════════════
# Integration: parse from VCF CSQ fields
# ═══════════════════════════════════════════════════════════════════════════


class TestVCFIntegration:
    """Integration tests parsing the demo VCF with in silico fields."""

    @pytest.fixture()
    def vcf_records(self):
        """Parse the demo VCF and return a list of (csq_field_dict, gene) tuples."""
        from scripts.intake.parse_annotation import parse_csq_header

        records = []
        fields = []
        with open(DEMO_VCF) as f:
            for line in f:
                line = line.strip()
                if line.startswith("##INFO=<ID=CSQ"):
                    fields = parse_csq_header(line)
                elif line.startswith("#"):
                    continue
                elif line and fields:
                    cols = line.split("\t")
                    info = cols[7]
                    csq_raw = ""
                    gene = ""
                    for part in info.split(";"):
                        if part.startswith("CSQ="):
                            csq_raw = part[4:]
                        elif part.startswith("Gene="):
                            gene = part[5:]
                    # Parse the first (best) transcript entry
                    if csq_raw:
                        values = csq_raw.split(",")[0].split("|")
                        entry = {}
                        for i, val in enumerate(values):
                            if i < len(fields):
                                entry[fields[i].lower()] = val
                        records.append((entry, gene))
        return records

    def test_tp53_high_revel_pp3_strong(self, vcf_records):
        """TP53 variant: REVEL=0.95 → PP3_Strong."""
        csq, gene = next((r for r in vcf_records if r[1] == "TP53"), (None, None))
        assert csq is not None
        scores = parse_in_silico_from_csq(csq)
        assert scores.revel == pytest.approx(0.95)
        evidence = generate_pp3_bp4(scores)
        assert evidence == ["PP3_Strong"]

    def test_brca2_moderate_revel_pp3_moderate(self, vcf_records):
        """BRCA2 variant: REVEL=0.70 → PP3_Moderate."""
        csq, gene = next((r for r in vcf_records if r[1] == "BRCA2"), (None, None))
        scores = parse_in_silico_from_csq(csq)
        assert scores.revel == pytest.approx(0.70)
        evidence = generate_pp3_bp4(scores)
        assert evidence == ["PP3_Moderate"]

    def test_apoe_low_revel_bp4(self, vcf_records):
        """APOE variant: REVEL=0.05 → BP4_Moderate."""
        csq, gene = next((r for r in vcf_records if r[1] == "APOE"), (None, None))
        scores = parse_in_silico_from_csq(csq)
        assert scores.revel == pytest.approx(0.05)
        evidence = generate_pp3_bp4(scores)
        assert evidence == ["BP4_Moderate"]

    def test_cyp2c19_splice_pp3_strong(self, vcf_records):
        """CYP2C19 splice variant: SpliceAI max=0.85 → PP3_Strong."""
        csq, gene = next((r for r in vcf_records if r[1] == "CYP2C19"), (None, None))
        scores = parse_in_silico_from_csq(csq)
        assert scores.spliceai_max == pytest.approx(0.85)
        evidence = generate_pp3_bp4(scores)
        assert evidence == ["PP3_Strong"]

    def test_mlh1_missing_scores(self, vcf_records):
        """MLH1 variant: all in silico fields empty → no evidence."""
        csq, gene = next((r for r in vcf_records if r[1] == "MLH1"), (None, None))
        scores = parse_in_silico_from_csq(csq)
        assert scores.revel is None
        assert scores.cadd_phred is None
        # Only SpliceAI present (all zeros)
        evidence = generate_pp3_bp4(scores)
        # SpliceAI zeros < 0.1 → BP4_Supporting
        assert evidence == ["BP4_Supporting"]

    def test_msh2_fallback_cadd_am(self, vcf_records):
        """MSH2 variant: no REVEL, CADD=32.5 + AM=likely_pathogenic → PP3_Supporting."""
        csq, gene = next((r for r in vcf_records if r[1] == "MSH2"), (None, None))
        scores = parse_in_silico_from_csq(csq)
        assert scores.revel is None
        assert scores.cadd_phred == pytest.approx(32.5)
        assert scores.alphamissense_class == "likely_pathogenic"
        evidence = generate_pp3_bp4(scores)
        assert evidence == ["PP3_Supporting"]

    def test_apc_low_revel_bp4(self, vcf_records):
        """APC variant: REVEL=0.01 → BP4_Strong."""
        csq, gene = next((r for r in vcf_records if r[1] == "APC"), (None, None))
        scores = parse_in_silico_from_csq(csq)
        assert scores.revel == pytest.approx(0.01)
        evidence = generate_pp3_bp4(scores)
        assert evidence == ["BP4_Strong"]

    def test_atm_splice_acceptor_high_spliceai(self, vcf_records):
        """ATM splice acceptor: SpliceAI DG=0.60 → PP3_Strong."""
        csq, gene = next((r for r in vcf_records if r[1] == "ATM"), (None, None))
        scores = parse_in_silico_from_csq(csq)
        assert scores.spliceai_max == pytest.approx(0.60)
        evidence = generate_pp3_bp4(scores)
        assert evidence == ["PP3_Strong"]
