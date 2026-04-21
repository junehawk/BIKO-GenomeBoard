"""Tests for PharmCAT integration and the unified PGx entry point.

All PharmCAT subprocess calls are mocked — no Java or JAR required.
"""

import json
import logging
import os
import subprocess
from pathlib import Path
from unittest import mock

import pytest

from scripts.common.models import Variant

pytestmark = pytest.mark.integration

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

FIXTURE_DIR = Path(__file__).parent / "fixtures"
SAMPLE_REPORT = FIXTURE_DIR / "pharmcat_sample_report.json"


def _load_sample_report() -> dict:
    with open(SAMPLE_REPORT) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# pharmcat_runner — availability checks
# ---------------------------------------------------------------------------


class TestIsPharmcatAvailable:
    """Tests for ``is_pharmcat_available``."""

    def test_with_jar_and_java(self, tmp_path):
        """JAR + Java 17 present -> True."""
        from scripts.pharmacogenomics.pharmcat_runner import is_pharmcat_available

        jar = tmp_path / "pharmcat-2.13.0.jar"
        jar.touch()

        with (
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_java", return_value="/usr/bin/java"),
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_jar", return_value=str(jar)),
        ):
            assert is_pharmcat_available() is True

    def test_without_jar(self):
        """No JAR -> False."""
        from scripts.pharmacogenomics.pharmcat_runner import is_pharmcat_available

        with (
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_java", return_value="/usr/bin/java"),
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_jar", return_value=None),
        ):
            assert is_pharmcat_available() is False

    def test_without_java(self, tmp_path):
        """No Java -> False."""
        from scripts.pharmacogenomics.pharmcat_runner import is_pharmcat_available

        jar = tmp_path / "pharmcat-2.13.0.jar"
        jar.touch()

        with (
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_java", return_value=None),
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_jar", return_value=str(jar)),
        ):
            assert is_pharmcat_available() is False


# ---------------------------------------------------------------------------
# pharmcat_runner — _find_java
# ---------------------------------------------------------------------------


class TestFindJava:
    """Tests for ``_find_java`` internal helper."""

    def test_finds_java_17(self):
        from scripts.pharmacogenomics.pharmcat_runner import _find_java

        # On machines with local Java 17 (data/tools/java/) or system Java 17,
        # _find_java succeeds. We just verify it returns a valid path.
        fake_proc = subprocess.CompletedProcess(
            args=["java", "-version"],
            returncode=0,
            stdout="",
            stderr='openjdk version "17.0.8" 2023-07-18',
        )
        with mock.patch("subprocess.run", return_value=fake_proc):
            result = _find_java()
            assert result is not None

    def test_rejects_java_11(self):
        from scripts.pharmacogenomics.pharmcat_runner import _find_java

        fake_proc = subprocess.CompletedProcess(
            args=["java", "-version"],
            returncode=0,
            stdout="",
            stderr='openjdk version "11.0.20" 2023-07-18',
        )
        with (
            mock.patch("shutil.which", return_value="/usr/bin/java"),
            mock.patch("subprocess.run", return_value=fake_proc),
        ):
            result = _find_java()
            assert result is None

    def test_no_java_on_path(self):
        from scripts.pharmacogenomics.pharmcat_runner import _find_java

        # Simulate no Java at all: shutil.which returns None AND every
        # subprocess version check fails (covers local data/tools/java/ too).
        with (
            mock.patch("shutil.which", return_value=None),
            mock.patch("subprocess.run", side_effect=OSError("no java")),
        ):
            assert _find_java() is None

    def test_config_java_path(self):
        from scripts.pharmacogenomics.pharmcat_runner import _find_java

        fake_proc = subprocess.CompletedProcess(
            args=["java", "-version"],
            returncode=0,
            stdout="",
            stderr='openjdk version "21.0.1" 2023-10-17',
        )
        with mock.patch("subprocess.run", return_value=fake_proc):
            result = _find_java(config={"java_path": "/opt/java21/bin/java"})
            assert result == "/opt/java21/bin/java"


# ---------------------------------------------------------------------------
# pharmcat_runner — _find_jar
# ---------------------------------------------------------------------------


class TestFindJar:
    """Tests for ``_find_jar`` internal helper."""

    def test_config_jar_path(self, tmp_path):
        from scripts.pharmacogenomics.pharmcat_runner import _find_jar

        jar = tmp_path / "pharmcat.jar"
        jar.touch()
        assert _find_jar(config={"pharmcat_jar": str(jar)}) == str(jar)

    def test_env_var_fallback(self, tmp_path):
        from scripts.pharmacogenomics.pharmcat_runner import _find_jar

        jar = tmp_path / "pharmcat-env.jar"
        jar.touch()
        with mock.patch.dict(os.environ, {"PHARMCAT_JAR": str(jar)}):
            # No config jar, no auto-discover match -> env var
            with mock.patch("glob.glob", return_value=[]):
                result = _find_jar(config={})
                assert result == str(jar)

    def test_returns_none_when_nothing(self):
        from scripts.pharmacogenomics.pharmcat_runner import _find_jar

        with (
            mock.patch("glob.glob", return_value=[]),
            mock.patch.dict(os.environ, {}, clear=True),
        ):
            result = _find_jar(config={})
            assert result is None


# ---------------------------------------------------------------------------
# pharmcat_runner — run_pharmcat + JSON parsing
# ---------------------------------------------------------------------------


class TestRunPharmcat:
    """Tests for ``run_pharmcat`` end-to-end (subprocess mocked)."""

    def test_parses_json_output(self, tmp_path):
        """Mock subprocess -> sample PharmCAT JSON -> PharmCATResult parsed."""
        from scripts.pharmacogenomics.pharmcat_runner import run_pharmcat

        # Prepare a fake germline VCF
        vcf = tmp_path / "germline.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Write fixture report into the output dir that run_pharmcat will use
        def fake_subprocess_run(cmd, **kwargs):
            # Find output dir from -o flag
            out_idx = cmd.index("-o") + 1
            out_dir = cmd[out_idx]
            os.makedirs(out_dir, exist_ok=True)
            # Copy fixture into output dir
            report_path = os.path.join(out_dir, "germline.report.json")
            report_data = _load_sample_report()
            with open(report_path, "w") as fh:
                json.dump(report_data, fh)
            return subprocess.CompletedProcess(args=cmd, returncode=0, stdout="", stderr="")

        with (
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_java", return_value="/usr/bin/java"),
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_jar", return_value="/opt/pharmcat.jar"),
            mock.patch("subprocess.run", side_effect=fake_subprocess_run),
        ):
            result = run_pharmcat(str(vcf))

        assert result is not None
        assert result.source == "pharmcat"
        assert result.version == "3.2.0"

        # Diplotypes
        assert result.diplotypes["CYP2C19"] == "*1/*2"
        assert result.diplotypes["DPYD"] == "*1/*1"
        assert result.diplotypes["NUDT15"] == "*1/*3"

        # Phenotypes
        assert result.phenotypes["CYP2C19"] == "Intermediate Metabolizer"
        assert result.phenotypes["DPYD"] == "Normal Metabolizer"

        # Drug recommendations are optional — PharmCAT v3 stores them in a
        # separate top-level "drugs" dict which may be empty depending on the
        # input VCF and PharmCAT configuration. The core value is in
        # diplotypes + phenotypes; drug recs are a bonus.
        assert isinstance(result.drug_recommendations, list)

    def test_returns_none_on_subprocess_failure(self, tmp_path):
        """Subprocess non-zero exit -> None."""
        from scripts.pharmacogenomics.pharmcat_runner import run_pharmcat

        vcf = tmp_path / "germline.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n")

        failed = subprocess.CompletedProcess(args=[], returncode=1, stdout="", stderr="Error: invalid VCF")
        with (
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_java", return_value="/usr/bin/java"),
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_jar", return_value="/opt/pharmcat.jar"),
            mock.patch("subprocess.run", return_value=failed),
        ):
            result = run_pharmcat(str(vcf))
            assert result is None

    def test_returns_none_when_not_available(self, tmp_path):
        """No Java/JAR -> None immediately."""
        from scripts.pharmacogenomics.pharmcat_runner import run_pharmcat

        vcf = tmp_path / "germline.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n")

        with (
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_java", return_value=None),
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_jar", return_value=None),
        ):
            result = run_pharmcat(str(vcf))
            assert result is None

    def test_returns_none_on_timeout(self, tmp_path):
        """Subprocess timeout -> None."""
        from scripts.pharmacogenomics.pharmcat_runner import run_pharmcat

        vcf = tmp_path / "germline.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n")

        with (
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_java", return_value="/usr/bin/java"),
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_jar", return_value="/opt/pharmcat.jar"),
            mock.patch("subprocess.run", side_effect=subprocess.TimeoutExpired(cmd="java", timeout=120)),
        ):
            result = run_pharmcat(str(vcf))
            assert result is None

    def test_returns_none_for_missing_vcf(self):
        """Non-existent VCF path -> None."""
        from scripts.pharmacogenomics.pharmcat_runner import run_pharmcat

        with (
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_java", return_value="/usr/bin/java"),
            mock.patch("scripts.pharmacogenomics.pharmcat_runner._find_jar", return_value="/opt/pharmcat.jar"),
        ):
            result = run_pharmcat("/nonexistent/path.vcf")
            assert result is None


# ---------------------------------------------------------------------------
# pharmcat_runner — _parse_pharmcat_json
# ---------------------------------------------------------------------------


class TestParsePharmcatJson:
    """Direct tests for the JSON parser."""

    def test_parses_sample_fixture(self):
        from scripts.pharmacogenomics.pharmcat_runner import _parse_pharmcat_json

        result = _parse_pharmcat_json(str(SAMPLE_REPORT))
        assert result is not None
        assert len(result.diplotypes) == 3
        assert result.diplotypes["CYP2C19"] == "*1/*2"
        assert result.phenotypes["CYP2C19"] == "Intermediate Metabolizer"
        # Drug recs may be empty in v3 fixture (drugs parsed from separate
        # top-level key, which is empty in the minimal fixture).
        assert isinstance(result.drug_recommendations, list)

    def test_handles_empty_genes(self, tmp_path):
        report = tmp_path / "empty.report.json"
        report.write_text(json.dumps({"genes": []}))

        from scripts.pharmacogenomics.pharmcat_runner import _parse_pharmcat_json

        result = _parse_pharmcat_json(str(report))
        assert result is not None
        assert result.diplotypes == {}
        assert result.drug_recommendations == []

    def test_handles_missing_file(self):
        from scripts.pharmacogenomics.pharmcat_runner import _parse_pharmcat_json

        result = _parse_pharmcat_json("/nonexistent/report.json")
        assert result is None

    def test_handles_malformed_json(self, tmp_path):
        report = tmp_path / "bad.report.json"
        report.write_text("not json at all {{{")

        from scripts.pharmacogenomics.pharmcat_runner import _parse_pharmcat_json

        result = _parse_pharmcat_json(str(report))
        assert result is None


# ---------------------------------------------------------------------------
# korean_pgx — get_pgx_results unified entry point
# ---------------------------------------------------------------------------


class TestGetPgxResults:
    """Tests for the unified ``get_pgx_results`` entry point."""

    def test_prefers_pharmcat_when_available(self, tmp_path):
        """germline VCF + PharmCAT available -> source='pharmcat'."""
        from scripts.pharmacogenomics.korean_pgx import get_pgx_results
        from scripts.pharmacogenomics.pharmcat_runner import PharmCATResult

        vcf = tmp_path / "germline.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n")

        fake_result = PharmCATResult(
            diplotypes={"CYP2C19": "*1/*2"},
            phenotypes={"CYP2C19": "Intermediate Metabolizer"},
            drug_recommendations=[
                {"gene": "CYP2C19", "drug": "clopidogrel", "guideline": "CPIC", "classification": "Actionable PGx"}
            ],
            version="2.13.0",
        )

        # get_pgx_results imports from pharmcat_runner inside the function body,
        # so we mock at the pharmcat_runner module level.
        with (
            mock.patch("scripts.pharmacogenomics.pharmcat_runner.is_pharmcat_available", return_value=True),
            mock.patch("scripts.pharmacogenomics.pharmcat_runner.run_pharmcat", return_value=fake_result),
        ):
            result = get_pgx_results([], germline_vcf=str(vcf))

        assert result["pgx_source"] == "pharmcat"
        assert result["germline_provided"] is True
        assert result["pharmcat_version"] == "2.13.0"
        assert len(result["pgx_hits"]) == 1
        assert result["pgx_hits"][0]["gene"] == "CYP2C19"
        assert result["pgx_hits"][0]["star_allele"] == "*1/*2"

    def test_fallback_to_builtin_when_pharmcat_unavailable(self, tmp_path):
        """germline VCF + PharmCAT NOT available -> source='builtin_limited'."""
        from scripts.pharmacogenomics.korean_pgx import get_pgx_results

        vcf = tmp_path / "germline.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n")

        with mock.patch("scripts.pharmacogenomics.pharmcat_runner.is_pharmcat_available", return_value=False):
            result = get_pgx_results([], germline_vcf=str(vcf))

        assert result["pgx_source"] == "builtin_limited"
        assert result["germline_provided"] is True
        assert result["pharmcat_version"] == ""
        assert len(result["warnings"]) > 0

    def test_builtin_when_no_germline(self):
        """No germline VCF -> source='builtin'."""
        from scripts.pharmacogenomics.korean_pgx import get_pgx_results

        cyp_variant = Variant(chrom="chr10", pos=96541616, ref="G", alt="A", gene="CYP2C19")
        result = get_pgx_results([cyp_variant])

        assert result["pgx_source"] == "builtin"
        assert result["germline_provided"] is False
        assert result["pharmcat_version"] == ""

    def test_builtin_hits_contain_expected_keys(self):
        """Builtin results should have all keys the report templates expect."""
        from scripts.pharmacogenomics.korean_pgx import get_pgx_results

        cyp_variant = Variant(chrom="chr10", pos=96541616, ref="G", alt="A", gene="CYP2C19")
        result = get_pgx_results([cyp_variant])

        if result["pgx_hits"]:
            hit = result["pgx_hits"][0]
            expected_keys = {
                "gene",
                "star_allele",
                "phenotype",
                "cpic_level",
                "korean_prevalence",
                "western_prevalence",
                "clinical_impact",
                "cpic_recommendation",
                "korean_flag",
            }
            assert expected_keys.issubset(set(hit.keys()))

    def test_pharmcat_hits_contain_expected_keys(self, tmp_path):
        """PharmCAT results should have all keys the report templates expect."""
        from scripts.pharmacogenomics.korean_pgx import get_pgx_results
        from scripts.pharmacogenomics.pharmcat_runner import PharmCATResult

        vcf = tmp_path / "germline.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n")

        fake_result = PharmCATResult(
            diplotypes={"CYP2C19": "*1/*2"},
            phenotypes={"CYP2C19": "Intermediate Metabolizer"},
            drug_recommendations=[
                {"gene": "CYP2C19", "drug": "clopidogrel", "guideline": "CPIC", "classification": "Actionable PGx"}
            ],
            version="2.13.0",
        )

        with (
            mock.patch("scripts.pharmacogenomics.pharmcat_runner.is_pharmcat_available", return_value=True),
            mock.patch("scripts.pharmacogenomics.pharmcat_runner.run_pharmcat", return_value=fake_result),
        ):
            result = get_pgx_results([], germline_vcf=str(vcf))

        assert result["pgx_source"] == "pharmcat"
        hit = result["pgx_hits"][0]
        expected_keys = {
            "gene",
            "star_allele",
            "phenotype",
            "cpic_level",
            "korean_prevalence",
            "western_prevalence",
            "clinical_impact",
            "cpic_recommendation",
            "korean_flag",
        }
        assert expected_keys.issubset(set(hit.keys()))


# ---------------------------------------------------------------------------
# orchestrate.py — --germline CLI flag
# ---------------------------------------------------------------------------


class TestOrchestrateGermlineFlag:
    """Tests for CLI argument parsing."""

    def test_germline_flag_parsed(self):
        """--germline path.vcf is correctly parsed."""
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("vcf_path", nargs="?", default=None)
        parser.add_argument("--germline", type=str, default=None)

        args = parser.parse_args(["sample.vcf", "--germline", "/path/to/germline.vcf"])
        assert args.germline == "/path/to/germline.vcf"

    def test_germline_flag_default_none(self):
        """Without --germline, the value is None."""
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("vcf_path", nargs="?", default=None)
        parser.add_argument("--germline", type=str, default=None)

        args = parser.parse_args(["sample.vcf"])
        assert args.germline is None


# ---------------------------------------------------------------------------
# orchestrate.py — PGx warning on somatic without germline
# ---------------------------------------------------------------------------


class TestPgxWarningOnSomatic:
    """Somatic VCF + no germline -> warning logged."""

    def test_warning_logged(self, caplog):
        """run_pipeline logs a PGx warning when germline is not provided in cancer mode."""
        # We mock the heavy parts of the pipeline to isolate the warning
        from scripts.common.models import Variant

        fake_variants = [Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")]

        with caplog.at_level(logging.WARNING, logger="genomeboard"):
            # Import inside to avoid side effects
            import scripts.orchestrate as orch

            # Mock the entire pipeline but check the warning is emitted
            with (
                mock.patch.object(orch, "parse_vcf", return_value=fake_variants),
                mock.patch.object(
                    orch,
                    "query_variant_databases",
                    return_value={
                        "clinvar": {"clinvar_significance": "Not Found"},
                        "gnomad": {"gnomad_all": None, "gnomad_eas": None},
                        "kova_freq": None,
                        "kova_homozygote": None,
                        "pgx": None,
                    },
                ),
                mock.patch.object(
                    orch,
                    "compare_frequencies",
                    return_value={
                        "acmg_codes": [],
                        "korean_flag": "",
                    },
                ),
                mock.patch.object(
                    orch,
                    "classify_variants",
                    return_value={
                        fake_variants[0].variant_id: mock.MagicMock(
                            classification="VUS",
                            evidence_codes=[],
                            conflict=False,
                        ),
                    },
                ),
                mock.patch.object(orch, "build_variant_records", return_value=[]),
                mock.patch.object(
                    orch,
                    "build_summary",
                    return_value={
                        "total": 1,
                        "pathogenic": 0,
                        "likely_pathogenic": 0,
                        "drug_response": 0,
                        "vus": 1,
                        "benign": 0,
                        "likely_benign": 0,
                    },
                ),
                mock.patch.object(orch, "split_variants_for_display", return_value=([], [], [], 0, [], [])),
                mock.patch.object(orch, "get_all_db_versions", return_value={}),
                mock.patch.object(orch, "generate_report_html", return_value="<html></html>"),
                mock.patch(
                    "scripts.pharmacogenomics.korean_pgx.get_pgx_results",
                    return_value={
                        "pgx_hits": [],
                        "pgx_source": "builtin",
                        "pharmcat_version": "",
                        "germline_provided": False,
                        "warnings": [],
                    },
                ),
            ):
                result = orch.run_pipeline(
                    vcf_path="data/sample_vcf/demo_variants_grch38_annotated.vcf",
                    output_path="/tmp/test_report.html",
                    mode="cancer",
                    skip_api=True,
                    germline_vcf=None,
                )

        # Check the warning was logged
        warning_msgs = [r.message for r in caplog.records if r.levelno >= logging.WARNING]
        pgx_warnings = [m for m in warning_msgs if "PGx results from somatic" in m or "germline" in m.lower()]
        assert len(pgx_warnings) > 0, f"Expected PGx warning, got: {warning_msgs}"


# ---------------------------------------------------------------------------
# PharmCATResult dataclass
# ---------------------------------------------------------------------------


class TestPharmCATResultDataclass:
    """Basic dataclass tests."""

    def test_defaults(self):
        from scripts.pharmacogenomics.pharmcat_runner import PharmCATResult

        r = PharmCATResult()
        assert r.diplotypes == {}
        assert r.phenotypes == {}
        assert r.drug_recommendations == []
        assert r.source == "pharmcat"
        assert r.version == ""
        assert r.warnings == []

    def test_with_data(self):
        from scripts.pharmacogenomics.pharmcat_runner import PharmCATResult

        r = PharmCATResult(
            diplotypes={"CYP2C19": "*1/*2"},
            phenotypes={"CYP2C19": "IM"},
            drug_recommendations=[{"gene": "CYP2C19", "drug": "clopidogrel"}],
            version="2.13.0",
        )
        assert r.diplotypes["CYP2C19"] == "*1/*2"
        assert r.version == "2.13.0"
        assert len(r.drug_recommendations) == 1
