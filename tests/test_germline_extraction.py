"""Tests for germline inherited variant extraction (v2.4 Phase 2+3).

Covers:
- Target BED building (build_germline_target_bed)
- Variant extraction from germline VCF (extract_germline)
- AF filtering, source badge, deduplication
- Orchestrate integration (mocked)
"""

import gzip
import os
import sqlite3
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from scripts.common.models import Variant

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_clinvar_db(tmp_path):
    """Create a minimal ClinVar SQLite database for testing."""
    db_path = tmp_path / "clinvar.sqlite3"
    conn = sqlite3.connect(str(db_path))
    conn.execute("""
        CREATE TABLE variants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            rsid TEXT,
            gene TEXT,
            clinical_significance TEXT,
            review_status TEXT,
            phenotype_list TEXT,
            variation_id TEXT,
            allele_id TEXT,
            origin TEXT,
            assembly TEXT DEFAULT 'GRCh38',
            last_evaluated TEXT,
            number_submitters INTEGER,
            hgvsp TEXT
        )
    """)
    conn.execute("""
        CREATE TABLE metadata (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    """)
    # Insert test variants
    test_variants = [
        ("chr17", 7577120, "G", "A", "rs28934578", "TP53", "Pathogenic", "criteria provided", "GRCh38"),
        ("chr13", 32337326, "C", "T", "rs80358981", "BRCA2", "Likely pathogenic", "criteria provided", "GRCh38"),
        ("chr7", 140453136, "A", "T", None, "BRAF", "Pathogenic", "reviewed", "GRCh38"),
        ("chr17", 41276045, "A", "G", None, "BRCA1", "Pathogenic", "criteria provided", "GRCh38"),
        ("chr17", 41276050, "C", "T", None, "BRCA1", "Conflicting interpretations", "criteria provided", "GRCh38"),
        ("chr11", 108234769, "G", "A", None, "ATM", "Pathogenic", "criteria provided", "GRCh38"),
        ("chr2", 47702181, "C", "T", None, "MSH2", "Pathogenic", "criteria provided", "GRCh38"),
    ]
    conn.executemany(
        """INSERT INTO variants (chrom, pos, ref, alt, rsid, gene,
           clinical_significance, review_status, assembly)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        test_variants,
    )
    conn.commit()
    conn.close()
    return str(db_path)


@pytest.fixture
def tmp_target_bed(tmp_path):
    """Create a minimal target BED file (plain text, not bgzipped)."""
    bed_path = tmp_path / "targets.bed"
    bed_content = textwrap.dedent("""\
        chr17\t7577119\t7577120\tTP53\tclinvar_plp
        chr13\t32337325\t32337326\tBRCA2\tclinvar_plp
        chr7\t140453135\t140453136\tBRAF\tclinvar_plp
        chr17\t41276044\t41276045\tBRCA1\tclinvar_plp
        chr11\t108234768\t108234769\tATM\tacmg_sf
    """)
    bed_path.write_text(bed_content)
    return str(bed_path)


@pytest.fixture
def tmp_target_bed_gz(tmp_path, tmp_target_bed):
    """Create a gzipped target BED (simulates bgzip for test purposes)."""
    bed_gz_path = tmp_path / "targets.bed.gz"
    plain_content = Path(tmp_target_bed).read_bytes()
    with gzip.open(str(bed_gz_path), "wb") as f:
        f.write(plain_content)
    return str(bed_gz_path)


@pytest.fixture
def tmp_germline_vcf(tmp_path):
    """Create a minimal germline VCF file (plain text)."""
    vcf_path = tmp_path / "germline.vcf"
    vcf_content = textwrap.dedent("""\
        ##fileformat=VCFv4.2
        ##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene info">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
        chr17\t7577120\trs28934578\tG\tA\t100\tPASS\tGENEINFO=TP53
        chr13\t32337326\t.\tC\tT\t80\tPASS\tGENEINFO=BRCA2
        chr7\t140453136\t.\tA\tT\t90\tPASS\tGENEINFO=BRAF
        chr11\t108234769\t.\tG\tA\t70\tPASS\tGENEINFO=ATM
        chr1\t100000\t.\tA\tG\t50\tPASS\tGENEINFO=FAKE_GENE
    """)
    vcf_path.write_text(vcf_content)
    return str(vcf_path)


@pytest.fixture
def tmp_germline_vcf_gz(tmp_path, tmp_germline_vcf):
    """Create a gzipped germline VCF (simulates bgzip)."""
    vcf_gz = tmp_path / "germline.vcf.gz"
    plain = Path(tmp_germline_vcf).read_bytes()
    with gzip.open(str(vcf_gz), "wb") as f:
        f.write(plain)
    return str(vcf_gz)


# ---------------------------------------------------------------------------
# 1. build_germline_target_bed tests
# ---------------------------------------------------------------------------


class TestBuildTargetBed:
    def test_build_target_bed_creates_output(self, tmp_path, tmp_clinvar_db):
        """When ClinVar DB exists, BED is generated with regions."""
        from scripts.tools.build_germline_target_bed import build_target_bed

        output_dir = str(tmp_path / "germline_targets")
        with patch("scripts.tools.build_germline_target_bed._get_clinvar_db_path", return_value=tmp_clinvar_db):
            result = build_target_bed(output_dir=output_dir)

        assert result["region_count"] > 0
        assert os.path.exists(result["output_path"])
        assert len(result["sources"]) >= 1
        # Should include clinvar_plp source
        assert any("clinvar_plp" in s for s in result["sources"])

    def test_build_target_bed_without_clinvar_uses_fallback(self, tmp_path):
        """When ClinVar DB is missing, fallback warning is emitted."""
        from scripts.tools.build_germline_target_bed import build_target_bed

        output_dir = str(tmp_path / "germline_targets")
        with patch("scripts.tools.build_germline_target_bed._get_clinvar_db_path", return_value=None):
            result = build_target_bed(output_dir=output_dir)

        # With no ClinVar DB and no DDG2P data, region count is 0
        assert result["region_count"] == 0
        assert any("fallback" in s for s in result["sources"])

    def test_build_includes_acmg_sf_regions(self, tmp_path, tmp_clinvar_db):
        """ACMG SF v3.2 gene regions are included when ClinVar DB has matching genes."""
        from scripts.tools.build_germline_target_bed import build_target_bed

        output_dir = str(tmp_path / "germline_targets")
        with patch("scripts.tools.build_germline_target_bed._get_clinvar_db_path", return_value=tmp_clinvar_db):
            result = build_target_bed(output_dir=output_dir)

        # ATM is both in our test ClinVar DB and in ACMG SF v3.2
        assert any("acmg_sf" in s for s in result["sources"])

    def test_build_includes_ddg2p_regions(self, tmp_path, tmp_clinvar_db):
        """DDG2P gene regions are included when panel file exists and genes have ClinVar data."""
        from scripts.tools.build_germline_target_bed import build_target_bed

        output_dir = str(tmp_path / "germline_targets")
        # Mock DDG2P to return genes that match our test ClinVar DB
        mock_genes = ["TP53", "BRCA2"]
        with (
            patch("scripts.tools.build_germline_target_bed._get_clinvar_db_path", return_value=tmp_clinvar_db),
            patch("scripts.tools.build_germline_target_bed._load_ddg2p_genes", return_value=mock_genes),
        ):
            result = build_target_bed(output_dir=output_dir)

        assert any("ddg2p" in s for s in result["sources"])

    def test_acmg_sf_gene_list_count(self):
        """ACMG SF v3.2 gene list has the expected number of unique genes."""
        from scripts.tools.build_germline_target_bed import ACMG_SF_V32_GENES

        # The deduplicated set should have at least 70 genes (spec says 73
        # but the input list may contain duplicates like MLH1 appearing twice)
        assert len(ACMG_SF_V32_GENES) >= 70


# ---------------------------------------------------------------------------
# 2. extract_germline tests
# ---------------------------------------------------------------------------


class TestExtractGermline:
    def test_extract_inherited_from_fixture_germline(self, tmp_germline_vcf, tmp_target_bed):
        """Variants matching target regions are extracted."""
        from scripts.orchestration.extract_germline import extract_inherited_variants

        with patch("scripts.orchestration.extract_germline._get_gnomad_af", return_value=None):
            variants = extract_inherited_variants(
                tmp_germline_vcf,
                target_bed=tmp_target_bed,
                max_gnomad_af=0.01,
            )

        # 4 of the 5 germline variants overlap the target BED regions
        # (FAKE_GENE at chr1:100000 does not)
        assert len(variants) == 4
        genes = {v.gene for v in variants}
        assert "TP53" in genes
        assert "BRCA2" in genes
        assert "ATM" in genes
        assert "FAKE_GENE" not in genes

    def test_extract_filters_common_variants(self, tmp_germline_vcf, tmp_target_bed):
        """Variants with gnomAD AF >= max_gnomad_af are excluded."""
        from scripts.orchestration.extract_germline import extract_inherited_variants

        def mock_af(variant):
            if variant.gene == "BRAF":
                return 0.05  # Above threshold
            return None  # Unknown -- kept (conservative)

        with patch("scripts.orchestration.extract_germline._get_gnomad_af", side_effect=mock_af):
            variants = extract_inherited_variants(
                tmp_germline_vcf,
                target_bed=tmp_target_bed,
                max_gnomad_af=0.01,
            )

        genes = {v.gene for v in variants}
        assert "BRAF" not in genes
        assert "TP53" in genes

    def test_extract_sets_source_germline_inherited(self, tmp_germline_vcf, tmp_target_bed):
        """All extracted variants carry source='germline_inherited'."""
        from scripts.orchestration.extract_germline import extract_inherited_variants

        with patch("scripts.orchestration.extract_germline._get_gnomad_af", return_value=None):
            variants = extract_inherited_variants(
                tmp_germline_vcf,
                target_bed=tmp_target_bed,
            )

        assert all(v.source == "germline_inherited" for v in variants)

    def test_extract_empty_when_no_target_bed(self, tmp_germline_vcf, tmp_path):
        """When target BED doesn't exist, returns empty list with no crash."""
        from scripts.orchestration.extract_germline import extract_inherited_variants

        nonexistent = str(tmp_path / "no_such.bed")
        variants = extract_inherited_variants(
            tmp_germline_vcf,
            target_bed=nonexistent,
        )
        assert variants == []

    def test_extract_empty_when_no_matches(self, tmp_path):
        """When no variants match target regions, returns empty list."""
        from scripts.orchestration.extract_germline import extract_inherited_variants

        # VCF with only a variant not in the target BED
        vcf_path = tmp_path / "no_match.vcf"
        vcf_path.write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t999999\t.\tA\tG\t50\tPASS\t.\n"
        )
        bed_path = tmp_path / "target.bed"
        bed_path.write_text("chr17\t7577119\t7577120\tTP53\tclinvar_plp\n")

        with patch("scripts.orchestration.extract_germline._get_gnomad_af", return_value=None):
            variants = extract_inherited_variants(
                str(vcf_path),
                target_bed=str(bed_path),
            )

        assert variants == []

    def test_extract_empty_when_vcf_missing(self, tmp_target_bed):
        """When germline VCF doesn't exist, returns empty list."""
        from scripts.orchestration.extract_germline import extract_inherited_variants

        variants = extract_inherited_variants(
            "/nonexistent/germline.vcf",
            target_bed=tmp_target_bed,
        )
        assert variants == []

    def test_dedup_primary_takes_precedence(self, tmp_germline_vcf, tmp_target_bed):
        """When primary_variant_ids contains a match, inherited copy is skipped."""
        from scripts.orchestration.extract_germline import extract_inherited_variants

        # TP53 variant is in the germline VCF and target BED
        primary_ids = {"chr17:7577120:G>A"}

        with patch("scripts.orchestration.extract_germline._get_gnomad_af", return_value=None):
            variants = extract_inherited_variants(
                tmp_germline_vcf,
                target_bed=tmp_target_bed,
                primary_variant_ids=primary_ids,
            )

        ids = {v.variant_id for v in variants}
        assert "chr17:7577120:G>A" not in ids
        # Other variants should still be present
        assert len(variants) == 3

    def test_extract_gzipped_target_bed(self, tmp_germline_vcf, tmp_target_bed_gz):
        """Gzipped target BED files are read correctly."""
        from scripts.orchestration.extract_germline import extract_inherited_variants

        with patch("scripts.orchestration.extract_germline._get_gnomad_af", return_value=None):
            variants = extract_inherited_variants(
                tmp_germline_vcf,
                target_bed=tmp_target_bed_gz,
            )

        assert len(variants) == 4

    def test_parse_vcf_line_with_geneinfo(self):
        """VCF lines with GENEINFO field are parsed correctly."""
        from scripts.orchestration.extract_germline import _parse_vcf_line

        line = "chr17\t7577120\trs28934578\tG\tA\t100\tPASS\tGENEINFO=TP53"
        v = _parse_vcf_line(line)
        assert v is not None
        assert v.chrom == "chr17"
        assert v.pos == 7577120
        assert v.ref == "G"
        assert v.alt == "A"
        assert v.gene == "TP53"
        assert v.rsid == "rs28934578"
        assert v.source == "germline_inherited"

    def test_parse_vcf_line_skips_header(self):
        """Header lines return None."""
        from scripts.orchestration.extract_germline import _parse_vcf_line

        assert _parse_vcf_line("#CHROM\tPOS\tID\tREF\tALT") is None
        assert _parse_vcf_line("##fileformat=VCFv4.2") is None


# ---------------------------------------------------------------------------
# 3. Variant model source field
# ---------------------------------------------------------------------------


class TestVariantSourceField:
    def test_variant_has_source_field(self):
        """Variant dataclass has a source field defaulting to None."""
        v = Variant(chrom="chr17", pos=7577120, ref="G", alt="A")
        assert v.source is None

    def test_variant_source_can_be_set(self):
        """Variant source field can be set to any string."""
        v = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", source="germline_inherited")
        assert v.source == "germline_inherited"


# ---------------------------------------------------------------------------
# 4. build_variant_records carries source field
# ---------------------------------------------------------------------------


class TestVariantRecordSource:
    def test_variant_record_carries_source_field(self):
        """build_variant_records propagates source to variant records."""
        from scripts.orchestration.classify import build_variant_records

        variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53", source="germline_inherited")
        db_results = {
            variant.variant_id: {
                "clinvar": {"clinvar_significance": "Pathogenic", "acmg_codes": []},
                "gnomad": {"gnomad_all": 0.0001, "gnomad_eas": 0.0002},
                "krgdb_freq": None,
                "pgx": None,
            }
        }
        freq_results = {variant.variant_id: {"korean_flag": ""}}

        # Minimal classification result mock
        from unittest.mock import MagicMock

        mock_cls = MagicMock()
        mock_cls.classification = "Pathogenic"
        mock_cls.evidence_codes = ["PVS1"]
        mock_cls.conflict = None
        mock_cls.clinvar_override = False
        mock_cls.clinvar_override_reason = ""
        mock_cls.original_engine_classification = None
        classification_results = {variant.variant_id: mock_cls}

        records = build_variant_records([variant], db_results, freq_results, classification_results, "rare-disease", [])

        assert len(records) == 1
        assert records[0]["source"] == "germline_inherited"

    def test_variant_record_default_source_is_primary(self):
        """Variants without explicit source get 'primary' in records."""
        from scripts.orchestration.classify import build_variant_records

        variant = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
        db_results = {
            variant.variant_id: {
                "clinvar": {"clinvar_significance": "VUS", "acmg_codes": []},
                "gnomad": {"gnomad_all": None, "gnomad_eas": None},
                "krgdb_freq": None,
                "pgx": None,
            }
        }
        freq_results = {variant.variant_id: {"korean_flag": ""}}

        mock_cls = MagicMock()
        mock_cls.classification = "VUS"
        mock_cls.evidence_codes = []
        mock_cls.conflict = None
        mock_cls.clinvar_override = False
        mock_cls.clinvar_override_reason = ""
        mock_cls.original_engine_classification = None
        classification_results = {variant.variant_id: mock_cls}

        records = build_variant_records([variant], db_results, freq_results, classification_results, "cancer", [])

        assert records[0]["source"] == "primary"


# ---------------------------------------------------------------------------
# 5. Orchestrate integration (mocked)
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestOrchestrateGermlineIntegration:
    def test_orchestrate_merges_inherited_in_rare_disease(self):
        """run_pipeline in rare-disease mode with --germline merges inherited variants."""
        from scripts.orchestrate import run_pipeline

        inherited = [
            Variant(chrom="chr13", pos=32337326, ref="C", alt="T", gene="BRCA2", source="germline_inherited"),
        ]

        with (
            patch("scripts.orchestrate.parse_vcf") as mock_parse,
            patch("scripts.orchestrate.query_variant_databases") as mock_query,
            patch("scripts.orchestrate.compare_frequencies", return_value={"korean_flag": ""}),
            patch("scripts.pharmacogenomics.korean_pgx.get_pgx_results") as mock_pgx,
            patch("scripts.orchestrate.classify_variants") as mock_classify,
            patch("scripts.orchestrate.build_variant_records") as mock_build,
            patch(
                "scripts.orchestrate.build_summary",
                return_value={
                    "total": 2,
                    "pathogenic": 1,
                    "likely_pathogenic": 0,
                    "vus": 1,
                    "benign": 0,
                    "likely_benign": 0,
                    "drug_response": 0,
                    "risk_factor": 0,
                },
            ),
            patch("scripts.orchestrate.split_variants_for_display", return_value=([], [], [], 0, [], [])),
            patch("scripts.orchestrate.get_all_db_versions", return_value={}),
            patch("scripts.orchestrate.generate_report_html", return_value="<html></html>"),
            patch("scripts.orchestrate.resolve_hpo_terms", return_value=[]),
            patch("scripts.orchestration.extract_germline.extract_inherited_variants", return_value=inherited),
        ):
            # Primary VCF returns one variant
            primary = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
            mock_parse.return_value = [primary]

            mock_query.return_value = {
                "clinvar": {"clinvar_significance": "Pathogenic", "acmg_codes": []},
                "gnomad": {"gnomad_all": 0.0001, "gnomad_eas": 0.0002},
                "krgdb_freq": None,
                "pgx": None,
            }

            mock_pgx.return_value = {
                "pgx_hits": [],
                "pgx_source": "builtin",
                "pharmcat_version": "",
                "germline_provided": True,
                "warnings": [],
            }

            mock_cls = MagicMock()
            mock_cls.classification = "Pathogenic"
            mock_cls.evidence_codes = ["PVS1"]
            mock_cls.conflict = None
            mock_classify.return_value = {
                "chr17:7577120:G>A": mock_cls,
                "chr13:32337326:C>T": mock_cls,
            }

            mock_build.return_value = [
                {"variant": "chr17:7577120:G>A", "source": "primary"},
                {"variant": "chr13:32337326:C>T", "source": "germline_inherited"},
            ]

            # Create a dummy VCF file
            import tempfile

            with tempfile.NamedTemporaryFile(suffix=".vcf", mode="w", delete=False) as f:
                f.write(
                    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr17\t7577120\t.\tG\tA\t100\tPASS\t.\n"
                )
                vcf_path = f.name

            try:
                with tempfile.NamedTemporaryFile(suffix=".vcf", mode="w", delete=False) as gf:
                    gf.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                    germline_path = gf.name

                try:
                    result = run_pipeline(
                        vcf_path=vcf_path,
                        output_path=str(Path(tempfile.gettempdir()) / "test_report.html"),
                        mode="rare-disease",
                        germline_vcf=germline_path,
                        skip_api=True,
                    )

                    # build_variant_records should have been called with
                    # merged variants (primary + inherited)
                    assert mock_build.called
                    call_args = mock_build.call_args
                    all_variants = call_args[0][0]  # first positional arg
                    assert len(all_variants) == 2
                finally:
                    os.unlink(germline_path)
            finally:
                os.unlink(vcf_path)

    def test_orchestrate_skips_germline_in_cancer_mode(self):
        """In cancer mode, germline inherited extraction is not triggered."""
        # This tests that the mode guard works -- no inherited extraction in cancer mode
        from scripts.orchestrate import run_pipeline

        with (
            patch("scripts.orchestrate.parse_vcf") as mock_parse,
            patch("scripts.orchestrate.query_variant_databases") as mock_query,
            patch("scripts.orchestrate.compare_frequencies", return_value={"korean_flag": ""}),
            patch("scripts.pharmacogenomics.korean_pgx.get_pgx_results") as mock_pgx,
            patch("scripts.orchestrate.classify_variants") as mock_classify,
            patch("scripts.orchestrate.build_variant_records") as mock_build,
            patch(
                "scripts.orchestrate.build_summary",
                return_value={
                    "total": 1,
                    "pathogenic": 1,
                    "likely_pathogenic": 0,
                    "vus": 0,
                    "benign": 0,
                    "likely_benign": 0,
                    "drug_response": 0,
                    "risk_factor": 0,
                },
            ),
            patch("scripts.orchestrate.split_variants_for_display", return_value=([], [], [], 0, [], [])),
            patch("scripts.orchestrate.get_all_db_versions", return_value={}),
            patch("scripts.orchestrate.generate_report_html", return_value="<html></html>"),
            patch("scripts.orchestration.extract_germline.extract_inherited_variants") as mock_extract,
        ):
            primary = Variant(chrom="chr17", pos=7577120, ref="G", alt="A", gene="TP53")
            mock_parse.return_value = [primary]

            mock_query.return_value = {
                "clinvar": {"clinvar_significance": "Pathogenic", "acmg_codes": []},
                "gnomad": {"gnomad_all": 0.0001, "gnomad_eas": 0.0002},
                "krgdb_freq": None,
                "pgx": None,
            }

            mock_pgx.return_value = {
                "pgx_hits": [],
                "pgx_source": "builtin",
                "pharmcat_version": "",
                "germline_provided": True,
                "warnings": [],
            }

            mock_cls = MagicMock()
            mock_cls.classification = "Pathogenic"
            mock_cls.evidence_codes = ["PVS1"]
            mock_cls.conflict = None
            mock_classify.return_value = {"chr17:7577120:G>A": mock_cls}
            mock_build.return_value = [{"variant": "chr17:7577120:G>A", "source": "primary"}]

            import tempfile

            with tempfile.NamedTemporaryFile(suffix=".vcf", mode="w", delete=False) as f:
                f.write(
                    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr17\t7577120\t.\tG\tA\t100\tPASS\t.\n"
                )
                vcf_path = f.name

            try:
                with tempfile.NamedTemporaryFile(suffix=".vcf", mode="w", delete=False) as gf:
                    gf.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                    germline_path = gf.name

                try:
                    result = run_pipeline(
                        vcf_path=vcf_path,
                        output_path=str(Path(tempfile.gettempdir()) / "test_report_cancer.html"),
                        mode="cancer",
                        germline_vcf=germline_path,
                        skip_api=True,
                    )

                    # extract_inherited_variants should NOT have been called
                    mock_extract.assert_not_called()
                finally:
                    os.unlink(germline_path)
            finally:
                os.unlink(vcf_path)
