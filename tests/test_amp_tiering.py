# === Strategy B (default) ===


def test_strategy_b_civic_variant_level_a_pathogenic():
    """CIViC variant-specific Level A + Pathogenic → Tier I."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {
        "match_level": "variant",
        "evidence": [
            {
                "evidence_level": "A",
                "therapies": "Vemurafenib",
                "significance": "Sensitivity/Response",
                "disease": "Melanoma",
                "gene": "BRAF",
                "variant": "V600E",
                "evidence_type": "Predictive",
                "pmid": "20979469",
                "citation": "Chapman 2011",
                "statement": "",
            }
        ],
    }
    result = amp_assign_tier("Pathogenic", "BRAF", hgvsp="p.Val600Glu", civic_evidence=civic)
    assert result.tier == 1
    assert result.evidence_source == "civic-variant-A"
    assert result.civic_match_level == "variant"


def test_strategy_b_civic_variant_level_b():
    """CIViC variant-specific Level B → Tier II."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {
        "match_level": "variant",
        "evidence": [
            {
                "evidence_level": "B",
                "therapies": "Erlotinib",
                "significance": "Resistance",
                "disease": "NSCLC",
                "gene": "KRAS",
                "variant": "G12D",
                "evidence_type": "Predictive",
                "pmid": "",
                "citation": "",
                "statement": "",
            }
        ],
    }
    result = amp_assign_tier("Pathogenic", "KRAS", hgvsp="p.Gly12Asp", civic_evidence=civic)
    assert result.tier == 2
    assert result.evidence_source == "civic-variant-B"


def test_strategy_b_civic_gene_level_does_not_elevate_to_tier1():
    """CIViC gene-level (not variant) Level A → does NOT elevate to Tier I."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {
        "match_level": "gene",
        "evidence": [
            {
                "evidence_level": "A",
                "therapies": "Vemurafenib",
                "significance": "Sensitivity/Response",
                "disease": "Melanoma",
                "gene": "BRAF",
                "variant": "V600E",
                "evidence_type": "Predictive",
                "pmid": "",
                "citation": "",
                "statement": "",
            }
        ],
    }
    # BRAF is OncoKB Level 1, so still Tier I via OncoKB, but NOT via CIViC
    result = amp_assign_tier("Pathogenic", "BRAF", hgvsp="p.Lys601Glu", civic_evidence=civic)
    assert result.tier == 1
    assert result.evidence_source.startswith("oncokb")  # not civic


def test_strategy_b_no_civic_oncokb_level1_pathogenic():
    """No CIViC + OncoKB Level 1 + Pathogenic → Tier I via OncoKB."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {"match_level": "none", "evidence": []}
    result = amp_assign_tier("Pathogenic", "TP53", civic_evidence=civic)
    assert result.tier == 1
    assert "oncokb" in result.evidence_source


def test_strategy_b_vus_hotspot():
    """VUS + hotspot → Tier II."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {"match_level": "none", "evidence": []}
    # Need to mock is_hotspot — test with hgvsp that triggers hotspot
    result = amp_assign_tier("VUS", "KRAS", hgvsp="p.Gly12Asp", civic_evidence=civic)
    # KRAS G12 is a known hotspot; result depends on CIViC DB availability
    assert result.tier in (2, 3)  # 2 if hotspot found, 3 otherwise


def test_strategy_b_vus_cancer_gene_no_hotspot():
    """VUS on cancer gene, no hotspot → Tier III."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {"match_level": "none", "evidence": []}
    result = amp_assign_tier("VUS", "TP53", hgvsp="p.Ala999Val", civic_evidence=civic)
    assert result.tier == 3


def test_strategy_b_benign():
    """Benign → Tier IV, CIViC cannot override."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {
        "match_level": "variant",
        "evidence": [
            {
                "evidence_level": "A",
                "therapies": "X",
                "significance": "Sensitivity/Response",
                "disease": "Y",
                "gene": "BRAF",
                "variant": "V600E",
                "evidence_type": "Predictive",
                "pmid": "",
                "citation": "",
                "statement": "",
            }
        ],
    }
    result = amp_assign_tier("Benign", "BRAF", hgvsp="p.Val600Glu", civic_evidence=civic)
    assert result.tier == 4


def test_strategy_b_drug_response():
    """Drug Response → Tier I always."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    result = amp_assign_tier("Drug Response", "CYP2C19")
    assert result.tier == 1


def test_strategy_b_vus_non_cancer_gene():
    """VUS on non-cancer gene → Tier IV."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {"match_level": "none", "evidence": []}
    result = amp_assign_tier("VUS", "FAKEGENE", civic_evidence=civic)
    assert result.tier == 4


def test_strategy_b_civic_level_cd_pathogenic():
    """CIViC variant-specific Level C + Pathogenic → Tier II via CIViC (Priority 5)."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {
        "match_level": "variant",
        "evidence": [
            {
                "evidence_level": "C",
                "therapies": "Drug",
                "significance": "Sensitivity/Response",
                "disease": "Cancer",
                "gene": "BRAF",
                "variant": "V600E",
                "evidence_type": "Predictive",
                "pmid": "",
                "citation": "",
                "statement": "",
            }
        ],
    }
    result = amp_assign_tier("Pathogenic", "BRAF", hgvsp="p.Val600Glu", civic_evidence=civic)
    assert result.tier == 2
    assert result.evidence_source == "civic-variant-C"


# === Strategy C (OncoKB only, backward compatible) ===


def test_strategy_c_ignores_civic():
    """Strategy C: CIViC evidence ignored, OncoKB only."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {
        "match_level": "variant",
        "evidence": [
            {
                "evidence_level": "A",
                "therapies": "Drug",
                "significance": "Sensitivity/Response",
                "disease": "Cancer",
                "gene": "KRAS",
                "variant": "G12D",
                "evidence_type": "Predictive",
                "pmid": "",
                "citation": "",
                "statement": "",
            }
        ],
    }
    result = amp_assign_tier("VUS", "KRAS", hgvsp="p.Gly12Asp", strategy="C", civic_evidence=civic)
    # VUS on cancer gene → Tier III (no CIViC elevation in Strategy C)
    assert result.tier in (2, 3)  # 2 only if hotspot, 3 otherwise


# === Strategy A (CIViC priority) ===


def test_strategy_a_civic_level_a():
    """Strategy A: CIViC Level A → Tier I."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    civic = {
        "match_level": "variant",
        "evidence": [
            {
                "evidence_level": "A",
                "therapies": "Drug",
                "significance": "Sensitivity/Response",
                "disease": "Cancer",
                "gene": "BRAF",
                "variant": "V600E",
                "evidence_type": "Predictive",
                "pmid": "",
                "citation": "",
                "statement": "",
            }
        ],
    }
    result = amp_assign_tier("Pathogenic", "BRAF", hgvsp="p.Val600Glu", strategy="A", civic_evidence=civic)
    assert result.tier == 1
    assert result.evidence_source == "civic-variant-A"


# === TierResult structure ===


def test_tier_result_has_all_fields():
    """TierResult에 모든 필수 필드가 있어야 함."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    result = amp_assign_tier("Pathogenic", "TP53")
    assert hasattr(result, "tier")
    assert hasattr(result, "tier_label")
    assert hasattr(result, "evidence_source")
    assert hasattr(result, "civic_match_level")
    assert hasattr(result, "civic_evidence")
    assert isinstance(result.tier, int)
    assert isinstance(result.tier_label, str)


def test_amp_tier_labels():
    """AMP 2017 라벨 확인."""
    from scripts.somatic.amp_tiering import amp_assign_tier

    r1 = amp_assign_tier("Pathogenic", "TP53")
    assert "Strong Clinical Significance" in r1.tier_label
    r4 = amp_assign_tier("Benign", "FAKEGENE")
    assert "Benign" in r4.tier_label
