# BIKO GenomeBoard — 변이 분류 및 Tiering 기준

**Version**: v2.5.1 (2026-04-21) · **Audience**: 변이 분류 결과를 검토·공유받는 연구자와 임상의

BIKO GenomeBoard는 **두 개의 독립된 분류 시스템**을 병렬로 적용합니다.

| 모드 | 시스템 | 카테고리 | 근거 가이드라인 |
|---|---|---|---|
| **Rare disease (germline)** | ACMG/AMP 2015 + ClinGen SVI updates | Pathogenic / Likely Pathogenic / VUS / Likely Benign / Benign | Richards et al. 2015 · ClinGen SVI 2018-2023 |
| **Cancer (somatic)** | AMP/ASCO/CAP 2017 | Tier I / II / III / IV | Li et al. 2017 |
| **CNV/SV** | ACMG 2020 | 5-class (L1=Pathogenic ... L5=Benign) | Riggs et al. 2020 |

모든 판정은 **결정론적 Python 엔진**에서 나옵니다. LLM은 판정 결과를 바꾸지 않으며, 문장으로 풀어 쓰기만 합니다 (curate-then-narrate 아키텍처).

---

## 1. Rare disease — ACMG 5-tier 분류

### 1.1 Evidence 코드와 강도

| Prefix | Strength | Direction | 의미 |
|---|---|---|---|
| **PVS1** | Very Strong | Pathogenic | Null variant (LoF) in LoF-intolerant gene |
| **PS1-4** | Strong | Pathogenic | PS1 aa 동일, PS2 de novo 확인, PS3 functional, PS4 case-control |
| **PM1-6** | Moderate | Pathogenic | PM1 hotspot, PM2 rare, PM3 in trans, PM4 length change, PM5 peer residue, PM6 de novo assumed |
| **PP1-5** | Supporting | Pathogenic | PP1 cosegregation, PP2 missense on LoF-constrained gene, PP3 in silico, PP4 phenotype match, PP5 reputable source |
| **BA1** | Stand-alone | Benign | AF > 5% |
| **BS1-4** | Strong | Benign | BS1 AF 1-5%, BS2 observed healthy, BS3 functional, BS4 non-segregation |
| **BP1-7** | Supporting | Benign | BP1 missense on truncating-mech gene, BP2 in trans with P, BP4 in silico, etc. |

**ClinGen SVI 2018 downgrade** (PMID 30311383): PS2를 PS2_Moderate로, PM6를 PM6_Supporting으로 강도 조정 (BIKO는 준수).

### 1.2 Combination rules (`data/acmg_rules.json`)

**Pathogenic** (8개 조합 중 하나 성립 시):
- `1 PVS + 1 PS` · `1 PVS + 2 PM` · `1 PVS + 1 PM + 1 PP` · `1 PVS + 2 PP` *¹*
- `2 PS` · `1 PS + 3 PM` · `1 PS + 2 PM + 2 PP` · `1 PS + 1 PM + 4 PP`

**Likely Pathogenic** (7개 조합):
- `1 PVS + 1 PM` · `1 PS + 1 PM` · `1 PS + 2 PM` · `1 PS + 2 PP`
- `3 PM` · `2 PM + 2 PP` · `1 PM + 4 PP`

**Benign**: `1 BA` · `2 BS` · · · **Likely Benign**: `1 BS + 1 BP` · `2 BP`

위 조합에 매치되지 않으면 **VUS**.

*¹ BIKO는 Richards 2015 원문 조합을 그대로 수록하고 있습니다. ClinGen SVI Bayesian 재해석(Tavtigian 2018, PMID 30192042)은 `1 PVS + 2 PP`만으로 P 도달에 의문을 제기하므로, 실제 LP/P 판정 시 `PP3_Strong` 여부 등을 함께 봐야 합니다.*

### 1.3 BIKO가 발화하는 evidence 코드와 임계값

| 코드 | 트리거 | 출처 |
|---|---|---|
| **BA1** | KOVA AF > 5% 또는 gnomAD AF > 5% | `thresholds.ba1` (config) |
| **BS1** | KOVA AF > 1% 또는 gnomAD AF > 1% *²* | `thresholds.bs1` |
| **PM2_Supporting** | KOVA + gnomAD 모두 AF ≤ 0.1% (또는 absent) *³* | `thresholds.pm2` |
| **PS1** | 동일 아미노산 변화가 ClinVar 2★+ Pathogenic | `scripts/classification/evidence_collector.py` |
| **PM1** | `data/pm1_hotspot_domains.json` domain 내 missense | BIKO curated (ACMG PM1 기준) |
| **PM5** | 동일 residue 다른 missense가 ClinVar Pathogenic | `clinvar_pathogenic_positions` 캐시 |
| **PS2_Moderate** | 부모 양쪽 sequencing으로 de novo 확인 (trio + PED) | v2.3 trio workflow |
| **PM6_Supporting** | De novo "assumed" (parent confirmation 없음) | PMID 30311383 |
| **PP3 / PP3_Moderate / PP3_Strong** | SpliceAI ≥ 0.2 / REVEL ≥ 0.644 → 현재 BIKO는 `PP3_Moderate` 발화. REVEL ≥ 0.932 또는 SpliceAI ≥ 0.5 → `PP3_Strong` *⁴* | `scripts/classification/in_silico.py` |
| **BP4 / BP4_Moderate / BP4_Strong** | REVEL ≤ 0.183 → `BP4_Moderate`, ≤ 0.016 → `BP4_Strong`, SpliceAI < 0.1 | 동일 |
| **PP5** | ClinVar Pathogenic with review status *⁵* | ClinVar |
| **BP6** | ClinVar Benign with review status *⁵* | ClinVar |

**In silico fallback**: REVEL이 없는 경우 CADD Phred ≥ 25 + AlphaMissense likely_pathogenic 조합을 PP3 대용으로 사용 (`scripts/classification/in_silico.py`).

*² BS1 1% cutoff는 common benign 범위의 generic threshold입니다. AR / AD / XL 질환 별 다른 penetrance-informed cutoff는 현재 설정되지 않습니다 (AR variant라도 AF ≥ 1%면 BS1 발화). 특수 질환 검토 시 수동 점검 필요.*

*³ BIKO의 PM2_Supporting threshold는 AF ≤ 0.001 (0.1%)로 `config.yaml::thresholds.pm2`에 설정되어 있습니다. **이는 ClinGen SVI 2020 권고 AF ≤ 0.0001 (0.01%)보다 10배 관대합니다.** 현재 설정의 rationale은 Korean pop에서 일부 rare variant가 gnomAD에 등재되지 않은 경우 소실을 방지하기 위함이나, false-positive PM2 risk가 있으므로 LP/P 판정은 PM2 외 다른 근거와 함께 교차검토해야 합니다. v2.6에서 AR/AD 기반 동적 cutoff 도입 검토 중.*

*⁴ **ClinGen SVI 2022 (Pejaver et al., PMID 36413997) 권고는 4-tier**: Supporting (REVEL ≥ 0.644), Moderate (≥ 0.773), Strong (≥ 0.932), Very Strong (≥ 0.994). **BIKO 현재 구현은 2-tier** (Moderate=0.644, Strong=0.932)로 Pejaver Supporting 임계값을 Moderate label로 사용 중 — 이는 실제 strength가 한 단계 과대평가됨을 의미합니다. v2.6에서 3-tier로 재정렬 예정 (Supporting=0.644 / Moderate=0.773 / Strong=0.932). 현재로선 `PP3_Moderate`로 보고된 값이 실제로는 Supporting 강도에 해당함에 유의하십시오.*

*⁵ ClinGen SVI 2018은 PP5/BP6를 deprecated로 권고합니다 (ClinVar entry를 별도 evidence로 인정하지 않음). BIKO는 현재 PP5/BP6를 evidence code로 유지하고 있으며, ClinVar "conflicting" 항목에 대해서는 §1.4 4-조건 gate 성립 시에만 override합니다. 즉 BIKO 정책: PP5/BP6 발화 O + ClinVar 신뢰 O (단, A4 경로로 좁게 override).*

### 1.4 ClinVar 충돌의 좁은 범위 override (A4 규칙)

ClinVar "conflicting interpretations" 항목을 BIKO 엔진이 LP/P로 덮어쓸 수 있는 단 하나의 경로는 **4조건 전원 성립** 시입니다 (`scripts/classification/acmg_engine.py::apply_hotspot_conflict_reconciliation`):

1. ClinVar review status = "conflicting"
2. BIKO 엔진이 독립적으로 LP 또는 P 도달
3. PM1 발화 (`pm1_hotspot_domains.json` 도메인 hotspot)
4. PM5 발화 (동일 residue 다른 missense가 ClinVar Pathogenic)

이 4가지가 모두 맞을 때만 override 되며, 이유가 `clinvar_override_reason`에 기록되어 리포트에 `.override-notice` 메크로로 표시됩니다. **이것이 BIKO가 공인 DB를 거스르는 유일한 경로**입니다.

---

## 2. Cancer somatic — AMP/ASCO/CAP 2017 Tiering

### 2.1 4-Tier 정의

| Tier | 의미 | 해당 예 |
|---|---|---|
| **Tier I** | Strong Clinical Significance | FDA-approved therapy 타겟, 강력한 predictive/prognostic 근거 |
| **Tier II** | Potential Clinical Significance | Clinical trial evidence, 다른 암종에서 확립된 변이 |
| **Tier III** | Unknown Clinical Significance | Cancer gene의 VUS, hotspot VUS |
| **Tier IV** | Benign or Likely Benign | Cancer 맥락에서 임상 의의 없음 (germline benign 포함) |

### 2.2 BIKO 통합 전략 — Strategy B (default)

BIKO는 AMP tiering에 3가지 전략을 지원합니다 (`amp_assign_tier(strategy=...)`):

- **Strategy A — CIViC priority only**: CIViC variant-specific evidence만으로 Tier I/II 지정. OncoKB는 fallback에만 사용. CIViC coverage가 완전한 gene-variant에만 권장.
- **Strategy B — Combined (default)**: CIViC variant-level 우선, OncoKB gene-level로 보강. Pathogenic/LP + CIViC Level C/D도 Tier II로 끌어올림. **현재 BIKO 기본값.**
- **Strategy C — OncoKB only**: CIViC 의존성 제거. Academic OncoKB key 없이도 동작하는 경로 (제한적).

기본(Strategy B) **결정 순서** (`scripts/somatic/amp_tiering.py::amp_assign_tier`):

```
1. Drug Response 분류 → Tier I (PGx actionable)
2. Risk Factor 분류 → Tier IV (cancer biomarker 아님)
3. Benign / Likely Benign → Tier IV (CIViC가 override 불가)
4. Pathogenic/LP + CIViC variant-specific Level A → Tier I
5. Pathogenic/LP + CIViC variant-specific Level B → Tier II
6. Pathogenic/LP + CIViC variant-specific Level C/D → Tier II
7. Pathogenic/LP + OncoKB gene Level 1/2 → Tier I
8. Pathogenic/LP + OncoKB gene (any level) → Tier II
9. Pathogenic/LP + non-cancer gene → Tier IV (incidental germline) *⁶*
10. VUS + cancer gene + cancerhotspots_v2 hotspot → Tier II (VUS 승격) *⁷*
11. VUS + cancer gene (비-hotspot) → Tier III
12. 그 외 → Tier IV
```

*⁶ AMP/ASCO/CAP 2017 원문은 germline incidental finding을 somatic tiering 범위 외로 권고합니다. BIKO는 일관성 위해 Tier IV로 포함하되 리포트에 "incidental germline — cancer relevance unclear"로 flag합니다.*

*⁷ BIKO는 Hotspot → Tier II **승격만** 수행하며, 반대로 hotspot이지만 functional evidence가 benign 방향인 경우의 **demotion 경로는 구현되지 않음**. 그런 변이는 ACMG classification 자체가 Benign/LB로 오면서 Tier IV가 되는 간접 경로만 존재.*

**CIViC Evidence Level**:
- A = FDA-approved
- B = Professional guidelines
- C = Case reports (multiple)
- D = Preclinical

**OncoKB Level**:
- 1 = FDA-recognized biomarker
- 2 = Standard care biomarker
- 3/4 = Investigational / Preclinical

### 2.3 Hotspot VUS 승격

`cancerhotspots_v2_single` DB에서 protein position이 hotspot으로 등재된 VUS는 Tier III에서 Tier II로 승격됩니다. 이는 Tier III의 signal 강화 용도이며, **ACMG PM1 (rare disease)와 다른 테이블**입니다 (PM1은 `pm1_hotspot_domains.json` 별도 사용).

### 2.4 TMB (Tumor Mutational Burden)

- 별도 지표 (ACMG/AMP tiering과 독립)
- Low: < 10 mut/Mb · Intermediate: 10-20 · High: ≥ 20 mut/Mb
- 패널 크기 `panel_size_mb` 기준 (WGS 기본 33 Mb)
- Immunotherapy 적격성 근거로만 사용

---

## 3. CNV / SV — ACMG 2020

`scripts/somatic/amp_tiering.py`와 별도 엔진. Section 2 (dosage sensitivity) + Section 3 (number of protein-coding genes) + Section 4 (inheritance) + Section 5 (case-control) 포인트를 합산해 5-tier 분류:
- **L1 (Pathogenic)** · **L2 (Likely Pathogenic)** · **L3 (VUS)** · **L4 (Likely Benign)** · **L5 (Benign)**

입력은 AnnotSV 결과 + Canvas/Manta VCF.

---

## 4. Clinical Board에 올라가는 변이 선택 기준

ACMG 분류와 **별개**로, AI Clinical Board가 검토할 변이를 선택하는 기준 (`scripts/clinical_board/variant_selector.py::_REASON_PRIORITY`).

### 4.1 통합 priority 표 (cancer + rare disease 공용)

정렬은 **단일 priority 테이블**에서 수행됩니다. Cancer/rare-disease에 따라 어느 reason이 활성화되는지가 다를 뿐, 동일한 숫자는 동일 순위입니다.

| Priority | Reason key | 조건 | Mode |
|---|---|---|---|
| **0** | `P_LP` | 분류가 Pathogenic / Likely Pathogenic (consequence gate 미적용) | 공용 |
| **1** | `Tier_I` | AMP Tier I | Cancer |
| **2** | `Tier_II` | AMP Tier II | Cancer |
| **3** | `Tier_III_hotspot` | VUS + cancer gene + cancerhotspots_v2 매치 | Cancer |
| **4** | `Tier_III_oncokb_gene` | VUS + OncoKB-등재 cancer gene (비-hotspot) | Cancer |
| **5** | `VUS_hotspot` | VUS + ACMG PM1 hotspot domain | 공용 |
| **5** | `VUS_denovo_neurodev` | De novo VUS + DDG2P neurodev gene 또는 constrained gene (pLI ≥ 0.9 AND missense-Z ≥ 3.09) | Rare |
| **5** | `VUS_denovo_splice` | De novo VUS + SpliceAI delta ≥ 0.2 | Rare |
| **6** | `VUS_MMR_Lynch` | MLH1/MSH2/MSH6/PMS2/EPCAM의 protein-impacting VUS *⁸* | Rare |
| **7** | `VUS_TSG_LoF` | TSG gene + LoF consequence VUS | Cancer |
| **8** | `VUS_indel_hotspot` | Hotspot 위치 VUS indel | Cancer |
| **9** | `VUS_HPO_match` | HPO-matched gene VUS | Rare |

**Priority 5 공유**: de novo neurodev, de novo splice, hotspot 모두 우선순위 동일 — 리포트 순서는 tiebreak 키(아래 4.3)로 결정.

*⁸ MMR/Lynch VUS가 de novo(5)보다 후순위인 이유: de novo 변이는 임상 판단에 즉시 중요한 새 사건이지만, MMR VUS는 Lynch syndrome 모니터링 대상 후보로 기록 목적이 강함. Phenotype-first 원칙(§5)과 정합.*

**Consequence gate (B1)**: Non-coding VUS (intronic, UTR, upstream, downstream, synonymous)은 `P_LP` 이외 모든 MUST/MAY reason에서 **제외**됩니다. Splice-region / synonymous는 SpliceAI delta ≥ 0.2일 때 예외적으로 통과.

### 4.3 Soft caps + 우선순위 정렬

선택 결과는 `(board_admitted, classification_rank, hpo_score, has_gene, variant_id)` 복합 키로 정렬되어 리포트 페이지 1에 배치됩니다 (`scripts/clinical_board/runner.py::_variant_sort_key`):

1. Board가 선택한 변이 먼저
2. 분류 순위: P > LP > Drug Response > Risk Factor > VUS > LB > B
3. HPO 점수 높은 순
4. Gene 있는 것 우선
5. 변이 ID (deterministic tiebreak)

---

## 5. Phenotype-first Primary Diagnosis 규칙 (v2.5)

Rare disease의 **Primary diagnosis** 선정에서 Board Chair는 다음 규칙을 의무적으로 따릅니다:

> **`PRIMARY DIAGNOSIS CANDIDATES (phenotype-matched)` 섹션에 변이가 있다면**, primary_diagnosis는 **반드시** 그 섹션에서 선택합니다. 섹션에 없는 Pathogenic/LP 변이는 `incidental finding` 또는 `secondary finding`으로 분리 보고합니다.

이유: VUS인 phenotype-matched gene의 de novo 변이가 phenotype-무관 LP carrier finding보다 clinically more relevant. LLM의 학습된 gene-disease 연관 prior (예: CFTR→Cystic Fibrosis)가 실제 제출 HPO를 override하는 실패 모드를 방지합니다.

---

## 6. 한국인 인구집단 빈도

BIKO v2.5는 **KOVA v7** (Korean Variant Archive, 43.3M variants) 단일 소스를 primary Korean AF로 사용합니다.

| 소스 | 역할 | BA1/BS1/PM2 판정 |
|---|---|---|
| **KOVA v7** | Primary Korean | Korean AF 사용 |
| **gnomAD EAS** | Secondary (비교용, Korean enrichment 계산) | 참고 |
| **gnomAD ALL** | Global reference | 참고 |

**Korean enrichment ratio** = `kova_af / gnomad_eas_af`. 이 값이 크면 한국인 특이 변이 flag.

**Homozygote count** (KOVA 전용 필드): AR disease 판정 시 healthy homozygote가 KOVA에 관측되면 BS2 근거로 사용 가능.

이전의 KRGDB / Korea4K / NARD2는 2026-04-21부로 제거되었습니다 (데이터 접근 불가 확인).

---

## 7. 리포트에서 이 정보를 확인하는 곳

| 섹션 | 내용 |
|---|---|
| **Variant Detail Card** | ACMG classification, evidence codes, HGVSc/p, KOVA/gnomAD AF, in silico scores |
| **Override Notice** (해당 시) | ClinVar conflict A4 override 이유 (PM1/PM5/PMID 포함) |
| **De novo Badge** | PS2 확인 / PM6 assumed / DDG2P neurodev / SpliceAI splice |
| **Tier Badge** (cancer) | Tier I-IV + CIViC level + OncoKB level |
| **Candidate Gene Ranking** | Board가 선택한 변이 상위 15개 (분류 + HPO score 정렬) |
| **Clinical Board Opinion** | 4명의 도메인 에이전트 + Board Chair 최종 의견 (LLM 생성) |
| **Treatment Options** (cancer) | Curated: OncoKB + CIViC curator에서만 나옴, `(curated_id, variant_key)` 쌍 검증됨 |

---

## 8. 한계와 주의사항

- **연구용 참고 자료**: BIKO 리포트는 CLIA 임상 보고가 아닙니다. 판독은 담당 연구자 또는 임상의의 독립 판단을 요합니다.
- **결정론적 엔진 ≠ 정답**: ACMG 가이드라인 자체가 발전 중이며, ClinGen Expert Panel 업데이트를 주기적으로 반영해야 합니다.
- **LLM narrative는 기록만**: Board 의견이 분류를 바꾸지 않으며, narrative scrubber가 curated evidence 외 약물 언급을 삭제합니다.
- **Homozygote zygosity 미파악**: 단일 pathogenic 변이로 AR 질환을 단정하지 마십시오. Board prompt에 경고 포함.
- **Phasing 처리**: BIKO는 trio-based phasing(PED의 parent genotype 활용)만 지원하며 read-based phasing은 미지원. Trans/cis 판정은 compound het 판정의 한계이므로 PM3/BP2는 현재 발화하지 않음.
- **KOVA v7 단일 의존 리스크**: Korean AF를 KOVA v7 하나로 집약해 version drift (향후 v8 업데이트 시 threshold drift), regional sub-population (특정 지역 cohort 과소대표) 위험이 있습니다. gnomAD EAS와 항상 교차확인 권장.
- **ACMG SF 보조 findings**: `incidental finding` 표기는 BIKO 고유 표현이며, 정식 "secondary findings" 리포팅은 별도 정책(ACMG SF v3.2 gene list + 환자 opt-in)에 따라 운영해야 합니다. BIKO 리포트는 ACMG SF 자동 스크리닝을 수행하지 않습니다.

---

## 9. 참고문헌

1. Richards S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. **Genet Med** 2015;17(5):405-24. PMID: 25741868
2. Li MM, et al. Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer. **J Mol Diagn** 2017;19(1):4-23. PMID: 27993330
3. Riggs ER, et al. Technical standards for the interpretation and reporting of constitutional copy-number variants. **Genet Med** 2020;22(2):245-257. PMID: 31690835
4. Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. **AJHG** 2022;109(12):2163-2177. PMID: 36413997
5. Brnich SE, et al. Recommendations for application of the functional evidence PS3/BS3 criterion using the ACMG/AMP framework. **Genome Med** 2019;12(1):3. PMID: 31892348

---

*문의: BIKO GenomeBoard 팀 (juneh@kisti.re.kr). 이 문서는 v2.5 (2026-04-21) 기준이며, 가이드라인 업데이트 시 재개정됩니다.*
