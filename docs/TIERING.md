# BIKO GenomeBoard — 변이 분류 및 Tiering 기준

**Version**: v2.5 (2026-04-21) · **Audience**: 변이 분류 결과를 검토·공유받는 연구자와 임상의

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
- `1 PVS + 1 PS` · `1 PVS + 2 PM` · `1 PVS + 1 PM + 1 PP` · `1 PVS + 2 PP`
- `2 PS` · `1 PS + 3 PM` · `1 PS + 2 PM + 2 PP` · `1 PS + 1 PM + 4 PP`

**Likely Pathogenic** (7개 조합):
- `1 PVS + 1 PM` · `1 PS + 1 PM` · `1 PS + 2 PM` · `1 PS + 2 PP`
- `3 PM` · `2 PM + 2 PP` · `1 PM + 4 PP`

**Benign**: `1 BA` · `2 BS` · · · **Likely Benign**: `1 BS + 1 BP` · `2 BP`

위 조합에 매치되지 않으면 **VUS**.

### 1.3 BIKO가 발화하는 evidence 코드와 임계값

| 코드 | 트리거 | 출처 |
|---|---|---|
| **BA1** | KOVA AF > 5% 또는 gnomAD AF > 5% | `thresholds.ba1` (config) |
| **BS1** | KOVA AF > 1% 또는 gnomAD AF > 1% | `thresholds.bs1` |
| **PM2_Supporting** | KOVA + gnomAD 모두 AF ≤ 0.1% (또는 absent) | `thresholds.pm2` · ClinGen SVI 2020 |
| **PS1** | 동일 아미노산 변화가 ClinVar 2★+ Pathogenic | `scripts/classification/evidence_collector.py` |
| **PM1** | `data/pm1_hotspot_domains.json` domain 내 missense | BIKO curated (ACMG PM1 기준) |
| **PM5** | 동일 residue 다른 missense가 ClinVar Pathogenic | `clinvar_pathogenic_positions` 캐시 |
| **PS2_Moderate** | 부모 양쪽 sequencing으로 de novo 확인 (trio + PED) | v2.3 trio workflow |
| **PM6_Supporting** | De novo "assumed" (parent confirmation 없음) | PMID 30311383 |
| **PP3** | REVEL > 0.644 (moderate) / > 0.932 (strong) · SpliceAI > 0.2 / > 0.5 | ClinGen SVI 2022 (Pejaver et al.) |
| **BP4** | REVEL < 0.183 / < 0.016 · SpliceAI < 0.1 | 동일 |
| **PP5** | ClinVar Pathogenic with review status | ClinVar |
| **BP6** | ClinVar Benign with review status | ClinVar |

**In silico fallback**: REVEL이 없는 경우 CADD Phred ≥ 25 + AlphaMissense likely_pathogenic 조합을 PP3 대용으로 사용 (`scripts/classification/in_silico.py`).

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

**결정 순서** (`scripts/somatic/amp_tiering.py::amp_assign_tier`):

```
1. Drug Response 분류 → Tier I (PGx actionable)
2. Risk Factor 분류 → Tier IV (cancer biomarker 아님)
3. Benign / Likely Benign → Tier IV (CIViC가 override 불가)
4. Pathogenic/LP + CIViC variant-specific Level A → Tier I
5. Pathogenic/LP + CIViC variant-specific Level B → Tier II
6. Pathogenic/LP + CIViC variant-specific Level C/D → Tier II
7. Pathogenic/LP + OncoKB gene Level 1/2 → Tier I
8. Pathogenic/LP + OncoKB gene (any level) → Tier II
9. Pathogenic/LP + non-cancer gene → Tier IV (incidental germline)
10. VUS + cancer gene + cancerhotspots_v2 hotspot → Tier II (VUS 승격)
11. VUS + cancer gene (비-hotspot) → Tier III
12. 그 외 → Tier IV
```

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

ACMG 분류와 **별개**로, AI Clinical Board가 검토할 변이를 선택하는 기준 (`scripts/clinical_board/variant_selector.py`):

### 4.1 Cancer mode — MUST reasons (무조건 포함)

| Reason key | 조건 | Priority |
|---|---|---|
| `P_LP` | 분류가 Pathogenic/Likely Pathogenic | 1 |
| `tier_1_2` | AMP Tier I/II | 2 |
| `hotspot_vus` | VUS on cancer gene + hotspot 매치 | 3 |
| `indel_hotspot_vus` | Hotspot 위치 VUS indel | 4 |
| `tsg_lof_vus` | TSG gene + LoF consequence VUS | 5 |
| `oncokb_gene_vus` | OncoKB level-가진 gene의 VUS | 6 |
| `drug_response` | PGx drug response | 7 |

### 4.2 Rare disease mode — MUST reasons

| Reason key | 조건 | Priority |
|---|---|---|
| `P_LP` | Pathogenic/Likely Pathogenic (consequence gate 미적용) | 1 |
| `VUS_hotspot` | VUS on known hotspot domain | 2 |
| `VUS_denovo_neurodev` | De novo VUS + DDG2P neurodev gene | 3 |
| `VUS_denovo_hpo_match` | De novo VUS + HPO-matched gene | 4 |
| `VUS_MMR_Lynch` | MLH1/MSH2/MSH6/PMS2/EPCAM의 protein-impacting VUS | 6 |
| `VUS_TSG_LoF` | LoF consequence on TSG | 7 |
| `VUS_HPO_match` | HPO-matched gene VUS | 8 |

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
- **Homozygote zygosity 미파악**: 단일 pathogenic 변이로 AR 질환을 단정하지 마십시오 (carrier status 가능). Board prompt에 경고 포함.
- **zygosity 정보 부재 시**: 일부 filter는 heterozygous 가정 — trio/quartet PED 입력 시만 정확한 de novo 판정.

---

## 9. 참고문헌

1. Richards S, et al. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. **Genet Med** 2015;17(5):405-24. PMID: 25741868
2. Li MM, et al. Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer. **J Mol Diagn** 2017;19(1):4-23. PMID: 27993330
3. Riggs ER, et al. Technical standards for the interpretation and reporting of constitutional copy-number variants. **Genet Med** 2020;22(2):245-257. PMID: 31690835
4. Pejaver V, et al. Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria. **AJHG** 2022;109(12):2163-2177. PMID: 36413997
5. Brnich SE, et al. Recommendations for application of the functional evidence PS3/BS3 criterion using the ACMG/AMP framework. **Genome Med** 2019;12(1):3. PMID: 31892348

---

*문의: BIKO GenomeBoard 팀 (juneh@kisti.re.kr). 이 문서는 v2.5 (2026-04-21) 기준이며, 가이드라인 업데이트 시 재개정됩니다.*
