# AI Clinical Board v2 — Design Spec

**Date:** 2026-04-14
**Status:** Approved
**Baseline:** v0.5.0 (tag)

---

## 1. Overview

AI Clinical Board v2는 v0.5.0의 AI Board를 3개 축으로 개선한다:

1. **Grounded Prompting** — 에이전트별 도메인 특화 컨텍스트 주입 + 모드별 에이전트 역할 분리
2. **임상 노트 입력** — 자유 텍스트 임상 정보를 옵셔널 입력으로 수용
3. **지식 베이스** — Board 판단 이력 축적 및 재활용 (하이브리드 Wiki + SQLite)

### 핵심 원칙

- **온프레미스 전용**: 외부 API 호출 없음. 모든 데이터는 로컬 DB에서 조회.
- **결정적 분류 불변**: ACMG/AMP 분류 엔진은 AI Board에 의해 변경되지 않음.
- **Hallucination 최소화**: LLM에게 "추론"을 시키지 않고 "정리"를 시킴. 사실 근거를 프롬프트에 주입.
- **일관성 우선**: 동일 변이 + 동일 맥락 → 동일 결과를 목표.

---

## 2. Phase 1: Grounded Prompting + 모드별 에이전트 분리

### 2.1 에이전트 구성 — 모드별 분리

| 슬롯 | Rare Disease 모드 | Cancer 모드 |
|------|------------------|-------------|
| Agent 1 | Variant Pathologist (단백질 기능 영향) | Therapeutic Target Analyst (druggable target, resistance mutation) |
| Agent 2 | Disease Geneticist (변이-질환 연관, 감별진단) | Tumor Genomics Specialist (driver vs passenger, clonal 의미) |
| Agent 3 | PGx Specialist (약물 대사, 동일) | PGx Specialist (항암제 대사 포함) |
| Agent 4 | Literature Analyst (임상 근거 종합) | Clinical Evidence Analyst (치료 반응, 가이드라인 참조) |
| Chair | Board Chair — 진단 종합 | Board Chair — 치료 전략 종합 |

### 2.2 프롬프트 구조 변경

**현재 (v0.5.0):**
```
[System Prompt] + [공통 브리핑 ~3K tokens]
→ 모든 에이전트가 동일한 입력
```

**v2 (토큰 예산 배분: MedGemma 27B = 8K context):**
```
[System Prompt ~1K] + [공통 브리핑 ~2K] + [Domain Sheet ~4K] + [Prior Knowledge ~1K]
→ 에이전트별로 다른 도메인 데이터, 각 구간에 hard char limit 적용
```

### 2.3 도메인 시트 — Rare Disease 모드

| 에이전트 | Domain Sheet 내용 |
|---------|-------------------|
| Variant Pathologist | ClinVar full entry (significance, review status, condition), protein domain annotation, in silico 상세 점수 + 판정 근거 (REVEL, CADD, AlphaMissense, **SpliceAI**), SIFT/PolyPhen raw scores |
| Disease Geneticist | OMIM genemap 전체 phenotype 목록 + inheritance pattern, GeneReviews 질환 설명 텍스트, HPO 매칭 상세 (gene-phenotype 연관 목록), Orphanet 유병률 |
| PGx Specialist | CPIC 가이드라인 텍스트 (drug-gene 권고사항), Korean vs Western prevalence 상세, star allele 해석 |
| Literature Analyst | CIViC evidence 전문 (drug, evidence level, direction, PMID), GeneReviews NBK/PMID 목록, ClinGen validity classification + date |

### 2.4 도메인 시트 — Cancer 모드

| 에이전트 | Domain Sheet 내용 |
|---------|-------------------|
| Therapeutic Target Analyst | ClinVar full entry, protein domain annotation, in silico 상세, CIViC variant-level drug evidence (drug, level, direction) |
| Tumor Genomics Specialist | OncoKB gene info (level, cancer type), CIViC gene summary, TMB 데이터, cancer hotspot 정보, **VAF (variant allele frequency), tumor purity estimate, co-occurring mutations** |
| PGx Specialist | CPIC 가이드라인 (항암제 포함), Korean vs Western prevalence, star allele 해석 |
| Clinical Evidence Analyst | CIViC evidence 전문 (drug, evidence level, PMID), NCCN 스타일 치료 가이드라인 텍스트 (지식 베이스에서 로드) |

### 2.5 Cancer Board Chair 출력 스키마

v0.5.0 `BoardOpinion` (진단 중심)과 별도로 `CancerBoardOpinion` (치료 중심)을 추가:

```python
@dataclass
class CancerBoardOpinion:
    therapeutic_implications: str          # 치료적 함의 (primary_diagnosis 대체)
    therapeutic_evidence: str              # 근거 요약
    treatment_options: List[dict]          # [{drug, evidence_level, line, PMID, notes, resistance_notes}]
    actionable_findings: List[str]         # 치료에 영향을 미치는 소견
    clinical_actions: List[str]            # 구체적 임상 조치
    immunotherapy_eligibility: str         # TMB 기반 면역치료 적합성 (MSI는 향후)
    agent_opinions: List[AgentOpinion]
    agent_consensus: str
    dissenting_opinions: List[str]
    monitoring_plan: List[str]             # 반응 모니터링 계획
    confidence: str
    disclaimer: str
```

Rare Disease 모드의 `BoardOpinion`은 현재 구조 유지.

### 2.6 Temperature 변경

모든 에이전트 temperature: 0.3 → **0.1** (일관성 강화).
temperature는 `config.yaml`의 `clinical_board.temperature`에서 읽고, `base.py`와 `board_chair.py`의 하드코딩된 값을 config 참조로 교체한다.

### 2.7 구현 대상 파일

**신규:**
- `scripts/clinical_board/domain_sheets.py` — 에이전트별 + 모드별 도메인 시트 빌더
- `scripts/clinical_board/agents/therapeutic_target.py` — Cancer Agent 1
- `scripts/clinical_board/agents/tumor_genomics.py` — Cancer Agent 2
- `scripts/clinical_board/agents/clinical_evidence.py` — Cancer Agent 4

**수정:**
- `scripts/clinical_board/agents/base.py` — `_build_prompt()`에 `domain_sheet` 파라미터 추가 (기본값 `""`), temperature를 config에서 읽도록 변경
- `scripts/clinical_board/agents/board_chair.py` — 모드별 프롬프트 + 출력 스키마 분기, temperature config 참조
- `scripts/clinical_board/models.py` — `CancerBoardOpinion` 추가, disclaimer를 고정 텍스트로 교체 (section 7.1, 영문/한글 모두 통일)
- `scripts/clinical_board/runner.py` — 모드에 따라 에이전트 셋 선택, 반환 타입을 `Optional[Union[BoardOpinion, CancerBoardOpinion]]`로 변경
- `scripts/clinical_board/render.py` — `isinstance` 분기로 BoardOpinion / CancerBoardOpinion 각각 렌더링
- `scripts/clinical_board/ollama_client.py` — generate() default temperature를 config 기반으로 변경
- `templates/cancer/report.html` — 치료 중심 AI Board 섹션

---

## 3. Phase 2: 임상 노트 입력

### 3.1 CLI 인터페이스

```bash
# 직접 텍스트
python scripts/orchestrate.py sample.vcf --clinical-note "55세 남성, 대장암 3기, FOLFOX 후 재발"

# 파일 경로
python scripts/orchestrate.py sample.vcf --clinical-note-file clinical_history.txt

# HPO와 함께 (둘 다 옵셔널, 독립적)
python scripts/orchestrate.py patient.vcf --mode rare-disease \
  --hpo HP:0001250,HP:0001263 \
  --clinical-note "6세 여아, 18개월부터 발달지연, 최근 경련 시작"
```

### 3.2 데이터 흐름

```
CLI 입력 (--clinical-note / --clinical-note-file)
  → orchestrate.py: report_data["clinical_note"]에 저장
  → case_briefing.py: 공통 브리핑에 == CLINICAL CONTEXT == 섹션 추가
  → 모든 에이전트가 공통 브리핑을 통해 임상 맥락 참조
  → 리포트 HTML에는 원문 포함하지 않음 (개인정보 보호)
```

### 3.3 제약 사항

- 임상 노트는 AI Board 에이전트의 **입력으로만** 사용. 결정적 분류 엔진에 영향 없음.
- 리포트에 원문 미포함 (로컬 처리이므로 외부 전송은 없지만 출력물에 환자 정보 잔류 방지).
- 최대 1500자. 초과 시 트리밍 + 경고 로그. 한글 멀티바이트 문자 경계에서 절단하지 않도록 처리.
- 한글/영문 모두 가능. 에이전트 프롬프트에 "Clinical note may be in Korean; interpret accordingly" 안내.
- 출력 언어는 `--board-lang` 설정을 따름 (입력 언어와 무관).
- `--clinical-note`와 `--clinical-note-file` 동시 제공 시 에러 (mutually exclusive).

### 3.4 구현 대상 파일

**수정:**
- `scripts/orchestrate.py` — `--clinical-note`, `--clinical-note-file` argparse 추가
- `scripts/clinical_board/case_briefing.py` — `_build_clinical_note_section()` 추가
- `scripts/clinical_board/agents/base.py` — 프롬프트에 multilingual note 안내 추가

---

## 4. Phase 3: 지식 베이스 (하이브리드 Wiki + SQLite)

### 4.1 디렉토리 구조

```
data/knowledge_base/
├── kb.sqlite3                      ← 구조화 데이터 (판단 이력, 통계, 검색)
├── index.md                        ← 전체 Wiki 인덱스 (자동 관리)
├── log.md                          ← 케이스 처리 이력 (append-only)
│
├── cancer/
│   ├── genes/                      ← flat, 유전자별 종합 해석
│   │   ├── TP53.md
│   │   ├── BRCA2.md
│   │   └── ...
│   └── treatments/                 ← 암종별 치료 가이드라인 (사전 구축)
│       ├── colorectal_cancer.md
│       ├── lung_cancer_nsclc.md
│       ├── breast_cancer.md
│       └── ...
│
└── rare-disease/
    ├── genes/                      ← flat, 유전자별 종합 해석
    │   ├── TP53.md
    │   ├── CFTR.md
    │   └── ...
    └── diagnoses/                  ← 질환별 변이-진단 패턴
        ├── li_fraumeni.md
        ├── cystic_fibrosis.md
        └── ...
```

### 4.2 SQLite 스키마 (`kb.sqlite3`)

```sql
CREATE TABLE board_decisions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id TEXT NOT NULL,
    mode TEXT NOT NULL,                        -- 'cancer' | 'rare-disease'
    date TEXT NOT NULL,                        -- ISO 8601
    gene TEXT NOT NULL,
    variant TEXT NOT NULL,                     -- chr:pos:ref:alt
    hgvsc TEXT,
    hgvsp TEXT,
    classification TEXT,                       -- ACMG classification
    board_diagnosis TEXT,                      -- Board의 primary diagnosis / therapeutic implication
    board_confidence TEXT,                     -- high / moderate / low
    clinical_context_summary TEXT,             -- 구조화된 요약만 저장 (e.g., "age:50s|sex:M|dx:lung_ca|family_hx:none") — 원문 미저장, 재식별 방지
    agent_consensus TEXT,                      -- unanimous / majority / split
    raw_opinion_json TEXT                      -- BoardOpinion / CancerBoardOpinion 전체 직렬화
);

CREATE INDEX idx_bd_variant ON board_decisions(gene, variant, mode);
CREATE INDEX idx_bd_gene ON board_decisions(gene, mode);
CREATE INDEX idx_bd_date ON board_decisions(date);
CREATE INDEX idx_bd_hgvsp ON board_decisions(gene, hgvsp, mode);

-- 변이별 집계 뷰
CREATE VIEW variant_stats AS
SELECT
    gene, variant, mode,
    COUNT(*) as total_cases,
    SUM(CASE WHEN board_confidence = 'high' THEN 1 ELSE 0 END) as high_confidence_count,
    GROUP_CONCAT(DISTINCT board_diagnosis) as diagnoses_seen,
    MAX(date) as last_seen
FROM board_decisions
GROUP BY gene, variant, mode;
```

### 4.3 Wiki 페이지 관리

Wiki 페이지는 **LLM이 생성/편집하지 않음**. 코드가 SQLite 데이터를 기반으로 마크다운을 자동 생성한다.

```
Board 실행 완료
  → SQLite에 board_decisions INSERT (변이별 1행씩)
  → 해당 유전자의 .md 파일을 코드가 재생성:
     - 유전자 기본 정보
     - 해당 유전자 변이들의 판단 이력 테이블 (상위 10개)
     - 맥락별 해석 패턴 요약
  → log.md에 케이스 append
  → index.md 업데이트
```

이점:
- Wiki 내용이 항상 SQLite와 동기화 (single source of truth = SQLite)
- LLM hallucination으로 Wiki가 오염될 위험 없음
- 사람이 .md 파일을 직접 읽고 검토 가능

### 4.4 새 케이스 분석 시 활용 흐름

```
새 케이스 (TP53 R175H 포함, cancer 모드)
  → kb_query.py: SQLite에서 variant_stats 조회
     "cancer 모드에서 4건 분석됨: 3건 somatic driver (high), 1건 resistance (moderate)"
  → 프롬프트의 [Prior Knowledge ~1K] 구간에 별도 섹션으로 주입:
     "== PRIOR BOARD KNOWLEDGE ==
      이전 4건에서 이 변이를 분석함:
      - Somatic driver, platinum sensitivity (3건, high confidence)
      - Resistance mutation context (1건, moderate confidence)
      ⚠ Prior knowledge는 참고 자료일 뿐입니다. 현재 케이스의 임상 맥락이 우선합니다.
      현재 케이스의 맥락을 고려하여 독립적으로 판단하시오."
  → 에이전트가 prior knowledge를 참고하되 anchoring bias 없이 분석
```

Prior knowledge는 Domain Sheet와 별개의 프롬프트 구간이다 (section 2.2의 `[Prior Knowledge ~1K]`).
`kb_query.py`가 텍스트를 생성하고, `runner.py`가 에이전트 호출 시 별도 인자로 전달한다.

Wiki `.md` 파일은 **사람이 검토하기 위한 용도**이며, 에이전트 프롬프트에는 주입되지 않는다.
에이전트에게는 SQLite 쿼리 결과를 정형화한 텍스트만 전달한다.

### 4.5 판단 충돌 처리

모든 판단을 보존한다. 충돌 시:

```
EGFR L858R (cancer 모드, SQLite 기록):

| 케이스 | 날짜 | 판단 | Confidence | 맥락 |
|--------|------|------|-----------|------|
| S001 | 2026-04-14 | Erlotinib sensitive | high | 비소세포폐암, 1차 치료 |
| S015 | 2026-04-20 | Osimertinib preferred | high | 비소세포폐암, T790M 동반 |

→ Prior knowledge로 주입 시:
  "이전 2건에서 이 변이를 분석함:
   - 1차 치료 맥락 → Erlotinib sensitive (1건, high)
   - T790M 동반 맥락 → Osimertinib preferred (1건, high)
   현재 케이스의 맥락을 고려하여 판단하시오."
```

cancer/rare-disease 모드별로 분리 검색하므로, 다른 모드의 판단이 섞이지 않는다. 같은 모드 내에서의 맥락 차이만 prior knowledge에 표시.

### 4.6 Cancer treatments/ 초기 데이터

NCCN 스타일 가이드라인 요약을 사전 구축한다. 공개 가이드라인 및 CIViC 데이터 기반.

초기 대상 (주요 암종 10개):
- Colorectal cancer
- Non-small cell lung cancer
- Breast cancer
- Prostate cancer
- Melanoma
- Renal cell carcinoma
- Hepatocellular carcinoma
- Pancreatic cancer
- Gastric cancer
- Ovarian cancer

각 파일은 주요 치료 라인, 표적 유전자/변이, 약물, evidence level을 포함.
모든 KB 가이드라인 파일에는 `source:` frontmatter 필드를 포함한다 (e.g., `source: CIViC + public guidelines`).
에이전트 프롬프트에 "KB 요약을 가이드라인 수준 근거로 인용하지 마시오" 지시를 포함.

### 4.7 구현 대상 파일

**신규:**
- `scripts/clinical_board/knowledge_base.py` — KB 매니저 (SQLite CRUD + Wiki 자동 생성)
- `scripts/clinical_board/kb_query.py` — 새 케이스용 prior knowledge 검색
- `scripts/db/build_kb_db.py` — SQLite 스키마 초기화
- `data/knowledge_base/cancer/treatments/*.md` — 주요 암종 10개 초기 가이드라인
- `data/knowledge_base/rare-disease/diagnoses/*.md` — 주요 희귀질환 초기 데이터 (선택)

**수정:**
- `scripts/clinical_board/runner.py` — 실행 전 KB 조회, 실행 후 KB 저장
- `scripts/clinical_board/runner.py` — prior knowledge를 별도 인자로 에이전트에 전달
- `config.yaml` — `knowledge_base.enabled`, `knowledge_base.path` 추가
- `.gitignore` — `data/knowledge_base/kb.sqlite3` 추가 (생성 데이터)

---

## 5. 전체 데이터 흐름

```
입력: VCF + (HPO) + (임상 노트)
  │
  ▼
파이프라인: 분류 엔진 (결정적, 변경 없음)
  │
  ▼
AI Board 실행 전:
  ├── KB에서 prior knowledge 검색 (SQLite 쿼리 → Wiki 로드)
  ├── 로컬 DB에서 에이전트별 도메인 시트 빌드
  └── 공통 브리핑 생성 (케이스 개요 + 임상 노트)
  │
  ▼
AI Board 실행:
  ├── 모드별 에이전트 4명 순차 실행 (temperature 0.1)
  │   각 에이전트: [System Prompt] + [Core Briefing] + [Domain Sheet] + [Prior Knowledge]
  └── Board Chair 종합 (모드별 다른 출력 스키마)
  │
  ▼
AI Board 실행 후:
  ├── SQLite에 변이별 판단 이력 저장
  ├── Wiki 페이지 자동 재생성 (해당 유전자)
  ├── log.md append
  └── 리포트에 Board Opinion 렌더링
```

---

## 6. 구현 순서

| 순서 | 작업 | 의존성 |
|------|------|--------|
| 1 | KB SQLite 스키마 + `build_kb_db.py` (trivial, 병렬 가능) | 없음 |
| 2 | `base.py` 수정 — `domain_sheet` 파라미터 (기본값 `""`), temperature 0.1 | 없음 |
| 3 | `domain_sheets.py` — Rare Disease 도메인 시트 빌더 | 2 |
| 4 | `domain_sheets.py` — Cancer 도메인 시트 빌더 | 2 |
| 5 | Cancer 에이전트 3개 신규 작성 | 2, 4 |
| 6 | `board_chair.py` 모드별 분기 + `CancerBoardOpinion` | 5 |
| 7 | `runner.py` 모드별 에이전트 셋 선택 | 6 |
| 8 | `render.py` — CancerBoardOpinion 렌더링 | 6 |
| 9 | 임상 노트 입력 (orchestrate.py + case_briefing.py) | 2 |
| 10 | KB 저장 (runner.py 수정) | 1, 7 |
| 11 | KB 조회 + domain_sheets.py 통합 | 1, 4 |
| 12 | treatments/ 초기 가이드라인 데이터 구축 | 1 |
| 13 | Wiki 자동 생성 로직 | 1 |
| 14 | 테스트 작성 (전 Phase, ~25-30개 신규) | 전체 |
| 15 | 통합 테스트 + 쇼케이스 리포트 재생성 | 14 |

---

## 7. 테스트 전략

- **도메인 시트 빌더**: mock DB 데이터 → 시트 생성 → 토큰 예산 내 확인
- **Cancer 에이전트**: mock Ollama → JSON 출력 스키마 검증
- **CancerBoardOpinion**: 필드 검증, 렌더링 HTML 확인 (치료 중심 레이아웃 포함)
- **임상 노트**: 길이 트리밍, 한글/영문 입력, 파일 입력
- **지식 베이스**: SQLite CRUD, variant_stats 뷰, Wiki 생성, prior knowledge 검색, empty KB 엣지 케이스 (DB 없음/테이블 없음/행 없음), cross-mode 격리 검증, Wiki 멱등성
- **모드 스위칭**: 모드별 에이전트 셋 정확 인스턴스화 검증
- **하위 호환성**: `domain_sheet=""` 기본값으로 기존 638 테스트 통과 확인, `BoardChair.synthesize()` mode 파라미터 기본값
- **통합**: 실제 MedGemma로 end-to-end 실행 (Ollama 가용 시)
- **예상 신규 테스트**: ~25-30개

### 7.1 Disclaimer 고정 텍스트

모든 Board 출력에 아래 고정 disclaimer를 포함한다 (LLM 생성 아님):

```
[AI-Generated] This opinion was generated by a local multi-specialist AI system
for research assistance only. It is NOT a diagnostic instrument. Final clinical
judgment must be made by the attending clinician. Variant classification (ACMG/AMP)
was performed by a deterministic engine and is not altered by this opinion.
AI Board is powered by Google MedGemma, which is not clinical-grade.
```

---

## 8. Config 변경

```yaml
clinical_board:
  enabled: false
  ollama_url: "http://localhost:11434"
  agent_model: "alibayram/medgemma:27b"
  chair_model: "alibayram/medgemma:27b"
  timeout: 300
  max_retries: 1
  language: "en"
  include_in_report: true
  temperature: 0.1                          # 신규: 기본 temperature
  clinical_note_max_chars: 1500             # 신규: 임상 노트 최대 길이

knowledge_base:
  enabled: true                             # 신규
  path: "data/knowledge_base"               # 신규
  prior_knowledge_max_tokens: 1000          # 신규: prior knowledge 토큰 예산
```

---

## 9. 범위 외 (Not in Scope)

- Vector DB / 임베딩 기반 검색 — 정확 매칭으로 충분
- PubMed API 연동 — 온프레미스 원칙 위반
- 임상 노트에서 HPO 자동 추출 — 정확도 이슈, 분류 엔진 영향
- LLM이 Wiki를 직접 편집 — hallucination 위험
- MSI 분석 — 별도 데이터 소스 필요 (향후)

---

## 10. 구현 현황 (Implementation Status) — 2026-04-14

Phase 1–3이 v2.1 쇼케이스 시점에 모두 통합되었다. 아래는 실제 머지된 커밋 요약이다 (메인 브랜치).

### 10.1 AI Board v2 Phase 1–3 (2026-04-10 ~ 2026-04-13)

| 영역 | 커밋 | 비고 |
|------|------|------|
| Base agent `domain_sheet` 파라미터 + temperature config | 상위 커밋 | 하위 호환 기본값 `""` |
| Rare/Cancer 도메인 시트 빌더 (`domain_sheets.py`) | 상위 커밋 | ClinVar / OncoKB / CIViC / HPO / OMIM / CPIC 주입 |
| Cancer 에이전트 3종 (`therapeutic_target.py`, `tumor_genomics.py`, `clinical_evidence.py`) | `8a0c984` | rare-disease 에이전트와 분리 |
| Board Chair 모드 분기 + `CancerBoardOpinion` | `fe4e90f` | 치료 중심 출력 스키마 |
| Runner 모드 기반 에이전트 셋 선택 | `0346585` | `Optional[Union[BoardOpinion, CancerBoardOpinion]]` |
| 임상 노트 입력 (`--clinical-note` / `--clinical-note-file`) | `02c7f96` | 1500자 캡, KR/EN, mutually exclusive |
| Runner KB 통합 (prior knowledge + 저장) | `820ca35` | SQLite 조회 + Wiki 자동 재생성 |
| Config: temperature, clinical_note_max_chars, knowledge_base | `932638d` | section 8 yaml 그대로 |
| CancerBoardOpinion 렌더링 (치료 중심 레이아웃) | `d3aa5c2` | cancer 템플릿 신규 섹션 |
| KB 스키마 자동 초기화 | `d1ffe2b` | 최초 사용 시 lazy init |
| Cancer orchestrate 로그 `therapeutic_implications` | `2210de9` | rare-disease 전용 로그와 분기 |
| Domain sheet 입력 보호 + KB gene-less 행 skip | `f17a37c` | 문자열 HPO/phenotype 허용 |

### 10.2 Variant Selector v1 (2026-04-14)

임상 기준 검토(`_workspace/variant-selector/00_clinical_review.md`, AMP 2017 PMID 27993330 / ACMG/AMP 2015 PMID 25741868 / ACMG SF v3.2 PMID 37347242 / Schwaederle 2015 / Harris 2016 / Van Allen 2014 근거)를 선행한 후 구현.

| 영역 | 커밋 | 비고 |
|------|------|------|
| `variant_selector.py` 모듈 + 단위 테스트 | `75ea590` | AMP 2017 / ACMG 2015 tiered filter |
| case_briefing + domain_sheets 통합 | `0e52578` | selector output을 브리핑·도메인 시트 입력으로 사용 |
| Runner 배선 + 선택 메타데이터 전파 | `af43c46` | `selection_reason`, `truncated`, `n_dropped`, `empty_reason` |
| `config.yaml` selector caps + criteria overrides | `3cbf518` | soft caps 30 / 20, VUS MAY 10 |
| `render.py` 선택 메타데이터 푸터 | `2fb7a3a` | "Pre-analytic filtering: N → M / Criteria: …" |
| No-findings placeholder (AgentOpinion 0건) | `e318c86` | EN / KO 고정 문구 |
| `scripts/rerender_report.py` | `f46fcbd` | 캐시 JSON → HTML 리렌더 |
| orchestrate `clinical_board` asdict dump | `f595159` | rerender roundtrip 가능 |
| v2.1 쇼케이스 재생성 (cancer / codegen / rare-disease) | `c00e54e`, `5f54708` | codegen: 777 raw → 2 board-presented (KRAS G12D LP + TP53 Tier II VUS) |

### 10.3 임상 검토에 반영된 설계 결정

clinical review에서 요구한 5가지 수정은 전부 반영되었다:

1. **Tier III hotspot/OncoKB MUST-include** — 구현됨. `selection_reason = "Tier_III_hotspot"`로 태그.
2. **"Protein-impacting" VUS 거부** — 구현됨. VUS 인정 기준은 hotspot / OncoKB oncogenic / CGC Tier 1 TSG LoF / indel hotspot으로 엄격화.
3. **HPO threshold + ACMG SF v3.2 제외** — 구현됨. 임상 기준은 operational threshold를 사용하며, SF v3.2 유전자는 opt-in 경로가 모델링되기 전까지 자동 포함 차단.
4. **Soft caps** — 구현됨. MUST-include는 절대 잘리지 않고, MAY-include overflow만 truncate. 메타데이터에 `truncated`/`n_dropped`.
5. **Top-3 fallback 제거** — 구현됨. 조건 불만족 시 `empty_reason` 기반의 honest empty selection. TMB·QC·방법론은 정상 출력.

### 10.4 쇼케이스 산출물

- `docs/showcase/sample_cancer_report_v2.html` — 10-variant demo (NSCLC stage III)
- `docs/showcase/sample_codegen_777_report_v2.html` — 777-variant WGS demo (pancreatic cancer stage IV), selector의 pre-analytic filtering 효과 시연
- `docs/showcase/sample_rare_disease_report_v2.html` — 5-variant demo (pediatric GDD)

세 HTML 모두 AI Board 섹션에 pre-analytic filtering 캡션과 (해당 시) no-findings placeholder가 렌더링된다.

### 10.5 v3 후보 (본 스펙 범위 외)

- ACMG SF v3.2 opt-in 경로 모델링 — 현재는 selector에서 silent exclude
- 동적 tumor-type × OncoKB level 매핑 (현재 Tier I/II는 global rule)
- variant selector UI/API 수준 노출 — 현재는 내부 파이프라인 전용
- MSI / fusion / clonal fraction 기반 추가 필터 — 데이터 소스 확보 후
