# Clinical Board — Multi-Agent Diagnostic Synthesis System

**Date:** 2026-04-13
**Status:** Approved
**Scope:** Local LLM 기반 다전문가 종합 진단 추론 시스템

---

## Overview

BIKO GenomeBoard의 결정적(deterministic) 분류 엔진 위에, Local LLM(MedGemma 27B + Gemma4 31B)을 활용한 다전문가 Clinical Board를 구축한다. 4명의 도메인 전문 Agent가 독립 분석 후, Board Chair가 종합하여 진단 의견서를 생성한다.

**핵심 원칙:**
- 변이 분류(ACMG/AMP)는 기존 결정적 엔진이 수행 — 변경 없음
- LLM은 분류 결과를 **입력으로만** 사용하여 종합 추론 수행
- 모든 LLM 생성 콘텐츠에 `[AI-Generated]` 태그 필수 표시
- 온프레미스 LLM(Ollama) 사용 — 환자 데이터 외부 전송 없음

---

## Architecture

```
[결정적 분류 엔진] → run_pipeline() 결과 (variant_records, classification, HPO 등)
        ↓
[Case Briefing] → 구조화된 케이스 요약 생성
        ↓
[4 Domain Agents — 병렬, MedGemma 27B]
  ├── Variant Pathologist (변이 병리)
  ├── Disease Geneticist (질환 유전학)
  ├── PGx Specialist (약물유전체)
  └── Literature Analyst (문헌 근거)
        ↓
[Board Chair — Gemma4 31B]
  → 4개 소견 종합 → 진단 의견서
        ↓
[Report Integration] → 기존 HTML 리포트에 "Clinical Board Opinion" 섹션 추가
```

---

## Components

### 1. Ollama Client (`scripts/clinical_board/ollama_client.py`)

Ollama REST API 클라이언트:
- `generate(model, prompt, system, format) -> str` — 텍스트 생성
- `generate_json(model, prompt, system, schema) -> dict` — 구조화된 JSON 출력
- Health check, timeout 관리, retry 로직
- 기본 endpoint: `http://localhost:11434`
- config.yaml에서 설정 가능: `clinical_board.ollama_url`, `clinical_board.timeout`

### 2. Case Briefing Builder (`scripts/clinical_board/case_briefing.py`)

run_pipeline() 결과를 LLM이 이해할 수 있는 구조화된 케이스 요약으로 변환:

```python
def build_case_briefing(report_data: dict, mode: str) -> str:
    """Convert pipeline results to a structured case briefing for LLM agents."""
```

포함 내용:
- 환자 정보 (sample_id, 분석 모드)
- 각 변이의 분류 결과, evidence codes, tier, ClinVar 정보
- In silico scores (REVEL, CADD, AlphaMissense, SpliceAI)
- 한국인 빈도 데이터 (KRGDB, Korea4K, NARD2, gnomAD)
- HPO 표현형 매칭 결과 (rare disease)
- OMIM 유전자-질환 매핑, 유전패턴
- PGx 결과
- SV/CNV 결과 (있는 경우)
- TMB 결과 (cancer)
- Gene knowledge 요약

### 3. Domain Agents (`scripts/clinical_board/agents/`)

각 Agent는 동일한 인터페이스를 따름:

```python
@dataclass
class AgentOpinion:
    agent_name: str           # 전문의 이름
    domain: str               # 전문 영역
    findings: list[dict]      # 주요 소견 [{finding, evidence, confidence}]
    recommendations: list[str] # 권고사항
    concerns: list[str]       # 우려/주의사항
    references: list[str]     # PMID 등 참고문헌
    confidence: str           # high/moderate/low

def analyze(case_briefing: str, model: str = "medgemma:27b") -> AgentOpinion
```

#### 3a. Variant Pathologist (`variant_pathologist.py`)
- **역할:** 각 변이의 기능적 영향 분석
- **분석:** 단백질 도메인 위치, 구조적 영향, conservation, 유사 위치 pathogenic 변이
- **핵심 질문:** "이 변이가 단백질 기능을 얼마나 손상시키는가?"
- **VUS 초점:** VUS 변이에 대해 재분류 가능성 평가 (추가 근거 제안)
- **모델:** MedGemma 27B

#### 3b. Disease Geneticist (`disease_geneticist.py`)
- **역할:** 변이-질환 연관성 및 감별진단
- **분석:** 유전패턴(AD/AR/XL) + zygosity 맥락, 표현형-유전형 상관, 감별진단 목록
- **핵심 질문:** "이 변이 조합이 어떤 질환을 시사하는가?"
- **종합 시각:** 개별 변이가 아닌 전체 변이 프로필을 종합적으로 해석
- **모델:** MedGemma 27B

#### 3c. PGx Specialist (`pgx_specialist.py`)
- **역할:** 약물유전체 상호작용 분석
- **분석:** PGx 변이의 약물 영향, 복합 PGx 효과, 한국인 특화 약물 반응
- **핵심 질문:** "이 환자의 약물 처방에서 주의할 점은?"
- **한국인 특화:** 한국인 고빈도 PGx 변이와 서양 가이드라인의 차이
- **모델:** MedGemma 27B
- **조건:** PGx 결과가 있을 때만 활성화

#### 3d. Literature Analyst (`literature_analyst.py`)
- **역할:** 변이별 최신 문헌 근거 종합
- **분석:** 기능연구, 케이스리포트, 메타분석 등 문헌 기반 evidence
- **핵심 질문:** "이 변이에 대해 어떤 임상 근거가 있는가?"
- **RAG 연동:** (향후) 로컬 PubMed abstract 벡터 DB 검색
- **현재:** LLM의 학습 데이터 기반 문헌 지식 활용
- **모델:** MedGemma 27B

### 4. Board Chair (`scripts/clinical_board/board_chair.py`)

4명의 전문의 소견을 종합하는 위원장:

```python
@dataclass
class BoardOpinion:
    primary_diagnosis: str          # 1차 진단
    differential_diagnoses: list    # 감별진단 목록 [{diagnosis, likelihood, evidence}]
    key_findings: list[str]         # 핵심 소견
    recommendations: list[str]     # 최종 권고사항
    agent_consensus: str           # 전문의 합의 수준 (unanimous/majority/split)
    dissenting_opinions: list[str] # 소수 의견 (있을 경우)
    follow_up: list[str]          # 추적 관찰 사항
    confidence: str               # high/moderate/low
    disclaimer: str               # AI 생성 면책 고지

def synthesize(case_briefing: str, agent_opinions: list[AgentOpinion],
               model: str = "gemma4:31b") -> BoardOpinion
```

- **모델:** Gemma4 31B (가장 강력한 추론)
- 의견 불일치 시 양쪽 논거를 명시하고 근거가 더 강한 쪽을 채택
- 모든 권고에 근거 출처 첨부
- 면책 고지 자동 포함

### 5. Board Runner (`scripts/clinical_board/runner.py`)

전체 Board 프로세스를 오케스트레이션:

```python
def run_clinical_board(report_data: dict, mode: str,
                       ollama_url: str = None) -> BoardOpinion:
    """Run the full Clinical Board analysis."""
    # 1. Build case briefing
    # 2. Run 4 domain agents (parallel via ThreadPoolExecutor)
    # 3. Run Board Chair synthesis
    # 4. Return structured opinion
```

### 6. Report Integration (`scripts/clinical_board/render.py`)

BoardOpinion을 HTML로 렌더링:
- 기존 리포트의 Methodology 섹션 앞에 "Clinical Board Opinion" 섹션 삽입
- 각 전문의 소견 요약 (접힘/펼침)
- 최종 진단 + 감별진단 테이블
- 권고사항 목록
- `[AI-Generated — 임상의 검토 필수]` 배너

---

## Pipeline Integration

### CLI
```bash
# Clinical Board 활성화
python scripts/orchestrate.py sample.vcf -o report.html --clinical-board

# 특정 모델 지정
python scripts/orchestrate.py sample.vcf -o report.html --clinical-board \
  --board-model gemma4:31b --agent-model medgemma:27b

# Board만 단독 실행 (기존 리포트 기반)
python scripts/clinical_board/runner.py --report report.json
```

### config.yaml
```yaml
clinical_board:
  enabled: false           # 기본 비활성화
  ollama_url: "http://localhost:11434"
  agent_model: "medgemma:27b"
  chair_model: "gemma4:31b"
  timeout: 120             # 에이전트당 최대 초
  max_retries: 2
  include_in_report: true  # HTML 리포트에 섹션 추가
```

### run_pipeline 통합

```python
# orchestrate.py run_pipeline() 마지막에 추가
if clinical_board:
    from scripts.clinical_board.runner import run_clinical_board
    board_opinion = run_clinical_board(report_data, mode)
    report_data["clinical_board"] = board_opinion
```

---

## File Structure

```
scripts/clinical_board/
├── __init__.py
├── ollama_client.py          — Ollama REST API 클라이언트
├── case_briefing.py          — 케이스 요약 빌더
├── models.py                 — AgentOpinion, BoardOpinion 데이터 모델
├── runner.py                 — Board 실행 오케스트레이터
├── render.py                 — HTML 렌더링
└── agents/
    ├── __init__.py
    ├── base.py               — BaseAgent 추상 클래스
    ├── variant_pathologist.py
    ├── disease_geneticist.py
    ├── pgx_specialist.py
    └── literature_analyst.py
    └── board_chair.py

templates/shared/
└── clinical_board_section.html  — Board Opinion HTML 템플릿

tests/
├── test_ollama_client.py
├── test_case_briefing.py
├── test_clinical_board_agents.py
├── test_board_chair.py
└── test_clinical_board_integration.py
```

---

## Testing Strategy

- **Unit tests:** mock Ollama API, 각 Agent의 prompt + output 파싱 검증
- **Integration tests:** 실제 Ollama 연동 테스트 (`@pytest.mark.ollama` 마커, CI에서 제외)
- **Golden test:** 샘플 VCF → 전체 Board 실행 → 출력 구조 검증
- 기존 599 tests는 영향 없음 (clinical_board는 독립 모듈)

---

## Implementation Phases

### Phase 1: Foundation
- Ollama client + health check
- Data models (AgentOpinion, BoardOpinion)
- Case briefing builder
- config.yaml 설정 추가

### Phase 2: Domain Agents
- BaseAgent 추상 클래스
- 4개 domain agent 구현 (system prompt + output parsing)
- 각 agent 단독 테스트

### Phase 3: Board Chair + Integration
- Board Chair synthesis logic
- Board runner orchestrator
- HTML 렌더링 + 템플릿
- orchestrate.py / CLI 통합

### Phase 4: Testing + Polish
- 전체 통합 테스트
- 샘플 리포트 생성 + 품질 검증
- 프롬프트 튜닝
