# GenomeBoard 시스템 아키텍처

## 시스템 개요

GenomeBoard는 **Paperclip-native 멀티 에이전트 컨실리움(consilium)**이다. 각 에이전트는 특정 전문 분야를 담당하며, LLM은 오케스트레이션과 리포트 작성에만 사용된다. 모든 DB 쿼리와 ACMG 분류는 결정론적(deterministic) Python 스크립트로 처리한다.

---

## 핵심 설계 원칙

> **LLM은 오케스트레이션과 리포트 작성에만 사용한다. 모든 DB 쿼리와 ACMG 분류는 결정론적 Python 스크립트가 담당한다.**

이 원칙은 분석 재현성과 신뢰성을 보장한다. LLM이 변이 데이터를 직접 "판단"하거나 "추측"하지 않는다.

---

## 조직도 (Org Chart)

```
Board (User)
     │
    CEO  ──────────────── Counselor
     │                       │
    CTO                   리포트 생성
     │
     ├── Intake Agent
     ├── Clinical Geneticist
     ├── Korean Pop Geneticist
     └── Pharmacogenomicist (Pharma)
```

### 에이전트 어댑터

| 에이전트 | 어댑터 | 모델 |
|---------|--------|------|
| CEO | claude-local | claude-sonnet-4 |
| CTO | claude-local | claude-sonnet-4 |
| Intake | process | — (Python 스크립트) |
| Clinical | process | — (Python 스크립트) |
| Korean Pop | process | — (Python 스크립트) |
| Pharma | codex-local | OpenAI Codex |
| Counselor | claude-local | claude-sonnet-4 |

---

## 데이터 플로우

```
VCF 파일 / 텍스트 입력
        │
        ▼
  [Intake Agent]
  - VCF/텍스트 파싱
  - 변이 정규화
  - 변이별 티켓 생성
        │
        ▼
      [CTO]
  - 에이전트별 태스크 배정
        │
   ┌────┴────┬──────────────┐
   ▼         ▼              ▼
[Clinical]  [Korean Pop]  [Pharma]
 ClinVar     KRGDB         PharmGKB
 조회         gnomAD EAS    CPIC
 ACMG codes  gnomAD ALL    PGx 등급
   │         ACMG codes      │
   └────┬────┘               │
        ▼                    │
      [CTO]  ◄───────────────┘
  ACMG 코드 수집
        │
        ▼
  [ACMG Engine]
  결정론적 분류
  (acmg_engine.py)
        │
        ▼
    [Counselor]
  결과 종합 + 한국인 특이 소견
        │
        ▼
    PDF 리포트
```

---

## 에이전트 런타임: Heartbeat 기반 폴링

Paperclip은 **heartbeat 기반 폴링** 방식으로 에이전트를 실행한다.

- CEO: `heartbeat_interval_seconds: 300` (5분)
- CTO: `heartbeat_interval_seconds: 120` (2분)
- 각 에이전트는 idle 상태에서 heartbeat마다 메시지 큐를 확인한다
- `idle_timeout_seconds: 1800` — 30분 무활동 시 에이전트 종료

> **중요**: LLM 에이전트 (CEO, CTO, Counselor)는 매 heartbeat마다 `/paperclip` 스킬을 invoke해야 한다. 이 스킬 없이는 Paperclip 메시지 큐를 폴링하거나 다른 에이전트에게 메시지를 전송할 수 없다.

이 구조는 항상 켜진 서버 없이도 비동기 멀티 에이전트 협업을 가능하게 한다.

---

## 파일 구조

```
genomeboard/
├── paperclip-config/          # Paperclip 에이전트 설정
│   ├── paperclip.manifest.json  # 회사 목표, 예산, 에이전트 계층 구조
│   └── agents/
│       ├── ceo.json
│       ├── cto.json
│       ├── intake.json
│       ├── clinical.json
│       ├── korean-pop.json
│       ├── pharma.json
│       └── counselor.json
│
├── scripts/                   # 결정론적 Python 스크립트
│   ├── common/
│   │   ├── models.py          # 공유 데이터 모델 (Variant, AcmgEvidence, ...)
│   │   ├── api_utils.py       # HTTP 유틸리티, rate limiting
│   │   └── paperclip_client.py  # Paperclip API 클라이언트
│   ├── intake/
│   │   ├── parse_vcf.py       # VCF 파싱 (cyvcf2)
│   │   └── parse_text.py      # 텍스트 변이 파싱
│   ├── clinical/
│   │   └── query_clinvar.py   # ClinVar E-utilities
│   ├── korean_pop/
│   │   ├── query_krgdb.py     # KRGDB 로컬 TSV 조회
│   │   ├── query_gnomad.py    # gnomAD GraphQL API
│   │   └── compare_freq.py    # 3단계 빈도 비교 + ACMG 코드 산출
│   ├── pharma/
│   │   ├── query_pharmgkb.py  # PharmGKB REST API
│   │   └── korean_pgx.py      # 한국인 PGx 5대 유전자 패턴 매칭
│   ├── classification/
│   │   └── acmg_engine.py     # 결정론적 ACMG 분류 엔진
│   └── counselor/
│       └── generate_pdf.py    # WeasyPrint PDF 생성
│
├── templates/                 # Jinja2 PDF 템플릿
├── data/                      # KRGDB TSV 등 로컬 데이터
├── tests/                     # pytest 테스트 스위트
├── docs/                      # 문서
├── requirements.txt
└── package.json
```

---

## ACMG 분류 엔진

`scripts/classification/acmg_engine.py`는 순수 결정론적 로직으로 동작한다.

- 입력: 각 에이전트가 산출한 `AcmgEvidence` 코드 목록
- 처리: ACMG/AMP 2015 가이드라인 규칙 적용
- 출력: `Pathogenic / Likely Pathogenic / VUS / Likely Benign / Benign`

LLM은 이 분류 결과를 사후에 해석하고 설명하는 역할만 한다.
