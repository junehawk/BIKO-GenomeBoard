---
name: genomeboard-orchestrate
description: "GenomeBoard 기능 개발 오케스트레이터. 새 기능 추가, 데이터 소스 통합, 파이프라인 수정, 리포트 개선, 리팩토링 등 모든 GenomeBoard 개발 작업을 에이전트 팀으로 조율. '구현', '통합', '추가', '개선', '개발', '만들어', '수정', 'implement', 'integrate', 'add', 'improve', 'develop', 'fix', 'refactor' 등 개발 작업 요청 시 반드시 이 스킬을 사용. GenomeBoard 프로젝트의 개발 계획을 수립하거나 현재 상태를 파악할 때도 사용."
---

# GenomeBoard Orchestrator

GenomeBoard 프로젝트의 모든 개발 작업을 에이전트 팀으로 조율하는 오케스트레이터.

## 실행 모드: Agent Team

## 에이전트 구성

| 팀원 | 에이전트 타입 | 역할 | 스킬 | 산출물 |
|------|-------------|------|------|--------|
| clinical-advisor | general-purpose | 임상 자문 | clinical-validation | 임상 의견서 |
| pipeline-dev | general-purpose | 파이프라인 개발 | genomeboard-conventions | Python 모듈 |
| db-dev | general-purpose | DB 통합 | genomeboard-conventions | build/query 모듈 |
| report-dev | general-purpose | 리포트 개발 | genomeboard-conventions | HTML 템플릿 |
| qa-engineer | general-purpose | QA/테스트 | genomeboard-conventions | 테스트 코드 + 검증 리포트 |

**모든 Agent 호출에 `model: "opus"` 필수.**

## 작업 유형별 팀 구성

모든 작업에 5명 전원이 필요하지는 않다. 작업 유형에 따라 최적 부분집합을 선택한다.

| 작업 유형 | 필수 에이전트 | 선택 에이전트 |
|----------|-------------|-------------|
| 데이터 소스 통합 (Korea4K, NARD2 등) | db-dev, pipeline-dev, qa-engineer | report-dev |
| 분류 로직 변경 (InterVar, in silico 등) | clinical-advisor, pipeline-dev, qa-engineer | — |
| 리포트 개선 (새 섹션, 레이아웃) | report-dev, qa-engineer | clinical-advisor |
| 전체 기능 (MSI, Trio, Fusion 등) | 전원 | — |
| 리팩토링/성능 | pipeline-dev, qa-engineer | — |
| 계획 수립/분석 | clinical-advisor | pipeline-dev |

## 워크플로우

### Phase 1: 준비 (리더 단독)

1. 사용자 요청을 분석하여 작업 유형 판별
2. `_workspace/` 폴더 생성 (이미 있으면 재사용)
3. 필요한 에이전트 부분집합 결정
4. 작업 분해: 각 에이전트의 구체적 태스크 정의

```bash
mkdir -p _workspace/
```

### Phase 2: 임상 사전 검토 (해당 시)

분류/tiering/임상적 결정이 필요한 작업인 경우:

1. clinical-advisor에게 사전 검토 요청 (서브에이전트로 실행)
2. 검토 결과를 `_workspace/00_clinical_review.md`에 저장
3. 이 결과를 이후 팀원들의 프롬프트에 포함

임상 검토가 불필요한 작업(리팩토링, UI 개선 등)은 이 Phase를 건너뛴다.

### Phase 3: 팀 구성 및 작업 배분

```
TeamCreate(
  team_name: "genomeboard-team",
  members: [선택된 에이전트들]
)

TaskCreate(tasks: [
  { title, description, assignee, depends_on }
])
```

**태스크 설계 원칙:**
- 에이전트당 1-3개 태스크 (명확한 산출물 기준)
- db-dev 태스크가 pipeline-dev보다 먼저 완료되도록 의존성 설정
- qa-engineer 태스크는 구현 태스크에 의존
- 각 태스크 description에 genomeboard-conventions 스킬 참조 지시 포함

### Phase 4: 구현 (팀 자율 조율)

팀원들이 자체 조율하며 작업 수행. 리더는:
- 진행 상황 모니터링 (TaskGet)
- 교착 상태 감지 시 개입 (SendMessage)
- 임상 자문 결과 공유 필요 시 중재

**파일 기반 데이터 전달 규칙:**
- 중간 산출물: `_workspace/{agent}_{artifact}.md`
- 최종 산출물: 프로젝트 소스 디렉토리에 직접 작성

### Phase 5: 검증

1. qa-engineer가 테스트 작성 + 실행
2. 전체 테스트 스위트 실행: `python -m pytest tests/ -x -q`
3. 테스트 리포트 기반으로 문제 발견 시 해당 에이전트에게 수정 요청
4. clinical-advisor의 최종 임상 검토 (해당 시)

### Phase 6: 정리

1. 모든 태스크 완료 확인
2. 팀 해체 (TeamDelete)
3. `_workspace/` 보존 (삭제 금지 — 감사 추적용)
4. 사용자에게 결과 요약 보고

## 데이터 흐름

```
[리더] ──TeamCreate──→ [팀원들]
                         │
  clinical-advisor ──→ _workspace/00_clinical_review.md
  db-dev ──→ scripts/db/build_*.py, query_*.py
  pipeline-dev ──→ scripts/**/*.py, orchestrate.py
  report-dev ──→ templates/**/*.html
  qa-engineer ──→ tests/test_*.py
                         │
[리더] ←──TaskGet/결과수집──┘
```

## 에러 핸들링

| 상황 | 대응 |
|------|------|
| 팀원 1명 실패 | 상태 확인 → 1회 재시도 → 실패 시 해당 결과 없이 진행, 리포트에 누락 명시 |
| 과반수 실패 | 사용자에게 알림, 진행 여부 확인 |
| 테스트 실패 | 실패 테스트 분석 → 해당 에이전트에게 수정 요청 → 재테스트 |
| 데이터 충돌 | 출처 병기, 삭제하지 않고 보존 |
| 임상 의견 불일치 | clinical-advisor 의견 우선, 사유 기록 |

## 테스트 시나리오

### 정상 흐름: Korea4K AF 데이터 소스 통합

1. 리더가 작업 분석 → db-dev + pipeline-dev + qa-engineer 선정
2. db-dev: build_korea4k_db.py + query_korea4k.py 작성
3. pipeline-dev: orchestrate.py에 Korea4K 조회 추가, compare_freq.py 확장
4. qa-engineer: test_korea4k.py 작성 + 전체 테스트 실행
5. 전체 통과 → 결과 보고

### 에러 흐름: DB 빌드 실패

1. db-dev가 Korea4K 데이터 형식 오류로 빌드 실패
2. db-dev가 리더에게 실패 보고 + 데이터 형식 분석 결과 공유
3. 리더가 사용자에게 데이터 형식 확인 요청
4. 수정된 데이터로 재빌드 → 계속 진행

## 개선 로드맵 참조

현재 프로젝트 상태와 향후 계획은 `references/improvement-roadmap.md`를 참조한다. 계획 수립이나 우선순위 결정 시 이 문서를 Read하여 참고한다.
