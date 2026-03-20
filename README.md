# GenomeBoard

**한국인 특화 유전체 변이 해석 AI 컨실리움** — Paperclip 멀티 에이전트 기반으로 VCF 입력부터 ACMG 분류, 약물유전체 스크리닝, 한국인 집단 빈도 분석까지 자동화된 PDF 리포트를 생성한다.

> **Research Use Only** — 이 소프트웨어는 연구 목적으로만 사용한다. 임상 진단이나 의료적 의사결정에 사용해서는 안 된다.

---

## 주요 기능

- **ACMG/AMP 2015 분류** — 결정론적 규칙 엔진으로 Pathogenic/VUS/Benign 판정
- **한국인 집단 빈도 분석** — KRGDB → gnomAD EAS → gnomAD ALL 3단계 비교
- **약물유전체(PGx) 스크리닝** — CYP2D6, CYP2C19, CYP2C9, HLA-B, NUDT15 한국인 특화 패턴
- **자동 PDF 리포트** — Genetic Counselor 에이전트가 한국어로 최종 리포트 생성

---

## 빠른 시작

```bash
# 1. 의존성 설치
pip install -r requirements.txt && pnpm install

# 2. 환경 변수 설정
cp .env.example .env  # ANTHROPIC_API_KEY, OPENAI_API_KEY 입력

# 3. Paperclip 실행
pnpm dev
```

Paperclip Board가 열리면 CEO에게 메시지를 보내 분석을 시작한다.

자세한 설치 방법은 [docs/SETUP.md](docs/SETUP.md)를 참조.

---

## 아키텍처 개요

GenomeBoard는 Paperclip 멀티 에이전트 프레임워크 위에서 동작한다. **LLM은 오케스트레이션과 리포트 작성만 담당하고, 모든 DB 쿼리와 ACMG 분류는 결정론적 Python 스크립트가 처리**한다.

```
Board (User) → CEO → CTO → [Intake | Clinical | Korean Pop | Pharma]
                              ↓ ACMG codes 수집
                         ACMG Engine (deterministic)
                              ↓
                         Counselor → PDF 리포트
```

상세 아키텍처는 [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)를 참조.

---

## 에이전트 설명

| 에이전트 | 역할 | 어댑터 |
|---------|------|--------|
| **CEO** | 사용자 소통, 분석 요청 접수, 예산 관리 | Claude |
| **CTO** | 분석 오케스트레이션, ACMG 코드 수집 및 판정 | Claude |
| **Intake** | VCF/텍스트 파싱, 변이 정규화 | Bash (Python) |
| **Clinical Geneticist** | ClinVar 조회, ACMG evidence 코드 산출 | Bash (Python) |
| **Korean Pop Geneticist** | KRGDB + gnomAD 빈도 비교, 한국인 특이 변이 플래그 | Bash (Python) |
| **Pharmacogenomicist** | PharmGKB/CPIC 조회, 한국인 PGx 5대 유전자 패턴 | Codex |
| **Genetic Counselor** | 전체 결과 종합, 한국어 PDF 리포트 생성 | Claude |

---

## 기술 스택

| 구분 | 기술 |
|------|------|
| 에이전트 프레임워크 | [Paperclip](https://paperclipai.com) |
| LLM | Claude Sonnet 4 (Anthropic), Codex (OpenAI) |
| VCF 파싱 | cyvcf2 |
| 변이 DB | ClinVar (E-utilities), gnomAD (GraphQL), KRGDB (로컬 TSV) |
| PGx DB | PharmGKB REST API, CPIC 가이드라인 |
| PDF 생성 | WeasyPrint + Jinja2 |
| 테스트 | pytest |
| 패키지 관리 | pip (Python), pnpm (Node.js) |

---

## 문서

| 문서 | 내용 |
|------|------|
| [docs/SETUP.md](docs/SETUP.md) | 설치 및 실행 가이드 |
| [docs/API_KEYS.md](docs/API_KEYS.md) | API 키 설정 레퍼런스 |
| [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) | 시스템 아키텍처 |
| [docs/KOREAN_STRATEGY.md](docs/KOREAN_STRATEGY.md) | 한국인 특화 분석 전략 |

---

## 라이선스

MIT License — 자세한 내용은 [LICENSE](LICENSE) 파일 참조.

---

> **Research Use Only**: GenomeBoard는 연구 및 교육 목적으로 개발된 소프트웨어다. 임상 진단, 치료 결정, 또는 의료적 의사결정에 사용해서는 안 된다. 모든 변이 해석은 반드시 자격을 갖춘 임상 유전학자 또는 의료 전문가의 검토를 거쳐야 한다.
