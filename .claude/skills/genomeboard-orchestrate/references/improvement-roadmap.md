# BIKO GenomeBoard Improvement Roadmap

**기준일: 2026-04-09** | **현재 상태: v1.0 — 463 tests, main branch**

## 목차

1. [완료된 기능](#완료된-기능)
2. [단기 계획](#단기-계획)
3. [중기 계획](#중기-계획)
4. [장기 계획](#장기-계획)
5. [신규 발견 개선 사항](#신규-발견-개선-사항)
6. [사용자 파이프라인 정보](#사용자-파이프라인-정보)

---

## 완료된 기능

| 기능 | 상태 | 핵심 파일 |
|------|------|----------|
| SNV/Indel 파이프라인 (VCF → ACMG 분류 → 리포트) | ✅ | orchestrate.py, acmg_engine.py |
| Cancer AMP/ASCO/CAP 2017 tiering (Strategy A/B/C) | ✅ | amp_tiering.py |
| CIViC variant-level evidence + match_level 추적 | ✅ | query_civic.py |
| CNV/SV via AnnotSV (ACMG CNV 2020 분류) | ✅ | parse_annotsv.py, models.py |
| TMB 계산 (nonsynonymous variants/Mb) | ✅ | tmb.py |
| Gene Knowledge Phase 2 (CPIC/CIViC/NCBI/GeneReviews/Orphanet/OMIM) | ✅ | build_gene_knowledge.py |
| 로컬 DB: ClinVar, gnomAD(tabix), CIViC, HPO, ClinGen, Orphanet, GeneReviews, OMIM | ✅ | scripts/db/ |
| VUS 필터링 + Tier 순서 정렬 | ✅ | orchestrate.py, templates/ |
| PGx 12개 유전자 + APOE 위험인자 | ✅ | korean_pgx.py |
| 배치 처리 (변이 중복제거 + 병렬 리포트 생성) | ✅ | orchestrate.py |
| Docker 지원 | ✅ | Dockerfile, docker-compose.yml |
| Cancer/Rare-disease 이중 리포트 모드 | ✅ | templates/cancer/, templates/rare-disease/ |
| CNV/SV 리포트 섹션 (Class 4-5 소견, Class 3 DS, 통계) | ✅ | shared/sv_section.html |
| SV VUS dosage sensitivity 필터링 | ✅ | models.py (is_dosage_sensitive) |
| 한국인 빈도 비교 (KRGDB vs gnomAD EAS vs ALL) | ✅ | compare_freq.py |

---

## 단기 계획 (바로 가능, 난이도 낮음)

### 1. OMIM genemap2.txt 통합
- **상태:** 승인 대기
- **내용:** 유전패턴(AD/AR/XL) 11개 정적 dict → genemap2.txt 전체 확장
- **작업:** query_omim.py 수정, build script 확장
- **에이전트:** db-dev + pipeline-dev
- **예상 소요:** 반나절

### 2. Korea4K AF 통합
- **내용:** korea4k.10000genomes.org AF 데이터. KRGDB 패턴과 동일 (TSV → 로컬 조회)
- **작업:** build_korea4k_db.py, query_korea4k.py, compare_freq.py 확장
- **에이전트:** db-dev + pipeline-dev + qa-engineer
- **의존성:** 없음

### 3. NARD2 통합
- **내용:** nard.macrogen.com 빈도 데이터
- **작업:** Korea4K와 동일 패턴
- **에이전트:** db-dev + pipeline-dev + qa-engineer
- **의존성:** Korea4K와 병렬 가능

---

## 중기 계획 (설계 필요, 난이도 중간)

### 4. In silico predictions 통합
- **내용:** REVEL, CADD, AlphaMissense, SpliceAI 점수 표시
- **소스:** VEP 주석에서 파싱 (사용자 파이프라인이 이미 ANNOVAR로 주석 중)
- **CADD:** offline GRCh38 설치 필요 (대용량)
- **작업:** VEP CSQ 필드 파싱 확장, 리포트 표시, PP3/BP4 evidence code 활용
- **에이전트:** pipeline-dev + report-dev + qa-engineer + clinical-advisor

### 5. MSI (Microsatellite Instability)
- **내용:** microsatellite loci DB + 알고리즘 (MSIsensor2/MANTIS 스타일)
- **임상 의미:** 면역치료 바이오마커 (pembrolizumab)
- **작업:** MSI 계산 엔진, microsatellite loci DB 구축, 리포트 섹션
- **에이전트:** 전원 (clinical-advisor 필수)

### 6. Family trio/quartet analysis
- **내용:** 부모 genotype 표시 (ref_read, alt_read, VAF), de novo variant 탐지
- **작업:** multi-sample VCF 파싱, de novo 판정 로직, 리포트 가계도 섹션
- **에이전트:** 전원 (rare disease 핵심 기능)

### 7. InterVar ACMG 자동 분류 통합
- **내용:** InterVar 결과 파싱 또는 자체 ACMG 규칙 강화
- **작업:** InterVar TSV 파싱 또는 acmg_engine.py 확장
- **에이전트:** clinical-advisor + pipeline-dev + qa-engineer

### 8. CNV/SV Phase 2 — ACMG CNV 2020 자체 분류 엔진
- **내용:** AnnotSV 없이 독자 ACMG CNV 분류
- **참조:** Riggs et al. 2020 가이드라인
- **작업:** CNV 분류 엔진, evidence scoring, 리포트 통합
- **에이전트:** 전원

---

## 장기 계획 (큰 작업)

### 9. CNV/SV Phase 3 — Fusion 전용 처리
- **내용:** Manta BND → gene fusion detection (ALK, ROS1, NTRK 등)
- **작업:** fusion caller 결과 파싱, 융합 유전자 DB, 리포트 fusion 섹션
- **에이전트:** 전원

### 10. v2.0 Report Results DB
- **내용:** PostgreSQL/SQLite DB로 리포트 결과 저장/검색/통계
- **스키마:** sample_id, 분석일, 변이 수, 분류 결과 요약, DB 버전
- **기능:** 변이별 분류 이력, 재분류 추적, 검색/필터/통계 API, 대시보드
- **에이전트:** pipeline-dev + db-dev + (프론트엔드)

---

## 신규 발견 개선 사항

코드베이스 분석을 통해 새로 식별된 개선 사항.

### A. orchestrate.py 리팩토링 (권장: 중기)
- **문제:** 1,200줄 모놀리스. run_pipeline, _query_variant_databases, _classify_variants, _build_variant_records 등이 단일 파일에 집중
- **제안:** 기능별 모듈 분리
  - `scripts/pipeline/query.py` — DB 쿼리 오케스트레이션
  - `scripts/pipeline/classify.py` — 분류 + tiering
  - `scripts/pipeline/assemble.py` — 레코드 조립
  - `scripts/orchestrate.py` — CLI + 진입점만 유지
- **에이전트:** pipeline-dev + qa-engineer
- **위험:** 대규모 리팩토링으로 인한 회귀 — 전체 테스트 통과 필수

### B. VCF 변이 정규화 강화 (권장: 중기)
- **문제:** 현재 텍스트 기반 VCF 파싱. left-align, multi-allelic decompose 미지원
- **제안:** 정규화 단계 추가 (pyvcf3 또는 cyvcf2 활용)
- **영향:** ClinVar/gnomAD 매칭 정확도 향상

### C. Config 스키마 검증 (권장: 단기)
- **문제:** config.yaml 유효성 검증 없음. 잘못된 키/타입 시 런타임 에러
- **제안:** pydantic 또는 jsonschema로 로드 시점 검증
- **에이전트:** pipeline-dev

### D. REST API 레이어 (권장: 장기)
- **문제:** CLI 전용. 웹 통합 시 subprocess 호출 필요
- **제안:** FastAPI/Flask 기반 REST API (run_pipeline을 HTTP endpoint로 노출)
- **의존성:** Report DB (v2.0)와 함께 진행

### E. 리포트 국제화 (권장: 장기)
- **문제:** 리포트가 영문 고정. 한국 임상 환경에서 한글 리포트 요구 가능
- **제안:** Jinja2 i18n 확장 또는 템플릿 이중화
- **에이전트:** report-dev

### F. 증분 분석 (권장: 중기)
- **문제:** 재분석 시 전체 변이를 처음부터 다시 쿼리
- **제안:** variant_cache.sqlite3 활용하여 변경된 변이만 재쿼리
- **에이전트:** pipeline-dev

### G. 모니터링/로깅 강화 (권장: 단기)
- **문제:** verbose 모드만 존재. 구조화된 로깅 없음
- **제안:** Python logging 모듈 + 구조화 로그 (JSON format)
- **에이전트:** pipeline-dev

---

## 사용자 파이프라인 정보

실제 사용 환경의 업스트림 파이프라인 구성 (BIKO GenomeBoard 입력 데이터 생성):

- **시퀀싱:** WGS 기반
- **Variant Callers:** SNV/Indel (Mutect2 등) + Canvas (CNV) + Manta (SV)
- **주석:** ANNOVAR + VEP (SNP/Indel), AnnotSV (CNV/SV)
- **주석 DB:** OMIM genemap2.txt, ClinVar, InterVar, REVEL, CADD, AlphaMissense, SpliceAI, gnomAD v4, NARD2, Korea4K AF
- **Gene info:** pLI/O-E/missense Z (gnomAD), ClinGen dosage sensitivity, MANE transcript
- **Family:** trio 기본, quartet/quintet 가능

---

## 우선순위 매트릭스

| 우선순위 | 항목 | 임상 가치 | 구현 난이도 | 의존성 |
|---------|------|----------|-----------|--------|
| 🔴 1 | OMIM genemap2.txt | 높음 | 낮음 | 승인 대기 |
| 🔴 2 | Korea4K AF + NARD2 | 높음 (한국인 특화) | 낮음 | 없음 |
| 🟠 3 | In silico predictions | 높음 | 중간 | VEP 파싱 확장 |
| 🟠 4 | Config 스키마 검증 | 중간 (안정성) | 낮음 | 없음 |
| 🟠 5 | 모니터링/로깅 | 중간 (운영) | 낮음 | 없음 |
| 🟡 6 | Family trio/quartet | 매우 높음 | 높음 | multi-sample VCF |
| 🟡 7 | InterVar ACMG | 높음 | 중간 | InterVar 설치 |
| 🟡 8 | MSI | 높음 (면역치료) | 높음 | loci DB |
| 🟡 9 | orchestrate.py 리팩토링 | 중간 (유지보수) | 중간 | 없음 (테스트 위험) |
| 🟡 10 | VCF 정규화 | 중간 | 중간 | 없음 |
| 🔵 11 | CNV/SV Phase 2 | 높음 | 높음 | ACMG CNV 2020 |
| 🔵 12 | 증분 분석 | 중간 (성능) | 중간 | 캐시 확장 |
| 🔵 13 | CNV/SV Phase 3 (Fusion) | 높음 | 매우 높음 | fusion DB |
| 🔵 14 | REST API | 중간 | 높음 | Report DB |
| 🔵 15 | Report DB (v2.0) | 높음 | 매우 높음 | 없음 |
| 🔵 16 | 리포트 국제화 | 낮음 | 중간 | 없음 |
