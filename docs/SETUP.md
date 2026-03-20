# GenomeBoard 설치 가이드

## 사전 요구사항 (Prerequisites)

| 도구 | 최소 버전 | 비고 |
|------|-----------|------|
| Python | 3.10+ | 가상환경 권장 |
| Node.js | 18+ | LTS 버전 권장 |
| pnpm | 8+ | `npm install -g pnpm` |
| bcftools | 1.17+ | VCF 전처리용 (선택) |

---

## 1. 저장소 클론

```bash
git clone https://github.com/your-org/genomeboard.git
cd genomeboard
```

---

## 2. 환경 변수 설정

```bash
cp .env.example .env
```

`.env` 파일을 열어 필요한 API 키를 입력한다. 각 키에 대한 설명은 [docs/API_KEYS.md](./API_KEYS.md)를 참조.

---

## 3. Python 환경 설정

```bash
# 가상환경 생성 (권장)
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate

# 의존성 설치
pip install -r requirements.txt
```

주요 Python 의존성:
- `cyvcf2` — VCF 파일 파싱
- `requests` — 외부 DB API 호출 (ClinVar, gnomAD, PharmGKB)
- `WeasyPrint` — PDF 리포트 생성
- `Jinja2` — PDF 템플릿 렌더링

---

## 4. Node.js 환경 설정 (Paperclip)

```bash
pnpm install
```

---

## 5. VCF 파일 전처리

GenomeBoard는 PASS 필터 변이만 처리한다. 분석 전 bcftools로 필터링을 권장한다.

### 기본 필터링

```bash
# PASS 변이 + Gene 어노테이션이 있는 변이만 추출
bcftools view -f PASS -i 'INFO/Gene!="."' input.vcf > filtered.vcf
```

### 단일 유전자 필터링

```bash
bcftools view -f PASS -i 'INFO/Gene="BRCA1"' input.vcf > brca1.vcf
```

### 다중 유전자 필터링

```bash
bcftools view -f PASS -i 'INFO/Gene="CYP2D6" || INFO/Gene="CYP2C19"' input.vcf > pgx.vcf
```

### BGZ 압축 및 인덱싱

```bash
bgzip filtered.vcf
tabix -p vcf filtered.vcf.gz
```

> **참고**: VCF는 GRCh38 기준을 권장한다. KRGDB 데이터가 GRCh38 좌표로 구축되어 있다.

---

## 6. 테스트 실행

```bash
# 전체 테스트 (verbose)
pytest tests/ -v

# 커버리지 포함
pytest tests/ --cov=scripts --cov-report=html

# 특정 모듈만
pytest tests/test_acmg_logic.py -v
pytest tests/test_korean_pgx.py -v
```

---

## 7. Paperclip 개발 서버 실행

```bash
pnpm dev
```

브라우저에서 Paperclip Board UI가 열리면 CEO 에이전트에게 메시지를 보내 분석을 시작할 수 있다.

---

## 8. 분석 요청 예시

Paperclip Board에서 CEO에게 다음과 같이 요청한다:

```
BRCA1 chr17:43057051 A>G 변이를 분석해주세요.
```

또는 VCF 파일 경로를 전달한다:

```
data/sample.vcf 파일의 변이를 분석해주세요.
```

---

## 문제 해결

**WeasyPrint 설치 실패 (macOS)**
```bash
brew install pango cairo gdk-pixbuf libffi
pip install WeasyPrint
```

**cyvcf2 설치 실패**
```bash
brew install htslib
pip install cyvcf2
```
