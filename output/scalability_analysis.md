# GenomeBoard Scalability Analysis

**Date:** 2026-03-24
**Purpose:** Throughput and scalability analysis for mass report generation (1K-10K+ samples)

---

## Current Bottleneck Analysis

| Pipeline Step | Time/variant | Type | Severity |
|--------------|-------------|------|----------|
| VCF parsing | 0.01ms | CPU | Negligible |
| **ClinVar API** | **400-1000ms** | **Network + rate limit** | **#1 bottleneck** |
| **gnomAD API** | **300-800ms** | **Network** | **#2 bottleneck** |
| KRGDB lookup | <0.001ms | Memory cache | Negligible |
| Freq comparison | <0.01ms | CPU | Negligible |
| PGx check | <0.01ms | Memory cache | Negligible |
| ACMG classification | <0.1ms | CPU | Negligible |
| **PDF generation** | **2-10s/report** | **CPU** | **#3 bottleneck** |

### Root Cause

GenomeBoard was designed as a single-sample clinical tool. The architecture at `orchestrate.py:174` processes variants sequentially with per-request HTTP connections (no session reuse), no response caching, and hard NCBI rate limits (3 req/s without key, 10/s with key).

---

## Throughput Projections

### Current Architecture (Sequential, No Caching, No API Key)

| Scenario | Variants | API calls | Rate-limited time | Report gen | **Total** |
|----------|----------|-----------|-------------------|------------|-----------|
| 1 sample, 100 variants | 100 | ~350 | ~117s (3/s) | 5s (PDF) | **~2 min** |
| 100 samples × 100 variants | 10,000 | ~35,000 | ~3.2hr (3/s) | 500s | **~3.4 hr** |
| 1,000 samples × 100 variants | 100,000 | ~350,000 | ~32hr (3/s) | 5,000s | **~34 hr** |
| 10,000 samples × 50 variants | 500,000 | ~1,750,000 | ~162hr (3/s) | 50,000s | **~7.3 days** |

### Optimized Architecture (Local DBs, Caching, Parallelism)

| Scenario | Unique variants | Cache hit | Local DB lookup | Report gen (8 workers) | **Total** |
|----------|----------------|-----------|-----------------|----------------------|-----------|
| 1 sample, 100 variants | 100 | 0% | ~0.1s | 5s | **~5s** |
| 100 samples × 100 variants | ~6,000 | 40% | ~0.6s | 63s | **~64s** |
| 1,000 samples × 100 variants | ~30,000 | 70% | ~3s | 625s | **~10 min** |
| 10,000 samples × 50 variants | ~80,000 | 84% | ~8s | 6,250s | **~105 min** |

---

## Optimization Priorities

### Priority 1: Local Database Materialization (HIGH impact, MEDIUM effort)

Replace live API calls with local database lookups.

- **ClinVar**: Download weekly XML dump (`variant_summary.txt.gz`, ~100MB). Load into SQLite indexed by `(chrom, pos, ref, alt)` and `rsid`. Lookup drops from 600ms to <1ms.
- **gnomAD**: Download sites-only VCF (~70GB genomes, ~20GB exomes). Pre-process into SQLite or tabix-indexed VCF.
- Implementation: `scripts/db/build_clinvar_db.py`, `scripts/db/build_gnomad_db.py`, config toggle `annotation.source: "local"` vs `"api"`

### Priority 2: Response Cache Layer (HIGH impact, LOW effort)

SQLite variant cache (`data/cache/variant_cache.sqlite3`).

- Key: `(chrom, pos, ref, alt)`, Value: JSON blob of ClinVar + gnomAD response
- TTL: 7 days (configurable)
- 60-85% cache hit rate for clinical panels
- Check cache before API calls

### Priority 3: Batch Orchestration (HIGH impact, MEDIUM effort)

`run_batch_pipeline()` function accepting directory of VCFs or manifest CSV.

- Phase 1: Parse all VCFs, deduplicate variants across samples
- Phase 2: Bulk annotate unique variants only (local DB or cached API)
- Phase 3: Per-sample assembly — look up pre-computed annotations, generate report
- Converts O(samples × variants) to O(unique_variants)

### Priority 4: Connection Pooling (MEDIUM impact, LOW effort)

Replace `requests.get()` with shared `requests.Session()`.

- Reuses TCP connections (HTTP keep-alive)
- Saves ~50-100ms per request in TLS overhead
- ~10-15% speedup on API calls

### Priority 5: Cross-Variant Parallelism (MEDIUM impact, MEDIUM effort)

Replace sequential variant loop with ThreadPoolExecutor + rate limiter.

- Semaphore-based rate limiting (10 concurrent for NCBI with API key)
- 5-10x speedup within rate limits
- For 100 variants: ~2min → ~20-30s

### Priority 6: Multiprocess PDF Generation (MEDIUM impact, LOW effort)

`ProcessPoolExecutor` for WeasyPrint.

- 8 workers: 1000 PDFs from ~83min to ~10min
- Each worker gets pre-computed report_data dict

### Priority 7: ClinVar Batch API (LOW-MEDIUM impact, LOW effort)

ClinVar esummary supports up to 500 UIDs per request.

- Batch esearch, then batch-fetch summaries
- Reduces HTTP requests from ~2N to ~N/250 + N/500

### Priority 8: Queue-Based Architecture (HIGH impact, HIGH effort)

For 10K+ samples/day production:

- Celery + Redis/RabbitMQ task queue
- Kubernetes Jobs for horizontal scaling
- Progress monitoring, retry logic, fault isolation

---

## Implementation Roadmap

### Phase A (1-2 days): Quick Wins
- NCBI API key in config (10x rate limit)
- `requests.Session()` connection pooling
- SQLite variant cache

### Phase B (3-5 days): Batch Mode
- `run_batch_pipeline()` with variant deduplication
- `--batch` CLI flag
- ThreadPoolExecutor + rate limiter for variants
- ProcessPoolExecutor for PDF generation

### Phase C (1-2 weeks): Local Databases
- ClinVar local DB ETL
- gnomAD local DB ETL
- Config toggle: `annotation.source: "local"` vs `"api"`
- Production KRGDB (5M+ lines)

### Phase D (2-4 weeks): Production Infrastructure
- Celery task queue
- Progress monitoring
- Scheduled DB refresh jobs
- Kubernetes deployment

---

## Trade-offs

| Option | Pros | Cons |
|--------|------|------|
| Local DB | 1000x faster; works offline | 70GB+ storage; weekly ETL updates |
| SQLite cache | Simple; huge win for panels | Cache staleness; first-seen variants still hit API |
| Batch orchestration | Eliminates redundant work | More complex code; partial failure handling |
| Connection pooling | Easy 10-15% speedup | Marginal if moving to local DBs |
| Async/concurrent | 5-10x speedup within rate limits | Still rate-limited; complexity |
| Multiprocess PDF | Linear speedup with cores | Memory overhead per worker |
| Celery queue | Production-grade scaling | Operational complexity; overkill for <1000 samples |

---

## Code References

- `orchestrate.py:174` — Sequential variant loop (primary bottleneck)
- `orchestrate.py:88` — ThreadPoolExecutor for per-variant parallelism
- `api_utils.py:18` — `requests.get()` without session reuse
- `query_clinvar.py:26-52` — Up to 4 HTTP requests per variant
- `query_gnomad.py:72-85` — Sequential dataset fallback
- `generate_pdf.py:106` — WeasyPrint single-threaded
- `config.yaml:20` — `ncbi_api_key: ""` (3 req/s limit)
