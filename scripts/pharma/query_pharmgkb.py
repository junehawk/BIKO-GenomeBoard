# scripts/pharma/query_pharmgkb.py
import time
from typing import Optional, Dict
from scripts.common.models import Variant
from scripts.common.api_utils import fetch_with_retry

PHARMGKB_API = "https://api.pharmgkb.org/v1/data"
_last_request_time = 0

def _rate_limit():
    global _last_request_time
    elapsed = time.time() - _last_request_time
    if elapsed < 0.5:
        time.sleep(0.5 - elapsed)
    _last_request_time = time.time()

def _fetch_pharmgkb(gene: str) -> Optional[dict]:
    _rate_limit()
    url = f"{PHARMGKB_API}/clinicalAnnotation"
    params = {"gene": gene, "view": "base"}
    return fetch_with_retry(url, params=params)

def query_pharmgkb(variant: Variant) -> Optional[Dict]:
    if not variant.gene:
        return None
    data = _fetch_pharmgkb(variant.gene)
    if not data:
        return None
    return {
        "agent": "pharmacogenomicist",
        "variant": variant.variant_id,
        "gene": variant.gene,
        "drug_name": data.get("name", ""),
        "recommendation": data.get("recommendation", ""),
    }


if __name__ == "__main__":
    import sys, json
    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: python -m scripts.pharma.query_pharmgkb 'chr10:96541616 G>A' [gene]"}))
        sys.exit(1)
    v = Variant.from_string(sys.argv[1])
    if len(sys.argv) > 2:
        v.gene = sys.argv[2]
    result = query_pharmgkb(v)
    print(json.dumps(result, indent=2))
