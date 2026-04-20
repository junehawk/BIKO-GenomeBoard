# scripts/pharma/query_pharmgkb.py
import threading
import time
from typing import Dict, Optional

from scripts.common.api_utils import fetch_with_retry
from scripts.common.config import get
from scripts.common.models import Variant

PHARMGKB_API = get("api.pharmgkb", "https://api.pharmgkb.org/v1/data")
_last_request_time = 0
_RATE_LOCK = threading.Lock()


def _rate_limit():
    global _last_request_time
    with _RATE_LOCK:
        elapsed = time.time() - _last_request_time
        rate_limit = get("api.pharmgkb_rate_limit", 0.5)
        if elapsed < rate_limit:
            time.sleep(rate_limit - elapsed)
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
    import json
    import sys

    if len(sys.argv) < 2:
        print(
            json.dumps(
                {"error": "Usage: python -m scripts.pharmacogenomics.query_pharmgkb 'chr10:96541616 G>A' [gene]"}
            )
        )
        sys.exit(1)
    v = Variant.from_string(sys.argv[1])
    if len(sys.argv) > 2:
        v.gene = sys.argv[2]
    result = query_pharmgkb(v)
    print(json.dumps(result, indent=2))
