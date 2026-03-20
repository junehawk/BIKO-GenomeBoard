# scripts/common/api_utils.py
import time
import requests
import requests.exceptions
from typing import Optional

def fetch_with_retry(
    url: str,
    params: Optional[dict] = None,
    headers: Optional[dict] = None,
    max_retries: int = 3,
    backoff_base: float = 1.0,
) -> Optional[dict]:
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, params=params, headers=headers, timeout=30)
            if 400 <= resp.status_code < 500:
                return None
            if resp.status_code >= 500:
                if attempt < max_retries - 1:
                    time.sleep(backoff_base * (2 ** attempt))
                continue
            return resp.json()
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            if attempt < max_retries - 1:
                time.sleep(backoff_base * (2 ** attempt))
    return None
