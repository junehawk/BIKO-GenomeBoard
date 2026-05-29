# scripts/common/api_utils.py
import threading
import time
from typing import Optional

import requests
import requests.exceptions

from scripts.common.config import get

_session = None
_session_lock = threading.Lock()


def _get_session() -> requests.Session:
    global _session
    with _session_lock:
        if _session is None:
            _session = requests.Session()
            _session.headers.update({"User-Agent": "BIKO GenomeBoard/0.9"})
        return _session


def fetch_with_retry(
    url: str,
    params: Optional[dict] = None,
    headers: Optional[dict] = None,
    max_retries: int = None,
    backoff_base: float = None,
    method: str = "GET",
    json_data: Optional[dict] = None,
) -> Optional[dict]:
    max_retries = max_retries if max_retries is not None else get("api.max_retries", 3)
    backoff_base = backoff_base if backoff_base is not None else get("api.backoff_base", 1.0)
    timeout = get("api.timeout", 30)
    session = _get_session()

    for attempt in range(max_retries):
        try:
            if method == "POST":
                resp = session.post(url, params=params, headers=headers, json=json_data, timeout=timeout)
            else:
                resp = session.get(url, params=params, headers=headers, timeout=timeout)

            if resp.status_code == 429:
                # Rate limited — retryable, not a hard 4xx. Honor Retry-After
                # when the server provides it, else fall back to exponential
                # backoff (v2.7 review ENRI-01).
                if attempt < max_retries - 1:
                    retry_after = resp.headers.get("Retry-After")
                    try:
                        delay = float(retry_after) if retry_after else backoff_base * (2**attempt)
                    except (TypeError, ValueError):
                        delay = backoff_base * (2**attempt)
                    time.sleep(delay)
                    continue
                return None
            if resp.status_code >= 500:
                if attempt < max_retries - 1:
                    time.sleep(backoff_base * (2**attempt))
                    continue
            elif resp.status_code >= 400:
                return None

            resp.raise_for_status()
            return resp.json()
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            if attempt < max_retries - 1:
                time.sleep(backoff_base * (2**attempt))
        except Exception:
            return None
    return None
