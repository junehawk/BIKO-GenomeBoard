# scripts/common/api_utils.py
import requests
import requests.exceptions
import time
import threading
from typing import Optional
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

            if resp.status_code >= 500:
                if attempt < max_retries - 1:
                    time.sleep(backoff_base * (2 ** attempt))
                    continue
            elif resp.status_code >= 400:
                return None

            resp.raise_for_status()
            return resp.json()
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            if attempt < max_retries - 1:
                time.sleep(backoff_base * (2 ** attempt))
        except Exception:
            return None
    return None
