# scripts/common/paperclip_client.py
import os
from typing import Optional, Dict, List
from scripts.common.api_utils import fetch_with_retry
import requests

PAPERCLIP_BASE = os.environ.get("PAPERCLIP_URL", "http://localhost:3100")

def _url(path: str) -> str:
    return f"{PAPERCLIP_BASE}/api{path}"

def create_issue(company_id: str, title: str, body: str, assignee: Optional[str] = None) -> Optional[Dict]:
    """Create a new issue/ticket in Paperclip."""
    payload = {"title": title, "body": body}
    if assignee:
        payload["assignee"] = assignee
    try:
        resp = requests.post(_url(f"/companies/{company_id}/issues"), json=payload, timeout=10)
        resp.raise_for_status()
        return resp.json()
    except Exception:
        return None

def add_comment(issue_id: str, body: str, author: Optional[str] = None) -> Optional[Dict]:
    """Add a comment to an existing issue."""
    payload = {"body": body}
    if author:
        payload["author"] = author
    try:
        resp = requests.post(_url(f"/issues/{issue_id}/comments"), json=payload, timeout=10)
        resp.raise_for_status()
        return resp.json()
    except Exception:
        return None

def update_issue_status(issue_id: str, status: str) -> Optional[Dict]:
    """Update issue status (open, in_progress, resolved, closed)."""
    try:
        resp = requests.patch(_url(f"/issues/{issue_id}"), json={"status": status}, timeout=10)
        resp.raise_for_status()
        return resp.json()
    except Exception:
        return None

def get_assigned_issues(employee_id: str) -> List[Dict]:
    """Get all issues assigned to an employee."""
    data = fetch_with_retry(_url(f"/employees/{employee_id}/issues"))
    return data if isinstance(data, list) else []
