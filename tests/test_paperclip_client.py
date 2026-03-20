# tests/test_paperclip_client.py
import pytest
from unittest.mock import MagicMock
from scripts.common.paperclip_client import create_issue, add_comment, update_issue_status, get_assigned_issues

def test_create_issue(mocker):
    mock_resp = MagicMock()
    mock_resp.status_code = 201
    mock_resp.json.return_value = {"id": "issue-1", "title": "Analyze TP53"}
    mocker.patch("requests.post", return_value=mock_resp)

    result = create_issue("genomeboard", "Analyze TP53", "Variant: chr17:7577120 G>A", assignee="cto")
    assert result is not None
    assert result["id"] == "issue-1"

def test_add_comment(mocker):
    mock_resp = MagicMock()
    mock_resp.status_code = 201
    mock_resp.json.return_value = {"id": "comment-1", "body": "ClinVar: Pathogenic"}
    mocker.patch("requests.post", return_value=mock_resp)

    result = add_comment("issue-1", "ClinVar: Pathogenic", author="clinical")
    assert result is not None
    assert result["body"] == "ClinVar: Pathogenic"

def test_update_issue_status(mocker):
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.json.return_value = {"id": "issue-1", "status": "resolved"}
    mocker.patch("requests.patch", return_value=mock_resp)

    result = update_issue_status("issue-1", "resolved")
    assert result is not None
    assert result["status"] == "resolved"

def test_get_assigned_issues(mocker):
    mocker.patch(
        "scripts.common.paperclip_client.fetch_with_retry",
        return_value=[{"id": "issue-1", "title": "Analyze TP53"}]
    )
    result = get_assigned_issues("cto")
    assert len(result) == 1
    assert result[0]["title"] == "Analyze TP53"

def test_create_issue_failure(mocker):
    mocker.patch("requests.post", side_effect=Exception("connection refused"))
    result = create_issue("genomeboard", "Test", "Test body")
    assert result is None
