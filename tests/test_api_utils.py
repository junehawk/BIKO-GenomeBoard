# tests/test_api_utils.py
from unittest.mock import MagicMock
from scripts.common.api_utils import fetch_with_retry


def test_fetch_success(mocker):
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.json.return_value = {"result": "ok"}
    mocker.patch("requests.get", return_value=mock_resp)

    result = fetch_with_retry("https://example.com/api")
    assert result == {"result": "ok"}


def test_fetch_retry_on_failure(mocker):
    fail_resp = MagicMock()
    fail_resp.status_code = 500
    fail_resp.raise_for_status.side_effect = Exception("500")

    ok_resp = MagicMock()
    ok_resp.status_code = 200
    ok_resp.json.return_value = {"result": "ok"}

    mocker.patch("requests.get", side_effect=[fail_resp, ok_resp])
    mocker.patch("time.sleep")  # skip actual sleep

    result = fetch_with_retry("https://example.com/api", max_retries=2)
    assert result == {"result": "ok"}


def test_fetch_all_retries_fail(mocker):
    fail_resp = MagicMock()
    fail_resp.status_code = 500
    fail_resp.raise_for_status.side_effect = Exception("500")

    mocker.patch("requests.get", return_value=fail_resp)
    mocker.patch("time.sleep")

    result = fetch_with_retry("https://example.com/api", max_retries=3)
    assert result is None
