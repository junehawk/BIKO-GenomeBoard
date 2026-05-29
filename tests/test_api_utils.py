# tests/test_api_utils.py
from unittest.mock import MagicMock

import scripts.common.api_utils as api_utils_mod
from scripts.common.api_utils import fetch_with_retry


def _mock_session(mocker, get_side_effect=None, get_return=None):
    """Helper: return a mock session whose .get() behaves as specified."""
    mock_sess = MagicMock()
    if get_side_effect is not None:
        mock_sess.get.side_effect = get_side_effect
    elif get_return is not None:
        mock_sess.get.return_value = get_return
    mocker.patch.object(api_utils_mod, "_get_session", return_value=mock_sess)
    return mock_sess


def test_fetch_success(mocker):
    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.json.return_value = {"result": "ok"}
    _mock_session(mocker, get_return=mock_resp)

    result = fetch_with_retry("https://example.com/api")
    assert result == {"result": "ok"}


def test_fetch_retry_on_failure(mocker):
    fail_resp = MagicMock()
    fail_resp.status_code = 500
    fail_resp.raise_for_status.side_effect = Exception("500")

    ok_resp = MagicMock()
    ok_resp.status_code = 200
    ok_resp.json.return_value = {"result": "ok"}

    _mock_session(mocker, get_side_effect=[fail_resp, ok_resp])
    mocker.patch("time.sleep")  # skip actual sleep

    result = fetch_with_retry("https://example.com/api", max_retries=2)
    assert result == {"result": "ok"}


def test_fetch_all_retries_fail(mocker):
    fail_resp = MagicMock()
    fail_resp.status_code = 500
    fail_resp.raise_for_status.side_effect = Exception("500")

    _mock_session(mocker, get_return=fail_resp)
    mocker.patch("time.sleep")

    result = fetch_with_retry("https://example.com/api", max_retries=3)
    assert result is None


def test_fetch_retries_on_429_honoring_retry_after(mocker):
    """429 is retryable (not a hard 4xx) and Retry-After is honored (ENRI-01)."""
    busy = MagicMock()
    busy.status_code = 429
    busy.headers = {"Retry-After": "0"}

    ok = MagicMock()
    ok.status_code = 200
    ok.json.return_value = {"result": "ok"}

    _mock_session(mocker, get_side_effect=[busy, ok])
    sleep_mock = mocker.patch("time.sleep")

    result = fetch_with_retry("https://example.com/api", max_retries=3)
    assert result == {"result": "ok"}
    sleep_mock.assert_called_once_with(0.0)  # Retry-After respected, not exponential


def test_fetch_429_exhausted_returns_none(mocker):
    """A persistent 429 (no Retry-After) eventually returns None, not on first hit."""
    busy = MagicMock()
    busy.status_code = 429
    busy.headers = {}
    sess = _mock_session(mocker, get_return=busy)
    mocker.patch("time.sleep")

    result = fetch_with_retry("https://example.com/api", max_retries=3)
    assert result is None
    assert sess.get.call_count == 3  # retried, not abandoned after the first 429
