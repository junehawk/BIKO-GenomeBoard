"""Tests for OllamaClient — all Ollama interactions are mocked."""
import json
from unittest.mock import patch, MagicMock

import pytest
import requests

from scripts.clinical_board.ollama_client import OllamaClient, _extract_json


# ── Fixtures ─────────────────────────────────────────────────────────────────

@pytest.fixture
def client():
    """OllamaClient with explicit settings (no config.yaml dependency)."""
    return OllamaClient(base_url="http://localhost:11434", timeout=30)


def _mock_response(json_data=None, status_code=200, text=""):
    """Create a mock requests.Response."""
    resp = MagicMock(spec=requests.Response)
    resp.status_code = status_code
    resp.text = text or json.dumps(json_data or {})
    resp.json.return_value = json_data or {}
    resp.raise_for_status.return_value = None
    if status_code >= 400:
        resp.raise_for_status.side_effect = requests.HTTPError(
            f"{status_code} Error", response=resp,
        )
    return resp


# ── is_available ─────────────────────────────────────────────────────────────

@patch("scripts.clinical_board.ollama_client.requests.get")
def test_is_available_success(mock_get, client):
    mock_get.return_value = _mock_response({"models": []}, 200)
    assert client.is_available() is True
    mock_get.assert_called_once()


@patch("scripts.clinical_board.ollama_client.requests.get")
def test_is_available_failure(mock_get, client):
    mock_get.side_effect = requests.ConnectionError("refused")
    assert client.is_available() is False


@patch("scripts.clinical_board.ollama_client.requests.get")
def test_is_available_timeout(mock_get, client):
    mock_get.side_effect = requests.Timeout("timed out")
    assert client.is_available() is False


# ── list_models / has_model ──────────────────────────────────────────────────

@patch("scripts.clinical_board.ollama_client.requests.get")
def test_list_models(mock_get, client):
    mock_get.return_value = _mock_response({
        "models": [
            {"name": "medgemma:27b"},
            {"name": "gemma4:31b"},
        ]
    })
    models = client.list_models()
    assert models == ["medgemma:27b", "gemma4:31b"]


@patch("scripts.clinical_board.ollama_client.requests.get")
def test_list_models_error(mock_get, client):
    mock_get.side_effect = requests.ConnectionError("refused")
    assert client.list_models() == []


@patch("scripts.clinical_board.ollama_client.requests.get")
def test_has_model_found(mock_get, client):
    mock_get.return_value = _mock_response({
        "models": [{"name": "medgemma:27b"}, {"name": "gemma4:31b"}]
    })
    assert client.has_model("medgemma:27b") is True


@patch("scripts.clinical_board.ollama_client.requests.get")
def test_has_model_not_found(mock_get, client):
    mock_get.return_value = _mock_response({
        "models": [{"name": "gemma4:31b"}]
    })
    assert client.has_model("medgemma:27b") is False


# ── generate ─────────────────────────────────────────────────────────────────

@patch("scripts.clinical_board.ollama_client.requests.post")
def test_generate_success(mock_post, client):
    mock_post.return_value = _mock_response({
        "response": "The TP53 variant is likely pathogenic.",
    })
    result = client.generate("medgemma:27b", "Analyze TP53 R175L")
    assert "TP53" in result
    assert "pathogenic" in result.lower()

    # Verify the request payload
    call_args = mock_post.call_args
    payload = call_args[1]["json"]
    assert payload["model"] == "medgemma:27b"
    assert payload["stream"] is False
    assert "prompt" in payload


@patch("scripts.clinical_board.ollama_client.requests.post")
def test_generate_with_system_prompt(mock_post, client):
    mock_post.return_value = _mock_response({"response": "OK"})
    client.generate("medgemma:27b", "prompt", system="You are a doctor.")
    payload = mock_post.call_args[1]["json"]
    assert payload["system"] == "You are a doctor."


@patch("scripts.clinical_board.ollama_client.time.sleep")  # skip backoff
@patch("scripts.clinical_board.ollama_client.requests.post")
def test_generate_timeout_retry(mock_post, mock_sleep, client):
    """Timeout on first attempt, success on retry."""
    mock_post.side_effect = [
        requests.Timeout("timed out"),
        _mock_response({"response": "retried OK"}),
    ]
    result = client.generate("medgemma:27b", "test", max_retries=1)
    assert result == "retried OK"
    assert mock_post.call_count == 2


@patch("scripts.clinical_board.ollama_client.time.sleep")
@patch("scripts.clinical_board.ollama_client.requests.post")
def test_generate_all_retries_fail(mock_post, mock_sleep, client):
    """All retries exhausted → error message returned."""
    mock_post.side_effect = requests.Timeout("timed out")
    result = client.generate("medgemma:27b", "test", max_retries=2)
    assert result.startswith("[Error]")
    assert "Timeout" in result
    assert mock_post.call_count == 3  # initial + 2 retries


@patch("scripts.clinical_board.ollama_client.requests.post")
def test_generate_connection_error(mock_post, client):
    """Connection refused returns error string, does not raise."""
    mock_post.side_effect = requests.ConnectionError("refused")
    result = client.generate("medgemma:27b", "test", max_retries=0)
    assert "[Error]" in result
    assert "Connection refused" in result


# ── generate_json ────────────────────────────────────────────────────────────

@patch("scripts.clinical_board.ollama_client.requests.post")
def test_generate_json_success(mock_post, client):
    expected = {"diagnosis": "Li-Fraumeni syndrome", "confidence": "high"}
    mock_post.return_value = _mock_response({
        "response": json.dumps(expected),
    })
    result = client.generate_json("medgemma:27b", "Analyze this case")
    assert result == expected

    # Verify format: json is requested
    payload = mock_post.call_args[1]["json"]
    assert payload["format"] == "json"


@patch("scripts.clinical_board.ollama_client.requests.post")
def test_generate_json_malformed(mock_post, client):
    """Non-JSON response with embedded JSON → extracted."""
    raw = 'Here is the analysis: {"diagnosis": "BRCA2-related", "confidence": "moderate"} end.'
    mock_post.return_value = _mock_response({"response": raw})
    result = client.generate_json("medgemma:27b", "test")
    assert result["diagnosis"] == "BRCA2-related"


@patch("scripts.clinical_board.ollama_client.requests.post")
def test_generate_json_completely_unparseable(mock_post, client):
    """Completely non-JSON response → empty dict."""
    mock_post.return_value = _mock_response({"response": "no json here at all"})
    result = client.generate_json("medgemma:27b", "test")
    assert result == {}


@patch("scripts.clinical_board.ollama_client.requests.post")
def test_generate_json_connection_error(mock_post, client):
    mock_post.side_effect = requests.ConnectionError("refused")
    result = client.generate_json("medgemma:27b", "test")
    assert result == {}


# ── _extract_json helper ─────────────────────────────────────────────────────

def test_extract_json_clean():
    assert _extract_json('{"a": 1}') == {"a": 1}


def test_extract_json_embedded():
    text = 'Here: {"key": "val"} done.'
    assert _extract_json(text) == {"key": "val"}


def test_extract_json_nested():
    text = '{"outer": {"inner": 1}}'
    assert _extract_json(text) == {"outer": {"inner": 1}}


def test_extract_json_no_json():
    assert _extract_json("plain text") == {}


def test_extract_json_empty():
    assert _extract_json("") == {}


# ── Default config ───────────────────────────────────────────────────────────

def test_default_config():
    """Client uses config defaults when no explicit args are given."""
    with patch("scripts.clinical_board.ollama_client.get") as mock_get:
        mock_get.side_effect = lambda key, default=None: {
            "clinical_board.ollama_url": "http://myhost:11434",
            "clinical_board.timeout": 60,
        }.get(key, default)
        c = OllamaClient()
        assert c.base_url == "http://myhost:11434"
        assert c.timeout == 60


def test_explicit_args_override_config():
    """Explicit constructor args take precedence over config."""
    c = OllamaClient(base_url="http://custom:1234", timeout=99)
    assert c.base_url == "http://custom:1234"
    assert c.timeout == 99
