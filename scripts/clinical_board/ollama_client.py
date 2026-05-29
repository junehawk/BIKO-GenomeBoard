"""Ollama REST API client for Clinical Board."""

import json
import logging
import time

import requests

from scripts.common.config import get

logger = logging.getLogger(__name__)

DEFAULT_URL = "http://localhost:11434"


class OllamaClient:
    """Simple client for the Ollama REST API.

    All methods handle connection failures gracefully — they return empty/error
    values rather than raising, so callers never crash when Ollama is offline.
    """

    def __init__(self, base_url: str = None, timeout: int = None):
        self.base_url = base_url or get("clinical_board.ollama_url", DEFAULT_URL)
        self.timeout = timeout or get("clinical_board.timeout", 120)
        # Deterministic sampling seed so the Board narrative is reproducible
        # run-to-run (a stable research-reference document). Set
        # clinical_board.seed: null in config to restore non-deterministic
        # sampling (v2.7 review — reproducibility).
        self.seed = get("clinical_board.seed", 42)

    def _sampling_options(self, **extra) -> dict:
        """Build the Ollama ``options`` dict, injecting the reproducibility seed."""
        options = dict(extra)
        if self.seed is not None:
            options["seed"] = self.seed
        return options

    # ── Health & Discovery ───────────────────────────────────────────────────

    def is_available(self) -> bool:
        """Check if Ollama server is running."""
        try:
            resp = requests.get(
                f"{self.base_url}/api/tags",
                timeout=5,
            )
            return resp.status_code == 200
        except (requests.ConnectionError, requests.Timeout, OSError):
            return False

    def list_models(self) -> list[str]:
        """List available models on the Ollama server."""
        try:
            resp = requests.get(
                f"{self.base_url}/api/tags",
                timeout=10,
            )
            resp.raise_for_status()
            data = resp.json()
            return [m["name"] for m in data.get("models", [])]
        except (requests.ConnectionError, requests.Timeout, OSError, requests.HTTPError, KeyError, ValueError):
            logger.debug("Failed to list Ollama models")
            return []

    def has_model(self, model: str) -> bool:
        """Check if a specific model is available."""
        models = self.list_models()
        # Exact match or prefix match (e.g., "medgemma:27b" matches "medgemma:27b")
        for m in models:
            if m == model or m.startswith(model + ":") or model.startswith(m.split(":")[0]):
                return True
        return model in models

    # ── Generation ───────────────────────────────────────────────────────────

    def generate(
        self,
        model: str,
        prompt: str,
        system: str = "",
        temperature: float = 0.3,
        max_retries: int = None,
    ) -> str:
        """Generate text completion.

        Returns the response text on success, or an error description string
        on failure (never raises).
        """
        if max_retries is None:
            max_retries = get("clinical_board.max_retries", 2)

        payload = {
            "model": model,
            "prompt": prompt,
            "stream": False,
            "options": self._sampling_options(temperature=temperature),
        }
        if system:
            payload["system"] = system

        last_error = ""
        for attempt in range(max_retries + 1):
            t0 = time.time()
            try:
                logger.debug(
                    "Ollama generate: model=%s prompt_len=%d attempt=%d",
                    model,
                    len(prompt),
                    attempt,
                )
                resp = requests.post(
                    f"{self.base_url}/api/generate",
                    json=payload,
                    timeout=self.timeout,
                )
                resp.raise_for_status()
                data = resp.json()
                response_text = data.get("response", "")
                elapsed = time.time() - t0
                logger.debug(
                    "Ollama response: model=%s response_len=%d elapsed=%.1fs",
                    model,
                    len(response_text),
                    elapsed,
                )
                return response_text

            except requests.Timeout:
                last_error = f"Timeout after {self.timeout}s (attempt {attempt + 1})"
                logger.warning("Ollama timeout: %s", last_error)
            except requests.ConnectionError:
                last_error = "Connection refused — is Ollama running?"
                logger.warning("Ollama connection error: %s", last_error)
            except (requests.HTTPError, ValueError, KeyError) as exc:
                last_error = f"Request error: {exc}"
                logger.warning("Ollama error: %s", last_error)

            # Exponential backoff before retry
            if attempt < max_retries:
                backoff = 2**attempt
                logger.debug("Retrying in %ds...", backoff)
                time.sleep(backoff)

        return f"[Error] {last_error}"

    def generate_json(
        self,
        model: str,
        prompt: str,
        system: str = "",
        temperature: float = 0.1,
        max_retries: int = None,
    ) -> dict:
        """Generate structured JSON output.

        Returns a parsed dict on success, or an empty dict on failure. Like
        :meth:`generate`, transient transport failures (timeout / connection /
        HTTP) are retried with exponential backoff — the entire Clinical Board
        synthesis depends on this call, so a single transient stall must not
        silently empty the report (v2.7 review BOAR-01). Truncation at the
        ``num_predict`` ceiling is detected and logged (BOAR-03).
        """
        if max_retries is None:
            max_retries = get("clinical_board.max_retries", 2)
        num_predict = get("clinical_board.num_predict", 2048)

        payload = {
            "model": model,
            "prompt": prompt,
            "stream": False,
            "format": "json",
            "options": self._sampling_options(temperature=temperature, num_predict=num_predict),
        }
        if system:
            payload["system"] = system

        last_error = ""
        for attempt in range(max_retries + 1):
            t0 = time.time()
            try:
                logger.debug(
                    "Ollama generate_json: model=%s prompt_len=%d attempt=%d",
                    model,
                    len(prompt),
                    attempt,
                )
                resp = requests.post(
                    f"{self.base_url}/api/generate",
                    json=payload,
                    timeout=self.timeout,
                )
                resp.raise_for_status()
                data = resp.json()
                response_text = data.get("response", "")
                elapsed = time.time() - t0
                logger.debug(
                    "Ollama JSON response: model=%s response_len=%d elapsed=%.1fs",
                    model,
                    len(response_text),
                    elapsed,
                )

                # Truncation detection: a response cut off at the token ceiling
                # almost always yields invalid/incomplete JSON. Surface it.
                if data.get("done_reason") == "length":
                    logger.warning(
                        "Ollama JSON for model=%s was truncated at num_predict=%d "
                        "(done_reason=length); JSON may be incomplete — raise "
                        "clinical_board.num_predict if this recurs",
                        model,
                        num_predict,
                    )

                # A successful HTTP response is authoritative — do not retry a
                # parse failure (it will not differ on the next attempt at this
                # temperature). Retry is reserved for transport errors below.
                try:
                    return json.loads(response_text)
                except json.JSONDecodeError:
                    return _extract_json(response_text)

            except requests.Timeout:
                last_error = f"Timeout after {self.timeout}s (attempt {attempt + 1})"
                logger.warning("Ollama JSON timeout: %s", last_error)
            except requests.ConnectionError:
                last_error = "Connection refused — is Ollama running?"
                logger.warning("Ollama JSON connection error: %s", last_error)
            except (requests.HTTPError, ValueError, KeyError) as exc:
                last_error = f"Request error: {exc}"
                logger.warning("Ollama JSON error: %s", last_error)

            if attempt < max_retries:
                backoff = 2**attempt
                logger.debug("Retrying generate_json in %ds...", backoff)
                time.sleep(backoff)

        return {}


def _extract_json(text: str) -> dict:
    """Best-effort extraction of a JSON object from free-form text.

    Looks for the first ``{...}`` block and attempts to parse it.
    Returns an empty dict if nothing parseable is found.
    """
    # Find the outermost { ... }
    start = text.find("{")
    if start == -1:
        return {}

    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                try:
                    return json.loads(text[start : i + 1])
                except json.JSONDecodeError:
                    return {}
    return {}
