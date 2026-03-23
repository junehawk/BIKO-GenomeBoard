"""Centralized configuration loader for GenomeBoard."""
import os
import threading
from pathlib import Path
from typing import Any, Optional

_config = None
_config_lock = threading.Lock()

def _find_project_root() -> Path:
    """Find project root by looking for config.yaml or known markers."""
    # Start from this file's location and go up
    current = Path(__file__).resolve().parent
    for _ in range(5):  # Max 5 levels up
        if (current / "config.yaml").exists():
            return current
        if (current / "scripts").exists() and (current / "data").exists():
            return current
        current = current.parent
    return Path.cwd()

def load_config(config_path: Optional[str] = None) -> dict:
    """Load config from YAML file with env var overrides."""
    global _config
    with _config_lock:
        if _config is not None and config_path is None:
            return _config

        project_root = _find_project_root()

        if config_path:
            path = Path(config_path)
        elif os.environ.get("GB_CONFIG_PATH"):
            path = Path(os.environ["GB_CONFIG_PATH"])
        else:
            path = project_root / "config.yaml"

        if path.exists():
            import yaml
            with open(path) as f:
                config = yaml.safe_load(f) or {}
        else:
            config = {}

        # Set defaults for missing values
        defaults = {
            "paths": {
                "krgdb": "data/krgdb_freq.tsv",
                "gene_knowledge": "data/gene_knowledge.json",
                "pgx_table": "data/korean_pgx_table.json",
                "acmg_rules": "data/acmg_rules.json",
                "templates": "templates",
            },
            "thresholds": {"ba1": 0.05, "bs1": 0.01, "pm2": 0.001},
            "api": {
                "clinvar_esearch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                "clinvar_esummary": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                "gnomad": "https://gnomad.broadinstitute.org/api",
                "gnomad_datasets": ["gnomad_r4", "gnomad_r2_1"],
                "pharmgkb": "https://api.pharmgkb.org/v1/data",
                "timeout": 30,
                "max_retries": 3,
                "backoff_base": 1.0,
                "pharmgkb_rate_limit": 0.5,
            },
            "report": {"default_mode": "cancer", "default_genome_build": "GRCh38"},
            "pgx": {
                "genes": ["CYP2D6", "CYP2C19", "CYP2C9", "HLA-B", "HLA-A", "NUDT15", "TPMT", "DPYD"],
                "risk_factor_genes": ["APOE"],
            },
        }

        # Deep merge defaults into config
        for section, values in defaults.items():
            if section not in config:
                config[section] = values
            elif isinstance(values, dict):
                for key, val in values.items():
                    config[section].setdefault(key, val)

        # Environment variable overrides
        if os.environ.get("NCBI_API_KEY"):
            config.setdefault("api", {})["ncbi_api_key"] = os.environ["NCBI_API_KEY"]

        # Resolve relative paths to absolute using project root
        config["_project_root"] = str(project_root)
        for key, val in config.get("paths", {}).items():
            if not os.path.isabs(val):
                config["paths"][key] = str(project_root / val)

        if config_path is None:
            _config = config
        return config

def get(key: str, default: Any = None) -> Any:
    """Get a config value using dot notation: get('thresholds.ba1')"""
    config = load_config()
    parts = key.split(".")
    current = config
    for part in parts:
        if isinstance(current, dict) and part in current:
            current = current[part]
        else:
            return default
    return current

def reset():
    """Reset cached config (useful for testing)."""
    global _config
    with _config_lock:
        _config = None
