"""Configuration loader/saver — merges defaults with user YAML overrides."""

import copy
from pathlib import Path

import yaml

from .config import DEFAULTS


def get_config_path() -> Path:
    """Path to user config YAML file."""
    return Path.home() / "mpox-luminex-qc-results" / "config.yaml"


def _deep_merge(base: dict, overrides: dict) -> dict:
    """Recursively merge overrides into base dict."""
    result = copy.deepcopy(base)
    for key, val in overrides.items():
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = _deep_merge(result[key], val)
        else:
            result[key] = copy.deepcopy(val)
    return result


def load_config() -> dict:
    """Load config: defaults overlaid with user YAML overrides.

    Returns the full merged config dict.
    """
    config = copy.deepcopy(DEFAULTS)
    config_path = get_config_path()
    if config_path.exists():
        try:
            user_config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            if user_config and isinstance(user_config, dict):
                config = _deep_merge(config, user_config)
        except Exception:
            pass  # If YAML is corrupt, fall back to defaults
    return config


def save_config(config: dict) -> None:
    """Save full config to YAML file."""
    config_path = get_config_path()
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        yaml.dump(config, default_flow_style=False, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )


def reset_config() -> None:
    """Delete user config file, reverting to defaults."""
    config_path = get_config_path()
    if config_path.exists():
        config_path.unlink()


def get_antigen_names(config: dict) -> list[str]:
    """Extract antigen name list from config."""
    return [a["name"] for a in config["panel"]["antigens"]]


def get_kit_control_names(config: dict) -> list[str]:
    """Extract kit control name list from config."""
    return [c["name"] for c in config["panel"]["kit_controls"]]


def get_qc_thresholds(config: dict) -> dict:
    """Extract QC thresholds as a flat dict."""
    return config["qc_thresholds"]
