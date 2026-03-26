"""Classify wells as PC standard, NC, or specimen from sample names."""

import re

import pandas as pd

from .config import PC_PATTERNS, NC_PATTERNS


def classify_wells(df: pd.DataFrame, config: dict | None = None) -> pd.DataFrame:
    """Add well_type and dilution columns based on sample_name.

    well_type: 'pc', 'nc', or 'specimen'
    dilution: numeric dilution denominator for PC wells, NaN otherwise
    """
    pc_pats = PC_PATTERNS
    nc_pats = NC_PATTERNS
    if config is not None:
        wc = config.get("well_classification", {})
        pc_pats = wc.get("pc_patterns", pc_pats)
        nc_pats = wc.get("nc_patterns", nc_pats)

    df = df.copy()
    df["well_type"] = df["sample_name"].apply(lambda n: _classify_sample(n, pc_pats, nc_pats))
    df["dilution"] = df["sample_name"].apply(_extract_dilution)
    return df


def _classify_sample(name: str, pc_patterns: list[str], nc_patterns: list[str]) -> str:
    name = name.strip()
    for pat in pc_patterns:
        if re.match(pat, name, re.IGNORECASE):
            return "pc"
    for pat in nc_patterns:
        if re.match(pat, name, re.IGNORECASE):
            return "nc"
    return "specimen"


def _extract_dilution(name: str) -> float:
    """Extract dilution denominator from PC sample names like 'PC 1:50' or 'ITM PC 50'."""
    # Match "PC 1:50" or "ITM PC 1:50"
    m = re.search(r"1:(\d+)", name)
    if m:
        return float(m.group(1))
    # Match "PC 50" or "ITM PC 50" (just the number after PC)
    m = re.search(r"PC\s+(\d+)", name, re.IGNORECASE)
    if m:
        return float(m.group(1))
    return float("nan")
