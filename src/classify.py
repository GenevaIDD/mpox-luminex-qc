"""Classify wells as PC standard, NC, or specimen from sample names."""

import re

import pandas as pd

def classify_wells(df: pd.DataFrame) -> pd.DataFrame:
    """Add well_type and dilution columns based on sample_name.

    well_type: 'pc', 'nc', or 'specimen'
    dilution: numeric dilution denominator for PC wells, NaN otherwise
    """
    df = df.copy()
    df["well_type"] = df["sample_name"].apply(_classify_sample)
    df["dilution"] = df["sample_name"].apply(_extract_dilution)
    return df


def _classify_sample(name: str) -> str:
    name = name.strip()
    if re.match(r"^(PC|ITM\s*PC)\s", name, re.IGNORECASE):
        return "pc"
    if re.match(r"^(NC|ITM\s*NC)", name, re.IGNORECASE):
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
