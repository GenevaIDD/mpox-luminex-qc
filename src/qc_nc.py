"""Negative control well monitoring."""

import pandas as pd

from .config import MPXV_ANTIGENS


def qc_nc_levels(df: pd.DataFrame, config: dict | None = None) -> pd.DataFrame:
    """Extract NC well MFI values per antigen.

    Returns DataFrame with columns [well, analyte, mfi] for NC wells,
    filtered to the antigens defined in config (or defaults).
    """
    antigens = MPXV_ANTIGENS
    if config is not None:
        antigens = [a["name"] for a in config.get("panel", {}).get("antigens", [])] or antigens

    nc = df[(df["well_type"] == "nc") & (df["analyte"].isin(antigens))].copy()
    return nc[["well", "sample_name", "analyte", "mfi"]].reset_index(drop=True)
