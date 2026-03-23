"""Negative control well monitoring."""

import pandas as pd

from .config import MPXV_ANTIGENS


def qc_nc_levels(df: pd.DataFrame) -> pd.DataFrame:
    """Extract NC well MFI values per antigen.

    Returns DataFrame with columns [well, analyte, mfi] for NC wells,
    filtered to the 8 mpox antigens.
    """
    nc = df[(df["well_type"] == "nc") & (df["analyte"].isin(MPXV_ANTIGENS))].copy()
    return nc[["well", "sample_name", "analyte", "mfi"]].reset_index(drop=True)
