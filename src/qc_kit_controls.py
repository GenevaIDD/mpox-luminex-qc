"""Kit control bead QC: ScG, FC, IC, NC bead."""

import numpy as np
import pandas as pd

from .config import NC_BEAD_MFI_MAX


def qc_kit_controls(df: pd.DataFrame) -> dict:
    """QC the 4 kit control beads across all wells.

    Returns dict with:
        nc_bead: DataFrame of per-well NC bead (09) MFI with flags (> NC_BEAD_MFI_MAX)
        scg: DataFrame of per-well ScG (10) MFI with flags (low outliers)
        fc: DataFrame of per-well FC (11) MFI with flags (low outliers)
        ic: DataFrame of per-well IC (12) MFI with flags (outliers)
        flagged_wells: set of wells with any kit control flag
    """
    results = {}
    flagged_wells = set()

    # NC bead (09): should be < NC_BEAD_MFI_MAX
    nc_bead = df[df["analyte"] == "09 NC"][["well", "sample_name", "mfi"]].copy()
    nc_bead["flag"] = nc_bead["mfi"] > NC_BEAD_MFI_MAX
    results["nc_bead"] = nc_bead
    flagged_wells.update(nc_bead[nc_bead["flag"]]["well"].tolist())

    # ScG (10): sample addition control — flag low outliers
    scg = df[df["analyte"] == "10 ScG"][["well", "sample_name", "mfi"]].copy()
    scg["flag"] = _flag_low_outliers(scg["mfi"])
    results["scg"] = scg
    flagged_wells.update(scg[scg["flag"]]["well"].tolist())

    # FC (11): fluorescent conjugate control — flag low outliers
    fc = df[df["analyte"] == "11 FC"][["well", "sample_name", "mfi"]].copy()
    fc["flag"] = _flag_low_outliers(fc["mfi"])
    results["fc"] = fc
    flagged_wells.update(fc[fc["flag"]]["well"].tolist())

    # IC (12): instrument control — flag any outliers
    ic = df[df["analyte"] == "12 IC"][["well", "sample_name", "mfi"]].copy()
    ic["flag"] = _flag_outliers(ic["mfi"])
    results["ic"] = ic
    flagged_wells.update(ic[ic["flag"]]["well"].tolist())

    results["flagged_wells"] = flagged_wells
    return results


def _flag_low_outliers(values: pd.Series, k: float = 3.0) -> pd.Series:
    """Flag values below median - k * MAD."""
    median = values.median()
    mad = np.median(np.abs(values - median))
    if mad == 0:
        mad = values.std() * 0.6745  # fallback
    threshold = median - k * mad
    return values < threshold


def _flag_outliers(values: pd.Series, k: float = 3.0) -> pd.Series:
    """Flag values beyond median +/- k * MAD."""
    median = values.median()
    mad = np.median(np.abs(values - median))
    if mad == 0:
        mad = values.std() * 0.6745
    return (values < median - k * mad) | (values > median + k * mad)
