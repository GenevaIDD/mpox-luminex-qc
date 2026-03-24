"""Kit control bead QC: NC bead, ScG, FC, IC."""

import pandas as pd

from .config import FC_MFI_RANGE, IC_MFI_RANGE, NC_BEAD_MFI_MAX, SCG_MFI_MIN


def qc_kit_controls(df: pd.DataFrame) -> dict:
    """QC the 4 kit control beads across all wells.

    Returns dict with:
        nc_bead: DataFrame of per-well NC bead MFI with flags (> NC_BEAD_MFI_MAX)
        scg: DataFrame of per-well ScG MFI with flags (< SCG_MFI_MIN)
        fc: DataFrame of per-well FC MFI with flags (outside FC_MFI_RANGE)
        ic: DataFrame of per-well IC MFI with flags (outside IC_MFI_RANGE)
        flagged_wells: set of wells with any kit control flag
    """
    results = {}
    flagged_wells = set()

    # NC bead: should be ≤ NC_BEAD_MFI_MAX (150)
    nc_bead = df[df["analyte"] == "NC"][["well", "sample_name", "mfi"]].copy()
    nc_bead["flag"] = nc_bead["mfi"] > NC_BEAD_MFI_MAX
    results["nc_bead"] = nc_bead
    flagged_wells.update(nc_bead[nc_bead["flag"]]["well"].tolist())

    # ScG: sample addition control — should be ≥ SCG_MFI_MIN (10,000)
    scg = df[df["analyte"] == "ScG"][["well", "sample_name", "mfi"]].copy()
    scg["flag"] = scg["mfi"] < SCG_MFI_MIN
    results["scg"] = scg
    flagged_wells.update(scg[scg["flag"]]["well"].tolist())

    # FC: fluorescent conjugate control — should be within FC_MFI_RANGE (2000–5000)
    fc = df[df["analyte"] == "FC"][["well", "sample_name", "mfi"]].copy()
    fc["flag"] = (fc["mfi"] < FC_MFI_RANGE[0]) | (fc["mfi"] > FC_MFI_RANGE[1])
    results["fc"] = fc
    flagged_wells.update(fc[fc["flag"]]["well"].tolist())

    # IC: instrument control — should be within IC_MFI_RANGE (2000–3000)
    ic = df[df["analyte"] == "IC"][["well", "sample_name", "mfi"]].copy()
    ic["flag"] = (ic["mfi"] < IC_MFI_RANGE[0]) | (ic["mfi"] > IC_MFI_RANGE[1])
    results["ic"] = ic
    flagged_wells.update(ic[ic["flag"]]["well"].tolist())

    results["flagged_wells"] = flagged_wells
    return results
