"""PC standard curve replicate variability QC."""

import numpy as np
import pandas as pd

from .config import MPXV_ANTIGENS, PC_CV_THRESHOLD


def qc_pc_replicates(df: pd.DataFrame, cv_threshold: float | None = None, config: dict | None = None) -> dict:
    """Compute CV for PC standard duplicates.

    Returns dict with:
        has_replicates: bool — whether duplicate PC wells were found
        replicate_cv: DataFrame with columns [analyte, dilution, rep1, rep2, mean, cv, flag]
        n_flagged: number of pairs exceeding cv_threshold
    """
    if cv_threshold is None:
        if config is not None:
            cv_threshold = config.get("qc_thresholds", {}).get("pc_cv_threshold", PC_CV_THRESHOLD)
        else:
            cv_threshold = PC_CV_THRESHOLD

    antigens = MPXV_ANTIGENS
    if config is not None:
        antigens = [a["name"] for a in config.get("panel", {}).get("antigens", [])] or antigens

    pc = df[(df["well_type"] == "pc") & (df["analyte"].isin(antigens))].copy()

    # Group by pool if available
    has_pool = "pc_pool" in pc.columns
    group_cols = (["pc_pool", "analyte", "dilution"] if has_pool
                  else ["analyte", "dilution"])

    # Count replicates per group
    rep_counts = pc.groupby(group_cols)["mfi"].count()
    has_replicates = (rep_counts > 1).any()

    if not has_replicates:
        return {"has_replicates": False, "replicate_cv": pd.DataFrame(), "n_flagged": 0}

    rows = []
    for keys, group in pc.groupby(group_cols):
        if has_pool:
            pool, analyte, dilution = keys
        else:
            analyte, dilution = keys
            pool = None
        values = group["mfi"].values
        if len(values) < 2:
            continue
        rep1, rep2 = values[0], values[1]
        mean_val = np.mean([rep1, rep2])
        cv = np.std([rep1, rep2], ddof=0) / mean_val if mean_val > 0 else np.nan
        row = {
            "analyte": analyte,
            "dilution": dilution,
            "rep1": rep1,
            "rep2": rep2,
            "mean": mean_val,
            "cv": cv,
            "flag": cv > cv_threshold if not np.isnan(cv) else False,
        }
        if has_pool:
            row["pc_pool"] = pool
        rows.append(row)

    cv_df = pd.DataFrame(rows)
    n_flagged = cv_df["flag"].sum() if len(cv_df) > 0 else 0

    return {
        "has_replicates": True,
        "replicate_cv": cv_df,
        "n_flagged": int(n_flagged),
    }
