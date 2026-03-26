"""Bead count QC — flag wells with insufficient bead counts."""

import pandas as pd

from .config import BEAD_COUNT_MIN


def qc_bead_counts(df: pd.DataFrame, min_count: int | None = None, config: dict | None = None) -> dict:
    """Check bead counts and flag well-analyte pairs below threshold.

    Returns dict with:
        flagged: DataFrame of well-analyte pairs with count < min_count
        by_well: DataFrame with median bead count per well (for heatmap)
        n_flagged: total number of flagged pairs
    """
    if min_count is None:
        if config is not None:
            min_count = config.get("qc_thresholds", {}).get("bead_count_min", BEAD_COUNT_MIN)
        else:
            min_count = BEAD_COUNT_MIN

    flagged = df[df["count"] < min_count][["well", "sample_name", "analyte", "count"]].copy()

    by_well = (
        df.groupby(["well", "sample_name", "well_type"])["count"]
        .median()
        .reset_index()
        .rename(columns={"count": "median_count"})
    )

    return {
        "flagged": flagged,
        "by_well": by_well,
        "n_flagged": len(flagged),
    }
