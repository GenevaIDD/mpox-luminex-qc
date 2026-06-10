"""Sample label QC — flag specimen wells that share the same sample ID."""

import pandas as pd


def qc_duplicate_labels(df: pd.DataFrame) -> dict:
    """Find specimen wells that share the same sample name.

    Duplicate sample IDs usually indicate a labelling/transcription error
    (e.g. the same ID accidentally entered for two different wells).

    Returns dict with:
        has_duplicates: bool
        duplicates: list of {"sample_name": str, "wells": [str, ...], "n": int}
                    sorted by sample_name, wells sorted naturally
    """
    specimens = df[df["well_type"] == "specimen"][["well", "sample_name"]].drop_duplicates()

    counts = specimens.groupby("sample_name")["well"].nunique()
    dup_names = sorted(counts[counts > 1].index.tolist())

    duplicates = []
    for name in dup_names:
        wells = sorted(specimens.loc[specimens["sample_name"] == name, "well"].tolist())
        duplicates.append({"sample_name": name, "wells": wells, "n": len(wells)})

    return {
        "has_duplicates": len(duplicates) > 0,
        "duplicates": duplicates,
    }
