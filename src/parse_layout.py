"""Parse optional plate layout xlsx files."""

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def read_plate_layout(path: str | Path) -> pd.DataFrame | None:
    """Read the 'Sample list' sheet from a plate layout xlsx.

    Returns DataFrame with columns: well, sample_id, visit_date
    Returns None if file doesn't exist or can't be parsed.
    """
    path = Path(path)
    if not path.exists():
        return None

    try:
        df = pd.read_excel(path, sheet_name="Sample list", header=None)
    except Exception as exc:
        logger.warning("Could not parse plate layout '%s': %s", path, exc)
        return None

    # First row is the header
    if df.empty:
        return None
    headers = df.iloc[0].tolist()
    df = df.iloc[1:].reset_index(drop=True)
    df.columns = [str(h).strip() for h in headers]

    # Standardize column names
    col_map = {}
    for col in df.columns:
        lower = col.lower()
        if "well_number" in lower:
            col_map[col] = "well"
        elif lower == "well":
            col_map[col] = "well"
        elif "sample" in lower:
            col_map[col] = "sample_id"
        elif "date" in lower or "dt_visit" in lower:
            col_map[col] = "visit_date"
        elif "dilution" in lower:
            col_map[col] = "dilution"
    df = df.rename(columns=col_map)

    keep = [c for c in ["well", "sample_id", "visit_date", "dilution"] if c in df.columns]
    if "well" not in keep:
        return None

    return df[keep].reset_index(drop=True)
