"""Parse xPONENT CSV files from MagPix Luminex instrument."""

import csv
import re
from pathlib import Path

import pandas as pd


def parse_xponent_csv(path: str | Path) -> dict:
    """Parse an xPONENT CSV file into metadata and data DataFrames.

    Returns dict with keys:
        metadata: dict of plate-level metadata
        mfi: DataFrame (long format) with columns [well, sample_name, analyte, mfi]
        counts: DataFrame (long format) with columns [well, sample_name, analyte, count]
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    rows = list(csv.reader(lines))

    metadata = _parse_metadata(rows)
    metadata["file"] = path.name

    block_starts = _find_datatype_blocks(rows)

    mfi_wide = _parse_data_block(rows, block_starts["Median"])
    counts_wide = _parse_data_block(rows, block_starts["Count"])

    analytes = [c for c in mfi_wide.columns if c not in ("well", "sample_name", "total_events")]

    mfi_long = _wide_to_long(mfi_wide, analytes, "mfi")
    counts_long = _wide_to_long(counts_wide, analytes, "count")

    merged = mfi_long.merge(counts_long, on=["well", "sample_name", "analyte"])

    return {"metadata": metadata, "data": merged}


def _parse_metadata(rows: list[list[str]]) -> dict:
    """Extract plate metadata from header rows."""
    meta = {}
    field_map = {
        "Batch": "batch",
        "Date": "run_date",
        "Operator": "operator",
        "ProtocolName": "protocol",
        "SN": "instrument_sn",
        "ProtocolDescription": "protocol_description",
        "BatchStartTime": "batch_start_time",
    }
    for row in rows:
        if not row or not row[0].strip('"'):
            continue
        key = row[0].strip('"')
        if key in field_map and len(row) > 1:
            meta[field_map[key]] = row[1].strip('"')
        if key == "Date" and len(row) > 2:
            meta["run_date"] = f"{row[1].strip('\"')} {row[2].strip('\"')}"
        if key == "Samples" and len(row) > 1:
            meta["n_samples"] = int(row[1].strip('"'))
            break
    meta["plate_id"] = _extract_plate_id(meta.get("batch", ""))
    return meta


def _extract_plate_id(batch: str) -> str:
    """Extract a short plate ID from the batch string."""
    # e.g. "A260323-MP1822-KV02-Plate01-12PlxMPXVHIg" -> "A260323-MP1822-KV02-Plate01"
    m = re.match(r"(.+-Plate\d+)", batch)
    if m:
        return m.group(1)
    return batch


def _find_datatype_blocks(rows: list[list[str]]) -> dict[str, int]:
    """Find the starting row index for each DataType block."""
    blocks = {}
    for i, row in enumerate(rows):
        if len(row) >= 2 and row[0].strip('"') == "DataType:":
            dtype = row[1].strip('"')
            blocks[dtype] = i
    return blocks


def _parse_data_block(rows: list[list[str]], block_start: int) -> pd.DataFrame:
    """Parse a well-level data block (Median, Net MFI, or Count) into a wide DataFrame."""
    # Header row is the line after "DataType:,..."
    header_row = rows[block_start + 1]
    # Columns: Location, Sample, <analyte1>, <analyte2>, ..., Total Events
    col_names = [c.strip('"') for c in header_row]

    data_rows = []
    for i in range(block_start + 2, len(rows)):
        row = rows[i]
        if not row or not row[0].strip('"'):
            break
        data_rows.append([c.strip('"') for c in row])

    df = pd.DataFrame(data_rows, columns=col_names)

    # Parse well from Location: "1(1,A1)" -> "A1"
    df["well"] = df["Location"].apply(_parse_well_from_location)
    df = df.rename(columns={"Sample": "sample_name"})

    # Convert numeric columns
    analyte_cols = [c for c in col_names if c not in ("Location", "Sample")]
    for col in analyte_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Normalize column name
    if "Total Events" in df.columns:
        df = df.rename(columns={"Total Events": "total_events"})

    df = df.drop(columns=["Location"])
    return df


def _parse_well_from_location(loc: str) -> str:
    """Extract well ID from location string: '1(1,A1)' -> 'A1'."""
    m = re.search(r"\d+\(\d+,([A-H]\d+)\)", loc)
    if m:
        return m.group(1)
    return loc


def _strip_bead_prefix(name: str) -> str:
    """Strip numeric prefix from bead names: '01 MVA Ag' -> 'MVA Ag'."""
    m = re.match(r"^\d+\s+(.+)$", name)
    return m.group(1) if m else name


def _wide_to_long(df: pd.DataFrame, analyte_cols: list[str], value_name: str) -> pd.DataFrame:
    """Pivot wide data block to long format."""
    id_cols = ["well", "sample_name"]
    long = df[id_cols + analyte_cols].melt(
        id_vars=id_cols, var_name="analyte", value_name=value_name
    )
    long["analyte"] = long["analyte"].apply(_strip_bead_prefix)
    return long
