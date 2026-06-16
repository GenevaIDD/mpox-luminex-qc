"""One-off script to anonymize real xPONENT CSVs into shareable example files.

Replaces instrument/operator metadata and specimen sample IDs (e.g. "1976 S1")
with generic placeholders, while leaving control well labels (VIG, INRB PC, NC)
and all numeric results untouched. Run once per source file; not part of the
regular pipeline.
"""

import csv
import io
import re
import sys
from pathlib import Path

SPECIMEN_RE = re.compile(r"^(\d{3,4}) (S\d+)$")


def anonymize(in_path: Path, out_path: Path, batch_name: str):
    text = in_path.read_text(encoding="utf-8-sig")
    reader = csv.reader(io.StringIO(text))
    rows = list(reader)

    # Build specimen ID mapping (in order of first appearance)
    id_map = {}
    next_id = 1
    for row in rows:
        if len(row) >= 2:
            m = SPECIMEN_RE.match(row[1])
            if m:
                prefix = m.group(1)
                if prefix not in id_map:
                    id_map[prefix] = f"EX{next_id:03d}"
                    next_id += 1

    out_rows = []
    for row in rows:
        if not row:
            out_rows.append(row)
            continue
        key = row[0]
        if key == "Date" and len(row) >= 3:
            row = ["Date", "01/01/2026", "12:00 PM"]
        elif key == "SN":
            row = ["SN", "MAGPX00000000"]
        elif key == "Batch":
            row = ["Batch", batch_name]
        elif key == "Operator":
            row = ["Operator", "EXAMPLE"]
        elif key == "ComputerName":
            row = ["ComputerName", "DESKTOP-EXAMPLE"]
        elif key in ("BatchStartTime", "BatchStopTime"):
            row = [key, "1/1/2026 1:00:00 PM"]
        elif key in ("Last CAL Calibration", "Last VER Verification", "Last Fluidics Test"):
            label, rest = row[1].split(" ", 1)
            row = [key, f"{label} 12/31/2025 10:00:00"]
        elif key == "C09108":
            row = list(row)
            row[0] = "C00000"
            row[1] = "12/31/2025"
            row[2] = "12/31/2025 10:00:00 AM"
            row[-1] = "MAGPX00000000"
        else:
            row = list(row)
            for idx in (0, 1):
                if idx < len(row):
                    m = SPECIMEN_RE.match(row[idx])
                    if m:
                        row[idx] = f"{id_map[m.group(1)]} {m.group(2)}"
        out_rows.append(row)

    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, quoting=csv.QUOTE_ALL)
        for row in out_rows:
            writer.writerow(row)


if __name__ == "__main__":
    base = Path("/Users/andrew/Library/CloudStorage/Dropbox/INRB-UNIGE-Serology/Raw-Data")
    out = Path("examples")
    out.mkdir(exist_ok=True)

    anonymize(
        base / "A020626-MP1822-GK02-P09-12PLXMPXVSERUMTS.csv",
        out / "EXAMPLE01-MP0000-GK00-P01-12PLXMPXVSERUMTS.csv",
        "EXAMPLE01-MP0000-GK00-P01-12PLXMPXVSERUMTS",
    )
    anonymize(
        base / "A040626-MP1822-GK02-P10-12PLXMPXVSERUMTS.csv",
        out / "EXAMPLE02-MP0000-GK00-P02-12PLXMPXVSERUMTS.csv",
        "EXAMPLE02-MP0000-GK00-P02-12PLXMPXVSERUMTS",
    )
    print("done")
