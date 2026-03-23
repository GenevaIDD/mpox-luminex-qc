# MPXV Luminex QC Tool — Specification
# Version 0.1 (2026-03-23)

## Overview

Standalone quality control tool for 12-plex MPXV (monkeypox virus) Luminex immunoassays from Tetracore run on a MagPix instrument. Processes xPONENT CSV exports, performs automated QC checks, fits standard curves, computes Relative Antibody Units (RAU), and generates interactive HTML reports. Distributed as a macOS `.app` (or Windows `.exe`) requiring no installation.

## Assay Panel

### Antigens (8 beads)

| Bead Region | Analyte    | Description                     |
|-------------|------------|---------------------------------|
| 01          | MVA Ag     | Modified Vaccinia Ankara antigen |
| 02          | VACV A33R  | Vaccinia virus A33R protein      |
| 03          | MPXV A35R  | Monkeypox virus A35R protein     |
| 04          | MPXV B6R   | Monkeypox virus B6R protein      |
| 05          | MPXV A27   | Monkeypox virus A27 protein      |
| 06          | MPXV E8L   | Monkeypox virus E8L protein      |
| 07          | MPXV H3L   | Monkeypox virus H3L protein      |
| 08          | MPXV M1R   | Monkeypox virus M1R protein      |

### Kit Control Beads (4 beads)

| Bead Region | Name | Purpose                                                    |
|-------------|------|------------------------------------------------------------|
| 09          | NC   | Negative control bead — monitors non-specific binding       |
| 10          | ScG  | Sample control (goat) — confirms sample was added to well   |
| 11          | FC   | Fluorescent conjugate control — confirms PE reporter added  |
| 12          | IC   | Instrument control — monitors instrument performance        |

## Plate Layout

Standard 96-well plate:

- **Columns 1-2**: PC (positive control) standard curve in duplicate (7-point, 2-fold serial dilution: 1:50, 1:100, 1:200, 1:400, 1:800, 1:1600, 1:3200)
- **Row H, columns 1-2**: NC (negative control) wells
- **Remaining wells**: Specimen wells (single dilution, no replicates)

Variations handled:
- Single PC replicate (no duplicate in column 2)
- Missing NC wells
- Variable number of specimen wells (up to 80 per plate)

## Input Files

### xPONENT CSV (required)

Multi-block CSV exported from the MagPix xPONENT software. Contains:
- Header metadata (batch ID, run date, operator, protocol, instrument serial number)
- `Median` data block: median fluorescent intensity (MFI) per well per analyte
- `Count` data block: bead count per well per analyte

Filename convention: `A{date}-{batch}-{plate}-{panel}_{timestamp}.csv`

### Plate Layout XLSX (optional)

Excel workbook with a "Sample list" sheet containing columns:
- `well_number`: sequential well number
- `well`: well position (e.g., A1, B2)
- `sample_id`: sample identifier
- `date` or `dt_visit`: visit date

When provided, enriches specimen results with sample IDs and visit dates.

## QC Checks

### 1. Bead Counts

- **Threshold**: minimum 30 beads per well-analyte pair
- **Output**: count of flagged pairs, flagged details table, plate heatmap of median bead counts per well
- **Interpretation**: low bead counts indicate aspiration issues or bead loss; affected results are unreliable

### 2. Standard Curve Fitting (4-Parameter Logistic)

Model: `y = d + (a - d) / (1 + (x / c)^b)`

| Parameter | Meaning                                |
|-----------|----------------------------------------|
| a         | Minimum asymptote (high dilution → low MFI) |
| b         | Hill slope                              |
| c         | Inflection point (EC50)                 |
| d         | Maximum asymptote (low dilution → high MFI) |

- Fit independently for each of the 8 antigens
- Uses mean of PC duplicates at each dilution point
- Fitted via `scipy.optimize.curve_fit` with bounds
- Reports fit success/failure and all 4 parameters per analyte

### 3. PC Replicate Variability

- **Metric**: Coefficient of Variation (CV) between duplicate PC wells at each dilution
- **Threshold**: CV ≤ 25%
- **Skipped**: when only a single PC replicate is present

### 4. Negative Control Monitoring

- Reports mean MFI per analyte across NC wells
- Tracked historically across plates for drift detection

### 5. Kit Control Bead QC

| Control | Check                          | Threshold / Method      |
|---------|-------------------------------|------------------------|
| NC bead | Non-specific binding          | MFI ≤ 200              |
| ScG     | Sample addition verification   | Not low outlier (median − 3×MAD) |
| FC      | Conjugate addition verification | Not low outlier (median − 3×MAD) |
| IC      | Instrument consistency         | Not outlier (median ± 3×MAD) |

## Computed Outputs

### Relative Antibody Units (RAU)

For each specimen well and antigen:
1. Invert the fitted 4PL standard curve: given MFI, solve for the equivalent dilution factor
2. Compute RAU = 1 / dilution_factor

**Interpretation**: higher RAU = more antibody. A specimen with RAU = 0.02 produces the same signal as a 1:50 dilution of the positive control standard.

RAU is `NaN` when:
- The 4PL fit failed for that analyte
- The specimen MFI falls outside the standard curve range

## Reports

### Per-Plate HTML QC Report

Self-contained HTML file with embedded interactive Plotly charts. Sections:

1. **Plate Overview** — metadata table, well count summary (total, PC, NC, specimen)
2. **Bead Counts** — flagged pairs count, flagged details table, plate heatmap
3. **Standard Curves** — 4PL parameter table, 2×4 grid of curve plots with historical overlays (grey) and specimen rug marks
4. **PC Replicate CV** — flagged pairs, CV bar chart by dilution (conditional: only when duplicates present)
5. **Negative Control Levels** — NC MFI bar chart by analyte
6. **Kit Control Beads** — flagged wells, 1×4 bar chart per control bead
7. **Specimen Results** — wide-format table with MFI and RAU per antigen

### Per-Plate Specimen CSV

Long-format CSV with columns: `well, sample_name, analyte, mfi, count, well_type, dilution, rau`

### Export All Data (Excel workbook)

Combined data across all processed plates. Sheets:

| Sheet                  | Contents                                          |
|------------------------|--------------------------------------------------|
| `specimens`            | All specimen results with `plate_id` column       |
| `standard_curve_params`| 4PL parameters (a, b, c, d) per plate per analyte |
| `standard_curve_data`  | Raw standard curve MFI data points                |
| `nc_levels`            | Negative control MFI per plate per analyte        |

## Historical Tracking

JSON files accumulate data across plates for trend monitoring:

- `std_curve_history.json` — standard curve data points (plate_id, run_date, analyte, dilution, mfi)
- `nc_history.json` — NC MFI per analyte per plate
- `fit_history.json` — 4PL fit coefficients per plate per analyte

History is deduplicated on plate_id — reprocessing a plate overwrites its previous history entry. Historical standard curves appear as grey lines behind the current plate's curves in the report.

## Application Architecture

### Web Interface (Flask)

- Local-only web server (binds to `127.0.0.1` on a random free port)
- Browser auto-opens on launch
- Upload form: one or more xPONENT CSVs + optional layout XLSX
- Past reports table with view/download links
- "Export All Data" button for combined Excel workbook
- "Quit Application" button for graceful shutdown

### Data Storage

All persistent data stored in `~/mpox-luminex-qc-results/`:

```
~/mpox-luminex-qc-results/
  reports/        # Generated HTML QC reports
  specimens/      # Per-plate specimen CSV files
  history/        # JSON history files (accumulate across plates)
  uploads/        # Temporary (cleaned after each pipeline run)
```

### Distribution

- **macOS**: `.app` bundle via PyInstaller (~128 MB)
- **Windows**: folder with `.exe` via PyInstaller (build on Windows machine)
- No Python installation required on target machine
- No internet connection required (all assets embedded)

## Technology Stack

| Component     | Technology              | Purpose                        |
|---------------|------------------------|--------------------------------|
| Language      | Python 3.11+           | Core logic                     |
| Data          | pandas ≥ 2.0           | Data manipulation              |
| Curve fitting | scipy ≥ 1.11           | 4PL fitting (curve_fit)        |
| Visualization | plotly ≥ 5.18          | Interactive HTML charts        |
| Templating    | Jinja2 ≥ 3.1           | HTML report generation         |
| Excel I/O     | openpyxl ≥ 3.1         | Layout reading, Excel export   |
| Web server    | Flask ≥ 3.0            | Local web UI                   |
| Packaging     | PyInstaller ≥ 6.0      | Standalone executable          |
