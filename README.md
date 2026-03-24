# MPXV Luminex QC Tool

<p align="center">
  <img src="gdd_antibody_icon_bottom_right.png" alt="MPXV Luminex QC" width="128">
</p>

<p align="center">
  <strong>Automated quality control for 12-plex MPXV Luminex immunoassays</strong><br>
  <em>Geneva Centre for Emerging Viral Diseases &middot; University of Geneva</em>
</p>

<p align="center">
  <a href="SPECIFICATION.md">Full Specification</a>
</p>

---

## Overview

Standalone desktop tool for QC of Tetracore 12-plex monkeypox (MPXV) Luminex immunoassays run on a MagPix instrument. Upload xPONENT CSV exports and get interactive HTML reports with:

- **4-Parameter Logistic (4PL) standard curve fitting** with automated quality checks
- **1/RAU quantification** (inverse Relative Antibody Units) for all specimen wells
- **Bead count validation** (minimum 30 beads per well-analyte)
- **PC replicate variability** monitoring (CV threshold)
- **Negative control tracking** across plates
- **Kit control bead QC** (NC, ScG, FC, IC)
- **Historical trend plots** for standard curves, PC MFI, and NC levels across plates
- **Excel export** of all data across plates

No Python installation or internet connection required &mdash; runs as a self-contained macOS `.app` or Windows `.exe`.

## Assay Panel

| Bead | Analyte | | Bead | Control |
|------|---------|---|------|---------|
| 01 | MVA Ag | | 09 | NC (negative control) |
| 02 | VACV A33R | | 10 | ScG (sample control) |
| 03 | MPXV A35R | | 11 | FC (fluorescent conjugate) |
| 04 | MPXV B6R | | 12 | IC (instrument control) |
| 05 | MPXV A27 | | | |
| 06 | MPXV E8L | | | |
| 07 | MPXV H3L | | | |
| 08 | MPXV M1R | | | |

## Quick Start

### Download

Download the latest release for your platform from the [Releases](https://github.com/GenevaIDD/mpox-luminex-qc/releases) page:

- **macOS**: `MPXV Luminex QC.app` (~128 MB)
- **Windows**: `MPXV-Luminex-QC/` folder with `.exe`

### Run

1. **macOS**: Double-click `MPXV Luminex QC.app`. If macOS blocks it, right-click > Open, or run in Terminal:
   ```bash
   xattr -cr "MPXV Luminex QC.app"
   open "MPXV Luminex QC.app"
   ```

2. **Windows**: Double-click `MPXV-Luminex-QC.exe` inside the extracted folder.

3. Your browser opens automatically to the upload page.

4. Upload one or more xPONENT CSV files (and optionally a plate layout XLSX), then click **Generate Report**.

### Input Files

| File | Required | Description |
|------|----------|-------------|
| xPONENT CSV | Yes | MagPix instrument export with Median and Count data blocks |
| Plate layout XLSX | No | Maps well positions to sample IDs and visit dates |

## QC Checks

### Standard Curve Fit Quality

Each of the 8 antigen standard curves is fit independently using a 4PL model:

```
y = d + (a - d) / (1 + (x / c)^b)
```

| Check | Criterion | Rationale |
|-------|-----------|-----------|
| R² | >= 0.95 | Goodness of fit |
| IC50 | Within tested dilution range (±3x) | Inflection point plausibility |
| Hill slope | 0.3 -- 5.0 | Prevents flat or step-function fits |
| Dynamic range | >= 3-fold | Adequate signal separation |

### Other QC Checks

- **Bead counts**: Flag well-analyte pairs with < 30 beads
- **PC replicate CV**: Flag duplicates with CV > 25%
- **Kit control beads**: NC < 150, ScG > 10000, FC 2000--5000, IC 2000--3000

## Output

### 1/RAU (Inverse Relative Antibody Units)

For each specimen, the 4PL standard curve is inverted to find the equivalent dilution factor, then:

```
1/RAU = 1 / dilution_factor
```

Higher 1/RAU = more antibody. Values are flagged when extrapolated beyond the observed standard curve range.

### Reports & Exports

| Output | Format | Description |
|--------|--------|-------------|
| QC Report | HTML | Interactive per-plate report with Plotly charts |
| Specimen data | CSV | Per-plate long-format results |
| All data export | XLSX | Combined workbook across all plates (specimens, fit params, NC levels) |

All outputs are stored in `~/mpox-luminex-qc-results/`.

## Development

### Prerequisites

- Python 3.11+
- [uv](https://github.com/astral-sh/uv) (recommended) or pip

### Setup

```bash
git clone https://github.com/GenevaIDD/mpox-luminex-qc.git
cd mpox-luminex-qc
uv sync
```

### Run in dev mode

```bash
uv run python -m src.main
```

### Build standalone app

**macOS:**
```bash
uv run python -m PyInstaller mpox-luminex-qc.spec --clean -y
codesign --force --deep --sign - "dist/MPXV Luminex QC.app"
```

**Windows:**
```bash
python -m PyInstaller mpox-luminex-qc-win.spec --clean -y
```

### Regenerate app icons

```bash
uv run python scripts/make_icon.py
```

## Tech Stack

| Component | Technology |
|-----------|-----------|
| Data | pandas |
| Curve fitting | scipy (curve_fit) |
| Visualization | plotly |
| Web UI | Flask |
| Reports | Jinja2 |
| Excel I/O | openpyxl |
| Packaging | PyInstaller |

## Contact

**Andrew Azman** &mdash; [andrew.azman@unige.ch](mailto:andrew.azman@unige.ch)

Geneva Centre for Emerging Viral Diseases, University of Geneva

## License

This project is developed for internal use by the Geneva Centre for Emerging Viral Diseases.
