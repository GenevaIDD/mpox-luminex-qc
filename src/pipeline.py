"""Main QC pipeline — orchestrates parsing, QC, and report generation."""

from pathlib import Path

import pandas as pd

from .parse_xponent import parse_xponent_csv
from .classify import classify_wells
from .qc_beads import qc_bead_counts
from .qc_standard_curve import fit_standard_curves, compute_concentrations, compute_net_mfi
from .qc_replicates import qc_pc_replicates
from .qc_nc import qc_nc_levels
from .qc_kit_controls import qc_kit_controls
from .qc_history import load_history, append_history, save_history
from .parse_layout import read_plate_layout
from .plate_summary import plate_summary
from .report import generate_report
from .settings import load_config


def run_pipeline(
    csv_path: str | Path,
    output_dir: str | Path | None = None,
    layout_path: str | Path | None = None,
    history_dir: str | Path | None = None,
    config: dict | None = None,
) -> Path:
    """Run the full QC pipeline on a single plate CSV.

    Args:
        csv_path: path to xPONENT CSV file
        output_dir: where to write report and CSVs (defaults to csv parent dir)
        layout_path: optional path to plate layout xlsx
        history_dir: where to store history JSON files (defaults to output_dir/history)
        config: optional config dict (from settings.load_config); loaded if not provided

    Returns path to the generated HTML report.
    """
    csv_path = Path(csv_path)
    if output_dir is None:
        output_dir = csv_path.parent
    output_dir = Path(output_dir)
    if history_dir is None:
        history_dir = output_dir / "history"

    # Load config if not provided
    if config is None:
        config = load_config()

    # 1. Parse CSV
    parsed = parse_xponent_csv(csv_path)
    metadata = parsed["metadata"]
    data = parsed["data"]

    # 2. Classify wells (using config patterns)
    data = classify_wells(data, config=config)

    # 3. Optional layout enrichment
    if layout_path:
        layout = read_plate_layout(layout_path)
        if layout is not None and "sample_id" in layout.columns:
            data = data.merge(layout, on="well", how="left", suffixes=("", "_layout"))

    # 4. QC: bead counts
    bead_qc = qc_bead_counts(data, config=config)

    # 5. QC: 4PL standard curves
    fits = fit_standard_curves(data, config=config)

    # 6. QC: PC replicate CV
    replicate_qc = qc_pc_replicates(data, config=config)

    # 7. QC: NC levels
    nc_levels = qc_nc_levels(data, config=config)

    # 8. QC: kit controls
    kit_ctrl = qc_kit_controls(data, config=config)

    # 9. Compute specimen AU and Net MFI
    specimen_results = compute_concentrations(data, fits)
    data = compute_net_mfi(data)
    # Merge net_mfi back into specimen_results
    net_mfi_map = data[data["well_type"] == "specimen"][["well", "analyte", "net_mfi"]]
    specimen_results = specimen_results.merge(net_mfi_map, on=["well", "analyte"], how="left")

    # 10. Plate summary
    summary = plate_summary(data)

    # 11. History — load, append, save
    std_history_path = Path(history_dir) / "std_curve_history.json"
    nc_history_path = Path(history_dir) / "nc_history.json"
    fit_history_path = Path(history_dir) / "fit_history.json"

    # Standard curve history
    history_std = load_history(std_history_path)
    new_std = _build_std_history(metadata, fits)
    if not new_std.empty:
        history_std = append_history(history_std, new_std, ["plate_id", "analyte", "dilution"])
        save_history(history_std, std_history_path)

    # NC history
    history_nc = load_history(nc_history_path)
    new_nc = _build_nc_history(metadata, nc_levels)
    if not new_nc.empty:
        history_nc = append_history(history_nc, new_nc, ["plate_id", "analyte"])
        save_history(history_nc, nc_history_path)

    # Fit history
    history_fit = load_history(fit_history_path)
    new_fit = _build_fit_history(metadata, fits)
    if not new_fit.empty:
        history_fit = append_history(history_fit, new_fit, ["plate_id", "analyte"])
        save_history(history_fit, fit_history_path)

    # 12. Generate report
    report_name = f"QC_{metadata['plate_id']}.html"
    report_path = output_dir / report_name

    generate_report(
        metadata=metadata,
        data=data,
        bead_qc=bead_qc,
        fits=fits,
        replicate_qc=replicate_qc,
        nc_levels=nc_levels,
        kit_controls=kit_ctrl,
        specimen_results=specimen_results,
        summary=summary,
        history_std=history_std,
        history_nc=history_nc,
        output_path=report_path,
    )

    # 13. Export specimen results CSV
    if not specimen_results.empty:
        csv_out = output_dir / f"specimens_{metadata['plate_id']}.csv"
        export_df = specimen_results.rename(columns={"rau": "AU"})
        # Add au_censored: none / left (below LLOQ) / right (above ULOQ)
        if "below_lloq" in export_df.columns and "above_uloq" in export_df.columns:
            export_df["au_censored"] = "none"
            export_df.loc[export_df["below_lloq"], "au_censored"] = "left"
            export_df.loc[export_df["above_uloq"], "au_censored"] = "right"
        export_df.to_csv(csv_out, index=False, encoding="utf-8")

    return report_path


def _build_std_history(metadata: dict, fits: dict) -> pd.DataFrame:
    """Build standard curve history entries from current plate fits."""
    rows = []
    for analyte, fit in fits.items():
        std_data = fit.get("std_data", pd.DataFrame())
        if std_data.empty:
            continue
        for _, r in std_data.iterrows():
            rows.append({
                "plate_id": metadata["plate_id"],
                "run_date": metadata.get("run_date", ""),
                "analyte": analyte,
                "dilution": r["dilution"],
                "mfi": r["mfi"],
            })
    return pd.DataFrame(rows)


def _build_nc_history(metadata: dict, nc_levels: pd.DataFrame) -> pd.DataFrame:
    """Build NC history entries."""
    if nc_levels.empty:
        return pd.DataFrame()
    nc_mean = nc_levels.groupby("analyte")["mfi"].mean().reset_index()
    nc_mean["plate_id"] = metadata["plate_id"]
    nc_mean["run_date"] = metadata.get("run_date", "")
    return nc_mean


def _build_fit_history(metadata: dict, fits: dict) -> pd.DataFrame:
    """Build fit coefficient history entries."""
    rows = []
    for analyte, fit in fits.items():
        row = {
            "plate_id": metadata["plate_id"],
            "run_date": metadata.get("run_date", ""),
            "analyte": analyte,
            "fit_ok": fit["fit_ok"],
        }
        if fit["params"]:
            a, b, c, d = fit["params"]
            row.update({"a": a, "b": b, "c": c, "d": d})
        rows.append(row)
    return pd.DataFrame(rows)
