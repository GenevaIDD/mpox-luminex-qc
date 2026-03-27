"""Generate self-contained HTML QC report."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from jinja2 import Environment, FileSystemLoader

from .config import (
    APP_VERSION, MPXV_ANTIGENS, MPXV_KIT_CONTROLS, RECOVERY_TOLERANCE,
    NC_BEAD_MFI_MAX, SCG_MFI_MIN, FC_MFI_RANGE, IC_MFI_RANGE,
)


def generate_report(
    metadata: dict,
    data: pd.DataFrame,
    bead_qc: dict,
    fits: dict,
    replicate_qc: dict,
    nc_levels: pd.DataFrame,
    kit_controls: dict,
    specimen_results: pd.DataFrame,
    summary: dict,
    history_std: pd.DataFrame | None = None,
    history_nc: pd.DataFrame | None = None,
    output_path: str | Path = "report.html",
) -> Path:
    """Generate the full QC report as a self-contained HTML file."""
    output_path = Path(output_path)

    current_plate_id = metadata.get("plate_id", "")

    # Build all the plotly figures
    figures = {}
    figures["bead_heatmap"] = _make_bead_heatmap(bead_qc["by_well"])
    figures["standard_curves"] = _make_standard_curve_plots(fits, specimen_results, history_std)
    if replicate_qc["has_replicates"]:
        figures["replicate_cv"] = _make_replicate_plot(replicate_qc["replicate_cv"])
    figures["pc_mfi_history"] = _make_pc_mfi_history(history_std, current_plate_id)
    figures["nc_levels"] = _make_nc_plot(nc_levels, history_nc)
    figures["nc_history"] = _make_nc_history(history_nc, current_plate_id)
    figures["kit_controls"] = _make_kit_control_plots(kit_controls)
    if not specimen_results.empty:
        figures["specimen_mfi"] = _make_specimen_distribution(specimen_results)

    # Convert figures to JSON for embedding
    figure_json = {}
    for key, fig in figures.items():
        if fig is not None:
            figure_json[key] = json.loads(plotly.io.to_json(fig))

    # Render template (PyInstaller-aware path resolution)
    if getattr(sys, "frozen", False):
        template_dir = Path(sys._MEIPASS) / "templates"
    else:
        template_dir = Path(__file__).parent.parent / "templates"
    env = Environment(loader=FileSystemLoader(str(template_dir)), autoescape=True)
    template = env.get_template("report.html")

    recovery_lo = int((1.0 - RECOVERY_TOLERANCE) * 100)
    recovery_hi = int((1.0 + RECOVERY_TOLERANCE) * 100)

    html = template.render(
        metadata=metadata,
        summary=summary,
        bead_qc=bead_qc,
        fits=_fits_to_table(fits),
        replicate_qc=replicate_qc,
        nc_levels=nc_levels.to_dict("records") if not nc_levels.empty else [],
        kit_controls=_kit_controls_to_tables(kit_controls),
        kit_control_refs=_kit_control_reference_table(),
        specimen_results=_specimen_to_table(specimen_results),
        extrapolation=_extrapolation_summary(specimen_results),
        figures=figure_json,
        plotly_js=plotly.offline.get_plotlyjs(),
        version=APP_VERSION,
        recovery_range=f"{recovery_lo}–{recovery_hi}%",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    return output_path


def _make_bead_heatmap(by_well: pd.DataFrame) -> go.Figure:
    """Plate layout plot of median bead counts with marker shapes by well type."""
    by_well = by_well.copy()
    by_well["row"] = by_well["well"].str[0]
    by_well["col"] = by_well["well"].str[1:].astype(int)
    # Map row letters to numeric for plotting (A=0, B=1, ...)
    row_labels = list("ABCDEFGH")
    by_well["row_num"] = by_well["row"].map({r: i for i, r in enumerate(row_labels)})

    # Well type styling
    type_style = {
        "pc":       {"symbol": "circle",         "name": "PC (Standard)"},
        "nc":       {"symbol": "x",              "name": "NC (Negative)"},
        "specimen": {"symbol": "square",         "name": "Specimen"},
    }

    # Assign discrete colors: red < 30, yellow 30–50, green > 50
    def _bead_color(count):
        if count < 30:
            return "#e74c3c"   # red
        elif count <= 50:
            return "#f39c12"   # yellow/orange
        else:
            return "#27ae60"   # green

    by_well["marker_color"] = by_well["median_count"].apply(_bead_color)

    fig = go.Figure()

    for wtype, style in type_style.items():
        subset = by_well[by_well["well_type"] == wtype]
        if subset.empty:
            continue
        fig.add_trace(go.Scatter(
            x=subset["col"],
            y=subset["row_num"],
            mode="markers+text",
            marker=dict(
                symbol=style["symbol"],
                size=22,
                color=subset["marker_color"].tolist(),
                line=dict(width=1, color="grey"),
            ),
            text=subset["median_count"].round(0).astype(int).astype(str),
            textposition="middle center",
            textfont=dict(size=8),
            name=style["name"],
            hovertemplate=(
                "Well %{customdata[0]}: %{customdata[1]:.0f} beads"
                "<br>Type: " + style["name"] + "<extra></extra>"
            ),
            customdata=list(zip(subset["well"], subset["median_count"])),
        ))

    # Add invisible traces for the color legend
    for label, color in [("< 30 beads", "#e74c3c"), ("30–50 beads", "#f39c12"), ("> 50 beads", "#27ae60")]:
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode="markers",
            marker=dict(size=10, color=color),
            name=label, showlegend=True,
        ))

    fig.update_layout(
        title="Median Bead Count by Well",
        xaxis=dict(
            title="Column", tickmode="array",
            tickvals=list(range(1, 13)), ticktext=[str(i) for i in range(1, 13)],
            range=[0.5, 12.5],
        ),
        yaxis=dict(
            title="Row", tickmode="array",
            tickvals=list(range(8)), ticktext=row_labels,
            autorange="reversed", range=[-0.5, 7.5],
        ),
        height=400, margin=dict(t=40, b=40),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


def _make_standard_curve_plots(
    fits: dict, specimens: pd.DataFrame, history: pd.DataFrame | None
) -> go.Figure:
    """Standard curve plots for all 8 antigens in a 2x4 grid."""
    fig = make_subplots(rows=2, cols=4, subplot_titles=MPXV_ANTIGENS,
                        horizontal_spacing=0.06, vertical_spacing=0.12)

    for i, analyte in enumerate(MPXV_ANTIGENS):
        row = i // 4 + 1
        col = i % 4 + 1
        fit = fits.get(analyte, {})

        # Historical curves (grey)
        if history is not None and not history.empty:
            hist_analyte = history[history["analyte"] == analyte]
            for pid in hist_analyte["plate_id"].unique():
                hdata = hist_analyte[hist_analyte["plate_id"] == pid].sort_values("dilution")
                fig.add_trace(go.Scatter(
                    x=hdata["dilution"], y=hdata["mfi"],
                    mode="lines", line=dict(color="lightgrey", width=1),
                    showlegend=False, hoverinfo="skip",
                ), row=row, col=col)

        # Current plate observed points — colored by recovery status
        std_data = fit.get("std_data", pd.DataFrame())
        mean_data = fit.get("mean_data", pd.DataFrame())
        obs_exp = fit.get("obs_exp")

        # Build a set of dilutions that are out of tolerance
        bad_dilutions = set()
        if obs_exp:
            for pt in obs_exp:
                if not pt["in_range"]:
                    bad_dilutions.add(pt["dilution"])

        if not std_data.empty:
            # In-range points (blue circles)
            ok_mask = ~std_data["dilution"].isin(bad_dilutions)
            if ok_mask.any():
                ok_data = std_data[ok_mask]
                fig.add_trace(go.Scatter(
                    x=ok_data["dilution"], y=ok_data["mfi"],
                    mode="markers", marker=dict(color="blue", size=6),
                    name="Observed", showlegend=(i == 0),
                    hovertemplate="1:%{x} → MFI %{y:.0f}<extra></extra>",
                ), row=row, col=col)
            # Out-of-range points (red triangles)
            if (~ok_mask).any():
                bad_data = std_data[~ok_mask]
                fig.add_trace(go.Scatter(
                    x=bad_data["dilution"], y=bad_data["mfi"],
                    mode="markers", marker=dict(
                        color="red", size=8, symbol="triangle-up",
                        line=dict(width=1, color="darkred"),
                    ),
                    name="Out of tolerance", showlegend=(i == 0 and len(bad_dilutions) > 0),
                    hovertemplate="1:%{x} → MFI %{y:.0f} (out of tolerance)<extra></extra>",
                ), row=row, col=col)

        # Dropped outlier point (X marker, grey)
        dropped = fit.get("dropped_point")
        if dropped:
            fig.add_trace(go.Scatter(
                x=[dropped["dilution"]], y=[dropped["mfi"]],
                mode="markers", marker=dict(
                    color="grey", size=10, symbol="x",
                    line=dict(width=2, color="black"),
                ),
                name="Dropped outlier", showlegend=(i == 0),
                hovertemplate="DROPPED: 1:%{x} → MFI %{y:.0f}<extra></extra>",
            ), row=row, col=col)

        # Fitted curve
        params = fit.get("params")
        if params is not None:
            # Extend fit line across the full dilution range shown on the plot
            x_min = min(std_data["dilution"].min(), 30) if not std_data.empty else 30
            x_max = max(std_data["dilution"].max(), 3200) * 2
            x_fit = np.logspace(np.log10(x_min), np.log10(x_max), 200)
            y_fit = four_pl_np(x_fit, *params)
            fig.add_trace(go.Scatter(
                x=x_fit, y=y_fit,
                mode="lines", line=dict(color="red", width=2),
                name="4PL Fit", showlegend=(i == 0),
                hoverinfo="skip",
            ), row=row, col=col)

        # Reportable range box (LLOQ–ULOQ)
        rr = fit.get("reportable_range")
        if rr and params is not None and rr.get("lloq_dilution") and rr.get("uloq_dilution"):
            lloq_dil = rr["lloq_dilution"]
            uloq_dil = rr["uloq_dilution"]
            # Compute MFI at the LLOQ and ULOQ dilutions from the 4PL
            mfi_at_lloq = float(four_pl_np(np.array([lloq_dil]), *params)[0])
            mfi_at_uloq = float(four_pl_np(np.array([uloq_dil]), *params)[0])
            # Draw green shaded box — on log axes, use raw data values (Plotly logs them)
            # Subplot axis refs: first subplot is "x"/"y", rest are "x2"/"y2", "x3"/"y3", etc.
            ax_idx = i + 1
            xref = "x" if ax_idx == 1 else f"x{ax_idx}"
            yref = "y" if ax_idx == 1 else f"y{ax_idx}"
            fig.add_shape(
                type="rect",
                x0=uloq_dil, x1=lloq_dil,
                y0=min(mfi_at_lloq, mfi_at_uloq),
                y1=max(mfi_at_lloq, mfi_at_uloq),
                xref=xref, yref=yref,
                fillcolor="rgba(0, 180, 0, 0.08)",
                line=dict(color="rgba(0, 150, 0, 0.4)", width=1, dash="dot"),
            )

        # Specimen rug on y-axis
        if not specimens.empty:
            spec_mfi = specimens[specimens["analyte"] == analyte]["mfi"]
            if not spec_mfi.empty:
                fig.add_trace(go.Scatter(
                    x=[30] * len(spec_mfi), y=spec_mfi,
                    mode="markers", marker=dict(color="green", size=3, symbol="line-ew", line_width=1),
                    name="Specimens", showlegend=(i == 0),
                    hoverinfo="skip",
                ), row=row, col=col)

        fig.update_xaxes(type="log", autorange="reversed", row=row, col=col)
        fig.update_yaxes(type="log", row=row, col=col)

    # Axis labels — only show on edge subplots to avoid clutter
    for col in range(1, 5):
        fig.update_xaxes(title_text="Dilution", row=2, col=col)
    for row in range(1, 3):
        fig.update_yaxes(title_text="MFI", row=row, col=1)

    fig.update_layout(height=600, title_text="PC Standard Curves (4PL)", margin=dict(t=60))
    return fig


def four_pl_np(x, a, b, c, d):
    """Numpy-compatible 4PL for plotting."""
    return d + (a - d) / (1.0 + (x / c) ** b)


def _make_replicate_plot(cv_df: pd.DataFrame) -> go.Figure:
    """Bar chart of replicate CVs by analyte and dilution."""
    fig = go.Figure()
    for analyte in cv_df["analyte"].unique():
        adata = cv_df[cv_df["analyte"] == analyte].sort_values("dilution")
        colors = ["red" if f else "steelblue" for f in adata["flag"]]
        fig.add_trace(go.Bar(
            x=[f"1:{int(d)}" for d in adata["dilution"]],
            y=adata["cv"],
            name=analyte,
            marker_color=colors,
            visible=True if analyte == cv_df["analyte"].unique()[0] else "legendonly",
        ))
    fig.update_layout(
        title="PC Replicate CV by Dilution",
        xaxis_title="Dilution", yaxis_title="CV",
        yaxis=dict(tickformat=".0%"),
        height=400,
        barmode="group",
    )
    fig.add_hline(y=0.25, line_dash="dash", line_color="red", annotation_text="25% threshold")
    return fig


def _plate_sort_and_label(history: pd.DataFrame) -> list[tuple[str, str]]:
    """Sort plates by run_date and create short labels (date + PLATEXX).

    Returns list of (plate_id, label) tuples sorted by date.
    """
    import re
    from datetime import datetime

    plate_info = history.groupby("plate_id").first().reset_index()
    if "run_date" in plate_info.columns:
        # Parse run_date and sort
        def _parse_date(d):
            for fmt in ("%m/%d/%Y %I:%M %p", "%m/%d/%Y %I:%M:%S %p", "%Y-%m-%d"):
                try:
                    return datetime.strptime(str(d).strip(), fmt)
                except (ValueError, TypeError):
                    continue
            return datetime.min
        plate_info["_sort_date"] = plate_info["run_date"].apply(_parse_date)
        plate_info = plate_info.sort_values("_sort_date")

    result = []
    for _, row in plate_info.iterrows():
        pid = row["plate_id"]
        # Extract PLATEXX from plate_id
        m = re.search(r"(PLATE\s*\d+)", pid, re.IGNORECASE)
        plate_part = m.group(1) if m else pid[-8:]
        # Extract short date
        date_part = ""
        if "run_date" in row and row["run_date"]:
            try:
                dt = row["_sort_date"] if "_sort_date" in row else None
                if dt and dt != datetime.min:
                    date_part = dt.strftime("%m/%d")
            except Exception:
                pass
        label = f"{date_part} {plate_part}".strip() if date_part else plate_part
        result.append((pid, label))
    return result


def _make_pc_mfi_history(history: pd.DataFrame | None, current_plate_id: str) -> go.Figure | None:
    """Panel plot: PC MFI by plate for each analyte (2x4 grid)."""
    if history is None or history.empty:
        return None

    antigens = [a for a in MPXV_ANTIGENS if a in history["analyte"].unique()]
    if not antigens:
        return None

    fig = make_subplots(rows=2, cols=4, subplot_titles=antigens[:8],
                        horizontal_spacing=0.06, vertical_spacing=0.25)

    plate_order = _plate_sort_and_label(history)
    plates = [pid for pid, _ in plate_order]
    plate_labels = [label for _, label in plate_order]

    for i, analyte in enumerate(antigens[:8]):
        row = i // 4 + 1
        col = i % 4 + 1
        adata = history[history["analyte"] == analyte]

        for dil in sorted(adata["dilution"].unique()):
            dil_data = adata[adata["dilution"] == dil]
            plate_mfis = []
            labels = []
            colors = []
            for j, pid in enumerate(plates):
                pdata = dil_data[dil_data["plate_id"] == pid]
                if not pdata.empty:
                    plate_mfis.append(pdata["mfi"].mean())
                    labels.append(plate_labels[j])
                    colors.append("red" if pid == current_plate_id else "steelblue")

            fig.add_trace(go.Scatter(
                x=labels, y=plate_mfis,
                mode="markers+lines",
                marker=dict(size=6, color=colors),
                line=dict(color="lightgrey", width=1),
                name=f"1:{int(dil)}", showlegend=(i == 0),
                hovertemplate="%{x}: MFI %{y:.0f}<extra>1:" + str(int(dil)) + "</extra>",
            ), row=row, col=col)

        fig.update_yaxes(type="log", row=row, col=col)

    fig.update_layout(
        height=800,
        title_text="PC Standard MFI Across Plates (by dilution)",
        margin=dict(t=60, b=120),
    )
    fig.update_xaxes(tickangle=-90, tickfont=dict(size=10))
    return fig


def _make_nc_plot(nc_levels: pd.DataFrame, history: pd.DataFrame | None) -> go.Figure:
    """NC MFI by analyte — current plate bar chart."""
    if nc_levels.empty:
        fig = go.Figure()
        fig.update_layout(title="No NC wells on this plate")
        return fig

    mean_nc = nc_levels.groupby("analyte")["mfi"].mean().reset_index()
    fig = go.Figure(go.Bar(
        x=mean_nc["analyte"], y=mean_nc["mfi"],
        marker_color="steelblue",
        hovertemplate="%{x}: %{y:.1f}<extra></extra>",
    ))
    fig.update_layout(
        title="Negative Control MFI by Analyte",
        xaxis_title="Analyte", yaxis_title="MFI",
        height=350,
    )
    return fig


def _make_nc_history(history_nc: pd.DataFrame | None, current_plate_id: str) -> go.Figure | None:
    """Panel plot: NC MFI by plate for each analyte (2x4 grid)."""
    if history_nc is None or history_nc.empty:
        return None

    antigens = [a for a in MPXV_ANTIGENS if a in history_nc["analyte"].unique()]
    if not antigens:
        return None

    fig = make_subplots(rows=2, cols=4, subplot_titles=antigens[:8],
                        horizontal_spacing=0.06, vertical_spacing=0.25)

    plate_order = _plate_sort_and_label(history_nc)
    plates = [pid for pid, _ in plate_order]
    plate_labels = [label for _, label in plate_order]

    for i, analyte in enumerate(antigens[:8]):
        row = i // 4 + 1
        col = i % 4 + 1
        adata = history_nc[history_nc["analyte"] == analyte]

        mfis = []
        labels = []
        colors = []
        for j, pid in enumerate(plates):
            pdata = adata[adata["plate_id"] == pid]
            if not pdata.empty:
                mfis.append(pdata["mfi"].mean())
                labels.append(plate_labels[j])
                colors.append("red" if pid == current_plate_id else "steelblue")

        fig.add_trace(go.Scatter(
            x=labels, y=mfis,
            mode="markers+lines",
            marker=dict(size=8, color=colors),
            line=dict(color="lightgrey", width=1),
            showlegend=False,
            hovertemplate="%{x}: MFI %{y:.1f}<extra></extra>",
        ), row=row, col=col)

    fig.update_layout(
        height=800,
        title_text="Negative Control MFI Across Plates",
        margin=dict(t=60, b=120),
    )
    fig.update_xaxes(tickangle=-90, tickfont=dict(size=10))
    return fig


def _make_kit_control_plots(kit_controls: dict) -> go.Figure:
    """Combined plot for kit control beads."""
    fig = make_subplots(rows=1, cols=4,
                        subplot_titles=["NC Bead", "ScG", "FC", "IC"])

    for i, (key, title) in enumerate([
        ("nc_bead", "NC Bead"), ("scg", "ScG"), ("fc", "FC"), ("ic", "IC")
    ]):
        cdata = kit_controls[key]
        colors = ["red" if f else "steelblue" for f in cdata["flag"]]
        fig.add_trace(go.Bar(
            x=cdata["well"], y=cdata["mfi"],
            marker_color=colors,
            showlegend=False,
            hovertemplate="Well %{x}: %{y:.0f}<extra></extra>",
        ), row=1, col=i+1)

    fig.update_yaxes(title_text="MFI")
    fig.update_layout(height=300, title_text="Kit Control Beads", margin=dict(t=60))
    return fig


def _make_specimen_distribution(specimens: pd.DataFrame) -> go.Figure:
    """MFI distribution histograms for specimens, faceted by antigen."""
    antigens_present = [a for a in MPXV_ANTIGENS if a in specimens["analyte"].unique()]
    n = len(antigens_present)
    if n == 0:
        return None

    fig = make_subplots(rows=2, cols=4, subplot_titles=antigens_present[:8])

    for i, analyte in enumerate(antigens_present[:8]):
        row = i // 4 + 1
        col = i % 4 + 1
        vals = specimens[specimens["analyte"] == analyte]["mfi"].dropna()
        # Use log(MFI) for histogram; filter out non-positive values
        log_vals = np.log(vals[vals > 0])
        fig.add_trace(go.Histogram(
            x=log_vals, nbinsx=20, marker_color="steelblue",
            showlegend=False,
            hovertemplate="ln(MFI) %{x:.2f}: count %{y}<extra></extra>",
        ), row=row, col=col)

    # Axis labels on edge subplots
    for col in range(1, 5):
        fig.update_xaxes(title_text="ln(MFI)", row=2, col=col)
    for row in range(1, 3):
        fig.update_yaxes(title_text="Count", row=row, col=1)

    fig.update_layout(height=500, title_text="Specimen MFI Distribution (ln scale)", margin=dict(t=60))
    return fig


def _fits_to_table(fits: dict) -> list[dict]:
    """Convert fits dict to a list of dicts for the template."""
    rows = []
    for analyte in MPXV_ANTIGENS:
        f = fits.get(analyte, {})
        row = {"analyte": analyte, "fit_ok": f.get("fit_ok", False)}
        if f.get("params"):
            a, b, c, d = f["params"]
            row.update({"a": f"{a:.1f}", "b": f"{b:.3f}", "c": f"{c:.1f}", "d": f"{d:.1f}"})
        else:
            row.update({"a": "-", "b": "-", "c": "-", "d": "-"})
        row["error"] = f.get("error", "")
        row["qc_warnings"] = f.get("qc_warnings", [])
        row["obs_exp"] = f.get("obs_exp")
        row["reportable_range"] = f.get("reportable_range")
        row["dropped_point"] = f.get("dropped_point")
        rows.append(row)
    return rows


def _kit_control_reference_table() -> list[dict]:
    """Return reference info for kit control beads for display in report."""
    return [
        {
            "name": "NC (Negative Control)",
            "purpose": "Non-specific background; should be near-zero MFI",
            "criterion": f"MFI ≤ {NC_BEAD_MFI_MAX}",
        },
        {
            "name": "ScG (Sheep anti-Goat)",
            "purpose": "Secondary antibody binding control; confirms conjugate is working",
            "criterion": f"MFI ≥ {SCG_MFI_MIN:,}",
        },
        {
            "name": "FC (Fluorescent Conjugate)",
            "purpose": "Fluorescent reference bead; checks detector calibration",
            "criterion": f"MFI {FC_MFI_RANGE[0]:,} – {FC_MFI_RANGE[1]:,}",
        },
        {
            "name": "IC (Instrument Control)",
            "purpose": "Internal bead consistency control; should be stable across plates",
            "criterion": f"MFI {IC_MFI_RANGE[0]:,} – {IC_MFI_RANGE[1]:,}",
        },
    ]


def _kit_controls_to_tables(kit_controls: dict) -> dict:
    """Convert kit control DataFrames to lists of dicts for the template."""
    result = {}
    for key in ["nc_bead", "scg", "fc", "ic"]:
        df = kit_controls[key]
        result[key] = df.to_dict("records")
    result["n_flagged"] = len(kit_controls["flagged_wells"])
    result["flagged_wells"] = sorted(kit_controls["flagged_wells"])
    return result


def _specimen_to_table(specimens: pd.DataFrame) -> list[dict]:
    """Pivot specimen results to wide format for the table."""
    if specimens.empty:
        return []
    # Pivot: one row per well, columns for each antigen's MFI and RAU
    antigens = [a for a in MPXV_ANTIGENS if a in specimens["analyte"].unique()]

    rows = []
    for well in specimens["well"].unique():
        wdata = specimens[specimens["well"] == well]
        sample_name = wdata["sample_name"].iloc[0]
        row = {"well": well, "sample_name": sample_name}
        for analyte in antigens:
            arow = wdata[wdata["analyte"] == analyte]
            if not arow.empty:
                row[f"{analyte}_mfi"] = int(round(arow["mfi"].iloc[0]))
                rau = arow["rau"].iloc[0]
                below_lloq = bool(arow["below_lloq"].iloc[0]) if "below_lloq" in arow.columns else False
                above_uloq = bool(arow["above_uloq"].iloc[0]) if "above_uloq" in arow.columns else False
                extrap = bool(arow["extrapolated"].iloc[0]) if "extrapolated" in arow.columns else False
                # au_censored: none / left (below LLOQ) / right (above ULOQ)
                if below_lloq:
                    censored = "left"
                elif above_uloq:
                    censored = "right"
                else:
                    censored = "none"
                rau_str = None
                if pd.notna(rau):
                    if extrap and not below_lloq and not above_uloq:
                        rau_str = f"{rau:.1f}*"  # extrapolated but within range
                    else:
                        rau_str = f"{rau:.1f}"
                row[f"{analyte}_AU"] = rau_str
                row[f"{analyte}_censored"] = censored
                net_mfi = arow["net_mfi"].iloc[0] if "net_mfi" in arow.columns else None
                row[f"{analyte}_net_mfi"] = int(round(net_mfi)) if pd.notna(net_mfi) and net_mfi is not None else None
        rows.append(row)
    return rows


def _extrapolation_summary(specimens: pd.DataFrame) -> dict:
    """Summarise how many specimen-analyte pairs are extrapolated."""
    if specimens.empty or "extrapolated" not in specimens.columns:
        return {"n_extrapolated": 0, "details": []}
    extrap = specimens[specimens["extrapolated"]].copy()
    n = len(extrap)
    # Summarise by analyte
    details = []
    if n > 0:
        for analyte in MPXV_ANTIGENS:
            aext = extrap[extrap["analyte"] == analyte]
            if not aext.empty:
                above = (aext["mfi"] > aext["mfi"].median()).sum()  # rough
                details.append({"analyte": analyte, "count": len(aext)})
    return {"n_extrapolated": n, "details": details}
