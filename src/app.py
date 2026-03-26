"""Flask web app for MPXV Luminex QC tool."""

import io
import json
import os
import signal
import sys
import traceback

import yaml
from datetime import datetime
from pathlib import Path

import pandas as pd
from flask import (
    Flask,
    flash,
    redirect,
    render_template,
    request,
    send_file,
    url_for,
)
from werkzeug.utils import secure_filename

from .config import APP_VERSION
from .pipeline import run_pipeline
from .settings import load_config, save_config, reset_config, get_config_path

# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def _get_base_path() -> Path:
    """Return the base path for bundled resources.

    PyInstaller sets sys._MEIPASS when running from a bundle.
    In dev, use the project root (parent of src/).
    """
    if getattr(sys, "frozen", False):
        return Path(sys._MEIPASS)
    return Path(__file__).parent.parent


def _get_results_dir() -> Path:
    """Persistent results directory in user's home."""
    d = Path.home() / "mpox-luminex-qc-results"
    for sub in ("reports", "specimens", "history", "uploads"):
        (d / sub).mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# App factory
# ---------------------------------------------------------------------------

def create_app() -> Flask:
    base = _get_base_path()
    results = _get_results_dir()

    app = Flask(
        __name__,
        template_folder=str(base / "templates" / "web"),
        static_folder=str(base / "static") if (base / "static").exists() else None,
    )
    app.secret_key = os.urandom(24)
    app.config["RESULTS_DIR"] = results
    app.config["MAX_CONTENT_LENGTH"] = 50 * 1024 * 1024  # 50 MB

    # ------------------------------------------------------------------
    # Routes
    # ------------------------------------------------------------------

    @app.route("/")
    def index():
        reports = _list_reports(results)
        return render_template("index.html", reports=reports, version=APP_VERSION)

    @app.route("/upload", methods=["POST"])
    def upload():
        csv_files = request.files.getlist("csv_files")
        layout_file = request.files.get("layout_file")

        # Validate
        csv_files = [f for f in csv_files if f and f.filename]
        if not csv_files:
            flash("Please select at least one CSV file.", "error")
            return redirect(url_for("index"))

        # Save optional layout file
        layout_path = None
        if layout_file and layout_file.filename:
            layout_name = secure_filename(layout_file.filename)
            layout_path = results / "uploads" / layout_name
            layout_file.save(layout_path)

        last_report = None
        for csv_file in csv_files:
            csv_name = secure_filename(csv_file.filename)
            csv_path = results / "uploads" / csv_name
            csv_file.save(csv_path)

            try:
                config = load_config()
                report_path = run_pipeline(
                    csv_path=csv_path,
                    output_dir=results / "reports",
                    layout_path=layout_path,
                    history_dir=results / "history",
                    config=config,
                )
                last_report = report_path

                # Move specimen CSV from reports/ to specimens/
                plate_id = report_path.stem.replace("QC_", "")
                spec_csv = results / "reports" / f"specimens_{plate_id}.csv"
                if spec_csv.exists():
                    spec_csv.rename(results / "specimens" / spec_csv.name)

                flash(f"Report generated: {report_path.name}", "success")

            except Exception as exc:
                traceback.print_exc()
                flash(f"Error processing {csv_name}: {exc}", "error")

            finally:
                # Clean up uploaded CSV
                if csv_path.exists():
                    csv_path.unlink()

        # Clean up layout file
        if layout_path and layout_path.exists():
            layout_path.unlink()

        # Redirect to the last generated report, or back to index
        if last_report and last_report.exists():
            return redirect(url_for("view_report", filename=last_report.name))
        return redirect(url_for("index"))

    @app.route("/report/<filename>")
    def view_report(filename):
        report_file = results / "reports" / secure_filename(filename)
        if not report_file.exists():
            flash("Report not found.", "error")
            return redirect(url_for("index"))
        return send_file(report_file)

    @app.route("/download/report/<filename>")
    def download_report(filename):
        report_file = results / "reports" / secure_filename(filename)
        if not report_file.exists():
            flash("Report not found.", "error")
            return redirect(url_for("index"))
        return send_file(report_file, as_attachment=True)

    @app.route("/download/specimens/<filename>")
    def download_specimens(filename):
        spec_file = results / "specimens" / secure_filename(filename)
        if not spec_file.exists():
            flash("Specimen file not found.", "error")
            return redirect(url_for("index"))
        return send_file(spec_file, as_attachment=True)

    @app.route("/export/all")
    def export_all():
        """Export all data to date as an Excel workbook.

        Sheets:
        - specimens: all specimen results combined across plates
        - standard_curve_params: 4PL fit parameters (a, b, c, d) per plate/analyte
        - standard_curve_data: raw standard curve MFI data points
        - nc_levels: negative control MFI per plate/analyte
        """
        history_dir = results / "history"
        specimens_dir = results / "specimens"

        buf = io.BytesIO()
        with pd.ExcelWriter(buf, engine="openpyxl") as writer:
            # All specimens
            spec_frames = []
            for csv_file in sorted(specimens_dir.glob("specimens_*.csv")):
                df = pd.read_csv(csv_file, encoding="utf-8")
                plate_id = csv_file.stem.replace("specimens_", "")
                df.insert(0, "plate_id", plate_id)
                spec_frames.append(df)
            if spec_frames:
                pd.concat(spec_frames, ignore_index=True).to_excel(
                    writer, sheet_name="specimens", index=False
                )

            # Fit history (standard curve parameters)
            fit_path = history_dir / "fit_history.json"
            if fit_path.exists():
                fit_data = pd.DataFrame(json.loads(fit_path.read_text(encoding="utf-8")))
                if not fit_data.empty:
                    fit_data.to_excel(
                        writer, sheet_name="standard_curve_params", index=False
                    )

            # Standard curve raw data
            std_path = history_dir / "std_curve_history.json"
            if std_path.exists():
                std_data = pd.DataFrame(json.loads(std_path.read_text(encoding="utf-8")))
                if not std_data.empty:
                    std_data.to_excel(
                        writer, sheet_name="standard_curve_data", index=False
                    )

            # NC levels
            nc_path = history_dir / "nc_history.json"
            if nc_path.exists():
                nc_data = pd.DataFrame(json.loads(nc_path.read_text(encoding="utf-8")))
                if not nc_data.empty:
                    nc_data.to_excel(
                        writer, sheet_name="nc_levels", index=False
                    )

        buf.seek(0)
        return send_file(
            buf,
            mimetype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            as_attachment=True,
            download_name="mpxv_luminex_all_data.xlsx",
        )

    @app.route("/delete/<plate_id>", methods=["POST"])
    def delete_plate(plate_id):
        """Delete a plate's report, specimen CSV, and history entries."""
        plate_id = secure_filename(plate_id)

        # Delete report HTML
        report_file = results / "reports" / f"QC_{plate_id}.html"
        if report_file.exists():
            report_file.unlink()

        # Delete specimen CSV
        spec_file = results / "specimens" / f"specimens_{plate_id}.csv"
        if spec_file.exists():
            spec_file.unlink()

        # Remove plate from history JSON files
        history_dir = results / "history"
        for hist_file in history_dir.glob("*.json"):
            try:
                data = json.loads(hist_file.read_text(encoding="utf-8"))
                filtered = [r for r in data if r.get("plate_id") != plate_id]
                if len(filtered) < len(data):
                    hist_file.write_text(json.dumps(filtered, indent=2), encoding="utf-8")
            except Exception:
                pass

        flash(f"Deleted plate {plate_id}.", "success")
        return redirect(url_for("index"))

    @app.route("/specification")
    def specification():
        """Serve the SPECIFICATION.md as a simple HTML page."""
        spec_path = base / "SPECIFICATION.md"
        if not spec_path.exists():
            flash("Specification file not found.", "error")
            return redirect(url_for("index"))
        content = spec_path.read_text(encoding="utf-8")
        # Simple markdown-to-HTML: render as preformatted with basic styling
        html = (
            '<!DOCTYPE html><html><head><meta charset="UTF-8">'
            '<title>MPXV Luminex QC — Specification</title>'
            '<style>'
            'body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;'
            ' max-width: 900px; margin: 0 auto; padding: 20px; color: #333; }'
            'pre { white-space: pre-wrap; word-wrap: break-word; font-family: inherit;'
            ' line-height: 1.7; font-size: 14px; }'
            'a.back { display: inline-block; margin-bottom: 16px; padding: 8px 16px;'
            ' background: #3498db; color: #fff; border-radius: 6px; text-decoration: none;'
            ' font-size: 13px; font-weight: 600; }'
            '</style></head><body>'
            '<a class="back" href="/">&larr; Back to Menu</a>'
            f'<pre>{content}</pre>'
            '</body></html>'
        )
        return html

    # ------------------------------------------------------------------
    # Settings routes
    # ------------------------------------------------------------------

    @app.route("/settings")
    def settings():
        config = load_config()
        return render_template(
            "settings.html",
            config=config,
            version=APP_VERSION,
            config_path=str(get_config_path()),
        )

    @app.route("/settings", methods=["POST"])
    def save_settings():
        config = load_config()

        # Assay info
        config["assay"]["name"] = request.form.get("assay_name", "").strip()
        config["assay"]["description"] = request.form.get("assay_description", "").strip()
        config["standard"]["bead_batch"] = request.form.get("bead_batch", "").strip()

        # Panel antigens
        for i, ag in enumerate(config["panel"]["antigens"]):
            name = request.form.get(f"ag_name_{i}", "").strip()
            region = request.form.get(f"ag_region_{i}", "0")
            if name:
                ag["name"] = name
            try:
                ag["bead_region"] = int(region)
            except ValueError:
                pass

        # Panel kit controls
        for i, kc in enumerate(config["panel"]["kit_controls"]):
            name = request.form.get(f"kc_name_{i}", "").strip()
            region = request.form.get(f"kc_region_{i}", "0")
            if name:
                kc["name"] = name
            try:
                kc["bead_region"] = int(region)
            except ValueError:
                pass

        # Well classification patterns
        pc_pats = request.form.get("pc_patterns", "")
        nc_pats = request.form.get("nc_patterns", "")
        config["well_classification"]["pc_patterns"] = [p.strip() for p in pc_pats.split(",") if p.strip()]
        config["well_classification"]["nc_patterns"] = [p.strip() for p in nc_pats.split(",") if p.strip()]

        # Standard curve
        dilution_mode = request.form.get("dilution_mode", "auto")
        if dilution_mode == "auto":
            config["standard"]["dilutions"] = "auto"
        else:
            manual = request.form.get("manual_dilutions", "")
            try:
                config["standard"]["dilutions"] = [int(d.strip()) for d in manual.split(",") if d.strip()]
            except ValueError:
                flash("Invalid dilution values. Use comma-separated integers.", "error")
                return redirect(url_for("settings"))

        # Specimen dilution
        try:
            config["specimens"]["default_dilution"] = int(request.form.get("specimen_default_dilution", 100))
        except ValueError:
            pass

        # QC thresholds
        qc = config["qc_thresholds"]
        for key in ("bead_count_min", "nc_bead_mfi_max", "scg_mfi_min"):
            try:
                qc[key] = int(request.form.get(key, qc[key]))
            except (ValueError, TypeError):
                pass
        for key in ("pc_cv_threshold", "recovery_tolerance"):
            try:
                qc[key] = float(request.form.get(key, qc[key]))
            except (ValueError, TypeError):
                pass
        # Range thresholds
        try:
            qc["fc_mfi_range"] = [
                int(request.form.get("fc_mfi_min", qc["fc_mfi_range"][0])),
                int(request.form.get("fc_mfi_max", qc["fc_mfi_range"][1])),
            ]
        except (ValueError, TypeError):
            pass
        try:
            qc["ic_mfi_range"] = [
                int(request.form.get("ic_mfi_min", qc["ic_mfi_range"][0])),
                int(request.form.get("ic_mfi_max", qc["ic_mfi_range"][1])),
            ]
        except (ValueError, TypeError):
            pass

        save_config(config)
        flash("Settings saved.", "success")
        return redirect(url_for("settings"))

    @app.route("/settings/reset", methods=["POST"])
    def reset_settings():
        reset_config()
        flash("Settings reset to defaults.", "success")
        return redirect(url_for("settings"))

    @app.route("/download/plate-layout-template")
    def download_plate_layout_template():
        """Generate and serve a blank plate layout XLSX template."""
        import openpyxl
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.title = "Plate Layout"
        ws.append(["well", "sample_id", "visit_date", "dilution"])
        # Pre-fill well IDs for a 96-well plate
        for row_letter in "ABCDEFGH":
            for col_num in range(1, 13):
                ws.append([f"{row_letter}{col_num}", "", "", ""])
        buf = io.BytesIO()
        wb.save(buf)
        buf.seek(0)
        return send_file(
            buf,
            mimetype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            as_attachment=True,
            download_name="plate_layout_template.xlsx",
        )

    @app.route("/settings/export")
    def export_config():
        config = load_config()
        buf = io.BytesIO()
        buf.write(yaml.dump(config, default_flow_style=False, sort_keys=False, allow_unicode=True).encode("utf-8"))
        buf.seek(0)
        return send_file(buf, mimetype="text/yaml", as_attachment=True, download_name="mpxv_luminex_config.yaml")

    @app.route("/settings/import", methods=["POST"])
    def import_config():
        config_file = request.files.get("config_file")
        if not config_file or not config_file.filename:
            flash("No file selected.", "error")
            return redirect(url_for("settings"))
        try:
            content = config_file.read().decode("utf-8")
            imported = yaml.safe_load(content)
            if not isinstance(imported, dict):
                raise ValueError("Invalid YAML structure")
            save_config(imported)
            flash("Configuration imported.", "success")
        except Exception as exc:
            flash(f"Import failed: {exc}", "error")
        return redirect(url_for("settings"))

    @app.route("/shutdown", methods=["POST"])
    def shutdown():
        os.kill(os.getpid(), signal.SIGINT)
        return "Shutting down...", 200

    return app


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _list_reports(results_dir: Path) -> list[dict]:
    """List past reports sorted by modification time (most recent first)."""
    reports_dir = results_dir / "reports"
    specimens_dir = results_dir / "specimens"
    reports = []

    for html_file in sorted(reports_dir.glob("QC_*.html"), key=os.path.getmtime, reverse=True):
        plate_id = html_file.stem.replace("QC_", "")
        mtime = datetime.fromtimestamp(html_file.stat().st_mtime)

        # Check for matching specimen CSV
        spec_csv = specimens_dir / f"specimens_{plate_id}.csv"

        reports.append({
            "plate_id": plate_id,
            "filename": html_file.name,
            "date": mtime.strftime("%Y-%m-%d %H:%M"),
            "specimen_csv": spec_csv.name if spec_csv.exists() else None,
        })

    return reports
