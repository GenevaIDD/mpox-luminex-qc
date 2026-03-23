"""Tkinter GUI for mpox Luminex QC tool."""

import sys
import threading
import traceback
from pathlib import Path

from .pipeline import run_pipeline


def main():
    """Launch the QC tool GUI."""
    import tkinter as tk
    from tkinter import filedialog, messagebox

    root = tk.Tk()
    root.title("MPXV Luminex QC Tool")
    root.geometry("520x380")
    root.resizable(False, False)

    # State
    csv_path = tk.StringVar(value="")
    layout_path = tk.StringVar(value="")
    output_dir = tk.StringVar(value="")
    status_text = tk.StringVar(value="Select a CSV file to begin.")

    # --- CSV file selection ---
    tk.Label(root, text="xPONENT CSV File:", anchor="w").pack(fill="x", padx=15, pady=(15, 0))
    csv_frame = tk.Frame(root)
    csv_frame.pack(fill="x", padx=15)
    csv_entry = tk.Entry(csv_frame, textvariable=csv_path, state="readonly", width=50)
    csv_entry.pack(side="left", fill="x", expand=True)

    def select_csv():
        path = filedialog.askopenfilename(
            title="Select xPONENT CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if path:
            csv_path.set(path)
            if not output_dir.get():
                output_dir.set(str(Path(path).parent))
            status_text.set("Ready to generate report.")

    tk.Button(csv_frame, text="Browse...", command=select_csv).pack(side="right", padx=(5, 0))

    # --- Layout file selection (optional) ---
    tk.Label(root, text="Plate Layout XLSX (optional):", anchor="w").pack(fill="x", padx=15, pady=(10, 0))
    layout_frame = tk.Frame(root)
    layout_frame.pack(fill="x", padx=15)
    tk.Entry(layout_frame, textvariable=layout_path, state="readonly", width=50).pack(side="left", fill="x", expand=True)

    def select_layout():
        path = filedialog.askopenfilename(
            title="Select Plate Layout XLSX",
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
        )
        if path:
            layout_path.set(path)

    tk.Button(layout_frame, text="Browse...", command=select_layout).pack(side="right", padx=(5, 0))

    # --- Output directory ---
    tk.Label(root, text="Output Folder:", anchor="w").pack(fill="x", padx=15, pady=(10, 0))
    out_frame = tk.Frame(root)
    out_frame.pack(fill="x", padx=15)
    tk.Entry(out_frame, textvariable=output_dir, state="readonly", width=50).pack(side="left", fill="x", expand=True)

    def select_output():
        path = filedialog.askdirectory(title="Select Output Folder")
        if path:
            output_dir.set(path)

    tk.Button(out_frame, text="Browse...", command=select_output).pack(side="right", padx=(5, 0))

    # --- Generate button ---
    report_path_result = [None]

    def do_generate():
        if not csv_path.get():
            messagebox.showwarning("No file", "Please select a CSV file first.")
            return

        generate_btn.config(state="disabled")
        status_text.set("Generating report...")
        root.update()

        def worker():
            try:
                out = Path(output_dir.get()) if output_dir.get() else Path(csv_path.get()).parent
                layout = layout_path.get() or None
                report_file = run_pipeline(csv_path.get(), output_dir=out, layout_path=layout)
                report_path_result[0] = report_file
                status_text.set(f"Report saved: {report_file.name}")
                open_btn.config(state="normal")
            except Exception as e:
                traceback.print_exc()
                status_text.set(f"Error: {e}")
            finally:
                generate_btn.config(state="normal")

        threading.Thread(target=worker, daemon=True).start()

    generate_btn = tk.Button(root, text="Generate Report", command=do_generate,
                             bg="#3498db", fg="white", font=("Arial", 12, "bold"),
                             padx=20, pady=8)
    generate_btn.pack(pady=20)

    # --- Open report button ---
    def open_report():
        if report_path_result[0]:
            import webbrowser
            webbrowser.open(report_path_result[0].as_uri())

    open_btn = tk.Button(root, text="Open Report in Browser", command=open_report, state="disabled")
    open_btn.pack()

    # --- Status ---
    tk.Label(root, textvariable=status_text, fg="#666", wraplength=480).pack(pady=(10, 15))

    root.mainloop()


if __name__ == "__main__":
    # Allow running as python -m src.main
    main()
