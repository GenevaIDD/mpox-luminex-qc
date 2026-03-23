"""4-parameter logistic (4PL) curve fitting for PC standard curves."""

import warnings

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from .config import MPXV_ANTIGENS


def four_pl(x, a, b, c, d):
    """4PL model: y = d + (a - d) / (1 + (x / c)^b)

    Parameters:
        a: minimum asymptote (response at infinite concentration)
        b: Hill slope
        c: inflection point (EC50)
        d: maximum asymptote (response at zero concentration)
    """
    return d + (a - d) / (1.0 + (x / c) ** b)


def invert_4pl(y, a, b, c, d):
    """Invert the 4PL to get x (dilution) from y (MFI).

    Returns NaN if the value is outside the curve range.
    """
    y = np.asarray(y, dtype=float)
    ratio = (a - d) / (y - d) - 1.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # ratio must be positive for real-valued result
        valid = ratio > 0
        result = np.full_like(y, np.nan)
        result[valid] = c * ratio[valid] ** (1.0 / b)
    return result


def fit_standard_curves(df: pd.DataFrame) -> dict:
    """Fit 4PL curves to PC standard data for each antigen.

    Args:
        df: DataFrame with columns [well, sample_name, analyte, mfi, well_type, dilution]
            filtered to well_type == 'pc'

    Returns dict keyed by analyte name, each value is a dict:
        params: (a, b, c, d) tuple or None if fit failed
        fit_ok: bool
        std_data: DataFrame of the standard curve points used
        error: error message if fit failed
    """
    pc = df[df["well_type"] == "pc"].copy()
    results = {}

    for analyte in MPXV_ANTIGENS:
        adata = pc[pc["analyte"] == analyte].copy()
        if adata.empty:
            results[analyte] = {"params": None, "fit_ok": False, "std_data": adata, "error": "No PC data"}
            continue

        # Average replicates at each dilution
        means = adata.groupby("dilution")["mfi"].mean().reset_index()
        means = means.sort_values("dilution")

        x = means["dilution"].values
        y = means["mfi"].values

        params, fit_ok, error, qc_warnings = _fit_one(x, y)
        results[analyte] = {
            "params": params,
            "fit_ok": fit_ok,
            "std_data": adata,
            "mean_data": means,
            "error": error,
            "qc_warnings": qc_warnings,
        }

    return results


def _fit_one(x, y):
    """Fit 4PL to a single analyte's standard curve.

    Returns (params, fit_ok, error, warnings) where:
    - params: (a, b, c, d) tuple or None
    - fit_ok: True only if fit converges AND passes quality checks
    - error: error message if fit failed or quality check failed
    - warnings: list of quality warnings (may be non-empty even if fit_ok)
    """
    # Initial guesses
    d_init = np.max(y)   # upper asymptote (low dilution = high MFI)
    a_init = np.min(y)   # lower asymptote (high dilution = low MFI)
    c_init = np.median(x)  # inflection point
    b_init = 1.0          # Hill slope

    p0 = [a_init, b_init, c_init, d_init]
    bounds = (
        [0, 0.1, 1, 0],          # lower bounds
        [np.inf, 10, 1e6, np.inf]  # upper bounds
    )

    try:
        popt, _ = curve_fit(
            four_pl, x, y, p0=p0, bounds=bounds, maxfev=10000
        )
    except Exception as e:
        return None, False, str(e), []

    a, b, c, d = popt

    # --- Fit quality checks ---
    qc_warnings = []

    # 1. R² (goodness of fit)
    y_pred = four_pl(x, *popt)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    if r_squared < 0.95:
        qc_warnings.append(f"R²={r_squared:.3f} (< 0.95)")

    # 2. EC50 within plausible range
    if c < 10 or c > 10000:
        qc_warnings.append(f"EC50={c:.1f} outside range [10, 10000]")

    # 3. Hill slope in reasonable range
    if b < 0.3 or b > 5.0:
        qc_warnings.append(f"Hill slope={b:.2f} outside range [0.3, 5.0]")

    # 4. Dynamic range (ratio of upper to lower asymptote)
    upper = max(a, d)
    lower = max(min(a, d), 1.0)  # floor at 1 to avoid division by zero
    dynamic_range = upper / lower
    if dynamic_range < 3.0:
        qc_warnings.append(f"Dynamic range={dynamic_range:.1f}x (< 3x)")

    fit_ok = len(qc_warnings) == 0
    error = "; ".join(qc_warnings) if qc_warnings else None

    return tuple(popt), fit_ok, error, qc_warnings


def compute_concentrations(df: pd.DataFrame, fits: dict) -> pd.DataFrame:
    """Apply 4PL inversion to compute RAU for specimen wells.

    Adds an 'rau' column — Relative Antibody Units, defined as
    1 / (interpolated dilution factor). Higher RAU = more antibody.
    """
    specimens = df[df["well_type"] == "specimen"].copy()
    specimens["rau"] = np.nan

    for analyte, fit_result in fits.items():
        if not fit_result["fit_ok"]:
            continue
        a, b, c, d = fit_result["params"]
        mask = specimens["analyte"] == analyte
        mfi_vals = specimens.loc[mask, "mfi"].values
        dilution_equiv = invert_4pl(mfi_vals, a, b, c, d)
        specimens.loc[mask, "rau"] = 1.0 / dilution_equiv

    return specimens
