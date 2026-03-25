"""4-parameter logistic (4PL) curve fitting for PC standard curves."""

import warnings

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from .config import MPXV_ANTIGENS, RECOVERY_TOLERANCE


def four_pl(x, a, b, c, d):
    """4PL model: y = d + (a - d) / (1 + (x / c)^b)

    Parameters:
        a: minimum asymptote (response at infinite concentration)
        b: Hill slope
        c: inflection point (IC50)
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

        params, fit_ok, error, qc_warnings = _fit_one(x, y, x_min=x.min(), x_max=x.max())

        # Obs/Exp recovery and reportable range
        obs_exp = None
        reportable_range = None
        if params is not None:
            obs_exp = _compute_obs_exp(x, y, params, tolerance=RECOVERY_TOLERANCE)
            reportable_range = _compute_reportable_range(x, y, params, tolerance=RECOVERY_TOLERANCE)

        results[analyte] = {
            "params": params,
            "fit_ok": fit_ok,
            "std_data": adata,
            "mean_data": means,
            "error": error,
            "qc_warnings": qc_warnings,
            "obs_exp": obs_exp,
            "reportable_range": reportable_range,
        }

    return results


def _fit_one(x, y, x_min=None, x_max=None):
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

    # 2. IC50 within the tested dilution range (with some margin)
    ic50_lo = (x_min or x.min()) / 3.0
    ic50_hi = (x_max or x.max()) * 3.0
    if c < ic50_lo or c > ic50_hi:
        qc_warnings.append(f"IC50={c:.1f} outside range [{ic50_lo:.0f}, {ic50_hi:.0f}]")

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


def _compute_obs_exp(x_expected, y_observed, params, tolerance=0.30):
    """Backcalculate concentrations from MFI and compute Obs/Exp recovery %.

    For each standard point, invert the 4PL to get the "observed" dilution
    from the measured MFI, then compare to the expected (known) dilution.

    Returns a list of dicts with keys: dilution, mfi, obs_dilution, recovery_pct, in_range.
    """
    lo = (1.0 - tolerance) * 100.0
    hi = (1.0 + tolerance) * 100.0
    a, b, c, d = params
    obs_dilution = invert_4pl(y_observed, a, b, c, d)
    results = []
    for i in range(len(x_expected)):
        expected = x_expected[i]
        observed = obs_dilution[i]
        if np.isnan(observed) or expected == 0:
            recovery = np.nan
        else:
            recovery = (observed / expected) * 100.0
        in_range = not np.isnan(recovery) and lo <= recovery <= hi
        results.append({
            "dilution": expected,
            "mfi": y_observed[i],
            "obs_dilution": observed if not np.isnan(observed) else None,
            "recovery_pct": round(recovery, 1) if not np.isnan(recovery) else None,
            "in_range": in_range,
        })
    return results


def _compute_reportable_range(x, y, params, tolerance=0.20):
    """Determine the reportable range (LLOQ to ULOQ) based on Obs/Exp recovery.

    The reportable range is the dilution range where backcalculated recovery
    is within ±tolerance (default 20%) of the expected value.

    Returns dict with lloq, uloq (as 1/RAU values), lloq_dilution, uloq_dilution.
    """
    a, b, c, d = params
    obs_exp = _compute_obs_exp(x, y, params)

    # Find dilutions where recovery is within range
    valid_dilutions = [
        r["dilution"] for r in obs_exp
        if r["recovery_pct"] is not None and (100 - tolerance * 100) <= r["recovery_pct"] <= (100 + tolerance * 100)
    ]

    if not valid_dilutions:
        return {"lloq": None, "uloq": None, "lloq_dilution": None, "uloq_dilution": None}

    lloq_dilution = max(valid_dilutions)  # highest dilution = lowest concentration = LLOQ
    uloq_dilution = min(valid_dilutions)  # lowest dilution = highest concentration = ULOQ

    return {
        "lloq": 1.0 / lloq_dilution if lloq_dilution > 0 else None,
        "uloq": 1.0 / uloq_dilution if uloq_dilution > 0 else None,
        "lloq_dilution": lloq_dilution,
        "uloq_dilution": uloq_dilution,
    }


def compute_concentrations(df: pd.DataFrame, fits: dict) -> pd.DataFrame:
    """Apply 4PL inversion to compute 1/RAU for specimen wells.

    Adds columns:
      'rau' — 1/RAU (inverse Relative Antibody Units), defined as
              1 / (interpolated dilution factor). Higher value = more antibody.
      'extrapolated' — True if the specimen MFI falls outside the observed
                       standard curve range (i.e. above max or below min PC MFI).
    """
    specimens = df[df["well_type"] == "specimen"].copy()
    specimens["rau"] = np.nan
    specimens["extrapolated"] = False
    specimens["below_lloq"] = False
    specimens["above_uloq"] = False

    for analyte, fit_result in fits.items():
        # Compute 1/RAU whenever we have valid fit params, even if QC checks failed
        if fit_result.get("params") is None:
            continue
        a, b, c, d = fit_result["params"]
        mask = specimens["analyte"] == analyte
        mfi_vals = specimens.loc[mask, "mfi"].values
        dilution_equiv = invert_4pl(mfi_vals, a, b, c, d)
        rau_vals = 1.0 / dilution_equiv
        specimens.loc[mask, "rau"] = rau_vals

        # Flag specimens outside the observed standard curve MFI range
        std_data = fit_result.get("mean_data")
        if std_data is not None and not std_data.empty:
            mfi_lo = std_data["mfi"].min()
            mfi_hi = std_data["mfi"].max()
            specimens.loc[mask, "extrapolated"] = (mfi_vals < mfi_lo) | (mfi_vals > mfi_hi)

        # Flag specimens outside reportable range
        rr = fit_result.get("reportable_range")
        if rr and rr["lloq"] is not None and rr["uloq"] is not None:
            specimens.loc[mask, "below_lloq"] = rau_vals < rr["lloq"]
            specimens.loc[mask, "above_uloq"] = rau_vals > rr["uloq"]

    return specimens
