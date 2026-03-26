"""Kit control bead QC: NC bead, ScG, FC, IC."""

import pandas as pd

from .config import FC_MFI_RANGE, IC_MFI_RANGE, NC_BEAD_MFI_MAX, SCG_MFI_MIN


def qc_kit_controls(df: pd.DataFrame, config: dict | None = None) -> dict:
    """QC the 4 kit control beads across all wells.

    Returns dict with:
        nc_bead: DataFrame of per-well NC bead MFI with flags
        scg: DataFrame of per-well ScG MFI with flags
        fc: DataFrame of per-well FC MFI with flags
        ic: DataFrame of per-well IC MFI with flags
        flagged_wells: set of wells with any kit control flag
        n_flagged: total number of flagged wells
    """
    # Load thresholds from config or use defaults
    if config is not None:
        qc = config.get("qc_thresholds", {})
        nc_max = qc.get("nc_bead_mfi_max", NC_BEAD_MFI_MAX)
        scg_min = qc.get("scg_mfi_min", SCG_MFI_MIN)
        fc_range = tuple(qc.get("fc_mfi_range", list(FC_MFI_RANGE)))
        ic_range = tuple(qc.get("ic_mfi_range", list(IC_MFI_RANGE)))
    else:
        nc_max = NC_BEAD_MFI_MAX
        scg_min = SCG_MFI_MIN
        fc_range = FC_MFI_RANGE
        ic_range = IC_MFI_RANGE

    # Resolve kit control names from config
    kc_names = {"NC": "NC", "ScG": "ScG", "FC": "FC", "IC": "IC"}
    if config is not None:
        kc_list = config.get("panel", {}).get("kit_controls", [])
        for kc in kc_list:
            name = kc.get("name", "")
            # Map by position/convention
            if name:
                for key in kc_names:
                    if name.upper().startswith(key.upper()):
                        kc_names[key] = name
                        break

    results = {}
    flagged_wells = set()

    # NC bead
    nc_bead = df[df["analyte"] == kc_names["NC"]][["well", "sample_name", "mfi"]].copy()
    nc_bead["flag"] = nc_bead["mfi"] > nc_max
    results["nc_bead"] = nc_bead
    flagged_wells.update(nc_bead[nc_bead["flag"]]["well"].tolist())

    # ScG
    scg = df[df["analyte"] == kc_names["ScG"]][["well", "sample_name", "mfi"]].copy()
    scg["flag"] = scg["mfi"] < scg_min
    results["scg"] = scg
    flagged_wells.update(scg[scg["flag"]]["well"].tolist())

    # FC
    fc = df[df["analyte"] == kc_names["FC"]][["well", "sample_name", "mfi"]].copy()
    fc["flag"] = (fc["mfi"] < fc_range[0]) | (fc["mfi"] > fc_range[1])
    results["fc"] = fc
    flagged_wells.update(fc[fc["flag"]]["well"].tolist())

    # IC
    ic = df[df["analyte"] == kc_names["IC"]][["well", "sample_name", "mfi"]].copy()
    ic["flag"] = (ic["mfi"] < ic_range[0]) | (ic["mfi"] > ic_range[1])
    results["ic"] = ic
    flagged_wells.update(ic[ic["flag"]]["well"].tolist())

    results["flagged_wells"] = flagged_wells
    results["n_flagged"] = len(flagged_wells)
    return results
