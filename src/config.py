"""Panel definition and default thresholds for 12-plex MPXV Luminex assay.

All configurable values have defaults here. User overrides are stored in
~/mpox-luminex-qc-results/config.yaml and loaded at runtime via settings.py.
"""

APP_VERSION = "0.3"

# --- Default configuration (also exported as DEFAULTS dict below) ---

MPXV_ANTIGENS = [
    "MVA Ag",
    "VACV A33R",
    "MPXV A35R",
    "MPXV B6R",
    "MPXV A27",
    "MPXV E8L",
    "MPXV H3L",
    "MPXV M1R",
]

MPXV_KIT_CONTROLS = [
    "NC",
    "ScG",
    "FC",
    "IC",
]

ALL_BEADS = MPXV_ANTIGENS + MPXV_KIT_CONTROLS

# QC thresholds
BEAD_COUNT_MIN = 30
PC_CV_THRESHOLD = 0.25
RECOVERY_TOLERANCE = 0.30  # ±30% Obs/Exp recovery for standard curve points

# Kit control bead thresholds
NC_BEAD_MFI_MAX = 150       # NC: flag if above 150
SCG_MFI_MIN = 5000          # ScG: flag if below 5,000
FC_MFI_RANGE = (2000, 5000) # FC: flag if outside 2000–5000
IC_MFI_RANGE = (1500, 3300) # IC: flag if outside 1500–3300

# Specimen default dilution
SPECIMEN_DEFAULT_DILUTION = 100

# Well classification patterns
PC_PATTERNS = [r"^PC\s", r"^ITM\s*PC\s"]
NC_PATTERNS = [r"^NC", r"^ITM\s*NC"]

# --- Structured defaults dict for settings.py ---

DEFAULTS = {
    "assay": {
        "name": "12-Plex MPXV HIgG",
        "description": "Tetracore FIA 12-Plex MPXV Human IgG Test",
    },
    "panel": {
        "antigens": [
            {"name": "MVA Ag", "bead_region": 15},
            {"name": "VACV A33R", "bead_region": 26},
            {"name": "MPXV A35R", "bead_region": 37},
            {"name": "MPXV B6R", "bead_region": 44},
            {"name": "MPXV A27", "bead_region": 55},
            {"name": "MPXV E8L", "bead_region": 62},
            {"name": "MPXV H3L", "bead_region": 73},
            {"name": "MPXV M1R", "bead_region": 78},
        ],
        "kit_controls": [
            {"name": "NC", "bead_region": 19},
            {"name": "ScG", "bead_region": 33},
            {"name": "FC", "bead_region": 48},
            {"name": "IC", "bead_region": 57},
        ],
    },
    "well_classification": {
        "pc_patterns": PC_PATTERNS,
        "nc_patterns": NC_PATTERNS,
    },
    "standard": {
        "dilutions": "auto",
        "bead_batch": "",
    },
    "specimens": {
        "default_dilution": SPECIMEN_DEFAULT_DILUTION,
    },
    "qc_thresholds": {
        "bead_count_min": BEAD_COUNT_MIN,
        "pc_cv_threshold": PC_CV_THRESHOLD,
        "recovery_tolerance": RECOVERY_TOLERANCE,
        "nc_bead_mfi_max": NC_BEAD_MFI_MAX,  # 150
        "scg_mfi_min": SCG_MFI_MIN,        # 5000
        "fc_mfi_range": list(FC_MFI_RANGE), # [2000, 5000]
        "ic_mfi_range": list(IC_MFI_RANGE), # [1500, 3300]
    },
}
