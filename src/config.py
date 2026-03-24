"""Panel definition and default thresholds for 12-plex MPXV Luminex assay."""

APP_VERSION = "0.2"

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

# Kit control bead thresholds
NC_BEAD_MFI_MAX = 150       # 09 NC: flag if above 150
SCG_MFI_MIN = 10000         # 10 ScG: flag if below 10,000
FC_MFI_RANGE = (2000, 5000) # 11 FC: flag if outside 2000–5000
IC_MFI_RANGE = (2000, 3000) # 12 IC: flag if outside 2000–3000
