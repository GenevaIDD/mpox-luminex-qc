"""Panel definition and default thresholds for 12-plex MPXV Luminex assay."""

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

PC_DILUTIONS = {
    "1:50": 50,
    "1:100": 100,
    "1:200": 200,
    "1:400": 400,
    "1:800": 800,
    "1:1600": 1600,
    "1:3200": 3200,
}

# QC thresholds
BEAD_COUNT_MIN = 30
PC_CV_THRESHOLD = 0.25
NC_BEAD_MFI_MAX = 200
