"""Panel definition and default thresholds for 12-plex MPXV Luminex assay."""

MPXV_ANTIGENS = [
    "01 MVA Ag",
    "02 VACV A33R",
    "03 MPXV A35R",
    "04 MPXV B6R",
    "05 MPXV A27",
    "06 MPXV E8L",
    "07 MPXV H3L",
    "08 MPXV M1R",
]

MPXV_KIT_CONTROLS = [
    "09 NC",
    "10 ScG",
    "11 FC",
    "12 IC",
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
BEAD_COUNT_MIN = 50
PC_CV_THRESHOLD = 0.25
NC_BEAD_MFI_MAX = 200
