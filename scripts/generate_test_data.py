"""Generate synthetic xPONENT CSV files for testing the QC pipeline.

Creates plates with:
- Plate 2: Normal plate with PC standard (duplicate), NC, and 80 specimens. Clean data.
- Plate 3: Plate with issues: some low bead counts, high replicate CV, elevated NC.
- Plate 4: Single PC replicate (no duplicate), NC dropped, 88 specimens.
"""

import csv
import io
import random
from pathlib import Path

import numpy as np

# 12-plex bead names
ANALYTES = [
    "01 MVA Ag", "02 VACV A33R", "03 MPXV A35R", "04 MPXV B6R",
    "05 MPXV A27", "06 MPXV E8L", "07 MPXV H3L", "08 MPXV M1R",
    "09 NC", "10 ScG", "11 FC", "12 IC",
]

PC_DILUTIONS = [50, 100, 200, 400, 800, 1600, 3200]

# Realistic 4PL parameters per antigen (d=max at low dilution, a=min at high dilution)
ANTIGEN_PARAMS = {
    "01 MVA Ag":    (20, 0.8, 200, 1500),
    "02 VACV A33R": (50, 0.9, 150, 30000),
    "03 MPXV A35R": (30, 0.85, 180, 28000),
    "04 MPXV B6R":  (25, 0.9, 100, 20000),
    "05 MPXV A27":  (30, 0.95, 120, 25000),
    "06 MPXV E8L":  (20, 0.85, 160, 20000),
    "07 MPXV H3L":  (25, 0.9, 140, 22000),
    "08 MPXV M1R":  (15, 1.0, 100, 8000),
}

# Kit control typical values
KIT_CONTROL_VALUES = {
    "09 NC": (10, 20),      # mean, sd — should be low
    "10 ScG": (18000, 2000), # high — sample control
    "11 FC": (3000, 300),    # fluorescent control
    "12 IC": (2300, 50),     # instrument control — very consistent
}


def four_pl(x, a, b, c, d):
    return d + (a - d) / (1 + (x / c) ** b)


def generate_pc_mfi(dilution, analyte, noise_pct=0.05):
    """Generate MFI for a PC standard well."""
    if analyte in ANTIGEN_PARAMS:
        a, b, c, d = ANTIGEN_PARAMS[analyte]
        mfi = four_pl(dilution, a, b, c, d)
        mfi *= (1 + np.random.normal(0, noise_pct))
        return max(5, round(mfi, 1))
    return generate_kit_control(analyte)


def generate_nc_mfi(analyte, elevated=False):
    """Generate MFI for a negative control well."""
    if analyte in ANTIGEN_PARAMS:
        base = np.random.uniform(8, 25)
        if elevated:
            base = np.random.uniform(150, 400)  # problematic
        return round(base, 1)
    return generate_kit_control(analyte)


def generate_specimen_mfi(analyte, serostatus="mixed"):
    """Generate specimen MFI. Serostatus: 'pos', 'neg', or 'mixed'."""
    if analyte in ANTIGEN_PARAMS:
        a, b, c, d = ANTIGEN_PARAMS[analyte]
        if serostatus == "neg":
            dilution_equiv = np.random.uniform(2000, 10000)
        elif serostatus == "pos":
            dilution_equiv = np.random.uniform(50, 500)
        else:
            dilution_equiv = np.random.lognormal(5.5, 1.5)
        mfi = four_pl(dilution_equiv, a, b, c, d)
        mfi *= (1 + np.random.normal(0, 0.08))
        return max(5, round(mfi, 1))
    return generate_kit_control(analyte)


def generate_kit_control(analyte, low_scg=False):
    """Generate kit control bead MFI."""
    mean, sd = KIT_CONTROL_VALUES[analyte]
    if low_scg and analyte == "10 ScG":
        return round(np.random.uniform(50, 200), 1)  # sample missing
    return round(max(5, np.random.normal(mean, sd)), 1)


def generate_bead_count(low=False):
    """Generate bead count for one well-analyte."""
    if low:
        return random.randint(15, 45)
    return random.randint(80, 180)


def build_plate_wells(config):
    """Build list of (well, sample_name) tuples for a plate."""
    wells = []
    row_letters = "ABCDEFGH"

    # PC standards (columns 1-2 or just 1)
    pc_cols = config.get("pc_columns", [1, 2])
    for col in pc_cols:
        for i, dil in enumerate(PC_DILUTIONS):
            well = f"{row_letters[i]}{col}"
            wells.append((well, f"PC 1:{dil}"))
        # NC in row H
        if config.get("include_nc", True):
            wells.append((f"H{col}", "NC"))

    # Specimens fill remaining wells
    specimen_col_start = max(pc_cols) + 1
    specimen_idx = 1
    for col in range(specimen_col_start, 13):
        for row in row_letters:
            well = f"{row}{col}"
            # Generate a realistic sample name
            subject = 1800 + random.randint(0, 99)
            sample_num = random.randint(1, 6)
            wells.append((well, f"{subject} S{sample_num}"))
            specimen_idx += 1

    return wells


def generate_plate_csv(plate_config):
    """Generate a full xPONENT CSV string for one plate."""
    wells = build_plate_wells(plate_config)
    batch = plate_config["batch"]
    plate_num = plate_config["plate_num"]
    date = plate_config.get("date", "03/25/2026")
    operator = plate_config.get("operator", "VENKAT")

    # Generate MFI and count data per well
    well_data = []
    for well, sample_name in wells:
        row_mfi = {}
        row_count = {}
        is_pc = sample_name.startswith("PC")
        is_nc = sample_name.startswith("NC")

        for analyte in ANALYTES:
            # Bead count
            low_beads = plate_config.get("low_bead_wells", set())
            row_count[analyte] = generate_bead_count(low=(well in low_beads))

            # MFI
            if is_pc:
                dil = int(sample_name.split(":")[1])
                noise = plate_config.get("pc_noise", 0.05)
                row_mfi[analyte] = generate_pc_mfi(dil, analyte, noise_pct=noise)
            elif is_nc:
                elevated = well in plate_config.get("elevated_nc_wells", set())
                row_mfi[analyte] = generate_nc_mfi(analyte, elevated=elevated)
            else:
                # Specimen
                low_scg = well in plate_config.get("missing_sample_wells", set())
                if analyte in KIT_CONTROL_VALUES:
                    row_mfi[analyte] = generate_kit_control(analyte, low_scg=low_scg)
                else:
                    status = random.choice(["pos", "neg", "mixed", "mixed", "mixed"])
                    row_mfi[analyte] = generate_specimen_mfi(analyte, serostatus=status)

        well_data.append((well, sample_name, row_mfi, row_count))

    # Build CSV
    lines = []

    # Header
    lines.append(f'"Program","xPONENT","","MAGPIX"')
    lines.append(f'"Build","4.3.309.1"')
    lines.append(f'"Date","{date}","3:00 PM"')
    lines.append(f'')
    lines.append(f'"SN","MAGPX18258722"')
    lines.append(f'"Batch","{batch}"')
    lines.append(f'"Version","1"')
    lines.append(f'"Operator","{operator}"')
    lines.append(f'"ComputerName","DESKTOP-3O8RKHL"')
    lines.append(f'"Country Code","7F"')
    lines.append(f'"ProtocolName","TetracoreFIA 12-Plex MPXV HIgG Test"')
    lines.append(f'"ProtocolVersion","1.0"')
    lines.append(f'"ProtocolDescription","Mpox HIgG Antibody test"')
    lines.append(f'"ProtocolDevelopingCompany","Tetracore Inc."')
    lines.append(f'"SampleWash","Off"')
    lines.append(f'"SampleVolume","50 uL"')
    lines.append(f'"BatchStartTime","{date} 3:22:19 PM"')
    lines.append(f'"BatchStopTime","{date} 3:51:22 PM"')
    lines.append(f'"BatchDescription","<None>"')
    lines.append(f'"ProtocolPlate","Name","Griener Black microclear 96 well plate","Type","96","Plates","1"')
    lines.append(f'"ProtocolMicrosphere","Map","BP 50 regions","Type","MagPlex","Count","12"')
    lines.append(f'"ProtocolAnalysis","Off"')
    lines.append(f'"NormBead","None"')
    lines.append(f'"ProtocolHeater","Off"')
    lines.append(f'')
    lines.append(f'"Most Recent Calibration and Verification Results:"')
    lines.append(f'"Last CAL Calibration","Passed 03/18/2026 13:43:33"')
    lines.append(f'"Last VER Verification","Passed 03/18/2026 13:44:35"')
    lines.append(f'"Last Fluidics Test","Passed 03/18/2026 13:45:55"')
    lines.append(f'""')
    lines.append(f'"CALInfo:"')
    lines.append(f'"Calibrator"')
    lines.append(f'"Lot","ExpirationDate","CalibrationTime","CL1Temp","CL2Temp","RP1LongTemp","RP1ShortTemp","CL1Current","CL2Current","RP1LongCurrent","RP1ShortCurrent","CL1Factor","CL2Factor","RP1LongFactor","RP1ShortFactor","Result","MachineSerialNo"')
    lines.append(f'"C09108","08/21/2026","3/18/2026 1:43:33 PM","31.7","31.7","31.7","31.7","229","232","338","338","0.007","0.00695","0.00391","0.11535","Pass","MAGPX18258722"')
    lines.append(f'')
    lines.append(f'')
    lines.append(f'"Samples","{len(well_data)}","Min Events","50","Per Bead"')
    lines.append(f'')
    lines.append(f'"Results"')
    lines.append(f'')

    # Median block
    _add_data_block(lines, "Median", well_data, lambda wd: wd[2])

    lines.append(f'')

    # Net MFI block (same as Median for this instrument)
    _add_data_block(lines, "Net MFI", well_data, lambda wd: wd[2])

    lines.append(f'')

    # Count block
    _add_data_block(lines, "Count", well_data, lambda wd: wd[3])

    # Trailing sections
    lines.append(f'')
    lines.append(f'"DataType:","Units"')
    lines.append(f'"Analyte:","01 MVA Ag","02 VACV A33R","03 MPXV A35R","04 MPXV B6R","05 MPXV A27","06 MPXV E8L","07 MPXV H3L","08 MPXV M1R","09 NC","10 ScG","11 FC","12 IC"')
    lines.append(f'"BeadID:","30","21","33","29","65","20","19","35","54","52","53","45"')
    lines.append(f'"Units:","","","","","","","","","","","",""')
    lines.append(f'')
    lines.append(f'"DataType:","Per Bead Count"')
    lines.append(f'"Analyte:","01 MVA Ag","02 VACV A33R","03 MPXV A35R","04 MPXV B6R","05 MPXV A27","06 MPXV E8L","07 MPXV H3L","08 MPXV M1R","09 NC","10 ScG","11 FC","12 IC"')
    lines.append(f'"BeadID:","30","21","33","29","65","20","19","35","54","52","53","45"')
    lines.append(f'"Per Bead:","50","50","50","50","50","50","50","50","50","50","50","50"')
    lines.append(f'')
    lines.append(f'"DataType:","Dilution Factor"')
    lines.append(f'"Location","Sample","Dilution Factor"')
    for i, (well, sample_name, _, _) in enumerate(well_data, 1):
        col = int(well[1:])
        row = well[0]
        lines.append(f'"{i}(1,{well})","{sample_name}","1"')
    lines.append(f'')
    lines.append(f'')
    lines.append(f'"DataType:","Analysis Types"')
    lines.append(f'"Analyte:","01 MVA Ag","02 VACV A33R","03 MPXV A35R","04 MPXV B6R","05 MPXV A27","06 MPXV E8L","07 MPXV H3L","08 MPXV M1R","09 NC","10 ScG","11 FC","12 IC"')
    lines.append(f'"AnalysisType","None","None","None","None","None","None","None","None","None","None","None","None"')
    lines.append(f'')
    lines.append(f'"DataType:","Audit Logs"')
    lines.append(f'"UserId","Date","Message"')
    lines.append(f'')
    lines.append(f'"DataType:","Warnings/Errors"')
    lines.append(f'"Location","Status","Message"')

    return "\n".join(lines)


def _add_data_block(lines, dtype, well_data, value_fn):
    """Add a data block (Median, Net MFI, or Count) to CSV lines."""
    lines.append(f'"DataType:","{dtype}"')
    header = '"Location","Sample",' + ",".join(f'"{a}"' for a in ANALYTES) + ',"Total Events"'
    lines.append(header)

    for i, (well, sample_name, mfi_dict, count_dict) in enumerate(well_data, 1):
        values = value_fn((well, sample_name, mfi_dict, count_dict))
        total_events = sum(count_dict.values()) if dtype != "Count" else sum(values.values())
        vals = ",".join(f'"{values[a]}"' for a in ANALYTES)
        lines.append(f'"{i}(1,{well})","{sample_name}",{vals},"{total_events}"')


def main():
    np.random.seed(42)
    random.seed(42)

    output_dir = Path(__file__).parent.parent / "data" / "raw"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Plate 2: Clean plate with specimens
    plate2 = generate_plate_csv({
        "batch": "A260325-MP1822-KV02-Plate02-12PlxMPXVHIg",
        "plate_num": 2,
        "date": "03/25/2026",
        "pc_columns": [1, 2],
        "include_nc": True,
        "pc_noise": 0.03,
    })
    (output_dir / "A260325-MP1822-KV02-Plate02-12PlxMPXVHIg_20260325_150000.csv").write_text(plate2)
    print("Generated Plate 2 (clean)")

    # Plate 3: Plate with QC issues
    plate3 = generate_plate_csv({
        "batch": "A260326-MP1822-KV02-Plate03-12PlxMPXVHIg",
        "plate_num": 3,
        "date": "03/26/2026",
        "pc_columns": [1, 2],
        "include_nc": True,
        "pc_noise": 0.15,  # high replicate variability
        "low_bead_wells": {"C5", "D5", "E7", "F7"},  # low bead counts
        "elevated_nc_wells": {"H1"},  # elevated NC
        "missing_sample_wells": {"A8", "B8"},  # low ScG — sample not added
    })
    (output_dir / "A260326-MP1822-KV02-Plate03-12PlxMPXVHIg_20260326_150000.csv").write_text(plate3)
    print("Generated Plate 3 (with QC issues)")

    # Plate 4: Single PC replicate, no NC
    plate4 = generate_plate_csv({
        "batch": "A260327-MP1822-KV02-Plate04-12PlxMPXVHIg",
        "plate_num": 4,
        "date": "03/27/2026",
        "pc_columns": [1],
        "include_nc": False,
        "pc_noise": 0.04,
    })
    (output_dir / "A260327-MP1822-KV02-Plate04-12PlxMPXVHIg_20260327_150000.csv").write_text(plate4)
    print("Generated Plate 4 (single PC replicate, no NC)")


if __name__ == "__main__":
    main()
