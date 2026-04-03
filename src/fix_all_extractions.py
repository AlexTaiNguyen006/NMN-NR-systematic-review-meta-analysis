#!/usr/bin/env python3
"""
Fix extraction errors for Dollerup_2018, Yoshino_2021, Katayoshi_2023, Morifuji_2024.

All corrections derived from published paper tables:
- Dollerup_2018: Table 3 (post-treatment blood biochemistry), Table 2 (insulin)
- Yoshino_2021: Table 1 (After)
- Katayoshi_2023: Table 2 (12 weeks)
- Morifuji_2024: Table 3 (body comp 12W), Table 6 (blood chemistry 12W)

Studies verified as clean (NOT changed):
- Conze_2019: all 6 outcomes match
- Remie_2020: all 11 outcomes match
- Bandi_2025: all 2 outcomes match
- Huang_2022: already corrected in prior pass
"""

import csv
import math
import os
import shutil
from datetime import datetime

DATA_FILE = os.path.join(os.path.dirname(__file__), '..', 'data', 'extraction', 'nma_input_long.csv')

def compute_se(sd1, n1, sd2, n2):
    return math.sqrt(sd1**2/n1 + sd2**2/n2)

def make_row(study, precursor, outcome, n_t, n_c, t_mean, t_sd, c_mean, c_sd,
             unit, design, orig_disp, orig_unit, notes):
    md = t_mean - c_mean
    se = compute_se(t_sd, n_t, c_sd, n_c)
    return {
        'study_id': study,
        'precursor': precursor,
        'outcome': outcome,
        'n_treatment': str(n_t),
        'n_control': str(n_c),
        'treat_mean': f"{t_mean:.4f}",
        'treat_sd': f"{t_sd:.4f}",
        'ctrl_mean': f"{c_mean:.4f}",
        'ctrl_sd': f"{c_sd:.4f}",
        'md': f"{md:.4f}",
        'se_md': f"{se:.4f}",
        'unit': unit,
        'design': design,
        'original_dispersion': orig_disp,
        'original_unit': orig_unit,
        'notes': notes
    }

def get_dollerup_corrections():
    """
    Dollerup 2018, Table 3 post-treatment (NR=treatment, Placebo=control)
    All SEM, n=20 per group. SD = SEM * sqrt(20) = SEM * 4.4721
    """
    n = 20
    s = math.sqrt(n)
    rows = []

    # FBG: NR Post=5.6±0.1, Pbo Post=5.8±0.1 mmol/L → mg/dL (×18.018)
    c = 18.018
    rows.append(make_row("Dollerup_2018", "NR", "FBG", n, n,
        5.6*c, 0.1*s*c, 5.8*c, 0.1*s*c,
        "mg/dL", "Parallel RCT", "SEM", "mmol/L", "Need ×18.018 for mg/dL"))

    # HbA1c: NR=5.6±0.1%, Pbo=5.8±0.1%
    rows.append(make_row("Dollerup_2018", "NR", "HbA1c", n, n,
        5.6, 0.1*s, 5.8, 0.1*s,
        "%", "Parallel RCT", "SEM", "%", "From Table 3"))

    # Fasting insulin (Table 2 Basal Post): NR=65.9±5.6, Pbo=68.8±5.6 pmol/L → uU/mL (/6.945)
    d = 6.945
    rows.append(make_row("Dollerup_2018", "NR", "fasting_insulin", n, n,
        65.9/d, 5.6*s/d, 68.8/d, 5.6*s/d,
        "uU/mL", "Parallel RCT", "SEM", "pmol/L", "Need /6.945 for uU/mL"))

    # HOMA-IR: derived from corrected glucose×insulin/22.5
    # NR: 5.6*(65.9/6.945)/22.5 = 2.3604
    # Pbo: 5.8*(68.8/6.945)/22.5 = 2.5538
    # Use baseline SEM for SD estimates since post-treatment HOMA-IR not directly tabled
    # Table 1: NR 2.5±0.2(SEM), Pbo 2.8±0.3(SEM)
    rows.append(make_row("Dollerup_2018", "NR", "HOMA-IR", n, n,
        5.6*(65.9/d)/22.5, 0.2*s, 5.8*(68.8/d)/22.5, 0.3*s,
        "", "Parallel RCT", "SEM", "", "Derived from glucose×insulin/22.5"))

    # TC: NR=5.3±0.2, Pbo=5.3±0.2 mmol/L → mg/dL (×38.67)
    c2 = 38.67
    rows.append(make_row("Dollerup_2018", "NR", "TC", n, n,
        5.3*c2, 0.2*s*c2, 5.3*c2, 0.2*s*c2,
        "mg/dL", "Parallel RCT", "SEM", "mmol/L", ""))

    # LDL: NR=3.3±0.2, Pbo=3.3±0.2
    rows.append(make_row("Dollerup_2018", "NR", "LDL", n, n,
        3.3*c2, 0.2*s*c2, 3.3*c2, 0.2*s*c2,
        "mg/dL", "Parallel RCT", "SEM", "mmol/L", ""))

    # HDL: NR=1.2±0.1, Pbo=1.3±0.1
    rows.append(make_row("Dollerup_2018", "NR", "HDL", n, n,
        1.2*c2, 0.1*s*c2, 1.3*c2, 0.1*s*c2,
        "mg/dL", "Parallel RCT", "SEM", "mmol/L", ""))

    # TG: NR=1.8±0.2, Pbo=1.5±0.1 mmol/L → mg/dL (×88.57)
    c3 = 88.57
    rows.append(make_row("Dollerup_2018", "NR", "TG", n, n,
        1.8*c3, 0.2*s*c3, 1.5*c3, 0.1*s*c3,
        "mg/dL", "Parallel RCT", "SEM", "mmol/L", "TG sig increased in NR group"))

    # ALT: NR=29.6±3.6, Pbo=32.0±3.4 U/L
    rows.append(make_row("Dollerup_2018", "NR", "ALT", n, n,
        29.6, 3.6*s, 32.0, 3.4*s,
        "U/L", "Parallel RCT", "SEM", "U/L", ""))

    # BMI: No post-treatment BMI in paper. Use baseline from Table 1.
    # NR=32.4±0.8(SEM), Pbo=33.3±1.0(SEM)
    rows.append(make_row("Dollerup_2018", "NR", "BMI", n, n,
        32.4, 0.8*s, 33.3, 1.0*s,
        "kg/m2", "Parallel RCT", "SEM", "kg/m2", "Baseline values; no post-treatment BMI in paper"))

    return rows


def get_yoshino_corrections():
    """
    Yoshino 2021, Table 1 "After" column
    NMN n=13, Placebo n=12. All SEM. SD = SEM * sqrt(n)
    """
    n_t = 13  # NMN
    n_c = 12  # Placebo
    st = math.sqrt(n_t)
    sc = math.sqrt(n_c)
    rows = []

    # FBG: NMN After=5.6±0.2, Pbo After=5.6±0.2 mmol/L → mg/dL
    c = 18.018
    rows.append(make_row("Yoshino_2021", "NMN", "FBG", n_t, n_c,
        5.6*c, 0.2*st*c, 5.6*c, 0.2*sc*c,
        "mg/dL", "Parallel RCT", "SEM", "mmol/L", "Need ×18.018 for mg/dL"))

    # HbA1c: NMN After=5.5±0.1, Pbo After=5.5±0.1 %
    rows.append(make_row("Yoshino_2021", "NMN", "HbA1c", n_t, n_c,
        5.5, 0.1*st, 5.5, 0.1*sc,
        "%", "Parallel RCT", "SEM", "%", "From Table 1 After"))

    # Fasting insulin: NMN After=15.8±2.7, Pbo After=17.2±2.5 μU/mL
    rows.append(make_row("Yoshino_2021", "NMN", "fasting_insulin", n_t, n_c,
        15.8, 2.7*st, 17.2, 2.5*sc,
        "uU/mL", "Parallel RCT", "SEM", "uU/mL", ""))

    # HDL: NMN After=1.30±0.10, Pbo After=1.34±0.07 mmol/L → mg/dL
    c2 = 38.67
    rows.append(make_row("Yoshino_2021", "NMN", "HDL", n_t, n_c,
        1.30*c2, 0.10*st*c2, 1.34*c2, 0.07*sc*c2,
        "mg/dL", "Parallel RCT", "SEM", "mmol/L", ""))

    # TG: NMN After=1.26±0.18, Pbo After=1.48±0.23 mmol/L → mg/dL
    c3 = 88.57
    rows.append(make_row("Yoshino_2021", "NMN", "TG", n_t, n_c,
        1.26*c3, 0.18*st*c3, 1.48*c3, 0.23*sc*c3,
        "mg/dL", "Parallel RCT", "SEM", "mmol/L", ""))

    # Body weight: NMN After=89±4, Pbo After=87±3 kg (SEM)
    rows.append(make_row("Yoshino_2021", "NMN", "body_weight", n_t, n_c,
        89.0, 4.0*st, 87.0, 3.0*sc,
        "kg", "Parallel RCT", "SEM", "kg", ""))

    # BMI: NMN After=33.8±1.5, Pbo After=33.3±0.9 kg/m2 (SEM)
    rows.append(make_row("Yoshino_2021", "NMN", "BMI", n_t, n_c,
        33.8, 1.5*st, 33.3, 0.9*sc,
        "kg/m2", "Parallel RCT", "SEM", "kg/m2", ""))

    # SBP: NMN After=125±5, Pbo After=131±4 mmHg (SEM)
    rows.append(make_row("Yoshino_2021", "NMN", "SBP", n_t, n_c,
        125.0, 5.0*st, 131.0, 4.0*sc,
        "mmHg", "Parallel RCT", "SEM", "mmHg", ""))

    # DBP: NMN After=75±3, Pbo After=77±3 mmHg (SEM)
    rows.append(make_row("Yoshino_2021", "NMN", "DBP", n_t, n_c,
        75.0, 3.0*st, 77.0, 3.0*sc,
        "mmHg", "Parallel RCT", "SEM", "mmHg", ""))

    return rows


def get_katayoshi_corrections():
    """
    Katayoshi 2023, Table 2, "12 weeks" column
    NMN n=17, Placebo n=17. Reports SD directly.
    """
    n = 17
    rows = []

    # FBG: NMN 12wk=89.0±5.6, Pbo 12wk=89.9±6.5 mg/dL (direct, no conversion needed)
    rows.append(make_row("Katayoshi_2023", "NMN", "FBG", n, n,
        89.0, 5.6, 89.9, 6.5,
        "mg/dL", "Parallel RCT", "SD", "mg/dL", ""))

    # HDL: NMN 12wk=73.5±16.2, Pbo 12wk=74.2±16.9 mg/dL
    rows.append(make_row("Katayoshi_2023", "NMN", "HDL", n, n,
        73.5, 16.2, 74.2, 16.9,
        "mg/dL", "Parallel RCT", "SD", "mg/dL", ""))

    # LDL: NMN 12wk=120.3±24.9, Pbo 12wk=121.3±47.6 mg/dL
    rows.append(make_row("Katayoshi_2023", "NMN", "LDL", n, n,
        120.3, 24.9, 121.3, 47.6,
        "mg/dL", "Parallel RCT", "SD", "mg/dL", ""))

    # TG: NMN 12wk=70.6±42.1, Pbo 12wk=83.1±46.6 mg/dL
    rows.append(make_row("Katayoshi_2023", "NMN", "TG", n, n,
        70.6, 42.1, 83.1, 46.6,
        "mg/dL", "Parallel RCT", "SD", "mg/dL", ""))

    # Body weight: NMN 12wk=61.9±18.6, Pbo 12wk=57.6±7.8 kg
    rows.append(make_row("Katayoshi_2023", "NMN", "body_weight", n, n,
        61.9, 18.6, 57.6, 7.8,
        "kg", "Parallel RCT", "SD", "kg", ""))

    # BMI: NMN 12wk=21.9±4.3, Pbo 12wk=21.5±2.2 kg/m2
    rows.append(make_row("Katayoshi_2023", "NMN", "BMI", n, n,
        21.9, 4.3, 21.5, 2.2,
        "kg/m2", "Parallel RCT", "SD", "kg/m2", ""))

    # SBP: NMN 12wk=119.2±18.2, Pbo 12wk=127.1±13.0 mmHg
    rows.append(make_row("Katayoshi_2023", "NMN", "SBP", n, n,
        119.2, 18.2, 127.1, 13.0,
        "mmHg", "Parallel RCT", "SD", "mmHg", ""))

    # DBP: NMN 12wk=70.8±15.5, Pbo 12wk=79.4±10.3 mmHg
    rows.append(make_row("Katayoshi_2023", "NMN", "DBP", n, n,
        70.8, 15.5, 79.4, 10.3,
        "mmHg", "Parallel RCT", "SD", "mmHg", ""))

    # ALT: NMN 12wk=17.4±8.3, Pbo 12wk=22.4±9.8 U/L
    rows.append(make_row("Katayoshi_2023", "NMN", "ALT", n, n,
        17.4, 8.3, 22.4, 9.8,
        "U/L", "Parallel RCT", "SD", "U/L", ""))

    # AST: NMN 12wk=20.8±5.4, Pbo 12wk=23.1±6.1 U/L
    rows.append(make_row("Katayoshi_2023", "NMN", "AST", n, n,
        20.8, 5.4, 23.1, 6.1,
        "U/L", "Parallel RCT", "SD", "U/L", ""))

    return rows


def get_morifuji_corrections():
    """
    Morifuji 2024, Table 3 (body composition 12W), Table 6 (blood chemistry 12W)
    NMN n=30, Placebo n=29. Reports SD directly.
    NOTE: CSV had groups swapped (n_treat=29, n_ctrl=30). Correcting to NMN=30, Pbo=29.
    """
    n_t = 30  # NMN
    n_c = 29  # Placebo
    rows = []

    # FBG: NMN 12W=5.35±0.40, Pbo 12W=5.19±0.36 mmol/L → mg/dL
    c = 18.018
    rows.append(make_row("Morifuji_2024", "NMN", "FBG", n_t, n_c,
        5.35*c, 0.40*c, 5.19*c, 0.36*c,
        "mg/dL", "Parallel RCT", "SD", "mmol/L", "Need ×18.018 for mg/dL"))

    # HbA1c: NMN 12W=5.46±0.27, Pbo 12W=5.47±0.30 %
    rows.append(make_row("Morifuji_2024", "NMN", "HbA1c", n_t, n_c,
        5.46, 0.27, 5.47, 0.30,
        "%", "Parallel RCT", "SD", "%", "From Table 6"))

    # TC: NMN 12W=5.66±1.00, Pbo 12W=5.69±0.88 mmol/L → mg/dL
    c2 = 38.67
    rows.append(make_row("Morifuji_2024", "NMN", "TC", n_t, n_c,
        5.66*c2, 1.00*c2, 5.69*c2, 0.88*c2,
        "mg/dL", "Parallel RCT", "SD", "mmol/L", ""))

    # HDL: NMN 12W=1.72±0.38, Pbo 12W=1.68±0.52 mmol/L → mg/dL
    rows.append(make_row("Morifuji_2024", "NMN", "HDL", n_t, n_c,
        1.72*c2, 0.38*c2, 1.68*c2, 0.52*c2,
        "mg/dL", "Parallel RCT", "SD", "mmol/L", ""))

    # LDL: NMN 12W=3.21±0.76, Pbo 12W=3.19±0.71 mmol/L → mg/dL
    rows.append(make_row("Morifuji_2024", "NMN", "LDL", n_t, n_c,
        3.21*c2, 0.76*c2, 3.19*c2, 0.71*c2,
        "mg/dL", "Parallel RCT", "SD", "mmol/L", ""))

    # TG: NMN 12W=1.04±0.55, Pbo 12W=1.26±0.74 mmol/L → mg/dL
    c3 = 88.57
    rows.append(make_row("Morifuji_2024", "NMN", "TG", n_t, n_c,
        1.04*c3, 0.55*c3, 1.26*c3, 0.74*c3,
        "mg/dL", "Parallel RCT", "SD", "mmol/L", ""))

    # Body weight: NMN 12W=58.2±10.6, Pbo 12W=59.4±11.9 kg
    rows.append(make_row("Morifuji_2024", "NMN", "body_weight", n_t, n_c,
        58.2, 10.6, 59.4, 11.9,
        "kg", "Parallel RCT", "SD", "kg", ""))

    # BMI: NMN 12W=22.2±2.7, Pbo 12W=22.6±3.7 kg/m2
    rows.append(make_row("Morifuji_2024", "NMN", "BMI", n_t, n_c,
        22.2, 2.7, 22.6, 3.7,
        "kg/m2", "Parallel RCT", "SD", "kg/m2", ""))

    # Fat mass: NMN 12W=16.6±5.0, Pbo 12W=17.1±6.4 kg
    rows.append(make_row("Morifuji_2024", "NMN", "fat_mass", n_t, n_c,
        16.6, 5.0, 17.1, 6.4,
        "kg", "Parallel RCT", "SD", "kg", ""))

    # Skeletal muscle: NMN 12W=22.3±4.7, Pbo 12W=22.7±4.9 kg
    rows.append(make_row("Morifuji_2024", "NMN", "skeletal_muscle", n_t, n_c,
        22.3, 4.7, 22.7, 4.9,
        "kg", "Parallel RCT", "SD", "kg", ""))

    # SBP: NMN 12W=122±15, Pbo 12W=124±19 mmHg
    rows.append(make_row("Morifuji_2024", "NMN", "SBP", n_t, n_c,
        122.0, 15.0, 124.0, 19.0,
        "mmHg", "Parallel RCT", "SD", "mmHg", ""))

    # DBP: NMN 12W=74±11, Pbo 12W=77±13 mmHg
    rows.append(make_row("Morifuji_2024", "NMN", "DBP", n_t, n_c,
        74.0, 11.0, 77.0, 13.0,
        "mmHg", "Parallel RCT", "SD", "mmHg", ""))

    # ALT: NMN 12W=18.4±5.6, Pbo 12W=16.4±5.3 U/L
    rows.append(make_row("Morifuji_2024", "NMN", "ALT", n_t, n_c,
        18.4, 5.6, 16.4, 5.3,
        "U/L", "Parallel RCT", "SD", "U/L", ""))

    # AST: NMN 12W=23.3±5.1, Pbo 12W=22.1±5.0 U/L
    rows.append(make_row("Morifuji_2024", "NMN", "AST", n_t, n_c,
        23.3, 5.1, 22.1, 5.0,
        "U/L", "Parallel RCT", "SD", "U/L", ""))

    return rows


def main():
    # Read current data
    with open(DATA_FILE, 'r') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        all_rows = list(reader)

    print(f"Read {len(all_rows)} rows from {DATA_FILE}")

    # Studies to fix
    studies_to_fix = {'Dollerup_2018', 'Yoshino_2021', 'Katayoshi_2023', 'Morifuji_2024'}

    # Separate rows
    keep_rows = [r for r in all_rows if r['study_id'] not in studies_to_fix]
    removed_rows = [r for r in all_rows if r['study_id'] in studies_to_fix]

    print(f"Keeping {len(keep_rows)} rows (unchanged studies)")
    print(f"Replacing {len(removed_rows)} rows from: {', '.join(sorted(studies_to_fix))}")

    # Get corrected rows
    corrections = []
    corrections.extend(get_dollerup_corrections())
    corrections.extend(get_yoshino_corrections())
    corrections.extend(get_katayoshi_corrections())
    corrections.extend(get_morifuji_corrections())

    print(f"Generated {len(corrections)} corrected rows:")
    for study in sorted(studies_to_fix):
        n = sum(1 for r in corrections if r['study_id'] == study)
        print(f"  {study}: {n} outcomes")

    # Create backup
    backup_name = DATA_FILE.replace('.csv', f'_BACKUP2_{datetime.now().strftime("%Y%m%d")}.csv')
    shutil.copy2(DATA_FILE, backup_name)
    print(f"\nBackup saved to: {os.path.basename(backup_name)}")

    # Merge: keep unchanged + new corrections, sorted by study_id then outcome
    final_rows = keep_rows + corrections

    # Sort by study_id, then outcome for consistent ordering
    final_rows.sort(key=lambda r: (r['study_id'], r['outcome']))

    # Write
    with open(DATA_FILE, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(final_rows)

    print(f"\nWrote {len(final_rows)} rows to {os.path.basename(DATA_FILE)}")
    print(f"  Net change: {len(final_rows) - len(all_rows)} rows")

    # Summary of changes
    print("\n=== CHANGE SUMMARY ===")
    for study in sorted(studies_to_fix):
        old = [r for r in removed_rows if r['study_id'] == study]
        new = [r for r in corrections if r['study_id'] == study]
        old_outcomes = set(r['outcome'] for r in old)
        new_outcomes = set(r['outcome'] for r in new)
        print(f"\n{study}:")
        print(f"  Old outcomes ({len(old)}): {', '.join(sorted(old_outcomes))}")
        print(f"  New outcomes ({len(new)}): {', '.join(sorted(new_outcomes))}")
        added = new_outcomes - old_outcomes
        removed = old_outcomes - new_outcomes
        if added:
            print(f"  ADDED: {', '.join(sorted(added))}")
        if removed:
            print(f"  REMOVED: {', '.join(sorted(removed))}")


if __name__ == '__main__':
    main()
