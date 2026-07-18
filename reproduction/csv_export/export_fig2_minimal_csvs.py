#!/usr/bin/env python3
"""
Export REAL simulation data as minimal CSV files for paper-style Fig. 2.

Sources (no synthetic values):
  IR spectra:
    second_rigorous_check/data/derived/2025-10-03/average_spectrum_dir_{0,3,5,7,10}.txt
    (dir index maps to epsilon via peak_analysis.txt; dirs {0,3,5,7,10} → lambda {0,0.042,0.070,0.098,0.141})
  F(k,t):
    final_production_run/data/derived/2026-01-29/cavity_coupling_*/master_fskt_ref*.txt
    (replica-averaged; ref1..ref14 correspond to waiting times t_w = 0..2600 ps)
  Relaxation tau:
    final_production_run/data/derived/2026-01-29/relaxation_vs_reference_time_data.txt
    (pre-computed 1/e crossing times for all 5 couplings × 15 ref_nums)

Outputs (to --out_dir, default ~/Desktop/):
  figure2a_irspectrum.csv
  figure2b_fkt_waitingtimes_lambda_{0p000,0p042,0p070,0p098,0p141}.csv  (5 files)
  figure3b_relaxation_vs_lambda.csv
  figure3c_relaxation_vs_tw.csv
"""

import argparse
import csv
import os

ARCHIVE     = "/media/extradrive/Trajectories/cavity_supercooled_archive"
FPR_DERIVED = os.path.join(ARCHIVE, "final_production_run/data/derived/2026-01-29")
SRC_DERIVED = os.path.join(ARCHIVE, "second_rigorous_check/data/derived/2025-10-03")

# Five coupling strengths matching Fig. 2
# (archive dir name, lambda_au label used in figure, CSV filename tag)
COUPLINGS = [
    ("cavity_coupling_0epos00_switch_200.0ps", 0.0,   "0p000"),
    ("cavity_coupling_3eneg04_switch_200.0ps", 0.042, "0p042"),
    ("cavity_coupling_5eneg04_switch_200.0ps", 0.070, "0p070"),
    ("cavity_coupling_7eneg04_switch_200.0ps", 0.098, "0p098"),
    ("cavity_coupling_1eneg03_switch_200.0ps", 0.141, "0p141"),
]

# archive epsilon → figure lambda (for relaxation lookup)
EPS_TO_LAM = {0.0: 0.0, 3e-4: 0.042, 5e-4: 0.070, 7e-4: 0.098, 1e-3: 0.141}
LAM_TO_EPS = {v: k for k, v in EPS_TO_LAM.items()}

# second_rigorous_check directory indices for the five lambdas
IR_DIRS = {0.0: 0, 0.042: 3, 0.070: 5, 0.098: 7, 0.141: 10}

# Frequency window for IR plot (cm⁻¹)
OMEGA_MIN = 1350.0
OMEGA_MAX = 2100.0

# Max lag time to write in F(k,t) files (ps)
LAG_MAX_PS = 1600.0

# ref indices to use for t_w (ref1..ref14 → t_w = 0..2600 ps)
# t_w = ref_time_ps - 200  where ref_time_ps = ref_num * 200
REF_NUMS = list(range(1, 15))  # 1..14 inclusive


# ── helpers ───────────────────────────────────────────────────────────────────

def read_data_lines(path):
    """Read a whitespace/tab-separated data file; skip # comments and text headers."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # skip column-name header rows (start with a letter)
            if line[0].isalpha() or line[0] == '_':
                continue
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                continue
    return rows


def write_csv(path, header_lines, col_names, rows):
    with open(path, 'w', newline='') as fh:
        for line in header_lines:
            fh.write('# ' + line + '\n')
        w = csv.writer(fh)
        w.writerow(col_names)
        w.writerows(rows)
    print(f"  Wrote {path}  ({len(rows)} data rows)")


# ── 1. IR spectra ─────────────────────────────────────────────────────────────

def build_ir_csv(out_dir):
    rows = []
    for _coupling_dir, lam, _tag in COUPLINGS:
        dir_idx = IR_DIRS[lam]
        spec_path = os.path.join(SRC_DERIVED, f"average_spectrum_dir_{dir_idx}.txt")
        if not os.path.exists(spec_path):
            raise FileNotFoundError(f"IR spectrum file missing: {spec_path}")
        data = read_data_lines(spec_path)
        count = 0
        for row in data:
            omega = row[0]
            intensity = row[1]
            if OMEGA_MIN <= omega <= OMEGA_MAX:
                rows.append([round(omega, 6), lam, round(intensity, 6)])
                count += 1
        print(f"    lambda={lam}: {count} spectral points (dir {dir_idx})")
    write_csv(
        os.path.join(out_dir, "figure2a_irspectrum.csv"),
        [
            "IR absorption spectrum: Average_Intensity vs Frequency for each coupling lambda",
            "Source: second_rigorous_check/data/derived/2025-10-03/average_spectrum_dir_*.txt",
            "Epsilon mapping: dir 0=lambda 0, dir 3=lambda 0.042, dir 5=lambda 0.070,"
            " dir 7=lambda 0.098, dir 10=lambda 0.141",
            f"Frequency window: {OMEGA_MIN}--{OMEGA_MAX} cm-1",
            "Units: omega_cm1 [cm-1]; lambda_au [a.u.]; n_alpha_au [Average_Intensity, arb. units]",
        ],
        ["omega_cm1", "lambda_au", "n_alpha_au"],
        rows,
    )


# ── 2. F(k,t) per coupling ────────────────────────────────────────────────────

def build_fkt_csv(coupling_dir, lam, tag, out_dir):
    base = os.path.join(FPR_DERIVED, coupling_dir)
    rows = []
    for ref_num in REF_NUMS:
        t_w = float((ref_num - 1) * 200)  # t_w = ref_time_ps - 200
        fkt_path = os.path.join(base, f"master_fskt_ref{ref_num}.txt")
        if not os.path.exists(fkt_path):
            print(f"    WARNING: missing {fkt_path} — skipping ref{ref_num}")
            continue
        data = read_data_lines(fkt_path)
        if not data:
            print(f"    WARNING: empty {fkt_path}")
            continue
        # Normalize by F(k, lag=0) of this ref (first data row)
        norm = data[0][1]
        if norm == 0.0:
            print(f"    WARNING: zero norm at ref{ref_num} — skipping")
            continue
        for row in data:
            lag_time, fskt = row[0], row[1]
            if lag_time > LAG_MAX_PS:
                break
            rows.append([t_w, round(lag_time, 3), round(fskt / norm, 8)])
    write_csv(
        os.path.join(out_dir, f"figure2b_fkt_waitingtimes_lambda_{tag}.csv"),
        [
            f"Normalized ISF phi_k(t-t_w; t_w) = F(k,t-t_w; t_w) / F(k,0; t_w) for lambda = {lam} a.u.",
            f"Source: {coupling_dir}/master_fskt_ref*.txt",
            "Normalization: divide each ref's F(k,t) by its own F(k, lag_time=0)",
            "t_w = (ref_num - 1) * 200 ps  (coupling switched on at 200 ps lab time)",
            f"Lag time truncated at {LAG_MAX_PS} ps to match Fig. 2(a) x-axis",
            "Units: t_w_ps [ps]; t_minus_tw_ps [ps]; phi_k_normalized [dimensionless, in (0,1]]",
        ],
        ["t_w_ps", "t_minus_tw_ps", "phi_k_normalized"],
        rows,
    )


# ── 3. Relaxation times ───────────────────────────────────────────────────────

def load_relaxation_data():
    """Return dict: (rounded_eps, ref_num) -> relaxation_time_ps."""
    path = os.path.join(FPR_DERIVED, "relaxation_vs_reference_time_data.txt")
    result = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            # columns: coupling_name  coupling_value  coupling_label  ref_num  ref_time_ps  relaxation_time_ps  normalization_value
            # coupling_name starts with 'cavity_', rest are floats/ints
            if not parts[0].startswith('cavity_'):
                continue
            if len(parts) < 6:
                continue
            try:
                eps     = float(parts[1])
                ref_num = int(parts[3])
                tau     = float(parts[5])
                result[(round(eps, 9), ref_num)] = tau
            except (ValueError, IndexError):
                continue
    return result


def build_relax_csvs(out_dir):
    relax = load_relaxation_data()
    print(f"    Loaded {len(relax)} (eps, ref_num) entries from relaxation file")

    rows_b = []   # sorted by (t_w, lambda) → Fig. 2(b)
    rows_c = []   # sorted by (lambda, t_w) → Fig. 2(c)

    for ref_num in REF_NUMS:
        t_w = float((ref_num - 1) * 200)
        eps0 = round(0.0, 9)
        tau_lam0 = relax.get((eps0, ref_num))
        if tau_lam0 is None:
            print(f"    WARNING: no lambda=0 data for ref{ref_num} — skipping t_w={t_w}")
            continue

        for _coupling_dir, lam, _tag in COUPLINGS:
            eps = round(LAM_TO_EPS[lam], 9)
            tau = relax.get((eps, ref_num))
            if tau is None:
                print(f"    WARNING: no data for lambda={lam} ref{ref_num}")
                continue
            tau_tilde = round(tau / tau_lam0, 8)
            rows_b.append([lam, t_w, tau_tilde])
            rows_c.append([lam, t_w, tau_tilde])

    rows_b.sort(key=lambda r: (r[1], r[0]))   # (t_w, lambda)
    rows_c.sort(key=lambda r: (r[0], r[1]))   # (lambda, t_w)

    common_header = [
        "Normalized relaxation time tau_s_tilde = tau_s(lambda, t_w) / tau_s(lambda=0, t_w)",
        "Source: relaxation_vs_reference_time_data.txt (1/e criterion on replica-averaged F(k,t))",
        "t_w = (ref_num - 1) * 200 ps  (cavity switched on at 200 ps lab time)",
        "Units: lambda_au [a.u.]; t_w_ps [ps]; tau_s_tilde [dimensionless]",
        "lambda=0 series is exactly 1.0 for all t_w by construction",
    ]

    write_csv(
        os.path.join(out_dir, "figure3b_relaxation_vs_lambda.csv"),
        ["Paper-style Fig. 2(b): tau_s_tilde vs lambda, colored by t_w"] + common_header,
        ["lambda_au", "t_w_ps", "tau_s_tilde"],
        rows_b,
    )
    write_csv(
        os.path.join(out_dir, "figure3c_relaxation_vs_tw.csv"),
        ["Paper-style Fig. 2(c): tau_s_tilde vs t_w, colored by lambda"] + common_header,
        ["lambda_au", "t_w_ps", "tau_s_tilde"],
        rows_c,
    )


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out_dir", default=os.path.expanduser("~/Desktop/"),
                        help="Output directory (default: ~/Desktop/)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    print(f"Writing 8 CSV files to: {args.out_dir}")

    print("\n[1/3] IR spectra...")
    build_ir_csv(args.out_dir)

    print("\n[2/3] F(k,t) per coupling...")
    for coupling_dir, lam, tag in COUPLINGS:
        print(f"  lambda={lam}  ({coupling_dir})")
        build_fkt_csv(coupling_dir, lam, tag, args.out_dir)

    print("\n[3/3] Relaxation times...")
    build_relax_csvs(args.out_dir)

    print("\nDone. 8 files written.")


if __name__ == "__main__":
    main()
