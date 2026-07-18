#!/usr/bin/env python3
"""
export_figure3_csv.py

Merges the three 1e-3 coupling time-series files into a single
figure3.csv that minimally reproduces panels (b) and (c) of the paper figure.

Panel (b) - potential energy deviations from pre-switch equilibrium:
  delta_V_harmonic_au       : ΔU_harmonic (Hartree, a.u.)
  delta_V_LJ_plus_Coulomb_au: Δ(U_LJ + U_Coulomb) (Hartree, a.u.)
  delta_V_total_au          : ΔU_total (Hartree, a.u.)
  delta_V_equilibrium_au    : 0.0 (reference line)

Panel (c) - fictive / kinetic temperatures:
  T_v_K  : vibrational fictive temperature from harmonic energy (K)
  T_s_K  : structural fictive temperature from LJ+Coulombic energy (K)
  T_k_K  : kinetic / bath temperature, constant 100 K
"""

import argparse
import csv
import sys

import numpy as np
from pathlib import Path

DEFAULT_BASE = Path("/media/extradrive/Trajectories/cavity_supercooled_archive"
                    "/final_production_run/data/derived/2026-01-29/time_series_output")


def load(path: Path) -> np.ndarray:
    """
    Load whitespace-separated numeric data, skipping any line that starts with
    '#' or whose first token is not a float (e.g. bare column-header rows).
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tokens = line.split()
            try:
                float(tokens[0])
            except ValueError:
                continue  # bare header row like "time_ps ..."
            rows.append([float(x) for x in tokens])
    return np.array(rows, dtype=np.float64)


def load_named_columns(path: Path) -> dict:
    """Load a whitespace table keyed by its bare (non-comment) header row.

    Robust to 5- vs 6-column potential-energy layouts so ``total_ha`` is never
    silently confused with a combined ``lj_coul_ha`` column.
    """
    header = None
    rows = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            try:
                float(tokens[0])
            except ValueError:
                header = tokens
                continue
            rows.append([float(x) for x in tokens])
    if header is None:
        raise ValueError(f"No column header row found in {path}")
    arr = np.array(rows, dtype=np.float64)
    if arr.shape[1] != len(header):
        raise ValueError(
            f"Header/data column mismatch in {path}: "
            f"{len(header)} names vs {arr.shape[1]} data columns"
        )
    return {name: arr[:, i] for i, name in enumerate(header)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Export figure3.csv from time-series files")
    parser.add_argument("--profile", default="paper", help="Dataset profile (default: paper)")
    parser.add_argument("--staged-root", type=Path, default=None, help="Override staged data root")
    parser.add_argument("--coupling", type=float, default=1e-3, help="Axis value (ε or λ)")
    parser.add_argument("--output", type=Path, default=Path.home() / "Desktop" / "figure3.csv")
    args = parser.parse_args()

    REPRO_ROOT = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(REPRO_ROOT / "shared"))
    from repro_bootstrap import setup_profile

    profile = setup_profile(args, default="paper")
    if profile.name != "paper":
        base = profile.time_series_dir
        tag = profile.entry_for_axis_value(args.coupling).epsilon_tag
    else:
        base = DEFAULT_BASE
        paper_tags = {0.0: "0epos00", 3e-4: "3eneg04", 5e-4: "5eneg04", 7e-4: "7eneg04", 1e-3: "1eneg03"}
        tag = paper_tags[args.coupling]

    pe_file = base / f"coupling_{tag}_averaged_potential_energy.txt"
    fict_file = base / f"coupling_{tag}_averaged_fictive_temperatures.txt"
    fts_file = base / f"coupling_{tag}_fictive_time_series.txt"
    output = args.output

    # ------------------------------------------------------------------ #
    # Load raw arrays. PE columns are read by header name (robust to the 5- vs
    # 6-column layouts), so ``total_ha`` is never confused with ``lj_coul_ha``.
    pe_cols = load_named_columns(pe_file)
    # FICT columns: t_ps | placeholder | T_v (harmonic)
    fict = load(fict_file)
    # FTS  columns: t_ps | T_s (LJ+Coul fictive) | relax_time | predicted_tau
    fts = load(fts_file)
    kt_file = base / f"coupling_{tag}_kinetic_temperature.txt"

    FICT_T, FICT_TV = 0, 2
    FTS_T, FTS_TS   = 0, 1

    # ------------------------------------------------------------------ #
    # Match rows by time value (all files share the same time grid;
    # use the PE file's time axis as master and nearest-match the others)
    t_master = pe_cols["time_ps"]

    def align(arr: np.ndarray, t_col: int) -> np.ndarray:
        """Return rows of arr reindexed to t_master via exact match (round to 6 dp)."""
        t_src = np.round(arr[:, t_col], 6)
        t_tgt = np.round(t_master, 6)
        idx = np.searchsorted(t_src, t_tgt)
        idx = np.clip(idx, 0, len(t_src) - 1)
        return arr[idx]

    fict_a = align(fict, FICT_T)
    fts_a  = align(fts,  FTS_T)

    # ------------------------------------------------------------------ #
    # Panel (b): subtract pre-switch baseline (mean over t < 200 ps)
    pre = t_master < 200.0
    if pre.sum() == 0:
        raise ValueError("No rows with t_ps < 200 ps – check time axis.")

    # Interpolate over the small number of NaN values (appear at the ~200 ps spike)
    def interp_nans(col: np.ndarray) -> np.ndarray:
        col = col.copy()
        nans = np.isnan(col)
        if nans.any():
            idx = np.arange(len(col))
            col[nans] = np.interp(idx[nans], idx[~nans], col[~nans])
        return col

    U_harm = interp_nans(pe_cols["harmonic_ha"])
    U_tot  = interp_nans(pe_cols["total_ha"])

    eq_harm  = np.mean(U_harm[pre])
    eq_tot   = np.mean(U_tot[pre])

    dV_harm = U_harm - eq_harm
    dV_tot  = U_tot - eq_tot

    # Combined LJ+Coulomb deviation: prefer an explicit combined column, else
    # sum the separate LJ and Coulomb columns (Coulomb is zero for these runs).
    if "lj_coul_ha" in pe_cols:
        U_ljc  = interp_nans(pe_cols["lj_coul_ha"])
        dV_ljc = U_ljc - np.mean(U_ljc[pre])
    else:
        U_lj   = interp_nans(pe_cols["lj_ha"])
        U_coul = interp_nans(pe_cols.get("coul_ha", np.zeros_like(U_lj)))
        dV_ljc = (U_lj - np.mean(U_lj[pre])) + (U_coul - np.mean(U_coul[pre]))

    # Panel (c)
    T_v = fict_a[:, FICT_TV]
    T_s = fts_a[:,  FTS_TS]
    if kt_file.is_file():
        kt_cols = load_named_columns(kt_file)
        T_k = kt_cols["kinetic_temp_K"]
        if len(T_k) != len(t_master):
            kt_stack = np.column_stack([kt_cols["time_ps"], T_k])
            T_k = align(kt_stack, 0)[:, 1]
    else:
        T_k = np.full(len(t_master), 100.0)

    # ------------------------------------------------------------------ #
    # Write CSV
    output.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "t_ps",
        "delta_V_harmonic_au",
        "delta_V_LJ_plus_Coulomb_au",
        "delta_V_total_au",
        "delta_V_equilibrium_au",
        "T_v_K",
        "T_s_K",
        "T_k_K",
    ]

    with open(output, "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(header)
        for i in range(len(t_master)):
            writer.writerow([
                f"{t_master[i]:.6f}",
                f"{dV_harm[i]:.8g}",
                f"{dV_ljc[i]:.8g}",
                f"{dV_tot[i]:.8g}",
                "0.0",
                f"{T_v[i]:.6f}",
                f"{T_s[i]:.6f}",
                "100.0",
            ])

    # ------------------------------------------------------------------ #
    # Sanity report
    n_pre  = pre.sum()
    n_post = (~pre).sum()
    print(f"Wrote {len(t_master)} rows  →  {output}")
    print(f"\n--- Pre-switch baseline check  (t < 200 ps, {n_pre} rows) ---")
    print(f"  mean dV_harmonic_au        : {dV_harm[pre].mean():.2e}  (expect ~0)")
    print(f"  mean dV_LJ+Coulomb_au      : {dV_ljc[pre].mean():.2e}  (expect ~0)")
    print(f"  mean dV_total_au           : {dV_tot[pre].mean():.2e}  (expect ~0)")
    print(f"\n--- Post-switch sign check     (t >= 200 ps, {n_post} rows) ---")
    print(f"  mean dV_harmonic_au        : {dV_harm[~pre].mean():.4f}  (expect > 0)")
    print(f"  mean dV_LJ+Coulomb_au      : {dV_ljc[~pre].mean():.4f}  (expect < 0)")
    print(f"\n--- Temperature range check ---")
    print(f"  T_v_K  : min={T_v.min():.1f}  max={T_v.max():.1f} K")
    print(f"  T_s_K  : min={T_s.min():.1f}  max={T_s.max():.1f} K")
    print(f"  T_k_K  : min={T_k.min():.1f}  max={T_k.max():.1f} K")
    print(f"\nFirst 3 data rows (columns):")
    print("  " + ",".join(header))
    for i in range(3):
        print(f"  {t_master[i]:.3f},{dV_harm[i]:.3e},{dV_ljc[i]:.3e},"
              f"{dV_tot[i]:.3e},0.0,{T_v[i]:.2f},{T_s[i]:.2f},100.0")


if __name__ == "__main__":
    main()
