#!/usr/bin/env python3
"""
Export figure4.csv: minimal tidy CSV to reproduce FIG. 4 of the cavity-supercooled paper.

Four panels in a single long-format CSV:
  panel=a  |  h_lambda(t) material time — both F(k,t)-measured and TN model
  panel=b  |  Phi_k(h) data collapse + stretched-exponential fit
  panel=c  |  tau_tilde_s vs lambda (a.u.) at different t_w
  panel=d  |  tau_tilde_s vs t_w (ps) for different coupling strengths

CSV columns
-----------
  panel         : a | b | c | d
  series        : measured_xi | tn_xi | collapse_data | collapse_fit |
                  tau_vs_lambda | tau_vs_tw
  epsilon_au    : coupling constant epsilon (a.u.)
  lambda_au     : lambda = epsilon / omega_c (dimensionless, a.u.)
  t_w_ps        : waiting time (ps), convention: lab_time - 200
                  (t_w=0 when coupling turns on at lab t=200 ps)
  ref_num       : master-F(k,t) reference index (-1 for fit rows)
  x             : panel-specific x-axis coordinate
  y             : panel-specific y-axis coordinate

Panel x/y semantics
-------------------
  a : x = time (ps after coupling turn-on),  y = xi(t) [material time]
  b : x = h (material-time increment),        y = Phi_k(h) = F(k, lag)/F(k,0)
  c : x = lambda_au,                           y = tau_tilde_s
  d : x = t_w_ps,                              y = tau_tilde_s

Physical constants
------------------
  omega_c  = 1560 cm^-1  =  7.1079e-3 a.u.
  lambda   = epsilon / omega_c
  Measured xi(t) = integral dt/tau_0.1 from master F(k,t) with ref_000 normalization
  and t_w = ref_num * switch_time_ps; TN xi uses integral dt/tau_fictive (F/F0=0.1, no ln10).
"""

import argparse
import sys
import os
from pathlib import Path
from glob import glob
import numpy as np
import csv
from datetime import datetime


def _make_linear_interp(x, y):
    """
    Return a callable f(xq) that linearly interpolates (x, y),
    clamping out-of-bounds queries to the nearest endpoint.
    Replaces scipy.interpolate.interp1d to avoid version issues.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    def f(xq):
        xq = np.asarray(xq, dtype=float)
        # np.interp does linear interpolation with clamping
        return np.interp(xq, x, y)

    return f

# ── physical constants ───────────────────────────────────────────────────────
OMEGA_C_CM_INV = 1560.0
CM_INV_TO_AU   = 4.556335e-6
OMEGA_C_AU     = OMEGA_C_CM_INV * CM_INV_TO_AU   # 7.10786e-3 a.u.

# ── coupling catalogue ───────────────────────────────────────────────────────
COUPLING_VALUES = [0.0, 3e-4, 5e-4, 7e-4, 1e-3]
COUPLING_NAMES  = {
    0.0:  '0epos00',
    3e-4: '3eneg04',
    5e-4: '5eneg04',
    7e-4: '7eneg04',
    1e-3: '1eneg03',
}
# Column pairs (t_col, xi_col) in material_time_all_couplings.txt (legacy reconstruction)
MT_COLS = {
    0.0:  (0, 1),
    3e-4: (2, 3),
    5e-4: (4, 5),
    7e-4: (6, 7),
    1e-3: (8, 9),
}

def lam(eps):
    """Convert epsilon (a.u.) to lambda = epsilon / omega_c."""
    return eps / OMEGA_C_AU


# ── I/O helpers ──────────────────────────────────────────────────────────────

def load_measured_mt(data_root, switch_time_ps: float = 200.0):
    """
    Load direct measured material times: integral dt/tau_0.1 from master F(k,t).

    Returns dict: axis value -> (time_ps array, xi array).
    Time is in ps after coupling turn-on (table origin = 0).
    """
    from direct_material_time import load_direct_measured_mt_all

    result = load_direct_measured_mt_all(
        Path(data_root),
        COUPLING_VALUES,
        COUPLING_NAMES,
        switch_time_ps=switch_time_ps,
    )
    if result:
        sample_len = next(iter(result.values()))[0].size
        print(f"  Loaded direct measured material time: {sample_len} grid points per coupling")
    else:
        print("  WARNING: no direct measured material time loaded")
    return result


def load_tn_mt(data_root, switch_time_ps: float = 200.0):
    """
    Load fictive/TN material times from coupling_*_fictive_time_series.txt.

    Integrates 1/tau_fictive from cavity turn-on with no 1/ln(10) scaling.
    Returns dict: axis value -> (time_ps after coupling, xi_tn arrays).
    """
    from direct_material_time import tn_material_time_from_series

    ts_dir = Path(data_root) / 'time_series_output'
    result: dict[float, tuple[np.ndarray, np.ndarray]] = {}

    for axis_value in COUPLING_VALUES:
        tag = COUPLING_NAMES[axis_value]
        fpath = ts_dir / f'coupling_{tag}_fictive_time_series.txt'
        if not fpath.exists():
            print(f"  WARNING: TN file not found: {fpath}")
            continue
        raw = np.loadtxt(fpath, skiprows=12)
        time_ps = raw[:, 0]
        tau_f = raw[:, 2]
        time_shifted, xi_shifted, _aging_rate = tn_material_time_from_series(
            time_ps,
            tau_f,
            switch_time_ps=switch_time_ps,
        )
        result[axis_value] = (time_shifted, xi_shifted)
        print(
            f"  TN mt {axis_value:g}: {len(time_shifted)} pts, max xi={xi_shifted.max():.3f}"
        )

    return result


def load_relaxation_table(data_root, alt_root=None):
    """
    Load relaxation_vs_reference_time_data.txt.
    Returns a sorted list of (coupling_value, ref_num, ref_time_ps, relaxation_time_ps).
    """
    fpath = data_root / 'relaxation_vs_reference_time_data.txt'
    if not fpath.is_file() and alt_root is not None:
        candidate = Path(alt_root) / 'relaxation_vs_reference_time_data.txt'
        if candidate.is_file():
            fpath = candidate
    if not fpath.is_file():
        raise FileNotFoundError(
            f"relaxation_vs_reference_time_data.txt not found under {data_root}"
            + (f" or {alt_root}" if alt_root else "")
        )

    rows = []
    with open(fpath, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if 'coupling_name' in line or 'ref_num' in line:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                if parts[0].startswith('cavity_coupling_'):
                    coupling_value = float(parts[1])
                    ref_num = int(parts[3])
                    ref_time_ps = float(parts[4])
                    relaxation_time = float(parts[5])
                else:
                    ref_num = int(parts[0])
                    ref_time_ps = float(parts[1])
                    coupling_value = float(parts[4])
                    relaxation_time = float(parts[6])
            except ValueError:
                continue
            rows.append((coupling_value, ref_num, ref_time_ps, relaxation_time))
    print(f"  Loaded relaxation table: {len(rows)} entries from {fpath}")
    return rows


def load_fskt(data_root, eps, ref_num, coupling_names=None):
    """
    Load a single master_fskt_ref file.
    Returns (lag_time, fskt) arrays, or (None, None) if file missing.
    """
    names = coupling_names or COUPLING_NAMES
    coupling_dir = data_root / f'cavity_coupling_{names[eps]}_switch_200.0ps'
    padded = f"{ref_num:03d}"
    candidates = sorted(
        p for p in coupling_dir.glob("master_fskt_ref*.txt")
        if "_sample_counts" not in p.name and padded in p.stem
    )
    if not candidates:
        legacy = coupling_dir / f'master_fskt_ref{ref_num}.txt'
        fpath = legacy if legacy.exists() else None
    else:
        fpath = candidates[0]
    if fpath is None or not fpath.exists():
        return None, None
    # Header comment count drifts (RF masters add a lag-grid note). Skip '#'
    # comments and an optional non-numeric column header instead of fixed skiprows.
    rows: list[list[float]] = []
    with open(fpath, encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.lower().startswith("lag_time"):
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            rows.append([float(parts[0]), float(parts[1])])
    if not rows:
        return None, None
    data = np.asarray(rows, dtype=float)
    return data[:, 0], data[:, 1]


def load_fit_params(data_root, alt_root=None):
    """
    Load stretched exponential fit results.
    Returns dict: epsilon -> beta.
    """
    fpath = None
    for root in (data_root, alt_root):
        if root is None:
            continue
        candidate = Path(root) / 'stretched_exponential_fit_results.txt'
        if candidate.is_file():
            fpath = candidate
            break
    if fpath is None:
        raise FileNotFoundError(
            f"stretched_exponential_fit_results.txt not found under {data_root}"
            + (f" or {alt_root}" if alt_root else "")
        )
    result = {}
    with open(fpath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if 'Coupling_au' in line:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            eps  = float(parts[0])
            beta = float(parts[4])
            result[eps] = beta
    return result


# ── panel builders ────────────────────────────────────────────────────────────

def build_panel_a(meas_mt, tn_mt, downsample=10):
    """
    Panel (a): material time h_lambda(t) — measured and TN.
    Yields CSV row dicts.
    """
    print("Building panel (a)...")
    for eps in COUPLING_VALUES:
        la = lam(eps)

        # Measured
        if eps in meas_mt:
            t, xi = meas_mt[eps]
            idx = slice(None, None, downsample)
            for tv, xv in zip(t[idx], xi[idx]):
                yield dict(panel='a', series='measured_xi',
                           epsilon_au=eps, lambda_au=la,
                           t_w_ps='', ref_num='',
                           x=round(float(tv), 4), y=round(float(xv), 6))

        # TN model
        if eps in tn_mt:
            t, xi = tn_mt[eps]
            idx = slice(None, None, downsample)
            for tv, xv in zip(t[idx], xi[idx]):
                yield dict(panel='a', series='tn_xi',
                           epsilon_au=eps, lambda_au=la,
                           t_w_ps='', ref_num='',
                           x=round(float(tv), 4), y=round(float(xv), 6))

    print("  Panel (a) done.")


def build_panel_b(data_root, meas_mt, relaxation_rows, fit_params,
                  downsample=5, phi_min=0.01):
    """
    Panel (b): data collapse Phi_k(h) and stretched-exponential fit.
    h = xi(t_w + lag) - xi(t_w)  using the measured ξ(t) for each coupling.
    """
    print("Building panel (b)...")

    # Build ref_num -> ref_time_ps lookup per coupling
    ref_lookup = {}   # (eps, ref_num) -> ref_time_ps
    for (eps, rn, rt, tau) in relaxation_rows:
        ref_lookup[(eps, rn)] = rt

    # Build interpolators for measured xi(t) — t in ps after coupling
    xi_interp = {}
    for eps in COUPLING_VALUES:
        if eps not in meas_mt:
            continue
        t, xi = meas_mt[eps]
        xi_interp[eps] = _make_linear_interp(t, xi)

    # Emit collapse data for selected refs (every 3rd ref from 1..13)
    sel_refs = list(range(1, 14, 3))   # [1, 4, 7, 10, 13]

    for eps in COUPLING_VALUES:
        if eps not in xi_interp:
            continue
        la      = lam(eps)
        xi_fn   = xi_interp[eps]

        for rn in sel_refs:
            lag, fskt = load_fskt(data_root, eps, rn)
            if lag is None:
                continue

            ref_t_lab = ref_lookup.get((eps, rn))
            if ref_t_lab is None:
                continue

            # t_w in material-time-table coordinates (0 = coupling start)
            t_w_mt  = max(0.0, ref_t_lab - 200.0)
            xi_at_tw = float(xi_fn(t_w_mt))

            f0 = fskt[0]
            if f0 <= 0:
                continue

            phi = fskt / f0

            # Map lag to h = xi(t_w_mt + lag) - xi(t_w_mt)
            h   = xi_fn(t_w_mt + lag) - xi_at_tw

            # Filter: h > 0, phi > phi_min, and finite
            mask = (h > 0) & (phi > phi_min) & np.isfinite(h) & np.isfinite(phi)
            h_sel   = h[mask][::downsample]
            phi_sel = phi[mask][::downsample]

            t_w_shifted = ref_t_lab - 200.0
            for hv, pv in zip(h_sel, phi_sel):
                yield dict(panel='b', series='collapse_data',
                           epsilon_au=eps, lambda_au=la,
                           t_w_ps=round(t_w_shifted, 1), ref_num=rn,
                           x=round(float(hv), 6), y=round(float(pv), 6))

    # Stretched-exponential fit curve (universal, epsilon=0, mean beta)
    betas = list(fit_params.values())
    beta_mean = float(np.mean(betas))
    h_fit  = np.linspace(0.0, 2.0, 200)
    phi_fit = np.exp(-h_fit ** beta_mean)
    for hv, pv in zip(h_fit, phi_fit):
        yield dict(panel='b', series='collapse_fit',
                   epsilon_au='', lambda_au='',
                   t_w_ps='', ref_num=-1,
                   x=round(float(hv), 6), y=round(float(pv), 6))

    print(f"  Panel (b) done. Fit curve uses beta_mean={beta_mean:.4f}")


def build_panels_cd(relaxation_rows):
    """
    Panels (c) and (d): tau_tilde_s = tau(eps, t_w) / tau(eps=0, t_w).
    t_w_ps convention: ref_time_ps - 200 (t_w=0 when coupling turns on).
    Excludes ref_num=0 (t_w before coupling).
    """
    print("Building panels (c) and (d)...")

    # Organise into nested dict: eps -> ref_num -> tau
    tau_dict = {}
    ref_times = {}   # ref_num -> ref_time_ps (use eps=0 as ground truth)
    for (eps, rn, rt, tau) in relaxation_rows:
        tau_dict.setdefault(eps, {})[rn] = tau
        ref_times[rn] = rt

    # Build tau_tilde_s table
    tau0 = tau_dict.get(0.0, {})

    for (eps, rn, rt, tau) in relaxation_rows:
        if rn == 0:
            continue   # before coupling, skip

        tau_ref = tau0.get(rn)
        if tau_ref is None or tau_ref <= 0:
            continue

        tau_tilde = tau / tau_ref
        t_w = rt - 200.0      # shifted waiting time
        la  = lam(eps)

        # Panel (c): x = lambda, y = tau_tilde
        yield dict(panel='c', series='tau_vs_lambda',
                   epsilon_au=eps, lambda_au=la,
                   t_w_ps=round(t_w, 1), ref_num=rn,
                   x=round(la, 6), y=round(tau_tilde, 6))

        # Panel (d): x = t_w, y = tau_tilde
        yield dict(panel='d', series='tau_vs_tw',
                   epsilon_au=eps, lambda_au=la,
                   t_w_ps=round(t_w, 1), ref_num=rn,
                   x=round(t_w, 1), y=round(tau_tilde, 6))

    print("  Panels (c) and (d) done.")


# ── verification ──────────────────────────────────────────────────────────────

def verify_csv(
    output_path,
    fit_params,
    meas_mt,
    relaxation_rows,
    *,
    expected_lam_max: float | None = None,
):
    """
    Run automated sanity checks on the written CSV.
    Pure Python/NumPy implementation (no pandas) for compatibility.
    Prints a PASS/FAIL report.
    """
    # Read CSV with Python's csv module (skip # comment lines)
    rows = []
    with open(output_path, 'r', encoding='utf-8') as fh:
        reader = csv.DictReader(row for row in fh if not row.startswith('#'))
        for r in reader:
            rows.append(r)

    def to_float(v):
        try:
            return float(v)
        except (ValueError, TypeError):
            return float('nan')

    errors = []

    # 1. Lambda max
    lam_vals = [to_float(r['lambda_au']) for r in rows if r['lambda_au'] != '']
    lam_max  = max(lam_vals)
    if expected_lam_max is None:
        expected_lam_max = lam(1e-3)
    status = 'PASS' if abs(lam_max - expected_lam_max) < 0.002 else 'FAIL'
    if status == 'FAIL':
        errors.append(f"lambda_max={lam_max:.4f} expected ~{expected_lam_max:.4f}")
    print(f"  CHECK lambda_max: {lam_max:.4f}  (expected ~{expected_lam_max:.4f})  {status}")

    # 2. Panel (a): both series per coupling
    for eps in COUPLING_VALUES:
        la = round(lam(eps), 6)
        pa_rows = [r for r in rows if r['panel'] == 'a']
        has_meas = any(r['series'] == 'measured_xi' and abs(to_float(r['lambda_au']) - la) < 1e-4
                       for r in pa_rows)
        has_tn   = any(r['series'] == 'tn_xi'       and abs(to_float(r['lambda_au']) - la) < 1e-4
                       for r in pa_rows)
        status = 'PASS' if has_meas and has_tn else 'FAIL'
        if not (has_meas and has_tn):
            errors.append(f"panel a missing series for eps={eps}")
        print(f"  CHECK panel(a) eps={eps:.0e}  measured={has_meas}  tn={has_tn}  {status}")

    # 3. Measured xi at t ~ 2000 ps
    for eps in COUPLING_VALUES:
        la = round(lam(eps), 6)
        sub = [r for r in rows
               if r['panel'] == 'a' and r['series'] == 'measured_xi'
               and abs(to_float(r['lambda_au']) - la) < 1e-4
               and 1900 <= to_float(r['x']) <= 2100]
        if sub:
            xi_val = float(np.mean([to_float(r['y']) for r in sub]))
            status = 'PASS' if 5.0 < xi_val < 17.0 else 'FAIL'
            if status == 'FAIL':
                errors.append(f"xi@2000ps out of range for eps={eps}: {xi_val:.3f}")
            print(f"  CHECK xi_meas@~2000ps eps={eps:.0e}: xi={xi_val:.3f}  {status}")

    # 4. tau_tilde_s ratios
    c_rows = [r for r in rows if r['panel'] == 'c' and r['series'] == 'tau_vs_lambda']
    if c_rows:
        max_tautilde = max(to_float(r['y']) for r in c_rows)
        status = 'PASS' if 1.5 < max_tautilde < 7.0 else 'FAIL'
        if status == 'FAIL':
            errors.append(f"max tau_tilde={max_tautilde:.2f} out of [1.5,7]")
        print(f"  CHECK max tau_tilde_s: {max_tautilde:.3f}  {status}")

    # 5. Beta
    betas = list(fit_params.values())
    print(f"  CHECK beta values: {[round(b, 4) for b in betas]}  mean={np.mean(betas):.4f}")
    beta_ok = all(0.45 < b < 0.70 for b in betas)
    if not beta_ok:
        errors.append("some beta values outside [0.45, 0.70]")
    print(f"  CHECK beta range 0.45-0.70: {'PASS' if beta_ok else 'FAIL'}")

    # 6. Row counts
    n_total = len(rows)
    n_a = sum(1 for r in rows if r['panel'] == 'a')
    n_b = sum(1 for r in rows if r['panel'] == 'b')
    n_c = sum(1 for r in rows if r['panel'] == 'c')
    n_d = sum(1 for r in rows if r['panel'] == 'd')
    print(f"  CHECK row counts: total={n_total}  a={n_a}  b={n_b}  c={n_c}  d={n_d}")

    if errors:
        print(f"\n  VERIFICATION FAILED: {len(errors)} error(s):")
        for e in errors:
            print(f"    - {e}")
        return False
    else:
        print("\n  ALL CHECKS PASSED")
        return True


# ── main ──────────────────────────────────────────────────────────────────────

def configure_profile_catalogue(profile):
    """Rebind module-level coupling catalogue from a dataset profile."""
    global COUPLING_VALUES, COUPLING_NAMES, MT_COLS, lam
    COUPLING_VALUES = [entry.axis_value for entry in profile.couplings]
    COUPLING_NAMES = profile.tag_map()
    MT_COLS = {
        value: (2 * idx, 2 * idx + 1)
        for idx, value in enumerate(COUPLING_VALUES)
    }
    if profile.axis == "lambda":
        def lam(value):  # noqa: ANN001
            return float(value)
    else:
        def lam(value):  # noqa: ANN001
            return float(value) / OMEGA_C_AU
    return profile


def main():
    parser = argparse.ArgumentParser(
        description="Export figure4.csv from derived simulation data."
    )
    parser.add_argument(
        '--profile',
        default='paper',
        help='Dataset profile (default: paper)',
    )
    parser.add_argument(
        '--staged-root',
        type=Path,
        default=None,
        help='Override staged data root from profile',
    )
    parser.add_argument(
        '--data_root',
        default=None,
        help='Path to derived data directory (overrides profile staged_root)',
    )
    parser.add_argument(
        '--output',
        default=os.path.expanduser('~/Desktop/figure4.csv'),
        help='Output CSV path'
    )
    parser.add_argument(
        '--downsample_a', type=int, default=10,
        help='Downsample step for panel (a) curves (default: 10)'
    )
    parser.add_argument(
        '--downsample_b', type=int, default=5,
        help='Downsample step for panel (b) lag grids (default: 5)'
    )
    args = parser.parse_args()

    REPRO_ROOT = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(REPRO_ROOT / "shared"))
    from repro_bootstrap import setup_profile

    profile = setup_profile(args, default="paper")
    if args.profile != "paper":
        configure_profile_catalogue(profile)

    if args.data_root is not None:
        data_root = Path(args.data_root)
    elif profile.name != "paper":
        data_root = profile.staged_root
    else:
        data_root = Path(
            '/media/extradrive/Trajectories/cavity_supercooled_archive/'
            'final_production_run/data/derived/2026-01-29'
        )
    output    = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("export_figure4_csv.py")
    print(f"  data_root : {data_root}")
    print(f"  output    : {output}")
    print("=" * 70)

    # Load all data
    print("\n[1/5] Loading measured material times...")
    meas_mt = load_measured_mt(data_root, switch_time_ps=profile.switch_time_ps)

    print("\n[2/5] Loading TN material times...")
    tn_mt = load_tn_mt(data_root, switch_time_ps=profile.switch_time_ps)

    print("\n[3/5] Loading relaxation table...")
    relax_rows = load_relaxation_table(
        data_root,
        alt_root=profile.figures_dir if args.profile != "paper" else None,
    )

    print("\n[4/5] Loading fit parameters...")
    fit_params = load_fit_params(data_root, alt_root=profile.figures_dir if args.profile != "paper" else None)
    for eps, beta in sorted(fit_params.items()):
        print(f"  eps={eps:.0e}  beta={beta:.4f}  lambda={lam(eps):.4f}")

    # Build rows
    print("\n[5/5] Assembling CSV rows...")
    rows = []
    rows.extend(build_panel_a(meas_mt, tn_mt, downsample=args.downsample_a))
    rows.extend(build_panel_b(data_root, meas_mt, relax_rows, fit_params,
                              downsample=args.downsample_b))
    rows.extend(build_panels_cd(relax_rows))

    # ── shared metadata strings ────────────────────────────────────────────
    beta_mean = float(np.mean(list(fit_params.values())))
    ts        = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    desktop   = output.parent

    # ── 1. Combined figure4.csv (all panels, long format) ─────────────────
    fieldnames = ['panel', 'series', 'epsilon_au', 'lambda_au',
                  't_w_ps', 'ref_num', 'x', 'y']

    header_lines = [
        f"# figure4.csv — data for FIG. 4 (cavity-supercooled aging paper)",
        f"# Generated: {ts}",
        f"# Data root: {data_root}",
        f"# Source files:",
        f"#   panel a measured : direct integral dt/tau_0.1 from {data_root}/cavity_coupling_*/master_fskt_ref*.txt",
        f"#   panel a TN       : {data_root}/time_series_output/coupling_*_fictive_time_series.txt",
        f"#   panel b collapse : {data_root}/cavity_coupling_*/master_fskt_ref*.txt",
        f"#   panel b fit      : {data_root}/stretched_exponential_fit_results.txt",
        f"#   panel c/d tau    : {data_root}/relaxation_vs_reference_time_data.txt",
        f"#",
        f"# Column notes:",
        f"#   t_w_ps     : waiting time (ps) = lab_time - 200",
        f"#                (t_w=0 corresponds to coupling turn-on at lab t=200 ps)",
        f"#   x (panel a): time (ps) after coupling turn-on",
        f"#   x (panel b): h = xi(t_w+lag) - xi(t_w)  [material time increment]",
        f"#   x (panel c): lambda (a.u.)",
        f"#   x (panel d): t_w (ps, same convention as t_w_ps)",
        f"#   y (panel a): xi(t)  [material time]",
        f"#   y (panel b): Phi_k(h) = F(k, lag) / F(k, 0)",
        f"#   y (panel c/d): tau_tilde_s = tau(lambda, t_w) / tau(lambda=0, t_w)",
        f"#",
        f"# Physical constants:",
        f"#   omega_c = {OMEGA_C_CM_INV} cm^-1 = {OMEGA_C_AU:.6e} a.u.",
        f"#   lambda  = epsilon / omega_c",
        f"#   xi(t)   = integral dt/tau_0.1 with F/F0=0.1 criterion (measured and TN)",
        f"#",
        f"# Fit parameters (panel b, Phi_k(h) = exp(-h^beta)):",
        f"#   " + "  ".join(f"eps={e:.0e} beta={b:.4f}" for e, b in sorted(fit_params.items())),
        f"#   beta_mean used for collapse_fit curve: {beta_mean:.4f}",
    ]

    with open(output, 'w', newline='', encoding='utf-8') as fh:
        for hl in header_lines:
            fh.write(hl + '\n')
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} rows to {output}")

    # ── 2. Per-panel split files ───────────────────────────────────────────
    print("\nWriting split panel files...")

    def write_split(path, comment_lines, field_names, split_rows):
        with open(path, 'w', newline='', encoding='utf-8') as fh:
            for cl in comment_lines:
                fh.write(cl + '\n')
            w = csv.DictWriter(fh, fieldnames=field_names, extrasaction='ignore')
            w.writeheader()
            w.writerows(split_rows)
        print(f"  {path.name}: {len(split_rows)} rows")

    # Panel a — measured (solid lines, from F(k,t))
    pa_meas = [r for r in rows if r['panel'] == 'a' and r['series'] == 'measured_xi']
    write_split(
        desktop / 'figure4_panel_a_measured.csv',
        [
            f"# figure4_panel_a_measured.csv — FIG. 4(a) measured material time h_lambda(t)",
            f"# Source: direct integral dt/tau_0.1 from master F(k,t) (F/F0=0.1 criterion)",
            f"# Source files: {data_root}/cavity_coupling_*/master_fskt_ref*.txt",
            f"# Generated: {ts}",
            f"# Columns: lambda_au — coupling strength lambda=epsilon/omega_c (a.u.)",
            f"#           t_ps     — time (ps) after coupling turn-on (lab t - 200 ps)",
            f"#           h_lambda — material time xi(t) (dimensionless)",
            f"# omega_c = {OMEGA_C_CM_INV} cm^-1 = {OMEGA_C_AU:.6e} a.u.",
        ],
        ['lambda_au', 't_ps', 'h_lambda'],
        [{'lambda_au': r['lambda_au'], 't_ps': r['x'], 'h_lambda': r['y']} for r in pa_meas],
    )

    # Panel a — TN model (dashed lines)
    pa_tn = [r for r in rows if r['panel'] == 'a' and r['series'] == 'tn_xi']
    write_split(
        desktop / 'figure4_panel_a_tn.csv',
        [
            f"# figure4_panel_a_tn.csv — FIG. 4(a) TN model material time h_{{lambda,TN}}(t)",
            f"# Source: fictive relaxation time integration (Tool-Narayanaswamy formalism)",
            f"# Source files: {data_root}/time_series_output/coupling_*_fictive_time_series.txt",
            f"# Generated: {ts}",
            f"# Columns: lambda_au    — coupling strength lambda=epsilon/omega_c (a.u.)",
            f"#           t_ps        — time (ps) after coupling turn-on (lab t - 200 ps)",
            f"#           h_lambda_TN — TN material time xi_TN(t) = integral(1/tau_f) since turn-on",
            f"# omega_c  = {OMEGA_C_CM_INV} cm^-1 = {OMEGA_C_AU:.6e} a.u.",
        ],
        ['lambda_au', 't_ps', 'h_lambda_TN'],
        [{'lambda_au': r['lambda_au'], 't_ps': r['x'], 'h_lambda_TN': r['y']} for r in pa_tn],
    )

    # Panel b — collapse data + fit
    pb = [r for r in rows if r['panel'] == 'b']
    pb_data = [r for r in pb if r['series'] == 'collapse_data']
    pb_fit  = [r for r in pb if r['series'] == 'collapse_fit']
    write_split(
        desktop / 'figure4_panel_b.csv',
        [
            f"# figure4_panel_b.csv — FIG. 4(b) ISF data collapse Phi_k(h) and fit",
            f"# Source: master F(k,t) files + measured material time interpolation",
            f"# Source files: {data_root}/cavity_coupling_*/master_fskt_ref*.txt",
            f"#               {data_root}/material_time_all_couplings.txt",
            f"#               {data_root}/stretched_exponential_fit_results.txt",
            f"# Generated: {ts}",
            f"# Columns: series     — collapse_data | collapse_fit",
            f"#           lambda_au  — coupling strength (empty for fit curve)",
            f"#           t_w_ps     — waiting time t_w = lab_time - 200 ps",
            f"#           ref_num    — reference frame index (-1 for fit)",
            f"#           h          — material time increment xi(t_w+lag) - xi(t_w)",
            f"#           Phi_k      — normalized ISF = F(k, lag) / F(k, 0)",
            f"# Fit: Phi_k(h) = exp(-h^beta), beta_mean = {beta_mean:.4f}",
            f"# Per-coupling betas: " + "  ".join(f"lam={lam(e):.3f} b={b:.4f}"
                                                   for e, b in sorted(fit_params.items())),
        ],
        ['series', 'lambda_au', 't_w_ps', 'ref_num', 'h', 'Phi_k'],
        [{'series': r['series'], 'lambda_au': r['lambda_au'], 't_w_ps': r['t_w_ps'],
          'ref_num': r['ref_num'], 'h': r['x'], 'Phi_k': r['y']} for r in pb],
    )

    # Panel c — tau_tilde vs lambda
    pc = [r for r in rows if r['panel'] == 'c']
    write_split(
        desktop / 'figure4_panel_c.csv',
        [
            f"# figure4_panel_c.csv — FIG. 4(c) normalized relaxation time vs coupling",
            f"# tau_tilde_s = tau(lambda, t_w) / tau(lambda=0, t_w)",
            f"# Source file: {data_root}/relaxation_vs_reference_time_data.txt",
            f"# Generated: {ts}",
            f"# Columns: t_w_ps      — waiting time (ps) = lab_time - 200 (colorbar variable)",
            f"#           lambda_au   — coupling strength lambda (a.u.) [x axis]",
            f"#           tau_tilde_s — normalized relaxation time (dimensionless) [y axis]",
        ],
        ['t_w_ps', 'lambda_au', 'tau_tilde_s'],
        [{'t_w_ps': r['t_w_ps'], 'lambda_au': r['lambda_au'], 'tau_tilde_s': r['y']} for r in pc],
    )

    # Panel d — tau_tilde vs t_w
    pd_rows = [r for r in rows if r['panel'] == 'd']
    write_split(
        desktop / 'figure4_panel_d.csv',
        [
            f"# figure4_panel_d.csv — FIG. 4(d) normalized relaxation time vs waiting time",
            f"# tau_tilde_s = tau(lambda, t_w) / tau(lambda=0, t_w)",
            f"# Source file: {data_root}/relaxation_vs_reference_time_data.txt",
            f"# Generated: {ts}",
            f"# Columns: lambda_au   — coupling strength lambda (a.u.) [legend variable]",
            f"#           t_w_ps     — waiting time (ps) = lab_time - 200 [x axis]",
            f"#           tau_tilde_s — normalized relaxation time (dimensionless) [y axis]",
        ],
        ['lambda_au', 't_w_ps', 'tau_tilde_s'],
        [{'lambda_au': r['lambda_au'], 't_w_ps': r['t_w_ps'], 'tau_tilde_s': r['y']} for r in pd_rows],
    )

    print(f"\nSplit files written to {desktop}/")

    # ── Verify combined file ───────────────────────────────────────────────
    print("\n── Verification ─────────────────────────────────────────────────────")
    ok = verify_csv(
        output,
        fit_params,
        meas_mt,
        relax_rows,
        expected_lam_max=max(entry.axis_value for entry in profile.couplings),
    )
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
