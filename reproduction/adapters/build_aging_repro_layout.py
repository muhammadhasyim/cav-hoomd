#!/usr/bin/env python3
"""Build paper-format staging layout for aging_weak_lambda reproduction scripts."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

REPRO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPRO_ROOT / "shared"))

from dataset_profile import (  # noqa: E402
    DatasetProfile,
    add_profile_args,
    load_profile,
    setup_profile,
    staged_coupling_dir,
)
from plot_fkt_analysis import find_relaxation_time, read_fkt_data  # noqa: E402

MIN_FKT_AGE_SEC = 120
MIN_FKT_DATA_LINES = 3
FICTIVE_HEADER_LINES = 12


def _require_master_fskt(run_dir: Path) -> None:
    masters = [
        p for p in run_dir.glob("master_fskt_ref*.txt")
        if "_sample_counts" not in p.name
    ]
    if not masters:
        raise FileNotFoundError(f"No master F(k,t) files in {run_dir}")


def _preflight_fkt_replicas(run_dir: Path, quarantine_name: str = "fkt_incomplete_quarantine") -> None:
    """Skip or quarantine in-flight replica FKT files (same policy as careful postproc)."""
    quarantine = run_dir / quarantine_name
    quarantine.mkdir(parents=True, exist_ok=True)
    now = time.time()
    kept = skipped = 0

    for path in sorted(run_dir.glob("prod-*_fkt_ref_*.txt")):
        try:
            age = now - path.stat().st_mtime
            if age < MIN_FKT_AGE_SEC:
                skipped += 1
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            if path.exists():
                shutil.move(str(path), str(quarantine / path.name))
            skipped += 1
            continue

        lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
        if len(lines) < MIN_FKT_DATA_LINES:
            shutil.move(str(path), str(quarantine / path.name))
            skipped += 1
            continue
        kept += 1

    print(f"  FKT preflight {run_dir.name}: kept {kept}, skipped/quarantined {skipped}")


def _write_symlinks(profile: DatasetProfile, dry_run: bool = False) -> list[Path]:
    profile.staged_root.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []

    for entry in profile.couplings:
        target = profile.data_root / entry.run_dir
        link = profile.staged_root / entry.staged_dir_name
        if not target.is_dir():
            raise FileNotFoundError(f"Run directory missing: {target}")
        _require_master_fskt(target)
        _preflight_fkt_replicas(target)

        rel_target = Path("..") / ".." / Path(entry.run_dir)
        if dry_run:
            print(f"[dry-run] symlink {link} -> {rel_target}")
            created.append(link)
            continue

        if link.is_symlink() or link.exists():
            if link.is_symlink():
                link.unlink()
            elif link.is_dir():
                raise FileExistsError(f"Refusing to replace non-symlink directory: {link}")

        link.symlink_to(rel_target, target_is_directory=True)
        print(f"  symlink {link.name} -> {rel_target}")
        created.append(link)

    return created


def _load_energy_csv(path: Path) -> dict[str, np.ndarray]:
    time_ps: list[float] = []
    harmonic: list[float] = []
    lj_coul: list[float] = []
    total: list[float] = []
    kinetic: list[float] | None = None

    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        has_kinetic = "kinetic_temp_K" in fieldnames
        if has_kinetic:
            kinetic = []
        for row in reader:
            time_ps.append(float(row["time_ps"]))
            harmonic.append(float(row["harmonic_ha"]))
            lj_coul.append(float(row["lj_coul_ha"]))
            total.append(float(row["system_total_ha"]))
            if has_kinetic and kinetic is not None:
                kinetic.append(float(row["kinetic_temp_K"]))

    result = {
        "time_ps": np.asarray(time_ps, dtype=float),
        "harmonic_ha": np.asarray(harmonic, dtype=float),
        "lj_coul_ha": np.asarray(lj_coul, dtype=float),
        "total_ha": np.asarray(total, dtype=float),
    }
    if kinetic is not None:
        result["kinetic_temp_K"] = np.asarray(kinetic, dtype=float)
    return result


def _load_fictive_csv(path: Path) -> dict[str, np.ndarray]:
    time_ps: list[float] = []
    harmonic: list[float] = []
    lj_coul: list[float] = []

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            if line.lower().startswith("time_ps"):
                continue
            parts = line.split(",")
            time_ps.append(float(parts[0]))
            harmonic.append(float(parts[1]))
            lj_coul.append(float(parts[2]))

    return {
        "time_ps": np.asarray(time_ps, dtype=float),
        "harmonic_fictive_K": np.asarray(harmonic, dtype=float),
        "lj_coul_fictive_K": np.asarray(lj_coul, dtype=float),
    }


def _write_potential_energy(path: Path, data: dict[str, np.ndarray]) -> None:
    """Write the paper 5-column layout: time, harmonic, lj, coulombic, total.

    The aging CSVs provide a combined LJ+Coulomb component, so it is written in
    the ``lj`` column with a zero ``coulombic`` column; downstream scripts sum
    the two (``U_lj + U_coul``) to recover the combined term while reading the
    system total from the final column.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# Averaged potential energy time series (generated for reproduction staging)\n")
        fh.write("# Columns: time_ps harmonic_ha lj_ha coul_ha total_ha\n")
        fh.write("time_ps harmonic_ha lj_ha coul_ha total_ha\n")
        for idx in range(len(data["time_ps"])):
            lj = data["lj_coul_ha"][idx]
            fh.write(
                f"{data['time_ps'][idx]:.6e} "
                f"{data['harmonic_ha'][idx]:.6e} "
                f"{lj:.6e} "
                f"0.000000e+00 "
                f"{data['total_ha'][idx]:.6e}\n"
            )


def _write_kinetic_temperature(path: Path, data: dict[str, np.ndarray]) -> None:
    """Write replica-averaged kinetic temperature time series."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# Averaged kinetic temperature time series (generated for reproduction staging)\n")
        fh.write("# Source: temperatures/kinetic from observables HDF5 via average_aging_energies.py\n")
        fh.write("# Columns: time_ps kinetic_temp_K\n")
        fh.write("time_ps kinetic_temp_K\n")
        for idx in range(len(data["time_ps"])):
            fh.write(
                f"{data['time_ps'][idx]:.6e} "
                f"{data['kinetic_temp_K'][idx]:.6e}\n"
            )


def _write_fictive_temperatures(path: Path, data: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# Averaged fictive temperatures (generated for reproduction staging)\n")
        fh.write("# Columns: time_ps placeholder harmonic_fictive_K\n")
        fh.write("time_ps placeholder harmonic_fictive_K\n")
        for idx in range(len(data["time_ps"])):
            fh.write(
                f"{data['time_ps'][idx]:.6e} "
                f"0.000000e+00 "
                f"{data['harmonic_fictive_K'][idx]:.6e}\n"
            )


def _relaxation_model(profile: DatasetProfile):
    from hoomd.cavitymd.controllers.diffeq import RelaxationTimeModel

    if profile.relaxation_equilibrium is None:
        raise ValueError("Profile missing relaxation_equilibrium path")
    return RelaxationTimeModel(str(profile.relaxation_equilibrium))


def _write_fictive_time_series(
    path: Path,
    fictive: dict[str, np.ndarray],
    model,
    tag: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    time_ps = fictive["time_ps"]
    t_struct = fictive["lj_coul_fictive_K"]
    tau = np.array([float(model.get_relaxation_time(float(temp))) for temp in t_struct], dtype=float)

    header = f"""# Fictive time series for coupling_{tag}
# Generated by build_aging_repro_layout.py
# Source: avg_fictive_temps + RelaxationTimeModel
# Columns:
#   time_ps (col 0)
#   structural_fictive_temperature_K (col 1)
#   fictive_relaxation_time_ps (col 2)
#   predicted_tau_alpha_ps (col 3)
#
# time_ps structural_fictive_temperature_K fictive_relaxation_time_ps predicted_tau_alpha_ps
"""
    lines = header.splitlines()
    while len(lines) < FICTIVE_HEADER_LINES:
        lines.append("#")
    lines[FICTIVE_HEADER_LINES - 1] = (
        "time_ps structural_fictive_temperature_K "
        "fictive_relaxation_time_ps predicted_tau_alpha_ps"
    )

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
        for idx in range(len(time_ps)):
            fh.write(
                f"{time_ps[idx]:.6e} "
                f"{t_struct[idx]:.6e} "
                f"{tau[idx]:.6e} "
                f"{tau[idx]:.6e}\n"
            )


def _write_time_series_output(profile: DatasetProfile, dry_run: bool = False) -> list[Path]:
    if profile.analysis_dir is None:
        raise ValueError("Profile missing analysis_dir")

    model = _relaxation_model(profile)
    out_dir = profile.time_series_dir
    written: list[Path] = []

    for entry in profile.couplings:
        if entry.analysis_tag is None:
            raise ValueError(f"Coupling {entry.id} missing analysis_tag")

        energy_csv = profile.analysis_dir / f"avg_energies_{entry.analysis_tag}.csv"
        fictive_csv = profile.analysis_dir / f"avg_fictive_temps_{entry.analysis_tag}.csv"
        if not energy_csv.is_file():
            raise FileNotFoundError(f"Missing energy CSV: {energy_csv}")
        if not fictive_csv.is_file():
            raise FileNotFoundError(f"Missing fictive CSV: {fictive_csv}")

        tag = entry.effective_tag
        pe_path = out_dir / f"coupling_{tag}_averaged_potential_energy.txt"
        kt_path = out_dir / f"coupling_{tag}_kinetic_temperature.txt"
        ft_path = out_dir / f"coupling_{tag}_averaged_fictive_temperatures.txt"
        fts_path = out_dir / f"coupling_{tag}_fictive_time_series.txt"

        if dry_run:
            print(f"[dry-run] would write {pe_path.name}, {kt_path.name}, {ft_path.name}, {fts_path.name}")
            written.extend([pe_path, kt_path, ft_path, fts_path])
            continue

        energy = _load_energy_csv(energy_csv)
        fictive = _load_fictive_csv(fictive_csv)
        _write_potential_energy(pe_path, energy)
        if "kinetic_temp_K" in energy:
            _write_kinetic_temperature(kt_path, energy)
        else:
            print(f"  WARNING: {energy_csv.name} lacks kinetic_temp_K; skipping kinetic staging")
        _write_fictive_temperatures(ft_path, fictive)
        _write_fictive_time_series(fts_path, fictive, model, tag)
        print(f"  time series for {entry.id} ({tag})")
        written.extend([pe_path, ft_path, fts_path])
        if "kinetic_temp_K" in energy:
            written.append(kt_path)

    return written


def _parse_ref_number(path: Path) -> int | None:
    match = re.search(r"ref[_]?(\d+)", path.name)
    if not match:
        return None
    return int(match.group(1))


def _write_relaxation_table(profile: DatasetProfile, dry_run: bool = False) -> Path:
    out_path = profile.staged_root / "relaxation_vs_reference_time_data.txt"
    rows: list[tuple[str, float, str, int, float, float, float]] = []

    for entry in profile.couplings:
        coupling_dir = (
            profile.data_root / entry.run_dir
            if dry_run
            else staged_coupling_dir(profile, entry.id)
        )
        coupling_name = entry.staged_dir_name
        axis_value = entry.axis_value
        label = f"{axis_value:g}" if profile.axis == "lambda" else f"{axis_value:.6e}"

        ref_files = sorted(
            p for p in coupling_dir.glob("master_fskt_ref*.txt")
            if "_sample_counts" not in p.name
        )
        if not ref_files:
            raise FileNotFoundError(f"No master F(k,t) in staged dir {coupling_dir}")

        norm_value = None
        ref0_candidates = [p for p in ref_files if _parse_ref_number(p) == 0]
        if ref0_candidates:
            t0, f0, _counts0 = read_fkt_data(ref0_candidates[0])
            if t0 is not None and f0 is not None:
                nonzero = f0[f0 != 0]
                if len(nonzero):
                    norm_value = float(nonzero[0])

        for ref_path in ref_files:
            ref_num = _parse_ref_number(ref_path)
            if ref_num is None:
                continue
            time_arr, fkt_arr, _counts = read_fkt_data(ref_path)
            if time_arr is None or fkt_arr is None:
                continue
            tau = find_relaxation_time(
                time_arr, fkt_arr, target_value=0.1, normalization_value=norm_value
            )
            if np.isnan(tau):
                continue
            ref_time_ps = float(ref_num * 200)
            rows.append(
                (
                    coupling_name,
                    axis_value,
                    label,
                    ref_num,
                    ref_time_ps,
                    float(tau),
                    float(norm_value) if norm_value is not None else float("nan"),
                )
            )

    if dry_run:
        print(f"[dry-run] would write relaxation table with {len(rows)} rows -> {out_path}")
        return out_path
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("# Relaxation Time vs Reference Time Data\n")
        fh.write("# Generated by build_aging_repro_layout.py (F/F0 = 0.1 criterion)\n")
        fh.write(
            "# Columns: coupling_name, coupling_value, coupling_label, "
            "ref_num, ref_time_ps, relaxation_time_ps, normalization_value\n#\n"
        )
        fh.write(
            f"{'coupling_name':<40} {'coupling_value':<15} {'coupling_label':<20} "
            f"{'ref_num':<8} {'ref_time_ps':<12} {'relaxation_time_ps':<18} "
            f"{'normalization_value':<18}\n"
        )
        for row in rows:
            fh.write(
                f"{row[0]:<40} {row[1]:<15.6e} {row[2]:<20} "
                f"{row[3]:<8} {row[4]:<12.1f} {row[5]:<18.6f} {row[6]:<18.6e}\n"
            )

    print(f"  relaxation table: {len(rows)} rows -> {out_path.name}")
    return out_path


def _write_material_time_all_couplings(profile: DatasetProfile, dry_run: bool = False) -> Path | None:
    """Optional: run corrected material-time reconstructor once per coupling."""
    sys.path.insert(0, str(REPRO_ROOT / "figure4_material_time"))
    from run_corrected_analysis import CorrectedMaterialTimeAnalyzer  # noqa: E402

    analyzer = CorrectedMaterialTimeAnalyzer(
        criterion_value=0.1,
        alpha_smooth=0.01,
        n_grid_points=50,
        verbose=False,
    )

    series: list[tuple[float, np.ndarray, np.ndarray]] = []
    for entry in profile.couplings:
        coupling_dir = staged_coupling_dir(profile, entry.id)
        if dry_run:
            print(f"[dry-run] material time for {entry.id}")
            continue
        fkt_files = analyzer.load_fkt_data(str(coupling_dir))
        normalized = analyzer.normalize_fkt_data(fkt_files)
        time_grid, xi_grid = analyzer.calculate_material_time(normalized)
        series.append((entry.axis_value, time_grid, xi_grid))

    out_path = profile.staged_root / "material_time_all_couplings.txt"
    if dry_run:
        print(f"[dry-run] would write {out_path}")
        return out_path

    if not series:
        print("  skipping material_time_all_couplings.txt (no series computed)")
        return None

    common_t = series[0][1]
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("# Material time functions for all couplings (staging adapter output)\n")
        fh.write("# Criterion: F(k,t) = 0.1, corrected least-squares reconstruction\n")
        fh.write("# Each column pair is (time_ps, xi) for one coupling.\n#\n")
        header_cols: list[str] = []
        for axis_value, _, _ in series:
            tag = f"{axis_value:g}".replace(".", "p")
            header_cols.extend([f"t_{tag}_ps", f"xi_{tag}"])
        fh.write(" ".join(header_cols) + "\n")
        for idx in range(len(common_t)):
            parts: list[str] = []
            for _, time_grid, xi in series:
                parts.append(f"{time_grid[idx]:.6e}")
                parts.append(f"{xi[idx]:.6e}")
            fh.write(" ".join(parts) + "\n")

    print(f"  wrote {out_path.name}")
    return out_path


def build_layout(
    profile: DatasetProfile,
    *,
    dry_run: bool = False,
    skip_material_time: bool = False,
) -> None:
    print(f"Building reproduction layout for profile={profile.name}")
    print(f"  staged_root: {profile.staged_root}")

    print("\n[1/4] Flat coupling symlinks")
    _write_symlinks(profile, dry_run=dry_run)

    print("\n[2/4] time_series_output/")
    _write_time_series_output(profile, dry_run=dry_run)

    print("\n[3/4] relaxation_vs_reference_time_data.txt")
    _write_relaxation_table(profile, dry_run=dry_run)

    if skip_material_time:
        print("\n[4/4] Skipping material_time_all_couplings.txt")
    else:
        print("\n[4/4] material_time_all_couplings.txt (optional)")
        try:
            _write_material_time_all_couplings(profile, dry_run=dry_run)
        except Exception as exc:
            print(f"  WARNING: material time export failed: {exc}")

    print("\nStaging complete.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_profile_args(parser, default="aging_weak_lambda")
    parser.add_argument("--dry-run", action="store_true", help="Print actions without writing files")
    parser.add_argument(
        "--skip-material-time",
        action="store_true",
        help="Do not emit material_time_all_couplings.txt",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    profile = setup_profile(args, default="aging_weak_lambda")
    build_layout(
        profile,
        dry_run=args.dry_run,
        skip_material_time=args.skip_material_time,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
