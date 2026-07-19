#!/usr/bin/env python3
"""Convert averaged aging energy CSVs to fictive temperatures and plot them."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    from cavitymd.controllers.empirical import EmpiricalTemperatureData
except ModuleNotFoundError:
    from hoomd.cavitymd.controllers.empirical import EmpiricalTemperatureData

SWITCH_TIME_PS = 200.0
BASELINE_TARGET_K = 100.0
DEFAULT_ANALYSIS_DIR = Path("/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda/analysis")
DEFAULT_EMPIRICAL = Path(
    "/scratch/mh7373/projects/cav-hoomd/reproduction/calibration/potential_energy_vs_T.txt"
)

LAMBDA_CSV: dict[float, str] = {
    0.0: "avg_energies_lam0p0.csv",
    0.01: "avg_energies_lam0p01.csv",
    0.016667: "avg_energies_lam0p016667.csv",
    0.023333: "avg_energies_lam0p023333.csv",
    0.03: "avg_energies_lam0p03.csv",
}


def load_empirical(path: Path, energy_component: str) -> EmpiricalTemperatureData:
    """Load an empirical energy -> temperature converter."""
    if energy_component == "lj":
        return _load_lj_only_empirical(path)
    return EmpiricalTemperatureData(str(path), energy_component=energy_component, create_plots=False)


def _load_lj_only_empirical(path: Path) -> EmpiricalTemperatureData:
    """Build a converter from the calibration LJ column only (RF proxy)."""
    import tempfile

    lines = path.read_text(encoding="utf-8").splitlines()
    header_idx = next(
        index
        for index, line in enumerate(lines)
        if line.strip().startswith("temperature")
    )
    header = lines[header_idx].split("\t")
    col_index = {name.strip(): idx for idx, name in enumerate(header)}
    temp_idx = col_index["temperature"]
    lj_idx = col_index["lj_hartree"]

    tmp_lines = ["temperature\tlj_hartree\tcoulombic_hartree"]
    for line in lines[header_idx + 1 :]:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) <= max(temp_idx, lj_idx):
            continue
        tmp_lines.append(f"{parts[temp_idx].strip()}\t{parts[lj_idx].strip()}\t0.0")

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix="_lj_only_calibration.txt",
        delete=False,
        encoding="utf-8",
    ) as handle:
        handle.write("\n".join(tmp_lines) + "\n")
        tmp_path = handle.name

    return EmpiricalTemperatureData(tmp_path, energy_component="lj_coulombic", create_plots=False)


def shift_baseline_to_target(
    time_ps: np.ndarray,
    values: np.ndarray,
    t_max_ps: float,
    target_k: float,
) -> tuple[np.ndarray, float, float]:
    """Shift a series so its mean on ``[min(t), t_max]`` equals ``target_k``."""
    time_arr = np.asarray(time_ps, dtype=float)
    value_arr = np.asarray(values, dtype=float)
    mask = time_arr <= t_max_ps
    if not np.any(mask):
        raise ValueError(f"no samples with time <= {t_max_ps} ps")
    pre_mean = float(np.mean(value_arr[mask]))
    offset = target_k - pre_mean
    return value_arr + offset, pre_mean, offset


def convert_energy_arrays(
    time_ps: np.ndarray,
    harmonic_ha: np.ndarray,
    lj_coul_ha: np.ndarray,
    empirical_harmonic: EmpiricalTemperatureData,
    empirical_lj_coul: EmpiricalTemperatureData,
) -> dict[str, np.ndarray]:
    """Map energy time series to empirical fictive temperatures."""
    harmonic_T = np.array(
        [empirical_harmonic.calculate_systemic_temperature(float(energy)) for energy in harmonic_ha],
        dtype=float,
    )
    lj_T = np.array(
        [empirical_lj_coul.calculate_systemic_temperature(float(energy)) for energy in lj_coul_ha],
        dtype=float,
    )
    return {
        "time_ps": np.asarray(time_ps, dtype=float),
        "harmonic_fictive_K": harmonic_T,
        "lj_coul_fictive_K": lj_T,
    }


def apply_pre_switch_baseline_shift(
    data: dict[str, np.ndarray],
    t_max_ps: float,
    target_k: float,
) -> dict[str, np.ndarray | float]:
    """Shift harmonic and LJ/Coul Tf so each pre-switch mean is ``target_k``."""
    harm_shifted, harm_mean, harm_offset = shift_baseline_to_target(
        data["time_ps"],
        data["harmonic_fictive_K"],
        t_max_ps,
        target_k,
    )
    lj_shifted, lj_mean, lj_offset = shift_baseline_to_target(
        data["time_ps"],
        data["lj_coul_fictive_K"],
        t_max_ps,
        target_k,
    )
    return {
        **data,
        "harmonic_fictive_K": harm_shifted,
        "lj_coul_fictive_K": lj_shifted,
        "harmonic_pre_mean_K": harm_mean,
        "harmonic_offset_K": harm_offset,
        "lj_coul_pre_mean_K": lj_mean,
        "lj_coul_offset_K": lj_offset,
    }


def convert_csv(
    energy_csv: Path,
    empirical_harmonic: EmpiricalTemperatureData,
    empirical_lj_coul: EmpiricalTemperatureData,
) -> dict[str, np.ndarray]:
    table = np.genfromtxt(energy_csv, delimiter=",", names=True)
    converted = convert_energy_arrays(
        table["time_ps"],
        table["harmonic_ha"],
        table["lj_coul_ha"],
        empirical_harmonic,
        empirical_lj_coul,
    )
    return apply_pre_switch_baseline_shift(converted, SWITCH_TIME_PS, BASELINE_TARGET_K)


def save_fictive_csv(path: Path, data: dict[str, np.ndarray | float]) -> None:
    header = (
        "time_ps,harmonic_fictive_K,lj_coul_fictive_K "
        f"# pre0-{SWITCH_TIME_PS:.0f}ps means shifted to {BASELINE_TARGET_K:.0f} K; "
        f"harm_offset={data['harmonic_offset_K']:.6f} "
        f"lj_offset={data['lj_coul_offset_K']:.6f}"
    )
    columns = np.column_stack(
        [data["time_ps"], data["harmonic_fictive_K"], data["lj_coul_fictive_K"]]
    )
    np.savetxt(path, columns, delimiter=",", header=header, comments="")


def make_plots(output_dir: Path, results: dict[float, dict[str, np.ndarray | float]]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for lam, data in sorted(results.items()):
        axes[0].plot(data["time_ps"], data["harmonic_fictive_K"], label=f"λ={lam:g}", alpha=0.85)
        axes[1].plot(data["time_ps"], data["lj_coul_fictive_K"], label=f"λ={lam:g}", alpha=0.85)
    for ax in axes:
        ax.axhline(BASELINE_TARGET_K, color="k", linestyle=":", alpha=0.4)
        ax.axvline(SWITCH_TIME_PS, color="k", linestyle="--", alpha=0.4)
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Harmonic Tf (K)")
    axes[1].set_ylabel("LJ + Coul Tf (K)")
    axes[1].set_xlabel("Time (ps)")
    axes[0].set_title(
        "Replica-averaged fictive temperatures — aging campaign "
        f"(each curve shifted so mean[0,{SWITCH_TIME_PS:.0f} ps] = {BASELINE_TARGET_K:.0f} K)"
    )
    axes[0].legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(output_dir / "avg_fictive_temperature_vs_time.png", dpi=150)
    plt.close(fig)

    fig2, axes2 = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    zoom = (150.0, 250.0)
    for lam, data in sorted(results.items()):
        mask = (data["time_ps"] >= zoom[0]) & (data["time_ps"] <= zoom[1])
        axes2[0].plot(data["time_ps"][mask], data["harmonic_fictive_K"][mask], alpha=0.85)
        axes2[1].plot(data["time_ps"][mask], data["lj_coul_fictive_K"][mask], alpha=0.85)
    axes2[0].set_ylabel("Harmonic Tf (K)")
    axes2[1].set_ylabel("LJ + Coul Tf (K)")
    axes2[1].set_xlabel("Time (ps)")
    fig2.suptitle(f"Fictive-temperature switch window ({zoom[0]:.0f}–{zoom[1]:.0f} ps)")
    fig2.tight_layout()
    fig2.savefig(output_dir / "avg_fictive_temperature_switch_zoom.png", dpi=150)
    plt.close(fig2)


def print_summary(results: dict[float, dict[str, np.ndarray | float]]) -> None:
    print("\n=== Pre-switch baseline shift (0–200 ps mean → 100 K) ===")
    for lam, data in sorted(results.items()):
        print(
            f"λ={lam:8.5f}  harm: mean={data['harmonic_pre_mean_K']:.2f} K  "
            f"offset={data['harmonic_offset_K']:+.2f} K  "
            f"lj: mean={data['lj_coul_pre_mean_K']:.2f} K  "
            f"offset={data['lj_coul_offset_K']:+.2f} K"
        )
    print("\n=== Switch response: ΔTf from t=150→250 ps ===")
    for lam, data in sorted(results.items()):
        post_idx = int(np.argmin(np.abs(data["time_ps"] - 250.0)))
        pre_idx = int(np.argmin(np.abs(data["time_ps"] - 150.0)))
        delta_harm = data["harmonic_fictive_K"][post_idx] - data["harmonic_fictive_K"][pre_idx]
        delta_lj = data["lj_coul_fictive_K"][post_idx] - data["lj_coul_fictive_K"][pre_idx]
        print(
            f"λ={lam:8.5f}  ΔT_harm={delta_harm:+.2f} K  ΔT_lj={delta_lj:+.2f} K  [OK]"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS_DIR)
    parser.add_argument("--empirical-data", type=Path, default=DEFAULT_EMPIRICAL)
    parser.add_argument(
        "--structural-energy-component",
        choices=("lj_coulombic", "lj"),
        default="lj_coulombic",
        help="Calibration channel for the lj_coul_ha CSV column (use lj for RF campaigns)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    analysis_dir = args.analysis_dir
    analysis_dir.mkdir(parents=True, exist_ok=True)
    empirical_harmonic = load_empirical(args.empirical_data, "harmonic")
    empirical_lj_coul = load_empirical(
        args.empirical_data,
        args.structural_energy_component,
    )

    results: dict[float, dict[str, np.ndarray | float]] = {}
    for lam, csv_name in LAMBDA_CSV.items():
        energy_csv = analysis_dir / csv_name
        if not energy_csv.is_file():
            print(f"λ={lam}: missing {csv_name}", flush=True)
            continue
        converted = convert_csv(energy_csv, empirical_harmonic, empirical_lj_coul)
        out_csv = analysis_dir / csv_name.replace("avg_energies", "avg_fictive_temps")
        save_fictive_csv(out_csv, converted)
        print(
            f"λ={lam}: wrote {out_csv.name} "
            f"(harm_offset={converted['harmonic_offset_K']:.3f} K, "
            f"lj_offset={converted['lj_coul_offset_K']:.3f} K)",
            flush=True,
        )
        results[lam] = converted

    if not results:
        raise SystemExit("No energy CSVs found to convert.")
    make_plots(analysis_dir, results)
    print_summary(results)
    print(f"Saved fictive-temperature CSVs and plots under {analysis_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
