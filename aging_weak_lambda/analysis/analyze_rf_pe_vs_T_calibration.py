#!/usr/bin/env python3
"""Analyze reaction-field PE-vs-T calibration HDF5 outputs."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

DEFAULT_ROOT = Path(__file__).resolve().parent / "pe_vs_T_calib_rf"
DEFAULT_DISCARD_PS = 1000.0
DEFAULT_PROD_END_PS = 2000.0
DEFAULT_BLOCK_PS = 1.0
DEFAULT_SAMPLE_PS = 0.01
TEMP_DIR_RE = re.compile(r"^temp_(?P<temp>\d+(?:\.\d+)?)K$")


def read_hdf5_sample_period_ps(h5_path: Path) -> float | None:
    """Read ``output_period_ps`` from observables HDF5 attrs when present."""
    with h5py.File(h5_path, "r") as handle:
        if "output_period_ps" in handle.attrs:
            return float(handle.attrs["output_period_ps"])
    return None


def resolve_sample_period_ps(
    sample_ps: float | None,
    *,
    hdf5_sample_ps: float | None,
) -> float:
    """
    Choose the HDF5 sample period used for block-size conversion.

    When ``output_period_ps`` is stored in the file, it overrides CLI defaults
    (legacy scripts assumed 0.01 ps while RF calib HDF5 uses 1.0 ps).
    """
    if hdf5_sample_ps is not None:
        return float(hdf5_sample_ps)
    if sample_ps is not None:
        return float(sample_ps)
    return DEFAULT_SAMPLE_PS


@dataclass(frozen=True)
class ProductionWindow:
    """Energy samples in the production window after equilibration discard."""

    time_ps: np.ndarray
    harmonic_ha: np.ndarray
    lj_ha: np.ndarray
    rf_coul_ha: np.ndarray
    t_kin_k: np.ndarray
    n_samples: int


@dataclass(frozen=True)
class TemperatureSummary:
    """Block-averaged calibration point for one temperature."""

    temperature_k: float
    harmonic_mean_ha: float
    harmonic_std_ha: float
    lj_mean_ha: float
    lj_std_ha: float
    coulombic_mean_ha: float
    coulombic_std_ha: float
    total_mean_ha: float
    total_std_ha: float
    t_kin_mean_k: float
    n_samples: int
    n_blocks: int
    t_end_ps: float


def extract_temperature_k_from_dirname(name: str) -> float:
    """Parse ``temp_300K`` style directory names."""
    match = TEMP_DIR_RE.match(name)
    if match is None:
        raise ValueError(f"cannot parse temperature from directory name: {name!r}")
    return float(match.group("temp"))


def block_average_series(
    values: np.ndarray,
    *,
    block_size: int,
) -> tuple[float, float, int]:
    """
    Compute mean and standard deviation of block means.

    Parameters
    ----------
    values
        1D array of samples (must be an integer multiple of block_size).
    block_size
        Number of consecutive samples per block.

    Returns
    -------
    mean
        Mean of block means.
    std
        Standard deviation of block means (0 for a single block).
    n_blocks
        Number of blocks used.
    """
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        raise ValueError("values must not be empty")
    if arr.size % block_size != 0:
        raise ValueError(
            f"values length {arr.size} must be a multiple of block_size {block_size}"
        )
    blocked = arr.reshape(-1, block_size)
    block_means = blocked.mean(axis=1)
    n_blocks = block_means.size
    mean = float(block_means.mean())
    std = float(block_means.std(ddof=1)) if n_blocks > 1 else 0.0
    return mean, std, n_blocks


def _trim_valid_time_series(
    time_ps: np.ndarray,
    *arrays: np.ndarray,
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Remove leading/trailing zero-padding from observable writer output."""
    time_arr = np.asarray(time_ps, dtype=float)
    if time_arr.size == 0:
        raise ValueError("empty time array")
    valid = time_arr > 0.0
    if not np.any(valid):
        raise ValueError("no positive simulation times in HDF5")
    first = int(np.argmax(valid))
    last = int(len(time_arr) - np.argmax(valid[::-1]))
    trimmed_time = time_arr[first:last]
    trimmed_arrays = [np.asarray(array, dtype=float)[first:last] for array in arrays]
    return trimmed_time, trimmed_arrays


def load_production_window(
    h5_path: Path,
    *,
    discard_ps: float,
    prod_end_ps: float,
) -> ProductionWindow:
    """Load harmonic, LJ, and RF Coulomb energies in the production window."""
    with h5py.File(h5_path, "r") as handle:
        time_raw = handle["time"][:]
        harmonic_raw = handle["energies/harmonic"][:]
        lj_raw = handle["energies/lj"][:]
        rf_raw = handle["energies/ewald_short"][:]
        if "temperature" in handle["energies"]:
            temp_raw = handle["energies/temperature"][:]
        elif "temperatures/kinetic" in handle:
            temp_raw = handle["temperatures/kinetic"][:]
        else:
            temp_raw = np.zeros_like(time_raw)

    time_ps, arrays = _trim_valid_time_series(
        time_raw,
        harmonic_raw,
        lj_raw,
        rf_raw,
        temp_raw,
    )
    harmonic, lj, rf_coul, t_kin = arrays
    mask = (time_ps > discard_ps) & (time_ps <= prod_end_ps + 1e-9)
    if not np.any(mask):
        raise ValueError(
            f"no samples in production window [{discard_ps}, {prod_end_ps}] ps "
            f"for {h5_path}"
        )
    return ProductionWindow(
        time_ps=time_ps[mask],
        harmonic_ha=harmonic[mask],
        lj_ha=lj[mask],
        rf_coul_ha=rf_coul[mask],
        t_kin_k=t_kin[mask],
        n_samples=int(mask.sum()),
    )


def summarize_temperature_run(
    temp_dir: Path,
    *,
    discard_ps: float,
    prod_end_ps: float,
    block_ps: float,
    sample_ps: float | None,
) -> TemperatureSummary:
    """Block-average one temperature directory."""
    h5_path = temp_dir / "no_cavity" / "observables_replica_0.h5"
    if not h5_path.is_file():
        raise FileNotFoundError(h5_path)
    effective_sample_ps = resolve_sample_period_ps(
        sample_ps,
        hdf5_sample_ps=read_hdf5_sample_period_ps(h5_path),
    )
    window = load_production_window(
        h5_path,
        discard_ps=discard_ps,
        prod_end_ps=prod_end_ps,
    )
    block_size = int(round(block_ps / effective_sample_ps))
    if block_size < 1:
        raise ValueError(
            f"block_ps={block_ps} must be >= sample_ps={effective_sample_ps}"
        )

    n_use = (window.n_samples // block_size) * block_size
    if n_use == 0:
        raise ValueError(f"insufficient production samples in {h5_path}")

    harmonic_mean, harmonic_std, n_blocks = block_average_series(
        window.harmonic_ha[:n_use],
        block_size=block_size,
    )
    lj_mean, lj_std, _ = block_average_series(window.lj_ha[:n_use], block_size=block_size)
    coul_mean, coul_std, _ = block_average_series(
        window.rf_coul_ha[:n_use],
        block_size=block_size,
    )
    total = window.harmonic_ha[:n_use] + window.lj_ha[:n_use] + window.rf_coul_ha[:n_use]
    total_mean, total_std, _ = block_average_series(total, block_size=block_size)
    t_kin_mean = float(np.mean(window.t_kin_k[:n_use])) if np.any(window.t_kin_k[:n_use]) else 0.0

    return TemperatureSummary(
        temperature_k=extract_temperature_k_from_dirname(temp_dir.name),
        harmonic_mean_ha=harmonic_mean,
        harmonic_std_ha=harmonic_std,
        lj_mean_ha=lj_mean,
        lj_std_ha=lj_std,
        coulombic_mean_ha=coul_mean,
        coulombic_std_ha=coul_std,
        total_mean_ha=total_mean,
        total_std_ha=total_std,
        t_kin_mean_k=t_kin_mean,
        n_samples=n_use,
        n_blocks=n_blocks,
        t_end_ps=float(window.time_ps[-1]),
    )


def discover_temperature_dirs(root: Path) -> list[Path]:
    """Return sorted temp_*K directories under root."""
    dirs = [path for path in root.iterdir() if path.is_dir() and TEMP_DIR_RE.match(path.name)]
    return sorted(dirs, key=lambda path: extract_temperature_k_from_dirname(path.name), reverse=True)


def write_calibration_table(
    summaries: list[TemperatureSummary],
    output_path: Path,
    *,
    block_ps: float = DEFAULT_BLOCK_PS,
    sample_ps: float | None = None,
) -> None:
    """Write tab-separated calibration table compatible with EmpiricalTemperatureData."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sample_label = f"{sample_ps:g}" if sample_ps is not None else "HDF5 output_period_ps"
    lines = [
        "# Potential Energy Components vs Temperature (reaction field, eps_rf=80)",
        "# All energies in Hartree (atomic units)",
        "# Production window: last 1 ns after 1 ns equilibration discard",
        f"# Block averaging: {block_ps:g} ps blocks from {sample_label} ps HDF5 samples",
        "# Columns:",
        "#   temperature: set-point temperature in Kelvin",
        "#   harmonic_hartree: block-averaged harmonic bond energy",
        "#   harmonic_energy_std: std of 1 ps block means",
        "#   lj_hartree: block-averaged Lennard-Jones energy",
        "#   lj_energy_std: std of 1 ps block means",
        "#   coulombic_hartree: block-averaged reaction-field Coulomb energy (ewald_short)",
        "#   coulombic_energy_std: std of 1 ps block means",
        "#   total_potential_energy_hartree: harmonic + lj + coulombic",
        "#   total_potential_energy_std: std of 1 ps block means",
        "#",
        "temperature\tharmonic_hartree\tharmonic_energy_std\tlj_hartree\t"
        "lj_energy_std\tcoulombic_hartree\tcoulombic_energy_std\t"
        "total_potential_energy_hartree\ttotal_potential_energy_std",
    ]
    for summary in sorted(summaries, key=lambda item: item.temperature_k, reverse=True):
        lines.append(
            f"{summary.temperature_k:18.2f}\t"
            f"{summary.harmonic_mean_ha:.8e}\t{summary.harmonic_std_ha:.8e}\t"
            f"{summary.lj_mean_ha:.8e}\t{summary.lj_std_ha:.8e}\t"
            f"{summary.coulombic_mean_ha:.8e}\t{summary.coulombic_std_ha:.8e}\t"
            f"{summary.total_mean_ha:.8e}\t{summary.total_std_ha:.8e}"
        )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_analyze_log(
    summaries: list[TemperatureSummary],
    log_path: Path,
    *,
    discard_ps: float,
    prod_end_ps: float,
) -> None:
    """Write human-readable per-temperature diagnostics."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"discard_ps={discard_ps} prod_end_ps={prod_end_ps}",
        "T_set(K)  T_kin(K)  harm(Ha)  lj(Ha)  rf_coul(Ha)  total(Ha)  n  t_end(ps)",
    ]
    for summary in sorted(summaries, key=lambda item: item.temperature_k, reverse=True):
        lines.append(
            f"T_set={summary.temperature_k:.0f}K  "
            f"T_kin={summary.t_kin_mean_k:.2f}K  "
            f"harm={summary.harmonic_mean_ha:.6f}  "
            f"lj={summary.lj_mean_ha:.6f}  "
            f"rf_coul={summary.coulombic_mean_ha:.6f}  "
            f"total={summary.total_mean_ha:.6f}  "
            f"n={summary.n_samples}  t_end={summary.t_end_ps:.1f}ps"
        )
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def analyze_calibration_root(
    root: Path,
    *,
    discard_ps: float = DEFAULT_DISCARD_PS,
    prod_end_ps: float = DEFAULT_PROD_END_PS,
    block_ps: float = DEFAULT_BLOCK_PS,
    sample_ps: float | None = None,
) -> list[TemperatureSummary]:
    """Analyze all temp_*K directories under root."""
    summaries: list[TemperatureSummary] = []
    for temp_dir in discover_temperature_dirs(root):
        summaries.append(
            summarize_temperature_run(
                temp_dir,
                discard_ps=discard_ps,
                prod_end_ps=prod_end_ps,
                block_ps=block_ps,
                sample_ps=sample_ps,
            )
        )
    return summaries


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="Calibration output root (default: pe_vs_T_calib_rf)",
    )
    parser.add_argument(
        "--discard-ps",
        type=float,
        default=DEFAULT_DISCARD_PS,
        help="Equilibration discard in ps (default: 1000)",
    )
    parser.add_argument(
        "--prod-end-ps",
        type=float,
        default=DEFAULT_PROD_END_PS,
        help="End of production window in ps (default: 2000)",
    )
    parser.add_argument(
        "--block-ps",
        type=float,
        default=DEFAULT_BLOCK_PS,
        help="Block size for block averaging in ps (default: 1.0)",
    )
    parser.add_argument(
        "--sample-ps",
        type=float,
        default=None,
        help=(
            "HDF5 sample period in ps (default: read output_period_ps from each HDF5; "
            f"fallback {DEFAULT_SAMPLE_PS} ps if missing)"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output table path (default: root/potential_energy_vs_T_rf.txt)",
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=None,
        help="Summary log path (default: root/analyze.log)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root
    output = args.output or (root / "potential_energy_vs_T_rf.txt")
    log_path = args.log or (root / "analyze.log")
    summaries = analyze_calibration_root(
        root,
        discard_ps=args.discard_ps,
        prod_end_ps=args.prod_end_ps,
        block_ps=args.block_ps,
        sample_ps=args.sample_ps,
    )
    if not summaries:
        raise SystemExit(f"no temp_*K directories found under {root}")
    effective_sample_ps = args.sample_ps
    if effective_sample_ps is None:
        temp_dirs = discover_temperature_dirs(root)
        if temp_dirs:
            h5_path = temp_dirs[0] / "no_cavity" / "observables_replica_0.h5"
            effective_sample_ps = read_hdf5_sample_period_ps(h5_path)
    write_calibration_table(
        summaries,
        output,
        block_ps=args.block_ps,
        sample_ps=effective_sample_ps,
    )
    write_analyze_log(
        summaries,
        log_path,
        discard_ps=args.discard_ps,
        prod_end_ps=args.prod_end_ps,
    )
    print(f"Wrote {output}")
    print(f"Wrote {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
