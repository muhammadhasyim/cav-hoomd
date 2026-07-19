#!/usr/bin/env python3
"""Average energy components over completed aging-campaign replicas."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Literal

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prune_timeseries import trim_padded_trajectory  # noqa: E402

DEFAULT_STRIDE = 1
LEGACY_STRIDE = 100
LEGACY_COMPLETE_MIN_BYTES = 1_380_000_000
PRUNED_COMPLETE_MIN_BYTES = 400_000
SWITCH_TIME_PS = 200.0
EARLY_COUPLING_PROBE_PS = 150.0
EARLY_COUPLING_THRESHOLD = 1.0e-4

LAMBDA_TAGS: dict[float, str] = {
    0.0: "lam0p0",
    0.01: "lam0p01",
    0.016667: "lam0p016667",
    0.023333: "lam0p023333",
    0.03: "lam0p03",
}


def find_run_dir(lam_dir: Path) -> Path | None:
    for child in sorted(lam_dir.iterdir()):
        if child.is_dir() and child.name.startswith("cavity_coupling_"):
            return child
    return None


def interp_nans(values: np.ndarray) -> np.ndarray:
    """Fill NaN entries by linear interpolation along the array index."""
    arr = np.asarray(values, dtype=float).copy()
    nans = np.isnan(arr)
    if not np.any(nans):
        return arr
    finites = np.flatnonzero(~nans)
    if finites.size == 0:
        return arr
    arr[nans] = np.interp(np.flatnonzero(nans), finites, arr[finites])
    return arr


def has_early_cavity_coupling(
    time_ps: np.ndarray,
    cavity_coupling: np.ndarray,
    probe_time_ps: float = EARLY_COUPLING_PROBE_PS,
    threshold: float = EARLY_COUPLING_THRESHOLD,
) -> bool:
    """Return True if cavity coupling is already on before the nominal switch."""
    idx = int(np.argmin(np.abs(time_ps - probe_time_ps)))
    return bool(np.mean(cavity_coupling[max(0, idx - 2) : idx + 3]) > threshold)


def align_switch_time(
    time_ps: np.ndarray,
    cavity_coupling: np.ndarray,
    switch_time_ps: float,
    *series: np.ndarray,
) -> tuple[tuple[np.ndarray, ...], bool]:
    """Shift early-on trajectories by ``switch_time_ps``."""
    shifted = has_early_cavity_coupling(time_ps, cavity_coupling)
    if shifted:
        time_ps = time_ps + switch_time_ps
    return (time_ps, *series), shifted


def interpolate_in_range(
    grid: np.ndarray,
    time_ps: np.ndarray,
    values: np.ndarray,
) -> np.ndarray:
    """Interpolate onto ``grid``, writing NaN outside the replica time span."""
    out = np.full(grid.shape, np.nan, dtype=float)
    mask = (grid >= time_ps.min()) & (grid <= time_ps.max())
    if not np.any(mask):
        return out
    out[mask] = np.interp(grid[mask], time_ps, values)
    return out


def average_trimmed_series(
    series_list: list[tuple[np.ndarray, ...]],
) -> dict[str, np.ndarray | int] | None:
    """Average trimmed replica series onto a common time grid."""
    usable = [series for series in series_list if series is not None and len(series[0]) > 1]
    if not usable:
        return None

    # Union time grid: keep the longest trajectory so late times from 2000 ps
    # replicas are not truncated to a shorter 1600 ps partner.
    grid = max((series[0] for series in usable), key=lambda time_ps: float(time_ps[-1])).copy()
    n_components = len(usable[0]) - 1
    layers: list[np.ndarray] = []
    for component_idx in range(n_components):
        stack = []
        for series in usable:
            time_ps = series[0]
            values = series[component_idx + 1]
            stack.append(interpolate_in_range(grid, time_ps, values))
        layers.append(np.vstack(stack))

    means = np.array([np.nanmean(layer, axis=0) for layer in layers])
    counts = np.sum(np.isfinite(layers[0]), axis=0).astype(int)
    stds = np.array(
        [
            np.nanstd(layer, axis=0, ddof=1) if layer.shape[0] > 1 else np.zeros(layer.shape[1])
            for layer in layers
        ]
    )
    # load_replica / align_switch_time order after time_ps:
    # harmonic, lj, coulombic, lj_coul, system_total, universe_total, kinetic_temp
    if means.shape[0] < 7:
        raise ValueError(
            f"expected 7 averaged components after time_ps, got {means.shape[0]}"
        )
    lj = means[1]
    coulombic = means[2]
    return {
        "time_ps": grid,
        "harmonic": means[0],
        "lj": lj,
        "coulombic": coulombic,
        "lj_coul": means[3],
        "system_total": means[4],
        "universe_total": means[5],
        "kinetic_temp": means[6],
        "harmonic_std": stds[0],
        "lj_std": stds[1],
        "coulombic_std": stds[2],
        "lj_coul_std": stds[3],
        "n_replicas": len(usable),
        "coverage": counts,
    }


def _complete_min_bytes(path: Path) -> int:
    with h5py.File(path, "r") as handle:
        period = float(handle.attrs.get("output_period_ps", 0.01))
    if period >= 0.5:
        return PRUNED_COMPLETE_MIN_BYTES
    return LEGACY_COMPLETE_MIN_BYTES


def _effective_stride(path: Path, stride: int) -> int:
    with h5py.File(path, "r") as handle:
        period = float(handle.attrs.get("output_period_ps", 0.01))
    if period >= 0.5:
        return max(1, stride // LEGACY_STRIDE) if stride >= LEGACY_STRIDE else 1
    return stride


def load_replica(
    path: Path,
    stride: int = DEFAULT_STRIDE,
    *,
    structural_energy: Literal["lj_coul", "lj_only"] = "lj_only",
) -> tuple[tuple[np.ndarray, ...], bool] | None:
    """Load, trim, and optionally switch-align one observables file.

    Parameters
    ----------
    path
        Observables HDF5 path for one replica.
    stride
        Output subsampling stride.
    structural_energy
        ``lj_only`` uses the LJ channel alone (default; RF Coulomb omitted).
        ``lj_coul`` averages LJ plus Ewald short/long (or RF in ewald_short).
    """
    try:
        stride = _effective_stride(path, stride)
        with h5py.File(path, "r") as handle:
            sl = slice(None, None, stride)
            time_ps = handle["time"][sl]
            harmonic = handle["energies/harmonic"][sl]
            lj = handle["energies/lj"][sl]
            ewald_short = handle["energies/ewald_short"][sl]
            ewald_long = handle["energies/ewald_long"][sl]
            system_total = handle["energies/system_total"][sl]
            universe_total = handle["energies/universe_total"][sl]
            kinetic_temp = handle["temperatures/kinetic"][sl]
            cavity_coupling = handle["energies/cavity_coupling"][sl]
    except OSError as exc:
        print(f"  skip {path.name}: {exc}", flush=True)
        return None

    lj_coul = lj + ewald_short + ewald_long
    coulombic = ewald_short + ewald_long
    if structural_energy == "lj_only":
        lj_coul = lj.copy()
        coulombic = np.zeros_like(lj)
    trimmed = trim_padded_trajectory(
        time_ps,
        harmonic,
        lj,
        coulombic,
        lj_coul,
        system_total,
        universe_total,
        kinetic_temp,
        cavity_coupling,
    )
    if trimmed is None:
        return None

    (
        time_ps,
        harmonic,
        lj,
        coulombic,
        lj_coul,
        system_total,
        universe_total,
        kinetic_temp,
        cavity_coupling,
    ) = trimmed
    aligned, shifted = align_switch_time(
        time_ps,
        cavity_coupling,
        SWITCH_TIME_PS,
        harmonic,
        lj,
        coulombic,
        lj_coul,
        system_total,
        universe_total,
        kinetic_temp,
    )
    return aligned, shifted


def average_lambda(
    run_dir: Path,
    stride: int = DEFAULT_STRIDE,
    max_replicas: int | None = None,
) -> dict[str, np.ndarray | int] | None:
    replica_paths = sorted(run_dir.glob("observables_replica_*.h5"))
    loaded: list[tuple[np.ndarray, ...]] = []
    n_shifted = 0
    for idx, replica_path in enumerate(replica_paths, start=1):
        if max_replicas is not None and len(loaded) >= max_replicas:
            break
        if replica_path.stat().st_size < _complete_min_bytes(replica_path):
            continue
        result = load_replica(replica_path, stride=stride)
        if result is None:
            continue
        series, shifted = result
        loaded.append(series)
        n_shifted += int(shifted)
        if idx % 50 == 0:
            print(
                f"  ... {idx}/{len(replica_paths)} "
                f"({len(loaded)} loaded, {n_shifted} switch-aligned)",
                flush=True,
            )

    averaged = average_trimmed_series(loaded)
    if averaged is None:
        return None
    averaged["n_shifted"] = n_shifted
    return averaged


def save_csv(output_dir: Path, lam: float, data: dict[str, np.ndarray | int]) -> None:
    columns = [
        data["time_ps"],
        data["harmonic"],
        data["lj"],
        data["coulombic"],
        data["lj_coul"],
        data["system_total"],
        data["universe_total"],
        data["kinetic_temp"],
    ]
    header = (
        "time_ps,harmonic_ha,lj_ha,coulombic_ha,lj_coul_ha,"
        "system_total_ha,universe_total_ha,kinetic_temp_K"
    )
    if "coverage" in data:
        columns.append(data["coverage"])
        header += ",coverage"
    if "coulombic_std" in data:
        columns.extend(
            [
                data["harmonic_std"],
                data["lj_std"],
                data["coulombic_std"],
                data["lj_coul_std"],
            ]
        )
        header += ",harmonic_std_ha,lj_std_ha,coulombic_std_ha,lj_coul_std_ha"
    tag = LAMBDA_TAGS[lam]
    out_path = output_dir / f"avg_energies_{tag}.csv"
    np.savetxt(
        out_path,
        np.column_stack(columns),
        delimiter=",",
        header=header,
        comments="",
    )


def make_plots(output_dir: Path, results: dict[float, dict[str, np.ndarray | int]]) -> None:
    labels = [
        ("harmonic", "Bond (harmonic) Ha"),
        ("lj_coul", "LJ + Coul Ha"),
        ("system_total", "System total Ha"),
    ]
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for ax, (key, ylabel) in zip(axes, labels, strict=True):
        for lam, data in sorted(results.items()):
            ax.plot(
                data["time_ps"],
                data[key],
                label=f"λ={lam:g} (n={data['n_replicas']})",
                alpha=0.85,
            )
        ax.axvline(SWITCH_TIME_PS, color="k", linestyle="--", alpha=0.4)
        ax.set_ylabel(ylabel)
    axes[-1].set_xlabel("Time (ps)")
    axes[0].set_title("Replica-averaged energies — weak-coupling aging campaign")
    axes[0].legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(output_dir / "avg_energy_components_vs_time.png", dpi=150)
    plt.close(fig)

    fig2, axes2 = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    zoom = (150.0, 250.0)
    for ax, (key, ylabel) in zip(axes2, labels, strict=True):
        for lam, data in sorted(results.items()):
            mask = (data["time_ps"] >= zoom[0]) & (data["time_ps"] <= zoom[1])
            ax.plot(data["time_ps"][mask], data[key][mask], alpha=0.85)
        ax.set_ylabel(ylabel)
        ax.set_xlim(zoom)
    axes2[-1].set_xlabel("Time (ps)")
    fig2.suptitle(f"Cavity switch window ({zoom[0]:.0f}–{zoom[1]:.0f} ps)")
    fig2.tight_layout()
    fig2.savefig(output_dir / "avg_energy_switch_zoom.png", dpi=150)
    plt.close(fig2)


def print_summary(results: dict[float, dict[str, np.ndarray | int]]) -> None:
    print("\n=== Switch response: ΔE from t=150→250 ps ===")
    for lam, data in sorted(results.items()):
        post_idx = int(np.argmin(np.abs(data["time_ps"] - 250.0)))
        pre_idx = int(np.argmin(np.abs(data["time_ps"] - 150.0)))
        delta_harmonic = data["harmonic"][post_idx] - data["harmonic"][pre_idx]
        delta_lj_coul = data["lj_coul"][post_idx] - data["lj_coul"][pre_idx]
        delta_system = data["system_total"][post_idx] - data["system_total"][pre_idx]
        print(
            f"λ={lam:8.5f} n={data['n_replicas']:<4d} shifted={data.get('n_shifted', 0):4d} "
            f"Δbond={delta_harmonic:+.5f} Δlj_coul={delta_lj_coul:+.5f} "
            f"Δsystem={delta_system:+.5f}  [OK]"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-base",
        type=Path,
        default=Path("/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda"),
    )
    parser.add_argument("--stride", type=int, default=DEFAULT_STRIDE)
    parser.add_argument("--max-replicas", type=int, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_base = args.output_base
    analysis_dir = output_base / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    results: dict[float, dict[str, np.ndarray | int]] = {}

    for lam, tag in LAMBDA_TAGS.items():
        dir_name = {
            "lam0p0": "lambda0",
            "lam0p01": "lambda0p01",
            "lam0p016667": "lambda0p016667",
            "lam0p023333": "lambda0p023333",
            "lam0p03": "lambda0p03",
        }[tag]
        lam_dir = output_base / dir_name
        run_dir = find_run_dir(lam_dir)
        if run_dir is None:
            print(f"λ={lam}: no run directory", flush=True)
            continue
        print(f"λ={lam}: scanning {run_dir}", flush=True)
        averaged = average_lambda(run_dir, stride=args.stride, max_replicas=args.max_replicas)
        if averaged is None:
            print(f"λ={lam}: no complete replicas", flush=True)
            continue
        print(
            f"λ={lam}: averaged {averaged['n_replicas']} replicas "
            f"({averaged.get('n_shifted', 0)} switch-aligned)",
            flush=True,
        )
        save_csv(analysis_dir, lam, averaged)
        results[lam] = averaged

    if not results:
        raise SystemExit("No completed replicas found.")
    make_plots(analysis_dir, results)
    print_summary(results)
    print(f"Saved outputs under {analysis_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
