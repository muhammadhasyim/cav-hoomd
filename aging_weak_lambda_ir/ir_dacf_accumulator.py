#!/usr/bin/env python3
"""Incremental dipole-autocorrelation accumulation for continuous IR production."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from analyze_ir_from_dipole import (
    CAMPAIGN_LAMBDAS,
    DEFAULT_TEMPERATURE_K,
    dipole_acf_xy_fft,
    ir_spectrum_dct_from_acf,
    ir_spectrum_mesa_fpe_from_dipole_xy,
    trim_monotonic_trajectory,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def average_acfs(acfs: list[np.ndarray]) -> tuple[np.ndarray, int]:
    """Average ACF arrays, truncating each to the shortest length."""
    if not acfs:
        return np.array([]), 0
    min_len = min(acf.size for acf in acfs)
    stacked = np.stack([acf[:min_len] for acf in acfs], axis=0)
    return np.mean(stacked, axis=0), stacked.shape[0]


def discover_run_observables(run_dir: Path) -> list[Path]:
    """Return observables HDF5 files under a single run directory."""
    return sorted(run_dir.glob("**/observables_replica_*.h5"))


def load_partial_dipole_hdf5(
    h5_path: Path,
    *,
    max_frames: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Load dipole trajectory, optionally limiting to the first ``max_frames``."""
    last_error: Exception | None = None
    for swmr in (True, False):
        try:
            with h5py.File(h5_path, "r", swmr=swmr) as handle:
                times = np.asarray(handle["time"][:], dtype=float)
                dipole = np.asarray(
                    handle["order_parameters/dipole/components"][:],
                    dtype=float,
                )
            break
        except OSError as exc:
            last_error = exc
            continue
    else:
        assert last_error is not None
        raise last_error

    if max_frames is not None and max_frames > 0:
        times = times[:max_frames]
        dipole = dipole[:max_frames]
    times, dipole = trim_monotonic_trajectory(times, dipole)
    return times, dipole


def acf_from_hdf5(h5_path: Path) -> tuple[np.ndarray, np.ndarray, int]:
    """Compute DACF and return (acf, times_ps, n_frames)."""
    times, dipole = load_partial_dipole_hdf5(h5_path)
    if dipole.shape[0] < 4:
        return np.array([]), times, dipole.shape[0]
    return dipole_acf_xy_fft(dipole), times, dipole.shape[0]


def spectrum_from_mean_acf(
    mean_acf: np.ndarray,
    dt_fs: float,
) -> tuple[np.ndarray, np.ndarray]:
    """DCT IR spectrum from ensemble-mean in-plane ACF."""
    return ir_spectrum_dct_from_acf(mean_acf, dt_fs)


def ensemble_mesa_spectrum_from_state(
    state: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    """Average MESA/FPE spectra from mu_x and mu_y for each stored run."""
    specs: list[np.ndarray] = []
    freqs_ref: np.ndarray | None = None
    for run in state.get("runs", {}).values():
        if not run.get("completed", False):
            continue
        h5_path = Path(run["h5_path"])
        if not h5_path.is_file():
            continue
        try:
            times, dipole = load_partial_dipole_hdf5(h5_path)
        except OSError:
            continue
        freqs, spec = ir_spectrum_mesa_fpe_from_dipole_xy(dipole, times)
        if freqs.size == 0:
            continue
        if freqs_ref is None:
            freqs_ref = freqs
        elif freqs_ref.shape != freqs.shape:
            continue
        specs.append(spec)
    if freqs_ref is None or not specs:
        return np.array([]), np.array([])
    return freqs_ref, np.mean(np.stack(specs, axis=0), axis=0)


def _empty_state(lam: float, tag: str) -> dict[str, Any]:
    return {
        "lam": lam,
        "tag": tag,
        "updated_at": utc_now(),
        "run_count": 0,
        "completed_runs": 0,
        "mean_acf_len": 0,
        "dt_fs": None,
        "runs": {},
    }


def _load_state(path: Path, lam: float, tag: str) -> dict[str, Any]:
    if not path.is_file():
        return _empty_state(lam, tag)
    state = json.loads(path.read_text(encoding="utf-8"))
    state.setdefault("runs", {})
    return state


def _save_state(path: Path, state: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    state["updated_at"] = utc_now()
    path.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")


def update_lambda_state(
    state_path: Path,
    *,
    run_id: str,
    h5_path: Path,
    lam: float,
    tag: str,
    replica: int,
    completed: bool,
) -> dict[str, Any]:
    """Update per-lambda accumulator state from one run's current HDF5 snapshot."""
    state = _load_state(state_path, lam, tag)
    acf, times, n_frames = acf_from_hdf5(h5_path)
    if acf.size == 0:
        return state

    times_fs = np.asarray(times, dtype=float) * 1000.0
    dt_fs = float(np.mean(np.diff(times_fs))) if times_fs.size > 1 else None
    if dt_fs is not None and dt_fs > 0.0:
        state["dt_fs"] = dt_fs

    state["runs"][run_id] = {
        "h5_path": str(h5_path),
        "replica": replica,
        "n_frames": n_frames,
        "duration_ps": float(times[-1]) if times.size else 0.0,
        "completed": completed,
        "updated_at": utc_now(),
        "acf": acf.tolist(),
    }
    state["run_count"] = len(state["runs"])
    state["completed_runs"] = sum(
        1 for run in state["runs"].values() if run.get("completed")
    )

    acfs = [np.asarray(run["acf"], dtype=float) for run in state["runs"].values()]
    mean_acf, _ = average_acfs(acfs)
    state["mean_acf_len"] = int(mean_acf.size)
    state["mean_acf"] = mean_acf.tolist()
    _save_state(state_path, state)
    return state


def write_dashboard(
    summary: dict[str, Any],
    text_path: Path,
    json_path: Path,
) -> None:
    """Write human-readable and JSON dashboard snapshots."""
    text_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    lines = [
        f"IR continuous dashboard updated {summary.get('updated_at', '')}",
        "",
        "Active runs:",
    ]
    active = summary.get("active_runs", [])
    if not active:
        lines.append("  (none)")
    else:
        for run in active:
            lines.append(
                "  gpu={gpu} lam={lam:g} tag={tag} replica={replica} "
                "run_id={run_id} frames={n_frames} duration_ps={duration_ps:.2f}".format(
                    gpu=run.get("gpu", "?"),
                    lam=run.get("lam", float("nan")),
                    tag=run.get("tag", "?"),
                    replica=run.get("replica", "?"),
                    run_id=run.get("run_id", "?"),
                    n_frames=run.get("n_frames", 0),
                    duration_ps=run.get("duration_ps", 0.0),
                )
            )

    lines.extend(["", "Lambda ensemble:"])
    for tag, info in sorted(summary.get("lambdas", {}).items()):
        lines.append(
            f"  lambda{tag}: runs={info.get('run_count', 0)} "
            f"completed={info.get('completed_runs', 0)} "
            f"frames_latest={info.get('n_frames_latest', 0)} "
            f"peak_cm1={info.get('peak_cm1', float('nan')):.1f}"
        )
    text_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@dataclass(frozen=True)
class ActiveRun:
    run_id: str
    run_dir: Path
    lam: float
    tag: str
    replica: int
    gpu: int
    started_at: str
    completed: bool = False


def load_active_runs(registry_path: Path) -> list[ActiveRun]:
    if not registry_path.is_file():
        return []
    runs: list[ActiveRun] = []
    for line in registry_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        row = json.loads(line)
        runs.append(
            ActiveRun(
                run_id=row["run_id"],
                run_dir=Path(row["run_dir"]),
                lam=float(row["lam"]),
                tag=row["tag"],
                replica=int(row["replica"]),
                gpu=int(row["gpu"]),
                started_at=row["started_at"],
                completed=bool(row.get("completed", False)),
            )
        )
    return runs


def append_active_run(registry_path: Path, run: ActiveRun) -> None:
    registry_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "run_id": run.run_id,
        "run_dir": str(run.run_dir),
        "lam": run.lam,
        "tag": run.tag,
        "replica": run.replica,
        "gpu": run.gpu,
        "started_at": run.started_at,
        "completed": run.completed,
    }
    with registry_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload) + "\n")


def mark_run_completed(registry_path: Path, run_id: str) -> None:
    runs = load_active_runs(registry_path)
    registry_path.write_text("", encoding="utf-8")
    for run in runs:
        completed = run.completed or run.run_id == run_id
        append_active_run(
            registry_path,
            ActiveRun(
                run_id=run.run_id,
                run_dir=run.run_dir,
                lam=run.lam,
                tag=run.tag,
                replica=run.replica,
                gpu=run.gpu,
                started_at=run.started_at,
                completed=completed,
            ),
        )


def refresh_accumulator(
    *,
    base_dir: Path,
    registry_path: Path,
    derived_dir: Path,
    method: str = "mesa",
    temperature_k: float = DEFAULT_TEMPERATURE_K,
) -> dict[str, Any]:
    """Poll active runs, update ACF states, and regenerate spectra/dashboard."""
    state_dir = derived_dir / "accumulator"
    state_dir.mkdir(parents=True, exist_ok=True)

    runs = load_active_runs(registry_path)
    active_runs: list[dict[str, Any]] = []
    lambda_summary: dict[str, dict[str, Any]] = {}

    for run in runs:
        h5_paths = discover_run_observables(run.run_dir)
        if not h5_paths:
            continue
        h5_path = h5_paths[0]
        try:
            state = update_lambda_state(
                state_dir / f"lambda{run.tag}_state.json",
                run_id=run.run_id,
                h5_path=h5_path,
                lam=run.lam,
                tag=run.tag,
                replica=run.replica,
                completed=run.completed,
            )
        except OSError:
            continue
        if not run.completed:
            active_runs.append(
                {
                    "run_id": run.run_id,
                    "gpu": run.gpu,
                    "lam": run.lam,
                    "tag": run.tag,
                    "replica": run.replica,
                    "n_frames": state["runs"][run.run_id]["n_frames"],
                    "duration_ps": state["runs"][run.run_id]["duration_ps"],
                }
            )

    for lam, tag in CAMPAIGN_LAMBDAS:
        state_path = state_dir / f"lambda{tag}_state.json"
        if not state_path.is_file():
            continue
        state = json.loads(state_path.read_text(encoding="utf-8"))
        mean_acf = np.asarray(state.get("mean_acf", []), dtype=float)
        dt_fs = state.get("dt_fs")
        peak_cm1 = float("nan")
        if mean_acf.size >= 4 and dt_fs is not None:
            if method == "mesa":
                freqs, spec = ensemble_mesa_spectrum_from_state(state)
            else:
                freqs, spec = spectrum_from_mean_acf(mean_acf, float(dt_fs))
            if freqs.size > 0:
                mask = (freqs >= 1200.0) & (freqs <= 2200.0)
                if np.any(mask):
                    peak_cm1 = float(freqs[mask][int(np.argmax(spec[mask]))])
                out_txt = derived_dir / f"ir_spectrum_lambda{tag}.txt"
                np.savetxt(
                    out_txt,
                    np.column_stack([freqs, spec]),
                    header="frequency_cm-1 spectrum_au",
                    comments="# ",
                )

        n_frames_latest = 0
        if state.get("runs"):
            n_frames_latest = max(
                int(run["n_frames"]) for run in state["runs"].values()
            )
        lambda_summary[tag] = {
            "lam": lam,
            "run_count": state.get("run_count", 0),
            "completed_runs": state.get("completed_runs", 0),
            "n_frames_latest": n_frames_latest,
            "mean_acf_len": state.get("mean_acf_len", 0),
            "peak_cm1": peak_cm1,
        }

    summary = {
        "updated_at": utc_now(),
        "base_dir": str(base_dir),
        "active_runs": active_runs,
        "lambdas": lambda_summary,
    }
    write_dashboard(
        summary,
        derived_dir / "dashboard.txt",
        derived_dir / "status.json",
    )
    return summary
