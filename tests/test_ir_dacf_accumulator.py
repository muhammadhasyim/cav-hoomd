"""Tests for incremental IR DACF accumulation."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_IR_DIR = _REPO_ROOT / "aging_weak_lambda_ir"
if str(_IR_DIR) not in sys.path:
    sys.path.insert(0, str(_IR_DIR))

from ir_dacf_accumulator import (  # noqa: E402
    average_acfs,
    discover_run_observables,
    load_partial_dipole_hdf5,
    update_lambda_state,
    write_dashboard,
)


def _write_observables_h5(
    path: Path,
    n_frames: int,
    dt_ps: float,
    freq_cm1: float = 1560.0,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    times = np.arange(n_frames, dtype=float) * dt_ps
    omega = 2.0 * np.pi * freq_cm1 * 2.99792458e10 * 1e-12
    dipole = np.zeros((n_frames, 3), dtype=float)
    dipole[:, 0] = np.cos(omega * times)
    with h5py.File(path, "w") as handle:
        handle.create_dataset("time", data=times)
        handle.create_dataset(
            "order_parameters/dipole/components",
            data=dipole,
        )


def test_average_acfs_equal_length() -> None:
    acf_a = np.array([1.0, 0.5, 0.25], dtype=float)
    acf_b = np.array([3.0, 1.5, 0.75], dtype=float)
    mean, count = average_acfs([acf_a, acf_b])
    assert count == 2
    np.testing.assert_allclose(mean, [2.0, 1.0, 0.5])


def test_average_acfs_truncates_to_shortest() -> None:
    acf_a = np.array([1.0, 0.5, 0.25, 0.1], dtype=float)
    acf_b = np.array([3.0, 1.5], dtype=float)
    mean, count = average_acfs([acf_a, acf_b])
    assert count == 2
    assert mean.shape == (2,)
    np.testing.assert_allclose(mean, [2.0, 1.0])


def test_update_lambda_state_tracks_multiple_runs(tmp_path: Path) -> None:
    state_path = tmp_path / "lambda0p01_state.json"
    h5_a = tmp_path / "run_a" / "observables_replica_0.h5"
    h5_b = tmp_path / "run_b" / "observables_replica_0.h5"
    _write_observables_h5(h5_a, n_frames=200, dt_ps=0.002)
    _write_observables_h5(h5_b, n_frames=200, dt_ps=0.002)

    state = update_lambda_state(
        state_path,
        run_id="run_a",
        h5_path=h5_a,
        lam=0.01,
        tag="0p01",
        replica=0,
        completed=False,
    )
    assert state["run_count"] == 1
    assert state["runs"]["run_a"]["n_frames"] == 200

    state = update_lambda_state(
        state_path,
        run_id="run_b",
        h5_path=h5_b,
        lam=0.01,
        tag="0p01",
        replica=1,
        completed=True,
    )
    assert state["run_count"] == 2
    assert state["completed_runs"] == 1
    assert state["mean_acf_len"] == 200


def test_discover_run_observables(tmp_path: Path) -> None:
    h5 = tmp_path / "coupling_100000e-07" / "observables_replica_0.h5"
    _write_observables_h5(h5, n_frames=10, dt_ps=0.002)
    found = discover_run_observables(tmp_path)
    assert found == [h5]


def test_load_partial_dipole_hdf5_respects_max_frames(tmp_path: Path) -> None:
    h5 = tmp_path / "observables_replica_0.h5"
    _write_observables_h5(h5, n_frames=50, dt_ps=0.002)
    times, dipole = load_partial_dipole_hdf5(h5, max_frames=20)
    assert times.shape == (20,)
    assert dipole.shape == (20, 3)


def test_write_dashboard_creates_summary(tmp_path: Path) -> None:
    summary = {
        "updated_at": "2026-07-19T10:00:00Z",
        "active_runs": [{"run_id": "r1", "gpu": 0, "lam": 0.01, "tag": "0p01"}],
        "lambdas": {
            "0p01": {
                "run_count": 2,
                "completed_runs": 1,
                "n_frames_latest": 1000,
                "peak_cm1": 1560.0,
            }
        },
    }
    out_txt = tmp_path / "dashboard.txt"
    out_json = tmp_path / "status.json"
    write_dashboard(summary, out_txt, out_json)
    assert out_txt.is_file()
    assert out_json.is_file()
    loaded = json.loads(out_json.read_text(encoding="utf-8"))
    assert loaded["lambdas"]["0p01"]["run_count"] == 2
