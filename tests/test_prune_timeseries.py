"""Tests for aging_weak_lambda time-series pruning utilities."""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "aging_weak_lambda"),
)

from prune_timeseries import (
    downsample_to_period,
    is_uniform_grid,
    prune_fkt_replica_file,
    prune_hdf5_observables,
    read_fkt_series,
    trim_padded_trajectory,
    verify_pruned_hdf5,
)


def _write_synthetic_hdf5(path: Path, n_valid: int = 1001, pad: int = 50) -> None:
    """Create a minimal observables HDF5 file at 0.01 ps sampling."""
    total = n_valid + pad
    time = np.zeros(total, dtype=np.float64)
    time[:2] = 0.0  # leading zeros before first write
    time[2:n_valid] = 0.001 + np.arange(n_valid - 2) * 0.01
    harmonic = np.zeros(total, dtype=np.float64)
    harmonic[2:n_valid] = -1.0 + 1.0e-4 * np.arange(n_valid - 2)
    lj = np.zeros(total, dtype=np.float64)
    lj[2:n_valid] = -1.4 + 2.0e-4 * np.arange(n_valid - 2)

    with h5py.File(path, "w", libver="latest") as handle:
        handle.attrs["format_version"] = "1.0"
        handle.attrs["output_period_ps"] = 0.01
        handle.create_group("energies")
        handle.create_group("temperatures")
        handle.create_dataset("time", data=time)
        handle["energies"].create_dataset("harmonic", data=harmonic)
        handle["energies"].create_dataset("lj", data=lj)
        handle["temperatures"].create_dataset("kinetic", data=100.0 + 0.01 * harmonic)


def test_trim_padded_trajectory_removes_leading_and_trailing() -> None:
    time = np.array([0.0, 0.0, 0.001, 0.011, 0.021, 0.0, 0.0])
    harmonic = np.array([0.0, 0.0, -1.0, -0.9, -0.8, 0.0, 0.0])
    lj = np.array([0.0, 0.0, -1.4, -1.3, -1.2, 0.0, 0.0])
    trimmed = trim_padded_trajectory(time, harmonic, lj)
    assert trimmed is not None
    t, h, l = trimmed
    np.testing.assert_allclose(t, [0.001, 0.011, 0.021])
    np.testing.assert_allclose(h, [-1.0, -0.9, -0.8])
    np.testing.assert_allclose(l, [-1.4, -1.3, -1.2])


def test_trim_padded_trajectory_returns_none_for_empty() -> None:
    assert trim_padded_trajectory(np.zeros(5), np.zeros(5)) is None


def test_downsample_to_period_selects_nearest_one_ps_grid() -> None:
    time = 0.001 + np.arange(100) * 0.01
    harmonic = time.copy()
    t_ds, h_ds = downsample_to_period(time, harmonic, target_period_ps=1.0, origin_ps=0.0)
    assert len(t_ds) == 2  # 0 ps and 1 ps within 0.001..0.991
    assert np.all(np.diff(t_ds) >= 0.99)
    assert h_ds.shape == t_ds.shape


def test_prune_hdf5_observables(tmp_path: Path) -> None:
    src = tmp_path / "observables_replica_0.h5"
    dst = tmp_path / "observables_replica_0_pruned.h5"
    _write_synthetic_hdf5(src, n_valid=1001, pad=50)
    stats = prune_hdf5_observables(src, dst, target_period_ps=1.0)
    assert stats["n_out"] < stats["n_in"]
    assert stats["target_period_ps"] == 1.0

    with h5py.File(dst, "r") as handle:
        assert handle.attrs["output_period_ps"] == 1.0
        time = handle["time"][:]
        assert time[-1] <= 10.0 + 1e-6
        assert len(time) == stats["n_out"]
        assert handle["energies/harmonic"].shape == time.shape
        assert handle["temperatures/kinetic"].shape == time.shape

    verify_pruned_hdf5(src, dst, target_period_ps=1.0)


def test_verify_pruned_hdf5_rejects_bad_grid(tmp_path: Path) -> None:
    src = tmp_path / "src.h5"
    bad = tmp_path / "bad.h5"
    _write_synthetic_hdf5(src)
    prune_hdf5_observables(src, bad, target_period_ps=1.0)
    with h5py.File(bad, "a") as handle:
        handle["time"][1] = handle["time"][1] + 0.5
    with pytest.raises(ValueError, match="uniform"):
        verify_pruned_hdf5(src, bad, target_period_ps=1.0)


def test_read_fkt_series_parses_interleaved_file(tmp_path: Path) -> None:
    path = tmp_path / "prod-0_fkt_ref_000.txt"
    path.write_text(
        "# F(k,t) correlation function\n"
        "# Reference time: 0.001 ps\n"
        "# lag_time_ps\tF(k,t)\n"
        "0.000000\t920.63913159\n"
        "445.264000\t-87.73328097\n"
        "1.000000\t587.44373806\n"
        "2.001000\t599.71878341\n",
        encoding="utf-8",
    )
    header, times, values = read_fkt_series(path)
    assert "Reference time" in header[1]
    assert len(times) == 4
    assert times[0] == pytest.approx(0.0)


def test_prune_fkt_replica_file_to_uniform_grid(tmp_path: Path) -> None:
    src = tmp_path / "prod-0_fkt_ref_000.txt"
    dst = tmp_path / "prod-0_fkt_ref_000_pruned.txt"
    lines = ["# F(k,t) correlation function\n", "# Reference time: 0.001 ps\n", "# lag_time_ps\tF(k,t)\n"]
    for lag in range(0, 11):
        lines.append(f"{float(lag):.6f}\t{100.0 - lag:.8f}\n")
    src.write_text("".join(lines), encoding="utf-8")
    stats = prune_fkt_replica_file(src, dst, target_period_ps=1.0)
    assert stats["n_out"] == 11
    header, times, values = read_fkt_series(dst)
    assert is_uniform_grid(times, period_ps=1.0)
    assert len(header) >= 3
    assert values[0] == pytest.approx(100.0)


def test_is_uniform_grid() -> None:
    assert is_uniform_grid(np.array([0.0, 1.0, 2.0, 3.0]), period_ps=1.0)
    assert not is_uniform_grid(np.array([0.0, 1.5, 2.0]), period_ps=1.0)
