"""Tests for RF calib sample_ps from HDF5 output_period_ps."""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "aging_weak_lambda" / "analysis"),
)

from analyze_rf_pe_vs_T_calibration import (  # noqa: E402
    read_hdf5_sample_period_ps,
    resolve_sample_period_ps,
    summarize_temperature_run,
)


def _write_noisy_rf_hdf5(
    path: Path,
    *,
    output_period_ps: float,
    rf_values: np.ndarray,
    discard_ps: float = 10.0,
) -> None:
    n_discard = int(round(discard_ps / output_period_ps))
    prod = np.asarray(rf_values, dtype=float)
    n_prod = prod.size
    time = np.concatenate(
        [
            np.arange(1, n_discard + 1) * output_period_ps,
            discard_ps + np.arange(1, n_prod + 1) * output_period_ps,
        ]
    )
    n_total = len(time)
    harmonic = np.full(n_total, 0.04)
    lj = np.full(n_total, -0.69)
    rf = np.zeros(n_total)
    rf[n_discard:] = prod
    temp = np.full(n_total, 100.0)

    with h5py.File(path, "w", libver="latest") as handle:
        handle.attrs["output_period_ps"] = output_period_ps
        handle.create_dataset("time", data=time)
        energies = handle.create_group("energies")
        energies.create_dataset("harmonic", data=harmonic)
        energies.create_dataset("lj", data=lj)
        energies.create_dataset("ewald_short", data=rf)
        energies.create_dataset("temperature", data=temp)


def test_read_hdf5_sample_period_ps_from_attrs(tmp_path: Path) -> None:
    path = tmp_path / "obs.h5"
    with h5py.File(path, "w") as handle:
        handle.attrs["output_period_ps"] = 1.0
        handle.create_dataset("time", data=[1.0, 2.0])
    assert read_hdf5_sample_period_ps(path) == pytest.approx(1.0)


def test_resolve_sample_period_ps_prefers_hdf5_over_default() -> None:
    assert resolve_sample_period_ps(0.01, hdf5_sample_ps=1.0) == pytest.approx(1.0)
    assert resolve_sample_period_ps(None, hdf5_sample_ps=1.0) == pytest.approx(1.0)
    assert resolve_sample_period_ps(0.5, hdf5_sample_ps=None) == pytest.approx(0.5)


def test_summarize_uses_hdf5_output_period_for_block_size(tmp_path: Path) -> None:
    """With output_period_ps=1.0 and block_ps=1.0, block_size=1 and std matches raw samples."""
    run_dir = tmp_path / "temp_100K" / "no_cavity"
    run_dir.mkdir(parents=True)
    rng = np.random.default_rng(0)
    rf_prod = rng.normal(loc=1.25, scale=0.40, size=100)
    _write_noisy_rf_hdf5(
        run_dir / "observables_replica_0.h5",
        output_period_ps=1.0,
        rf_values=rf_prod,
    )
    summary = summarize_temperature_run(
        run_dir.parent,
        discard_ps=10.0,
        prod_end_ps=110.0,
        block_ps=1.0,
        sample_ps=None,
    )
    assert summary.n_blocks == 100
    assert summary.coulombic_std_ha == pytest.approx(float(np.std(rf_prod, ddof=1)), rel=1e-6)
