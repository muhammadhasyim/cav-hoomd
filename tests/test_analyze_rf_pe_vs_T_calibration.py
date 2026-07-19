"""Tests for RF potential-energy vs temperature calibration analysis."""

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
    block_average_series,
    extract_temperature_k_from_dirname,
    load_production_window,
    summarize_temperature_run,
    write_calibration_table,
)


def test_block_average_series_known_mean_and_std() -> None:
    """Block means of a constant series have zero spread."""
    values = np.full(200, 2.5, dtype=float)
    mean, std, n_blocks = block_average_series(values, block_size=100)
    assert mean == pytest.approx(2.5)
    assert std == pytest.approx(0.0)
    assert n_blocks == 2


def test_block_average_series_linear_trend() -> None:
    """Two blocks [0..99] and [100..199] give block means 49.5 and 149.5."""
    values = np.arange(200, dtype=float)
    mean, std, n_blocks = block_average_series(values, block_size=100)
    assert mean == pytest.approx(99.5)
    assert std == pytest.approx(70.71067811865476)
    assert n_blocks == 2


def test_block_average_series_rejects_incomplete_tail() -> None:
    with pytest.raises(ValueError, match="multiple of block_size"):
        block_average_series(np.arange(105, dtype=float), block_size=100)


def test_extract_temperature_k_from_dirname() -> None:
    assert extract_temperature_k_from_dirname("temp_300K") == 300.0
    assert extract_temperature_k_from_dirname("temp_80K") == 80.0


def _write_calibration_hdf5(
    path: Path,
    *,
    n_prod_points: int = 200,
    sample_ps: float = 0.01,
    discard_ps: float = 10.0,
    harmonic_prod: float = 0.05,
    lj_prod: float = -0.70,
    rf_prod: float = 1.25,
) -> None:
    """Minimal HDF5 mimicking 2 ns calib output with short discard for tests."""
    n_discard = int(round(discard_ps / sample_ps))
    time = np.concatenate(
        [
            np.arange(1, n_discard + 1) * sample_ps,
            discard_ps + np.arange(1, n_prod_points + 1) * sample_ps,
        ]
    )
    n_total = len(time)
    harmonic = np.zeros(n_total, dtype=float)
    lj = np.zeros(n_total, dtype=float)
    rf = np.zeros(n_total, dtype=float)
    temp = np.zeros(n_total, dtype=float)
    harmonic[:n_discard] = 0.20
    lj[:n_discard] = -0.30
    rf[:n_discard] = 0.50
    temp[:n_discard] = 250.0
    harmonic[n_discard:] = harmonic_prod
    lj[n_discard:] = lj_prod
    rf[n_discard:] = rf_prod
    temp[n_discard:] = 100.0

    with h5py.File(path, "w", libver="latest") as handle:
        handle.attrs["output_period_ps"] = sample_ps
        handle.create_dataset("time", data=time)
        energies = handle.create_group("energies")
        energies.create_dataset("harmonic", data=harmonic)
        energies.create_dataset("lj", data=lj)
        energies.create_dataset("ewald_short", data=rf)
        energies.create_dataset("ewald_long", data=np.zeros(n_total))
        energies.create_dataset("temperature", data=temp)


def test_load_production_window_applies_discard(tmp_path: Path) -> None:
    h5_path = tmp_path / "observables_replica_0.h5"
    _write_calibration_hdf5(h5_path, discard_ps=10.0, n_prod_points=200)
    window = load_production_window(h5_path, discard_ps=10.0, prod_end_ps=12.0)
    assert window.n_samples == 200
    assert window.time_ps[0] == pytest.approx(10.01, abs=1e-9)
    assert window.harmonic_ha[0] == pytest.approx(0.05)
    assert window.lj_ha[0] == pytest.approx(-0.70)
    assert window.rf_coul_ha[0] == pytest.approx(1.25)


def test_summarize_temperature_run(tmp_path: Path) -> None:
    run_dir = tmp_path / "temp_100K" / "no_cavity"
    run_dir.mkdir(parents=True)
    h5_path = run_dir / "observables_replica_0.h5"
    _write_calibration_hdf5(
        h5_path,
        discard_ps=10.0,
        n_prod_points=200,
        harmonic_prod=0.04,
        lj_prod=-0.69,
        rf_prod=1.28,
    )
    summary = summarize_temperature_run(
        run_dir.parent,
        discard_ps=10.0,
        prod_end_ps=12.0,
        block_ps=1.0,
        sample_ps=0.01,
    )
    assert summary.temperature_k == 100.0
    assert summary.harmonic_mean_ha == pytest.approx(0.04)
    assert summary.lj_mean_ha == pytest.approx(-0.69)
    assert summary.coulombic_mean_ha == pytest.approx(1.28)
    assert summary.t_kin_mean_k == pytest.approx(100.0)
    assert summary.n_blocks == 2


def test_write_calibration_table(tmp_path: Path) -> None:
    run_root = tmp_path / "pe_vs_T_calib_rf"
    for temp_k in (300.0, 100.0):
        run_dir = run_root / f"temp_{int(temp_k)}K" / "no_cavity"
        run_dir.mkdir(parents=True)
        _write_calibration_hdf5(run_dir / "observables_replica_0.h5", discard_ps=10.0)

    out_path = run_root / "potential_energy_vs_T_rf.txt"
    summaries = [
        summarize_temperature_run(
            run_root / f"temp_{int(temp_k)}K",
            discard_ps=10.0,
            prod_end_ps=12.0,
            block_ps=1.0,
            sample_ps=0.01,
        )
        for temp_k in (300.0, 100.0)
    ]
    write_calibration_table(summaries, out_path)
    text = out_path.read_text(encoding="utf-8")
    assert "temperature\tharmonic_hartree" in text
    assert "300.00" in text
    assert "100.00" in text
