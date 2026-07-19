"""Tests for aging campaign resume I/O helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from examples.slurm.aging_resume_io import (
    infer_hdf5_write_index,
    load_fkt_state,
    max_lag_for_reference,
    normalize_fkt_series,
    parse_fkt_reference_index,
    parse_fkt_reference_time,
    save_fkt_state,
)


def test_parse_fkt_reference_time(tmp_path: Path) -> None:
    path = tmp_path / "prod-0_fkt_ref_000.txt"
    path.write_text(
        "# F(k,t) correlation function\n"
        "# Reference time: 1404.000 ps\n"
        "# lag_time_ps\tF(k,t)\n"
        "0.000000\t1.000000\n",
        encoding="utf-8",
    )
    assert parse_fkt_reference_time(path) == 1404.0


def test_parse_fkt_reference_index() -> None:
    assert parse_fkt_reference_index(Path("prod-12_fkt_ref_007.txt")) == 7


def test_max_lag_for_reference() -> None:
    assert max_lag_for_reference(reference_time_ps=1400.0, runtime_ps=2000.0) == 600.0


def test_normalize_fkt_series_clips_and_normalizes() -> None:
    lag_times = np.array([0.0, 100.0, 700.0])
    values = np.array([2.0, 1.0, 0.5])
    clipped_lags, normalized = normalize_fkt_series(
        lag_times,
        values,
        reference_time_ps=1400.0,
        runtime_ps=2000.0,
    )
    assert np.allclose(clipped_lags, [0.0, 100.0])
    assert np.allclose(normalized, [1.0, 0.5])


def test_fkt_state_round_trip(tmp_path: Path) -> None:
    references = [
        {
            "timestep": 10,
            "time_ps": 1.0,
            "rhok_real": np.array([1.0, 2.0]),
            "rhok_imag": np.array([0.5, -0.5]),
        }
    ]
    path = tmp_path / "prod-0_fkt_state.npz"
    save_fkt_state(
        path,
        references=references,
        last_reference_time=1.0,
        last_output_time=0.0,
        kmag=1.0,
        reference_interval_ps=200.0,
    )
    loaded = load_fkt_state(path)
    assert len(loaded["references"]) == 1
    assert loaded["references"][0]["time_ps"] == 1.0
    assert loaded["last_reference_time"] == 1.0


def test_infer_hdf5_write_index_ignores_zero_padding() -> None:
    times = np.array([0.0, 1.0, 2.0, 0.0, 0.0])
    assert infer_hdf5_write_index(times) == 3
