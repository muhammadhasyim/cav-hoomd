"""Tests for Figure 3 display-only smoothing (science CSV unchanged)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPRO_ROOT = Path(__file__).resolve().parents[1] / "reproduction"
sys.path.insert(0, str(REPRO_ROOT / "figure3_energy"))

from plot_figure3 import (  # noqa: E402
    apply_display_smoothing,
    load_figure3_data,
    moving_average,
)


def _write_staged_pe(path: Path, *, n: int = 50) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    time_ps = np.arange(1, n + 1, dtype=float)
    harmonic = np.full(n, 0.04)
    lj = np.full(n, -0.69)
    coul = 1.2 + rng.normal(0, 0.4, size=n)
    total = harmonic + lj + coul
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# test staged PE\n")
        fh.write("time_ps harmonic_ha lj_ha coul_ha total_ha\n")
        for idx in range(n):
            fh.write(
                f"{time_ps[idx]:.6e} {harmonic[idx]:.6e} "
                f"{lj[idx]:.6e} {coul[idx]:.6e} {total[idx]:.6e}\n"
            )


def _write_staged_companion_files(base_dir: Path, tag: str, *, n: int = 50) -> None:
    time_ps = np.arange(1, n + 1, dtype=float)
    fict = np.column_stack([time_ps, np.full(n, 100.0), np.full(n, 105.0), np.full(n, 110.0)])
    fts = np.column_stack([time_ps, np.full(n, 102.0), np.full(n, 0.5), np.full(n, 0.5)])
    np.savetxt(
        base_dir / f"coupling_{tag}_averaged_fictive_temperatures.txt",
        fict,
        fmt="%.6e",
    )
    np.savetxt(
        base_dir / f"coupling_{tag}_fictive_time_series.txt",
        fts,
        fmt="%.6e",
    )


def test_moving_average_reduces_high_frequency_noise() -> None:
    x = np.array([0.0, 1.0, -1.0, 1.0, -1.0, 0.0], dtype=float)
    smoothed = moving_average(x, window=3)
    assert smoothed.shape == x.shape
    assert float(np.std(smoothed)) < float(np.std(x))


def test_apply_display_smoothing_leaves_harmonic_unsmoothed_by_default() -> None:
    data = {
        "dV_harm": np.linspace(0.0, 1.0, 20),
        "dV_ljc": np.random.default_rng(0).normal(size=20),
        "dV_tot": np.random.default_rng(1).normal(size=20),
        "dV_lj": None,
        "dV_coul": None,
    }
    display = apply_display_smoothing(
        data,
        smooth_window=5,
        smooth_window_total=5,
        smooth_harmonic=False,
    )
    assert np.allclose(display["dV_harm"], data["dV_harm"])
    assert not np.allclose(display["dV_ljc"], data["dV_ljc"])


def test_load_figure3_data_unchanged_by_display_smoothing(tmp_path: Path) -> None:
    """Staged science series is read once; display smoothing is a separate step."""
    tag = "0epos00"
    base = tmp_path / "time_series_output"
    base.mkdir()
    _write_staged_pe(base / f"coupling_{tag}_averaged_potential_energy.txt")
    _write_staged_companion_files(base, tag)

    class _Entry:
        epsilon_tag = tag

    class _Profile:
        def entry_for_axis_value(self, _coupling: float) -> _Entry:
            return _Entry()

    data = load_figure3_data(base, 0.0, profile=_Profile())
    raw_ljc = data["dV_ljc"].copy()
    display = apply_display_smoothing(data, smooth_window=7, smooth_window_total=7)
    assert np.allclose(data["dV_ljc"], raw_ljc)
    assert float(np.std(display["dV_ljc"])) <= float(np.std(raw_ljc)) + 1e-12

    staged_text = (base / f"coupling_{tag}_averaged_potential_energy.txt").read_text()
    data_reread = load_figure3_data(base, 0.0, profile=_Profile())
    assert np.allclose(data_reread["dV_ljc"], raw_ljc)
    assert "time_ps harmonic_ha" in staged_text
