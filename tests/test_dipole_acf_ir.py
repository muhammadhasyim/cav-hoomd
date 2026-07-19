"""Tests for in-plane dipole ACF and IR spectrum conventions."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest

_REPO = Path(__file__).resolve().parents[1]
_IR = _REPO / "aging_weak_lambda_ir"
if str(_IR) not in sys.path:
    sys.path.insert(0, str(_IR))

from analyze_ir_from_dipole import (  # noqa: E402
    dipole_acf_xy_direct,
    dipole_acf_xy_fft,
    ir_spectrum_dct_from_acf,
    ir_spectrum_mesa_fpe_from_dipole_xy,
)

C_LIGHT_CM_PS = 2.99792458e10


def _cosine_dipole(
    freq_cm1: float,
    duration_ps: float,
    dt_ps: float,
    dims: tuple[int, ...] = (0,),
) -> tuple[np.ndarray, np.ndarray]:
    times = np.arange(0.0, duration_ps, dt_ps, dtype=float)
    dipole = np.zeros((times.size, 3), dtype=float)
    omega = 2.0 * np.pi * freq_cm1 * C_LIGHT_CM_PS * 1e-12
    for dim in dims:
        dipole[:, dim] = np.cos(omega * times)
    return times, dipole


def test_fft_and_direct_xy_acf_agree_for_synthetic_signal() -> None:
    times, dipole = _cosine_dipole(1560.0, duration_ps=20.0, dt_ps=0.002, dims=(0, 1))
    acf_fft = dipole_acf_xy_fft(dipole)
    acf_direct = dipole_acf_xy_direct(dipole)
    n = min(acf_fft.size, acf_direct.size, 500)
    np.testing.assert_allclose(acf_fft[:n], acf_direct[:n], rtol=1e-10, atol=1e-10)


def test_xy_acf_excludes_z_contribution() -> None:
    times, dipole = _cosine_dipole(800.0, duration_ps=10.0, dt_ps=0.002, dims=(2,))
    acf = dipole_acf_xy_fft(dipole)
    assert acf[0] == pytest.approx(0.0, abs=1e-12)


def test_dct_on_acf_resolves_in_plane_mode() -> None:
    times, dipole = _cosine_dipole(1560.0, duration_ps=50.0, dt_ps=0.002, dims=(0,))
    acf = dipole_acf_xy_fft(dipole)
    dt_fs = float(np.mean(np.diff(times)) * 1000.0)
    freqs, spec = ir_spectrum_dct_from_acf(acf, dt_fs)
    mask = (freqs > 1000.0) & (freqs < 2000.0)
    peak = float(freqs[mask][int(np.argmax(spec[mask]))])
    assert abs(peak - 1560.0) < 40.0


def test_mesa_sums_x_and_y_component_spectra() -> None:
    pytest.importorskip("memspectrum")
    times, dipole = _cosine_dipole(1560.0, duration_ps=50.0, dt_ps=0.002, dims=(0, 1))
    freqs, spec = ir_spectrum_mesa_fpe_from_dipole_xy(dipole, times)
    mask = (freqs > 1000.0) & (freqs < 2000.0)
    peak = float(freqs[mask][int(np.argmax(spec[mask]))])
    assert abs(peak - 1560.0) < 40.0


def test_fdr_tracker_uses_unwrapped_positions() -> None:
    from hoomd.cavitymd.analysis.fdr import DipoleMomentFDRTracker

    tracker = DipoleMomentFDRTracker(time_tracker=SimpleNamespace(elapsed_time=0.0))
    simulation = MagicMock()
    snap = MagicMock()
    snap.particles.position = np.array([[0.1, 0.0, 0.0], [0.0, 0.0, 0.0]])
    snap.particles.image = np.array([[1, 0, 0], [0, 0, 0]])
    snap.particles.charge = np.array([1.0, 0.0])
    snap.global_box = SimpleNamespace(L=np.array([10.0, 10.0, 10.0]))
    simulation.state.cpu_local_snapshot.__enter__.return_value = snap
    simulation.state.cpu_local_snapshot.__exit__.return_value = None

    tracker.attach(simulation)
    dipole = tracker._calculate_total_dipole_moment()
    assert dipole[0] == pytest.approx(10.1)
