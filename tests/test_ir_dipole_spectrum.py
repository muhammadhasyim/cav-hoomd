"""Unit tests for IR spectrum from fine-sampled dipole trajectories."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
_IR_DIR = _REPO_ROOT / "aging_weak_lambda_ir"
if str(_IR_DIR) not in sys.path:
    sys.path.insert(0, str(_IR_DIR))

from analyze_ir_from_dipole import (  # noqa: E402
    dipole_autocorrelation_fftconvolve,
    ir_spectrum_dct_from_dipole,
    ir_spectrum_from_dipole,
    ir_spectrum_mesa_fpe_from_dipole,
)

C_LIGHT_CM_PS = 2.99792458e10  # cm/s; lag in ps -> frequency in cm^-1


def _synthetic_dipole_cosine(
    freq_cm1: float,
    duration_ps: float,
    dt_ps: float,
    amplitude: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Dipole along x: mu(t) = A cos(2 pi nu t) with t in ps, nu in cm^-1."""
    times = np.arange(0.0, duration_ps, dt_ps, dtype=float)
    omega_rad_ps = 2.0 * np.pi * freq_cm1 * C_LIGHT_CM_PS * 1e-12
    dipole = np.zeros((times.size, 3), dtype=float)
    dipole[:, 0] = amplitude * np.cos(omega_rad_ps * times)
    return times, dipole


def test_acf_truncation_uses_early_fraction_only() -> None:
    dt_ps = 0.002
    times, dipole = _synthetic_dipole_cosine(1560.0, duration_ps=10.0, dt_ps=dt_ps)
    acf = dipole_autocorrelation_fftconvolve(dipole)
    assert acf.shape[0] == dipole.shape[0]
    freqs_short, _, _ = ir_spectrum_from_dipole(dipole, times, acf_fraction=0.01)
    freqs_full, _, _ = ir_spectrum_from_dipole(dipole, times, acf_fraction=1.0)
    assert freqs_full.size > freqs_short.size


def test_ir_spectrum_resolves_1560_cm1_peak_at_2fs_sampling() -> None:
    dt_ps = 0.002
    times, dipole = _synthetic_dipole_cosine(1560.0, duration_ps=50.0, dt_ps=dt_ps)
    freqs, spec, method = ir_spectrum_from_dipole(dipole, times, acf_fraction=1.0)
    assert method == "dct"
    assert freqs.size > 0
    mask = (freqs > 1000.0) & (freqs < 2000.0)
    assert np.any(mask)
    idx_peak = int(np.argmax(spec[mask]))
    freqs_band = freqs[mask]
    assert abs(freqs_band[idx_peak] - 1560.0) < 30.0


def test_ir_spectrum_resolves_4800_cm1_peak_at_2fs_sampling() -> None:
    dt_ps = 0.002
    times, dipole = _synthetic_dipole_cosine(4800.0, duration_ps=50.0, dt_ps=dt_ps)
    freqs, spec, method = ir_spectrum_from_dipole(dipole, times, acf_fraction=1.0)
    assert method == "dct"
    nyquist = 0.5 / (dt_ps * 1e-12) / C_LIGHT_CM_PS
    assert nyquist > 4800.0
    mask = (freqs > 4500.0) & (freqs < 5100.0)
    assert np.any(mask)
    idx_peak = int(np.argmax(spec[mask]))
    freqs_band = freqs[mask]
    assert abs(freqs_band[idx_peak] - 4800.0) < 50.0


def test_mesa_fpe_resolves_1560_cm1_peak_at_2fs_sampling() -> None:
    pytest = __import__("pytest")
    try:
        import memspectrum  # noqa: F401
    except ImportError:
        pytest.skip("memspectrum not installed")

    dt_ps = 0.002
    times, dipole = _synthetic_dipole_cosine(1560.0, duration_ps=50.0, dt_ps=dt_ps)
    freqs, spec = ir_spectrum_mesa_fpe_from_dipole(dipole, times)
    assert freqs.size > 0
    mask = (freqs > 1000.0) & (freqs < 2000.0)
    assert np.any(mask)
    idx_peak = int(np.argmax(spec[mask]))
    freqs_band = freqs[mask]
    assert abs(freqs_band[idx_peak] - 1560.0) < 30.0


def test_auto_prefers_mesa_when_dct_invalid(monkeypatch) -> None:
    pytest = __import__("pytest")
    try:
        import memspectrum  # noqa: F401
    except ImportError:
        pytest.skip("memspectrum not installed")

    dt_ps = 0.002
    times, dipole = _synthetic_dipole_cosine(1560.0, duration_ps=50.0, dt_ps=dt_ps)

    def _bad_dct(*args, **kwargs):
        return np.array([]), np.array([])

    monkeypatch.setattr(
        "analyze_ir_from_dipole.ir_spectrum_dct_from_dipole",
        _bad_dct,
    )
    freqs, spec, method = ir_spectrum_from_dipole(
        dipole,
        times,
        method="auto",
        acf_fraction=1.0,
    )
    assert method == "mesa"
    assert freqs.size > 0
    assert spec.size > 0
