"""Tests for converting 1/e equilibrium tau(T) tables to F/F0=0.1."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    REPO_ROOT
    / "reproduction"
    / "calibration"
    / "convert_relaxation_table_to_f01.py"
)

import sys

sys.path.insert(0, str(SCRIPT.parent))

from convert_relaxation_table_to_f01 import (  # noqa: E402
    convert_tau_1e_to_f01,
    convert_table,
    kww_1e_to_f01_factor,
)


def test_kww_factor_at_beta_055() -> None:
    factor = kww_1e_to_f01_factor(0.55)
    assert factor == pytest.approx((-np.log(0.1)) ** (1.0 / 0.55), rel=1e-12)
    assert 4.0 < factor < 5.0


def test_convert_tau_scales_monotonically() -> None:
    tau_1e = np.array([10.0, 30.0, 100.0])
    tau_01 = convert_tau_1e_to_f01(tau_1e, beta=0.55)
    assert np.all(tau_01 > tau_1e)
    assert np.all(np.diff(tau_01) > 0.0)
    # At beta=0.55, Phi=exp(-(t/tau)^beta): t(0.1)/t(1/e) is the factor
    assert tau_01[1] / tau_1e[1] == pytest.approx(kww_1e_to_f01_factor(0.55))


def test_convert_table_preserves_temperature_grid(tmp_path: Path) -> None:
    src = tmp_path / "tau_1e.txt"
    src.write_text(
        "# Relaxation time analysis from F(k,t) using 1/e method\n"
        "# Temperature(K) 1/T(1/K) tau_relax(ps) F_initial F_final_norm decay_extent method success\n"
        "100.00 0.010000 29.0 700.0 0.0 1.0 interpolated True\n"
        "120.00 0.008333 10.0 700.0 0.0 1.0 interpolated True\n",
        encoding="utf-8",
    )
    dst = tmp_path / "tau_f01.txt"
    rows = convert_table(src, dst, beta=0.55)
    assert len(rows) == 2
    assert rows[0][0] == pytest.approx(100.0)
    assert rows[0][2] == pytest.approx(29.0 * kww_1e_to_f01_factor(0.55))
    text = dst.read_text(encoding="utf-8")
    assert "F/F0 = 0.1" in text
    assert "1/e" in text
    assert "beta=0.55" in text


def test_anchor_matches_reference_temperature(tmp_path: Path) -> None:
    src = tmp_path / "tau_1e.txt"
    src.write_text(
        "100.00 0.010000 29.0 700.0 0.0 1.0 interpolated True\n"
        "120.00 0.008333 10.0 700.0 0.0 1.0 interpolated True\n",
        encoding="utf-8",
    )
    dst = tmp_path / "tau_f01.txt"
    rows = convert_table(
        src,
        dst,
        beta=0.55,
        anchor_temperature_k=100.0,
        anchor_tau_ps=177.0,
    )
    assert rows[0][2] == pytest.approx(177.0, rel=1e-9)
    # Relative T-dependence preserved
    assert rows[1][2] / rows[0][2] == pytest.approx(10.0 / 29.0, rel=1e-9)
