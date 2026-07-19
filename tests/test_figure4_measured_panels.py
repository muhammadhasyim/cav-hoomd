"""Tests for measured F(k,t) panels (c)/(d) in Figure 4."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "reproduction" / "figure4_material_time"))
sys.path.insert(0, str(REPO_ROOT / "reproduction" / "csv_export"))

from plot_figure4 import build_measured_tau_tilde_tables  # noqa: E402
from export_figure4_csv import build_panels_cd  # noqa: E402


def _fake_relaxation_rows() -> list[tuple[float, int, float, float]]:
    """Minimal table: lambda=0 and 0.01, refs 1-9, tw up to 1600."""
    rows: list[tuple[float, int, float, float]] = []
    for lam in (0.0, 0.01):
        for ref in range(1, 10):
            ref_time = ref * 200.0
            tau = 180.0 + ref * 5.0 + lam * 100.0
            rows.append((lam, ref, ref_time, tau))
    return rows


def test_build_measured_tau_tilde_tables_reaches_tw_1600() -> None:
    rows = _fake_relaxation_rows()
    panel_c, panel_d = build_measured_tau_tilde_tables(rows)
    assert panel_c
    assert panel_d
    assert 1600.0 in panel_c
    assert 0.0 in panel_d
    assert float(panel_d[0.0]["t_w"][-1]) == pytest.approx(1600.0)
    assert float(panel_d[0.01]["t_w"][-1]) == pytest.approx(1600.0)


def test_build_panels_cd_and_measured_tables_agree_on_max_tw() -> None:
    rows = _fake_relaxation_rows()
    cd_rows = list(build_panels_cd(rows))
    panel_c, panel_d = build_measured_tau_tilde_tables(rows)
    csv_max_tw = max(float(r["x"]) for r in cd_rows if r["panel"] == "d")
    plot_max_tw = max(float(panel_d[lam]["t_w"][-1]) for lam in panel_d)
    assert csv_max_tw == pytest.approx(1600.0)
    assert plot_max_tw == pytest.approx(1600.0)
    assert len(panel_c[1600.0]["lambda"]) >= 2


def test_weak_lambda_profile_max_coupling() -> None:
    """CSV verifier should accept lambda=0.03 for the weak-lambda preliminary profile."""
    sys.path.insert(0, str(REPO_ROOT / "reproduction" / "shared"))
    from dataset_profile import load_profile  # noqa: WPS433

    profile = load_profile("aging_weak_lambda_preliminary")
    lam_max = max(entry.axis_value for entry in profile.couplings)
    assert lam_max == pytest.approx(0.03)
    assert abs(lam_max - 0.03) < 0.002
