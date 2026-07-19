"""Tests for Figure 4 TN-based panel assembly."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "reproduction" / "shared"))
sys.path.insert(0, str(REPO_ROOT / "reproduction" / "figure2_3_relaxation"))
sys.path.insert(0, str(REPO_ROOT / "reproduction" / "figure4_material_time"))

from plot_figure4 import (  # noqa: E402
    DEFAULT_BETA,
    build_tn_tau_tilde_tables,
    load_panel_b_csv,
    stretched_exponential_guide,
)
from plot_fictive_three_panel_analysis import FictiveThreePanelAnalyzer  # noqa: E402
from dataset_profile import load_profile  # noqa: E402


STAGED = REPO_ROOT / "aging_weak_lambda" / "derived" / "repro_layout"
FIGURES = REPO_ROOT / "aging_weak_lambda" / "derived" / "repro_figures"


@pytest.fixture(scope="module")
def analyzer() -> FictiveThreePanelAnalyzer:
    if not STAGED.is_dir():
        pytest.skip(f"staged root missing: {STAGED}")
    profile = load_profile("aging_weak_lambda")
    return FictiveThreePanelAnalyzer(
        base_dir=str(STAGED),
        use_latex=False,
        profile=profile,
    )


def test_load_panel_b_csv_groups_trajectories() -> None:
    path = FIGURES / "figure4_panel_b.csv"
    if not path.is_file():
        pytest.skip("figure4_panel_b.csv missing")
    data = load_panel_b_csv(path)
    assert data
    assert data[0]["h"].size == data[0]["phi"].size
    assert np.all(data[0]["phi"] <= 1.05)


def test_tn_tau_tilde_tables_are_smooth_and_normalized(analyzer: FictiveThreePanelAnalyzer) -> None:
    panel_c, panel_d = build_tn_tau_tilde_tables(analyzer)
    assert panel_c
    assert panel_d
    # lambda=0 baseline should sit at 1
    assert 0.0 in panel_d
    assert np.allclose(panel_d[0.0]["tau_tilde"], 1.0, atol=1e-6)
    # curves should be denser / smoother than sparse measured tables
    sample_tw = next(iter(panel_c))
    assert panel_c[sample_tw]["lambda"].size >= 3
    # no exploding values for weak-lambda campaign
    for data in panel_c.values():
        assert np.all(np.isfinite(data["tau_tilde"]))
        assert float(np.nanmax(data["tau_tilde"])) < 10.0


def test_default_beta_matches_weak_lambda_collapse() -> None:
    # Paper uses 0.55; this campaign's F=0.1 collapse overlays near ~0.72.
    assert DEFAULT_BETA == pytest.approx(0.72)


def test_stretched_exponential_guide_matches_f01_definition() -> None:
    """At h=1, guide must equal 0.1 for any beta (F/F0=0.1 criterion)."""
    assert stretched_exponential_guide(np.array([1.0]), 0.55)[0] == pytest.approx(0.1)
    assert stretched_exponential_guide(np.array([1.0]), 0.73)[0] == pytest.approx(0.1)
    assert stretched_exponential_guide(np.array([0.0]), 0.55)[0] == pytest.approx(1.0)
