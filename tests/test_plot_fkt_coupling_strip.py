"""Tests for paper-style F(k,t) coupling-strip helpers and layout."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SHARED = _REPO_ROOT / "reproduction" / "shared"
if str(_SHARED) not in sys.path:
    sys.path.insert(0, str(_SHARED))

from plot_fkt_analysis import (  # noqa: E402
    _apply_measured_relaxation_rcparams,
    colorbar_ticks_every,
    plot_fkt_by_coupling,
    waiting_time_ps_from_ref,
)


def test_waiting_time_ps_from_ref_uses_200_ps_interval():
    assert waiting_time_ps_from_ref(0) == pytest.approx(0.0)
    assert waiting_time_ps_from_ref(1) == pytest.approx(200.0)
    assert waiting_time_ps_from_ref(12) == pytest.approx(2400.0)


def test_colorbar_ticks_every_400_ps():
    ticks = colorbar_ticks_every(0.0, 2400.0, step=400.0)
    assert ticks == [0.0, 400.0, 800.0, 1200.0, 1600.0, 2000.0, 2400.0]


def test_colorbar_ticks_clips_to_range():
    ticks = colorbar_ticks_every(100.0, 900.0, step=400.0)
    assert ticks == [400.0, 800.0]


def test_apply_measured_relaxation_rcparams_uses_computer_modern():
    _apply_measured_relaxation_rcparams(use_latex=False)
    assert plt.rcParams["mathtext.fontset"] == "cm"
    assert "cmr10" in plt.rcParams["font.family"] or plt.rcParams["font.family"] == "cmr10" or (
        isinstance(plt.rcParams["font.family"], list)
        and "cmr10" in plt.rcParams["font.family"]
    ) or plt.rcParams["font.family"] in ("cmr10", "serif")


def test_plot_fkt_by_coupling_strip_layout_and_outputs(tmp_path):
    """Synthetic 1×N strip: shared y, viridis tw, PNG+PDF, CM math labels."""
    t = np.linspace(0.0, 2000.0, 101)
    data_dict = {}
    coupling_info = {}
    for lam, name in [(0.0, "lambda_0"), (0.01, "lambda_0p01"), (0.03, "lambda_0p03")]:
        folder = {}
        for ref in range(3):
            # Simple exponential decay; slower for larger tw.
            tau = 200.0 + 50.0 * ref
            fkt = np.exp(-t / tau)
            folder[ref] = (t, fkt, None)
        data_dict[name] = folder
        coupling_info[name] = (lam, f"{lam:g}")

    # Bypass paper epsilon whitelist.
    import plot_fkt_analysis as pfa

    old = pfa._USE_PROFILE_COUPLINGS
    pfa._USE_PROFILE_COUPLINGS = True
    try:
        plot_fkt_by_coupling(data_dict, coupling_info, output_dir=str(tmp_path))
    finally:
        pfa._USE_PROFILE_COUPLINGS = old

    png = tmp_path / "fkt_by_coupling_filtered.png"
    pdf = tmp_path / "fkt_by_coupling_filtered.pdf"
    assert png.is_file()
    assert pdf.is_file()
    assert png.stat().st_size > 1000
    assert pdf.stat().st_size > 1000

    # CM / mathtext still configured after plot.
    assert plt.rcParams["mathtext.fontset"] == "cm" or plt.rcParams["text.usetex"] is True
