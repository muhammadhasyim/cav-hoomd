"""Tests for per-reference lambda=0 tau normalization."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SHARED = _REPO_ROOT / "reproduction" / "shared"
if str(_SHARED) not in sys.path:
    sys.path.insert(0, str(_SHARED))

from plot_fkt_analysis import (  # noqa: E402
    _format_lambda_tick_label,
    build_lambda0_tau_baseline_by_ref,
    merge_taus_with_kww_fallback,
    normalize_tau_by_lambda0_per_ref,
)


def test_build_lambda0_tau_baseline_by_ref_skips_invalid_entries():
    raw_taus = {
        ("lambda0", 1): 100.0,
        ("lambda0", 2): 0.0,
        ("lambda0", 3): -5.0,
        ("lambda0.03", 1): 300.0,
    }
    baseline = build_lambda0_tau_baseline_by_ref(
        raw_taus, ref_coupling_name="lambda0", refs=[1, 2, 3, 4]
    )
    assert baseline == {1: 100.0}


def test_normalize_tau_lambda0_is_unity():
    ref_tau_by_ref = {1: 100.0, 2: 120.0}
    assert normalize_tau_by_lambda0_per_ref(100.0, 1, ref_tau_by_ref) == pytest.approx(1.0)
    assert normalize_tau_by_lambda0_per_ref(120.0, 2, ref_tau_by_ref) == pytest.approx(1.0)


def test_normalize_tau_finite_ratio_for_nonzero_coupling():
    ref_tau_by_ref = {1: 100.0}
    assert normalize_tau_by_lambda0_per_ref(250.0, 1, ref_tau_by_ref) == pytest.approx(2.5)


def test_normalize_tau_missing_ref_returns_nan():
    ref_tau_by_ref = {1: 100.0}
    assert math.isnan(normalize_tau_by_lambda0_per_ref(200.0, 99, ref_tau_by_ref))


def test_normalize_tau_nonfinite_input_returns_nan():
    ref_tau_by_ref = {1: 100.0}
    assert math.isnan(normalize_tau_by_lambda0_per_ref(float("nan"), 1, ref_tau_by_ref))


def test_format_lambda_tick_label_three_sig_figs():
    assert _format_lambda_tick_label(0.0, "0") == r"$0$"
    assert _format_lambda_tick_label(0.01, "0.01") == r"$0.01$"
    assert _format_lambda_tick_label(0.016667, "0.016667") == r"$0.0167$"
    assert _format_lambda_tick_label(0.023333, "0.023333") == r"$0.0233$"


def test_merge_taus_with_kww_fallback_prefers_direct_crossing():
    raw = {("lam", 1): 100.0}
    kww = {
        ("lam", 1): {"tau_s": 999.0, "method": "kww_shared_beta"},
        ("lam", 2): {"tau_s": 200.0, "method": "kww_shared_beta"},
        ("lam", 3): {"tau_s": float("nan"), "method": "failed"},
    }
    merged = merge_taus_with_kww_fallback(raw, kww)
    assert merged[("lam", 1)] == pytest.approx(100.0)
    assert merged[("lam", 2)] == pytest.approx(200.0)
    assert ("lam", 3) not in merged
