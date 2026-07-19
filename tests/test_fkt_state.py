"""Tests for F(k,t) state serialization."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from cavitymd.analysis.fkt_state import load_fkt_state, save_fkt_state


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
    assert loaded["references"][0]["time_ps"] == 1.0
