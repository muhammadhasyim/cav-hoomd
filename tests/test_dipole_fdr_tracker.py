"""Tests for DipoleMomentFDRTracker attach behavior."""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from hoomd.cavitymd.analysis.fdr import DipoleMomentFDRTracker


class _DummyTimeTracker:
    elapsed_time = 0.0


def test_dipole_fdr_tracker_uses_attached_simulation() -> None:
    tracker = DipoleMomentFDRTracker(time_tracker=_DummyTimeTracker())
    simulation = MagicMock()
    snap = MagicMock()
    snap.particles.position = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    snap.particles.image = np.zeros((2, 3), dtype=int)
    snap.particles.charge = np.array([1.0, 0.0])
    snap.global_box = SimpleNamespace(L=np.array([10.0, 10.0, 10.0]))
    simulation.state.cpu_local_snapshot.__enter__.return_value = snap
    simulation.state.cpu_local_snapshot.__exit__.return_value = None

    tracker.attach(simulation)
    dipole = tracker._calculate_total_dipole_moment()

    assert dipole.shape == (3,)
    assert dipole[0] == pytest.approx(1.0)


def test_dipole_fdr_tracker_requires_attachment() -> None:
    tracker = DipoleMomentFDRTracker(time_tracker=_DummyTimeTracker())
    with pytest.raises(AttributeError, match="not attached"):
        tracker._calculate_total_dipole_moment()
