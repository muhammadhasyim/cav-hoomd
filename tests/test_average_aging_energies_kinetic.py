"""Regression: kinetic temperature must not be aliased to universe_total."""

from __future__ import annotations

import numpy as np
import pytest

from aging_weak_lambda.analysis.average_aging_energies import average_trimmed_series


def test_average_trimmed_series_maps_kinetic_temperature_not_universe_energy() -> None:
    """load_replica order ends with kinetic_temp; that slot must stay in Kelvin."""
    time_ps = np.array([0.0, 1.0, 2.0], dtype=float)
    harmonic = np.array([0.04, 0.04, 0.05], dtype=float)
    lj = np.array([-0.7, -0.7, -0.69], dtype=float)
    coulombic = np.zeros_like(time_ps)
    lj_coul = lj.copy()
    system_total = np.array([0.1, 0.1, 0.12], dtype=float)
    universe_total = np.array([0.8, 0.9, 1.0], dtype=float)
    kinetic_temp = np.array([99.5, 100.2, 100.8], dtype=float)

    averaged = average_trimmed_series(
        [
            (
                time_ps,
                harmonic,
                lj,
                coulombic,
                lj_coul,
                system_total,
                universe_total,
                kinetic_temp,
            )
        ]
    )
    assert averaged is not None
    assert np.allclose(averaged["kinetic_temp"], kinetic_temp)
    assert np.allclose(averaged["universe_total"], universe_total)
    assert np.allclose(averaged["lj_coul"], lj_coul)
    assert float(np.mean(averaged["kinetic_temp"])) == pytest.approx(100.1666667, abs=1e-5)
    # Guard against the old off-by-one that wrote universe_total into kinetic_temp.
    assert float(np.mean(averaged["kinetic_temp"])) > 50.0
