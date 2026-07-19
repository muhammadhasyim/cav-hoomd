"""Tests for empirical Tf mapping with reaction-field PE-vs-T tables."""

from __future__ import annotations

from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RF_TABLE = (
    REPO_ROOT
    / "aging_weak_lambda"
    / "analysis"
    / "pe_vs_T_calib_rf"
    / "potential_energy_vs_T_rf.txt"
)


@pytest.mark.skipif(not RF_TABLE.is_file(), reason="RF calibration table not built yet")
def test_rf_lj_coulombic_recovers_calibration_temperatures() -> None:
    """Inverted RF U(T) must still invert to ~set-point T (interp, not T^{3/5})."""
    try:
        from hoomd.cavitymd.controllers.empirical import EmpiricalTemperatureData
    except ModuleNotFoundError:
        from cavitymd.controllers.empirical import EmpiricalTemperatureData

    emp = EmpiricalTemperatureData(
        str(RF_TABLE),
        energy_component="lj_coulombic",
        create_plots=False,
    )
    # Table row at 100 K: lj + coulombic
    energy_100k = -6.86101720e-01 + 1.22271306e00
    tf = emp.calculate_systemic_temperature(energy_100k)
    assert tf == pytest.approx(100.0, abs=2.0)

    energy_300k = -4.19664493e-01 + 2.30678922e-01
    tf_hot = emp.calculate_systemic_temperature(energy_300k)
    assert tf_hot == pytest.approx(300.0, abs=5.0)
