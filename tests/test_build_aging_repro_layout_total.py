"""Tests for potential-energy staging (Total = PE, not KE+PE)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from reproduction.adapters.build_aging_repro_layout import (
    _load_energy_csv,
    _write_potential_energy,
)


def _write_energy_csv(
    path: Path,
    *,
    harmonic: np.ndarray,
    lj_coul: np.ndarray,
    system_total: np.ndarray,
    lj: np.ndarray | None = None,
    coulombic: np.ndarray | None = None,
) -> None:
    n = len(harmonic)
    time_ps = np.arange(n, dtype=float)
    header = "time_ps,harmonic_ha,lj_coul_ha,system_total_ha,universe_total_ha,kinetic_temp_K"
    if lj is not None and coulombic is not None:
        header = (
            "time_ps,harmonic_ha,lj_ha,coulombic_ha,lj_coul_ha,"
            "system_total_ha,universe_total_ha,kinetic_temp_K"
        )
        cols = np.column_stack(
            [time_ps, harmonic, lj, coulombic, lj_coul, system_total, system_total, np.full(n, 100.0)]
        )
    else:
        cols = np.column_stack(
            [time_ps, harmonic, lj_coul, system_total, system_total, np.full(n, 100.0)]
        )
    np.savetxt(path, cols, delimiter=",", header=header, comments="")


def test_load_energy_csv_total_is_potential_not_system_total(tmp_path: Path) -> None:
    """total_ha must be harmonic + lj_coul, never system_total (KE+PE)."""
    harmonic = np.array([0.04, 0.05], dtype=float)
    lj_coul = np.array([-1.0, -0.9], dtype=float)
    system_total = np.array([0.5, 0.6], dtype=float)  # includes KE offset
    csv_path = tmp_path / "avg.csv"
    _write_energy_csv(
        csv_path,
        harmonic=harmonic,
        lj_coul=lj_coul,
        system_total=system_total,
    )

    loaded = _load_energy_csv(csv_path)
    expected_total = harmonic + lj_coul
    assert np.allclose(loaded["total_ha"], expected_total)
    assert not np.allclose(loaded["total_ha"], system_total)


def test_load_energy_csv_reads_separate_lj_and_coulombic(tmp_path: Path) -> None:
    harmonic = np.array([0.04], dtype=float)
    lj = np.array([-0.69], dtype=float)
    coul = np.array([1.20], dtype=float)
    lj_coul = lj + coul
    system_total = np.array([2.0], dtype=float)
    csv_path = tmp_path / "avg.csv"
    _write_energy_csv(
        csv_path,
        harmonic=harmonic,
        lj_coul=lj_coul,
        system_total=system_total,
        lj=lj,
        coulombic=coul,
    )

    loaded = _load_energy_csv(csv_path)
    assert float(loaded["lj_ha"][0]) == pytest.approx(-0.69)
    assert float(loaded["coulombic_ha"][0]) == pytest.approx(1.20)
    assert np.allclose(loaded["total_ha"], harmonic + lj + coul)


def test_write_potential_energy_writes_pe_total_and_split_columns(tmp_path: Path) -> None:
    harmonic = np.array([0.04, 0.05], dtype=float)
    lj = np.array([-0.69, -0.68], dtype=float)
    coul = np.array([1.20, 1.18], dtype=float)
    data = {
        "time_ps": np.array([1.0, 2.0]),
        "harmonic_ha": harmonic,
        "lj_ha": lj,
        "coulombic_ha": coul,
        "lj_coul_ha": lj + coul,
        "total_ha": harmonic + lj + coul,
    }
    out = tmp_path / "pe.txt"
    _write_potential_energy(out, data)
    text = out.read_text(encoding="utf-8")
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    assert lines[0].startswith("time_ps")
    row = [float(x) for x in lines[1].split()]
    assert row[1] == pytest.approx(0.04)  # harmonic
    assert row[2] == pytest.approx(-0.69)  # lj
    assert row[3] == pytest.approx(1.20)  # coul
    assert row[4] == pytest.approx(0.04 + (-0.69) + 1.20)  # total PE
