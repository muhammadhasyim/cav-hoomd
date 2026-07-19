"""Tests for averaging all available 1600 ps and ~2000 ps replicas together."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from aging_weak_lambda.analysis.average_aging_energies import average_trimmed_series
from aging_weak_lambda.analysis.preliminary_repro_layout import list_union_complete_replicas
from examples.slurm.aging_campaign_status import fkt_file_path, write_protocol_marker
from reproduction.figure2_fkt.process_fskt_only import choose_detected_max_lag
from tests.test_aging_campaign_status import TEST_PROVENANCE, _write_complete_replica


def _write_fkt(path: Path, reference_time_ps: float, max_lag_ps: float) -> None:
    rows = "0.000000\t1.000000\n"
    if max_lag_ps > 0.0:
        rows += f"{max_lag_ps:.6f}\t0.100000\n"
    path.write_text(
        "# F(k,t) correlation function\n"
        f"# Reference time: {reference_time_ps:.3f} ps\n"
        "# lag_time_ps\tF(k,t)\n"
        f"{rows}",
        encoding="utf-8",
    )


def test_list_union_complete_replicas_includes_primary_and_near_target(
    tmp_path: Path,
) -> None:
    """1600/8 and 1999/10 replicas both contribute; incomplete ones do not."""
    lam = 0.01
    directory = tmp_path / "run"
    directory.mkdir()
    min_bytes = 5_000_000

    _write_complete_replica(
        directory,
        1,
        min_bytes=min_bytes,
        runtime_ps=1600.0,
        reference_interval_ps=200.0,
        max_references=8,
    )
    write_protocol_marker(
        directory,
        replica=1,
        lam=lam,
        switch_time_ps=200.0,
        seed=2,
        provenance=TEST_PROVENANCE,
    )

    # Near-target extension replica: scientifically complete F(k,t), but no
    # protocol marker (common for extension jobs). Must still enter the union.
    _write_complete_replica(
        directory,
        2,
        min_bytes=min_bytes,
        runtime_ps=1999.0,
        reference_interval_ps=200.0,
        max_references=10,
    )

    (directory / "observables_replica_3.h5").write_bytes(b"x" * min_bytes)
    for ref_index in range(4):
        _write_fkt(
            fkt_file_path(directory, 3, ref_index),
            ref_index * 200.0,
            1600.0 - ref_index * 200.0,
        )

    replicas = list_union_complete_replicas(
        directory,
        lam,
        n_replicas=10,
        switch_time_ps=200.0,
    )
    assert replicas == [1, 2]


def test_list_union_accepts_reference_time_drift_on_target_tier(
    tmp_path: Path,
) -> None:
    """Target tier tolerates ~5 ps header drift on late references."""
    directory = tmp_path / "run"
    directory.mkdir()
    min_bytes = 5_000_000
    (directory / "observables_replica_0.h5").write_bytes(b"x" * min_bytes)
    runtime = 1999.0
    for ref_index in range(10):
        # Production late refs drift to ~1805 ps vs policy 1800 ps.
        header_ref = 1805.0 if ref_index == 9 else float(ref_index * 200)
        _write_fkt(
            fkt_file_path(directory, 0, ref_index),
            header_ref,
            runtime - header_ref,
        )

    replicas = list_union_complete_replicas(
        directory,
        0.0,
        n_replicas=5,
        switch_time_ps=200.0,
    )
    assert replicas == [0]


def test_average_trimmed_series_uses_longest_time_grid() -> None:
    """Late-time points keep the longer replicas instead of truncating to the first."""
    short_t = np.array([0.0, 1.0, 2.0], dtype=float)
    long_t = np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=float)
    short = (
        short_t,
        np.ones_like(short_t),
        np.ones_like(short_t),
        np.ones_like(short_t) * 2,
        np.ones_like(short_t) * 3,
        np.ones_like(short_t) * 4,
        np.ones_like(short_t) * 5,
        np.ones_like(short_t) * 100,
    )
    long = (
        long_t,
        np.ones_like(long_t) * 10,
        np.ones_like(long_t) * 20,
        np.ones_like(long_t) * 30,
        np.ones_like(long_t) * 40,
        np.ones_like(long_t) * 50,
        np.ones_like(long_t) * 60,
        np.ones_like(long_t) * 100,
    )
    averaged = average_trimmed_series([short, long])
    assert averaged is not None
    assert float(averaged["time_ps"][-1]) == pytest.approx(4.0)
    assert int(averaged["coverage"][0]) == 2
    assert int(averaged["coverage"][-1]) == 1
    assert float(averaged["harmonic"][-1]) == pytest.approx(10.0)
    # Component order: harm, lj, coul, lj_coul, system, universe, kinetic
    assert float(averaged["lj_coul"][-1]) == pytest.approx(40.0)
    assert float(averaged["system_total"][-1]) == pytest.approx(50.0)
    assert float(averaged["universe_total"][-1]) == pytest.approx(60.0)
    assert float(averaged["kinetic_temp"][-1]) == pytest.approx(100.0)
    assert float(averaged["kinetic_temp"][0]) == pytest.approx(100.0)


def test_choose_detected_max_lag_union_vs_conservative() -> None:
    times = [100.0, 400.0, 250.0]
    assert choose_detected_max_lag(times, mode="min") == pytest.approx(100.0)
    assert choose_detected_max_lag(times, mode="max") == pytest.approx(400.0)
