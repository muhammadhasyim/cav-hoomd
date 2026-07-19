"""Tests for aging campaign protocol-marker backfill."""

from __future__ import annotations

from pathlib import Path

from examples.slurm.aging_campaign_status import (
    is_complete_replica,
    protocol_marker_path,
    run_dir,
)
from examples.slurm.backfill_aging_protocol_markers import (
    backfill_protocol_markers,
    find_protocol_backfill_candidates,
)
from tests.test_aging_campaign_status import TEST_PROVENANCE, _write_complete_replica


def test_find_protocol_backfill_candidates_detects_missing_marker(
    tmp_path: Path,
) -> None:
    """Candidates should include data-complete replicas without markers."""
    lam = 0.01
    directory = run_dir(tmp_path, lam, switch_time_ps=200.0)
    directory.mkdir(parents=True)
    _write_complete_replica(directory, 3, min_bytes=1024)

    candidates = find_protocol_backfill_candidates(
        tmp_path,
        (lam,),
        n_replicas=10,
        runtime_ps=14.0,
        reference_interval_ps=1.0,
        max_references=13,
        switch_time_ps=200.0,
        min_bytes=1024,
    )
    assert candidates == [(lam, 3)]


def test_backfill_protocol_markers_writes_valid_marker(tmp_path: Path) -> None:
    """Backfill should promote protocol-only partials to complete."""
    lam = 0.03
    directory = run_dir(tmp_path, lam, switch_time_ps=200.0)
    directory.mkdir(parents=True)
    _write_complete_replica(directory, 5, min_bytes=1024)

    written = backfill_protocol_markers(
        tmp_path,
        (lam,),
        provenance=TEST_PROVENANCE,
        n_replicas=10,
        runtime_ps=14.0,
        reference_interval_ps=1.0,
        max_references=13,
        switch_time_ps=200.0,
        min_bytes=1024,
        dry_run=False,
    )
    assert written == [(lam, 5)]
    assert protocol_marker_path(directory, 5).is_file()
    assert is_complete_replica(
        directory,
        5,
        min_bytes=1024,
        runtime_ps=14.0,
        reference_interval_ps=1.0,
        max_references=13,
        require_step_protocol=True,
        expected_lam=lam,
        expected_switch_time_ps=200.0,
    )


def test_backfill_protocol_markers_dry_run_is_non_mutating(tmp_path: Path) -> None:
    """Dry-run should report candidates without writing markers."""
    lam = 0.01
    directory = run_dir(tmp_path, lam, switch_time_ps=200.0)
    directory.mkdir(parents=True)
    _write_complete_replica(directory, 2, min_bytes=1024)

    written = backfill_protocol_markers(
        tmp_path,
        (lam,),
        provenance=TEST_PROVENANCE,
        n_replicas=10,
        runtime_ps=14.0,
        reference_interval_ps=1.0,
        max_references=13,
        switch_time_ps=200.0,
        min_bytes=1024,
        dry_run=True,
    )
    assert written == [(lam, 2)]
    assert not protocol_marker_path(directory, 2).is_file()
