"""Tests for aging campaign part 2 submit watcher."""

from __future__ import annotations

from pathlib import Path

from examples.slurm.watch_submit_aging_part2 import (
    SubmitGate,
    part2_already_submitted,
)


def test_submit_gate_allows_when_under_limit() -> None:
    """Part 2 should submit when projected total stays below the QOS cap."""
    gate = SubmitGate(
        active_jobs=760,
        part2_tasks=1226,
        max_submit_jobs=2000,
        safety_margin=8,
    )
    assert gate.projected_total == 1986
    assert gate.allowed is True


def test_submit_gate_blocks_when_over_limit() -> None:
    """Part 2 should wait when projected total exceeds the QOS cap."""
    gate = SubmitGate(
        active_jobs=800,
        part2_tasks=1226,
        max_submit_jobs=2000,
        safety_margin=8,
    )
    assert gate.projected_total == 2026
    assert gate.allowed is False


def test_part2_already_submitted_marker(tmp_path: Path) -> None:
    """Marker file should short-circuit duplicate submissions."""
    submitted_log = tmp_path / "submitted.txt"
    marker = tmp_path / "part2.marker"
    marker.write_text("job_id=123\n", encoding="utf-8")
    assert part2_already_submitted(
        submitted_log=submitted_log,
        marker_file=marker,
    )


def test_part2_already_submitted_log(tmp_path: Path) -> None:
    """Submitted log entry should short-circuit duplicate submissions."""
    submitted_log = tmp_path / "submitted.txt"
    submitted_log.write_text(
        "14243455:packed-r1-defer24-part2:1226:walltime=01:00:00\n",
        encoding="utf-8",
    )
    marker = tmp_path / "part2.marker"
    assert part2_already_submitted(
        submitted_log=submitted_log,
        marker_file=marker,
    )
