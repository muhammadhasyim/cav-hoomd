"""Tests for SLURM watcher helpers and strong-adaptive submit script."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from examples.slurm.slurm_job_watch import jobs_fully_complete, parse_marker_job_id


def test_parse_marker_job_id_reads_job_id() -> None:
    assert parse_marker_job_id("job_id=14250221\n") == "14250221"


def test_parse_marker_job_id_rejects_missing_value() -> None:
    with pytest.raises(ValueError, match="job_id"):
        parse_marker_job_id("submitted=true\n")


def test_jobs_fully_complete_true_when_all_counts_zero(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_count(job_id: str, *, user: str | None = None) -> int:
        return 0

    monkeypatch.setattr(
        "examples.slurm.slurm_job_watch.count_active_job_tasks",
        fake_count,
    )
    assert jobs_fully_complete(["1", "2"]) is True


def test_jobs_fully_complete_false_when_any_active(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    active = {"14250221": 0, "14260000": 3}

    def fake_count(job_id: str, *, user: str | None = None) -> int:
        return active[job_id]

    monkeypatch.setattr(
        "examples.slurm.slurm_job_watch.count_active_job_tasks",
        fake_count,
    )
    assert jobs_fully_complete(["14250221", "14260000"]) is False


def test_submit_strong_adaptive_plan_only_uses_weak_lambda_walltime() -> None:
    repository = Path(__file__).resolve().parents[1]
    script = repository / "examples" / "slurm" / "submit_aging_strong_adaptive.sh"
    plan_id = "test-plan-only"
    plan_dir = repository / "examples" / "slurm" / "generated" / f"aging_packed_{plan_id}"

    result = subprocess.run(
        ["bash", str(script), "--plan-only", "--plan-id", plan_id],
        cwd=repository / "examples" / "slurm",
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    try:
        sbatch = (plan_dir / "job_part1.sbatch").read_text(encoding="utf-8")
        assert "#SBATCH --time=02:00:00" in sbatch
        assert "#SBATCH --array=" in sbatch
        assert "%24" in sbatch
        assert "#SBATCH --cpus-per-task=8" in sbatch
        assert "#SBATCH --mem=32G" in sbatch
        assert "--adaptive-timestep" in sbatch
        assert "--error-tolerance 5.0" in sbatch
        assert "--initial-fraction 1e-5" in sbatch
        assert "--time-constant-ps 50.0" in sbatch
        assert "--fixed-timestep" not in sbatch
    finally:
        if plan_dir.is_dir():
            shutil.rmtree(plan_dir)
