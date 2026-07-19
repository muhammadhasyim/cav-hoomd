"""Helpers for polling SLURM job completion in aging campaign watchers."""

from __future__ import annotations

import re
import subprocess


def parse_marker_job_id(marker_text: str) -> str:
    """Return a SLURM job ID from a watcher marker file body."""
    match = re.search(r"(?:^|\n)job_id=(\d+)(?:\n|$)", marker_text.strip())
    if match is None:
        raise ValueError("marker does not contain job_id=<digits>")
    return match.group(1)


def count_active_job_tasks(
    job_id: str,
    *,
    user: str | None = None,
) -> int:
    """Count PENDING/RUNNING array tasks for one SLURM job ID."""
    command = [
        "squeue",
        "-j",
        job_id,
        "-h",
        "-t",
        "PENDING,RUNNING",
        "-r",
    ]
    if user is not None:
        command[1:1] = ["-u", user]
    result = subprocess.run(
        command,
        check=True,
        capture_output=True,
        text=True,
    )
    return len([line for line in result.stdout.splitlines() if line.strip()])


def jobs_fully_complete(job_ids: list[str], *, user: str | None = None) -> bool:
    """Return True when none of the job IDs have active array tasks."""
    return all(
        count_active_job_tasks(job_id, user=user) == 0 for job_id in job_ids
    )
