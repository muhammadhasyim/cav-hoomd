#!/usr/bin/env python3
"""Poll SLURM queue depth and submit aging campaign part 2 when QOS allows."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_MAX_SUBMIT_JOBS = 2000
DEFAULT_PART2_TASKS = 1226
DEFAULT_POLL_INTERVAL_S = 300
DEFAULT_SAFETY_MARGIN = 8


@dataclass(frozen=True)
class SubmitGate:
    """Decision inputs for whether part 2 can be submitted."""

    active_jobs: int
    part2_tasks: int
    max_submit_jobs: int
    safety_margin: int

    @property
    def projected_total(self) -> int:
        """Return active jobs plus the part 2 array size."""
        return self.active_jobs + self.part2_tasks

    @property
    def allowed(self) -> bool:
        """Return True when submitting part 2 stays under the QOS cap."""
        return self.projected_total <= (self.max_submit_jobs - self.safety_margin)


def count_active_jobs(*, user: str) -> int:
    """Count pending and running SLURM jobs for one user.

    Each array task is counted separately, matching gpu48 MaxSubmitPU behavior.
    """
    result = subprocess.run(
        ["squeue", "-u", user, "-h", "-t", "PENDING,RUNNING"],
        check=True,
        capture_output=True,
        text=True,
    )
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    return len(lines)


def part2_already_submitted(
    *,
    submitted_log: Path,
    marker_file: Path,
) -> bool:
    """Return True when part 2 was already launched or recorded."""
    if marker_file.is_file():
        return True
    if not submitted_log.is_file():
        return False
    text = submitted_log.read_text(encoding="utf-8")
    return "packed-r1-defer24-part2" in text


class SubmitResult:
    """Outcome of one part 2 submission attempt."""

    def __init__(
        self,
        *,
        status: str,
        job_id: str | None = None,
        message: str = "",
    ) -> None:
        self.status = status
        self.job_id = job_id
        self.message = message


def submit_part2(*, sbatch_file: Path) -> SubmitResult:
    """Submit part 2, returning success, quota-blocked, or hard failure."""
    result = subprocess.run(
        ["sbatch", "--parsable", str(sbatch_file)],
        capture_output=True,
        text=True,
    )
    stdout = result.stdout.strip()
    stderr = (result.stderr or "").strip()
    if result.returncode == 0:
        job_id = stdout.split(";")[0]
        if not job_id.isdigit():
            return SubmitResult(
                status="error",
                message=f"unexpected sbatch output: {stdout!r}",
            )
        return SubmitResult(status="submitted", job_id=job_id)

    combined = f"{stdout}\n{stderr}".strip()
    if "QOSMaxSubmitJobPerUserLimit" in combined:
        return SubmitResult(status="quota_blocked", message=combined)
    return SubmitResult(status="error", message=combined or f"exit {result.returncode}")


def append_submission_record(
    *,
    submitted_log: Path,
    job_id: str,
    part2_tasks: int,
    walltime: str,
) -> None:
    """Append one submission record to the campaign log."""
    submitted_log.parent.mkdir(parents=True, exist_ok=True)
    with submitted_log.open("a", encoding="utf-8") as handle:
        handle.write(
            f"{job_id}:packed-r1-defer24-part2:{part2_tasks}:walltime={walltime}\n"
        )


def write_marker(*, marker_file: Path, job_id: str) -> None:
    """Persist the successful part 2 submission."""
    marker_file.parent.mkdir(parents=True, exist_ok=True)
    marker_file.write_text(
        f"job_id={job_id}\n"
        f"submitted_at={datetime.now(timezone.utc).isoformat()}\n",
        encoding="utf-8",
    )


def log_message(log_file: Path, message: str) -> None:
    """Append one timestamped line to the watcher log."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with log_file.open("a", encoding="utf-8") as handle:
        handle.write(f"{timestamp} {message}\n")


def run_watcher(args: argparse.Namespace) -> int:
    """Poll until part 2 submits successfully or the watcher is stopped."""
    user = args.user or os.environ.get("USER", "")
    if not user:
        raise ValueError("cannot determine SLURM user")

    if part2_already_submitted(
        submitted_log=args.submitted_log,
        marker_file=args.marker_file,
    ):
        log_message(args.log_file, "part 2 already submitted; exiting")
        return 0

    deadline = time.monotonic() + args.max_runtime_s
    while True:
        active_jobs = count_active_jobs(user=user)
        gate = SubmitGate(
            active_jobs=active_jobs,
            part2_tasks=args.part2_tasks,
            max_submit_jobs=args.max_submit_jobs,
            safety_margin=args.safety_margin,
        )
        log_message(
            args.log_file,
            (
                f"active_jobs={active_jobs} projected_total={gate.projected_total} "
                f"limit={args.max_submit_jobs - args.safety_margin} "
                f"squeue_gate={gate.allowed}"
            ),
        )

        outcome = submit_part2(sbatch_file=args.sbatch_file)
        if outcome.status == "submitted" and outcome.job_id is not None:
            append_submission_record(
                submitted_log=args.submitted_log,
                job_id=outcome.job_id,
                part2_tasks=args.part2_tasks,
                walltime=args.walltime,
            )
            write_marker(marker_file=args.marker_file, job_id=outcome.job_id)
            log_message(
                args.log_file,
                f"submitted part 2 as job {outcome.job_id}; watcher exiting",
            )
            return 0
        if outcome.status == "quota_blocked":
            log_message(
                args.log_file,
                (
                    "submit blocked by QOSMaxSubmitJobPerUserLimit; "
                    "waiting for part 1 tasks to finish"
                ),
            )
        else:
            log_message(args.log_file, f"submit failed: {outcome.message}")

        if time.monotonic() >= deadline:
            log_message(args.log_file, "max runtime reached without submitting part 2")
            return 1

        time.sleep(args.poll_interval_s)


def build_parser() -> argparse.ArgumentParser:
    """Build the watcher CLI."""
    script_dir = Path(__file__).resolve().parent
    default_plan_dir = (
        script_dir / "generated" / "aging_packed_20260719T051155Z-r1-defer24"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sbatch-file",
        type=Path,
        default=default_plan_dir / "job_part2.sbatch",
    )
    parser.add_argument(
        "--submitted-log",
        type=Path,
        default=script_dir / "submitted_aging_campaign.txt",
    )
    parser.add_argument(
        "--marker-file",
        type=Path,
        default=default_plan_dir / "part2_submitted.marker",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=script_dir / "logs" / "watch_submit_part2.log",
    )
    parser.add_argument("--user", type=str, default=None)
    parser.add_argument(
        "--max-submit-jobs",
        type=int,
        default=DEFAULT_MAX_SUBMIT_JOBS,
    )
    parser.add_argument(
        "--part2-tasks",
        type=int,
        default=DEFAULT_PART2_TASKS,
    )
    parser.add_argument(
        "--safety-margin",
        type=int,
        default=DEFAULT_SAFETY_MARGIN,
    )
    parser.add_argument(
        "--poll-interval-s",
        type=int,
        default=DEFAULT_POLL_INTERVAL_S,
    )
    parser.add_argument(
        "--max-runtime-s",
        type=int,
        default=7 * 24 * 3600,
        help="Stop after this many seconds if part 2 never submits.",
    )
    parser.add_argument(
        "--walltime",
        type=str,
        default="01:00:00",
    )
    return parser


def main() -> int:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()
    if args.part2_tasks <= 0:
        parser.error("--part2-tasks must be positive")
    if args.max_submit_jobs <= 0:
        parser.error("--max-submit-jobs must be positive")
    if args.poll_interval_s <= 0:
        parser.error("--poll-interval-s must be positive")
    if not args.sbatch_file.is_file():
        parser.error(f"sbatch file not found: {args.sbatch_file}")
    return run_watcher(args)


if __name__ == "__main__":
    sys.exit(main())
