#!/usr/bin/env python3
"""Hourly watcher that bundles incomplete aging replicas into retry job arrays."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Protocol, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from examples.slurm.aging_campaign_status import (
        COMPLETE_MIN_BYTES,
        DEFAULT_FKT_MAX_REFERENCES,
        DEFAULT_FKT_REFERENCE_INTERVAL_PS,
        DEFAULT_RUNTIME_PS,
        DEFAULT_SWITCH_TIME_PS,
        fkt_file_path,
        is_complete_replica,
        lambda_to_tag,
        protocol_marker_path,
        replica_gsd_path,
        replica_h5_path,
        run_dir,
        scan_lambda_run,
    )
else:
    from .aging_campaign_status import (
        COMPLETE_MIN_BYTES,
        DEFAULT_FKT_MAX_REFERENCES,
        DEFAULT_FKT_REFERENCE_INTERVAL_PS,
        DEFAULT_RUNTIME_PS,
        DEFAULT_SWITCH_TIME_PS,
        fkt_file_path,
        is_complete_replica,
        lambda_to_tag,
        protocol_marker_path,
        replica_gsd_path,
        replica_h5_path,
        run_dir,
        scan_lambda_run,
    )


DEFAULT_LAMBDAS = (0.0, 0.01, 0.016667, 0.023333, 0.03)
DEFAULT_TARGET_VALID = 500
DEFAULT_POLL_INTERVAL_S = 3600
DEFAULT_MAX_SUBMIT_JOBS = 2000
DEFAULT_SAFETY_MARGIN = 8
DEFAULT_ARRAY_CONCURRENCY = 24
DEFAULT_JOB_NAME = "cav-aging-retry"
DEFAULT_PARTITION = "l40s_public,h200_public"
DEFAULT_ACCOUNT = "torch_pr_283_chemistry"
DEFAULT_WALLTIME = "01:00:00"


class CommandRunner(Protocol):
    """Minimal subprocess interface for tests."""

    def __call__(
        self,
        args: Sequence[str],
        *,
        check: bool,
        capture_output: bool,
        text: bool,
    ) -> subprocess.CompletedProcess[str]:
        """Run one subprocess command."""


@dataclass(frozen=True)
class RetryCandidate:
    """One incomplete replica eligible for retry."""

    lambda_tag: str
    lam: float
    replica: int


@dataclass(frozen=True)
class RetrySubmitGate:
    """Decision inputs for whether a retry array can be submitted."""

    active_jobs: int
    retry_tasks: int
    max_submit_jobs: int
    safety_margin: int

    @property
    def projected_total(self) -> int:
        """Return active jobs plus the retry array size."""
        return self.active_jobs + self.retry_tasks

    @property
    def allowed(self) -> bool:
        """Return True when submitting stays under the QOS cap."""
        return self.projected_total <= (self.max_submit_jobs - self.safety_margin)


def replica_was_attempted(run_directory: Path, replica: int) -> bool:
    """Return True when a replica has any on-disk attempt artifacts."""
    if replica_h5_path(run_directory, replica).is_file():
        return True
    if protocol_marker_path(run_directory, replica).is_file():
        return True
    if replica_gsd_path(run_directory, replica).is_file():
        return True
    return any(
        fkt_file_path(run_directory, replica, reference_index).is_file()
        for reference_index in range(DEFAULT_FKT_MAX_REFERENCES)
    )


def collect_retry_candidates(
    output_base: Path,
    lambdas: Sequence[float],
    *,
    target_valid: int,
    n_replicas: int,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    runtime_ps: float = DEFAULT_RUNTIME_PS,
    reference_interval_ps: float = DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    max_references: int = DEFAULT_FKT_MAX_REFERENCES,
    min_bytes: int = COMPLETE_MIN_BYTES,
) -> tuple[RetryCandidate, ...]:
    """Collect incomplete replicas that were attempted but not yet valid.

    Parameters
    ----------
    output_base
        Campaign output root.
    lambdas
        Coupling values to scan.
    target_valid
        Maximum replica ID domain for the production campaign.
    n_replicas
        Upper scan bound passed to ``scan_lambda_run``.

    Returns
    -------
    tuple[RetryCandidate, ...]
        Deterministically sorted retry rows.
    """
    if target_valid <= 0:
        raise ValueError("target_valid must be positive")
    if n_replicas <= 0:
        raise ValueError("n_replicas must be positive")

    candidates: list[RetryCandidate] = []
    for lam in lambdas:
        require_step_protocol = lam != 0.0
        scan = scan_lambda_run(
            output_base=output_base,
            lam=lam,
            n_replicas=n_replicas,
            switch_time_ps=switch_time_ps,
            min_bytes=min_bytes,
            require_fkt=True,
            runtime_ps=runtime_ps,
            reference_interval_ps=reference_interval_ps,
            max_references=max_references,
            require_step_protocol=require_step_protocol,
            expected_lam=lam if require_step_protocol else None,
            expected_switch_time_ps=switch_time_ps if require_step_protocol else None,
        )
        directory = scan.run_dir
        retry_replicas = set(scan.partial_replicas)
        for replica in scan.missing_replicas:
            if replica >= target_valid:
                continue
            if replica_was_attempted(directory, replica):
                retry_replicas.add(replica)
        for replica in sorted(retry_replicas):
            if is_complete_replica(
                directory,
                replica,
                min_bytes=min_bytes,
                runtime_ps=runtime_ps,
                reference_interval_ps=reference_interval_ps,
                max_references=max_references,
                require_step_protocol=require_step_protocol,
                expected_lam=lam if require_step_protocol else None,
                expected_switch_time_ps=switch_time_ps if require_step_protocol else None,
            ):
                continue
            candidates.append(
                RetryCandidate(
                    lambda_tag=lambda_to_tag(lam),
                    lam=lam,
                    replica=replica,
                )
            )

    candidates.sort(key=lambda item: (item.lam, item.replica))
    return tuple(candidates)


def build_retry_manifest_lines(candidates: Sequence[RetryCandidate]) -> list[str]:
    """Build one packed-task row per retry candidate."""
    return [
        "\t".join([candidate.lambda_tag, f"{candidate.lam:g}", str(candidate.replica), ""])
        for candidate in candidates
    ]


def write_retry_manifest(
    path: Path,
    candidates: Sequence[RetryCandidate],
) -> None:
    """Write a retry manifest with one replica per array task."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = build_retry_manifest_lines(candidates)
    path.write_text("".join(f"{line}\n" for line in lines), encoding="utf-8")


def parse_manifest_row(row: str) -> tuple[str, float, int]:
    """Decode one manifest row to ``(lambda_tag, lam, replica)``."""
    parts = row.split("\t")
    if len(parts) < 3:
        raise ValueError(f"invalid manifest row: {row!r}")
    return parts[0], float(parts[1]), int(parts[2])


def manifest_path_for_cycle(plan_dir: Path, cycle_id: str) -> Path:
    """Return the manifest path for one watcher cycle."""
    return plan_dir / f"manifest_retry_{cycle_id}.tsv"


def sbatch_path_for_cycle(plan_dir: Path, cycle_id: str) -> Path:
    """Return the sbatch path for one watcher cycle."""
    return plan_dir / f"job_retry_{cycle_id}.sbatch"


def render_retry_sbatch(
    *,
    manifest_path: Path,
    n_tasks: int,
    array_concurrency: int,
    job_name: str,
    partition: str,
    account: str,
    output_base: Path,
    examples_dir: Path,
    input_gsd: Path,
    log_dir: Path,
    walltime: str,
) -> str:
    """Render one retry sbatch script for a fixed-size array."""
    if n_tasks <= 0:
        raise ValueError("n_tasks must be positive")
    if array_concurrency <= 0:
        raise ValueError("array_concurrency must be positive")

    last_index = n_tasks - 1
    repo_root = Path(__file__).resolve().parents[2]
    runner = repo_root / "examples" / "slurm" / "run_packed_aging_task.py"
    conda_sh = Path("/scratch/mh7373/miniforge3/etc/profile.d/conda.sh")
    return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --account={account}
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time={walltime}
#SBATCH --array=0-{last_index}%{array_concurrency}
#SBATCH --output={log_dir}/{job_name}_%A_%a.out
#SBATCH --error={log_dir}/{job_name}_%A_%a.err
#SBATCH --mail-type=NONE

set -euo pipefail
if [[ "$(hostname -s)" == "ga015" ]]; then
  echo "Refusing ga015: CUDA init failures observed on this node" >&2
  exit 1
fi
source {conda_sh}
conda activate hoomd
export OMP_NUM_THREADS=4

python {runner} \\
    --manifest {manifest_path} \\
    --task-index "${{SLURM_ARRAY_TASK_ID}}" \\
    --output-base {output_base} \\
    --examples-dir {examples_dir} \\
    --input-gsd {input_gsd} \\
    --log-dir {log_dir}
"""


def count_active_jobs(*, user: str, runner: CommandRunner | None = None) -> int:
    """Count pending and running SLURM jobs for one user."""
    run = runner or subprocess.run
    result = run(
        ["squeue", "-u", user, "-h", "-t", "PENDING,RUNNING"],
        check=True,
        capture_output=True,
        text=True,
    )
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    return len(lines)


def retry_job_in_flight(
    *,
    user: str,
    job_name_prefix: str,
    runner: CommandRunner | None = None,
) -> bool:
    """Return True when a retry array is already pending or running."""
    run = runner or subprocess.run
    result = run(
        [
            "squeue",
            "-u",
            user,
            "-h",
            "-t",
            "PENDING,RUNNING",
            "-o",
            "%T %j",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        _state, name = parts
        if name == job_name_prefix or name.startswith(f"{job_name_prefix}-"):
            return True
    return False


class SubmitResult:
    """Outcome of one retry submission attempt."""

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


def submit_retry_array(*, sbatch_file: Path) -> SubmitResult:
    """Submit one retry array job."""
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
    if "QOSMaxSubmitJobPerUserLimit" in combined or "QOSMaxGRESPerUser" in combined:
        return SubmitResult(status="quota_blocked", message=combined)
    return SubmitResult(status="error", message=combined or f"exit {result.returncode}")


def log_message(log_file: Path, message: str) -> None:
    """Append one timestamped line to the watcher log."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with log_file.open("a", encoding="utf-8") as handle:
        handle.write(f"{timestamp} {message}\n")


def append_submission_record(
    *,
    submitted_log: Path,
    job_id: str,
    n_tasks: int,
    walltime: str,
    cycle_id: str,
) -> None:
    """Append one retry submission record to the campaign log."""
    submitted_log.parent.mkdir(parents=True, exist_ok=True)
    with submitted_log.open("a", encoding="utf-8") as handle:
        handle.write(
            f"{job_id}:packed-r1-defer24-retry:{n_tasks}:walltime={walltime}:cycle={cycle_id}\n"
        )


def write_state(
    *,
    state_file: Path,
    cycle_id: str,
    job_id: str,
    n_tasks: int,
    manifest_path: Path,
    sbatch_path: Path,
) -> None:
    """Persist metadata for the latest retry submission."""
    payload = {
        "cycle_id": cycle_id,
        "job_id": job_id,
        "n_tasks": n_tasks,
        "manifest_path": str(manifest_path),
        "sbatch_path": str(sbatch_path),
        "submitted_at": datetime.now(timezone.utc).isoformat(),
    }
    state_file.parent.mkdir(parents=True, exist_ok=True)
    state_file.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run_retry_cycle(args: argparse.Namespace) -> str:
    """Scan once and submit a retry array when appropriate.

    Returns
    -------
    str
        One of ``submitted``, ``nothing_to_retry``, ``retry_in_flight``,
        ``quota_blocked``, or ``error``.
    """
    user = args.user or os.environ.get("USER", "")
    if not user:
        raise ValueError("cannot determine SLURM user")

    candidates = collect_retry_candidates(
        args.output_base,
        args.lambdas,
        target_valid=args.target_valid,
        n_replicas=args.scan_replicas,
        switch_time_ps=args.switch_time_ps,
        runtime_ps=args.runtime_ps,
        reference_interval_ps=args.reference_interval_ps,
        max_references=args.max_references,
        min_bytes=args.min_bytes,
    )
    log_message(args.log_file, f"retry_candidates={len(candidates)}")
    if not candidates:
        return "nothing_to_retry"

    if retry_job_in_flight(user=user, job_name_prefix=args.job_name):
        log_message(args.log_file, "retry array already pending or running")
        return "retry_in_flight"

    active_jobs = count_active_jobs(user=user)
    gate = RetrySubmitGate(
        active_jobs=active_jobs,
        retry_tasks=len(candidates),
        max_submit_jobs=args.max_submit_jobs,
        safety_margin=args.safety_margin,
    )
    log_message(
        args.log_file,
        (
            f"active_jobs={active_jobs} retry_tasks={len(candidates)} "
            f"projected_total={gate.projected_total} qos_allowed={gate.allowed}"
        ),
    )
    if not gate.allowed:
        log_message(args.log_file, "retry blocked by projected QOS submit cap")
        return "quota_blocked"

    cycle_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest_path = manifest_path_for_cycle(args.plan_dir, cycle_id)
    sbatch_path = sbatch_path_for_cycle(args.plan_dir, cycle_id)
    write_retry_manifest(manifest_path, candidates)
    sbatch_path.write_text(
        render_retry_sbatch(
            manifest_path=manifest_path,
            n_tasks=len(candidates),
            array_concurrency=args.array_concurrency,
            job_name=args.job_name,
            partition=args.partition,
            account=args.account,
            output_base=args.output_base,
            examples_dir=args.examples_dir,
            input_gsd=args.input_gsd,
            log_dir=args.log_dir,
            walltime=args.walltime,
        ),
        encoding="utf-8",
    )
    sbatch_path.chmod(0o750)

    if args.dry_run:
        log_message(
            args.log_file,
            f"dry-run wrote manifest={manifest_path} sbatch={sbatch_path}",
        )
        return "dry_run"

    outcome = submit_retry_array(sbatch_file=sbatch_path)
    if outcome.status == "submitted" and outcome.job_id is not None:
        append_submission_record(
            submitted_log=args.submitted_log,
            job_id=outcome.job_id,
            n_tasks=len(candidates),
            walltime=args.walltime,
            cycle_id=cycle_id,
        )
        write_state(
            state_file=args.state_file,
            cycle_id=cycle_id,
            job_id=outcome.job_id,
            n_tasks=len(candidates),
            manifest_path=manifest_path,
            sbatch_path=sbatch_path,
        )
        log_message(
            args.log_file,
            f"submitted retry job {outcome.job_id} tasks={len(candidates)}",
        )
        return "submitted"
    if outcome.status == "quota_blocked":
        log_message(args.log_file, "sbatch blocked by QOS quota")
        return "quota_blocked"

    log_message(args.log_file, f"submit failed: {outcome.message}")
    return "error"


def run_watcher(args: argparse.Namespace) -> int:
    """Poll hourly until stopped, submitting retry arrays as needed."""
    deadline = time.monotonic() + args.max_runtime_s
    while True:
        status = run_retry_cycle(args)
        if status == "error" and args.fail_on_error:
            return 1
        if time.monotonic() >= deadline:
            log_message(args.log_file, "max runtime reached; watcher exiting")
            return 0
        time.sleep(args.poll_interval_s)


def _parse_lambdas(raw: str) -> tuple[float, ...]:
    values = tuple(float(item.strip()) for item in raw.split(",") if item.strip())
    if not values:
        raise ValueError("at least one lambda value is required")
    return values


def build_parser() -> argparse.ArgumentParser:
    """Build the retry watcher CLI."""
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[1]
    default_plan_dir = (
        script_dir / "generated" / "aging_packed_20260719T051155Z-r1-defer24"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=default_plan_dir)
    parser.add_argument(
        "--output-base",
        type=Path,
        default=repo_root / "aging_weak_lambda",
    )
    parser.add_argument(
        "--examples-dir",
        type=Path,
        default=repo_root / "examples",
    )
    parser.add_argument(
        "--input-gsd",
        type=Path,
        default=repo_root
        / "aging_weak_lambda_ic_library"
        / "init-500-from-ic-library.gsd",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=script_dir / "logs" / "packed",
    )
    parser.add_argument(
        "--submitted-log",
        type=Path,
        default=script_dir / "submitted_aging_campaign.txt",
    )
    parser.add_argument(
        "--state-file",
        type=Path,
        default=default_plan_dir / "retry_watcher_state.json",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=script_dir / "logs" / "watch_retry_failed.log",
    )
    parser.add_argument("--user", type=str, default=None)
    parser.add_argument("--job-name", type=str, default=DEFAULT_JOB_NAME)
    parser.add_argument("--partition", type=str, default=DEFAULT_PARTITION)
    parser.add_argument("--account", type=str, default=DEFAULT_ACCOUNT)
    parser.add_argument("--walltime", type=str, default=DEFAULT_WALLTIME)
    parser.add_argument(
        "--lambdas",
        type=_parse_lambdas,
        default=",".join(str(value) for value in DEFAULT_LAMBDAS),
    )
    parser.add_argument("--target-valid", type=int, default=DEFAULT_TARGET_VALID)
    parser.add_argument("--scan-replicas", type=int, default=DEFAULT_TARGET_VALID)
    parser.add_argument(
        "--switch-time-ps",
        type=float,
        default=DEFAULT_SWITCH_TIME_PS,
    )
    parser.add_argument("--runtime-ps", type=float, default=DEFAULT_RUNTIME_PS)
    parser.add_argument(
        "--reference-interval-ps",
        type=float,
        default=DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    )
    parser.add_argument(
        "--max-references",
        type=int,
        default=DEFAULT_FKT_MAX_REFERENCES,
    )
    parser.add_argument("--min-bytes", type=int, default=COMPLETE_MIN_BYTES)
    parser.add_argument(
        "--max-submit-jobs",
        type=int,
        default=DEFAULT_MAX_SUBMIT_JOBS,
    )
    parser.add_argument(
        "--safety-margin",
        type=int,
        default=DEFAULT_SAFETY_MARGIN,
    )
    parser.add_argument(
        "--array-concurrency",
        type=int,
        default=DEFAULT_ARRAY_CONCURRENCY,
    )
    parser.add_argument(
        "--poll-interval-s",
        type=int,
        default=DEFAULT_POLL_INTERVAL_S,
    )
    parser.add_argument(
        "--max-runtime-s",
        type=int,
        default=14 * 24 * 3600,
        help="Stop after this many seconds.",
    )
    parser.add_argument(
        "--once",
        action="store_true",
        help="Run one scan/submit cycle and exit.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write manifest and sbatch but do not call sbatch.",
    )
    parser.add_argument(
        "--fail-on-error",
        action="store_true",
        help="Exit with status 1 when a cycle returns error.",
    )
    return parser


def main() -> int:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()
    if args.target_valid <= 0:
        parser.error("--target-valid must be positive")
    if args.scan_replicas <= 0:
        parser.error("--scan-replicas must be positive")
    if args.poll_interval_s <= 0:
        parser.error("--poll-interval-s must be positive")
    if args.array_concurrency <= 0:
        parser.error("--array-concurrency must be positive")
    if args.once:
        status = run_retry_cycle(args)
        if status == "error" and args.fail_on_error:
            return 1
        return 0
    return run_watcher(args)


if __name__ == "__main__":
    sys.exit(main())
