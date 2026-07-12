"""Run one manifest row containing up to two concurrent aging replicas."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
import fcntl
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from typing import Iterator, Mapping, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from examples.slurm.aging_campaign_status import (
        cleanup_replica_artifacts,
        is_complete_replica,
        lambda_to_tag,
        run_dir,
    )
else:
    from .aging_campaign_status import (
        cleanup_replica_artifacts,
        is_complete_replica,
        lambda_to_tag,
        run_dir,
    )

DEFAULT_OUTPUT_BASE = Path(
    "/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda"
)
CAMPAIGN_LAMBDAS = (0.0, 0.01, 0.016667, 0.023333, 0.03)
MAX_REPLICA_ID = 999


@dataclass(frozen=True)
class ManifestTask:
    """One same-coupling row from the packed campaign manifest."""

    lambda_tag: str
    lam: float
    replicas: tuple[int, ...]

    def __post_init__(self) -> None:
        """Reject unsafe task rows before any outputs are touched."""
        if not self.lambda_tag:
            raise ValueError("lambda tag cannot be empty")
        if not math.isfinite(self.lam) or self.lam < 0.0:
            raise ValueError("lambda value must be finite and nonnegative")
        if not 1 <= len(self.replicas) <= 2:
            raise ValueError("a manifest task requires one or two replicas")
        if len(set(self.replicas)) != len(self.replicas):
            raise ValueError("a manifest task cannot repeat a replica")
        if any(not 0 <= replica <= MAX_REPLICA_ID for replica in self.replicas):
            raise ValueError("replica IDs must be in the campaign domain 0..999")
        if self.lam not in CAMPAIGN_LAMBDAS:
            raise ValueError("lambda value is not in the aging campaign")
        if self.lambda_tag != lambda_to_tag(self.lam):
            raise ValueError("lambda tag does not match the lambda value")


@dataclass(frozen=True)
class ReplicaProcessSpec:
    """Subprocess command, working directory, environment, and log paths."""

    replica: int
    command: tuple[str, ...]
    cwd: Path
    stdout_path: Path
    stderr_path: Path
    environment: Mapping[str, str]


def read_manifest_task(path: Path, task_index: int) -> ManifestTask:
    """Read and validate one zero-based manifest row.

    Parameters
    ----------
    path
        Tab-separated packed campaign manifest.
    task_index
        Zero-based SLURM array index.

    Returns
    -------
    ManifestTask
        Validated coupling and one or two replica IDs.
    """
    if task_index < 0:
        raise IndexError("task index must be nonnegative")
    try:
        rows = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise ValueError(f"cannot read manifest {path}: {error}") from error
    if task_index >= len(rows):
        raise IndexError(
            f"task index {task_index} is outside manifest size {len(rows)}"
        )

    fields = rows[task_index].split("\t")
    if len(fields) != 4:
        raise ValueError("manifest rows must contain exactly four fields")
    lambda_tag, lambda_text, first_text, second_text = fields
    if not first_text:
        raise ValueError("manifest row is missing its first replica")
    try:
        lam = float(lambda_text)
        replicas = [int(first_text)]
        if second_text:
            replicas.append(int(second_text))
    except ValueError as error:
        raise ValueError("manifest row contains a nonnumeric value") from error
    return ManifestTask(lambda_tag, lam, tuple(replicas))


def child_environment(
    parent_environment: Mapping[str, str],
) -> dict[str, str]:
    """Return an environment that forces explicit replica selection.

    ``05_advanced_run.py`` otherwise replaces ``--replicas`` with the enclosing
    SLURM array index. Removing array identity lets both child simulations use
    their distinct manifest replica IDs while retaining the allocation job ID.

    Parameters
    ----------
    parent_environment
        Environment inherited by the packed task.

    Returns
    -------
    dict[str, str]
        Copy without SLURM array identity variables.
    """
    environment = dict(parent_environment)
    environment.pop("SLURM_ARRAY_TASK_ID", None)
    environment.pop("SLURM_ARRAY_JOB_ID", None)
    return environment


@contextmanager
def replica_locks(
    run_directory: Path,
    replicas: Sequence[int],
) -> Iterator[None]:
    """Hold exclusive per-replica locks in deterministic ID order.

    Parameters
    ----------
    run_directory
        Coupling output directory containing lock files.
    replicas
        Replica IDs protected from revalidation through post-run validation.

    Yields
    ------
    None
        Control while every requested replica lock is held.

    Side Effects
    ------------
    Creates persistent zero-length lock files. Kernel locks are released when
    the context exits, including on cancellation or failure.
    """
    ordered = tuple(sorted(replicas))
    if not ordered or len(set(ordered)) != len(ordered):
        raise ValueError("replica locks require distinct replica IDs")
    if any(not 0 <= replica <= MAX_REPLICA_ID for replica in ordered):
        raise ValueError("replica IDs must be in the campaign domain 0..999")

    run_directory.mkdir(parents=True, exist_ok=True)
    handles = []
    try:
        for replica in ordered:
            handle = (run_directory / f".replica_{replica}.lock").open("a+")
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            handles.append(handle)
        yield
    finally:
        for handle in reversed(handles):
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
            handle.close()


def _request_process_termination(
    processes: Sequence[subprocess.Popen[bytes]],
) -> None:
    """Request termination without blocking inside a signal handler."""
    for process in processes:
        if process.poll() is None:
            process.terminate()


def _terminate_processes(processes: Sequence[subprocess.Popen[bytes]]) -> None:
    """Terminate all active child processes and wait briefly for exit."""
    _request_process_termination(processes)
    for process in processes:
        if process.poll() is None:
            try:
                process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()


def run_concurrent_processes(
    specifications: Sequence[ReplicaProcessSpec],
    *,
    termination_timeout_seconds: float = 30.0,
) -> int:
    """Start up to two replica processes concurrently and aggregate status.

    Parameters
    ----------
    specifications
        One or two distinct replica process definitions.
    termination_timeout_seconds
        Grace period after SIGINT or SIGTERM before active children are killed.

    Returns
    -------
    int
        Zero only when every child exits successfully. A forwarded termination
        signal returns the conventional ``128 + signal`` status.

    Side Effects
    ------------
    Creates separate standard-output and standard-error logs for every replica.
    """
    if not 1 <= len(specifications) <= 2:
        raise ValueError("one or two process specifications are required")
    if (
        not math.isfinite(termination_timeout_seconds)
        or termination_timeout_seconds <= 0.0
    ):
        raise ValueError("termination timeout must be finite and positive")
    replicas = [specification.replica for specification in specifications]
    if len(set(replicas)) != len(replicas):
        raise ValueError("process specifications repeat a replica")
    output_paths = [
        path
        for specification in specifications
        for path in (specification.stdout_path, specification.stderr_path)
    ]
    if len(set(output_paths)) != len(output_paths):
        raise ValueError("each child requires distinct output and error logs")

    processes: list[subprocess.Popen[bytes]] = []
    handles: list[object] = []
    received_signal: list[int | None] = [None]
    previous_handlers: dict[int, signal.Handlers] = {}

    def handle_signal(signum: int, _frame: object) -> None:
        received_signal[0] = signum
        _request_process_termination(processes)

    try:
        for signum in (signal.SIGINT, signal.SIGTERM):
            previous_handlers[signum] = signal.getsignal(signum)
            signal.signal(signum, handle_signal)

        for specification in specifications:
            if received_signal[0] is not None:
                break
            specification.stdout_path.parent.mkdir(parents=True, exist_ok=True)
            specification.stderr_path.parent.mkdir(parents=True, exist_ok=True)
            stdout_handle = specification.stdout_path.open("wb")
            stderr_handle = specification.stderr_path.open("wb")
            handles.extend((stdout_handle, stderr_handle))
            process = subprocess.Popen(
                specification.command,
                cwd=specification.cwd,
                env=dict(specification.environment),
                stdout=stdout_handle,
                stderr=stderr_handle,
            )
            processes.append(process)
            if received_signal[0] is not None:
                _request_process_termination((process,))
                break

        return_codes: list[int | None] = [None] * len(processes)
        termination_deadline: float | None = None
        while any(return_code is None for return_code in return_codes):
            for index, process in enumerate(processes):
                if return_codes[index] is None:
                    return_codes[index] = process.poll()
            if received_signal[0] is not None:
                if termination_deadline is None:
                    termination_deadline = (
                        time.monotonic() + termination_timeout_seconds
                    )
                if time.monotonic() >= termination_deadline:
                    for process in processes:
                        if process.poll() is None:
                            process.kill()
                    for index, process in enumerate(processes):
                        if return_codes[index] is None:
                            return_codes[index] = process.wait()
            if any(return_code is None for return_code in return_codes):
                time.sleep(0.05)
    except BaseException:
        _terminate_processes(processes)
        raise
    finally:
        for handle in handles:
            handle.close()
        for signum, previous_handler in previous_handlers.items():
            signal.signal(signum, previous_handler)

    if received_signal[0] is not None:
        return 128 + received_signal[0]
    return 0 if all(return_code == 0 for return_code in return_codes) else 1


def build_simulation_command(
    *,
    python_executable: Path,
    simulation_script: Path,
    replica: int,
    lam: float,
    input_gsd: Path,
) -> tuple[str, ...]:
    """Build the fixed aging-protocol command for one replica.

    Parameters use the established campaign units: temperature in kelvin,
    frequency in inverse centimeters, runtime in picoseconds, and timestep in
    femtoseconds.
    """
    return (
        str(python_executable),
        str(simulation_script),
        "--molecular-bath",
        "bussi",
        "--cavity-bath",
        "langevin",
        "--coupling",
        f"{lam:g}",
        "--temperature",
        "100.0",
        "--frequency",
        "1560.0",
        "--runtime",
        "2500.0",
        "--switch-time",
        "200.0",
        "--input-gsd",
        str(input_gsd),
        "--frame",
        "-1",
        "--device",
        "GPU",
        "--gpu-id",
        "0",
        "--fixed-timestep",
        "--timestep",
        "1.0",
        "--enable-energy-tracker",
        "--energy-output-period-ps",
        "1.0",
        "--enable-fkt",
        "--fkt-kmag",
        "1.0085",
        "--fkt-wavevectors",
        "50",
        "--fkt-ref-interval",
        "200.0",
        "--fkt-max-refs",
        "13",
        "--fkt-output-period-ps",
        "1.0",
        "--gsd-output-period-ps",
        "999999",
        "--console-output-period-ps",
        "100.0",
        "--replicas",
        str(replica),
        "--seed",
        str(replica + 1),
    )


def _parse_args() -> argparse.Namespace:
    """Parse packed-task command-line arguments."""
    examples_directory = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT_BASE)
    parser.add_argument("--examples-dir", type=Path, default=examples_directory)
    parser.add_argument(
        "--input-gsd",
        type=Path,
        default=examples_directory / "init-0.gsd",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "logs" / "packed",
    )
    parser.add_argument(
        "--python-executable",
        type=Path,
        default=Path(sys.executable),
    )
    return parser.parse_args()


def _task_index(argument_value: int | None) -> int:
    """Resolve a task index from CLI input or the SLURM environment."""
    if argument_value is not None:
        return argument_value
    environment_value = os.environ.get("SLURM_ARRAY_TASK_ID")
    if environment_value is None:
        raise ValueError(
            "--task-index is required outside a SLURM array allocation"
        )
    return int(environment_value)


def main() -> int:
    """Validate, clean, and concurrently execute one packed manifest row."""
    args = _parse_args()
    task_index = _task_index(args.task_index)
    task = read_manifest_task(args.manifest, task_index)
    coupling_directory = args.output_base / f"lambda{task.lambda_tag}"
    coupling_directory.mkdir(parents=True, exist_ok=True)
    run_directory = run_dir(args.output_base, task.lam)
    simulation_script = args.examples_dir / "05_advanced_run.py"
    job_id = os.environ.get("SLURM_JOB_ID", "local")
    array_id = os.environ.get("SLURM_ARRAY_TASK_ID", str(task_index))
    environment = child_environment(os.environ)
    specifications: list[ReplicaProcessSpec] = []

    with replica_locks(run_directory, task.replicas):
        for replica in task.replicas:
            if is_complete_replica(run_directory, replica):
                print(
                    f"Skipping scientifically complete lambda={task.lam:g} "
                    f"replica={replica}",
                    flush=True,
                )
                continue
            deleted = cleanup_replica_artifacts(run_directory, replica)
            if deleted:
                print(
                    f"Cleaned {len(deleted)} invalid artifacts for "
                    f"lambda={task.lam:g} replica={replica}",
                    flush=True,
                )
            log_stem = (
                f"cav-aging-lam{task.lambda_tag}_{job_id}_{array_id}"
                f"_replica_{replica}"
            )
            specifications.append(
                ReplicaProcessSpec(
                    replica=replica,
                    command=build_simulation_command(
                        python_executable=args.python_executable,
                        simulation_script=simulation_script,
                        replica=replica,
                        lam=task.lam,
                        input_gsd=args.input_gsd,
                    ),
                    cwd=coupling_directory,
                    stdout_path=args.log_dir / f"{log_stem}.out",
                    stderr_path=args.log_dir / f"{log_stem}.err",
                    environment=environment,
                )
            )

        if not specifications:
            print("All manifest replicas are already complete", flush=True)
            return 0

        result = run_concurrent_processes(specifications)
        incomplete = [
            specification.replica
            for specification in specifications
            if not is_complete_replica(run_directory, specification.replica)
        ]
        if incomplete:
            print(
                "Robust post-run validation failed for replicas "
                f"{incomplete}",
                file=sys.stderr,
                flush=True,
            )
        if result >= 128:
            return result
        return 0 if result == 0 and not incomplete else 1


if __name__ == "__main__":
    raise SystemExit(main())
