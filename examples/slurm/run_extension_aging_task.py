"""Run one extension manifest row for checkpointed 1600→2000 ps aging."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from examples.slurm.aging_campaign_config import (
        CAMPAIGN_EXTENSION_RUNTIME_PS,
        CAMPAIGN_EXTENSION_WALLTIME,
        CAMPAIGN_FKT_KMAG,
        CAMPAIGN_FKT_MAX_REFS_TARGET,
        CAMPAIGN_FKT_REFERENCE_INTERVAL_PS,
        CAMPAIGN_HDF5_PERIOD_PS,
        CAMPAIGN_PRIMARY_RUNTIME_PS,
        CAMPAIGN_SWITCH_TIME_PS,
        CAMPAIGN_TARGET_RUNTIME_PS,
        CHECKPOINT_GSD_TEMPLATE,
        FKT_STATE_TEMPLATE,
    )
    from examples.slurm.aging_campaign_status import (
        COMPLETE_MIN_BYTES,
        is_complete_replica,
        needs_extension_replica,
        replica_checkpoint_path,
        replica_fkt_state_path,
        run_dir,
        write_protocol_marker,
    )
    from examples.slurm.run_packed_aging_task import (
        CAMPAIGN_ELECTROSTATICS,
        CAMPAIGN_COULOMB_RCUT,
        CAMPAIGN_EPS_RF,
        CAMPAIGN_LAMBDAS,
        child_environment,
        collect_simulator_provenance,
        read_manifest_task,
        replica_locks,
        run_concurrent_processes,
        ReplicaProcessSpec,
    )
else:
    from .aging_campaign_config import (
        CAMPAIGN_EXTENSION_RUNTIME_PS,
        CAMPAIGN_EXTENSION_WALLTIME,
        CAMPAIGN_FKT_KMAG,
        CAMPAIGN_FKT_MAX_REFS_TARGET,
        CAMPAIGN_FKT_REFERENCE_INTERVAL_PS,
        CAMPAIGN_HDF5_PERIOD_PS,
        CAMPAIGN_PRIMARY_RUNTIME_PS,
        CAMPAIGN_SWITCH_TIME_PS,
        CAMPAIGN_TARGET_RUNTIME_PS,
        CHECKPOINT_GSD_TEMPLATE,
        FKT_STATE_TEMPLATE,
    )
    from .aging_campaign_status import (
        COMPLETE_MIN_BYTES,
        is_complete_replica,
        needs_extension_replica,
        replica_checkpoint_path,
        replica_fkt_state_path,
        run_dir,
        write_protocol_marker,
    )
    from .run_packed_aging_task import (
        CAMPAIGN_ELECTROSTATICS,
        CAMPAIGN_COULOMB_RCUT,
        CAMPAIGN_EPS_RF,
        CAMPAIGN_LAMBDAS,
        child_environment,
        collect_simulator_provenance,
        read_manifest_task,
        replica_locks,
        run_concurrent_processes,
        ReplicaProcessSpec,
    )

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_BASE = REPOSITORY_ROOT / "aging_weak_lambda"


def build_extension_command(
    *,
    python_executable: Path,
    simulation_script: Path,
    replica: int,
    lam: float,
    checkpoint_gsd: Path,
    fkt_state_path: Path,
) -> tuple[str, ...]:
    """Build the checkpoint-resume command for one extension replica."""
    return (
        str(python_executable),
        str(simulation_script),
        "--molecular-bath",
        "bussi",
        "--cavity-bath",
        "langevin",
        "--lambda-coupling",
        f"{lam:g}",
        "--coupling-type",
        "step",
        "--temperature",
        "100.0",
        "--frequency",
        "1560.0",
        "--runtime",
        f"{CAMPAIGN_EXTENSION_RUNTIME_PS:g}",
        "--time-offset-ps",
        f"{CAMPAIGN_PRIMARY_RUNTIME_PS:g}",
        "--switch-time",
        f"{CAMPAIGN_SWITCH_TIME_PS:g}",
        "--input-gsd",
        str(checkpoint_gsd),
        "--frame",
        "-1",
        "--no-restart-velocities",
        "--device",
        "GPU",
        "--gpu-id",
        "0",
        "--fixed-timestep",
        "--timestep",
        "1.0",
        "--enable-energy-tracker",
        "--energy-output-period-ps",
        f"{CAMPAIGN_HDF5_PERIOD_PS:g}",
        "--hdf5-output-period-ps",
        f"{CAMPAIGN_HDF5_PERIOD_PS:g}",
        "--resume-hdf5-append",
        "--enable-fkt",
        "--fkt-kmag",
        f"{CAMPAIGN_FKT_KMAG:g}",
        "--fkt-wavevectors",
        "50",
        "--fkt-ref-interval",
        f"{CAMPAIGN_FKT_REFERENCE_INTERVAL_PS:g}",
        "--fkt-max-refs",
        str(CAMPAIGN_FKT_MAX_REFS_TARGET),
        "--fkt-output-period-ps",
        f"{CAMPAIGN_HDF5_PERIOD_PS:g}",
        "--resume-fkt-state-file",
        str(fkt_state_path),
        "--electrostatics",
        CAMPAIGN_ELECTROSTATICS,
        "--eps-rf",
        f"{CAMPAIGN_EPS_RF:g}",
        "--coulomb-rcut",
        f"{CAMPAIGN_COULOMB_RCUT:g}",
        "--disable-gsd",
        "--console-output-period-ps",
        "100.0",
        "--replicas",
        str(replica),
        "--seed",
        str(replica + 1),
    )


def _parse_args() -> argparse.Namespace:
    examples_directory = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT_BASE)
    parser.add_argument("--examples-dir", type=Path, default=examples_directory)
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "logs" / "extension",
    )
    parser.add_argument(
        "--python-executable",
        type=Path,
        default=Path(sys.executable),
    )
    return parser.parse_args()


def _task_index(argument_value: int | None) -> int:
    if argument_value is not None:
        return argument_value
    environment_value = os.environ.get("SLURM_ARRAY_TASK_ID")
    if environment_value is None:
        raise ValueError(
            "--task-index is required outside a SLURM array allocation"
        )
    return int(environment_value)


def main() -> int:
    args = _parse_args()
    task_index = _task_index(args.task_index)
    task = read_manifest_task(args.manifest, task_index)
    if task.lam not in CAMPAIGN_LAMBDAS:
        raise ValueError("manifest row uses a non-campaign coupling")
    coupling_directory = args.output_base / f"lambda{task.lambda_tag}"
    coupling_directory.mkdir(parents=True, exist_ok=True)
    run_directory = run_dir(args.output_base, task.lam)
    simulation_script = args.examples_dir / "05_advanced_run.py"
    job_id = os.environ.get("SLURM_JOB_ID", "local")
    array_id = os.environ.get("SLURM_ARRAY_TASK_ID", str(task_index))
    environment = child_environment(os.environ)
    completion_options = {
        "runtime_ps": CAMPAIGN_TARGET_RUNTIME_PS,
        "max_references": CAMPAIGN_FKT_MAX_REFS_TARGET,
        "reference_interval_ps": CAMPAIGN_FKT_REFERENCE_INTERVAL_PS,
        "min_bytes": COMPLETE_MIN_BYTES,
        "require_step_protocol": task.lam != 0.0,
        "expected_lam": task.lam,
        "expected_switch_time_ps": CAMPAIGN_SWITCH_TIME_PS,
    }
    specifications: list[ReplicaProcessSpec] = []

    with replica_locks(run_directory, task.replicas):
        for replica in task.replicas:
            if is_complete_replica(run_directory, replica, **completion_options):
                print(
                    f"Skipping already-extended lambda={task.lam:g} replica={replica}",
                    flush=True,
                )
                continue
            if not needs_extension_replica(run_directory, replica, **completion_options):
                print(
                    f"Skipping non-extension lambda={task.lam:g} replica={replica}",
                    flush=True,
                )
                continue
            checkpoint = replica_checkpoint_path(run_directory, replica)
            fkt_state = replica_fkt_state_path(run_directory, replica)
            if not checkpoint.is_file():
                print(
                    f"Missing checkpoint for lambda={task.lam:g} replica={replica}",
                    file=sys.stderr,
                    flush=True,
                )
                continue
            log_stem = (
                f"cav-aging-ext-lam{task.lambda_tag}_{job_id}_{array_id}"
                f"_replica_{replica}"
            )
            specifications.append(
                ReplicaProcessSpec(
                    replica=replica,
                    command=build_extension_command(
                        python_executable=args.python_executable,
                        simulation_script=simulation_script,
                        replica=replica,
                        lam=task.lam,
                        checkpoint_gsd=checkpoint,
                        fkt_state_path=fkt_state if fkt_state.is_file() else fkt_state,
                    ),
                    cwd=coupling_directory,
                    stdout_path=args.log_dir / f"{log_stem}.out",
                    stderr_path=args.log_dir / f"{log_stem}.err",
                    environment=environment,
                )
            )

        if not specifications:
            print("No extension replicas to run in this manifest row", flush=True)
            return 0

        result = run_concurrent_processes(specifications)
        provenance = collect_simulator_provenance(
            python_executable=args.python_executable,
            simulation_script=simulation_script,
        )
        for specification in specifications:
            if is_complete_replica(
                run_directory,
                specification.replica,
                **completion_options,
            ):
                write_protocol_marker(
                    run_directory,
                    specification.replica,
                    lam=task.lam,
                    switch_time_ps=CAMPAIGN_SWITCH_TIME_PS,
                    seed=specification.replica + 1,
                    provenance=provenance,
                )
        if result >= 128:
            return result
        incomplete = [
            specification.replica
            for specification in specifications
            if not is_complete_replica(
                run_directory,
                specification.replica,
                **completion_options,
            )
        ]
        return 0 if result == 0 and not incomplete else 1


if __name__ == "__main__":
    raise SystemExit(main())
