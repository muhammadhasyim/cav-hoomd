"""Tests for failure-safe concurrent packed aging tasks."""

from __future__ import annotations

import fcntl
import os
from pathlib import Path
import subprocess
import sys
import time
from contextlib import contextmanager

import pytest

from examples.slurm import run_packed_aging_task as packed_runner
from examples.slurm.run_packed_aging_task import (
    ManifestTask,
    ReplicaProcessSpec,
    build_simulation_command,
    child_environment,
    collect_simulator_provenance,
    read_manifest_task,
    replica_locks,
    run_concurrent_processes,
)


def _process_spec(
    tmp_path: Path,
    replica: int,
    *,
    exit_code: int = 0,
) -> ReplicaProcessSpec:
    """Return a deterministic child process writing distinct output."""
    command = (
        sys.executable,
        "-c",
        (
            "import sys,time;"
            f"print('stdout-{replica}', flush=True);"
            f"print('stderr-{replica}', file=sys.stderr, flush=True);"
            "time.sleep(0.1);"
            f"raise SystemExit({exit_code})"
        ),
    )
    return ReplicaProcessSpec(
        replica=replica,
        command=command,
        cwd=tmp_path,
        stdout_path=tmp_path / f"replica-{replica}.out",
        stderr_path=tmp_path / f"replica-{replica}.err",
        environment=os.environ.copy(),
    )


def test_read_manifest_task_reads_pair(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("0p03\t0.03\t4\t5\n", encoding="utf-8")

    assert read_manifest_task(manifest, 0) == ManifestTask(
        lambda_tag="0p03",
        lam=0.03,
        replicas=(4, 5),
    )


def test_read_manifest_task_reads_singleton(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("0p03\t0.03\t4\t\n", encoding="utf-8")

    assert read_manifest_task(manifest, 0).replicas == (4,)


def test_read_manifest_task_rejects_out_of_range_index(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("0p03\t0.03\t4\t5\n", encoding="utf-8")

    with pytest.raises(IndexError, match="task index"):
        read_manifest_task(manifest, 1)


def test_read_manifest_task_rejects_out_of_domain_replica(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("0p03\t0.03\t1000\t\n", encoding="utf-8")

    with pytest.raises(ValueError, match="0..999"):
        read_manifest_task(manifest, 0)


def test_read_manifest_task_rejects_unknown_coupling(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("0p02\t0.02\t4\t\n", encoding="utf-8")

    with pytest.raises(ValueError, match="not in the aging campaign"):
        read_manifest_task(manifest, 0)


def test_run_concurrent_processes_succeeds_and_writes_distinct_logs(
    tmp_path: Path,
) -> None:
    specs = [_process_spec(tmp_path, 4), _process_spec(tmp_path, 5)]

    result = run_concurrent_processes(specs)

    assert result == 0
    assert specs[0].stdout_path.read_text(encoding="utf-8") == "stdout-4\n"
    assert specs[1].stdout_path.read_text(encoding="utf-8") == "stdout-5\n"
    assert specs[0].stderr_path.read_text(encoding="utf-8") == "stderr-4\n"
    assert specs[1].stderr_path.read_text(encoding="utf-8") == "stderr-5\n"


def test_run_concurrent_processes_returns_nonzero_if_one_child_fails(
    tmp_path: Path,
) -> None:
    successful = _process_spec(tmp_path, 4)
    failed = _process_spec(tmp_path, 5, exit_code=7)

    result = run_concurrent_processes([successful, failed])

    assert result != 0
    assert successful.stdout_path.exists()
    assert failed.stdout_path.exists()


def test_run_concurrent_processes_waits_for_successful_sibling(
    tmp_path: Path,
) -> None:
    successful = ReplicaProcessSpec(
        replica=4,
        command=(
            sys.executable,
            "-c",
            "import time;time.sleep(0.2);print('finished', flush=True)",
        ),
        cwd=tmp_path,
        stdout_path=tmp_path / "successful.out",
        stderr_path=tmp_path / "successful.err",
        environment=os.environ.copy(),
    )
    failed = _process_spec(tmp_path, 5, exit_code=7)

    result = run_concurrent_processes([successful, failed])

    assert result != 0
    assert successful.stdout_path.read_text(encoding="utf-8") == "finished\n"


def test_run_concurrent_processes_returns_nonzero_if_both_children_fail(
    tmp_path: Path,
) -> None:
    first = _process_spec(tmp_path, 4, exit_code=3)
    second = _process_spec(tmp_path, 5, exit_code=7)

    assert run_concurrent_processes([first, second]) != 0


def test_run_concurrent_processes_supports_singleton(tmp_path: Path) -> None:
    singleton = _process_spec(tmp_path, 4)

    assert run_concurrent_processes([singleton]) == 0
    assert singleton.stdout_path.exists()


def test_run_concurrent_processes_launches_children_concurrently(
    tmp_path: Path,
) -> None:
    first_ready = tmp_path / "first.ready"
    second_ready = tmp_path / "second.ready"

    def handshake_spec(
        replica: int,
        own_ready: Path,
        peer_ready: Path,
    ) -> ReplicaProcessSpec:
        code = (
            "import pathlib,time;"
            f"pathlib.Path({str(own_ready)!r}).write_text('ready');"
            "deadline=time.monotonic()+2;"
            f"peer=pathlib.Path({str(peer_ready)!r});"
            "\nwhile not peer.exists() and time.monotonic()<deadline:"
            "\n time.sleep(0.01)"
            "\nraise SystemExit(0 if peer.exists() else 9)"
        )
        return ReplicaProcessSpec(
            replica=replica,
            command=(sys.executable, "-c", code),
            cwd=tmp_path,
            stdout_path=tmp_path / f"{replica}.out",
            stderr_path=tmp_path / f"{replica}.err",
            environment=os.environ.copy(),
        )

    result = run_concurrent_processes(
        [
            handshake_spec(4, first_ready, second_ready),
            handshake_spec(5, second_ready, first_ready),
        ]
    )

    assert result == 0


def test_run_concurrent_processes_forwards_termination_signal(
    tmp_path: Path,
) -> None:
    repository = Path(__file__).parents[1]
    first_pid = tmp_path / "first.pid"
    second_pid = tmp_path / "second.pid"
    code = f"""
import os
from pathlib import Path
import sys
from examples.slurm.run_packed_aging_task import (
    ReplicaProcessSpec,
    run_concurrent_processes,
)
base = Path({str(tmp_path)!r})
specs = []
rows = (
    (4, Path({str(first_pid)!r})),
    (5, Path({str(second_pid)!r})),
)
for replica, pid_path in rows:
    child_code = (
        'import os,pathlib,time;'
        + 'pathlib.Path(' + repr(str(pid_path)) + ')'
        + '.write_text(str(os.getpid()));time.sleep(30)'
    )
    child = (
        sys.executable,
        '-c',
        child_code,
    )
    specs.append(
        ReplicaProcessSpec(
            replica,
            child,
            base,
            base / f'{{replica}}.out',
            base / f'{{replica}}.err',
            os.environ.copy(),
        )
    )
raise SystemExit(run_concurrent_processes(specs))
"""
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(repository)
    driver = subprocess.Popen(
        [sys.executable, "-c", code],
        cwd=repository,
        env=environment,
    )
    deadline = time.monotonic() + 5.0
    while (
        (not first_pid.exists() or not second_pid.exists())
        and time.monotonic() < deadline
    ):
        time.sleep(0.02)
    assert first_pid.exists() and second_pid.exists()

    driver.terminate()
    assert driver.wait(timeout=5) == 128 + 15

    for path in (first_pid, second_pid):
        pid = int(path.read_text(encoding="utf-8"))
        with pytest.raises(ProcessLookupError):
            os.kill(pid, 0)


def test_run_concurrent_processes_kills_child_ignoring_termination(
    tmp_path: Path,
) -> None:
    repository = Path(__file__).parents[1]
    child_pid = tmp_path / "ignoring.pid"
    code = f"""
import os
from pathlib import Path
import sys
from examples.slurm.run_packed_aging_task import (
    ReplicaProcessSpec,
    run_concurrent_processes,
)
base = Path({str(tmp_path)!r})
pid_path = Path({str(child_pid)!r})
child_code = (
    'import os,pathlib,signal,time;'
    + 'signal.signal(signal.SIGTERM, signal.SIG_IGN);'
    + 'pathlib.Path(' + repr(str(pid_path)) + ')'
    + '.write_text(str(os.getpid()));time.sleep(30)'
)
spec = ReplicaProcessSpec(
    4,
    (sys.executable, '-c', child_code),
    base,
    base / '4.out',
    base / '4.err',
    os.environ.copy(),
)
raise SystemExit(
    run_concurrent_processes([spec], termination_timeout_seconds=0.2)
)
"""
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(repository)
    driver = subprocess.Popen(
        [sys.executable, "-c", code],
        cwd=repository,
        env=environment,
    )
    deadline = time.monotonic() + 5.0
    while not child_pid.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert child_pid.exists()

    driver.terminate()
    assert driver.wait(timeout=5) == 128 + 15
    pid = int(child_pid.read_text(encoding="utf-8"))
    with pytest.raises(ProcessLookupError):
        os.kill(pid, 0)


def test_child_environment_removes_array_task_identity() -> None:
    parent = {
        "PATH": "/bin",
        "SLURM_JOB_ID": "123",
        "SLURM_ARRAY_JOB_ID": "120",
        "SLURM_ARRAY_TASK_ID": "7",
    }

    environment = child_environment(parent)

    assert environment["SLURM_JOB_ID"] == "123"
    assert "SLURM_ARRAY_TASK_ID" not in environment
    assert "SLURM_ARRAY_JOB_ID" not in environment


def test_build_simulation_command_sets_stable_seed_and_replica(
    tmp_path: Path,
) -> None:
    command = build_simulation_command(
        python_executable=Path(sys.executable),
        simulation_script=tmp_path / "05_advanced_run.py",
        replica=17,
        lam=0.03,
        input_gsd=tmp_path / "init-0.gsd",
    )

    assert command[command.index("--replicas") + 1] == "17"
    assert command[command.index("--seed") + 1] == "18"
    assert command[command.index("--gpu-id") + 1] == "0"
    assert command[command.index("--coupling-type") + 1] == "step"


def test_collect_simulator_provenance_hashes_executed_sources() -> None:
    repository = Path(__file__).parents[1]
    provenance = collect_simulator_provenance(
        python_executable=Path(sys.executable),
        simulation_script=repository / "examples" / "05_advanced_run.py",
    )

    assert Path(provenance["python_executable"]).is_absolute()
    assert Path(provenance["simulation_script"]).is_file()
    assert Path(provenance["cavitymd_core"]).is_file()
    assert len(provenance["simulation_script_sha256"]) == 64
    assert len(provenance["cavitymd_core_sha256"]) == 64


def test_replica_locks_exclude_duplicate_writer(tmp_path: Path) -> None:
    with replica_locks(tmp_path, [4, 5]):
        duplicate_handle = (tmp_path / ".replica_4.lock").open("a+")
        try:
            with pytest.raises(BlockingIOError):
                fcntl.flock(
                    duplicate_handle.fileno(),
                    fcntl.LOCK_EX | fcntl.LOCK_NB,
                )
        finally:
            duplicate_handle.close()

    released_handle = (tmp_path / ".replica_4.lock").open("a+")
    try:
        fcntl.flock(
            released_handle.fileno(),
            fcntl.LOCK_EX | fcntl.LOCK_NB,
        )
    finally:
        released_handle.close()


def test_main_skips_complete_replicas_without_cleanup(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("0p03\t0.03\t4\t5\n", encoding="utf-8")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_packed_aging_task.py",
            "--manifest",
            str(manifest),
            "--task-index",
            "0",
            "--output-base",
            str(tmp_path / "outputs"),
            "--log-dir",
            str(tmp_path / "logs"),
        ],
    )
    monkeypatch.setattr(
        packed_runner,
        "is_complete_replica",
        lambda _directory, _replica, **_kwargs: True,
    )
    monkeypatch.setattr(
        packed_runner,
        "cleanup_replica_artifacts",
        lambda *_args, **_kwargs: pytest.fail("cleanup must not run"),
    )
    monkeypatch.setattr(
        packed_runner,
        "run_concurrent_processes",
        lambda _specs: pytest.fail("children must not run"),
    )

    assert packed_runner.main() == 0


def test_main_postvalidates_both_children_after_process_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("0p03\t0.03\t4\t5\n", encoding="utf-8")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_packed_aging_task.py",
            "--manifest",
            str(manifest),
            "--task-index",
            "0",
            "--output-base",
            str(tmp_path / "outputs"),
            "--log-dir",
            str(tmp_path / "logs"),
        ],
    )
    validation_calls = {4: 0, 5: 0}

    validation_protocol_flags: list[bool] = []

    def fake_validation(
        _directory: Path,
        replica: int,
        **kwargs: object,
    ) -> bool:
        validation_calls[replica] += 1
        validation_protocol_flags.append(
            bool(kwargs.get("require_step_protocol"))
        )
        return validation_calls[replica] > 1 and replica == 4

    cleaned: list[int] = []
    marked: list[tuple[int, float, float, int, dict[str, str]]] = []
    events: list[str] = []

    @contextmanager
    def fake_locks(_directory: Path, replicas: tuple[int, ...]):
        events.append(f"lock:{replicas}")
        yield
        events.append("unlock")

    monkeypatch.setattr(
        packed_runner,
        "is_complete_replica",
        fake_validation,
    )
    monkeypatch.setattr(
        packed_runner,
        "cleanup_replica_artifacts",
        lambda _directory, replica: cleaned.append(replica) or [],
    )
    monkeypatch.setattr(
        packed_runner,
        "write_protocol_marker",
        lambda _directory, replica, *, lam, switch_time_ps, seed, provenance: (
            marked.append((replica, lam, switch_time_ps, seed, provenance))
        ),
    )
    provenance = {
        "python_executable": "/env/bin/python",
        "simulation_script": "/repo/examples/05_advanced_run.py",
        "simulation_script_sha256": "a" * 64,
        "cavitymd_core": "/env/core.py",
        "cavitymd_core_sha256": "b" * 64,
    }
    monkeypatch.setattr(
        packed_runner,
        "collect_simulator_provenance",
        lambda **_kwargs: provenance,
        raising=False,
    )
    monkeypatch.setattr(packed_runner, "replica_locks", fake_locks)

    def fake_run(_specs: list[ReplicaProcessSpec]) -> int:
        events.append("run")
        return 1

    monkeypatch.setattr(
        packed_runner,
        "run_concurrent_processes",
        fake_run,
    )

    assert packed_runner.main() == 1
    assert cleaned == [4, 5]
    assert marked == [(4, 0.03, 200.0, 5, provenance)]
    assert validation_calls == {4: 3, 5: 3}
    assert validation_protocol_flags == [
        True,
        True,
        False,
        False,
        True,
        True,
    ]
    assert events == ["lock:(4, 5)", "run", "unlock"]
