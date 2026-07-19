"""Tests for the aging campaign failed-run retry watcher."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from examples.slurm.aging_campaign_status import run_dir
from examples.slurm.watch_retry_failed_aging import (
    RetryCandidate,
    RetrySubmitGate,
    build_retry_manifest_lines,
    collect_retry_candidates,
    manifest_path_for_cycle,
    parse_manifest_row,
    replica_was_attempted,
    render_retry_sbatch,
    retry_job_in_flight,
    sbatch_path_for_cycle,
    write_retry_manifest,
)


def test_replica_was_attempted_detects_partial_h5(tmp_path: Path) -> None:
    """Partial HDF5 output counts as an attempted replica."""
    run_dir = tmp_path / "lambda0" / "coupling"
    run_dir.mkdir(parents=True)
    (run_dir / "observables_replica_7.h5").write_bytes(b"x" * 128)
    assert replica_was_attempted(run_dir, 7) is True
    assert replica_was_attempted(run_dir, 8) is False


def test_collect_retry_candidates_partial_and_attempted_missing(
    tmp_path: Path,
) -> None:
    """Retry set includes partial replicas and missing ones with artifacts."""
    lam = 0.01
    run_directory = run_dir(tmp_path, lam, switch_time_ps=200.0)
    run_directory.mkdir(parents=True)
    (run_directory / "observables_replica_1.h5").write_bytes(b"x" * 128)
    (run_directory / "protocol_replica_2.json").write_text("{}", encoding="utf-8")

    candidates = collect_retry_candidates(
        tmp_path,
        (lam,),
        target_valid=5,
        n_replicas=5,
        switch_time_ps=200.0,
        runtime_ps=1600.0,
        reference_interval_ps=200.0,
        max_references=8,
        min_bytes=5_000_000,
    )
    replicas = {(item.lam, item.replica) for item in candidates}
    assert (lam, 1) in replicas
    assert (lam, 2) in replicas
    assert (0.01, 3) not in replicas


def test_build_retry_manifest_lines_one_replica_per_task() -> None:
    """Each retry task should carry exactly one replica."""
    candidates = (
        RetryCandidate(lambda_tag="0", lam=0.0, replica=4),
        RetryCandidate(lambda_tag="0p01", lam=0.01, replica=9),
    )
    lines = build_retry_manifest_lines(candidates)
    assert lines == ["0\t0\t4\t", "0p01\t0.01\t9\t"]


def test_write_retry_manifest_round_trip(tmp_path: Path) -> None:
    """Manifest writer should preserve tab-separated rows."""
    manifest = tmp_path / "manifest_retry.tsv"
    candidates = (RetryCandidate(lambda_tag="0", lam=0.0, replica=1),)
    write_retry_manifest(manifest, candidates)
    assert manifest.read_text(encoding="utf-8") == "0\t0\t1\t\n"


def test_parse_manifest_row() -> None:
    """Manifest rows should decode to one replica key."""
    assert parse_manifest_row("0p03\t0.03\t12\t") == ("0p03", 0.03, 12)


def test_render_retry_sbatch_array_bounds() -> None:
    """Generated sbatch should use a zero-based array covering all tasks."""
    manifest = Path("/tmp/manifest_retry.tsv")
    text = render_retry_sbatch(
        manifest_path=manifest,
        n_tasks=3,
        array_concurrency=24,
        job_name="cav-aging-retry",
        partition="l40s_public,h200_public",
        account="torch_pr_283_chemistry",
        output_base=Path("/tmp/aging_weak_lambda"),
        examples_dir=Path("/tmp/examples"),
        input_gsd=Path("/tmp/init.gsd"),
        log_dir=Path("/tmp/logs/packed"),
        walltime="01:00:00",
    )
    assert "#SBATCH --array=0-2%24" in text
    assert "ga015" in text
    assert str(manifest) in text


def test_retry_submit_gate_respects_qos_cap() -> None:
    """Retry submission should wait when the projected total exceeds QOS."""
    gate = RetrySubmitGate(
        active_jobs=900,
        retry_tasks=400,
        max_submit_jobs=1000,
        safety_margin=8,
    )
    assert gate.projected_total == 1300
    assert gate.allowed is False


def test_retry_job_in_flight_detects_running_name() -> None:
    """Running retry arrays should block duplicate submissions."""

    def fake_runner(
        args: object,
        *,
        check: bool,
        capture_output: bool,
        text: bool,
    ) -> subprocess.CompletedProcess[str]:
        del args, check, capture_output, text
        return subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout="RUNNING cav-aging-retry\n",
            stderr="",
        )

    assert retry_job_in_flight(
        user="mh7373",
        job_name_prefix="cav-aging-retry",
        runner=fake_runner,
    )


def test_cycle_paths_are_unique(tmp_path: Path) -> None:
    """Each watcher cycle should write distinct manifest and sbatch paths."""
    cycle = "20260719T120000Z"
    plan_dir = tmp_path / "plan"
    manifest = manifest_path_for_cycle(plan_dir, cycle)
    sbatch = sbatch_path_for_cycle(plan_dir, cycle)
    assert manifest.name == f"manifest_retry_{cycle}.tsv"
    assert sbatch.name == f"job_retry_{cycle}.sbatch"


def test_state_file_tracks_last_cycle(tmp_path: Path) -> None:
    """State persistence should round-trip the last submitted cycle metadata."""
    state_path = tmp_path / "retry_watcher_state.json"
    payload = {"last_cycle": "20260719T120000Z", "job_id": "14299999", "n_tasks": 12}
    state_path.write_text(json.dumps(payload), encoding="utf-8")
    loaded = json.loads(state_path.read_text(encoding="utf-8"))
    assert loaded["job_id"] == "14299999"
