"""Tests for deterministic packed aging-campaign planning."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest

from examples.slurm.aging_campaign_planner import (
    PackedReplicaGroup,
    build_campaign_plan,
    campaign_plan_summary,
    interleave_groups,
    pair_replicas,
    select_replicas_to_target,
    write_manifest,
)
from examples.slurm.aging_campaign_status import fkt_file_path, run_dir


def _write_complete_replica(
    directory: Path,
    replica: int,
    *,
    min_bytes: int,
    runtime_ps: float,
    reference_interval_ps: float,
    max_references: int,
) -> None:
    """Write a compact scientifically complete planner fixture."""
    (directory / f"observables_replica_{replica}.h5").write_bytes(
        b"x" * min_bytes
    )
    for reference_index in range(max_references):
        reference_time_ps = reference_index * reference_interval_ps
        max_lag_ps = runtime_ps - reference_time_ps
        fkt_file_path(directory, replica, reference_index).write_text(
            "# F(k,t) correlation function\n"
            f"# Reference time: {reference_time_ps:.3f} ps\n"
            "# lag_time_ps\tF(k,t)\n"
            f"0.0\t1.0\n{max_lag_ps:.3f}\t0.1\n",
            encoding="utf-8",
        )


def test_select_replicas_to_target_uses_lowest_unique_candidates() -> None:
    selected = select_replicas_to_target(
        valid_replicas=[10, 11],
        candidate_replicas=[4, 2, 3, 2, 5],
        target_valid=5,
    )

    assert selected == (2, 3, 4)


def test_select_replicas_to_target_returns_empty_at_target() -> None:
    selected = select_replicas_to_target(
        valid_replicas=range(500),
        candidate_replicas=range(500, 1000),
        target_valid=500,
    )

    assert selected == ()


def test_select_replicas_to_target_returns_empty_above_target() -> None:
    selected = select_replicas_to_target(
        valid_replicas=range(680),
        candidate_replicas=range(680, 1000),
        target_valid=500,
    )

    assert selected == ()


def test_select_replicas_to_target_rejects_insufficient_candidates() -> None:
    with pytest.raises(ValueError, match="insufficient"):
        select_replicas_to_target(
            valid_replicas=[0],
            candidate_replicas=[1],
            target_valid=3,
        )


def test_pair_replicas_groups_even_count_with_same_lambda() -> None:
    groups = pair_replicas(0.03, [2, 3, 4, 5], group_size=2)

    assert groups == (
        PackedReplicaGroup("0p03", 0.03, (2, 3)),
        PackedReplicaGroup("0p03", 0.03, (4, 5)),
    )


def test_pair_replicas_keeps_odd_final_singleton() -> None:
    groups = pair_replicas(0.023333, [2, 3, 4], group_size=2)

    assert groups[-1].replicas == (4,)
    assert all(group.lam == 0.023333 for group in groups)


def test_pair_replicas_rejects_duplicate_ids() -> None:
    with pytest.raises(ValueError, match="duplicate"):
        pair_replicas(0.03, [2, 2], group_size=2)


def test_pair_replicas_rejects_more_than_two_per_task() -> None:
    with pytest.raises(ValueError, match="at most two"):
        pair_replicas(0.03, [1, 2, 3], group_size=3)


def test_interleave_groups_round_robins_couplings() -> None:
    first = pair_replicas(0.016667, [0, 1, 2, 3], group_size=2)
    second = pair_replicas(0.023333, [4, 5], group_size=2)
    third = pair_replicas(0.03, [6, 7, 8, 9], group_size=2)

    interleaved = interleave_groups([first, second, third])

    assert [group.lam for group in interleaved] == [
        0.016667,
        0.023333,
        0.03,
        0.016667,
        0.03,
    ]


def test_write_manifest_serializes_pairs_and_singletons(
    tmp_path: Path,
) -> None:
    groups = (
        PackedReplicaGroup("0p016667", 0.016667, (2, 3)),
        PackedReplicaGroup("0p03", 0.03, (4,)),
    )
    manifest = tmp_path / "aging.tsv"

    write_manifest(groups, manifest)

    assert manifest.read_text(encoding="utf-8") == (
        "0p016667\t0.016667\t2\t3\n"
        "0p03\t0.03\t4\t\n"
    )


def test_write_manifest_is_deterministic(tmp_path: Path) -> None:
    groups = pair_replicas(0.03, [5, 2, 4, 3], group_size=2)
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"

    write_manifest(groups, first)
    write_manifest(groups, second)

    assert first.read_bytes() == second.read_bytes()


def test_build_campaign_plan_counts_valid_outputs_and_selects_exact_top_up(
    tmp_path: Path,
) -> None:
    lam = 0.03
    directory = run_dir(tmp_path, lam)
    directory.mkdir(parents=True)
    for replica in (0, 1):
        _write_complete_replica(
            directory,
            replica,
            min_bytes=16,
            runtime_ps=4.0,
            reference_interval_ps=1.0,
            max_references=3,
        )

    plan = build_campaign_plan(
        output_base=tmp_path,
        lambdas=[lam],
        target_valid=3,
        max_replica_id=4,
        replicas_per_task=2,
        min_bytes=16,
        runtime_ps=4.0,
        reference_interval_ps=1.0,
        max_references=3,
    )

    lambda_plan = plan.lambda_plans[0]
    assert lambda_plan.scan.complete_replicas == (0, 1)
    assert lambda_plan.selected_replicas == (2,)
    assert plan.groups[0].replicas == (2,)


def test_build_campaign_plan_rejects_duplicate_lambdas(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="duplicate coupling"):
        build_campaign_plan(
            output_base=tmp_path,
            lambdas=[0.03, 0.03],
            target_valid=1,
            max_replica_id=2,
            min_bytes=16,
            runtime_ps=4.0,
            reference_interval_ps=1.0,
            max_references=3,
        )


@pytest.mark.parametrize(
    "invalid_lambda",
    [float("nan"), float("inf"), float("-inf")],
)
def test_build_campaign_plan_rejects_nonfinite_lambda(
    tmp_path: Path,
    invalid_lambda: float,
) -> None:
    with pytest.raises(ValueError, match="finite"):
        build_campaign_plan(
            output_base=tmp_path,
            lambdas=[invalid_lambda],
            target_valid=1,
            max_replica_id=2,
            min_bytes=16,
            runtime_ps=4.0,
            reference_interval_ps=1.0,
            max_references=3,
        )


def test_campaign_plan_summary_reports_valid_selected_and_groups(
    tmp_path: Path,
) -> None:
    plan = build_campaign_plan(
        output_base=tmp_path,
        lambdas=[0.03],
        target_valid=2,
        max_replica_id=2,
        min_bytes=16,
        runtime_ps=4.0,
        reference_interval_ps=1.0,
        max_references=3,
    )

    assert campaign_plan_summary(plan) == [
        {
            "lam": 0.03,
            "lambda_tag": "0p03",
            "run_dir": str(run_dir(tmp_path, 0.03)),
            "valid": 0,
            "invalid": 3,
            "target_valid": 2,
            "selected": 2,
            "selected_replicas": [0, 1],
            "groups": 1,
            "pairs": 1,
            "singletons": 0,
        }
    ]


def test_planner_supports_direct_cli_invocation() -> None:
    script = (
        Path(__file__).parents[1]
        / "examples"
        / "slurm"
        / "aging_campaign_planner.py"
    )

    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "--target-valid" in result.stdout
