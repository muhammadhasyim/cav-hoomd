"""Plan packed SLURM tasks for the weak-coupling aging campaign."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from examples.slurm.aging_campaign_status import (
        COMPLETE_MIN_BYTES,
        DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
        DEFAULT_FKT_MAX_REFERENCES,
        DEFAULT_FKT_REFERENCE_INTERVAL_PS,
        DEFAULT_RUNTIME_PS,
        DEFAULT_SWITCH_TIME_PS,
        LambdaScanResult,
        lambda_to_tag,
        scan_lambda_run,
    )
else:
    from .aging_campaign_status import (
        COMPLETE_MIN_BYTES,
        DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
        DEFAULT_FKT_MAX_REFERENCES,
        DEFAULT_FKT_REFERENCE_INTERVAL_PS,
        DEFAULT_RUNTIME_PS,
        DEFAULT_SWITCH_TIME_PS,
        LambdaScanResult,
        lambda_to_tag,
        scan_lambda_run,
    )

DEFAULT_OUTPUT_BASE = Path(
    "/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda"
)
DEFAULT_LAMBDAS = (0.0, 0.01, 0.016667, 0.023333, 0.03)
DEFAULT_TARGET_VALID = 500
DEFAULT_MAX_REPLICA_ID = 999
DEFAULT_REPLICAS_PER_TASK = 2


@dataclass(frozen=True)
class PackedReplicaGroup:
    """One array task containing replicas from a single coupling value."""

    lambda_tag: str
    lam: float
    replicas: tuple[int, ...]

    def __post_init__(self) -> None:
        """Validate IDs when constructing an immutable task group."""
        if not self.replicas:
            raise ValueError("a packed group requires at least one replica")
        if len(set(self.replicas)) != len(self.replicas):
            raise ValueError("a packed group cannot contain duplicate replicas")
        if any(replica < 0 for replica in self.replicas):
            raise ValueError("replica IDs must be nonnegative")


@dataclass(frozen=True)
class LambdaPlan:
    """Robust scan and exact top-up plan for one coupling value."""

    scan: LambdaScanResult
    target_valid: int
    selected_replicas: tuple[int, ...]
    groups: tuple[PackedReplicaGroup, ...]

    @property
    def singleton_count(self) -> int:
        """Return the number of groups containing one replica."""
        return sum(len(group.replicas) == 1 for group in self.groups)


@dataclass(frozen=True)
class CampaignPlan:
    """Complete packed plan for all requested coupling values."""

    lambda_plans: tuple[LambdaPlan, ...]
    groups: tuple[PackedReplicaGroup, ...]


def _validated_unique_ids(
    replica_ids: Iterable[int],
    *,
    description: str,
) -> tuple[int, ...]:
    """Return sorted unique nonnegative replica IDs.

    Parameters
    ----------
    replica_ids
        Replica IDs to validate.
    description
        Name used in validation error messages.

    Returns
    -------
    tuple[int, ...]
        Sorted unique IDs.
    """
    values = tuple(replica_ids)
    if any(replica < 0 for replica in values):
        raise ValueError(f"{description} contain a negative replica ID")
    return tuple(sorted(set(values)))


def select_replicas_to_target(
    valid_replicas: Iterable[int],
    candidate_replicas: Iterable[int],
    *,
    target_valid: int = DEFAULT_TARGET_VALID,
) -> tuple[int, ...]:
    """Select the lowest invalid IDs needed to reach a valid target.

    Parameters
    ----------
    valid_replicas
        Scientifically complete replica IDs.
    candidate_replicas
        Missing or invalid IDs eligible for execution.
    target_valid
        Desired total number of valid replicas.

    Returns
    -------
    tuple[int, ...]
        Exact deterministic top-up selection.

    Raises
    ------
    ValueError
        If inputs are invalid or candidates cannot reach the target.
    """
    if target_valid < 0:
        raise ValueError("target_valid must be nonnegative")
    valid = _validated_unique_ids(valid_replicas, description="valid replicas")
    candidates = _validated_unique_ids(
        candidate_replicas,
        description="candidate replicas",
    )
    valid_set = set(valid)
    available = tuple(replica for replica in candidates if replica not in valid_set)
    needed = max(0, target_valid - len(valid))
    if len(available) < needed:
        raise ValueError(
            "insufficient candidate replicas to reach the valid target"
        )
    return available[:needed]


def pair_replicas(
    lam: float,
    replica_ids: Sequence[int],
    *,
    group_size: int = DEFAULT_REPLICAS_PER_TASK,
) -> tuple[PackedReplicaGroup, ...]:
    """Group sorted replica IDs for concurrent execution.

    Parameters
    ----------
    lam
        Dimensionless campaign coupling value.
    replica_ids
        IDs selected for this coupling.
    group_size
        Maximum replicas assigned to one array task.

    Returns
    -------
    tuple[PackedReplicaGroup, ...]
        Same-coupling groups, including a final singleton when necessary.
    """
    if not 1 <= group_size <= DEFAULT_REPLICAS_PER_TASK:
        raise ValueError("group_size must assign at most two replicas")
    values = tuple(replica_ids)
    if len(set(values)) != len(values):
        raise ValueError("replica IDs contain a duplicate")
    ordered = _validated_unique_ids(values, description="replica IDs")
    tag = lambda_to_tag(lam)
    return tuple(
        PackedReplicaGroup(tag, lam, ordered[index : index + group_size])
        for index in range(0, len(ordered), group_size)
    )


def interleave_groups(
    groups_by_lambda: Sequence[Sequence[PackedReplicaGroup]],
) -> tuple[PackedReplicaGroup, ...]:
    """Round-robin task groups while preserving each coupling's order.

    Parameters
    ----------
    groups_by_lambda
        Ordered group sequences, one sequence per coupling value.

    Returns
    -------
    tuple[PackedReplicaGroup, ...]
        Deterministically interleaved groups.
    """
    queues = [deque(groups) for groups in groups_by_lambda if groups]
    interleaved: list[PackedReplicaGroup] = []
    while queues:
        remaining: list[deque[PackedReplicaGroup]] = []
        for queue in queues:
            interleaved.append(queue.popleft())
            if queue:
                remaining.append(queue)
        queues = remaining
    return tuple(interleaved)


def write_manifest(
    groups: Sequence[PackedReplicaGroup],
    path: Path,
) -> None:
    """Write deterministic tab-separated packed-task rows.

    Parameters
    ----------
    groups
        Interleaved task groups.
    path
        Destination manifest path.

    Side Effects
    ------------
    Creates parent directories and replaces the manifest atomically enough for
    submission planning on the shared filesystem.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for group in groups:
        if len(group.replicas) > DEFAULT_REPLICAS_PER_TASK:
            raise ValueError("manifest groups may contain at most two replicas")
        replicas = [str(replica) for replica in group.replicas]
        replicas.extend([""] * (DEFAULT_REPLICAS_PER_TASK - len(replicas)))
        lines.append(
            "\t".join(
                [group.lambda_tag, f"{group.lam:g}", *replicas]
            )
        )
    content = "".join(f"{line}\n" for line in lines)
    path.write_text(content, encoding="utf-8")


def build_campaign_plan(
    output_base: Path,
    lambdas: Sequence[float] = DEFAULT_LAMBDAS,
    *,
    target_valid: int = DEFAULT_TARGET_VALID,
    max_replica_id: int = DEFAULT_MAX_REPLICA_ID,
    replicas_per_task: int = DEFAULT_REPLICAS_PER_TASK,
    min_bytes: int = COMPLETE_MIN_BYTES,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    runtime_ps: float = DEFAULT_RUNTIME_PS,
    reference_interval_ps: float = DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    max_references: int = DEFAULT_FKT_MAX_REFERENCES,
    max_lag_tolerance_ps: float = DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
) -> CampaignPlan:
    """Build a robust, exact top-up plan for the requested couplings.

    Parameters are expressed in campaign units: times are picoseconds and
    coupling values are dimensionless lambda values.

    Returns
    -------
    CampaignPlan
        Per-coupling scans and one globally interleaved task sequence.
    """
    if max_replica_id < 0:
        raise ValueError("max_replica_id must be nonnegative")
    if any(not math.isfinite(lam) for lam in lambdas):
        raise ValueError("coupling values must be finite")
    if len(set(lambdas)) != len(lambdas):
        raise ValueError("duplicate coupling values are not allowed")
    lambda_plans: list[LambdaPlan] = []
    grouped: list[tuple[PackedReplicaGroup, ...]] = []

    for lam in lambdas:
        scan = scan_lambda_run(
            output_base=output_base,
            lam=lam,
            n_replicas=max_replica_id + 1,
            switch_time_ps=switch_time_ps,
            min_bytes=min_bytes,
            require_fkt=True,
            runtime_ps=runtime_ps,
            reference_interval_ps=reference_interval_ps,
            max_references=max_references,
            max_lag_tolerance_ps=max_lag_tolerance_ps,
        )
        selected = select_replicas_to_target(
            scan.complete_replicas,
            scan.missing_replicas,
            target_valid=target_valid,
        )
        groups = pair_replicas(
            lam,
            selected,
            group_size=replicas_per_task,
        )
        lambda_plans.append(
            LambdaPlan(
                scan=scan,
                target_valid=target_valid,
                selected_replicas=selected,
                groups=groups,
            )
        )
        grouped.append(groups)

    return CampaignPlan(
        lambda_plans=tuple(lambda_plans),
        groups=interleave_groups(grouped),
    )


def campaign_plan_summary(plan: CampaignPlan) -> list[dict[str, object]]:
    """Return a JSON-serializable per-coupling plan summary."""
    return [
        {
            "lam": lambda_plan.scan.lam,
            "lambda_tag": lambda_to_tag(lambda_plan.scan.lam),
            "run_dir": str(lambda_plan.scan.run_dir),
            "valid": lambda_plan.scan.n_complete,
            "invalid": lambda_plan.scan.n_missing,
            "target_valid": lambda_plan.target_valid,
            "selected": len(lambda_plan.selected_replicas),
            "selected_replicas": list(lambda_plan.selected_replicas),
            "groups": len(lambda_plan.groups),
            "pairs": len(lambda_plan.groups) - lambda_plan.singleton_count,
            "singletons": lambda_plan.singleton_count,
        }
        for lambda_plan in plan.lambda_plans
    ]


def _parse_args() -> argparse.Namespace:
    """Parse planner command-line options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT_BASE)
    parser.add_argument(
        "--lambdas",
        nargs="+",
        type=float,
        default=list(DEFAULT_LAMBDAS),
    )
    parser.add_argument("--target-valid", type=int, default=DEFAULT_TARGET_VALID)
    parser.add_argument(
        "--max-replica-id",
        type=int,
        default=DEFAULT_MAX_REPLICA_ID,
    )
    parser.add_argument(
        "--replicas-per-task",
        type=int,
        default=DEFAULT_REPLICAS_PER_TASK,
    )
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()


def main() -> int:
    """Build a campaign plan, optionally write its manifest, and report it."""
    args = _parse_args()
    plan = build_campaign_plan(
        output_base=args.output_base,
        lambdas=args.lambdas,
        target_valid=args.target_valid,
        max_replica_id=args.max_replica_id,
        replicas_per_task=args.replicas_per_task,
    )
    if args.manifest is not None:
        write_manifest(plan.groups, args.manifest)

    summary = campaign_plan_summary(plan)
    if args.json:
        print(
            json.dumps(
                {
                    "couplings": summary,
                    "total_groups": len(plan.groups),
                    "manifest": (
                        str(args.manifest)
                        if args.manifest is not None
                        else None
                    ),
                },
                indent=2,
            )
        )
    else:
        for item in summary:
            print(
                f"lam={item['lam']}: valid={item['valid']} "
                f"invalid={item['invalid']} selected={item['selected']} "
                f"pairs={item['pairs']} singletons={item['singletons']} "
                f"groups={item['groups']}"
            )
            if item["selected_replicas"]:
                print(
                    "  selected_replicas="
                    + ",".join(
                        str(replica)
                        for replica in item["selected_replicas"]
                    )
                )
        print(f"total_groups={len(plan.groups)}")
        if args.manifest is not None:
            print(f"manifest={args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
