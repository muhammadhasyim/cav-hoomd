"""Plan packed SLURM tasks for 1600→2000 ps aging extensions."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from examples.slurm.aging_campaign_planner import (
        DEFAULT_LAMBDAS,
        DEFAULT_MAX_REPLICA_ID,
        DEFAULT_REPLICAS_PER_TASK,
        DEFAULT_TARGET_VALID,
        CampaignPlan,
        LambdaPlan,
        PackedReplicaGroup,
        campaign_plan_summary,
        interleave_groups,
        pair_replicas,
        write_manifest,
    )
    from examples.slurm.aging_campaign_status import (
        COMPLETE_MIN_BYTES,
        DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
        DEFAULT_FKT_MAX_REFERENCES,
        DEFAULT_FKT_REFERENCE_INTERVAL_PS,
        DEFAULT_RUNTIME_PS,
        DEFAULT_SWITCH_TIME_PS,
        lambda_to_tag,
        scan_lambda_run,
    )
else:
    from .aging_campaign_planner import (
        DEFAULT_LAMBDAS,
        DEFAULT_MAX_REPLICA_ID,
        DEFAULT_REPLICAS_PER_TASK,
        DEFAULT_TARGET_VALID,
        CampaignPlan,
        LambdaPlan,
        PackedReplicaGroup,
        campaign_plan_summary,
        interleave_groups,
        pair_replicas,
        write_manifest,
    )
    from .aging_campaign_status import (
        COMPLETE_MIN_BYTES,
        DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
        DEFAULT_FKT_MAX_REFERENCES,
        DEFAULT_FKT_REFERENCE_INTERVAL_PS,
        DEFAULT_RUNTIME_PS,
        DEFAULT_SWITCH_TIME_PS,
        lambda_to_tag,
        scan_lambda_run,
    )

DEFAULT_OUTPUT_BASE = Path(
    "/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda"
)


@dataclass(frozen=True)
class ExtensionCampaignPlan:
    """Extension-only plan for replicas that finished the primary runtime."""

    lambda_plans: tuple[LambdaPlan, ...]
    groups: tuple[PackedReplicaGroup, ...]


def build_extension_plan(
    output_base: Path,
    lambdas: Sequence[float] = DEFAULT_LAMBDAS,
    *,
    max_replica_id: int = DEFAULT_MAX_REPLICA_ID,
    replicas_per_task: int = DEFAULT_REPLICAS_PER_TASK,
    min_bytes: int = COMPLETE_MIN_BYTES,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    runtime_ps: float = DEFAULT_RUNTIME_PS,
    reference_interval_ps: float = DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    max_references: int = DEFAULT_FKT_MAX_REFERENCES,
    max_lag_tolerance_ps: float = DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
) -> ExtensionCampaignPlan:
    """Build a manifest for replicas that need a checkpointed extension."""
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
            require_step_protocol=lam != 0.0,
            expected_lam=lam,
            expected_switch_time_ps=switch_time_ps,
        )
        selected = scan.extension_replicas
        groups = pair_replicas(
            lam,
            selected,
            group_size=replicas_per_task,
        )
        lambda_plans.append(
            LambdaPlan(
                scan=scan,
                target_valid=len(selected),
                selected_replicas=selected,
                groups=groups,
            )
        )
        grouped.append(groups)

    return ExtensionCampaignPlan(
        lambda_plans=tuple(lambda_plans),
        groups=interleave_groups(grouped),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT_BASE)
    parser.add_argument("--lambdas", type=float, nargs="+", default=list(DEFAULT_LAMBDAS))
    parser.add_argument("--max-replica-id", type=int, default=DEFAULT_MAX_REPLICA_ID)
    parser.add_argument("--replicas-per-task", type=int, default=DEFAULT_REPLICAS_PER_TASK)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    plan = build_extension_plan(
        output_base=args.output_base,
        lambdas=args.lambdas,
        max_replica_id=args.max_replica_id,
        replicas_per_task=args.replicas_per_task,
    )
    write_manifest(plan.groups, args.manifest)
    summary = campaign_plan_summary(
        CampaignPlan(lambda_plans=plan.lambda_plans, groups=plan.groups)
    )
    if args.json:
        payload = {
            "couplings": summary,
            "total_groups": len(plan.groups),
            "manifest": str(args.manifest),
        }
        print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
