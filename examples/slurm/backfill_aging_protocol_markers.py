#!/usr/bin/env python3
"""Backfill step-coupling protocol markers for data-complete aging replicas."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Mapping, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from examples.slurm.aging_campaign_status import (
        COMPLETE_MIN_BYTES,
        DEFAULT_FKT_MAX_REFERENCES,
        DEFAULT_FKT_REFERENCE_INTERVAL_PS,
        DEFAULT_RUNTIME_PS,
        DEFAULT_SWITCH_TIME_PS,
        is_complete_replica,
        protocol_marker_path,
        run_dir,
        write_protocol_marker,
    )
    from examples.slurm.run_packed_aging_task import collect_simulator_provenance
else:
    from .aging_campaign_status import (
        COMPLETE_MIN_BYTES,
        DEFAULT_FKT_MAX_REFERENCES,
        DEFAULT_FKT_REFERENCE_INTERVAL_PS,
        DEFAULT_RUNTIME_PS,
        DEFAULT_SWITCH_TIME_PS,
        is_complete_replica,
        protocol_marker_path,
        run_dir,
        write_protocol_marker,
    )
    from .run_packed_aging_task import collect_simulator_provenance

DEFAULT_LAMBDAS = (0.0, 0.01, 0.016667, 0.023333, 0.03)
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def find_protocol_backfill_candidates(
    output_base: Path,
    lambdas: Sequence[float],
    *,
    n_replicas: int,
    runtime_ps: float = DEFAULT_RUNTIME_PS,
    reference_interval_ps: float = DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    max_references: int = DEFAULT_FKT_MAX_REFERENCES,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    min_bytes: int = COMPLETE_MIN_BYTES,
) -> list[tuple[float, int]]:
    """Return replicas that are data-complete but missing a valid protocol marker."""
    candidates: list[tuple[float, int]] = []
    completion_kwargs = {
        "min_bytes": min_bytes,
        "runtime_ps": runtime_ps,
        "reference_interval_ps": reference_interval_ps,
        "max_references": max_references,
    }
    for lam in lambdas:
        if lam == 0.0:
            continue
        directory = run_dir(output_base, lam, switch_time_ps=switch_time_ps)
        for replica in range(n_replicas):
            if not is_complete_replica(
                directory,
                replica,
                require_step_protocol=False,
                **completion_kwargs,
            ):
                continue
            if protocol_marker_path(directory, replica).is_file():
                continue
            if is_complete_replica(
                directory,
                replica,
                require_step_protocol=True,
                expected_lam=lam,
                expected_switch_time_ps=switch_time_ps,
                **completion_kwargs,
            ):
                continue
            candidates.append((lam, replica))
    candidates.sort()
    return candidates


def backfill_protocol_markers(
    output_base: Path,
    lambdas: Sequence[float],
    *,
    provenance: Mapping[str, str],
    n_replicas: int,
    runtime_ps: float = DEFAULT_RUNTIME_PS,
    reference_interval_ps: float = DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    max_references: int = DEFAULT_FKT_MAX_REFERENCES,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    min_bytes: int = COMPLETE_MIN_BYTES,
    dry_run: bool = False,
) -> list[tuple[float, int]]:
    """Write missing protocol markers for data-complete replicas."""
    candidates = find_protocol_backfill_candidates(
        output_base,
        lambdas,
        n_replicas=n_replicas,
        runtime_ps=runtime_ps,
        reference_interval_ps=reference_interval_ps,
        max_references=max_references,
        switch_time_ps=switch_time_ps,
        min_bytes=min_bytes,
    )
    if dry_run:
        return candidates

    written: list[tuple[float, int]] = []
    for lam, replica in candidates:
        directory = run_dir(output_base, lam, switch_time_ps=switch_time_ps)
        write_protocol_marker(
            directory,
            replica,
            lam=lam,
            switch_time_ps=switch_time_ps,
            seed=replica + 1,
            provenance=provenance,
        )
        written.append((lam, replica))
    return written


def _parse_lambdas(raw: str) -> tuple[float, ...]:
    values = tuple(float(item.strip()) for item in raw.split(",") if item.strip())
    if not values:
        raise ValueError("at least one lambda value is required")
    return values


def build_parser() -> argparse.ArgumentParser:
    """Build the backfill CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-base",
        type=Path,
        default=REPOSITORY_ROOT / "aging_weak_lambda",
    )
    parser.add_argument(
        "--examples-dir",
        type=Path,
        default=REPOSITORY_ROOT / "examples",
    )
    parser.add_argument(
        "--python-executable",
        type=Path,
        default=Path(sys.executable),
    )
    parser.add_argument(
        "--lambdas",
        type=_parse_lambdas,
        default=",".join(str(value) for value in DEFAULT_LAMBDAS),
    )
    parser.add_argument("--scan-replicas", type=int, default=500)
    parser.add_argument("--switch-time-ps", type=float, default=DEFAULT_SWITCH_TIME_PS)
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
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main() -> int:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()
    if args.scan_replicas <= 0:
        parser.error("--scan-replicas must be positive")

    provenance = collect_simulator_provenance(
        python_executable=args.python_executable,
        simulation_script=args.examples_dir / "05_advanced_run.py",
    )
    written = backfill_protocol_markers(
        args.output_base,
        args.lambdas,
        provenance=provenance,
        n_replicas=args.scan_replicas,
        runtime_ps=args.runtime_ps,
        reference_interval_ps=args.reference_interval_ps,
        max_references=args.max_references,
        switch_time_ps=args.switch_time_ps,
        min_bytes=args.min_bytes,
        dry_run=args.dry_run,
    )
    action = "would backfill" if args.dry_run else "backfilled"
    print(f"{action} {len(written)} protocol markers")
    for lam, replica in written:
        print(f"  lam={lam:g} replica={replica}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
