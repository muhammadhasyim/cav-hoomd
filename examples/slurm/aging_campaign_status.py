"""Utilities for weak-coupling aging campaign resume and cleanup."""

from __future__ import annotations

import argparse
import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

COMPLETE_MIN_BYTES = 1_380_000_000
DEFAULT_SWITCH_TIME_PS = 200.0
DEFAULT_N_REPLICAS = 1000
REPLICA_H5_PATTERN = re.compile(r"observables_replica_(\d+)\.h5$")


@dataclass(frozen=True)
class LambdaScanResult:
    """Replica inventory for one lambda value."""

    lam: float
    run_dir: Path
    complete_replicas: tuple[int, ...]
    partial_replicas: tuple[int, ...]
    missing_replicas: tuple[int, ...]

    @property
    def n_complete(self) -> int:
        return len(self.complete_replicas)

    @property
    def n_partial(self) -> int:
        return len(self.partial_replicas)

    @property
    def n_missing(self) -> int:
        return len(self.missing_replicas)


def lambda_to_tag(lam: float) -> str:
    """Convert a lambda value to the directory tag used by the campaign."""
    if lam == 0.0:
        return "0"
    return str(lam).replace(".", "p")


def coupling_dir_name(coupling: float, switch_time_ps: float) -> str:
    """Match the directory naming logic in examples/05_advanced_run.py."""
    coupling_str = f"{coupling:.0e}".replace("-", "neg").replace("+", "pos")
    return f"cavity_coupling_{coupling_str}_switch_{switch_time_ps}ps"


def run_dir(
    output_base: Path,
    lam: float,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
) -> Path:
    """Return the expected run directory for one lambda value."""
    return output_base / f"lambda{lambda_to_tag(lam)}" / coupling_dir_name(lam, switch_time_ps)


def is_complete_h5(path: Path, min_bytes: int = COMPLETE_MIN_BYTES) -> bool:
    """Return True when an observables HDF5 file looks like a finished replica."""
    return path.is_file() and path.stat().st_size >= min_bytes


def replica_h5_path(run_directory: Path, replica: int) -> Path:
    """Return the observables HDF5 path for one replica."""
    return run_directory / f"observables_replica_{replica}.h5"


def replica_gsd_path(run_directory: Path, replica: int) -> Path:
    """Return the production GSD path for one replica."""
    return run_directory / f"prod-{replica}.gsd"


def slurm_array_spec(replica_ids: list[int]) -> str:
    """Convert replica IDs into a compact SLURM --array specification."""
    if not replica_ids:
        return ""

    ranges: list[str] = []
    start = prev = replica_ids[0]
    for replica in replica_ids[1:]:
        if replica == prev + 1:
            prev = replica
            continue
        ranges.append(f"{start}-{prev}" if start != prev else str(start))
        start = prev = replica
    ranges.append(f"{start}-{prev}" if start != prev else str(start))
    return ",".join(ranges)


def scan_lambda_run(
    output_base: Path,
    lam: float,
    n_replicas: int = DEFAULT_N_REPLICAS,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    min_bytes: int = COMPLETE_MIN_BYTES,
) -> LambdaScanResult:
    """Scan one lambda directory and classify replica outputs."""
    directory = run_dir(output_base, lam, switch_time_ps)
    complete: list[int] = []
    partial: list[int] = []

    if directory.is_dir():
        for path in directory.glob("observables_replica_*.h5"):
            match = REPLICA_H5_PATTERN.match(path.name)
            if match is None:
                continue
            replica = int(match.group(1))
            if is_complete_h5(path, min_bytes=min_bytes):
                complete.append(replica)
            else:
                partial.append(replica)

    complete.sort()
    partial.sort()
    missing = [replica for replica in range(n_replicas) if replica not in complete]

    return LambdaScanResult(
        lam=lam,
        run_dir=directory,
        complete_replicas=tuple(complete),
        partial_replicas=tuple(partial),
        missing_replicas=tuple(missing),
    )


def cleanup_partial_replicas(
    scan: LambdaScanResult,
    dry_run: bool = False,
) -> list[Path]:
    """Delete partial/failed replica artifacts while preserving complete outputs."""
    deleted: list[Path] = []

    for replica in scan.partial_replicas:
        for path in (
            replica_h5_path(scan.run_dir, replica),
            replica_gsd_path(scan.run_dir, replica),
        ):
            if not path.exists():
                continue
            deleted.append(path)
            if not dry_run:
                path.unlink()

    return deleted


def cleanup_stray_run_dirs(
    output_base: Path,
    lam: float,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    dry_run: bool = False,
) -> list[Path]:
    """Remove run directories under a lambda folder that do not match the campaign."""
    lam_dir = output_base / f"lambda{lambda_to_tag(lam)}"
    expected = coupling_dir_name(lam, switch_time_ps)
    removed: list[Path] = []

    if not lam_dir.is_dir():
        return removed

    for child in lam_dir.iterdir():
        if not child.is_dir() or child.name == expected:
            continue
        removed.append(child)
        if not dry_run:
            shutil.rmtree(child)

    return removed


def scan_campaign(
    output_base: Path,
    lambdas: list[float],
    n_replicas: int = DEFAULT_N_REPLICAS,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
) -> list[LambdaScanResult]:
    """Scan all lambda values in the campaign."""
    return [
        scan_lambda_run(
            output_base=output_base,
            lam=lam,
            n_replicas=n_replicas,
            switch_time_ps=switch_time_ps,
        )
        for lam in lambdas
    ]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-base",
        type=Path,
        default=Path("/scratch/mh7373/projects/cav-hoomd/aging_weak_lambda"),
    )
    parser.add_argument(
        "--lambdas",
        nargs="+",
        type=float,
        default=[0.0, 0.01, 0.016667, 0.023333, 0.03],
    )
    parser.add_argument("--n-replicas", type=int, default=DEFAULT_N_REPLICAS)
    parser.add_argument(
        "--switch-time-ps",
        type=float,
        default=DEFAULT_SWITCH_TIME_PS,
    )
    parser.add_argument("--cleanup", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    scans = scan_campaign(
        output_base=args.output_base,
        lambdas=args.lambdas,
        n_replicas=args.n_replicas,
        switch_time_ps=args.switch_time_ps,
    )

    summary = []
    for scan in scans:
        deleted: list[str] = []
        stray: list[str] = []
        if args.cleanup:
            deleted = [str(path) for path in cleanup_partial_replicas(scan, dry_run=args.dry_run)]
            stray = [
                str(path)
                for path in cleanup_stray_run_dirs(
                    args.output_base,
                    scan.lam,
                    switch_time_ps=args.switch_time_ps,
                    dry_run=args.dry_run,
                )
            ]

        summary.append(
            {
                "lam": scan.lam,
                "run_dir": str(scan.run_dir),
                "complete": scan.n_complete,
                "partial": scan.n_partial,
                "missing": scan.n_missing,
                "missing_replicas": list(scan.missing_replicas),
                "slurm_array": slurm_array_spec(list(scan.missing_replicas)),
                "deleted_files": deleted,
                "deleted_dirs": stray,
            }
        )

    if args.json:
        print(json.dumps(summary, indent=2))
    else:
        for item in summary:
            print(
                f"lam={item['lam']}: complete={item['complete']} "
                f"partial={item['partial']} missing={item['missing']}"
            )
            if item["slurm_array"]:
                print(f"  array={item['slurm_array'][:120]}{'...' if len(item['slurm_array']) > 120 else ''}")
            if args.cleanup:
                print(f"  deleted_files={len(item['deleted_files'])} deleted_dirs={len(item['deleted_dirs'])}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
