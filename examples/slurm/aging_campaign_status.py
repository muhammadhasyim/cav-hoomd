"""Utilities for weak-coupling aging campaign resume and cleanup."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

COMPLETE_MIN_BYTES = 1_380_000_000
DEFAULT_SWITCH_TIME_PS = 200.0
DEFAULT_N_REPLICAS = 1000
DEFAULT_RUNTIME_PS = 2500.0
DEFAULT_FKT_REFERENCE_INTERVAL_PS = 200.0
DEFAULT_FKT_MAX_REFERENCES = 13
DEFAULT_FKT_MAX_LAG_TOLERANCE_PS = 2.0
REPLICA_H5_PATTERN = re.compile(r"observables_replica_(\d+)\.h5$")
REFERENCE_TIME_PATTERN = re.compile(
    r"Reference time:\s*([0-9.eE+-]+)\s*ps"
)


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


def fkt_file_path(
    run_directory: Path,
    replica: int,
    reference_index: int,
) -> Path:
    """Return one structural-correlation output path.

    Parameters
    ----------
    run_directory
        Directory containing output for one coupling value.
    replica
        Replica identifier in the campaign ID domain.
    reference_index
        Zero-based F(k,t) reference index.

    Returns
    -------
    pathlib.Path
        Path using the production ``prod-<replica>_fkt_ref_<index>.txt``
        naming convention.
    """
    return run_directory / f"prod-{replica}_fkt_ref_{reference_index:03d}.txt"


def validate_fkt_file(
    path: Path,
    runtime_ps: float,
    expected_reference_time_ps: float,
    max_lag_tolerance_ps: float = DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
) -> bool:
    """Validate one two-column F(k,t) output file.

    The check rejects malformed, non-finite, duplicate, and out-of-order rows.
    It also requires a reference-time header and final lag coverage consistent
    with the 2500 ps aging protocol.

    Parameters
    ----------
    path
        F(k,t) text file to validate.
    runtime_ps
        Requested simulation runtime in picoseconds.
    expected_reference_time_ps
        Reference time implied by the file's reference index, in picoseconds.
    max_lag_tolerance_ps
        Maximum absolute final-lag discrepancy in picoseconds.

    Returns
    -------
    bool
        ``True`` only when the complete file satisfies every check.

    Side Effects
    ------------
    The file is read but never modified.
    """
    if not path.is_file():
        return False
    if (
        not math.isfinite(runtime_ps)
        or not math.isfinite(expected_reference_time_ps)
        or not math.isfinite(max_lag_tolerance_ps)
        or runtime_ps <= 0.0
        or expected_reference_time_ps < 0.0
        or expected_reference_time_ps > runtime_ps
        or max_lag_tolerance_ps < 0.0
    ):
        return False

    reference_time_ps: float | None = None
    previous_lag: float | None = None
    row_count = 0

    try:
        with path.open("r", encoding="utf-8", errors="strict") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith("#"):
                    match = REFERENCE_TIME_PATTERN.search(stripped)
                    if match is not None:
                        reference_time = float(match.group(1))
                        if (
                            not math.isfinite(reference_time)
                            or reference_time < 0.0
                            or reference_time > runtime_ps
                        ):
                            return False
                        if (
                            reference_time_ps is not None
                            and reference_time != reference_time_ps
                        ):
                            return False
                        reference_time_ps = reference_time
                    continue

                fields = stripped.split()
                if len(fields) != 2:
                    return False
                lag_time_ps, correlation = (float(value) for value in fields)
                if not (
                    math.isfinite(lag_time_ps)
                    and math.isfinite(correlation)
                    and lag_time_ps >= 0.0
                ):
                    return False
                if previous_lag is not None and lag_time_ps <= previous_lag:
                    return False
                previous_lag = lag_time_ps
                row_count += 1
    except (OSError, UnicodeError, ValueError):
        return False

    if reference_time_ps is None or row_count < 2 or previous_lag is None:
        return False
    if (
        abs(reference_time_ps - expected_reference_time_ps)
        > max_lag_tolerance_ps
    ):
        return False
    expected_max_lag_ps = runtime_ps - reference_time_ps
    if expected_max_lag_ps < 0.0:
        return False
    return abs(previous_lag - expected_max_lag_ps) <= max_lag_tolerance_ps


def is_complete_replica(
    run_directory: Path,
    replica: int,
    *,
    min_bytes: int = COMPLETE_MIN_BYTES,
    runtime_ps: float = DEFAULT_RUNTIME_PS,
    reference_interval_ps: float = DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    max_references: int = DEFAULT_FKT_MAX_REFERENCES,
    max_lag_tolerance_ps: float = DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
) -> bool:
    """Return whether a replica is complete for relaxation analysis.

    Parameters
    ----------
    run_directory
        Directory containing one coupling value's output.
    replica
        Replica identifier.
    min_bytes
        Minimum observables HDF5 size in bytes.
    runtime_ps
        Requested simulation runtime in picoseconds.
    reference_interval_ps
        Time between F(k,t) references in picoseconds.
    max_references
        Required number of reference files.
    max_lag_tolerance_ps
        Allowed final-lag discrepancy in picoseconds.

    Returns
    -------
    bool
        ``True`` when the HDF5 guard and every F(k,t) file pass.
    """
    if (
        not math.isfinite(runtime_ps)
        or not math.isfinite(reference_interval_ps)
        or not math.isfinite(max_lag_tolerance_ps)
        or runtime_ps <= 0.0
        or reference_interval_ps <= 0.0
        or max_lag_tolerance_ps < 0.0
        or max_references <= 0
    ):
        return False
    if not is_complete_h5(
        replica_h5_path(run_directory, replica),
        min_bytes=min_bytes,
    ):
        return False

    for reference_index in range(max_references):
        reference_time_ps = reference_index * reference_interval_ps
        expected_max_lag_ps = runtime_ps - reference_time_ps
        if expected_max_lag_ps < 0.0:
            return False
        if not validate_fkt_file(
            fkt_file_path(run_directory, replica, reference_index),
            runtime_ps=runtime_ps,
            expected_reference_time_ps=reference_time_ps,
            max_lag_tolerance_ps=max_lag_tolerance_ps,
        ):
            return False
    return True


def cleanup_replica_artifacts(
    run_directory: Path,
    replica: int,
    *,
    dry_run: bool = False,
) -> list[Path]:
    """Remove restart-sensitive artifacts for one invalid replica.

    Parameters
    ----------
    run_directory
        Directory containing one coupling value's output.
    replica
        Replica identifier to clean.
    dry_run
        When ``True``, report paths without deleting them.

    Returns
    -------
    list[pathlib.Path]
        Existing artifacts selected for deletion.

    Side Effects
    ------------
    Unless ``dry_run`` is set, removes the replica HDF5, GSD, and all of its
    F(k,t) reference files. Other replicas are untouched.
    """
    paths = [
        replica_h5_path(run_directory, replica),
        replica_gsd_path(run_directory, replica),
        *sorted(run_directory.glob(f"prod-{replica}_fkt_ref_*.txt")),
    ]
    existing_paths = [path for path in paths if path.exists()]
    if not dry_run:
        for path in existing_paths:
            path.unlink()
    return existing_paths


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
    require_fkt: bool = False,
    runtime_ps: float = DEFAULT_RUNTIME_PS,
    reference_interval_ps: float = DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    max_references: int = DEFAULT_FKT_MAX_REFERENCES,
    max_lag_tolerance_ps: float = DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
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
            if not 0 <= replica < n_replicas:
                continue
            if path != replica_h5_path(directory, replica):
                continue
            complete_replica = is_complete_h5(path, min_bytes=min_bytes)
            if require_fkt:
                complete_replica = is_complete_replica(
                    directory,
                    replica,
                    min_bytes=min_bytes,
                    runtime_ps=runtime_ps,
                    reference_interval_ps=reference_interval_ps,
                    max_references=max_references,
                    max_lag_tolerance_ps=max_lag_tolerance_ps,
                )
            if complete_replica:
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
        deleted.extend(
            cleanup_replica_artifacts(
                scan.run_dir,
                replica,
                dry_run=dry_run,
            )
        )

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
