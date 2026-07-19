"""Utilities for weak-coupling aging campaign resume and cleanup."""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from examples.slurm.aging_campaign_config import (
    CAMPAIGN_EXTENSION_RUNTIME_PS,
    CAMPAIGN_FKT_MAX_REFS_PRIMARY,
    CAMPAIGN_FKT_MAX_REFS_TARGET,
    CAMPAIGN_FKT_REFERENCE_INTERVAL_PS,
    CAMPAIGN_PRIMARY_RUNTIME_PS,
    CAMPAIGN_SWITCH_TIME_PS,
    CAMPAIGN_TARGET_RUNTIME_PS,
    CHECKPOINT_GSD_TEMPLATE,
    FKT_STATE_TEMPLATE,
)

# Target runtime at 1 ps HDF5 sampling; guard against truncated files.
COMPLETE_MIN_BYTES = 5_000_000
DEFAULT_SWITCH_TIME_PS = CAMPAIGN_SWITCH_TIME_PS
DEFAULT_N_REPLICAS = 1000
DEFAULT_RUNTIME_PS = CAMPAIGN_TARGET_RUNTIME_PS
DEFAULT_FKT_REFERENCE_INTERVAL_PS = CAMPAIGN_FKT_REFERENCE_INTERVAL_PS
DEFAULT_FKT_MAX_REFERENCES = CAMPAIGN_FKT_MAX_REFS_TARGET
DEFAULT_FKT_MAX_LAG_TOLERANCE_PS = 4.0
DEFAULT_PRIMARY_RUNTIME_PS = CAMPAIGN_PRIMARY_RUNTIME_PS
DEFAULT_EXTENSION_RUNTIME_PS = CAMPAIGN_EXTENSION_RUNTIME_PS
REPLICA_H5_PATTERN = re.compile(r"observables_replica_(\d+)\.h5$")
REFERENCE_TIME_PATTERN = re.compile(
    r"Reference time:\s*([0-9.eE+-]+)\s*ps"
)
PROVENANCE_PATH_KEYS = (
    "python_executable",
    "simulation_script",
    "cavitymd_core",
)
PROVENANCE_HASH_KEYS = (
    "simulation_script_sha256",
    "cavitymd_core_sha256",
)


@dataclass(frozen=True)
class LambdaScanResult:
    """Replica inventory for one lambda value."""

    lam: float
    run_dir: Path
    complete_replicas: tuple[int, ...]
    partial_replicas: tuple[int, ...]
    missing_replicas: tuple[int, ...]
    extension_replicas: tuple[int, ...] = ()
    preserved_primary_replicas: tuple[int, ...] = ()

    @property
    def n_complete(self) -> int:
        return len(self.complete_replicas)

    @property
    def n_partial(self) -> int:
        return len(self.partial_replicas)

    @property
    def n_missing(self) -> int:
        return len(self.missing_replicas)

    @property
    def fresh_run_replicas(self) -> tuple[int, ...]:
        """Replicas that still need a full target-runtime production run."""
        return tuple(
            replica
            for replica in (*self.partial_replicas, *self.missing_replicas)
            if replica not in self.extension_replicas
            and replica not in self.preserved_primary_replicas
        )


def lambda_to_tag(lam: float) -> str:
    """Convert a lambda value to the directory tag used by the campaign."""
    if lam == 0.0:
        return "0"
    return str(lam).replace(".", "p")


def _format_coupling_strength_for_runner(
    coupling: float,
    num_digits: int = 6,
) -> str:
    """Mirror ``cavitymd.utils.format_coupling_strength`` for directory names."""
    if coupling == 0.0:
        mantissa_str = "0" * num_digits
        return f"coupling_{mantissa_str}e+00"

    sign_prefix = ""
    abs_coupling = abs(coupling)
    if coupling < 0:
        sign_prefix = "neg"

    scientific_str = f"{abs_coupling:.{num_digits - 1}e}"
    mantissa_part, exp_part = scientific_str.split("e")
    mantissa = float(mantissa_part)
    exponent = int(exp_part)
    mantissa_scaled = int(mantissa * (10 ** (num_digits - 1)))
    mantissa_str = f"{mantissa_scaled:0{num_digits}d}"
    adjusted_exponent = exponent - (num_digits - 1)
    exp_str = f"{adjusted_exponent:+03d}"
    return f"coupling_{sign_prefix}{mantissa_str}e{exp_str}"


def coupling_dir_name(coupling: float, switch_time_ps: float) -> str:
    """Match the directory naming logic in examples/05_advanced_run.py."""
    return (
        f"{_format_coupling_strength_for_runner(coupling)}"
        f"_switch_{switch_time_ps}ps"
    )


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


def replica_checkpoint_path(run_directory: Path, replica: int) -> Path:
    """Return the single-frame checkpoint GSD used for runtime extension."""
    return run_directory / CHECKPOINT_GSD_TEMPLATE.format(replica=replica)


def replica_fkt_state_path(
    run_directory: Path,
    replica: int,
    *,
    job_name: str = "prod",
) -> Path:
    """Return the serialized F(k,t) reference state for resume."""
    return run_directory / FKT_STATE_TEMPLATE.format(job_name=job_name, replica=replica)


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


def protocol_marker_path(run_directory: Path, replica: int) -> Path:
    """Return the verified coupling-protocol marker path for one replica."""
    return run_directory / f"protocol_replica_{replica}.json"


def write_protocol_marker(
    run_directory: Path,
    replica: int,
    *,
    lam: float,
    switch_time_ps: float,
    seed: int,
    provenance: Mapping[str, str],
) -> Path:
    """Atomically record the step-coupling protocol used by a replica.

    Parameters
    ----------
    run_directory
        Coupling output directory.
    replica
        Replica identifier.
    lam
        Dimensionless target coupling.
    switch_time_ps
        Step activation time in picoseconds.
    seed
        Deterministic HOOMD random seed.
    provenance
        Executable/module paths and SHA-256 hashes used for the simulation.

    Returns
    -------
    pathlib.Path
        Final marker path.
    """
    if (
        not math.isfinite(lam)
        or lam < 0.0
        or not math.isfinite(switch_time_ps)
        or switch_time_ps < 0.0
        or seed <= 0
        or not _valid_provenance(provenance)
    ):
        raise ValueError("protocol marker values are not physically valid")
    run_directory.mkdir(parents=True, exist_ok=True)
    path = protocol_marker_path(run_directory, replica)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    payload = {
        "format_version": 1,
        "coupling_variant_type": "step",
        "lambda_coupling": lam,
        "switch_time_ps": switch_time_ps,
        "seed": seed,
        "provenance": dict(provenance),
    }
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)
    return path


def validate_protocol_marker(
    run_directory: Path,
    replica: int,
    *,
    expected_lam: float,
    expected_switch_time_ps: float,
    expected_seed: int,
) -> bool:
    """Validate a replica's recorded step-coupling protocol.

    Returns ``False`` for absent, malformed, constant-coupling, or mismatched
    markers. The marker supplies protocol provenance not present in legacy HDF5
    files.
    """
    if (
        not math.isfinite(expected_lam)
        or expected_lam < 0.0
        or not math.isfinite(expected_switch_time_ps)
        or expected_switch_time_ps < 0.0
        or expected_seed <= 0
    ):
        return False
    try:
        payload = json.loads(
            protocol_marker_path(run_directory, replica).read_text(
                encoding="utf-8"
            )
        )
        seed_value = payload.get("seed") if isinstance(payload, dict) else None
        provenance = (
            payload.get("provenance") if isinstance(payload, dict) else None
        )
        return (
            isinstance(payload, dict)
            and payload.get("format_version") == 1
            and payload.get("coupling_variant_type") == "step"
            and math.isclose(
                float(payload.get("lambda_coupling")),
                expected_lam,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            and math.isclose(
                float(payload.get("switch_time_ps")),
                expected_switch_time_ps,
                rel_tol=0.0,
                abs_tol=1.0e-9,
            )
            and isinstance(seed_value, int)
            and not isinstance(seed_value, bool)
            and seed_value == expected_seed
            and isinstance(provenance, dict)
            and _valid_provenance(provenance)
        )
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        return False


def _valid_provenance(provenance: Mapping[str, object]) -> bool:
    """Return whether executable provenance is complete and well formed."""
    for key in PROVENANCE_PATH_KEYS:
        value = provenance.get(key)
        if (
            not isinstance(value, str)
            or not value
            or not Path(value).is_absolute()
        ):
            return False
    for key in PROVENANCE_HASH_KEYS:
        value = provenance.get(key)
        if (
            not isinstance(value, str)
            or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)
        ):
            return False
    return True


def validate_fkt_file(
    path: Path,
    runtime_ps: float,
    expected_reference_time_ps: float,
    max_lag_tolerance_ps: float = DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
) -> bool:
    """Validate one two-column F(k,t) output file.

    The check rejects malformed, non-finite, duplicate, and out-of-order rows.
    It also requires a reference-time header and final lag coverage consistent
    with the campaign aging protocol runtime.

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

    if reference_time_ps is None:
        return False
    expected_max_lag_ps = runtime_ps - reference_time_ps
    if expected_max_lag_ps < 0.0:
        return False
    min_rows = (
        1
        if (runtime_ps - expected_reference_time_ps) <= max_lag_tolerance_ps
        else 2
    )
    if row_count < min_rows or previous_lag is None:
        return False
    if (
        abs(reference_time_ps - expected_reference_time_ps)
        > max_lag_tolerance_ps
    ):
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
    require_step_protocol: bool = False,
    expected_lam: float | None = None,
    expected_switch_time_ps: float | None = None,
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
    require_step_protocol
        Require a verified step-coupling provenance marker.
    expected_lam
        Expected dimensionless coupling when protocol validation is required.
    expected_switch_time_ps
        Expected coupling activation time in picoseconds.

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
    if require_step_protocol:
        if expected_lam is None or expected_switch_time_ps is None:
            return False
        if not validate_protocol_marker(
            run_directory,
            replica,
            expected_lam=expected_lam,
            expected_switch_time_ps=expected_switch_time_ps,
            expected_seed=replica + 1,
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


def is_complete_at_runtime(
    run_directory: Path,
    replica: int,
    runtime_ps: float,
    *,
    min_bytes: int = COMPLETE_MIN_BYTES,
    reference_interval_ps: float = DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    max_references: int | None = None,
    max_lag_tolerance_ps: float = DEFAULT_FKT_MAX_LAG_TOLERANCE_PS,
    require_step_protocol: bool = False,
    expected_lam: float | None = None,
    expected_switch_time_ps: float | None = None,
) -> bool:
    """Return whether a replica satisfies completion at an explicit runtime."""
    if max_references is None:
        max_references = int(math.floor(runtime_ps / reference_interval_ps)) + 1
    return is_complete_replica(
        run_directory,
        replica,
        min_bytes=min_bytes,
        runtime_ps=runtime_ps,
        reference_interval_ps=reference_interval_ps,
        max_references=max_references,
        max_lag_tolerance_ps=max_lag_tolerance_ps,
        require_step_protocol=require_step_protocol,
        expected_lam=expected_lam,
        expected_switch_time_ps=expected_switch_time_ps,
    )


def _runtime_completion_kwargs(
    **completion_kwargs: object,
) -> dict[str, object]:
    """Drop runtime aliases so explicit runtime positional args are not duplicated."""
    return {
        key: value
        for key, value in completion_kwargs.items()
        if key not in {"runtime_ps", "target_runtime_ps"}
    }


def _primary_completion_kwargs(
    **completion_kwargs: object,
) -> dict[str, object]:
    """Kwargs for primary-runtime checks with an explicit ``max_references``."""
    return {
        key: value
        for key, value in _runtime_completion_kwargs(**completion_kwargs).items()
        if key != "max_references"
    }


def needs_extension_replica(
    run_directory: Path,
    replica: int,
    *,
    primary_runtime_ps: float = DEFAULT_PRIMARY_RUNTIME_PS,
    target_runtime_ps: float = DEFAULT_RUNTIME_PS,
    **completion_kwargs: object,
) -> bool:
    """Return True when a replica finished the primary run and can extend."""
    target_kwargs = _runtime_completion_kwargs(**completion_kwargs)
    primary_kwargs = _primary_completion_kwargs(**completion_kwargs)
    checkpoint = replica_checkpoint_path(run_directory, replica)
    if not checkpoint.is_file():
        return False
    if is_complete_at_runtime(
        run_directory,
        replica,
        target_runtime_ps,
        **target_kwargs,  # type: ignore[arg-type]
    ):
        return False
    return is_complete_at_runtime(
        run_directory,
        replica,
        primary_runtime_ps,
        max_references=CAMPAIGN_FKT_MAX_REFS_PRIMARY,
        **primary_kwargs,  # type: ignore[arg-type]
    )


def replica_needs_fresh_run(
    run_directory: Path,
    replica: int,
    *,
    target_runtime_ps: float = DEFAULT_RUNTIME_PS,
    **completion_kwargs: object,
) -> bool:
    """Return True when a replica still needs a full target-runtime production run."""
    target_kwargs = _runtime_completion_kwargs(**completion_kwargs)
    primary_kwargs = _primary_completion_kwargs(**completion_kwargs)
    if is_complete_at_runtime(
        run_directory,
        replica,
        target_runtime_ps,
        **target_kwargs,  # type: ignore[arg-type]
    ):
        return False
    if needs_extension_replica(
        run_directory,
        replica,
        target_runtime_ps=target_runtime_ps,
        **completion_kwargs,  # type: ignore[arg-type]
    ):
        return False
    if is_complete_at_runtime(
        run_directory,
        replica,
        DEFAULT_PRIMARY_RUNTIME_PS,
        max_references=CAMPAIGN_FKT_MAX_REFS_PRIMARY,
        **primary_kwargs,  # type: ignore[arg-type]
    ):
        # Primary complete but no checkpoint: cannot extend safely.
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
        replica_checkpoint_path(run_directory, replica),
        replica_fkt_state_path(run_directory, replica),
        protocol_marker_path(run_directory, replica),
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
    require_step_protocol: bool = False,
    expected_lam: float | None = None,
    expected_switch_time_ps: float | None = None,
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
                    require_step_protocol=require_step_protocol,
                    expected_lam=expected_lam,
                    expected_switch_time_ps=expected_switch_time_ps,
                )
            if complete_replica:
                complete.append(replica)
            else:
                partial.append(replica)

    complete.sort()
    partial.sort()
    missing = [replica for replica in range(n_replicas) if replica not in complete]
    extension: list[int] = []
    preserved_primary: list[int] = []
    completion_kwargs = {
        "min_bytes": min_bytes,
        "runtime_ps": runtime_ps,
        "reference_interval_ps": reference_interval_ps,
        "max_references": max_references,
        "max_lag_tolerance_ps": max_lag_tolerance_ps,
        "require_step_protocol": require_step_protocol,
        "expected_lam": expected_lam,
        "expected_switch_time_ps": expected_switch_time_ps,
    }
    status_kwargs = {
        "min_bytes": min_bytes,
        "reference_interval_ps": reference_interval_ps,
        "max_lag_tolerance_ps": max_lag_tolerance_ps,
        "require_step_protocol": require_step_protocol,
        "expected_lam": expected_lam,
        "expected_switch_time_ps": expected_switch_time_ps,
    }
    for replica in range(n_replicas):
        if needs_extension_replica(directory, replica, **completion_kwargs):
            extension.append(replica)
        elif is_complete_at_runtime(
            directory,
            replica,
            DEFAULT_PRIMARY_RUNTIME_PS,
            max_references=CAMPAIGN_FKT_MAX_REFS_PRIMARY,
            **status_kwargs,
        ):
            preserved_primary.append(replica)

    return LambdaScanResult(
        lam=lam,
        run_dir=directory,
        complete_replicas=tuple(complete),
        partial_replicas=tuple(partial),
        missing_replicas=tuple(missing),
        extension_replicas=tuple(sorted(extension)),
        preserved_primary_replicas=tuple(sorted(preserved_primary)),
    )


def cleanup_partial_replicas(
    scan: LambdaScanResult,
    dry_run: bool = False,
) -> list[Path]:
    """Delete partial/failed replica artifacts while preserving complete outputs."""
    deleted: list[Path] = []

    for replica in scan.partial_replicas:
        if replica in scan.extension_replicas:
            continue
        if replica in scan.preserved_primary_replicas:
            continue
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


def cleanup_stray_gsd_files(
    run_directory: Path,
    *,
    dry_run: bool = False,
) -> list[Path]:
    """Delete ``prod-*.gsd`` trajectory files without touching HDF5 or F(k,t).

    Parameters
    ----------
    run_directory
        Coupling output directory.
    dry_run
        When ``True``, report paths without deleting them.

    Returns
    -------
    list[pathlib.Path]
        Existing GSD files selected for deletion.
    """
    if not run_directory.is_dir():
        return []
    paths = sorted(run_directory.glob("prod-*.gsd"))
    if not dry_run:
        for path in paths:
            path.unlink(missing_ok=True)
    return paths


def cleanup_stale_replica_locks(
    run_directory: Path,
    *,
    dry_run: bool = False,
) -> list[Path]:
    """Remove leftover ``.replica_*.lock`` files from crashed packed tasks.

    Parameters
    ----------
    run_directory
        Coupling output directory.
    dry_run
        When ``True``, report paths without deleting them.

    Returns
    -------
    list[pathlib.Path]
        Lock files selected for deletion.
    """
    if not run_directory.is_dir():
        return []
    paths = sorted(run_directory.glob(".replica_*.lock"))
    if not dry_run:
        for path in paths:
            path.unlink(missing_ok=True)
    return paths


def _hdf5_needs_clear(path: Path) -> bool:
    """Return True when a read-only open fails due to an inconsistent HDF5 file."""
    try:
        import h5py
    except ImportError:
        return False
    try:
        with h5py.File(path, "r"):
            return False
    except OSError as exc:
        message = str(exc).lower()
        return (
            "already open for write" in message
            or "file consistency" in message
            or "file is already open" in message
            or "unable to synchronously open" in message
        )


def _resolve_h5clear() -> str:
    """Return the ``h5clear`` executable path from PATH or the active env bin."""
    found = shutil.which("h5clear")
    if found:
        return found
    candidate = Path(sys.executable).resolve().parent / "h5clear"
    if candidate.is_file():
        return str(candidate)
    raise FileNotFoundError(
        "h5clear not found on PATH or next to the active Python executable"
    )


def clear_hdf5_consistency_locks(
    run_directory: Path,
    *,
    dry_run: bool = False,
) -> list[Path]:
    """Run ``h5clear -s`` on HDF5 files that fail a read-only open.

    Parameters
    ----------
    run_directory
        Coupling output directory containing ``observables_replica_*.h5``.
    dry_run
        When ``True``, report candidates without clearing them.

    Returns
    -------
    list[pathlib.Path]
        Files that needed (or would need) ``h5clear``.
    """
    if not run_directory.is_dir():
        return []
    cleared: list[Path] = []
    h5clear = None
    for path in sorted(run_directory.glob("observables_replica_*.h5")):
        if not _hdf5_needs_clear(path):
            continue
        cleared.append(path)
        if dry_run:
            continue
        if h5clear is None:
            h5clear = _resolve_h5clear()
        result = subprocess.run(
            [h5clear, "-s", str(path)],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"h5clear failed for {path}: {result.stderr.strip() or result.stdout.strip()}"
            )
    return cleared


def run_preflight_cleanup(
    output_base: Path,
    lambdas: list[float],
    *,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    dry_run: bool = False,
) -> list[dict[str, object]]:
    """Clear locked HDF5 files, stray GSDs, and stale replica locks.

    Does not delete scientifically complete HDF5 or F(k,t) data.
    """
    summaries: list[dict[str, object]] = []
    for lam in lambdas:
        directory = run_dir(output_base, lam, switch_time_ps)
        cleared = clear_hdf5_consistency_locks(directory, dry_run=dry_run)
        gsd_removed = cleanup_stray_gsd_files(directory, dry_run=dry_run)
        locks_removed = cleanup_stale_replica_locks(directory, dry_run=dry_run)
        summaries.append(
            {
                "lam": lam,
                "run_dir": str(directory),
                "h5_cleared": [str(path) for path in cleared],
                "gsd_removed": [str(path) for path in gsd_removed],
                "locks_removed": [str(path) for path in locks_removed],
            }
        )
    return summaries


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
    parser.add_argument(
        "--preflight",
        action="store_true",
        help=(
            "Clear locked HDF5 files (h5clear -s), delete prod-*.gsd only, "
            "and remove stale .replica_*.lock files without deleting science data"
        ),
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()

    if args.preflight:
        preflight = run_preflight_cleanup(
            output_base=args.output_base,
            lambdas=args.lambdas,
            switch_time_ps=args.switch_time_ps,
            dry_run=args.dry_run,
        )
        if args.json:
            print(json.dumps(preflight, indent=2))
        else:
            for item in preflight:
                print(
                    f"lam={item['lam']}: h5_cleared={len(item['h5_cleared'])} "
                    f"gsd_removed={len(item['gsd_removed'])} "
                    f"locks_removed={len(item['locks_removed'])}"
                )
        return 0

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
