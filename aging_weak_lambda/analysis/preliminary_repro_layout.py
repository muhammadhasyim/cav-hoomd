"""Build preliminary reproduction layout from completed campaign replicas."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from examples.slurm.aging_campaign_config import (  # noqa: E402
    CAMPAIGN_FKT_MAX_REFS_PRIMARY,
    CAMPAIGN_FKT_MAX_REFS_TARGET,
    CAMPAIGN_PRIMARY_RUNTIME_PS,
)
from examples.slurm.aging_campaign_status import (  # noqa: E402
    COMPLETE_MIN_BYTES,
    DEFAULT_FKT_MAX_REFERENCES,
    DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    DEFAULT_RUNTIME_PS,
    DEFAULT_SWITCH_TIME_PS,
    is_complete_replica,
    run_dir,
)
from examples.slurm.run_packed_aging_task import CAMPAIGN_LAMBDAS  # noqa: E402

LAMBDA_EPSILON_TAGS: Mapping[float, str] = {
    0.0: "0epos00",
    0.01: "1eneg02",
    0.016667: "2eneg02",
    0.023333: "2eneg02_l2333",
    0.03: "3eneg02",
}

FKT_PROCESS_SCRIPT = (
    REPO_ROOT / "reproduction" / "figure2_fkt" / "process_fskt_only.py"
)


def staged_coupling_dir_name(epsilon_tag: str, switch_time_ps: float) -> str:
    """Return the paper-style staged directory name for one coupling."""
    return f"cavity_coupling_{epsilon_tag}_switch_{switch_time_ps}ps"


def completion_options(
    *,
    runtime_ps: float = DEFAULT_RUNTIME_PS,
    reference_interval_ps: float = DEFAULT_FKT_REFERENCE_INTERVAL_PS,
    max_references: int = DEFAULT_FKT_MAX_REFERENCES,
    min_bytes: int = COMPLETE_MIN_BYTES,
) -> dict[str, float | int]:
    """Return shared kwargs for :func:`is_complete_replica`."""
    return {
        "runtime_ps": runtime_ps,
        "reference_interval_ps": reference_interval_ps,
        "max_references": max_references,
        "min_bytes": min_bytes,
    }


def list_complete_replicas(
    run_directory: Path,
    lam: float,
    *,
    n_replicas: int,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    **completion_kwargs: float | int,
) -> list[int]:
    """Return replica IDs that pass the campaign completion policy."""
    options = {
        **completion_kwargs,
        "require_step_protocol": lam != 0.0,
        "expected_lam": lam if lam != 0.0 else None,
        "expected_switch_time_ps": switch_time_ps if lam != 0.0 else None,
    }
    complete: list[int] = []
    for replica in range(n_replicas):
        if is_complete_replica(run_directory, replica, **options):
            complete.append(replica)
    return complete


def list_union_complete_replicas(
    run_directory: Path,
    lam: float,
    *,
    n_replicas: int,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    primary_runtime_ps: float = CAMPAIGN_PRIMARY_RUNTIME_PS,
    primary_max_references: int = CAMPAIGN_FKT_MAX_REFS_PRIMARY,
    target_runtime_ps: float = 1999.0,
    target_max_references: int = CAMPAIGN_FKT_MAX_REFS_TARGET,
    target_max_lag_tolerance_ps: float = 5.0,
    min_bytes: int = COMPLETE_MIN_BYTES,
) -> list[int]:
    """Return replicas complete at primary (1600/8) and/or near-target (~2000/10).

    Production 2000 ps jobs currently end at ``t_end ≈ 1999`` ps, so the target
    tier uses ``target_runtime_ps=1999`` by default. Reference times also drift by
    ~1 ps per reference (ref 9 headers near 1805 ps vs policy 1800 ps), so the
    target tier allows ``target_max_lag_tolerance_ps=5`` by default.

    Near-target extension replicas are accepted for averaging even when the
    step-protocol provenance marker is missing; primary-complete replicas still
    require the marker when ``lam != 0``. Late-time averages then use every
    replica that has data, with fewer contributors at the longest lags.
    """
    primary = set(
        list_complete_replicas(
            run_directory,
            lam,
            n_replicas=n_replicas,
            switch_time_ps=switch_time_ps,
            runtime_ps=primary_runtime_ps,
            max_references=primary_max_references,
            min_bytes=min_bytes,
        )
    )
    # Target tier: scientific F(k,t)/HDF5 completeness only (no protocol marker).
    target: set[int] = set()
    for replica in range(n_replicas):
        if is_complete_replica(
            run_directory,
            replica,
            min_bytes=min_bytes,
            runtime_ps=target_runtime_ps,
            max_references=target_max_references,
            max_lag_tolerance_ps=target_max_lag_tolerance_ps,
            require_step_protocol=False,
        ):
            target.add(replica)
    return sorted(primary | target)


def symlink_complete_fkt_files(
    source_run_dir: Path,
    staging_dir: Path,
    replicas: Sequence[int],
) -> int:
    """Symlink replica F(k,t) files into a staging directory."""
    staging_dir.mkdir(parents=True, exist_ok=True)
    linked = 0
    for replica in replicas:
        pattern = f"prod-{replica}_fkt_ref_*.txt"
        for source_path in sorted(source_run_dir.glob(pattern)):
            destination = staging_dir / source_path.name
            if destination.exists() or destination.is_symlink():
                destination.unlink()
            destination.symlink_to(source_path.resolve())
            linked += 1
    return linked


def average_fkt_in_directory(
    staging_dir: Path,
    *,
    python_executable: Path,
    max_time: float = 2500.0,
    dt: float = 1.0,
    runtime_ps: float | None = None,
    lag_extent_mode: str = "max",
) -> None:
    """Average symlinked replica F(k,t) files into master files."""
    for path in staging_dir.glob("master_fskt_ref*.txt"):
        path.unlink()
    for path in staging_dir.glob("master_fskt_ref*_sample_counts.txt"):
        path.unlink()

    env = dict(os.environ)
    pythonpath = str(REPO_ROOT)
    existing = env.get("PYTHONPATH")
    if existing:
        pythonpath = f"{pythonpath}:{existing}"
    env["PYTHONPATH"] = pythonpath
    command = [
        str(python_executable),
        str(FKT_PROCESS_SCRIPT),
        "--exp_dir",
        str(staging_dir),
        "--fskt_dt",
        str(dt),
        "--max_time",
        str(max_time),
        "--lag-extent-mode",
        lag_extent_mode,
    ]
    if runtime_ps is not None:
        command.extend(["--runtime-ps", f"{runtime_ps:g}"])
    subprocess.run(
        command,
        check=True,
        env=env,
        cwd=str(REPO_ROOT),
    )

    counts_dir = staging_dir / "fkt_sample_counts"
    counts_dir.mkdir(exist_ok=True)
    for path in staging_dir.glob("master_fskt_ref*_sample_counts.txt"):
        destination = counts_dir / path.name
        if destination.exists():
            destination.unlink()
        path.rename(destination)


def build_preliminary_layout(
    output_base: Path,
    staged_root: Path,
    *,
    python_executable: Path,
    n_replicas: int = 500,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    max_time: float = 2500.0,
    dt: float = 1.0,
    dry_run: bool = False,
    runtime_ps: float | None = None,
    max_references: int | None = None,
    use_all_available: bool = False,
) -> dict[float, dict[str, int | str]]:
    """Stage and average F(k,t) masters for all complete replicas per lambda.

    When ``use_all_available`` is True, include both primary-complete (1600 ps /
    8 refs) and near-target-complete (~1999 ps / 10 refs) replicas, average with
    per-lag sample counts, and do not clip lags to the primary runtime.
    """
    completion_kwargs = completion_options(
        **{
            key: value
            for key, value in {
                "runtime_ps": runtime_ps,
                "max_references": max_references,
            }.items()
            if value is not None
        }
    )
    summary: dict[float, dict[str, int | str]] = {}

    if staged_root.exists() and not dry_run:
        for child in staged_root.iterdir():
            if child.is_dir() and child.name.startswith("cavity_coupling_"):
                shutil.rmtree(child)

    for lam in CAMPAIGN_LAMBDAS:
        epsilon_tag = LAMBDA_EPSILON_TAGS[lam]
        source_run_dir = run_dir(output_base, lam, switch_time_ps=switch_time_ps)
        staging_dir = staged_root / staged_coupling_dir_name(
            epsilon_tag, switch_time_ps
        )
        if not source_run_dir.is_dir():
            summary[lam] = {
                "staged_dir": staging_dir.name,
                "complete_replicas": 0,
                "linked_fkt_files": 0,
            }
            continue

        if use_all_available:
            replicas = list_union_complete_replicas(
                source_run_dir,
                lam,
                n_replicas=n_replicas,
                switch_time_ps=switch_time_ps,
            )
        else:
            replicas = list_complete_replicas(
                source_run_dir,
                lam,
                n_replicas=n_replicas,
                switch_time_ps=switch_time_ps,
                **completion_kwargs,
            )
        if dry_run:
            refs_estimate = (
                CAMPAIGN_FKT_MAX_REFS_TARGET
                if use_all_available
                else int(completion_kwargs.get("max_references", DEFAULT_FKT_MAX_REFERENCES))
            )
            summary[lam] = {
                "staged_dir": staging_dir.name,
                "complete_replicas": len(replicas),
                "linked_fkt_files": len(replicas) * refs_estimate,
            }
            continue

        linked = symlink_complete_fkt_files(source_run_dir, staging_dir, replicas)
        if replicas:
            # Union mode: do not clip to primary runtime; each file contributes
            # only where it has data (lower N at late lags is expected).
            average_runtime_ps: float | None
            if use_all_available:
                average_runtime_ps = None
            elif runtime_ps is not None:
                average_runtime_ps = float(runtime_ps)
            else:
                average_runtime_ps = float(completion_kwargs["runtime_ps"])
            average_fkt_in_directory(
                staging_dir,
                python_executable=python_executable,
                max_time=max_time,
                dt=dt,
                runtime_ps=average_runtime_ps,
                # Always union lag coverage: a few short replicas must not
                # truncate the master grid (historical bug for λ=0 / 0.01).
                lag_extent_mode="max",
            )
        summary[lam] = {
            "staged_dir": staging_dir.name,
            "complete_replicas": len(replicas),
            "linked_fkt_files": linked,
        }

    return summary


def average_campaign_energies(
    output_base: Path,
    analysis_dir: Path,
    *,
    n_replicas: int = 500,
    switch_time_ps: float = DEFAULT_SWITCH_TIME_PS,
    stride: int = 1,
    structural_energy: str = "lj_only",
    runtime_ps: float | None = None,
    max_references: int | None = None,
    use_all_available: bool = False,
) -> dict[float, dict[str, np.ndarray | int]]:
    """Average observables HDF5 energies over scientifically complete replicas.

    Default ``lj_only`` stores LJ in the ``lj_coul_ha`` CSV column (RF Coulomb
    omitted) for fictive mapping against the ``lj_hartree`` column of
    ``pe_vs_T_calib_rf/potential_energy_vs_T_rf.txt``. Use ``lj_coul`` to
    average LJ + ``ewald_short`` instead.

    When ``use_all_available`` is True, average the union of primary-complete and
    near-target-complete replicas onto the longest time grid.
    """
    campaign_root = Path(__file__).resolve().parents[1]
    if str(campaign_root) not in sys.path:
        sys.path.insert(0, str(campaign_root))

    from aging_weak_lambda.analysis.average_aging_energies import (  # noqa: WPS433
        average_trimmed_series,
        load_replica,
        save_csv,
    )

    analysis_dir.mkdir(parents=True, exist_ok=True)
    completion_kwargs = completion_options(
        **{
            key: value
            for key, value in {
                "runtime_ps": runtime_ps,
                "max_references": max_references,
            }.items()
            if value is not None
        }
    )
    results: dict[float, dict[str, np.ndarray | int]] = {}

    for lam in CAMPAIGN_LAMBDAS:
        run_directory = run_dir(output_base, lam, switch_time_ps=switch_time_ps)
        if not run_directory.is_dir():
            continue

        if use_all_available:
            replicas = list_union_complete_replicas(
                run_directory,
                lam,
                n_replicas=n_replicas,
                switch_time_ps=switch_time_ps,
            )
        else:
            replicas = list_complete_replicas(
                run_directory,
                lam,
                n_replicas=n_replicas,
                switch_time_ps=switch_time_ps,
                **completion_kwargs,
            )
        loaded: list[tuple[np.ndarray, ...]] = []
        n_shifted = 0
        for replica in replicas:
            replica_path = run_directory / f"observables_replica_{replica}.h5"
            result = load_replica(
                replica_path,
                stride=stride,
                structural_energy=structural_energy,  # type: ignore[arg-type]
            )
            if result is None:
                continue
            series, shifted = result
            loaded.append(series)
            n_shifted += int(shifted)

        averaged = average_trimmed_series(loaded)
        if averaged is None:
            continue
        averaged["n_shifted"] = n_shifted
        save_csv(analysis_dir, lam, averaged)
        results[lam] = averaged

    return results
