"""Tests for aging campaign resume/cleanup helpers."""

from __future__ import annotations

from pathlib import Path

from examples.slurm.aging_campaign_status import (
    COMPLETE_MIN_BYTES,
    cleanup_replica_artifacts,
    cleanup_partial_replicas,
    cleanup_stray_run_dirs,
    coupling_dir_name,
    fkt_file_path,
    is_complete_h5,
    is_complete_replica,
    lambda_to_tag,
    run_dir,
    scan_lambda_run,
    slurm_array_spec,
    validate_fkt_file,
)


def _write_fkt_file(path: Path, reference_time_ps: float, max_lag_ps: float) -> None:
    """Write a minimal, valid two-column F(k,t) file."""
    path.write_text(
        "# F(k,t) correlation function\n"
        f"# Reference time: {reference_time_ps:.3f} ps\n"
        "# lag_time_ps\tF(k,t)\n"
        f"0.000000\t1.000000\n{max_lag_ps:.6f}\t0.100000\n",
        encoding="utf-8",
    )


def _write_complete_replica(
    directory: Path,
    replica: int,
    *,
    min_bytes: int,
    runtime_ps: float = 14.0,
    reference_interval_ps: float = 1.0,
    max_references: int = 13,
) -> None:
    """Write all artifacts required by the scientific completion policy."""
    (directory / f"observables_replica_{replica}.h5").write_bytes(
        b"x" * min_bytes
    )
    for ref_index in range(max_references):
        reference_time_ps = ref_index * reference_interval_ps
        _write_fkt_file(
            fkt_file_path(directory, replica, ref_index),
            reference_time_ps,
            runtime_ps - reference_time_ps,
        )


def test_lambda_to_tag() -> None:
    assert lambda_to_tag(0.0) == "0"
    assert lambda_to_tag(0.01) == "0p01"


def test_coupling_dir_name_matches_runner() -> None:
    assert coupling_dir_name(0.0, 200.0) == "cavity_coupling_0epos00_switch_200.0ps"
    assert coupling_dir_name(0.01, 200.0) == "cavity_coupling_1eneg02_switch_200.0ps"
    assert coupling_dir_name(0.03, 200.0) == "cavity_coupling_3eneg02_switch_200.0ps"


def test_run_dir() -> None:
    base = Path("/tmp/aging")
    assert run_dir(base, 0.01) == base / "lambda0p01" / "cavity_coupling_1eneg02_switch_200.0ps"


def test_is_complete_h5(tmp_path: Path) -> None:
    complete = tmp_path / "observables_replica_0.h5"
    partial = tmp_path / "observables_replica_1.h5"
    min_bytes = 1024
    complete.write_bytes(b"x" * min_bytes)
    partial.write_bytes(b"x" * 512)
    assert is_complete_h5(complete, min_bytes=min_bytes)
    assert not is_complete_h5(partial, min_bytes=min_bytes)


def test_validate_fkt_file_accepts_complete_ordered_data(tmp_path: Path) -> None:
    path = tmp_path / "prod-0_fkt_ref_000.txt"
    _write_fkt_file(path, reference_time_ps=0.0, max_lag_ps=10.0)

    assert validate_fkt_file(path, expected_max_lag_ps=10.0)


def test_validate_fkt_file_rejects_malformed_row(tmp_path: Path) -> None:
    path = tmp_path / "prod-0_fkt_ref_000.txt"
    path.write_text(
        "# Reference time: 0.0 ps\n0.0 1.0\nbroken\n10.0 0.1\n",
        encoding="utf-8",
    )

    assert not validate_fkt_file(path, expected_max_lag_ps=10.0)


def test_validate_fkt_file_rejects_duplicate_lag_time(tmp_path: Path) -> None:
    path = tmp_path / "prod-0_fkt_ref_000.txt"
    path.write_text(
        "# Reference time: 0.0 ps\n0.0 1.0\n1.0 0.5\n1.0 0.4\n10.0 0.1\n",
        encoding="utf-8",
    )

    assert not validate_fkt_file(path, expected_max_lag_ps=10.0)


def test_validate_fkt_file_rejects_out_of_order_lag_time(
    tmp_path: Path,
) -> None:
    path = tmp_path / "prod-0_fkt_ref_000.txt"
    path.write_text(
        "# Reference time: 0.0 ps\n0.0 1.0\n2.0 0.5\n1.0 0.4\n10.0 0.1\n",
        encoding="utf-8",
    )

    assert not validate_fkt_file(path, expected_max_lag_ps=10.0)


def test_validate_fkt_file_rejects_short_coverage(tmp_path: Path) -> None:
    path = tmp_path / "prod-0_fkt_ref_000.txt"
    _write_fkt_file(path, reference_time_ps=0.0, max_lag_ps=7.0)

    assert not validate_fkt_file(
        path,
        expected_max_lag_ps=10.0,
        max_lag_tolerance_ps=2.0,
    )


def test_is_complete_replica_requires_all_fkt_references(
    tmp_path: Path,
) -> None:
    min_bytes = 1024
    _write_complete_replica(tmp_path, 7, min_bytes=min_bytes)

    assert is_complete_replica(
        tmp_path,
        7,
        min_bytes=min_bytes,
        runtime_ps=14.0,
        reference_interval_ps=1.0,
        max_references=13,
    )

    fkt_file_path(tmp_path, 7, 12).unlink()
    assert not is_complete_replica(
        tmp_path,
        7,
        min_bytes=min_bytes,
        runtime_ps=14.0,
        reference_interval_ps=1.0,
        max_references=13,
    )


def test_slurm_array_spec() -> None:
    assert slurm_array_spec([0, 1, 2, 5, 7, 8]) == "0-2,5,7-8"
    assert slurm_array_spec([]) == ""


def test_scan_lambda_run(tmp_path: Path) -> None:
    lam = 0.01
    min_bytes = 1024
    directory = run_dir(tmp_path, lam)
    directory.mkdir(parents=True)

    complete_path = directory / "observables_replica_0.h5"
    partial_path = directory / "observables_replica_1.h5"
    complete_path.write_bytes(b"x" * min_bytes)
    partial_path.write_bytes(b"x" * 512)

    scan = scan_lambda_run(tmp_path, lam, n_replicas=4, min_bytes=min_bytes)
    assert scan.complete_replicas == (0,)
    assert scan.partial_replicas == (1,)
    assert scan.missing_replicas == (2, 3)


def test_cleanup_partial_replicas(tmp_path: Path) -> None:
    lam = 0.01
    directory = run_dir(tmp_path, lam)
    directory.mkdir(parents=True)

    partial_h5 = directory / "observables_replica_3.h5"
    partial_gsd = directory / "prod-3.gsd"
    partial_h5.write_bytes(b"x" * 512)
    partial_gsd.write_bytes(b"gsd")

    scan = scan_lambda_run(tmp_path, lam, n_replicas=5)
    deleted = cleanup_partial_replicas(scan, dry_run=False)

    assert partial_h5 in deleted
    assert partial_gsd in deleted
    assert not partial_h5.exists()
    assert not partial_gsd.exists()


def test_cleanup_replica_artifacts_removes_fkt_files(tmp_path: Path) -> None:
    _write_complete_replica(tmp_path, 4, min_bytes=1024)
    gsd_path = tmp_path / "prod-4.gsd"
    gsd_path.write_bytes(b"gsd")

    deleted = cleanup_replica_artifacts(tmp_path, 4)

    assert tmp_path / "observables_replica_4.h5" in deleted
    assert gsd_path in deleted
    assert fkt_file_path(tmp_path, 4, 0) in deleted
    assert not list(tmp_path.glob("prod-4_fkt_ref_*.txt"))


def test_cleanup_stray_run_dirs(tmp_path: Path) -> None:
    lam = 0.0
    expected = run_dir(tmp_path, lam)
    stray = expected.parent / "cavity_coupling_1eneg02_switch_200.0ps"
    expected.mkdir(parents=True)
    stray.mkdir(parents=True)
    (stray / "observables_replica_0.h5").write_bytes(b"x")

    removed = cleanup_stray_run_dirs(tmp_path, lam, dry_run=False)
    assert stray in removed
    assert not stray.exists()
    assert expected.exists()
