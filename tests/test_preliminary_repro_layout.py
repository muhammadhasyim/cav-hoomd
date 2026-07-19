"""Tests for preliminary reproduction layout helpers."""

from __future__ import annotations

from pathlib import Path

import inspect

from aging_weak_lambda.analysis.preliminary_repro_layout import (
    average_campaign_energies,
    list_complete_replicas,
    staged_coupling_dir_name,
    symlink_complete_fkt_files,
)
from examples.slurm.aging_campaign_status import run_dir, write_protocol_marker
from tests.test_aging_campaign_status import TEST_PROVENANCE, _write_complete_replica


def test_staged_coupling_dir_name_matches_repro_profile() -> None:
    assert staged_coupling_dir_name("1eneg02", 200.0) == (
        "cavity_coupling_1eneg02_switch_200.0ps"
    )


def test_average_campaign_energies_defaults_to_lj_only() -> None:
    """RF Coulomb is omitted; structural CSV column carries LJ only."""
    default = inspect.signature(average_campaign_energies).parameters[
        "structural_energy"
    ].default
    assert default == "lj_only"


def test_list_complete_replicas_respects_protocol_requirement(tmp_path: Path) -> None:
    lam = 0.01
    directory = run_dir(tmp_path, lam, switch_time_ps=200.0)
    directory.mkdir(parents=True)
    _write_complete_replica(directory, 4, min_bytes=1024)

    without_marker = list_complete_replicas(
        directory,
        lam,
        n_replicas=10,
        runtime_ps=14.0,
        reference_interval_ps=1.0,
        max_references=13,
        min_bytes=1024,
    )
    assert without_marker == []

    write_protocol_marker(
        directory,
        replica=4,
        lam=lam,
        switch_time_ps=200.0,
        seed=5,
        provenance=TEST_PROVENANCE,
    )
    with_marker = list_complete_replicas(
        directory,
        lam,
        n_replicas=10,
        runtime_ps=14.0,
        reference_interval_ps=1.0,
        max_references=13,
        min_bytes=1024,
    )
    assert with_marker == [4]


def test_symlink_complete_fkt_files_links_only_requested_replicas(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source"
    staging = tmp_path / "staging"
    source.mkdir()
    (source / "prod-1_fkt_ref_000.txt").write_text("# test\n0 1\n", encoding="utf-8")
    (source / "prod-2_fkt_ref_000.txt").write_text("# test\n0 1\n", encoding="utf-8")

    linked = symlink_complete_fkt_files(source, staging, [1])

    assert linked == 1
    assert (staging / "prod-1_fkt_ref_000.txt").is_symlink()
    assert not (staging / "prod-2_fkt_ref_000.txt").exists()
