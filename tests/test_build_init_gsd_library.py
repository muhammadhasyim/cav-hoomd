"""Tests for building multi-replica init GSD from IC library trajectories."""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import gsd.hoomd
import numpy as np
import pytest

sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "aging_weak_lambda_ic_library"),
)

from build_init_gsd_library import (
    add_cavity_particle_at_origin,
    build_init_gsd,
    min_frames_required,
    source_mapping,
    validate_source_trajectories,
)


MOLECULAR_TEMPLATE = (
    Path(__file__).resolve().parents[1]
    / "aging_weak_lambda_ic_library"
    / "init_molecular_only.gsd"
)


def _write_molecular_trajectory(path: Path, n_frames: int, seed: int) -> None:
    """Write a minimal molecular GSD trajectory by cloning the library seed."""
    del seed
    with gsd.hoomd.open(MOLECULAR_TEMPLATE, "r") as template_traj:
        template = template_traj[0]
    with gsd.hoomd.open(path, "w") as traj:
        for frame_idx in range(n_frames):
            frame = copy.deepcopy(template)
            frame.log = {"Time/elapsed_ps": np.array([frame_idx * 200.0])}
            frame.particles.position[0, 0] += frame_idx * 1.0e-6
            traj.append(frame)


def test_source_mapping_round_robin() -> None:
    assert source_mapping(0, 8) == (0, 0)
    assert source_mapping(7, 8) == (7, 0)
    assert source_mapping(8, 8) == (0, 1)
    assert source_mapping(499, 8) == (3, 62)


def test_min_frames_required() -> None:
    assert min_frames_required(500, 8) == 63
    assert min_frames_required(10, 2) == 5


def test_add_cavity_particle_at_origin() -> None:
    with gsd.hoomd.open(MOLECULAR_TEMPLATE, "r") as traj:
        snap = copy.deepcopy(traj[0])

    add_cavity_particle_at_origin(snap)

    assert snap.particles.N == 501
    assert list(snap.particles.types) == ["O", "N", "L"]
    assert snap.particles.typeid[-1] == 2
    np.testing.assert_allclose(snap.particles.position[-1], [0.0, 0.0, 0.0])


def test_validate_source_trajectories_detects_short_stream(tmp_path: Path) -> None:
    base = tmp_path / "ic_library"
    prod0 = base / "replicas/replica_0/no_cavity/prod-0.gsd"
    prod1 = base / "replicas/replica_1/no_cavity/prod-1.gsd"
    prod0.parent.mkdir(parents=True)
    prod1.parent.mkdir(parents=True)
    _write_molecular_trajectory(prod0, n_frames=5, seed=0)
    _write_molecular_trajectory(prod1, n_frames=3, seed=1)

    with pytest.raises(RuntimeError, match="replica_1"):
        validate_source_trajectories(base, n_source_replicas=2, n_output_replicas=10)


def test_build_init_gsd_round_robin(tmp_path: Path) -> None:
    base = tmp_path / "ic_library"
    n_source = 2
    n_output = 10
    n_frames = min_frames_required(n_output, n_source)
    for replica in range(n_source):
        prod = base / f"replicas/replica_{replica}/no_cavity/prod-{replica}.gsd"
        prod.parent.mkdir(parents=True)
        _write_molecular_trajectory(prod, n_frames=n_frames, seed=replica)

    output = tmp_path / "init-test.gsd"
    manifest = build_init_gsd(
        base=base,
        output_path=output,
        n_output_replicas=n_output,
        n_source_replicas=n_source,
        add_cavity=True,
    )

    assert output.is_file()
    assert manifest["n_output_replicas"] == n_output
    with gsd.hoomd.open(output, "r") as traj:
        assert len(traj) == n_output
        assert traj[0].particles.N == 501
        assert list(traj[0].particles.types) == ["O", "N", "L"]

        src_replica, src_frame = source_mapping(5, n_source)
        with gsd.hoomd.open(
            base / f"replicas/replica_{src_replica}/no_cavity/prod-{src_replica}.gsd",
            "r",
        ) as src:
            expected_pos = src[src_frame].particles.position[0]
        np.testing.assert_allclose(traj[5].particles.position[0], expected_pos)
