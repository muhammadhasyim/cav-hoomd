"""Tests for study configuration helpers."""

import importlib.util
from pathlib import Path

import pytest


def _load_module(module_path: Path):
    spec = importlib.util.spec_from_file_location(module_path.stem, module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


_STUDIES_MODULE = _load_module(
    Path(__file__).resolve().parents[1] / "src" / "cavitymd" / "studies.py"
)
compute_replica_count = _STUDIES_MODULE.compute_replica_count


def test_compute_replica_count_scaling():
    """Scale replicas inversely with system size."""
    replicas = compute_replica_count(base_molecules=250, base_replicas=500, target_molecules=10000)
    assert replicas == 13


def test_compute_replica_count_minimum():
    """Respect minimum replica count."""
    replicas = compute_replica_count(
        base_molecules=250,
        base_replicas=1,
        target_molecules=100000,
        min_replicas=5,
    )
    assert replicas == 5


def test_compute_replica_count_maximum():
    """Respect maximum replica count."""
    replicas = compute_replica_count(
        base_molecules=250,
        base_replicas=500,
        target_molecules=50,
        max_replicas=200,
    )
    assert replicas == 200


@pytest.mark.parametrize(
    "base_molecules, base_replicas, target_molecules",
    [
        (0, 500, 250),
        (250, 0, 250),
        (250, 500, 0),
    ],
)
def test_compute_replica_count_invalid_inputs(base_molecules, base_replicas, target_molecules):
    """Reject non-positive inputs."""
    with pytest.raises(ValueError):
        compute_replica_count(base_molecules, base_replicas, target_molecules)


def test_compute_replica_count_invalid_minimum():
    """Reject invalid minimum replica count."""
    with pytest.raises(ValueError):
        compute_replica_count(250, 500, 250, min_replicas=0)


def test_compute_replica_count_invalid_maximum():
    """Reject invalid maximum replica count."""
    with pytest.raises(ValueError):
        compute_replica_count(250, 500, 250, max_replicas=0)
