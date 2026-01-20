"""Study configuration helpers."""

from __future__ import annotations

import math
from typing import Optional


def compute_replica_count(
    base_molecules: int,
    base_replicas: int,
    target_molecules: int,
    min_replicas: int = 1,
    max_replicas: Optional[int] = None,
) -> int:
    """
    Compute replica count for finite-size scaling averages.

    The default scaling is inverse with system size so that the number of
    trajectories decreases as N increases while keeping statistical effort
    comparable to the base system.

    Parameters
    ----------
    base_molecules : int
        Reference number of molecules (e.g., 250).
    base_replicas : int
        Reference number of trajectories for base_molecules (e.g., 500).
    target_molecules : int
        Target system size for which to compute replicas.
    min_replicas : int, optional
        Minimum number of replicas to run (default: 1).
    max_replicas : Optional[int], optional
        Optional cap on replicas.

    Returns
    -------
    int
        Number of replicas to run for target_molecules.

    Raises
    ------
    ValueError
        If any input is not positive or min_replicas is invalid.
    """
    if base_molecules <= 0 or base_replicas <= 0 or target_molecules <= 0:
        raise ValueError("base_molecules, base_replicas, and target_molecules must be positive")
    if min_replicas <= 0:
        raise ValueError("min_replicas must be positive")
    if max_replicas is not None and max_replicas <= 0:
        raise ValueError("max_replicas must be positive when provided")

    scaling = base_replicas * (base_molecules / target_molecules)
    replicas = int(math.ceil(scaling))
    replicas = max(replicas, min_replicas)
    if max_replicas is not None:
        replicas = min(replicas, max_replicas)
    return replicas
