"""Utilities for validating Maxwell-Boltzmann velocity distributions."""

from __future__ import annotations

from pathlib import Path
from typing import Tuple

import numpy as np
import scipy.stats


def compute_single_molecule_coupling(
    collective_coupling: float,
    num_molecules: int,
) -> float:
    """
    Compute single-molecule coupling lambda for fixed collective coupling.

    The collective coupling scales as g = lambda * sqrt(N), so:
    lambda = g / sqrt(N).

    Parameters
    ----------
    collective_coupling : float
        Target collective coupling g in atomic units.
    num_molecules : int
        Number of molecules N (must be positive).

    Returns
    -------
    float
        Single-molecule coupling lambda in atomic units.

    Raises
    ------
    ValueError
        If num_molecules is not positive or collective_coupling is negative.
    """
    if num_molecules <= 0:
        raise ValueError("num_molecules must be positive")
    if collective_coupling < 0.0:
        raise ValueError("collective_coupling must be non-negative")
    return float(collective_coupling) / np.sqrt(float(num_molecules))


def expected_velocity_component_variance(kT: float, mass: float) -> float:
    """
    Compute the expected variance for a single velocity component.

    For Maxwell-Boltzmann statistics, each Cartesian component satisfies:
    <v_i^2> = kT / m.

    Parameters
    ----------
    kT : float
        Thermal energy in Hartree (atomic units).
    mass : float
        Particle mass in atomic units (electron mass units).

    Returns
    -------
    float
        Expected variance <v_i^2> in (Bohr / atomic time)^2.

    Raises
    ------
    ValueError
        If kT is not positive or mass is not positive.
    """
    if kT <= 0.0:
        raise ValueError("kT must be positive")
    if mass <= 0.0:
        raise ValueError("mass must be positive")
    return float(kT) / float(mass)


def ks_test_velocity_component(
    velocity_components: np.ndarray,
    mass: float,
    kT: float,
) -> Tuple[float, float]:
    """
    Perform a Kolmogorov-Smirnov test against a Maxwell-Boltzmann component.

    Parameters
    ----------
    velocity_components : np.ndarray
        Sampled velocity components (1D array) in atomic units.
    mass : float
        Particle mass in atomic units (electron mass units).
    kT : float
        Thermal energy in Hartree (atomic units).

    Returns
    -------
    Tuple[float, float]
        (KS statistic, p-value).

    Raises
    ------
    ValueError
        If inputs are invalid or the sample is empty.
    """
    if velocity_components is None or len(velocity_components) == 0:
        raise ValueError("velocity_components must be a non-empty array")
    variance = expected_velocity_component_variance(kT, mass)
    sigma = np.sqrt(variance)
    statistic, pvalue = scipy.stats.kstest(
        velocity_components,
        scipy.stats.norm(loc=0.0, scale=sigma).cdf,
    )
    return float(statistic), float(pvalue)


def resolve_gsd_path(input_path: str, base_dir: Path) -> str:
    """
    Resolve a GSD file path to an absolute path.

    Parameters
    ----------
    input_path : str
        Path to the GSD file (absolute or relative).
    base_dir : Path
        Base directory to resolve relative paths against.

    Returns
    -------
    str
        Absolute path to the GSD file.
    """
    input_candidate = Path(input_path)
    if input_candidate.is_absolute():
        return str(input_candidate)
    return str((base_dir / input_candidate).resolve())
