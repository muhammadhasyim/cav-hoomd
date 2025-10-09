# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Utility functions and constants for cavity molecular dynamics simulations."""

from typing import Union, List, Tuple
import numpy as np


class PhysicalConstants:
    r"""
    Physical constants and unit conversions for cavity molecular dynamics simulations.
    
    This class provides fundamental physical constants in atomic units and
    conversion functions between different unit systems commonly used in
    cavity-coupled molecular dynamics simulations.
    
    **Atomic Units Convention:**
    
    All calculations in cavity MD use atomic units where:
    - Length: 1 a.u. = 1 Bohr = 0.529177 Å
    - Energy: 1 a.u. = 1 Hartree = 27.2114 eV = 219474.63 cm⁻¹
    - Time: 1 a.u. = 2.418884 × 10⁻¹⁷ s = 0.02418884 ps
    - Mass: 1 a.u. = electron mass = 9.109 × 10⁻³¹ kg
    
    **Frequency Conversions:**
    
    Cavity frequencies are typically given in cm⁻¹ (wavenumbers) but must be
    converted to atomic units for simulation:
    
    .. math::
        
        \omega_{\text{a.u.}} = \frac{\tilde{\nu}_{\text{cm}^{-1}}}{219474.63}
    
    where :math:`\tilde{\nu}` is the wavenumber in cm⁻¹.
    
    **Thermostat Time Constants:**
    
    Thermostat damping coefficients relate to time constants via:
    
    .. math::
        
        \gamma = \frac{1}{\tau}
    
    where :math:`\tau` is the relaxation time constant and :math:`\gamma` is
    the damping coefficient used in Langevin dynamics.
    
    Attributes
    ----------
    HARTREE_TO_CM_MINUS1 : float
        Conversion factor from Hartree to cm⁻¹ (219474.63)
    KB_HARTREE_PER_K : float  
        Boltzmann constant in Hartree per Kelvin (3.167 × 10⁻⁶)
    ENERGY_JOULES : float
        Hartree to Joules conversion (4.35974 × 10⁻¹⁸)
    LENGTH_METERS : float
        Bohr to meters conversion (5.29177 × 10⁻¹¹)
    MASS_KG : float
        Electron mass in kg (9.10938 × 10⁻³¹)
    TIME_SECONDS : float
        Atomic time unit to seconds (2.418884 × 10⁻¹⁷)
    TIME_PS_CONVERSION : float
        Atomic time unit to picoseconds (2.418884 × 10⁻⁵)
    
    Examples
    --------
    **Frequency conversion:**
    
    >>> # Convert 2000 cm⁻¹ to atomic units
    >>> freq_au = 2000.0 / PhysicalConstants.HARTREE_TO_CM_MINUS1
    >>> print(f"2000 cm⁻¹ = {freq_au:.6f} a.u.")
    
    **Time conversions:**
    
    >>> # Convert 5 ps to atomic units
    >>> time_au = PhysicalConstants.ps_to_atomic_units(5.0)
    >>> print(f"5 ps = {time_au:.1f} a.u.")
    
    **Thermostat setup:**
    
    >>> # Calculate damping for 5 ps time constant
    >>> gamma = PhysicalConstants.gamma_from_tau_ps(5.0)
    >>> print(f"τ = 5 ps → γ = {gamma:.6f} a.u.⁻¹")
    
    Notes
    -----
    These constants are fundamental to ensuring proper unit consistency
    between HOOMD-blue (which uses reduced units) and the cavity-molecular coupling
    formulation (which uses atomic units).
    """
    
    # Energy conversions
    HARTREE_TO_CM_MINUS1: float = 219474.63
    """Conversion factor from Hartree to cm⁻¹ wavenumbers"""
    
    KB_HARTREE_PER_K: float = 3.167e-6  # Boltzmann constant in Hartree/K
    """Boltzmann constant in Hartree per Kelvin"""
    
    ENERGY_JOULES: float = 4.35974e-18  # Hartree to Joules
    """Hartree to Joules conversion factor"""
    
    # Length conversions  
    LENGTH_METERS: float = 5.29177210544e-11  # Bohr to meters
    """Bohr radius to meters conversion factor"""
    
    # Mass conversions
    MASS_KG: float = 9.1093837139e-31  # Electron mass in kg
    """Electron mass in kg"""
    
    AMU_TO_ELECTRON_MASS: float = 1822.888  # Atomic mass unit to electron mass
    """Conversion factor from AMU to electron masses (used in atomic units)"""
    
    # Time conversions
    TIME_SECONDS: float = 2.418884e-17  # Atomic time unit to seconds
    """Atomic time unit to seconds conversion factor"""
    
    TIME_PS_CONVERSION: float = 2.418884e-5  # a.u. to picoseconds (corrected from 0.02418884)
    """Atomic time unit to picoseconds conversion factor"""
    
    @classmethod
    def ps_to_atomic_units(cls, time_ps: float) -> float:
        """
        Convert time from picoseconds to atomic units.
        
        Parameters
        ----------
        time_ps : float
            Time in picoseconds
            
        Returns
        -------
        float
            Time in atomic units
            
        Examples
        --------
        >>> time_au = PhysicalConstants.ps_to_atomic_units(1.0)
        >>> print(f"1 ps = {time_au:.1f} a.u.")
        """
        return time_ps / cls.TIME_PS_CONVERSION
    
    @classmethod
    def atomic_units_to_ps(cls, time_au: float) -> float:
        """
        Convert time from atomic units to picoseconds.
        
        Parameters
        ----------
        time_au : float
            Time in atomic units
            
        Returns
        -------
        float
            Time in picoseconds
            
        Examples
        --------
        >>> time_ps = PhysicalConstants.atomic_units_to_ps(41341.4)
        >>> print(f"41341.4 a.u. = {time_ps:.1f} ps")
        """
        return time_au * cls.TIME_PS_CONVERSION
    
    @classmethod
    def fs_to_atomic_units(cls, time_fs: float) -> float:
        """
        Convert time from femtoseconds to atomic units.
        
        Parameters
        ----------
        time_fs : float
            Time in femtoseconds
            
        Returns
        -------
        float
            Time in atomic units
            
        Examples
        --------
        >>> time_au = PhysicalConstants.fs_to_atomic_units(1.0)
        >>> print(f"1 fs = {time_au:.1f} a.u.")
        """
        return time_fs / (cls.TIME_PS_CONVERSION * 1000.0)
    
    @classmethod
    def gamma_from_tau_ps(cls, tau_ps: float) -> float:
        r"""
        Calculate damping coefficient from time constant in picoseconds.
        
        For Langevin dynamics, the damping coefficient relates to the time
        constant via:
        
        .. math::
            
            \gamma = \frac{1}{\tau}
        
        where :math:`\tau` is the characteristic relaxation time.
        
        Parameters
        ----------
        tau_ps : float
            Time constant in picoseconds (must be positive)
            
        Returns
        -------
        float
            Damping coefficient :math:`\gamma` in atomic units (inverse time)
            
        Raises
        ------
        ValueError
            If tau_ps is not positive
            
        Examples
        --------
        >>> # 5 ps relaxation time
        >>> gamma = PhysicalConstants.gamma_from_tau_ps(5.0)
        >>> print(f"γ = {gamma:.6f} a.u.⁻¹")
        
        Notes
        -----
        In the overdamped limit (τ → 0), γ → ∞, which is unphysical for
        Langevin dynamics. Use Brownian dynamics instead for very fast
        relaxation.
        """
        if tau_ps <= 0.0:
            raise ValueError(
                f"ERROR: tau_ps must be positive, got {tau_ps} ps.\n"
                f"For Langevin dynamics, gamma = 1/tau, so tau must be > 0.\n"
                f"For overdamped dynamics (tau → 0), use Brownian dynamics instead."
            )
        tau_au = cls.ps_to_atomic_units(tau_ps)
        return 1.0 / tau_au


def unwrap_positions(positions: np.ndarray, 
                    images: np.ndarray, 
                    box_lengths: np.ndarray) -> np.ndarray:
    r"""
    Unwrap particle positions across periodic boundary conditions.
    
    Computes the actual positions of particles that may have been wrapped
    across periodic boundaries during simulation. This is essential for
    calculating properties like dipole moments that depend on absolute
    positions.
    
    **Mathematical Definition:**
    
    The unwrapped position is:
    
    .. math::
        
        \vec{R}_{\text{unwrapped}} = \vec{r}_{\text{wrapped}} + \vec{n} \cdot \vec{L}
    
    where:
    - :math:`\vec{r}_{\text{wrapped}}` is the wrapped position in [0, L)
    - :math:`\vec{n}` is the integer image vector
    - :math:`\vec{L}` is the box length vector
    
    Parameters
    ----------
    positions : np.ndarray
        Wrapped particle positions, shape (N, 3) where N is number of particles
    images : np.ndarray  
        Integer image vectors, shape (N, 3), indicating how many times each
        particle has crossed each boundary
    box_lengths : np.ndarray
        Simulation box lengths in each dimension, shape (3,)
        
    Returns
    -------
    np.ndarray
        Unwrapped positions, shape (N, 3), in the same units as input positions
        
    Examples
    --------
    >>> import numpy as np
    >>> positions = np.array([[0.1, 0.2, 0.3], [9.9, 8.8, 7.7]])
    >>> images = np.array([[0, 0, 0], [1, 1, 1]])  
    >>> box_lengths = np.array([10.0, 10.0, 10.0])
    >>> unwrapped = unwrap_positions(positions, images, box_lengths)
    >>> print(unwrapped)
    [[ 0.1  0.2  0.3]
     [19.9 18.8 17.7]]
    
    Notes
    -----
    This function is critical for computing molecular dipole moments, as
    wrapped positions can lead to incorrect values when molecules are
    split across boundaries.
    
    The unwrapping is performed independently for each spatial dimension.
    """
    return positions + images * box_lengths


# =============================================================================
# JOB MANAGEMENT UTILITIES
# =============================================================================

import os


def get_slurm_info() -> Tuple[Union[int, None], Union[int, None]]:
    """
    Extract SLURM array job information from environment variables.
    
    Returns
    -------
    Tuple[Union[int, None], Union[int, None]]
        (task_id, job_id) where:
        - task_id: SLURM array task ID (SLURM_ARRAY_TASK_ID)
        - job_id: SLURM job ID (SLURM_JOB_ID)
        - Both are None if not running under SLURM
        
    Examples
    --------
    >>> task_id, job_id = get_slurm_info()
    >>> if task_id is not None:
    ...     print(f"Running SLURM array task {task_id} in job {job_id}")
    ... else:
    ...     print("Not running under SLURM")
    """
    import os
    
    task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    job_id = os.environ.get('SLURM_JOB_ID')
    
    if task_id is not None:
        task_id = int(task_id)
    if job_id is not None:
        job_id = int(job_id)
        
    return task_id, job_id


def parse_replicas(replicas_str: Union[str, None]) -> List[int]:
    """
    Parse replica specification string into list of replica numbers.
    
    Supports both comma-separated lists and range specifications.
    
    Parameters
    ----------
    replicas_str : Union[str, None]
        Replica specification string. Formats supported:
        - "1,2,3,4,5": Comma-separated list
        - "1-5": Range specification (inclusive)
        - None: Returns [1] (single replica)
        
    Returns
    -------
    List[int]
        List of replica numbers
        
    Raises
    ------
    ValueError
        If the format is invalid or numbers are out of range
        
    Examples
    --------
    >>> parse_replicas("1,2,3")
    [1, 2, 3]
    >>> parse_replicas("1-5") 
    [1, 2, 3, 4, 5]
    >>> parse_replicas(None)
    [1]
    """
    if replicas_str is None:
        return [1]
    
    replicas_str = replicas_str.strip()
    
    if '-' in replicas_str and ',' not in replicas_str:
        # Range format: "1-5"
        try:
            start, end = map(int, replicas_str.split('-'))
            if start > end or start < 0:
                raise ValueError(f"Invalid range: {replicas_str}")
            return list(range(start, end + 1))
        except ValueError as e:
            raise ValueError(f"Invalid range format '{replicas_str}': {e}")
    else:
        # Comma-separated format: "1,2,3,4,5"
        try:
            replicas = [int(x.strip()) for x in replicas_str.split(',')]
            if any(r < 0 for r in replicas):
                raise ValueError("Replica numbers must be non-negative")
            return sorted(set(replicas))  # Remove duplicates and sort
        except ValueError as e:
            raise ValueError(f"Invalid replica list format '{replicas_str}': {e}") 