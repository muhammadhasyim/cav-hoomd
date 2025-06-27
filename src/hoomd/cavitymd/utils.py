# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Utility functions and constants for cavity molecular dynamics simulations."""

import numpy as np


class PhysicalConstants:
    """Physical constants and unit conversions for cavity MD simulations."""
    
    HARTREE_TO_CM_MINUS1 = 219474.63
    KB_HARTREE_PER_K = 3.167e-6  # Boltzmann constant in Hartree/K
    ENERGY_JOULES = 4.35974e-18  # Hartree to Joules
    LENGTH_METERS = 5.29177210544e-11  # Bohr to meters
    MASS_KG = 9.1093837139e-31  # Electron mass in kg
    TIME_SECONDS = 2.418884e-17  # Atomic time unit to seconds
    TIME_PS_CONVERSION = 2.418884e-5  # a.u. to picoseconds (corrected from 0.02418884)
    
    @classmethod
    def ps_to_atomic_units(cls, time_ps):
        """
        Convert time from picoseconds to atomic units.
        
        Args:
            time_ps: Time in picoseconds
            
        Returns:
            Time in atomic units
        """
        return time_ps / cls.TIME_PS_CONVERSION
    
    @classmethod
    def atomic_units_to_ps(cls, time_au):
        """
        Convert time from atomic units to picoseconds.
        
        Args:
            time_au: Time in atomic units
            
        Returns:
            Time in picoseconds
        """
        return time_au * cls.TIME_PS_CONVERSION
    
    @classmethod
    def gamma_from_tau_ps(cls, tau_ps):
        """
        Calculate gamma (damping coefficient) from time constant in picoseconds.
        For Langevin dynamics: gamma = 1/tau
        
        Args:
            tau_ps: Time constant in picoseconds
            
        Returns:
            Gamma in atomic units (inverse time)
        """
        if tau_ps <= 0.0:
            raise ValueError(
                f"ERROR: tau_ps must be positive, got {tau_ps} ps.\n"
                f"For Langevin dynamics, gamma = 1/tau, so tau must be > 0.\n"
                f"For overdamped dynamics (tau → 0), use Brownian dynamics instead."
            )
        tau_au = cls.ps_to_atomic_units(tau_ps)
        return 1.0 / tau_au


def unwrap_positions(positions, images, box_lengths):
    """
    Unwrap particle positions across periodic boundaries.
    
    Args:
        positions: Array of wrapped positions (N x 3)
        images: Array of image flags (N x 3)
        box_lengths: Array of box dimensions (3,)
        
    Returns:
        Array of unwrapped positions (N x 3)
    """
    # Convert inputs to numpy arrays if they aren't already
    pos = np.asarray(positions)
    img = np.asarray(images)
    box = np.asarray(box_lengths)
    
    # Unwrap by adding box lengths multiplied by image flags
    return pos + img * box[None, :] 