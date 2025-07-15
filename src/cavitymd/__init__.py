# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Cavity molecular dynamics simulation components.

This package provides components for cavity molecular dynamics simulations.
To use the components, import from the specific submodules:

- cavitymd.forces: Force implementations
- cavitymd.utils: Utility functions and constants
- cavitymd.analysis: Analysis and tracking tools
- cavitymd.simulation: Simulation framework components
- cavitymd.updaters: Particle updaters
- cavitymd.variants: Parameter variants

Example:
    from cavitymd.forces import CavityForce
    from cavitymd.utils import PhysicalConstants
    from cavitymd.analysis import EnergyTracker
"""

# Version compatibility check
import hoomd
_REQUIRED_HOOMD_VERSION = "5.2.0"

if hoomd.version.version != _REQUIRED_HOOMD_VERSION:
    raise ImportError(
        f"HOOMD-blue version mismatch: "
        f"CavityMD requires exactly version {_REQUIRED_HOOMD_VERSION}, "
        f"but found version {hoomd.version.version}. "
        f"Please install HOOMD-blue {_REQUIRED_HOOMD_VERSION} from source."
    ) 