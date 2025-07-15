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
import os
_REQUIRED_HOOMD_VERSION = "5.2.0"

# Check if we're running on Read the Docs
on_rtd = (os.environ.get('READTHEDOCS') == 'True' or 
          os.environ.get('RTD_ENV_NAME') is not None or
          'readthedocs' in os.environ.get('HOSTNAME', '').lower())

if not on_rtd:
    # Only enforce strict version checking when not building documentation
    if hoomd.version.version != _REQUIRED_HOOMD_VERSION:
        raise ImportError(
            f"HOOMD-blue version mismatch: "
            f"CavityMD requires exactly version {_REQUIRED_HOOMD_VERSION}, "
            f"but found version {hoomd.version.version}. "
            f"Please install HOOMD-blue {_REQUIRED_HOOMD_VERSION} from source."
        )
else:
    # On RTD, allow compatible versions but warn about version differences
    current_version = hoomd.version.version
    if current_version != _REQUIRED_HOOMD_VERSION:
        print(f"WARNING: HOOMD version {current_version} detected on Read the Docs")
        print(f"         CavityMD is tested with version {_REQUIRED_HOOMD_VERSION}")
        print(f"         Documentation build proceeding with available version") 