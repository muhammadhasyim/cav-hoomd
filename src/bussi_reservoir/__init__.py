# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Bussi thermostat implementation with reservoir energy tracking."""

# Version compatibility check
import hoomd
import warnings
_REQUIRED_HOOMD_VERSION = "5.2.0"
_ALLOWED_HOOMD_VERSIONS = ["5.2.0", "5.4.0"]

if hoomd.version.version not in _ALLOWED_HOOMD_VERSIONS:
    raise ImportError(
        f"HOOMD-blue version mismatch: "
        f"BussiReservoir requires version {_REQUIRED_HOOMD_VERSION}, "
        f"but found version {hoomd.version.version}. "
        f"Please install HOOMD-blue {_REQUIRED_HOOMD_VERSION} from source."
    )
elif hoomd.version.version != _REQUIRED_HOOMD_VERSION:
    warnings.warn(
        f"Using HOOMD {hoomd.version.version} instead of tested version {_REQUIRED_HOOMD_VERSION}."
    )

from .thermostats import BussiReservoir

__all__ = ['BussiReservoir'] 