# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Bussi thermostat implementation with reservoir energy tracking."""

# Version compatibility check
import hoomd
_REQUIRED_HOOMD_VERSION = "5.2.0"

if hoomd.version.version != _REQUIRED_HOOMD_VERSION:
    raise ImportError(
        f"HOOMD-blue version mismatch: "
        f"BussiReservoir requires exactly version {_REQUIRED_HOOMD_VERSION}, "
        f"but found version {hoomd.version.version}. "
        f"Please install HOOMD-blue {_REQUIRED_HOOMD_VERSION} from source."
    )

from .thermostats import BussiReservoir

__all__ = ['BussiReservoir'] 