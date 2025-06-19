# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Extended Bussi thermostat with reservoir energy tracking.

This HOOMD-blue plugin provides an extended version of the Bussi stochastic
velocity rescaling thermostat that tracks the energy dumped into the thermal
reservoir over time.
"""

from . import thermostats
from .version import version as __version__

__all__ = ['thermostats', '__version__']
