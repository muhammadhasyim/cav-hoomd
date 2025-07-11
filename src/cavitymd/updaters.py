# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Cavity particle displacer."""

import hoomd
from hoomd.operation import Updater
from hoomd.variant import Variant
from hoomd.trigger import Trigger, On

# Try to import C++ compiled module
try:
    from . import _cavitymd
    _cpp_available = True
except ImportError:
    _cpp_available = False

class CavityParticleDisplacer(Updater):
    """Displaces the cavity particle to its new equilibrium position.

    Args:
        trigger (hoomd.trigger.Trigger): The trigger to activate this updater.
        couplstr (hoomd.variant.Variant): The coupling strength variant.
        omegac (float): The cavity frequency.
        phmass (float): The photon mass.
    """

    def __init__(self, trigger, couplstr, omegac, phmass=1.0):
        super().__init__(trigger)
        
        # Store parameters for later use when attached to simulation
        self._couplstr = couplstr
        self._omegac = omegac
        self._phmass = phmass
        self._cpp_obj = None

        if not _cpp_available:
            raise RuntimeError("C++ implementation of CavityParticleDisplacer not available.")

    def _attach_hook(self):
        """Called when the updater is attached to a simulation."""
        # Now we can create the C++ object since we have access to the simulation
        if _cpp_available and self._cpp_obj is None:
            self._cpp_obj = _cavitymd.CavityParticleDisplacer(
                self._simulation.state._cpp_sys_def,
                self.trigger,
                self._couplstr,
                self._omegac,
                self._phmass
            ) 