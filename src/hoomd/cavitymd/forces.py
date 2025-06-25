# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Unified force interfaces with C++/Python fallback support."""

import hoomd
import warnings
import numpy as np
from hoomd.logging import log

# Try to import C++ compiled module
try:
    from . import _cavitymd  # C++ compiled module
    _cpp_available = True
except ImportError:
    _cpp_available = False

from .cavity_force_python import CavityForcePython


class CavityForce(hoomd.md.force.Force):
    """
    Cavity force with automatic C++/Python fallback.
    
    Implements the cavity-molecule interaction force from the Hamiltonian:
    H = (1/2) * K * q² + g * q · d + (g²/2K) * d²
    
    where q is the cavity mode position, d is the molecular dipole moment,
    g is the coupling strength, and K = phmass * omegac².
    
    Parameters:
    -----------
    kvector : array_like
        Cavity mode wave vector (currently not used but kept for compatibility)
    couplstr : float
        Coupling strength g in atomic units
    omegac : float
        Cavity frequency in atomic units (Hartree)
    phmass : float, optional
        Photon mass, determines K = phmass * omegac² (default: 1.0)
    force_python : bool, optional
        Force use of Python implementation even if C++ is available (default: False)
    """
    
    def __init__(self, kvector, couplstr, omegac, phmass=1.0, force_python=False):
        # Initialize the base class FIRST - this creates empty _param_dict and _typeparam_dict
        super().__init__()
        
        # Now set up parameter dictionaries using the proper HOOMD methods
        param_dict = hoomd.data.parameterdicts.ParameterDict(
            kvector=hoomd.data.typeconverter.to_type_converter([float, float, float]),
            couplstr=float,
            omegac=float,
            phmass=float,
            force_python=bool
        )
        param_dict['kvector'] = list(kvector)
        param_dict['couplstr'] = couplstr  
        param_dict['omegac'] = omegac
        param_dict['phmass'] = phmass
        param_dict['force_python'] = force_python
        
        # Update the existing _param_dict (don't replace it)
        self._param_dict.update(param_dict)
        
        # Store parameters for easy access
        self.kvector = np.array(kvector)
        self.couplstr = couplstr
        self.omegac = omegac
        self.phmass = phmass
        
        # Determine which implementation to use
        if force_python or not _cpp_available:
            if not force_python and not _cpp_available:
                warnings.warn(
                    "C++ cavity force implementation not available, falling back to Python. "
                    "For better performance, compile the C++ module.",
                    UserWarning
                )
            
            # Use Python implementation
            self._force_impl = CavityForcePython(
                kvector=kvector,
                couplstr=couplstr,
                omegac=omegac,
                phmass=phmass
            )
            self._implementation = "python"
            
        else:
            # Use C++ implementation (will be initialized during _attach_hook)
            self._force_impl = None
            self._implementation = "cpp"
        
        print(f"CavityForce initialized using {self._implementation} implementation")
    
    def _attach_hook(self):
        """Called when force is attached to simulation"""
        if self._implementation == "cpp" and self._force_impl is None:
            # Initialize C++ force implementation now that we have system definition
            try:
                # Check if we're running on GPU and GPU implementation is available
                device = self._simulation.device
                if hasattr(device, 'gpu_ids') or 'GPU' in str(type(device)):
                    # Try GPU implementation first
                    try:
                        if hasattr(_cavitymd, 'CavityForceComputeGPU'):
                            self._force_impl = _cavitymd.CavityForceComputeGPU(
                                self._simulation.state._cpp_sys_def,
                                self.omegac,
                                self.couplstr,
                                self.phmass
                            )
                            self._implementation = "cuda"
                            print(f"CUDA CavityForceComputeGPU initialized successfully")
                        else:
                            raise AttributeError("GPU implementation not available")
                    except Exception as gpu_error:
                        print(f"GPU implementation failed ({gpu_error}), falling back to CPU")
                        # Fall back to CPU implementation
                        self._force_impl = _cavitymd.CavityForceCompute(
                            self._simulation.state._cpp_sys_def,
                            self.omegac,
                            self.couplstr,
                            self.phmass
                        )
                        self._implementation = "cpp"
                        print(f"CPU CavityForceCompute initialized successfully")
                else:
                    # Use CPU implementation for CPU device
                    self._force_impl = _cavitymd.CavityForceCompute(
                        self._simulation.state._cpp_sys_def,
                        self.omegac,
                        self.couplstr,
                        self.phmass
                    )
                    print(f"CPU CavityForceCompute initialized successfully")
                
                # Set the C++ object for HOOMD's Force interface
                self._cpp_obj = self._force_impl
                
            except Exception as e:
                warnings.warn(
                    f"Failed to initialize C++ cavity force ({e}), falling back to Python",
                    UserWarning
                )
                # Fallback to Python implementation
                self._force_impl = CavityForcePython(
                    kvector=self.kvector,
                    couplstr=self.couplstr,
                    omegac=self.omegac,
                    phmass=self.phmass
                )
                self._implementation = "python_fallback"
        
        # For Python implementation, set up as a Custom force
        if self._implementation in ["python", "python_fallback"]:
            # Initialize the Python implementation with simulation state
            if hasattr(self._force_impl, '_simulation') and self._force_impl._simulation is None:
                self._force_impl._simulation = self._simulation
            # Create a custom force compute object
            self._cpp_obj = hoomd.md._md.CustomForceCompute(
                self._simulation.state._cpp_sys_def, 
                self.set_forces, 
                False  # aniso=False
            )
        
        # Call parent attach hook LAST
        super()._attach_hook()
        
        # Call implementation attach hook if needed  
        if hasattr(self._force_impl, '_attach_hook'):
            self._force_impl._attach_hook()
    
    @property
    def implementation(self):
        """Return the current implementation being used ('cpp', 'cuda', or 'python')."""
        return self._implementation
    
    @log(requires_run=True)
    def harmonic_energy(self):
        """Harmonic oscillator energy component: (1/2) * K * q²"""
        if self._implementation in ["cpp", "cuda"]:
            return self._force_impl.getHarmonicEnergy() if self._force_impl else 0.0
        else:
            return getattr(self._force_impl, 'harmonic_energy', 0.0)
    
    @log(requires_run=True)
    def coupling_energy(self):
        """Coupling interaction energy component: g * (q · d)"""
        if self._implementation in ["cpp", "cuda"]:
            return self._force_impl.getCouplingEnergy() if self._force_impl else 0.0
        else:
            return getattr(self._force_impl, 'coupling_energy', 0.0)
    
    @log(requires_run=True)
    def dipole_self_energy(self):
        """Dipole self-energy component: (g²/2K) * d²"""
        if self._implementation in ["cpp", "cuda"]:
            return self._force_impl.getDipoleSelfEnergy() if self._force_impl else 0.0
        else:
            return getattr(self._force_impl, 'dipole_self_energy', 0.0)
    
    @log(requires_run=True)
    def total_cavity_energy(self):
        """Total cavity energy: sum of all energy components"""
        return self.harmonic_energy + self.coupling_energy + self.dipole_self_energy
    
    @property
    def forces(self):
        """Get forces array (for HOOMD compatibility)"""
        if hasattr(self._force_impl, 'forces'):
            return self._force_impl.forces
        else:
            # For Python implementation, forces are managed by HOOMD's Custom force
            return None
    
    def _detach_hook(self):
        """Called when force is detached from simulation"""
        if hasattr(self._force_impl, '_detach_hook'):
            self._force_impl._detach_hook()
        super()._detach_hook()
    
    def set_forces(self, timestep):
        """For Python implementation compatibility"""
        if hasattr(self._force_impl, 'set_forces'):
            return self._force_impl.set_forces(timestep)
    
    # Removed __getattr__ method temporarily to debug _typeparam_dict issue 