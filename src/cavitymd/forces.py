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
    
    Now supports time-varying coupling strength and dissipation using hoomd.variant.
    
    Parameters:
    -----------
    kvector : array_like
        Cavity mode wave vector (currently not used but kept for compatibility)
    couplstr : float or hoomd.variant.Variant
        Coupling strength g in atomic units (can be time-varying)
    omegac : float
        Cavity frequency in atomic units (Hartree)
    phmass : float, optional
        Photon mass, determines K = phmass * omegac² (default: 1.0)
    dissipation : float or hoomd.variant.Variant, optional
        Dissipation rate in atomic units (default: 0.0, can be time-varying)
    force_python : bool, optional
        Force use of Python implementation even if C++ is available (default: False)
    """
    
    def __init__(self, kvector, couplstr, omegac, phmass=1.0, dissipation=0.0, force_python=False):
        # Initialize the base class FIRST - this creates empty _param_dict and _typeparam_dict
        super().__init__()
        
        # Handle variant parameters - convert to variants if needed
        if isinstance(couplstr, (int, float)):
            from hoomd.variant import Constant
            self.couplstr_variant = Constant(float(couplstr))
            self.uses_variant_coupling = False
        elif hasattr(couplstr, '__call__'):  # It's a variant
            from hoomd.variant import Constant
            self.couplstr_variant = couplstr
            self.uses_variant_coupling = not isinstance(couplstr, Constant)
        else:
            raise ValueError("couplstr must be a number or hoomd.variant.Variant")
            
        if isinstance(dissipation, (int, float)):
            from hoomd.variant import Constant
            self.dissipation_variant = Constant(float(dissipation))
            self.uses_variant_dissipation = False
        elif hasattr(dissipation, '__call__'):  # It's a variant
            from hoomd.variant import Constant
            self.dissipation_variant = dissipation
            self.uses_variant_dissipation = not isinstance(dissipation, Constant)
        else:
            # Default to a constant zero variant if dissipation is None or not provided
            from hoomd.variant import Constant
            self.dissipation_variant = Constant(0.0)
            self.uses_variant_dissipation = False

        # Now set up parameter dictionaries using the proper HOOMD methods
        param_dict = hoomd.data.parameterdicts.ParameterDict(
            kvector=hoomd.data.typeconverter.to_type_converter([float, float, float]),
            couplstr=hoomd.variant.Variant,
            omegac=float,
            phmass=float,
            dissipation=hoomd.variant.Variant,
            force_python=bool
        )
        param_dict['kvector'] = list(kvector)
        param_dict['couplstr'] = self.couplstr_variant
        param_dict['omegac'] = omegac
        param_dict['phmass'] = phmass
        param_dict['dissipation'] = self.dissipation_variant
        param_dict['force_python'] = force_python
        
        # Update the existing _param_dict (don't replace it)
        self._param_dict.update(param_dict)
        
        # Store parameters for easy access (no longer need initial values here)
        self.kvector = np.array(kvector)
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
                phmass=phmass,
                dissipation=dissipation
            )
            self._implementation = "python"
            
        else:
            # Use C++ implementation (will be initialized during _attach_hook)
            self._force_impl = None
            self._implementation = "cpp"
        
        print(f"CavityForce initialized using {self._implementation} implementation")
        if self.uses_variant_coupling:
            print(f"  Using time-varying coupling strength")
        else:
            print(f"  Using constant coupling strength: {self.couplstr_variant.value}")
        if self.uses_variant_dissipation:
            print(f"  Using time-varying dissipation")
        else:
            print(f"  Using constant dissipation: {self.dissipation_variant.value}")
    
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
                                self.couplstr_variant,
                                self.phmass,
                                self.dissipation_variant
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
                            self.couplstr_variant,
                            self.phmass,
                            self.dissipation_variant
                        )
                        self._implementation = "cpp"
                        print(f"CPU CavityForceCompute initialized successfully")
                else:
                    # Use CPU implementation for CPU device
                    self._force_impl = _cavitymd.CavityForceCompute(
                        self._simulation.state._cpp_sys_def,
                        self.omegac,
                        self.couplstr_variant,
                        self.phmass,
                        self.dissipation_variant
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
                    phmass=self.phmass,
                    dissipation=self.dissipation
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
    
    def set_forces(self, timestep):
        """
        Compute forces (Python implementation only).
        
        This method is now only used by the pure Python fallback implementation.
        """
        # Call the actual force computation
        if self._implementation in ["python", "python_fallback"]:
            self._force_impl.set_forces(timestep)

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
    def energy(self):
        """Override HOOMD's base Force.energy property to return sum of components instead of per-particle energies"""
        return self.total_cavity_energy
    
    @property
    def forces(self):
        """Get forces array (for HOOMD compatibility)"""
        if hasattr(self._force_impl, 'forces'):
            return self._force_impl.forces
        else:
            # For Python implementation, forces are managed by HOOMD's Custom force
            return None

    def set_params(self, **kwargs):
        """Set force parameters at runtime."""
        for key, value in kwargs.items():
            if key in self._param_dict:
                self._param_dict[key] = value
                setattr(self, key, value)
                # Update implementation if needed
                if self._force_impl and hasattr(self._force_impl, 'setParams'):
                    self._force_impl.setParams(self.omegac, self.couplstr, self.phmass, self.dissipation)
            else:
                raise ValueError(f"Unknown parameter: {key}") 