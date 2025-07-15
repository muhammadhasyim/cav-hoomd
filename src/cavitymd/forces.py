# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Unified force interfaces with C++/Python fallback support."""

from typing import Union, Optional, Dict, Any, List
import hoomd
import warnings
import numpy as np
from hoomd.logging import log
from hoomd.variant import Constant

# Try to import C++ compiled module
try:
    from . import _cavitymd  # C++ compiled module
    _cpp_available = True
except ImportError:
    _cpp_available = False


class CavityForce(hoomd.md.force.Force):
    r"""
    Cavity-molecule interaction force with automatic C++/Python fallback.
    
    Implements the cavity-molecule interaction force from the single-mode cavity Hamiltonian:
    
    .. math::
        
        H = E_{\text{molecular}}^{(0)} + \frac{1}{2} K \tilde{q}_{0,\lambda}^2 + 
            \tilde{\varepsilon}_{0,\lambda} \tilde{q}_{0,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda} + 
            \frac{\tilde{\varepsilon}_{0,\lambda}^2}{2K} \left(\sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}\right)^2
    
    where:
    
    - :math:`\tilde{q}_{0,\lambda}` is the normalized cavity mode coordinate
    - :math:`K = m_{0,\lambda}\omega_{0,\lambda}^2` is the cavity spring constant  
    - :math:`\tilde{\varepsilon}_{0,\lambda}` is the effective coupling strength
    - :math:`d_{ng,\lambda}` is the dipole moment component of molecule n in direction λ
    - :math:`N_{\text{sub}}` is the number of molecules in the simulation cell
    
    **Time-Varying Coupling Support:**
    
    The coupling strength can be time-dependent: :math:`\tilde{\varepsilon}_{0,\lambda}(t) = g(t) \tilde{\varepsilon}_{0,\lambda}^{(0)}`
    
    **Forces Derived from the Hamiltonian:**
    
    Nuclear forces:
    
    .. math::
        
        F_{nj} = -\tilde{\varepsilon}_{0,\lambda}(t) \tilde{q}_{0,\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}} + 
                 \frac{\tilde{\varepsilon}_{0,\lambda}(t)^2}{K} \sum_{l=1}^{N_{\text{sub}}} d_{lg,\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}
    
    Cavity mode force:
    
    .. math::
        
        F_{\text{cavity}} = -K \tilde{q}_{0,\lambda} - \tilde{\varepsilon}_{0,\lambda}(t) \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda} - \gamma(t) \dot{\tilde{q}}_{0,\lambda}
    
    where :math:`\gamma(t)` is the optional time-dependent dissipation coefficient.
    
    **Energy Components:**
    
    - **Harmonic energy**: :math:`E_{\text{harmonic}} = \frac{1}{2} K \tilde{q}_{0,\lambda}^2`
    - **Coupling energy**: :math:`E_{\text{coupling}} = \tilde{\varepsilon}_{0,\lambda}(t) \tilde{q}_{0,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}`
    - **Dipole self-energy**: :math:`E_{\text{dipole}} = \frac{\tilde{\varepsilon}_{0,\lambda}(t)^2}{2K} \left(\sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}\right)^2`
    
    Parameters
    ----------
    kvector : array_like
        Cavity mode wave vector (shape: (3,)). Currently not used but kept for compatibility 
        with future multimode implementations.
    couplstr : float or hoomd.variant.Variant
        Coupling strength :math:`\tilde{\varepsilon}_{0,\lambda}` in atomic units. 
        Can be time-varying using hoomd.variant objects for dynamic coupling experiments.
    omegac : float
        Cavity frequency :math:`\omega_{0,\lambda}` in atomic units (Hartree).
        Related to wavenumber: :math:`\omega = 2\pi c \tilde{\nu}` where :math:`\tilde{\nu}` is in cm⁻¹.
    phmass : float, optional
        Photon effective mass :math:`m_{0,\lambda}` in atomic units. 
        Determines spring constant: :math:`K = m_{0,\lambda}\omega_{0,\lambda}^2`. 
        Default: 1.0 a.u.
    dissipation : float or hoomd.variant.Variant, optional
        Dissipation coefficient :math:`\gamma` in atomic units (inverse time). 
        Can be time-varying for dynamic damping experiments. Default: 0.0 (no damping).
    
    Attributes
    ----------
    implementation : str
        Current implementation being used ('cpp', 'cuda', or 'python')
    uses_variant_coupling : bool
        Whether time-varying coupling is enabled
    uses_variant_dissipation : bool  
        Whether time-varying dissipation is enabled
    
    Methods
    -------
    harmonic_energy() : float
        Returns harmonic oscillator energy component
    coupling_energy() : float
        Returns coupling interaction energy component  
    dipole_self_energy() : float
        Returns dipole self-energy component
    total_cavity_energy() : float
        Returns sum of all cavity energy components
    set_params(**kwargs) : None
        Update force parameters at runtime
    
    Examples
    --------
    **Constant coupling:**
    
    >>> import numpy as np
    >>> from hoomd.cavitymd.forces import CavityForce
    >>> cavity_force = CavityForce(
    ...     kvector=np.array([0, 0, 1]),
    ...     couplstr=0.001,  # 1e-3 a.u.
    ...     omegac=0.00913,  # 2000 cm⁻¹ in a.u.
    ...     phmass=1.0
    ... )
    
    **Time-varying coupling:**
    
    >>> from hoomd.cavitymd.variants import StepVariant
    >>> coupling_variant = StepVariant(
    ...     target_value=0.001,
    ...     switch_time_ps=1.0,
    ...     time_tracker=time_tracker
    ... )
    >>> cavity_force = CavityForce(
    ...     kvector=np.array([0, 0, 1]),
    ...     couplstr=coupling_variant,
    ...     omegac=0.00913,
    ...     dissipation=0.0001  # Add damping
    ... )
    
    **Energy monitoring:**
    
    >>> # Access energy components during simulation
    >>> harmonic = cavity_force.harmonic_energy
    >>> coupling = cavity_force.coupling_energy  
    >>> dipole = cavity_force.dipole_self_energy
    >>> total = cavity_force.total_cavity_energy
    
    Notes
    -----
    - The cavity particle must have type name 'L' (typically typeid=2)
    - Only x,y components of cavity mode and molecular dipoles contribute to coupling
    - Automatically selects optimal implementation: CUDA > C++ > Python fallback
    - Energy conservation is maintained during time-varying coupling transitions
    - Compatible with both finite-q and q=0 cavity modes
    
    See Also
    --------
    hoomd.cavitymd.variants.StepVariant : For time-varying coupling
    hoomd.cavitymd.analysis.EnergyTracker : For monitoring energy conservation
    hoomd.cavitymd.simulation.CavityMDSimulation : For complete simulation setup
    
    References
    ----------
    For the underlying theory, see the cavity molecular dynamics formulation in the 
    theory documentation, particularly the single-mode approximation and time-varying 
    coupling sections.
    """
    
    def __init__(self, 
                 kvector: Union[List[float], np.ndarray], 
                 couplstr: Union[float, hoomd.variant.Variant], 
                 omegac: float, 
                 phmass: float = 1.0, 
                 dissipation: Union[float, hoomd.variant.Variant] = 0.0) -> None:
        # Initialize the base class FIRST - this creates empty _param_dict and _typeparam_dict
        super().__init__()
        
        # Handle variant parameters - convert to variants if needed
        if isinstance(couplstr, (int, float)):
            self.couplstr_variant = Constant(float(couplstr))
            self.uses_variant_coupling = False
        elif hasattr(couplstr, '__call__'):  # It's a variant
            self.couplstr_variant = couplstr
            self.uses_variant_coupling = not isinstance(couplstr, Constant)
        else:
            raise ValueError("couplstr must be a number or hoomd.variant.Variant")
            
        if isinstance(dissipation, (int, float)):
            self.dissipation_variant = Constant(float(dissipation))
            self.uses_variant_dissipation = False
        elif hasattr(dissipation, '__call__'):  # It's a variant
            self.dissipation_variant = dissipation
            self.uses_variant_dissipation = not isinstance(dissipation, Constant)
        else:
            # Default to a constant zero variant if dissipation is None or not provided
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
    def implementation(self) -> str:
        """Return the current implementation being used ('cpp', 'cuda', or 'python')."""
        return self._implementation
    
    @log(requires_run=True)
    def harmonic_energy(self) -> float:
        r"""
        Harmonic oscillator energy component of the cavity mode.
        
        Returns the harmonic energy contribution:
        
        .. math::
            
            E_{\text{harmonic}} = \frac{1}{2} K \tilde{q}_{0,\lambda}^2
        
        where :math:`K = m_{0,\lambda}\omega_{0,\lambda}^2` is the cavity spring constant
        and :math:`\tilde{q}_{0,\lambda}` is the normalized cavity mode coordinate.
        
        Returns
        -------
        float
            Harmonic energy in atomic units (Hartree)
            
        Notes
        -----
        This energy component represents the "photonic" contribution to the total
        cavity energy and corresponds to the kinetic + potential energy of a
        quantum harmonic oscillator in the classical limit.
        """
        if self._implementation in ["cpp", "cuda"]:
            return self._force_impl.getHarmonicEnergy() if self._force_impl else 0.0
        else:
            return getattr(self._force_impl, 'harmonic_energy', 0.0)
    
    @log(requires_run=True)
    def coupling_energy(self) -> float:
        r"""
        Coupling interaction energy between cavity mode and molecular dipoles.
        
        Returns the coupling energy contribution:
        
        .. math::
            
            E_{\text{coupling}} = \tilde{\varepsilon}_{0,\lambda}(t) \tilde{q}_{0,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}
        
        where:
        - :math:`\tilde{\varepsilon}_{0,\lambda}(t)` is the (possibly time-varying) coupling strength
        - :math:`\tilde{q}_{0,\lambda}` is the cavity mode coordinate
        - :math:`d_{ng,\lambda}` are the molecular dipole moment components
        
        Returns
        -------
        float
            Coupling energy in atomic units (Hartree)
            
        Notes
        -----
        This energy component represents the direct interaction between the cavity
        electromagnetic field and the molecular charge distributions. It can be
        positive or negative depending on the relative orientation of the cavity
        field and molecular dipoles.
        """
        if self._implementation in ["cpp", "cuda"]:
            return self._force_impl.getCouplingEnergy() if self._force_impl else 0.0
        else:
            return getattr(self._force_impl, 'coupling_energy', 0.0)
    
    @log(requires_run=True)
    def dipole_self_energy(self) -> float:
        r"""
        Dipole self-energy correction from cavity-induced interactions.
        
        Returns the dipole self-energy contribution:
        
        .. math::
            
            E_{\text{dipole}} = \frac{\tilde{\varepsilon}_{0,\lambda}(t)^2}{2K} \left(\sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}\right)^2
        
        where the summation represents the total molecular dipole moment squared,
        scaled by the square of the coupling strength and inversely by the cavity
        spring constant.
        
        Returns
        -------
        float
            Dipole self-energy in atomic units (Hartree)
            
        Notes
        -----
        This energy component arises from the modification of molecular interactions
        due to the cavity environment. It represents the energy cost/benefit of
        having aligned dipoles in the presence of the cavity mode and is always
        positive (stabilizing).
        """
        if self._implementation in ["cpp", "cuda"]:
            return self._force_impl.getDipoleSelfEnergy() if self._force_impl else 0.0
        else:
            return getattr(self._force_impl, 'dipole_self_energy', 0.0)
    
    @log(requires_run=True)
    def total_cavity_energy(self) -> float:
        r"""
        Total cavity energy as sum of all components.
        
        Returns the complete cavity contribution to the system energy:
        
        .. math::
            
            E_{\text{total}} = E_{\text{harmonic}} + E_{\text{coupling}} + E_{\text{dipole}}
        
        This represents the total modification to the system energy due to
        cavity-molecule coupling.
        
        Returns
        -------
        float
            Total cavity energy in atomic units (Hartree)
            
        Notes
        -----
        This total energy is conserved during time-varying coupling transitions,
        ensuring proper energy conservation in dynamic experiments.
        """
        return self.harmonic_energy + self.coupling_energy + self.dipole_self_energy
    
    @property
    def energy(self) -> float:
        """Override HOOMD's base Force.energy property to return sum of components instead of per-particle energies"""
        return self.total_cavity_energy
    
    @property
    def forces(self) -> Optional[np.ndarray]:
        """Get forces array (for HOOMD compatibility)"""
        if hasattr(self._force_impl, 'forces'):
            return self._force_impl.forces
        else:
            # For Python implementation, forces are managed by HOOMD's Custom force
            return None

    def set_params(self, **kwargs: Any) -> None:
        """
        Set force parameters at runtime.
        
        Parameters
        ----------
        **kwargs : Any
            Force parameters to update. Valid keys include:
            - 'omegac': Cavity frequency
            - 'couplstr': Coupling strength (constant or variant)
            - 'phmass': Photon mass
            - 'dissipation': Dissipation coefficient (constant or variant)
            
        Raises
        ------
        ValueError
            If unknown parameter key is provided
            
        Examples
        --------
        >>> cavity_force.set_params(omegac=0.01, couplstr=0.002)
        """
        for key, value in kwargs.items():
            if key in self._param_dict:
                self._param_dict[key] = value
                setattr(self, key, value)
                # Update implementation if needed
                if self._force_impl and hasattr(self._force_impl, 'setParams'):
                    self._force_impl.setParams(self.omegac, self.couplstr, self.phmass, self.dissipation)
            else:
                raise ValueError(f"Unknown parameter: {key}") 