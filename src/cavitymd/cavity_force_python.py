# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Pure Python implementation of cavity force computation."""

from typing import Union, List, Optional
import hoomd
import numpy as np
from .utils import unwrap_positions


class CavityForcePython(hoomd.md.force.Custom):
    r"""
    Pure Python implementation of the cavity-molecule interaction force.
    
    This class serves as a fallback when the optimized C++/CUDA implementations
    are not available. It implements the same physics as the compiled versions
    but with Python performance characteristics.
    
    **Mathematical Implementation:**
    
    Implements the single-mode cavity Hamiltonian:
    
    .. math::
        
        H = E_{\text{molecular}}^{(0)} + \frac{1}{2} K \tilde{q}_{0,\lambda}^2 + 
            \tilde{\varepsilon}_{0,\lambda}(t) \tilde{q}_{0,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda} + 
            \frac{\tilde{\varepsilon}_{0,\lambda}(t)^2}{2K} \left(\sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}\right)^2
    
    **Forces Computed:**
    
    **Molecular forces:**
    
    .. math::
        
        F_{nj} = -\tilde{\varepsilon}_{0,\lambda}(t) \tilde{q}_{0,\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}} + 
                 \frac{\tilde{\varepsilon}_{0,\lambda}(t)^2}{K} \sum_{l=1}^{N_{\text{sub}}} d_{lg,\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}
    
    **Cavity force:**
    
    .. math::
        
        F_{\text{cavity}} = -K \tilde{q}_{0,\lambda} - \tilde{\varepsilon}_{0,\lambda}(t) \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda} - \gamma(t) \dot{\tilde{q}}_{0,\lambda}
    
    **Time-Varying Parameter Support:**
    
    Supports time-dependent coupling :math:`\tilde{\varepsilon}_{0,\lambda}(t)` and 
    dissipation :math:`\gamma(t)` using hoomd.variant objects for dynamic experiments.
    
    **Implementation Details:**
    
    - Uses unwrapped positions for accurate dipole moment calculation
    - Only considers x,y components of cavity mode and molecular dipoles  
    - Cavity particle identified by typeid=2 (type name 'L')
    - Evaluates variants at each timestep for time-varying parameters
    - Maintains energy conservation during parameter changes
    
    Parameters
    ----------
    kvector : array_like
        Cavity mode wave vector (shape: (3,)). Currently unused but kept for compatibility.
    couplstr : float or hoomd.variant.Variant
        Coupling strength :math:`\tilde{\varepsilon}_{0,\lambda}` in atomic units.
        Can be time-varying using hoomd.variant objects.
    omegac : float
        Cavity frequency :math:`\omega_{0,\lambda}` in atomic units (Hartree).
    phmass : float, optional
        Photon effective mass :math:`m_{0,\lambda}` in atomic units. Default: 1.0
    dissipation : float or hoomd.variant.Variant, optional
        Dissipation coefficient :math:`\gamma` in atomic units. Default: 0.0
        
    Attributes
    ----------
    K : float
        Cavity spring constant :math:`K = m_{0,\lambda}\omega_{0,\lambda}^2`
    harmonic_energy : float
        Current harmonic energy component
    coupling_energy : float
        Current coupling energy component
    dipole_self_energy : float
        Current dipole self-energy component
    total_cavity_energy : float
        Total cavity energy (sum of components)
        
    Methods
    -------
    set_forces(timestep) : None
        Compute and set forces for current timestep
        
    Examples
    --------
    **Basic usage:**
    
    >>> from hoomd.cavitymd.cavity_force_python import CavityForcePython
    >>> 
    >>> # Create Python cavity force
    >>> cavity_force = CavityForcePython(
    ...     kvector=[0, 0, 1],
    ...     couplstr=0.001,
    ...     omegac=0.00913,  # 2000 cm⁻¹ in a.u.
    ...     phmass=1.0
    ... )
    
    **With time-varying parameters:**
    
    >>> from hoomd.cavitymd.variants import StepVariant
    >>> coupling_variant = StepVariant(
    ...     target_value=0.001,
    ...     switch_time_ps=1.0,
    ...     time_tracker=time_tracker
    ... )
    >>> 
    >>> cavity_force = CavityForcePython(
    ...     kvector=[0, 0, 1],
    ...     couplstr=coupling_variant,
    ...     omegac=0.00913,
    ...     dissipation=0.0001
    ... )
    
    Notes
    -----
    - Performance is significantly slower than C++/CUDA implementations
    - Used automatically as fallback when compiled extensions are unavailable
    - Maintains identical physics and energy conservation properties
    - Compatible with all variant types for time-dependent parameters
    - Suitable for testing, debugging, and small system calculations
    
    See Also
    --------
    hoomd.cavitymd.forces.CavityForce : Main cavity force interface with auto-fallback
    hoomd.cavitymd.utils.unwrap_positions : For position unwrapping
    hoomd.cavitymd.variants.StepVariant : For time-varying parameters
    """
    
    def __init__(self, 
                 kvector: Union[List[float], np.ndarray], 
                 couplstr: Union[float, hoomd.variant.Variant], 
                 omegac: float, 
                 phmass: float = 1.0, 
                 dissipation: Union[float, hoomd.variant.Variant] = 0.0) -> None:
        """
        Initialize the Python cavity force implementation.
        
        Parameters
        ----------
        kvector : array_like
            Cavity mode wave vector (currently unused)
        couplstr : float or hoomd.variant.Variant
            Coupling strength in atomic units
        omegac : float
            Cavity frequency in atomic units
        phmass : float, optional
            Photon mass in atomic units. Default: 1.0
        dissipation : float or hoomd.variant.Variant, optional
            Dissipation coefficient in atomic units. Default: 0.0
        """
        super().__init__(aniso=False)
        
        # Store parameters
        self.kvector: np.ndarray = np.array(kvector)
        self.couplstr: Union[float, hoomd.variant.Variant] = couplstr
        self.omegac: float = omegac
        self.phmass: float = phmass
        self.dissipation: Union[float, hoomd.variant.Variant] = dissipation
        self.K: float = phmass * omegac**2  # Cavity spring constant
        
        # Initialize energy components
        self.harmonic_energy: float = 0.0
        self.coupling_energy: float = 0.0
        self.dipole_self_energy: float = 0.0
        
        print(f"CavityForcePython initialized:")
        # Handle both regular floats and variant objects for printing
        if hasattr(self.couplstr, '__call__'):
            # It's a variant object, get its initial value
            try:
                coupling_val = self.couplstr(0)  # Get value at timestep 0
                print(f"  Coupling strength: {coupling_val:.6f} a.u. (variant)")
            except:
                print(f"  Coupling strength: {self.couplstr} a.u. (variant)")
        else:
            print(f"  Coupling strength: {self.couplstr:.6f} a.u.")
        
        # Handle dissipation parameter similarly
        if hasattr(self.dissipation, '__call__'):
            try:
                dissipation_val = self.dissipation(0)  # Get value at timestep 0
                print(f"  Dissipation rate: {dissipation_val:.6f} a.u. (variant)")
            except:
                print(f"  Dissipation rate: {self.dissipation} a.u. (variant)")
        else:
            print(f"  Dissipation rate: {self.dissipation:.6f} a.u.")
            
        print(f"  Cavity frequency: {self.omegac:.6f} a.u.")
        print(f"  Photon mass: {self.phmass:.6f} a.u.")
        print(f"  Spring constant K: {self.K:.6f} a.u.")
    
    @property
    def total_cavity_energy(self) -> float:
        """
        Total cavity energy as sum of all components.
        
        Returns
        -------
        float
            Total cavity energy in atomic units
        """
        return self.harmonic_energy + self.coupling_energy + self.dipole_self_energy
    
    def set_forces(self, timestep: int) -> None:
        r"""
        Compute cavity forces and potential energy for the current timestep.
        
        This method is called by HOOMD at each timestep to compute forces.
        Implements the full cavity-molecule interaction Hamiltonian with
        optional time-varying coupling and dissipation.
        
        **Force Calculation Steps:**
        
        1. **Evaluate time-varying parameters** at current timestep
        2. **Locate cavity particle** (typeid=2, type name 'L')
        3. **Unwrap particle positions** across periodic boundaries
        4. **Calculate molecular dipole moment** (x,y components only)
        5. **Compute energy components** from cavity Hamiltonian
        6. **Calculate and apply forces** to all particles
        
        **Mathematical Implementation:**
        
        Molecular forces include coupling and self-energy terms:
        
        .. math::
            
            F_{i} = -g(t) q_i \left[\vec{q}_{\text{cavity,xy}} + \frac{g(t)}{K} \vec{d}_{\text{xy}}\right]
        
        Cavity force includes harmonic, coupling, and dissipation:
        
        .. math::
            
            \vec{F}_{\text{cavity}} = -K \vec{q}_{\text{cavity}} - g(t) \vec{d}_{\text{xy}} - \gamma(t) \vec{v}_{\text{cavity}}
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep number
            
        Notes
        -----
        - Only x,y components of cavity mode and dipoles contribute
        - Cavity particle must have typeid=2 (type name 'L')
        - Energy components are updated for logging/monitoring
        - Forces are applied through HOOMD's custom force interface
        - Handles missing cavity particle gracefully with warnings
        """
        with self._state.cpu_local_snapshot as snap:
            try:
                # Evaluate variant parameters at current timestep
                if hasattr(self.couplstr, '__call__'):
                    coupling_strength = self.couplstr(timestep)
                else:
                    coupling_strength = self.couplstr
                    
                if hasattr(self.dissipation, '__call__'):
                    dissipation_rate = self.dissipation(timestep)
                else:
                    dissipation_rate = self.dissipation
                
                # Find cavity particle (type 'L', should be typeid=2)
                cavity_mask = snap.particles.typeid == 2
                cavity_indices = np.where(cavity_mask)[0]
                
                if len(cavity_indices) == 0:
                    print("Warning: No cavity particle found (typeid=2, type 'L')")
                    self._zero_forces_and_energy()
                    return
                elif len(cavity_indices) > 1:
                    print(f"Warning: Multiple cavity particles found ({len(cavity_indices)}), using first one")
                
                cavity_idx = cavity_indices[0]
                
                # Get unwrapped positions
                box_lengths = np.array([
                    snap.global_box.L[0],
                    snap.global_box.L[1], 
                    snap.global_box.L[2]
                ])
                
                unwrapped_positions = unwrap_positions(
                    snap.particles.position,
                    snap.particles.image,
                    box_lengths
                )
                
                # Calculate molecular dipole moment (only x,y components)
                dipole_moment = np.dot(snap.particles.charge, unwrapped_positions)
                dipole_xy = dipole_moment[:2]  # Only x,y components
                
                # Get cavity particle properties
                cavity_position = unwrapped_positions[cavity_idx][:2]  # Only x,y
                cavity_velocity = snap.particles.velocity[cavity_idx][:2]  # Only x,y
                
                # Calculate energy components
                self.harmonic_energy = 0.5 * self.K * np.dot(cavity_position, cavity_position)
                self.coupling_energy = coupling_strength * np.dot(cavity_position, dipole_xy)
                self.dipole_self_energy = 0.5 * (coupling_strength**2 / self.K) * np.dot(dipole_xy, dipole_xy)
                
                # Compute forces
                with self.cpu_local_force_arrays as arrays:
                    # Initialize arrays
                    arrays.force[:] = 0.0
                    arrays.potential_energy[:] = 0.0
                    
                    # DO NOT assign energy to particle potential energy - prevents double-counting
                    # Energy is accessed directly through force object attributes
                    arrays.potential_energy[cavity_idx] = 0.0
                    
                    # Force on molecular particles: F_i = -g * q_i * [q_xy + (g/K) * d_xy]
                    # Only x,y components contribute
                    force_factor = cavity_position + (coupling_strength / self.K) * dipole_xy
                    num_particles = len(snap.particles.position)
                    for i in range(num_particles):
                        if i != cavity_idx:  # Skip cavity particle
                            charge = snap.particles.charge[i]
                            force_molecular = -coupling_strength * charge * force_factor
                            arrays.force[i] = force_molecular
                    
                    # Force on cavity particle: F_cavity = -K * q - g * d_xy - γ * v_cavity
                    # Add dissipation term: -γ * v_cavity (friction force)
                    force_cavity = (-self.K * cavity_position - 
                                   coupling_strength * dipole_xy - 
                                   dissipation_rate * cavity_velocity)
                    arrays.force[cavity_idx] = force_cavity
                    
            except Exception as e:
                print(f"Error in cavity force calculation: {e}")
                self._zero_forces_and_energy()
    
    def _zero_forces_and_energy(self) -> None:
        """
        Zero all forces and energy components (error handling).
        
        Called when force calculation fails to ensure simulation stability.
        """
        # Zero energy components
        self.harmonic_energy = 0.0
        self.coupling_energy = 0.0 
        self.dipole_self_energy = 0.0
        
        # Zero forces if arrays are available
        try:
            with self.cpu_local_force_arrays as arrays:
                arrays.force[:] = 0.0
                arrays.potential_energy[:] = 0.0
        except:
            pass  # Arrays not available yet 