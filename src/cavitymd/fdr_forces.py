# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""External perturbation forces for FDR measurements."""

from typing import Union, List, Optional
import hoomd
import numpy as np
from hoomd.logging import log

# Try to import C++ compiled module
try:
    from hoomd.cavitymd import _cavitymd  # C++ compiled module
    _cpp_available = True
except ImportError:
    try:
        from .. import _cavitymd  # Alternative relative import
        _cpp_available = True
    except ImportError:
        _cpp_available = False


class PerturbationForce(hoomd.md.force.Force):
    r"""
    External perturbation force for FDR measurements.
    
    Implements the external perturbation force from the potential:
    
    .. math::
        
        U_{\text{ext}}(\mathbf{r}) = -\text{sign} \cdot h_0 \cos(\mathbf{k} \cdot \mathbf{r})
    
    The resulting force is:
    
    .. math::
        
        \mathbf{F}_{\text{pert}} = \text{sign} \cdot h_0 \mathbf{k} \sin(\mathbf{k} \cdot \mathbf{r})
    
    where:
    
    - :math:`h_0` is the perturbation amplitude (small for linear response)
    - :math:`\mathbf{k}` is the user-specified wave vector  
    - :math:`\text{sign}` is +1 for (+) clone or -1 for (-) clone
    
    **For FDR measurements:**
    
    Two clones are created from the same initial state at waiting time :math:`t_w`:
    
    - **Clone (+)**: Uses sign = +1, applies +h₀ perturbation
    - **Clone (−)**: Uses sign = -1, applies -h₀ perturbation
    
    The integrated response is computed as:
    
    .. math::
        
        \chi_k(t, t_w) = \frac{\rho_k^{(+)}(t) - \rho_k^{(-)}(t)}{2 N h_0}
    
    **Important notes:**
    
    - Only applied to molecular particles (cavity particle excluded)
    - Requires small h₀ for linear response regime
    - Compatible with both CPU and GPU execution
    
    Parameters
    ----------
    kvector : array_like
        Wave vector for perturbation (shape: (3,)). User specifies based on 
        system's structure factor S(k) or desired length scale.
    amplitude : float
        Perturbation amplitude h₀. Should be small (typically 1e-6 to 1e-4)
        to ensure linear response regime.
    sign : float, optional
        Sign of perturbation: +1 for (+) clone, -1 for (-) clone. Default: +1.
    start_time_ps : float, optional
        Time in picoseconds to start applying the perturbation. Default: 0.0
    stop_time_ps : float, optional  
        Time in picoseconds to stop applying the perturbation. If None, runs indefinitely. Default: None
    time_tracker : ElapsedTimeTracker, optional
        Time tracker for accurate timing control. If None, timing control disabled. Default: None
    
    Attributes
    ----------
    k_magnitude : float
        Magnitude of the k-vector |k|
    perturbation_energy : float
        Total perturbation energy (logged automatically)
    
    Methods
    -------
    set_params(**kwargs) : None
        Update force parameters at runtime
    enable() : None
        Enable the perturbation force
    disable() : None
        Disable the perturbation force
        
    Examples
    --------
    **Basic perturbation force:**
    
    >>> import numpy as np
    >>> from hoomd.cavitymd.fdr_forces import PerturbationForce
    >>> 
    >>> # Create perturbation for (+) clone
    >>> pert_plus = PerturbationForce(
    ...     kvector=np.array([1.0, 0.0, 0.0]),  # k in x-direction
    ...     amplitude=1e-6,                     # Small for linearity
    ...     sign=+1.0                           # (+) clone
    ... )
    
    **For FDR clone system:**
    
    >>> # Create both perturbation forces
    >>> k_vector = np.array([2.5, 0.0, 0.0])  # User chooses based on S(k)
    >>> h0 = 1e-6                             # Small amplitude
    >>> 
    >>> pert_plus = PerturbationForce(kvector=k_vector, amplitude=h0, sign=+1.0)
    >>> pert_minus = PerturbationForce(kvector=k_vector, amplitude=h0, sign=-1.0)
    
    **Runtime parameter updates:**
    
    >>> # Change amplitude for linearity testing
    >>> pert_plus.amplitude = 1e-5
    >>> 
    >>> # Disable temporarily
    >>> pert_plus.disable()
    >>> # ... run some steps ...
    >>> pert_plus.enable()
    
    Notes
    -----
    - The perturbation is only applied to molecular particles, not cavity particles
    - For FDR measurements, both (+) and (-) clones should use identical parameters
      except for the sign
    - Monitor perturbation_energy to verify energy conservation when needed
    - Automatically selects CPU or GPU implementation based on device availability
    
    See Also
    --------
    hoomd.cavitymd.fdr_manager : For complete FDR measurement workflow
    hoomd.cavitymd.analysis.IntegratedResponseTracker : For response calculation
    
    References
    ----------
    For the underlying FDR theory and implementation details, see:
    - Gnan et al., "Predicting the Effective Temperature of a Glass" (2009)
    - Literature on fluctuation-dissipation relation violations in glassy systems
    """
    
    def __init__(self, 
                 kvector: Union[List[float], np.ndarray], 
                 amplitude: float, 
                 sign: float = 1.0,
                 start_time_ps: float = 0.0,
                 stop_time_ps: Union[float, None] = None,
                 time_tracker: Union[object, None] = None) -> None:
        
        # Initialize the base class
        super().__init__()
        
        # Process k-vector input
        self._kvector = np.array(kvector, dtype=float)
        if self._kvector.shape != (3,):
            raise ValueError("kvector must be a 3-element array")
            
        self._amplitude = float(amplitude)
        self._sign = float(sign)
        
        # Compute k-vector properties
        self.k_magnitude = np.linalg.norm(self._kvector)
        
        # Convert to HOOMD format
        kvector_hoomd = hoomd.data.typeconverter.to_type_converter([float, float, float])(self._kvector)
        
        # Parameter validation
        if self.k_magnitude == 0:
            raise ValueError("k-vector cannot be zero")
        if abs(self._sign) != 1.0:
            self._sign = 1.0 if self._sign > 0 else -1.0
            
        # Create parameter dict
        param_dict = hoomd.data.parameterdicts.ParameterDict(
            kvector=hoomd.data.typeconverter.to_type_converter([float, float, float]),
            amplitude=float,
            sign=float
        )
        param_dict['kvector'] = kvector_hoomd
        param_dict['amplitude'] = self._amplitude
        param_dict['sign'] = self._sign
        
        # Store parameters
        self._param_dict = param_dict
        
        # Timing control parameters
        self.start_time_ps = float(start_time_ps)
        self.stop_time_ps = float(stop_time_ps) if stop_time_ps is not None else None
        self.time_tracker = time_tracker
        self._timing_enabled = time_tracker is not None
        self._is_active = False  # Track if currently applying force
        self._has_started = False
        self._has_stopped = False
        
        # Validate timing
        if self.stop_time_ps is not None and self.stop_time_ps <= self.start_time_ps:
            raise ValueError("stop_time_ps must be greater than start_time_ps")
        
        # Log creation
        print(f"PerturbationForce created:")
        print(f"  k-vector: [{self._kvector[0]:.3f}, {self._kvector[1]:.3f}, {self._kvector[2]:.3f}]")
        print(f"  |k| = {self.k_magnitude:.3f}")
        print(f"  amplitude = {self._amplitude:.2e}")
        print(f"  sign = {self._sign:+.0f}")
        if self._timing_enabled:
            print(f"  start_time_ps: {self.start_time_ps}")
            if self.stop_time_ps is not None:
                print(f"  stop_time_ps: {self.stop_time_ps}")
            else:
                print(f"  stop_time_ps: None (runs indefinitely)")
        else:
            print(f"  timing control: disabled (no time_tracker)")
        
    def _attach_hook(self):
        """Initialize the C++ compute object when attached to simulation."""
        # Create appropriate C++ object based on device
        if isinstance(self._simulation.device, hoomd.device.CPU):
            if not _cpp_available:
                raise RuntimeError("C++ module not available for PerturbationForce")
            cpp_class = _cavitymd.PerturbationForceCompute
        else:
            if not _cpp_available:
                raise RuntimeError("C++ module not available for PerturbationForce")
            cpp_class = _cavitymd.PerturbationForceComputeGPU
            
        # Convert k-vector to required format
        kvector_scalar3 = hoomd._hoomd.make_scalar3(float(self._kvector[0]), 
                                                     float(self._kvector[1]), 
                                                     float(self._kvector[2]))
        
        # Create C++ object
        self._cpp_obj = cpp_class(self._simulation.state._cpp_sys_def,
                                  kvector_scalar3,
                                  self._amplitude,
                                  self._sign)
        
        # If timing is enabled, start disabled and let timing control handle enabling
        if self._timing_enabled:
            self._cpp_obj.setEnabled(False)
        
        super()._attach_hook()
        
    def _update_timing_control(self):
        """Update force enabling based on timing control."""
        if not self._timing_enabled or not hasattr(self, '_cpp_obj'):
            return
            
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should start
        if current_time_ps >= self.start_time_ps and not self._has_started:
            if not self._is_active:
                print(f"🟢 PerturbationForce turned ON at t = {current_time_ps:.3f} ps")
                self._cpp_obj.setEnabled(True)
                self._is_active = True
                self._has_started = True
        
        # Check if we should stop  
        if (self.stop_time_ps is not None and 
            current_time_ps >= self.stop_time_ps and 
            not self._has_stopped):
            if self._is_active:
                print(f"🔴 PerturbationForce turned OFF at t = {current_time_ps:.3f} ps")
                self._cpp_obj.setEnabled(False)
                self._is_active = False
                self._has_stopped = True
                
    @property
    def amplitude(self):
        """float: Perturbation amplitude h₀."""
        return self._amplitude
        
    @amplitude.setter
    def amplitude(self, value):
        self._amplitude = float(value)
        self._param_dict['amplitude'] = self._amplitude
        if hasattr(self, '_cpp_obj'):
            kvector_scalar3 = (float(self._kvector[0]), float(self._kvector[1]), float(self._kvector[2]))
            self._cpp_obj.setParams(kvector_scalar3, self._amplitude, self._sign)
            
    @property
    def sign(self):
        """float: Perturbation sign (+1 or -1)."""
        return self._sign
        
    @sign.setter  
    def sign(self, value):
        self._sign = 1.0 if value > 0 else -1.0
        self._param_dict['sign'] = self._sign
        if hasattr(self, '_cpp_obj'):
            kvector_scalar3 = (float(self._kvector[0]), float(self._kvector[1]), float(self._kvector[2]))
            self._cpp_obj.setParams(kvector_scalar3, self._amplitude, self._sign)
            
    @property
    def kvector(self):
        """numpy.ndarray: Wave vector for perturbation."""
        return self._kvector.copy()
        
    @kvector.setter
    def kvector(self, value):
        new_kvector = np.array(value, dtype=float)
        if new_kvector.shape != (3,):
            raise ValueError("kvector must be a 3-element array")
        if np.linalg.norm(new_kvector) == 0:
            raise ValueError("k-vector cannot be zero")
            
        self._kvector = new_kvector
        self.k_magnitude = np.linalg.norm(self._kvector)
        kvector_hoomd = hoomd.data.typeconverter.to_type_converter([float, float, float])(self._kvector)
        self._param_dict['kvector'] = kvector_hoomd
        
        if hasattr(self, '_cpp_obj'):
            kvector_scalar3 = (float(self._kvector[0]), float(self._kvector[1]), float(self._kvector[2]))
            self._cpp_obj.setParams(kvector_scalar3, self._amplitude, self._sign)
            
    def enable(self):
        """Enable the perturbation force."""
        if hasattr(self, '_cpp_obj'):
            self._cpp_obj.setEnabled(True)
            
    def disable(self):
        """Disable the perturbation force."""
        if hasattr(self, '_cpp_obj'):
            self._cpp_obj.setEnabled(False)
            
    @log(category="scalar")
    def perturbation_energy(self):
        """float: Total perturbation energy."""
        if hasattr(self, '_cpp_obj'):
            return self._cpp_obj.getPerturbationEnergy()
        return 0.0
    
    # Timing control properties
    @property 
    def has_started(self) -> bool:
        """Whether the perturbation has started."""
        return self._has_started
        
    @property
    def has_stopped(self) -> bool:
        """Whether the perturbation has been stopped."""
        return self._has_stopped
        
    @property 
    def is_active(self) -> bool:
        """Whether the perturbation is currently active."""
        return self._is_active


class PerturbationTimingUpdater(hoomd.custom.Action):
    """Timing control updater for PerturbationForce objects."""
    
    def __init__(self, perturbation_forces: List[PerturbationForce]):
        """
        Initialize timing updater.
        
        Parameters
        ---------- 
        perturbation_forces : List[PerturbationForce]
            List of PerturbationForce objects to control
        """
        self.perturbation_forces = perturbation_forces
        
    def act(self, timestep):
        """Update timing control for all registered forces."""
        for force in self.perturbation_forces:
            force._update_timing_control()


class DipoleResponseForce(hoomd.md.force.Force):
    r"""
    Electric field force for dipole moment FDR measurements.
    
    Implements a static electric field that couples to the total dipole moment:
    
    .. math::
        
        U_{\text{dipole}}(\{\mathbf{r}_i\}) = -\text{sign} \cdot E_0 \sum_i q_i (\mathbf{E} \cdot \mathbf{r}_i)
    
    The resulting force on each particle is:
    
    .. math::
        
        \mathbf{F}_i = \text{sign} \cdot E_0 q_i \mathbf{E}
    
    where:
    
    - :math:`E_0` is the electric field strength (small for linear response)
    - :math:`\mathbf{E}` is the unit electric field direction vector
    - :math:`q_i` is the charge on particle i
    - :math:`\text{sign}` is +1 for (+) clone or -1 for (-) clone
    
    **For Dipole Moment FDR:**
    
    This force enables measurement of the dipole moment susceptibility by applying
    a conjugate field to the total dipole moment. Two clones are created:
    
    - **Clone (+)**: Uses sign = +1, applies +E₀ field
    - **Clone (−)**: Uses sign = -1, applies -E₀ field
    
    The dipole susceptibility is computed as:
    
    .. math::
        
        \chi_{\mu}(t, t_w) = \frac{\langle\mu^{(+)}(t)\rangle - \langle\mu^{(-)}(t)\rangle}{2 E_0}
    
    **Physics Notes:**
    
    - Applied to all charged particles (molecular particles)
    - Cavity particle typically excluded (set to zero charge or neutral)
    - Field strength E₀ should be small for linear response
    - Related to dielectric response and polarizability
    
    **Comparison to Cavity Force:**
    
    This is like a static version of the cavity force but:
    - Acts on charges uniformly (not position-dependent)
    - Couples to total dipole moment (not local density)
    - Used for equilibrium fluctuation-response relations
    
    Parameters
    ----------
    field_vector : array_like
        Electric field direction vector (shape: (3,)). Will be normalized internally.
    field_strength : float
        Electric field strength E₀. Should be small (typically 1e-6 to 1e-4)
        to ensure linear response regime.
    sign : float, optional
        Sign of field: +1 for (+) clone, -1 for (-) clone. Default: +1.
    exclude_cavity : bool, optional
        Whether to exclude cavity particles from the field. Default: True
    
    Examples
    --------
    **Basic dipole response force:**
    
    >>> from hoomd.cavitymd.fdr_forces import DipoleResponseForce
    >>> 
    >>> # Apply electric field in z-direction
    >>> field = DipoleResponseForce(
    ...     field_vector=[0, 0, 1],
    ...     field_strength=1e-5,
    ...     sign=+1  # (+) clone
    ... )
    >>> 
    >>> # Add to simulation
    >>> sim.operations.integrator.forces.append(field)
    
    **Fork-and-clone FDR measurement:**
    
    >>> # Create two clones with opposite signs
    >>> field_plus = DipoleResponseForce(
    ...     field_vector=[0, 0, 1],
    ...     field_strength=1e-5,
    ...     sign=+1
    ... )
    >>> field_minus = DipoleResponseForce(
    ...     field_vector=[0, 0, 1],
    ...     field_strength=1e-5,
    ...     sign=-1
    ... )
    
    Notes
    -----
    - GPU-accelerated C++/CUDA backend for high performance
    - Automatic GPU/CPU selection based on execution context
    - Energy logging and force magnitude tracking included
    - Compatible with existing FDR workflow infrastructure
    - Used in conjunction with DipoleMomentFDRTracker for complete FDR analysis
    """
    
    def __init__(self,
                 field_vector: Union[List[float], np.ndarray],
                 field_strength: float,
                 sign: float = 1.0,
                 exclude_cavity: bool = True):
        
        super().__init__()
        
        # Validate and store field parameters
        field_vector = np.array(field_vector, dtype=np.float64)
        if len(field_vector) != 3:
            raise ValueError("field_vector must have exactly 3 components")
        
        # Normalize field vector
        field_magnitude = np.linalg.norm(field_vector)
        if field_magnitude < 1e-12:
            raise ValueError("field_vector cannot be zero vector")
        field_vector = field_vector / field_magnitude
        
        # Store parameters for later C++ object creation
        self.field_vector = field_vector
        self.field_strength = float(field_strength)
        self.sign = float(sign)
        self.exclude_cavity = bool(exclude_cavity)
        
        # C++ object will be created in _attach_hook when simulation is available
        self._cpp_obj = None
    
    def _attach_hook(self):
        """Attach the C++ compute to the simulation."""
        if isinstance(self._simulation.device, hoomd.device.CPU):
            # Use CPU implementation
            if hasattr(_cavitymd, 'DipoleResponseForceCompute'):
                self._cpp_obj = _cavitymd.DipoleResponseForceCompute(
                    self._simulation.state._cpp_sys_def,
                    (float(self.field_vector[0]), float(self.field_vector[1]), float(self.field_vector[2])),
                    self.field_strength,
                    self.sign,
                    self.exclude_cavity
                )
        else:
            # Use GPU implementation
            if hasattr(_cavitymd, 'DipoleResponseForceComputeGPU'):
                self._cpp_obj = _cavitymd.DipoleResponseForceComputeGPU(
                    self._simulation.state._cpp_sys_def,
                    (float(self.field_vector[0]), float(self.field_vector[1]), float(self.field_vector[2])),
                    self.field_strength,
                    self.sign,
                    self.exclude_cavity
                )
        
        # Enable the force by default
        if self._cpp_obj is not None:
            self._cpp_obj.setEnabled(True)
            
        super()._attach_hook()
    
    @log(category="scalar", default=True)
    def energy(self):
        """Get total electric field energy."""
        if self._cpp_obj is not None:
            return self._cpp_obj.getElectricFieldEnergy()
        return 0.0
    
    def setParams(self, field_vector=None, field_strength=None, sign=None, exclude_cavity=None):
        """Update force parameters.
        
        Parameters
        ----------
        field_vector : array_like, optional
            New electric field direction
        field_strength : float, optional  
            New electric field strength
        sign : float, optional
            New field sign
        exclude_cavity : bool, optional
            New exclude cavity setting
        """
        # Update stored parameters
        if field_vector is not None:
            field_vector = np.array(field_vector, dtype=np.float64)
            field_magnitude = np.linalg.norm(field_vector)
            if field_magnitude > 1e-12:
                self.field_vector = field_vector / field_magnitude
            else:
                raise ValueError("field_vector cannot be zero vector")
        
        if field_strength is not None:
            self.field_strength = float(field_strength)
        
        if sign is not None:
            self.sign = float(sign)
            
        if exclude_cavity is not None:
            self.exclude_cavity = bool(exclude_cavity)
        
        # Update C++ object
        if self._cpp_obj is not None:
            field_vector_scalar3 = (float(self.field_vector[0]), float(self.field_vector[1]), float(self.field_vector[2]))
            self._cpp_obj.setParams(field_vector_scalar3, self.field_strength, self.sign, self.exclude_cavity)
    
    def setEnabled(self, enabled: bool):
        """Enable or disable the electric field.
        
        Parameters
        ----------
        enabled : bool
            Whether to enable the electric field
        """
        if self._cpp_obj is not None:
            self._cpp_obj.setEnabled(enabled)
    
    def getEnabled(self) -> bool:
        """Check if the electric field is enabled.
        
        Returns
        -------
        bool
            Whether the electric field is enabled
        """
        if self._cpp_obj is not None:
            return self._cpp_obj.getEnabled()
        return False
