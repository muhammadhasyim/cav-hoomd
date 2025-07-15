# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Custom variant classes for cavity molecular dynamics."""

from typing import Union, Optional
import hoomd.variant
import numpy as np
from .utils import PhysicalConstants
from .analysis import ElapsedTimeTracker


class StepVariant(hoomd.variant.Variant):
    r"""
    Step function variant for instantaneous coupling switching.
    
    Implements a discontinuous step function that switches from 0.0 to a target
    value at a specified switch time, enabling study of time-varying cavity coupling:
    
    .. math::
        
        g(t) = \begin{cases}
        0 & \text{if } t < t_{\text{switch}} \\
        g_{\text{target}} & \text{if } t \geq t_{\text{switch}}
        \end{cases}
    
    This creates a time-dependent coupling strength:
    
    .. math::
        
        \tilde{\varepsilon}_{0,\lambda}(t) = g(t) \tilde{\varepsilon}_{0,\lambda}^{(0)}
    
    where :math:`\tilde{\varepsilon}_{0,\lambda}^{(0)}` is the base coupling strength.
    
    **Physical Interpretation:**
    
    The step variant simulates experimental scenarios where cavity coupling is
    suddenly activated, such as:
    
    - Pump-probe experiments with laser activation
    - Sudden cavity tuning into molecular resonance  
    - Dynamic switching protocols
    - Non-equilibrium cavity-coupling studies
    
    **Energy Conservation:**
    
    During switching, total energy is conserved. However, energy redistribution
    occurs between molecular and cavity modes according to the new coupling strength.
    
    Parameters
    ----------
    target_value : float
        The target value :math:`g_{\text{target}}` to switch to at switch_time
    switch_time_ps : float
        The switch time :math:`t_{\text{switch}}` in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing in adaptive timestep simulations
        
    Attributes
    ----------
    target_value : float
        Target coupling value after switching
    switch_time_ps : float
        Switch time in picoseconds
    current_value : float
        Current value of the variant (0.0 before switch, target_value after)
    has_switched : bool
        Whether the switching has occurred
        
    Examples
    --------
    **Basic step variant:**
    
    >>> from hoomd.cavitymd.variants import StepVariant
    >>> from hoomd.cavitymd.analysis import ElapsedTimeTracker
    >>> 
    >>> # Create time tracker
    >>> time_tracker = ElapsedTimeTracker(sim, runtime_ps=10.0)
    >>> 
    >>> # Create step variant: switch at 1 ps to coupling 0.001
    >>> coupling_variant = StepVariant(
    ...     target_value=0.001,
    ...     switch_time_ps=1.0,
    ...     time_tracker=time_tracker
    ... )
    >>> 
    >>> # Use in cavity force
    >>> cavity_force = CavityForce(
    ...     kvector=[0, 0, 1],
    ...     couplstr=coupling_variant,
    ...     omegac=0.00913  # 2000 cm⁻¹
    ... )
    
    **Time-varying dissipation:**
    
    >>> # Also switch dissipation simultaneously
    >>> dissipation_variant = StepVariant(
    ...     target_value=0.0001,  # Add damping after switch
    ...     switch_time_ps=1.0,
    ...     time_tracker=time_tracker
    ... )
    >>> 
    >>> cavity_force = CavityForce(
    ...     kvector=[0, 0, 1],
    ...     couplstr=coupling_variant,
    ...     omegac=0.00913,
    ...     dissipation=dissipation_variant
    ... )
    
    Notes
    -----
    - Requires ElapsedTimeTracker for accurate timing in adaptive timestep simulations
    - The switch is instantaneous (discontinuous) which may cause transients
    - Compatible with both coupling strength and dissipation parameters
    - Switching preserves total system energy but redistributes it between modes
    
    See Also
    --------
    ConstantVariant : For time-independent parameters
    hoomd.cavitymd.analysis.ElapsedTimeTracker : For time tracking
    hoomd.cavitymd.forces.CavityForce : For cavity force with time-varying coupling
    
    References
    ----------
    For the theoretical framework of time-varying coupling, see the theory
    documentation section on "Time-Varying Coupling Theory".
    """
    
    def __init__(self, 
                 target_value: float, 
                 switch_time_ps: float, 
                 time_tracker: ElapsedTimeTracker) -> None:
        """
        Initialize the step variant.
        
        Parameters
        ----------
        target_value : float
            Value to switch to at switch_time_ps
        switch_time_ps : float  
            Switch time in picoseconds
        time_tracker : ElapsedTimeTracker
            Time tracker for accurate timing
        """
        hoomd.variant.Variant.__init__(self)
        self.target_value: float = float(target_value)
        self.switch_time_ps: float = float(switch_time_ps)
        self.time_tracker: ElapsedTimeTracker = time_tracker
        
        # Store for debugging/logging
        self.current_value: float = 0.0
        self._has_switched: bool = False
        
        print(f"StepVariant initialized:")
        print(f"  Target value: {self.target_value}")
        print(f"  Switch time: {self.switch_time_ps} ps")
        print(f"  Using ElapsedTimeTracker for accurate timing")
    
    def __call__(self, timestep: int) -> float:
        """
        Evaluate the step function at the given timestep.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep (unused, time comes from time_tracker)
            
        Returns
        -------
        float
            Current variant value (0.0 before switch, target_value after)
        """
        # Get current elapsed time from time tracker
        current_time_ps = self.time_tracker.elapsed_time
        
        if current_time_ps >= self.switch_time_ps:
            if not self._has_switched:
                print(f"StepVariant: Switching at t = {current_time_ps:.6f} ps "
                      f"(target: {self.switch_time_ps:.6f} ps)")
                self._has_switched = True
            self.current_value = self.target_value
            return self.target_value
        else:
            self.current_value = 0.0
            return 0.0
    
    def _min(self) -> float:
        """Return the minimum value of the variant."""
        return 0.0
    
    def _max(self) -> float:
        """Return the maximum value of the variant."""
        return self.target_value
    
    @property
    def has_switched(self) -> bool:
        """Whether the variant has switched to the target value."""
        return self._has_switched


class ConstantVariant(hoomd.variant.Variant):
    """
    Simple constant variant for backward compatibility.
    
    Provides a time-independent parameter value, equivalent to using a 
    plain float but maintaining consistency with the variant interface.
    
    .. math::
        
        g(t) = g_{\text{constant}}
    
    This is useful for maintaining uniform interfaces when some parameters
    are time-varying and others are constant.
    
    Parameters
    ----------
    value : float
        The constant value to return
        
    Examples
    --------
    >>> from hoomd.cavitymd.variants import ConstantVariant
    >>> constant_coupling = ConstantVariant(0.001)
    >>> print(constant_coupling(1000))  # Always returns 0.001
    0.001
    
    Notes
    -----
    While you can use plain floats for constant parameters, ConstantVariant
    provides explicit typing and consistent interface with time-varying variants.
    """
    
    def __init__(self, value: float) -> None:
        """
        Initialize the constant variant.
        
        Parameters
        ----------
        value : float
            The constant value to return
        """
        hoomd.variant.Variant.__init__(self)
        self.value: float = float(value)
    
    def __call__(self, timestep: int) -> float:
        """
        Return the constant value.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep (ignored)
            
        Returns
        -------
        float
            The constant value
        """
        return self.value
    
    def _min(self) -> float:
        """Return the minimum value (constant)."""
        return self.value
    
    def _max(self) -> float:
        """Return the maximum value (constant)."""
        return self.value 