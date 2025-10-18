# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Custom variant classes for cavity molecular dynamics."""

from typing import Union, Optional
import hoomd.variant
import numpy as np
from ..utils import PhysicalConstants
from ..analysis import ElapsedTimeTracker


class StepVariant(hoomd.variant.Variant):
    r"""
    Step function variant for instantaneous coupling switching with optional exponential decay.
    
    Implements a discontinuous step function that switches from 0.0 to a target
    value at a specified switch time, with optional exponential decay after switching:
    
    **Without decay:**
    
    .. math::
        
        g(t) = \begin{cases}
        0 & \text{if } t < t_{\text{switch}} \\
        g_{\text{target}} & \text{if } t \geq t_{\text{switch}}
        \end{cases}
    
    **With decay:**
    
    .. math::
        
        g(t) = \begin{cases}
        0 & \text{if } t < t_{\text{switch}} \\
        g_{\text{target}} \exp\left(-\frac{t - t_{\text{switch}}}{\tau_{\text{decay}}}\right) & \text{if } t \geq t_{\text{switch}}
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
    - Dynamic switching protocols with decay
    - Non-equilibrium cavity-coupling studies
    - Cavity mode decoherence and damping effects
    
    **Energy Conservation:**
    
    During switching, total energy is conserved. However, energy redistribution
    occurs between molecular and cavity modes according to the new coupling strength.
    With decay, energy is gradually dissipated from the cavity mode.
    
    Parameters
    ----------
    target_value : float
        The target value :math:`g_{\text{target}}` to switch to at switch_time
    switch_time_ps : float
        The switch time :math:`t_{\text{switch}}` in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing in adaptive timestep simulations
    decay_time_constant_ps : float, optional
        The exponential decay time constant :math:`\tau_{\text{decay}}` in picoseconds.
        If None, no decay occurs (standard step function). If provided, the coupling
        decays exponentially after switching with time constant τ_decay.
        
    Attributes
    ----------
    target_value : float
        Target coupling value after switching
    switch_time_ps : float
        Switch time in picoseconds
    decay_time_constant_ps : float or None
        Decay time constant in picoseconds, or None for no decay
    current_value : float
        Current value of the variant (0.0 before switch, decaying after switch)
    has_switched : bool
        Whether the switching has occurred
        
    Examples
    --------
    **Basic step variant (no decay):**
    
    >>> from hoomd.cavitymd.coupling import StepVariant
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
    
    **Step variant with exponential decay:**
    
    >>> # Switch at 1 ps, then decay with 2 ps time constant
    >>> decaying_variant = StepVariant(
    ...     target_value=0.001,
    ...     switch_time_ps=1.0,
    ...     time_tracker=time_tracker,
    ...     decay_time_constant_ps=2.0  # Decay with τ = 2 ps
    ... )
    >>> 
    >>> # Use in cavity force
    >>> cavity_force = CavityForce(
    ...     kvector=[0, 0, 1],
    ...     couplstr=decaying_variant,
    ...     omegac=0.00913  # 2000 cm⁻¹
    ... )
    
    **Time-varying dissipation with decay:**
    
    >>> # Also switch dissipation simultaneously with different decay
    >>> dissipation_variant = StepVariant(
    ...     target_value=0.0001,  # Add damping after switch
    ...     switch_time_ps=1.0,
    ...     time_tracker=time_tracker,
    ...     decay_time_constant_ps=5.0  # Longer decay for dissipation
    ... )
    >>> 
    >>> cavity_force = CavityForce(
    ...     kvector=[0, 0, 1],
    ...     couplstr=decaying_variant,
    ...     omegac=0.00913,
    ...     dissipation=dissipation_variant
    ... )
    
    Notes
    -----
    - Requires ElapsedTimeTracker for accurate timing in adaptive timestep simulations
    - The switch is instantaneous (discontinuous) which may cause transients
    - Exponential decay is continuous after switching
    - Compatible with both coupling strength and dissipation parameters
    - Decay time constant should be chosen based on physical cavity decoherence times
    - For adaptive timestepping, timing accuracy depends on ElapsedTimeTracker
    
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
                 time_tracker: ElapsedTimeTracker,
                 decay_time_constant_ps: Optional[float] = None,
                 turn_off_time_ps: Optional[float] = None,
                 simulation=None) -> None:
        """
        Initialize the step variant with optional exponential decay.
        
        Parameters
        ----------
        target_value : float
            Value to switch to at switch_time_ps
        switch_time_ps : float  
            Switch time in picoseconds
        time_tracker : ElapsedTimeTracker
            Time tracker for accurate timing
        decay_time_constant_ps : float, optional
            Exponential decay time constant in picoseconds. If None, no decay occurs.
        turn_off_time_ps : float, optional
            Time to turn off the coupling (force to zero) in picoseconds. Default: None (never turn off)
        """
        hoomd.variant.Variant.__init__(self)
        self.target_value: float = float(target_value)
        self.switch_time_ps: float = float(switch_time_ps)
        self.time_tracker: ElapsedTimeTracker = time_tracker
        self.decay_time_constant_ps: Optional[float] = (
            float(decay_time_constant_ps) if decay_time_constant_ps is not None else None
        )
        self.turn_off_time_ps: Optional[float] = (
            float(turn_off_time_ps) if turn_off_time_ps is not None else None
        )
        self.simulation = simulation
        
        # Store for debugging/logging
        self.current_value: float = 0.0
        self._has_switched: bool = False
        self._has_stopped: bool = False
        
        print(f"StepVariant initialized:")
        print(f"  Target value: {self.target_value}")
        print(f"  Switch time: {self.switch_time_ps} ps")
        if self.decay_time_constant_ps is not None:
            print(f"  Decay time constant: {self.decay_time_constant_ps} ps")
        else:
            print(f"  No decay (standard step function)")
        if self.turn_off_time_ps is not None:
            print(f"  Turn-off time: {self.turn_off_time_ps} ps")
        else:
            print(f"  No turn-off time (runs indefinitely)")
        print(f"  Using ElapsedTimeTracker for accurate timing")
    
    def __call__(self, timestep: int) -> float:
        """
        Evaluate the step function with optional decay at the given timestep.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep (unused, time comes from time_tracker)
            
        Returns
        -------
        float
            Current variant value (0.0 before switch, target_value or decaying after)
        """
        # Check for auto-stop signal first
        if self._check_auto_stop_signal():
            self.current_value = 0.0
            return 0.0
        
        # Get current elapsed time from time tracker
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should turn off
        if self.turn_off_time_ps is not None and current_time_ps >= self.turn_off_time_ps:
            if not self._has_stopped:
                print(f"StepVariant: Turned off at t = {current_time_ps:.6f} ps")
                self._has_stopped = True
            self.current_value = 0.0
            return 0.0
        
        if current_time_ps >= self.switch_time_ps:
            if not self._has_switched:
                print(f"StepVariant: Switching at t = {current_time_ps:.6f} ps "
                      f"(target: {self.switch_time_ps:.6f} ps)")
                self._has_switched = True
            
            if self.decay_time_constant_ps is not None:
                # Exponential decay: ε₀ * exp(-(t - t_switch)/τ_decay)
                time_since_switch = current_time_ps - self.switch_time_ps
                decay_factor = np.exp(-time_since_switch / self.decay_time_constant_ps)
                self.current_value = self.target_value * decay_factor
                return self.current_value
            else:
                # Standard step function (no decay)
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
    
    @property
    def has_stopped(self) -> bool:
        """Whether the variant has been turned off."""
        return self._has_stopped
    
    def _check_auto_stop_signal(self) -> bool:
        """
        Check if the auto-stop coupling signal has been set.
        
        Returns
        -------
        bool
            True if coupling should be disabled (return 0), False if normal operation
        """
        if self.simulation is None:
            return False
            
        try:
            # Check if auto-stop signal has been set
            if hasattr(self.simulation, 'auto_stop_coupling_signal'):
                if self.simulation.auto_stop_coupling_signal:
                    return True
                    
            # Also check HOOMD simulation object if available
            hoomd_sim = getattr(self.simulation, 'sim', self.simulation)
            if hasattr(hoomd_sim, 'auto_stop_coupling_signal'):
                if hoomd_sim.auto_stop_coupling_signal:
                    return True
                    
            # Check user_data for the signal
            if hasattr(self.simulation, 'state') and hasattr(self.simulation.state, 'user_data'):
                # Skip user_data check - not available in this HOOMD version
                pass  # if self.simulation.state.user_data.get('auto_stop_coupling_disabled', False):
                return True
                    
        except Exception:
            pass  # Continue with normal operation if any check fails
            
        return False


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
    >>> from hoomd.cavitymd.coupling import ConstantVariant
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


