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


class PeriodicVariant(hoomd.variant.Variant):
    r"""
    Periodic sinusoidal variant for oscillating coupling constant protocols.
    
    Implements a continuous sinusoidal function that oscillates the coupling
    constant with a specified period and amplitude:
    
    .. math::
        
        g(t) = A \frac{1 + \sin\left(\frac{2\pi t}{T} + \phi\right)}{2}
    
    where:
    - :math:`A` is the amplitude (maximum coupling strength)
    - :math:`T` is the period in picoseconds  
    - :math:`\phi` is the phase offset in radians
    - :math:`t` is the elapsed time in picoseconds
    
    This creates a time-dependent coupling strength that oscillates between 0 and A:
    
    .. math::
        
        \tilde{\varepsilon}_{0,\lambda}(t) = g(t) \tilde{\varepsilon}_{0,\lambda}^{(0)}
    
    where :math:`\tilde{\varepsilon}_{0,\lambda}^{(0)}` is the base coupling strength scale.
    
    **Physical Interpretation:**
    
    The periodic variant simulates experimental scenarios with oscillating 
    cavity coupling, such as:
    
    - Pump-probe experiments with periodic driving
    - Dynamic cavity tuning with sinusoidal modulation
    - AC cavity-molecule coupling protocols  
    - Periodic Rabi oscillations in cavity QED
    - Studying dynamic phase transitions under periodic driving
    - Cavity mode frequency modulation effects
    
    **Energy Conservation:**
    
    The periodic driving adds energy to the system at each oscillation cycle.
    This represents external driving of the cavity mode, which can lead to
    heating unless balanced by dissipation. Energy input rate depends on
    the amplitude and frequency of oscillation.
    
    Parameters
    ----------
    amplitude : float
        The amplitude :math:`A` of the oscillation (maximum coupling strength)
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing in adaptive timestep simulations
    period_ps : float, optional
        The oscillation period :math:`T` in picoseconds. Default: 1.0 ps
    phase_offset : float, optional
        Phase offset :math:`\phi` in radians. Default: 0.0 (starts at zero)
    start_time_ps : float, optional
        Time in ps to start the oscillation. Before this time, returns 0.0.
        Default: 0.0 (oscillates from simulation start)
    simulation : hoomd.Simulation, optional
        Reference to simulation object for auto-stop communication. Default: None
        
    Attributes
    ----------
    amplitude : float
        Oscillation amplitude (maximum coupling strength)
    period_ps : float
        Oscillation period in picoseconds
    phase_offset : float
        Phase offset in radians
    start_time_ps : float
        Start time for oscillation in picoseconds
    frequency_hz : float
        Oscillation frequency in Hz (derived from period)
    angular_frequency : float
        Angular frequency in rad/ps
    current_value : float
        Current value of the variant
        
    Examples
    --------
    **Basic periodic coupling (1 ps period):**
    
    >>> from hoomd.cavitymd.variants import PeriodicVariant
    >>> from hoomd.cavitymd.analysis import ElapsedTimeTracker
    >>> 
    >>> # Create time tracker
    >>> time_tracker = ElapsedTimeTracker(sim, runtime_ps=10.0)
    >>> 
    >>> # Create periodic variant: oscillate with amplitude 0.001, period 1 ps
    >>> periodic_coupling = PeriodicVariant(
    ...     amplitude=0.001,
    ...     time_tracker=time_tracker,
    ...     period_ps=1.0
    ... )
    
    **Periodic coupling with custom period and phase:**
    
    >>> # Oscillate with 2 ps period, starting at maximum (π/2 phase offset)
    >>> periodic_coupling = PeriodicVariant(
    ...     amplitude=0.001,
    ...     time_tracker=time_tracker,
    ...     period_ps=2.0,
    ...     phase_offset=np.pi/2  # Start at maximum
    ... )
    >>> 
    >>> # Use in cavity force
    >>> cavity_force = CavityForce(
    ...     kvector=[0, 0, 1],
    ...     couplstr=periodic_coupling,
    ...     omegac=0.00913  # 2000 cm⁻¹
    ... )
    
    **High-frequency oscillation (0.1 ps period):**
    
    >>> # Fast oscillation for studying high-frequency effects
    >>> fast_coupling = PeriodicVariant(
    ...     amplitude=0.0005,
    ...     time_tracker=time_tracker,
    ...     period_ps=0.1  # 10 THz frequency
    ... )
    
    **Delayed periodic coupling:**
    
    >>> # Start oscillation after 2 ps of equilibration
    >>> delayed_coupling = PeriodicVariant(
    ...     amplitude=0.001,
    ...     time_tracker=time_tracker,
    ...     period_ps=1.0,
    ...     start_time_ps=2.0  # Start oscillating at t=2 ps
    ... )
    
    Notes
    -----
    - Requires ElapsedTimeTracker for accurate timing in adaptive timestep simulations
    - The oscillation is continuous and smooth (no discontinuities)
    - Energy is continuously added to the system during oscillation
    - Compatible with both coupling strength and dissipation parameters
    - Period should be chosen based on relevant physical timescales
    - For adaptive timestepping, timing accuracy depends on ElapsedTimeTracker
    - Very high frequencies may be limited by timestep resolution
    - Auto-stop functionality can disable coupling when temperatures converge
    
    See Also
    --------
    StepVariant : For discontinuous switching protocols
    ConstantVariant : For time-independent parameters
    hoomd.cavitymd.analysis.ElapsedTimeTracker : For time tracking
    hoomd.cavitymd.forces.CavityForce : For cavity force with time-varying coupling
    
    References
    ----------
    For the theoretical framework of time-varying coupling, see the theory
    documentation section on "Time-Varying Coupling Theory" and periodic
    driving in cavity QED systems.
    """
    
    def __init__(self, 
                 amplitude: float, 
                 time_tracker: ElapsedTimeTracker,
                 period_ps: float = 1.0,
                 phase_offset: float = 0.0,
                 start_time_ps: float = 0.0,
                 stop_time_ps: Optional[float] = None,
                 simulation=None) -> None:
        """
        Initialize the periodic variant.
        
        Parameters
        ----------
        amplitude : float
            Oscillation amplitude (maximum coupling strength)
        time_tracker : ElapsedTimeTracker
            Time tracker for accurate timing
        period_ps : float, optional
            Oscillation period in picoseconds. Default: 1.0 ps
        phase_offset : float, optional
            Phase offset in radians. Default: 0.0
        start_time_ps : float, optional
            Start time for oscillation in picoseconds. Default: 0.0
        stop_time_ps : float, optional
            Time to stop the oscillation in picoseconds. Default: None (never stop)
        simulation : hoomd.Simulation, optional
            Reference to simulation object for auto-stop communication. Default: None
        """
        hoomd.variant.Variant.__init__(self)
        self.amplitude: float = float(amplitude)
        self.time_tracker: ElapsedTimeTracker = time_tracker
        self.period_ps: float = float(period_ps)
        self.phase_offset: float = float(phase_offset)
        self.start_time_ps: float = float(start_time_ps)
        self.stop_time_ps: Optional[float] = (
            float(stop_time_ps) if stop_time_ps is not None else None
        )
        self.simulation = simulation
        
        # Derived quantities
        self.frequency_hz: float = 1.0 / (self.period_ps * 1e-12)  # Convert ps to Hz
        self.angular_frequency: float = 2.0 * np.pi / self.period_ps  # rad/ps
        
        # State variables
        self.current_value: float = 0.0
        self._has_started: bool = False
        self._has_stopped: bool = False
        
        print(f"PeriodicVariant initialized:")
        print(f"  Amplitude: {self.amplitude}")
        print(f"  Period: {self.period_ps} ps")
        print(f"  Frequency: {self.frequency_hz:.3e} Hz")
        print(f"  Phase offset: {self.phase_offset:.3f} rad")
        print(f"  Start time: {self.start_time_ps} ps")
        if self.stop_time_ps is not None:
            print(f"  Stop time: {self.stop_time_ps} ps")
        else:
            print(f"  No stop time (runs indefinitely)")
        print(f"  Using ElapsedTimeTracker for accurate timing")
    
    def __call__(self, timestep: int) -> float:
        """
        Evaluate the periodic function at the given timestep.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep (unused, time comes from time_tracker)
            
        Returns
        -------
        float
            Current variant value (0.0 before start, oscillating after)
        """
        # Check for auto-stop signal first
        if self._check_auto_stop_signal():
            self.current_value = 0.0
            return 0.0
        
        # Get current elapsed time from time tracker
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should stop
        if self.stop_time_ps is not None and current_time_ps >= self.stop_time_ps:
            if not self._has_stopped:
                print(f"PeriodicVariant: Stopped at t = {current_time_ps:.6f} ps")
                self._has_stopped = True
            self.current_value = 0.0
            return 0.0
        
        if current_time_ps >= self.start_time_ps:
            if not self._has_started:
                print(f"PeriodicVariant: Oscillation started at t = {current_time_ps:.6f} ps")
                self._has_started = True
            
            # Calculate time since oscillation started
            time_since_start = current_time_ps - self.start_time_ps
            
            # Periodic function: A * (1 + sin(2π*t/T + φ)) / 2
            # This ensures coupling varies from 0 to amplitude (never negative)
            phase = self.angular_frequency * time_since_start + self.phase_offset
            self.current_value = self.amplitude * (1 + np.sin(phase)) / 2
            
            return self.current_value
        else:
            self.current_value = 0.0
            return 0.0
    
    def _check_auto_stop_signal(self) -> bool:
        """
        Check if the auto-stop coupling signal has been set by empirical feedback.
        
        Returns
        -------
        bool
            True if coupling should be disabled (return 0), False if normal operation
        """
        if self.simulation is None:
            return False
            
        try:
            # Check if auto-stop signal has been set by empirical feedback
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
                #     return True
                    
        except Exception:
            pass  # Continue with normal operation if any check fails
            
        return False
    
    def _min(self) -> float:
        """Return the minimum value of the variant."""
        return 0.0
    
    def _max(self) -> float:
        """Return the maximum value of the variant."""
        return self.amplitude
    
    @property
    def has_started(self) -> bool:
        """Whether the periodic oscillation has started."""
        return self._has_started
    
    @property
    def has_stopped(self) -> bool:
        """Whether the periodic oscillation has been stopped."""
        return self._has_stopped


class ExponentialDecayVariant(hoomd.variant.Variant):
    r"""
    Exponential decay variant for time-dependent coupling with optional turn-off.
    
    Implements exponential decay starting from an initial value at time zero,
    with optional turn-on/turn-off timing control:
    
    .. math::
        
        g(t) = \begin{cases}
        0 & \text{if } t < t_{\text{on}} \\
        A \exp\left(-\frac{t - t_{\text{on}}}{\tau}\right) & \text{if } t_{\text{on}} \leq t < t_{\text{off}} \\
        0 & \text{if } t \geq t_{\text{off}}
        \end{cases}
    
    where:
    - :math:`A` is the initial amplitude
    - :math:`\tau` is the decay time constant in picoseconds
    - :math:`t_{\text{on}}` is the turn-on time (default: 0.0)
    - :math:`t_{\text{off}}` is the turn-off time (default: never)
    
    **Physical Interpretation:**
    
    This variant simulates physical processes with exponential relaxation:
    
    - Cavity mode decoherence and damping
    - Exponential energy dissipation
    - Time-dependent coupling strength decay
    - Pump-probe recovery dynamics
    - Temperature-dependent coupling relaxation
    
    **Energy Conservation:**
    
    The exponential decay naturally conserves energy by gradually reducing
    the coupling strength. Energy initially stored in the cavity mode is
    slowly released back to the molecular system.
    
    Parameters
    ----------
    amplitude : float
        Initial amplitude :math:`A` of the coupling
    decay_time_constant_ps : float
        Exponential decay time constant :math:`\tau` in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing in adaptive timestep simulations
    turn_on_time_ps : float, optional
        Time to start the decay in picoseconds. Default: 0.0
    turn_off_time_ps : float, optional
        Time to stop the decay (force to zero) in picoseconds. Default: None (never turn off)
        
    Attributes
    ----------
    amplitude : float
        Initial amplitude of the decay
    decay_time_constant_ps : float
        Decay time constant in picoseconds
    turn_on_time_ps : float
        Turn-on time in picoseconds
    turn_off_time_ps : float or None
        Turn-off time in picoseconds, or None for no turn-off
    current_value : float
        Current value of the variant
    has_started : bool
        Whether the decay has started
    has_stopped : bool
        Whether the decay has been turned off
        
    Examples
    --------
    **Basic exponential decay (starts at t=0):**
    
    >>> from hoomd.cavitymd.variants import ExponentialDecayVariant
    >>> from hoomd.cavitymd.analysis import ElapsedTimeTracker
    >>> 
    >>> # Create time tracker
    >>> time_tracker = ElapsedTimeTracker(sim, runtime_ps=10.0)
    >>> 
    >>> # Create exponential decay: amplitude 0.001, decay 2 ps
    >>> decay_coupling = ExponentialDecayVariant(
    ...     amplitude=0.001,
    ...     decay_time_constant_ps=2.0,
    ...     time_tracker=time_tracker
    ... )
    
    **Delayed exponential decay:**
    
    >>> # Start decay at 1 ps
    >>> delayed_decay = ExponentialDecayVariant(
    ...     amplitude=0.001,
    ...     decay_time_constant_ps=2.0,
    ...     time_tracker=time_tracker,
    ...     turn_on_time_ps=1.0
    ... )
    
    **Finite-duration exponential decay:**
    
    >>> # Decay from 1 ps to 5 ps, then stop
    >>> finite_decay = ExponentialDecayVariant(
    ...     amplitude=0.001,
    ...     decay_time_constant_ps=2.0,
    ...     time_tracker=time_tracker,
    ...     turn_on_time_ps=1.0,
    ...     turn_off_time_ps=5.0
    ... )
    
    Notes
    -----
    - Requires ElapsedTimeTracker for accurate timing
    - The decay is continuous and smooth
    - Compatible with both coupling strength and dissipation parameters
    - Decay time constant should be based on physical decoherence times
    - Turn-off is instantaneous (discontinuous)
    """
    
    def __init__(self,
                 amplitude: float,
                 decay_time_constant_ps: float,
                 time_tracker: ElapsedTimeTracker,
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None) -> None:
        """
        Initialize the exponential decay variant.
        
        Parameters
        ----------
        amplitude : float
            Initial amplitude of the decay
        decay_time_constant_ps : float
            Exponential decay time constant in picoseconds  
        time_tracker : ElapsedTimeTracker
            Time tracker for accurate timing
        turn_on_time_ps : float, optional
            Time to start the decay in picoseconds. Default: 0.0
        turn_off_time_ps : float, optional
            Time to stop the decay in picoseconds. Default: None (never turn off)
        """
        hoomd.variant.Variant.__init__(self)
        self.amplitude: float = float(amplitude)
        self.decay_time_constant_ps: float = float(decay_time_constant_ps)
        self.time_tracker: ElapsedTimeTracker = time_tracker
        self.turn_on_time_ps: float = float(turn_on_time_ps)
        self.turn_off_time_ps: Optional[float] = (
            float(turn_off_time_ps) if turn_off_time_ps is not None else None
        )
        
        # State variables
        self.current_value: float = 0.0
        self._has_started: bool = False
        self._has_stopped: bool = False
        
        print(f"ExponentialDecayVariant initialized:")
        print(f"  Amplitude: {self.amplitude}")
        print(f"  Decay time constant: {self.decay_time_constant_ps} ps")
        print(f"  Turn-on time: {self.turn_on_time_ps} ps")
        if self.turn_off_time_ps is not None:
            print(f"  Turn-off time: {self.turn_off_time_ps} ps")
        else:
            print(f"  No turn-off time (decays indefinitely)")
    
    def __call__(self, timestep: int) -> float:
        """
        Evaluate the exponential decay function at the given timestep.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep (unused, time comes from time_tracker)
            
        Returns
        -------
        float
            Current variant value
        """
        # Get current elapsed time from time tracker
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should turn off
        if self.turn_off_time_ps is not None and current_time_ps >= self.turn_off_time_ps:
            if not self._has_stopped:
                print(f"ExponentialDecayVariant: Turned off at t = {current_time_ps:.6f} ps")
                self._has_stopped = True
            self.current_value = 0.0
            return 0.0
        
        # Check if we should start
        if current_time_ps >= self.turn_on_time_ps:
            if not self._has_started:
                print(f"ExponentialDecayVariant: Started at t = {current_time_ps:.6f} ps")
                self._has_started = True
            
            # Calculate time since start
            time_since_start = current_time_ps - self.turn_on_time_ps
            
            # Exponential decay: A * exp(-t/τ)
            self.current_value = self.amplitude * np.exp(-time_since_start / self.decay_time_constant_ps)
            return self.current_value
        else:
            self.current_value = 0.0
            return 0.0
    
    def _min(self) -> float:
        """Return the minimum value of the variant."""
        return 0.0
    
    def _max(self) -> float:
        """Return the maximum value of the variant."""
        return self.amplitude
    
    @property
    def has_started(self) -> bool:
        """Whether the exponential decay has started."""
        return self._has_started
    
    @property
    def has_stopped(self) -> bool:
        """Whether the exponential decay has been turned off."""
        return self._has_stopped


class SquareWaveVariant(hoomd.variant.Variant):
    r"""
    Square wave variant for periodic on/off coupling protocols.
    
    Implements a square wave that switches between 0 and a target amplitude
    with configurable duty cycle and timing:
    
    .. math::
        
        g(t) = \begin{cases}
        0 & \text{if } t < t_{\text{start}} \\
        A & \text{if } t_{\text{start}} \leq t < t_{\text{stop}} \text{ and } \text{wave is high} \\
        0 & \text{if } t_{\text{start}} \leq t < t_{\text{stop}} \text{ and } \text{wave is low} \\
        0 & \text{if } t \geq t_{\text{stop}}
        \end{cases}
    
    The square wave pattern within the active period is:
    
    .. math::
        
        \text{phase}(t) = \frac{2\pi (t - t_{\text{start}})}{T} + \phi
        
        g_{\text{wave}}(t) = \begin{cases}
        A & \text{if } \sin(\text{phase}(t)) \geq 0 \text{ and duty cycle conditions met} \\
        0 & \text{otherwise}
        \end{cases}
    
    where:
    - :math:`A` is the amplitude when "on"
    - :math:`T` is the period in picoseconds
    - :math:`\phi` is the phase offset in radians
    - Duty cycle controls the fraction of time the wave is "high"
    
    **Physical Interpretation:**
    
    Square wave coupling simulates experimental scenarios with discrete on/off control:
    
    - Pulsed laser experiments with defined on/off periods
    - Switching cavity resonance in and out of molecular resonance
    - Digital control protocols for cavity coupling
    - Interrupted coupling experiments
    - Pump-probe with controlled dark periods
    
    **Energy Conservation:**
    
    Energy input occurs only during "on" periods. The discrete switching
    can cause transients, but total energy is conserved across switching events.
    
    Parameters
    ----------
    amplitude : float
        Amplitude :math:`A` when the square wave is "high"
    period_ps : float
        Period :math:`T` of the square wave in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing in adaptive timestep simulations
    duty_cycle : float, optional
        Fraction of each period that the wave is "high" (0.0 to 1.0). Default: 0.5
    phase_offset : float, optional
        Phase offset :math:`\phi` in radians. Default: 0.0
    start_time_ps : float, optional
        Time to start the square wave in picoseconds. Default: 0.0
    stop_time_ps : float, optional
        Time to stop the square wave in picoseconds. Default: None (never stop)
        
    Attributes
    ----------
    amplitude : float
        Square wave amplitude
    period_ps : float
        Period in picoseconds
    duty_cycle : float
        Duty cycle (fraction of time "high")
    phase_offset : float
        Phase offset in radians
    start_time_ps : float
        Start time in picoseconds
    stop_time_ps : float or None
        Stop time in picoseconds, or None for no stop
    current_value : float
        Current value of the variant
    is_high : bool
        Whether the square wave is currently "high"
    has_started : bool
        Whether the square wave has started
    has_stopped : bool
        Whether the square wave has been stopped
        
    Examples
    --------
    **Basic square wave (50% duty cycle):**
    
    >>> from hoomd.cavitymd.variants import SquareWaveVariant
    >>> from hoomd.cavitymd.analysis import ElapsedTimeTracker
    >>> 
    >>> # Create time tracker
    >>> time_tracker = ElapsedTimeTracker(sim, runtime_ps=10.0)
    >>> 
    >>> # Create square wave: amplitude 0.001, period 2 ps, 50% duty cycle
    >>> square_coupling = SquareWaveVariant(
    ...     amplitude=0.001,
    ...     period_ps=2.0,
    ...     time_tracker=time_tracker
    ... )
    
    **Custom duty cycle (75% on, 25% off):**
    
    >>> # Square wave with 75% duty cycle
    >>> high_duty_coupling = SquareWaveVariant(
    ...     amplitude=0.001,
    ...     period_ps=2.0,
    ...     time_tracker=time_tracker,
    ...     duty_cycle=0.75
    ... )
    
    **Delayed finite-duration square wave:**
    
    >>> # Start at 1 ps, stop at 10 ps
    >>> finite_square = SquareWaveVariant(
    ...     amplitude=0.001,
    ...     period_ps=1.0,
    ...     time_tracker=time_tracker,
    ...     start_time_ps=1.0,
    ...     stop_time_ps=10.0,
    ...     duty_cycle=0.3  # 30% on, 70% off
    ... )
    
    Notes
    -----
    - Requires ElapsedTimeTracker for accurate timing
    - The switching is instantaneous (discontinuous)
    - Can cause transients due to sudden coupling changes
    - Compatible with both coupling strength and dissipation parameters
    - Period should be chosen based on relevant physical timescales
    - Duty cycle of 0.0 means always off, 1.0 means always on
    """
    
    def __init__(self,
                 amplitude: float,
                 period_ps: float,
                 time_tracker: ElapsedTimeTracker,
                 duty_cycle: float = 0.5,
                 phase_offset: float = 0.0,
                 start_time_ps: float = 0.0,
                 stop_time_ps: Optional[float] = None,
                 simulation=None) -> None:
        """
        Initialize the square wave variant.
        
        Parameters
        ----------
        amplitude : float
            Square wave amplitude
        period_ps : float
            Period of the square wave in picoseconds
        time_tracker : ElapsedTimeTracker
            Time tracker for accurate timing
        duty_cycle : float, optional
            Fraction of each period that the wave is "high". Default: 0.5
        phase_offset : float, optional
            Phase offset in radians. Default: 0.0
        start_time_ps : float, optional
            Time to start the square wave in picoseconds. Default: 0.0
        stop_time_ps : float, optional
            Time to stop the square wave in picoseconds. Default: None
        """
        hoomd.variant.Variant.__init__(self)
        self.amplitude: float = float(amplitude)
        self.period_ps: float = float(period_ps)
        self.time_tracker: ElapsedTimeTracker = time_tracker
        self.duty_cycle: float = float(duty_cycle)
        self.phase_offset: float = float(phase_offset)
        self.start_time_ps: float = float(start_time_ps)
        self.stop_time_ps: Optional[float] = (
            float(stop_time_ps) if stop_time_ps is not None else None
        )
        self.simulation = simulation
        
        # Validate duty cycle
        if not 0.0 <= self.duty_cycle <= 1.0:
            raise ValueError(f"Duty cycle must be between 0.0 and 1.0, got {self.duty_cycle}")
        
        # Derived quantities
        self.angular_frequency: float = 2.0 * np.pi / self.period_ps  # rad/ps
        
        # State variables
        self.current_value: float = 0.0
        self.is_high: bool = False
        self._has_started: bool = False
        self._has_stopped: bool = False
        
        print(f"SquareWaveVariant initialized:")
        print(f"  Amplitude: {self.amplitude}")
        print(f"  Period: {self.period_ps} ps")
        print(f"  Duty cycle: {self.duty_cycle:.1%}")
        print(f"  Phase offset: {self.phase_offset:.3f} rad")
        print(f"  Start time: {self.start_time_ps} ps")
        if self.stop_time_ps is not None:
            print(f"  Stop time: {self.stop_time_ps} ps")
        else:
            print(f"  No stop time (runs indefinitely)")
    
    def __call__(self, timestep: int) -> float:
        """
        Evaluate the square wave function at the given timestep.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep (unused, time comes from time_tracker)
            
        Returns
        -------
        float
            Current variant value (0.0 or amplitude)
        """
        # Check for auto-stop signal first
        if self._check_auto_stop_signal():
            self.current_value = 0.0
            self.is_high = False
            return 0.0
        
        # Get current elapsed time from time tracker
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should stop
        if self.stop_time_ps is not None and current_time_ps >= self.stop_time_ps:
            if not self._has_stopped:
                print(f"SquareWaveVariant: Stopped at t = {current_time_ps:.6f} ps")
                self._has_stopped = True
            self.current_value = 0.0
            self.is_high = False
            return 0.0
        
        # Check if we should start
        if current_time_ps >= self.start_time_ps:
            if not self._has_started:
                print(f"SquareWaveVariant: Started at t = {current_time_ps:.6f} ps")
                self._has_started = True
            
            # Calculate time since start
            time_since_start = current_time_ps - self.start_time_ps
            
            # Calculate phase within the current period
            phase = self.angular_frequency * time_since_start + self.phase_offset
            phase_mod = phase % (2 * np.pi)  # Normalize to [0, 2π]
            
            # Determine if we're in the "high" part of the duty cycle
            # High period starts at phase_offset and lasts for duty_cycle * 2π
            high_phase_duration = self.duty_cycle * 2 * np.pi
            
            if phase_mod < high_phase_duration:
                self.is_high = True
                self.current_value = self.amplitude
            else:
                self.is_high = False
                self.current_value = 0.0
            
            return self.current_value
        else:
            self.current_value = 0.0
            self.is_high = False
            return 0.0
    
    def _min(self) -> float:
        """Return the minimum value of the variant."""
        return 0.0
    
    def _max(self) -> float:
        """Return the maximum value of the variant."""
        return self.amplitude
    
    @property
    def has_started(self) -> bool:
        """Whether the square wave has started."""
        return self._has_started
    
    @property
    def has_stopped(self) -> bool:
        """Whether the square wave has been stopped."""
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
                #     return True
                    
        except Exception:
            pass  # Continue with normal operation if any check fails
            
        return False


class DecayingSquareWaveVariant(hoomd.variant.Variant):
    r"""
    Decaying square wave variant for coupling protocols with amplitude reduction over time.
    
    Implements a square wave that switches between 0 and a decaying amplitude,
    where the amplitude decreases by a fixed percentage after each complete period:
    
    .. math::
        
        A_n = A_0 \times (1 - r)^n
    
    where:
    - :math:`A_0` is the initial amplitude
    - :math:`r` is the decay rate per period (0.0 to 1.0)
    - :math:`n` is the number of completed periods
    
    **Physical Interpretation:**
    
    This variant simulates experimental scenarios where cavity coupling strength 
    weakens over time due to decoherence, laser power decreases, or material 
    properties change under repeated excitation.
    
    Parameters
    ----------
    initial_amplitude : float
        Initial amplitude when the square wave is "high"
    period_ps : float
        Period of the square wave in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing in adaptive timestep simulations
    decay_rate_per_period : float, optional
        Fractional amplitude reduction after each period (0.0 to 1.0). Default: 0.0
    duty_cycle : float, optional
        Fraction of each period that the wave is "high" (0.0 to 1.0). Default: 0.5
    phase_offset : float, optional
        Phase offset in radians. Default: 0.0
    start_time_ps : float, optional
        Time to start the square wave in picoseconds. Default: 0.0
    stop_time_ps : float, optional
        Time to stop the square wave in picoseconds. Default: None (never stop)
    minimum_amplitude : float, optional
        Minimum amplitude below which the wave stops. Default: 1e-6
        
    Examples
    --------
    **Decaying square wave (10% decay per period):**
    
    >>> from hoomd.cavitymd.variants import DecayingSquareWaveVariant
    >>> from hoomd.cavitymd.analysis import ElapsedTimeTracker
    >>> 
    >>> # Create time tracker
    >>> time_tracker = ElapsedTimeTracker(sim, runtime_ps=20.0)
    >>> 
    >>> # Create decaying square wave: 10% amplitude reduction per period
    >>> decaying_square = DecayingSquareWaveVariant(
    ...     initial_amplitude=0.001,
    ...     period_ps=2.0,
    ...     time_tracker=time_tracker,
    ...     decay_rate_per_period=0.10  # 10% decay per period
    ... )
    """
    
    def __init__(self,
                 initial_amplitude: float,
                 period_ps: float,
                 time_tracker: ElapsedTimeTracker,
                 decay_rate_per_period: float = 0.0,
                 duty_cycle: float = 0.5,
                 phase_offset: float = 0.0,
                 start_time_ps: float = 0.0,
                 stop_time_ps: Optional[float] = None,
                 minimum_amplitude: float = 1e-6) -> None:
        
        hoomd.variant.Variant.__init__(self)
        
        # Core parameters
        self.initial_amplitude: float = float(initial_amplitude)
        self.current_amplitude: float = float(initial_amplitude)
        self.period_ps: float = float(period_ps)
        self.time_tracker: ElapsedTimeTracker = time_tracker
        self.decay_rate_per_period: float = float(decay_rate_per_period)
        self.duty_cycle: float = float(duty_cycle)
        self.phase_offset: float = float(phase_offset)
        self.start_time_ps: float = float(start_time_ps)
        self.stop_time_ps: Optional[float] = (
            float(stop_time_ps) if stop_time_ps is not None else None
        )
        self.minimum_amplitude: float = float(minimum_amplitude)
        
        # Validate parameters
        if not 0.0 <= self.decay_rate_per_period <= 1.0:
            raise ValueError(f"Decay rate must be between 0.0 and 1.0, got {self.decay_rate_per_period}")
        if not 0.0 <= self.duty_cycle <= 1.0:
            raise ValueError(f"Duty cycle must be between 0.0 and 1.0, got {self.duty_cycle}")
        if self.minimum_amplitude < 0.0:
            raise ValueError(f"Minimum amplitude must be >= 0.0, got {self.minimum_amplitude}")
        
        # Derived quantities
        self.angular_frequency: float = 2.0 * np.pi / self.period_ps  # rad/ps
        
        # State variables
        self.current_value: float = 0.0
        self.completed_periods: int = 0
        self.is_high: bool = False
        self._has_started: bool = False
        self._has_stopped: bool = False
        self._amplitude_too_small: bool = False
        
        print(f"DecayingSquareWaveVariant initialized:")
        print(f"  Initial amplitude: {self.initial_amplitude}")
        print(f"  Period: {self.period_ps} ps")
        print(f"  Decay rate: {self.decay_rate_per_period:.1%} per period")
        print(f"  Duty cycle: {self.duty_cycle:.1%}")
    
    def __call__(self, timestep: int) -> float:
        """Evaluate the decaying square wave function at the given timestep."""
        # Get current elapsed time from time tracker
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should stop due to time limit
        if self.stop_time_ps is not None and current_time_ps >= self.stop_time_ps:
            if not self._has_stopped:
                print(f"DecayingSquareWaveVariant: Stopped at t = {current_time_ps:.6f} ps (time limit)")
                self._has_stopped = True
            self.current_value = 0.0
            self.is_high = False
            return 0.0
        
        # Check if amplitude is too small to continue
        if self.current_amplitude < self.minimum_amplitude:
            if not self._amplitude_too_small:
                print(f"DecayingSquareWaveVariant: Stopped at t = {current_time_ps:.6f} ps "
                      f"(amplitude {self.current_amplitude:.2e} < minimum {self.minimum_amplitude:.2e})")
                self._amplitude_too_small = True
            self.current_value = 0.0
            self.is_high = False
            return 0.0
        
        # Check if we should start
        if current_time_ps >= self.start_time_ps:
            if not self._has_started:
                print(f"DecayingSquareWaveVariant: Started at t = {current_time_ps:.6f} ps")
                self._has_started = True
            
            # Calculate time since start
            time_since_start = current_time_ps - self.start_time_ps
            
            # Calculate which period we're in
            current_period = int(time_since_start / self.period_ps)
            
            # Update amplitude if we've completed additional periods
            if current_period > self.completed_periods:
                periods_completed = current_period - self.completed_periods
                for _ in range(periods_completed):
                    self.current_amplitude *= (1.0 - self.decay_rate_per_period)
                
                print(f"DecayingSquareWaveVariant: Period {current_period} completed, "
                      f"amplitude decayed to {self.current_amplitude:.6e}")
                self.completed_periods = current_period
                
                # Check if amplitude is now too small
                if self.current_amplitude < self.minimum_amplitude:
                    return self.__call__(timestep)  # Re-evaluate with new amplitude
            
            # Calculate phase within the current period
            phase = self.angular_frequency * time_since_start + self.phase_offset
            phase_mod = phase % (2 * np.pi)  # Normalize to [0, 2π]
            
            # Determine if we're in the "high" part of the duty cycle
            high_phase_duration = self.duty_cycle * 2 * np.pi
            
            if phase_mod < high_phase_duration:
                self.is_high = True
                self.current_value = self.current_amplitude
            else:
                self.is_high = False
                self.current_value = 0.0
            
            return self.current_value
        else:
            self.current_value = 0.0
            self.is_high = False
            return 0.0
    
    def _min(self) -> float:
        """Return the minimum value of the variant."""
        return 0.0
    
    def _max(self) -> float:
        """Return the maximum value of the variant."""
        return self.initial_amplitude
    
    @property
    def has_started(self) -> bool:
        """Whether the square wave has started."""
        return self._has_started
    
    @property
    def has_stopped(self) -> bool:
        """Whether the square wave has been stopped."""
        return self._has_stopped or self._amplitude_too_small
    
    @property
    def decay_factor(self) -> float:
        """Current decay factor relative to initial amplitude."""
        return self.current_amplitude / self.initial_amplitude if self.initial_amplitude > 0 else 0.0


class AdaptiveSquareWaveVariant(hoomd.variant.Variant):
    r"""
    Adaptive square wave variant that adjusts amplitude based on harmonic equipartition temperature.
    
    Implements an adaptive feedback control scheme that modifies the square wave amplitude
    based on the measured harmonic equipartition temperature at the end of each duty cycle:
    
    .. math::
        
        g_{\text{next}} = g_{\text{target}} \times \sqrt{\frac{T_{\text{target}}}{T_{\text{harmonic}}}}
    
    where:
    - :math:`g_{\text{target}}` is the reference coupling strength (from --coupling)
    - :math:`T_{\text{target}}` is the target temperature (from controller)
    - :math:`T_{\text{harmonic}}` is the measured harmonic equipartition temperature
    
    **Algorithm:**
    
    1. **Measure** :math:`T_{\text{harmonic}}` at the end of each duty cycle
    2. **Calculate** adaptive coupling: :math:`g_{\text{next}} = g_{\text{target}} \sqrt{T_{\text{target}} / T_{\text{harmonic}}}`
    3. **Apply** :math:`g_{\text{next}}` as the amplitude for the next square wave cycle
    4. **Repeat** this feedback loop throughout the simulation
    
    **Physical Interpretation:**
    
    This creates a self-regulating coupling control:
    - If :math:`T_{\text{harmonic}} > T_{\text{target}}`, then :math:`g_{\text{next}} < g_{\text{target}}` (reduce coupling)
    - If :math:`T_{\text{harmonic}} < T_{\text{target}}`, then :math:`g_{\text{next}} > g_{\text{target}}` (increase coupling)
    - The square root provides gentler adjustments than linear scaling
    - Uses harmonic temperature which reflects vibrational energy state
    
    **Temperature Sources:**
    
    The target temperature :math:`T_{\text{target}}` comes from either:
    - Quench controller: ``--quench-target-temperature``
    - Gradient descent controller: ``--gd-target-temperature``
    
    Parameters
    ----------
    target_coupling : float
        Reference coupling strength :math:`g_{\text{target}}` (from --coupling argument)
    target_temperature : float
        Target temperature :math:`T_{\text{target}}` in Kelvin (from controller)
    period_ps : float
        Period of the square wave in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing in adaptive timestep simulations
    temperature_tracker : object
        Temperature tracker object that provides harmonic_equipartition_temperature
    duty_cycle : float, optional
        Fraction of each period that the wave is "high" (0.0 to 1.0). Default: 0.5
    phase_offset : float, optional
        Phase offset in radians. Default: 0.0
    start_time_ps : float, optional
        Time to start the square wave in picoseconds. Default: 0.0
    stop_time_ps : float, optional
        Time to stop the square wave in picoseconds. Default: None (never stop)
    min_amplitude : float, optional
        Minimum allowed amplitude (safety limit). Default: 1e-8
    max_amplitude : float, optional
        Maximum allowed amplitude (safety limit). Default: 1e-1
        
    Examples
    --------
    **Adaptive square wave with quench controller:**
    
    >>> # Target temperature from quench controller
    >>> adaptive_coupling = AdaptiveSquareWaveVariant(
    ...     target_coupling=1e-3,
    ...     target_temperature=10.0,  # From --quench-target-temperature
    ...     period_ps=5.0,
    ...     time_tracker=time_tracker,
    ...     temperature_tracker=temp_tracker,
    ...     duty_cycle=0.75
    ... )
    
    **Adaptive square wave with gradient descent controller:**
    
    >>> # Target temperature from GD controller
    >>> adaptive_coupling = AdaptiveSquareWaveVariant(
    ...     target_coupling=1e-3,
    ...     target_temperature=25.0,  # From --gd-target-temperature
    ...     period_ps=8.0,
    ...     time_tracker=time_tracker,
    ...     temperature_tracker=temp_tracker
    ... )
    
    Notes
    -----
    - Requires a temperature tracker that provides harmonic_equipartition_temperature
    - Amplitude is updated only at the end of each duty cycle (not continuously)
    - Safety limits prevent extreme amplitude values
    - Compatible with both quench and gradient descent controllers
    - The square root scaling provides stable feedback control
    """
    
    def __init__(self,
                 target_coupling: float,
                 target_temperature: float,
                 period_ps: float,
                 time_tracker: ElapsedTimeTracker,
                 temperature_tracker,
                 duty_cycle: float = 0.5,
                 phase_offset: float = 0.0,
                 start_time_ps: float = 0.0,
                 stop_time_ps: Optional[float] = None,
                 min_amplitude: float = 1e-8,
                 max_amplitude: float = 1e-1,
                 simulation=None) -> None:
        """
        Initialize the adaptive square wave variant.
        
        Parameters
        ----------
        target_coupling : float
            Reference coupling strength (g_target)
        target_temperature : float
            Target temperature in Kelvin (T_target)
        period_ps : float
            Period of the square wave in picoseconds
        time_tracker : ElapsedTimeTracker
            Time tracker for accurate timing
        temperature_tracker : object
            Temperature tracker that provides harmonic_equipartition_temperature
        duty_cycle : float, optional
            Fraction of each period that the wave is "high". Default: 0.5
        phase_offset : float, optional
            Phase offset in radians. Default: 0.0
        start_time_ps : float, optional
            Time to start the square wave in picoseconds. Default: 0.0
        stop_time_ps : float, optional
            Time to stop the square wave in picoseconds. Default: None
        min_amplitude : float, optional
            Minimum allowed amplitude (safety limit). Default: 1e-8
        max_amplitude : float, optional
            Maximum allowed amplitude (safety limit). Default: 1e-1
        """
        hoomd.variant.Variant.__init__(self)
        
        # Core parameters
        self.target_coupling: float = float(target_coupling)
        self.target_temperature: float = float(target_temperature)
        self.current_amplitude: float = float(target_coupling)  # Start with target
        self.period_ps: float = float(period_ps)
        self.time_tracker: ElapsedTimeTracker = time_tracker
        self.temperature_tracker = temperature_tracker
        self.duty_cycle: float = float(duty_cycle)
        self.phase_offset: float = float(phase_offset)
        self.start_time_ps: float = float(start_time_ps)
        self.stop_time_ps: Optional[float] = (
            float(stop_time_ps) if stop_time_ps is not None else None
        )
        self.min_amplitude: float = float(min_amplitude)
        self.max_amplitude: float = float(max_amplitude)
        self.simulation = simulation
        
        # Validate parameters
        if not 0.0 <= self.duty_cycle <= 1.0:
            raise ValueError(f"Duty cycle must be between 0.0 and 1.0, got {self.duty_cycle}")
        if self.min_amplitude < 0.0:
            raise ValueError(f"Minimum amplitude must be >= 0.0, got {self.min_amplitude}")
        if self.max_amplitude <= self.min_amplitude:
            raise ValueError(f"Maximum amplitude must be > minimum amplitude")
        if self.target_temperature <= 0.0:
            raise ValueError(f"Target temperature must be > 0.0, got {self.target_temperature}")
        
        # Derived quantities
        self.angular_frequency: float = 2.0 * np.pi / self.period_ps  # rad/ps
        
        # State variables
        self.current_value: float = 0.0
        self.completed_periods: int = 0
        self.is_high: bool = False
        self._has_started: bool = False
        self._has_stopped: bool = False
        self._last_duty_cycle_end: float = -1.0
        self._amplitude_history: list = []
        
        print(f"AdaptiveSquareWaveVariant initialized:")
        print(f"  Target coupling (g_target): {self.target_coupling:.6e}")
        print(f"  Target temperature (T_target): {self.target_temperature:.1f} K")
        print(f"  Period: {self.period_ps} ps")
        print(f"  Duty cycle: {self.duty_cycle:.1%}")
        print(f"  Phase offset: {self.phase_offset:.3f} rad")
        print(f"  Start time: {self.start_time_ps} ps")
        print(f"  Amplitude limits: [{self.min_amplitude:.2e}, {self.max_amplitude:.2e}]")
        if self.stop_time_ps is not None:
            print(f"  Stop time: {self.stop_time_ps} ps")
        else:
            print(f"  No stop time (runs indefinitely)")
    
    def __call__(self, timestep: int) -> float:
        """
        Evaluate the adaptive square wave function at the given timestep.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep (unused, time comes from time_tracker)
            
        Returns
        -------
        float
            Current variant value (target_coupling or current_amplitude)
        """
        # Get current elapsed time from time tracker
        current_time_ps = self.time_tracker.elapsed_time
        
        # Return 0 before start time
        if current_time_ps < self.start_time_ps:
            self.current_value = 0.0
            self.is_high = False
            return 0.0
            
        # Check for auto-stop signal or time limit
        if self._check_auto_stop_signal() or (
            self.stop_time_ps is not None and current_time_ps >= self.stop_time_ps
        ):
            if not self._has_stopped:
                print(f"AdaptiveSquareWaveVariant: Auto-stopped at t = {current_time_ps:.6f} ps")
                print(f"  Transitioning to zero coupling: 0.0")
                self._has_stopped = True
            self.current_value = 0.0  # Return zero coupling for steady state
            self.is_high = False
            return 0.0
            
        # Main coupling logic after start time
        if not self._has_started:
            print(f"AdaptiveSquareWaveVariant: Started at t = {current_time_ps:.6f} ps")
            self._has_started = True
        
        # Calculate time since start and phase
        time_since_start = current_time_ps - self.start_time_ps
        phase = self.angular_frequency * time_since_start + self.phase_offset
        phase_mod = phase % (2 * np.pi)  # Normalize to [0, 2π]
        
        # Determine if we're in the "high" part of the duty cycle
        high_phase_duration = self.duty_cycle * 2 * np.pi
        
        # Check if we just completed a duty cycle (high period ended)
        was_high = self.is_high
        self.is_high = phase_mod < high_phase_duration
        
        # Update amplitude at the end of each duty cycle (transition from high to low)
        # or if this is the first cycle and we haven't updated yet
        should_update_amplitude = False
        
        if was_high and not self.is_high:
            # Just transitioned from high to low - end of duty cycle
            should_update_amplitude = True
        elif not hasattr(self, '_last_amplitude_update_time'):
            # First time through - initialize
            should_update_amplitude = True
            self._last_amplitude_update_time = current_time_ps
        
        if should_update_amplitude:
            self._update_amplitude_from_temperature()
            self._last_amplitude_update_time = current_time_ps
        
        # Set current value based on current state
        if self.is_high:
            self.current_value = self.current_amplitude
        else:
            self.current_value = 0.0  # Square wave goes to ZERO during low phase
        
        # Debug output every few cycles and state changes
        current_period = int(time_since_start / self.period_ps)
        if current_period != getattr(self, '_last_debug_period', -1):
            self._last_debug_period = current_period
            print(f"AdaptiveSquareWave: t={current_time_ps:.2f}ps, period={current_period}, "
                  f"high={self.is_high}, amplitude={self.current_amplitude:.2e}, value={self.current_value:.2e}, "
                  f"phase={phase:.3f}, phase_mod={phase_mod:.3f}, high_dur={high_phase_duration:.3f}")
        
        # Also debug state transitions
        if hasattr(self, '_last_is_high') and self._last_is_high != self.is_high:
            print(f"AdaptiveSquareWave STATE CHANGE: t={current_time_ps:.2f}ps, "
                  f"{'HIGH→LOW' if self._last_is_high else 'LOW→HIGH'}, value: {self.current_value:.2e}")
        self._last_is_high = self.is_high
        
        return self.current_value
    
    def _update_amplitude_from_temperature(self) -> None:
        """
        Update the amplitude based on bath temperature measurement.
        
        Implements: ε_new = ε_target * (T_initial / T_bath)^0.5
        where:
        - ε_target is the target coupling (from --coupling)
        - T_initial is the initial temperature (from --temperature)  
        - T_bath is the current bath temperature
        """
        try:
            # Get current bath temperature from simulation
            T_bath = self._get_current_bath_temperature()
            
            # Skip if temperature is not available or invalid
            if T_bath is None or T_bath <= 0.0:
                print(f"AdaptiveSquareWaveVariant: Warning - invalid bath temperature: {T_bath}")
                return
            
            # Calculate new amplitude: ε_new = ε_target * (T_initial / T_bath)^0.5
            ratio = self.target_temperature / T_bath  # T_initial / T_bath
            new_amplitude = self.target_coupling * np.sqrt(ratio)
            
            # Apply safety limits
            new_amplitude = np.clip(new_amplitude, self.min_amplitude, self.max_amplitude)
            
            # Store amplitude history for analysis
            self._amplitude_history.append({
                'time_ps': self.time_tracker.elapsed_time,
                'T_bath': T_bath,
                'T_initial': self.target_temperature,
                'ratio': ratio,
                'old_amplitude': self.current_amplitude,
                'new_amplitude': new_amplitude
            })
            
            # Update amplitude
            old_amplitude = self.current_amplitude
            self.current_amplitude = new_amplitude
            
            print(f"AdaptiveSquareWaveVariant: Amplitude updated at t = {self.time_tracker.elapsed_time:.3f} ps")
            print(f"  T_bath = {T_bath:.2f} K, T_initial = {self.target_temperature:.2f} K")
            print(f"  Ratio = {ratio:.4f}, sqrt(ratio) = {np.sqrt(ratio):.4f}")
            print(f"  Amplitude: {old_amplitude:.6e} → {new_amplitude:.6e}")
            
        except Exception as e:
            print(f"AdaptiveSquareWaveVariant: Error updating amplitude - {e}")
    
    def _get_current_bath_temperature(self) -> float:
        """
        Get the current bath temperature from the simulation.
        
        Returns
        -------
        float
            Current bath temperature in Kelvin, or None if not available
        """
        try:
            # Try to get temperature from simulation object
            if self.simulation is None:
                print("AdaptiveSquareWaveVariant: No simulation object available")
                return None
            
            # Get HOOMD simulation object
            hoomd_sim = getattr(self.simulation, 'sim', self.simulation)
            if hoomd_sim is None:
                print("AdaptiveSquareWaveVariant: No HOOMD simulation object available")
                return None
            
            # Try to get molecular thermostat temperature (primary choice)
            molecular_thermostat = getattr(self.simulation, 'molecular_thermostat_obj', None)
            if molecular_thermostat is not None:
                try:
                    # Handle different thermostat types
                    if hasattr(molecular_thermostat, 'kT'):
                        kT = molecular_thermostat.kT
                        # Handle hoomd.variant.Constant objects
                        if hasattr(kT, '__call__'):
                            kT_value = kT(0)  # Call with dummy timestep
                        else:
                            kT_value = float(kT)
                        # Convert from kT to temperature (kT is in energy units, divide by kB)
                        from .utils import PhysicalConstants
                        T_K = kT_value / PhysicalConstants.KB_HARTREE_PER_K
                        return T_K
                    elif hasattr(molecular_thermostat, 'thermostat') and hasattr(molecular_thermostat.thermostat, 'kT'):
                        # Handle nested thermostats (e.g., Bussi within ConstantVolume)
                        kT = molecular_thermostat.thermostat.kT
                        if hasattr(kT, '__call__'):
                            kT_value = kT(0)
                        else:
                            kT_value = float(kT)
                        from .utils import PhysicalConstants
                        T_K = kT_value / PhysicalConstants.KB_HARTREE_PER_K
                        return T_K
                except Exception as e:
                    print(f"AdaptiveSquareWaveVariant: Error getting molecular thermostat temperature: {e}")
            
            # Fallback: try cavity thermostat
            cavity_thermostat = getattr(self.simulation, 'cavity_thermostat_obj', None)
            if cavity_thermostat is not None:
                try:
                    if hasattr(cavity_thermostat, 'kT'):
                        kT = cavity_thermostat.kT
                        if hasattr(kT, '__call__'):
                            kT_value = kT(0)
                        else:
                            kT_value = float(kT)
                        from .utils import PhysicalConstants
                        T_K = kT_value / PhysicalConstants.KB_HARTREE_PER_K
                        return T_K
                except Exception as e:
                    print(f"AdaptiveSquareWaveVariant: Error getting cavity thermostat temperature: {e}")
            
            # Final fallback: use target temperature (no adaptation)
            print("AdaptiveSquareWaveVariant: Warning - could not get bath temperature, using target temperature")
            return self.target_temperature
            
        except Exception as e:
            print(f"AdaptiveSquareWaveVariant: Error in _get_current_bath_temperature: {e}")
            return self.target_temperature
    
    def _min(self) -> float:
        """Return the minimum value of the variant."""
        return 0.0
    
    def _max(self) -> float:
        """Return the maximum value of the variant."""
        return self.max_amplitude
    
    @property
    def has_started(self) -> bool:
        """Whether the square wave has started."""
        return self._has_started
    
    @property
    def has_stopped(self) -> bool:
        """Whether the square wave has been stopped."""
        return self._has_stopped
    
    @property
    def amplitude_history(self) -> list:
        """History of amplitude adjustments for analysis."""
        return self._amplitude_history.copy()
    
    @property
    def feedback_ratio(self) -> float:
        """Current feedback ratio (T_initial / T_bath)."""
        try:
            T_bath = self._get_current_bath_temperature()
            if T_bath is not None and T_bath > 0.0:
                return self.target_temperature / T_bath
            return 1.0
        except:
            return 1.0
    
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
                if getattr(self.simulation, 'auto_stop_coupling_signal', False):
                    return True
                    
            # Also check HOOMD simulation object if available
            hoomd_sim = getattr(self.simulation, 'sim', self.simulation)
            if hasattr(hoomd_sim, 'auto_stop_coupling_signal'):
                if getattr(hoomd_sim, 'auto_stop_coupling_signal', False):
                    return True
                    
        except Exception:
            pass  # Continue with normal operation if any check fails
            
        return False


class ExponentialWaveVariant(hoomd.variant.Variant):
    r"""
    Exponential wave variant for periodic coupling protocols with exponential decay profile.
    
    Implements a periodic exponential wave where each period starts at an amplitude
    and decays exponentially to zero with a specified time constant:
    
    .. math::
        
        g(t) = \begin{cases}
        0 & \text{if } t < t_{\text{start}} \\
        A \cdot e^{-\frac{t_{\text{period}}(t)}{\tau}} & \text{if } t_{\text{start}} \leq t < t_{\text{stop}} \text{ and within period} \\
        0 & \text{if } t \geq t_{\text{stop}}
        \end{cases}
    
    where :math:`t_{\text{period}}(t)` is the time within the current period:
    
    .. math::
        
        t_{\text{period}}(t) = (t - t_{\text{start}}) \bmod T
    
    and:
    - :math:`A` is the initial amplitude for each exponential pulse
    - :math:`\tau` is the exponential decay time constant in picoseconds  
    - :math:`T` is the period in picoseconds
    
    **Physical Interpretation:**
    
    Exponential wave coupling simulates experimental scenarios with:
    
    - Pulsed laser experiments with exponential decay envelopes
    - Cavity field decay due to losses during each pulse
    - Repeated exponential excitation protocols
    - Natural exponential relaxation processes
    - Pump-probe experiments with decaying pump intensity
    
    **Energy Conservation:**
    
    Energy input follows the exponential profile within each period, providing
    smooth and physically realistic coupling dynamics without abrupt discontinuities.
    
    Parameters
    ----------
    amplitude : float
        Initial amplitude :math:`A` at the start of each exponential pulse
    period_ps : float
        Period :math:`T` of the exponential wave in picoseconds
    tau_ps : float
        Exponential decay time constant :math:`\tau` in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing in adaptive timestep simulations
    start_time_ps : float, optional
        Time to start the exponential wave in picoseconds. Default: 0.0
    stop_time_ps : float, optional
        Time to stop the exponential wave in picoseconds. Default: None (never stop)
        
    Attributes
    ----------
    amplitude : float
        Exponential wave amplitude
    period_ps : float
        Period in picoseconds
    tau_ps : float
        Decay time constant in picoseconds
    start_time_ps : float
        Start time in picoseconds
    stop_time_ps : float or None
        Stop time in picoseconds, or None for no stop
    current_value : float
        Current value of the variant
    current_period : int
        Current period number (0-indexed)
    has_started : bool
        Whether the exponential wave has started
    has_stopped : bool
        Whether the exponential wave has been stopped
        
    Examples
    --------
    **Basic exponential wave:**
    
    >>> from hoomd.cavitymd.variants import ExponentialWaveVariant
    >>> from hoomd.cavitymd.analysis import ElapsedTimeTracker
    >>> 
    >>> # Create time tracker
    >>> time_tracker = ElapsedTimeTracker(sim, runtime_ps=10.0)
    >>> 
    >>> # Create exponential wave: amplitude 0.001, period 2 ps, tau 0.5 ps
    >>> exp_coupling = ExponentialWaveVariant(
    ...     amplitude=0.001,
    ...     period_ps=2.0,
    ...     tau_ps=0.5,
    ...     time_tracker=time_tracker
    ... )
    
    **Fast decay exponential wave:**
    
    >>> # Exponential wave with fast decay (tau = 0.1 * period)
    >>> fast_decay_coupling = ExponentialWaveVariant(
    ...     amplitude=0.002,
    ...     period_ps=5.0,
    ...     tau_ps=0.5,  # Fast decay
    ...     time_tracker=time_tracker
    ... )
    
    **Delayed finite-duration exponential wave:**
    
    >>> # Start at 1 ps, stop at 10 ps
    >>> finite_exp = ExponentialWaveVariant(
    ...     amplitude=0.001,
    ...     period_ps=1.5,
    ...     tau_ps=0.3,
    ...     time_tracker=time_tracker,
    ...     start_time_ps=1.0,
    ...     stop_time_ps=10.0
    ... )
    
    Notes
    -----
    - Requires ElapsedTimeTracker for accurate timing
    - Each period starts with full amplitude and decays exponentially
    - Provides smooth, continuous coupling profile (no discontinuities)
    - Compatible with both coupling strength and dissipation parameters
    - Period should be chosen relative to the decay time constant for desired behavior
    - When tau_ps >> period_ps, the decay is slow (nearly constant per period)
    - When tau_ps << period_ps, the decay is fast (sharp exponential pulses)
    """
    
    def __init__(self,
                 amplitude: float,
                 period_ps: float,
                 tau_ps: float,
                 time_tracker: ElapsedTimeTracker,
                 start_time_ps: float = 0.0,
                 stop_time_ps: Optional[float] = None,
                 adaptive: bool = False,
                 temperature_tracker = None,
                 simulation=None) -> None:
        """
        Initialize the exponential wave variant.
        
        Parameters
        ----------
        amplitude : float
            Exponential wave amplitude (target amplitude for adaptive mode)
        period_ps : float
            Period of the exponential wave in picoseconds
        tau_ps : float
            Exponential decay time constant in picoseconds
        time_tracker : ElapsedTimeTracker
            Time tracker for accurate timing
        start_time_ps : float, optional
            Time to start the exponential wave in picoseconds. Default: 0.0
        stop_time_ps : float, optional
            Time to stop the exponential wave in picoseconds. Default: None
        adaptive : bool, optional
            Enable adaptive amplitude scaling based on T_bath/T_harmonic. Default: False
        temperature_tracker : object, optional
            Temperature tracker for adaptive mode (required if adaptive=True)
        """
        hoomd.variant.Variant.__init__(self)
        self.target_amplitude: float = float(amplitude)  # Store as target for adaptive mode
        self.current_amplitude: float = float(amplitude)  # Current working amplitude
        self.period_ps: float = float(period_ps)
        self.tau_ps: float = float(tau_ps)
        self.time_tracker: ElapsedTimeTracker = time_tracker
        self.start_time_ps: float = float(start_time_ps)
        self.stop_time_ps: Optional[float] = (
            float(stop_time_ps) if stop_time_ps is not None else None
        )
        self.adaptive: bool = adaptive
        self.temperature_tracker = temperature_tracker
        self.simulation = simulation
        
        # Validate parameters
        if self.target_amplitude < 0:
            raise ValueError(f"Amplitude must be non-negative, got {self.target_amplitude}")
        if self.period_ps <= 0:
            raise ValueError(f"Period must be positive, got {self.period_ps}")
        if self.tau_ps <= 0:
            raise ValueError(f"Decay time constant must be positive, got {self.tau_ps}")
        
        # Validate adaptive mode requirements
        if self.adaptive and self.temperature_tracker is None:
            raise ValueError("Adaptive exponential wave requires temperature_tracker")
        
        # State variables
        self.current_value: float = 0.0
        self.current_period: int = 0
        self._has_started: bool = False
        self._has_stopped: bool = False
        self._amplitude_history: list = []  # For tracking adaptive behavior
        
        print(f"ExponentialWaveVariant initialized:")
        print(f"  Target amplitude: {self.target_amplitude}")
        print(f"  Period: {self.period_ps} ps")
        print(f"  Decay time constant: {self.tau_ps} ps")
        print(f"  Start time: {self.start_time_ps} ps")
        if self.stop_time_ps is not None:
            print(f"  Stop time: {self.stop_time_ps} ps")
        else:
            print(f"  No stop time (runs indefinitely)")
        if self.adaptive:
            print(f"  Adaptive mode: ENABLED (T_bath/T_harmonic scaling)")
        else:
            print(f"  Adaptive mode: DISABLED (fixed amplitude)")
    
    def __call__(self, timestep: int) -> float:
        """
        Evaluate the exponential wave function at the given timestep.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep (unused, time comes from time_tracker)
            
        Returns
        -------
        float
            Current variant value (exponential decay within period, with adaptive scaling)
        """
        # Check for auto-stop signal first
        if self._check_auto_stop_signal():
            self.current_value = 0.0
            return 0.0
        
        # Get current elapsed time from time tracker
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should stop
        if self.stop_time_ps is not None and current_time_ps >= self.stop_time_ps:
            if not self._has_stopped:
                print(f"ExponentialWaveVariant: Stopped at t = {current_time_ps:.6f} ps")
                self._has_stopped = True
            self.current_value = 0.0
            return 0.0
        
        # Check if we should start
        if current_time_ps >= self.start_time_ps:
            if not self._has_started:
                print(f"ExponentialWaveVariant: Started at t = {current_time_ps:.6f} ps")
                self._has_started = True
            
            # Calculate time since start
            time_since_start = current_time_ps - self.start_time_ps
            
            # Calculate which period we're in and time within current period
            current_period = int(time_since_start / self.period_ps)
            time_in_period = time_since_start % self.period_ps
            
            # Update amplitude for adaptive mode (at start of each new period)
            if self.adaptive and current_period != self.current_period:
                self._update_adaptive_amplitude()
                self.current_period = current_period
            
            # Calculate exponential decay within the current period
            # g(t) = A_current * exp(-t_period / tau)
            amplitude_to_use = self.current_amplitude if self.adaptive else self.target_amplitude
            self.current_value = amplitude_to_use * np.exp(-time_in_period / self.tau_ps)
            
            return self.current_value
        else:
            self.current_value = 0.0
            return 0.0
    
    def _min(self) -> float:
        """Return the minimum value of the variant."""
        # Minimum approaches 0 as t -> infinity within each period
        return 0.0
    
    def _max(self) -> float:
        """Return the maximum value of the variant."""
        return self.amplitude
    
    @property
    def has_started(self) -> bool:
        """Whether the exponential wave has started."""
        return self._has_started
    
    @property
    def has_stopped(self) -> bool:
        """Whether the exponential wave has been stopped."""
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
                #     return True
                    
        except Exception:
            pass  # Continue with normal operation if any check fails
            
        return False
    def _update_adaptive_amplitude(self) -> None:
        """
        Update the amplitude based on T_bath/T_harmonic scaling.
        
        Implements energy conservation: 0.5*coupling_target^2*T_bath = 0.5*coupling_new^2*T_harmonic
        Therefore: coupling_new = coupling_target * sqrt(T_bath / T_harmonic)
        """
        if not self.adaptive or self.temperature_tracker is None:
            return
            
        try:
            # Get current bath temperature (molecular thermostat temperature)
            T_bath = self._get_bath_temperature()
            
            # Get current harmonic equipartition temperature  
            T_harmonic = self._get_harmonic_temperature()
            
            # Skip if temperatures are not available or invalid
            if T_bath is None or T_harmonic is None or T_bath <= 0.0 or T_harmonic <= 0.0:
                print(f"ExponentialWaveVariant: Warning - invalid temperatures: T_bath={T_bath}, T_harmonic={T_harmonic}")
                return
            
            # Calculate new amplitude: coupling_new = coupling_target * sqrt(T_bath / T_harmonic)
            ratio = T_bath / T_harmonic
            new_amplitude = self.target_amplitude * np.sqrt(ratio)
            
            # Store amplitude history for analysis
            self._amplitude_history.append({
                "time_ps": self.time_tracker.elapsed_time,
                "period": self.current_period,
                "T_bath": T_bath,
                "T_harmonic": T_harmonic,
                "ratio": ratio,
                "old_amplitude": self.current_amplitude,
                "new_amplitude": new_amplitude
            })
            
            # Update amplitude
            old_amplitude = self.current_amplitude
            self.current_amplitude = new_amplitude
            
            print(f"ExponentialWaveVariant: Adaptive amplitude updated at t = {self.time_tracker.elapsed_time:.3f} ps (period {self.current_period})")
            print(f"  T_bath = {T_bath:.2f} K, T_harmonic = {T_harmonic:.2f} K")
            print(f"  Ratio = {ratio:.4f}, sqrt(ratio) = {np.sqrt(ratio):.4f}")
            print(f"  Amplitude: {old_amplitude:.6e} → {new_amplitude:.6e}")
            
        except Exception as e:
            print(f"ExponentialWaveVariant: Error updating adaptive amplitude - {e}")
    
    def _get_bath_temperature(self) -> Optional[float]:
        """Get current bath temperature from molecular thermostat."""
        try:
            if self.simulation is None:
                return None
                
            # Get molecular thermostat temperature (representing T_bath)
            molecular_thermostat = getattr(self.simulation, "molecular_thermostat_obj", None)
            if molecular_thermostat is not None:
                if hasattr(molecular_thermostat, "kT"):
                    kT = molecular_thermostat.kT
                    # Handle hoomd.variant.Constant objects
                    if hasattr(kT, "__call__"):
                        kT_value = kT(0)  # Call with dummy timestep
                    else:
                        kT_value = float(kT)
                    # Convert from kT to temperature (kT is in energy units, divide by kB)
                    from .utils import PhysicalConstants
                    T_K = kT_value / PhysicalConstants.KB_HARTREE_PER_K
                    return T_K
                elif hasattr(molecular_thermostat, "thermostat") and hasattr(molecular_thermostat.thermostat, "kT"):
                    # Handle nested thermostats (e.g., Bussi within ConstantVolume)
                    kT = molecular_thermostat.thermostat.kT
                    if hasattr(kT, "__call__"):
                        kT_value = kT(0)
                    else:
                        kT_value = float(kT)
                    from .utils import PhysicalConstants
                    T_K = kT_value / PhysicalConstants.KB_HARTREE_PER_K
                    return T_K
                    
        except Exception as e:
            print(f"ExponentialWaveVariant: Error getting bath temperature: {e}")
            
        return None
    
    def _get_harmonic_temperature(self) -> Optional[float]:
        """Get current harmonic equipartition temperature."""
        try:
            if self.temperature_tracker is None:
                return None
                
            # Get harmonic equipartition temperature
            if hasattr(self.temperature_tracker, "harmonic_equipartition_temperature"):
                return self.temperature_tracker.harmonic_equipartition_temperature
            elif hasattr(self.temperature_tracker, "get_harmonic_equipartition_temperature"):
                return self.temperature_tracker.get_harmonic_equipartition_temperature()
                
        except Exception as e:
            print(f"ExponentialWaveVariant: Error getting harmonic temperature: {e}")
            
        return None
    
    @property
    def amplitude_history(self) -> list:
        """Return history of adaptive amplitude adjustments."""
        return self._amplitude_history.copy()
    
    @property
    def feedback_ratio(self) -> float:
        """Current feedback ratio (T_bath / T_harmonic)."""
        try:
            T_bath = self._get_bath_temperature()
            T_harmonic = self._get_harmonic_temperature()
            if T_bath is not None and T_harmonic is not None and T_harmonic > 0.0:
                return T_bath / T_harmonic
            return 1.0
        except:
            return 1.0

