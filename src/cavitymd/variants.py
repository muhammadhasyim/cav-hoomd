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
                 turn_off_time_ps: Optional[float] = None) -> None:
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
                if self.simulation.state.user_data.get('auto_stop_coupling_disabled', False):
                    return True
                    
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
                 stop_time_ps: Optional[float] = None) -> None:
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