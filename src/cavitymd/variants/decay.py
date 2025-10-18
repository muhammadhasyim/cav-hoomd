# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Custom variant classes for cavity molecular dynamics."""

from typing import Union, Optional
import hoomd.variant
import numpy as np
from ..utils import PhysicalConstants
from ..analysis import ElapsedTimeTracker


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
    
    >>> from hoomd.cavitymd.coupling import ExponentialDecayVariant
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
    
    >>> from hoomd.cavitymd.coupling import DecayingSquareWaveVariant
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


