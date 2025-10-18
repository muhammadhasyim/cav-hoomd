# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Custom variant classes for cavity molecular dynamics."""

from typing import Union, Optional
import hoomd.variant
import numpy as np
from ..utils import PhysicalConstants
from ..analysis import ElapsedTimeTracker


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
    
    >>> from hoomd.cavitymd.coupling import PeriodicVariant
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
    
    >>> from hoomd.cavitymd.coupling import SquareWaveVariant
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


