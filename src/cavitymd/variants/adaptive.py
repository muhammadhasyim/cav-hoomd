# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Custom variant classes for cavity molecular dynamics."""

from typing import Union, Optional
import hoomd.variant
import numpy as np
from ..utils import PhysicalConstants
from ..analysis import ElapsedTimeTracker


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
                        from ..utils import PhysicalConstants
                        T_K = kT_value / PhysicalConstants.KB_HARTREE_PER_K
                        return T_K
                    elif hasattr(molecular_thermostat, 'thermostat') and hasattr(molecular_thermostat.thermostat, 'kT'):
                        # Handle nested thermostats (e.g., Bussi within ConstantVolume)
                        kT = molecular_thermostat.thermostat.kT
                        if hasattr(kT, '__call__'):
                            kT_value = kT(0)
                        else:
                            kT_value = float(kT)
                        from ..utils import PhysicalConstants
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
                        from ..utils import PhysicalConstants
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
    
    >>> from hoomd.cavitymd.coupling import ExponentialWaveVariant
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
                    from ..utils import PhysicalConstants
                    T_K = kT_value / PhysicalConstants.KB_HARTREE_PER_K
                    return T_K
                elif hasattr(molecular_thermostat, "thermostat") and hasattr(molecular_thermostat.thermostat, "kT"):
                    # Handle nested thermostats (e.g., Bussi within ConstantVolume)
                    kT = molecular_thermostat.thermostat.kT
                    if hasattr(kT, "__call__"):
                        kT_value = kT(0)
                    else:
                        kT_value = float(kT)
                    from ..utils import PhysicalConstants
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

