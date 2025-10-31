#!/usr/bin/env python3
"""
Classical PID Temperature Controller for Cavity MD Simulations

This module implements a standard PID controller for temperature regulation:
    u(t) = Kp * (e(t) + (1/Ti) * ∫e(t)dt + Td * de(t)/dt)

where u(t) = T_bath(t) is the bath temperature command.

Features:
- Multiple temperature signal choices (lj_coulombic, harmonic_fictive, kinetic)
- Two operating modes: setpoint tracking or self-loop
- Auto-tuning via FOPDT step response and Ziegler-Nichols rules
- Manual gain specification
- Anti-windup for integral term
- Filtered derivative for noise rejection
- Output constraints and rate limiting

Author: CavityMD Development Team
Date: 2025-10-30
"""

import numpy as np
import hoomd
from typing import Optional
from scipy.optimize import curve_fit

from ..utils import PhysicalConstants


class PIDControl(hoomd.custom.Action):
    """
    Classical PID temperature controller for cavity MD simulations.
    
    Implements the standard PID control law:
        u(t) = Kp * (e(t) + (1/Ti) * ∫e(t)dt + Td * de(t)/dt)
    
    where u(t) = T_bath(t) is the bath temperature.
    
    Parameters
    ----------
    temperature_tracker : TemperatureTracker
        Temperature tracker for signal measurements
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    simulation : CavityMDSimulation
        Reference to simulation object
    molecular_thermostat : hoomd thermostat, optional
        Molecular thermostat to control
    cavity_thermostat : hoomd thermostat, optional
        Cavity thermostat to control
    signal_choice : str
        Temperature signal to control: 'lj_coulombic', 'harmonic_fictive', or 'kinetic'
    target_temperature : float
        Setpoint temperature in K (ignored if self_loop=True)
    self_loop : bool
        If True, setpoint = T_bath(t) (useful for scheduled coupling)
    Kp : float, optional
        Proportional gain (manual tuning)
    Ti : float, optional
        Integral time constant in ps (manual tuning)
    Td : float, optional
        Derivative time constant in ps (manual tuning)
    auto_tune : bool
        Enable auto-tuning via step response (default: True if gains not specified)
    auto_tune_step_size : float
        Step change size for auto-tuning in K
    auto_tune_duration_ps : float
        Duration to observe step response in ps
    turn_on_time_ps : float
        Time to activate controller in ps
    turn_off_time_ps : float, optional
        Time to deactivate controller in ps
    update_interval_ps : float
        Controller update period in ps
    apply_to : str
        Apply control to 'molecular', 'cavity', or 'both' thermostats
    T_min : float
        Minimum bath temperature in K
    T_max : float, optional
        Maximum bath temperature in K
    rate_limit_K_per_ps : float, optional
        Maximum rate of temperature change in K/ps
    derivative_filter_N : float
        Derivative filter parameter (tau_f = Td/N)
    enable_anti_windup : bool
        Enable anti-windup for integral term
    output_file : str
        CSV output file path
    console_output_period_ps : float
        Console output period in ps
    empirical_data_file : str, optional
        Empirical data file for lj_coulombic or harmonic_fictive signals
        
    Examples
    --------
    **Auto-tuning mode with kinetic temperature:**
    
    >>> pid = PIDControl(
    ...     temperature_tracker=temp_tracker,
    ...     time_tracker=time_tracker,
    ...     simulation=sim,
    ...     molecular_thermostat=mol_thermo,
    ...     signal_choice='kinetic',
    ...     target_temperature=100.0,
    ...     auto_tune=True
    ... )
    
    **Manual tuning mode:**
    
    >>> pid = PIDControl(
    ...     temperature_tracker=temp_tracker,
    ...     time_tracker=time_tracker,
    ...     simulation=sim,
    ...     molecular_thermostat=mol_thermo,
    ...     signal_choice='lj_coulombic',
    ...     target_temperature=150.0,
    ...     Kp=0.5,
    ...     Ti=20.0,
    ...     Td=5.0
    ... )
    
    **Self-loop mode for scheduled coupling:**
    
    >>> pid = PIDControl(
    ...     temperature_tracker=temp_tracker,
    ...     time_tracker=time_tracker,
    ...     simulation=sim,
    ...     molecular_thermostat=mol_thermo,
    ...     signal_choice='lj_coulombic',
    ...     self_loop=True,
    ...     Kp=0.3,
    ...     Ti=15.0,
    ...     Td=3.0
    ... )
    """
    
    def __init__(self,
                 temperature_tracker,
                 time_tracker,
                 simulation,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 signal_choice: str = 'lj_coulombic',
                 target_temperature: float = 100.0,
                 self_loop: bool = False,
                 Kp: Optional[float] = None,
                 Ti: Optional[float] = None,
                 Td: Optional[float] = None,
                 auto_tune: bool = True,
                 auto_tune_step_size: float = 20.0,
                 auto_tune_duration_ps: float = 50.0,
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 apply_to: str = 'both',
                 T_min: float = 0.1,
                 T_max: Optional[float] = None,
                 rate_limit_K_per_ps: Optional[float] = None,
                 derivative_filter_N: float = 10.0,
                 enable_anti_windup: bool = True,
                 output_file: str = 'pid_control.csv',
                 console_output_period_ps: float = 1.0,
                 empirical_data_file: Optional[str] = None):
        
        super().__init__()
        
        # Store references
        self.temperature_tracker = temperature_tracker
        self.time_tracker = time_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        
        # Signal selection
        valid_signals = ['lj_coulombic', 'harmonic_fictive', 'kinetic']
        if signal_choice not in valid_signals:
            raise ValueError(f"signal_choice must be one of {valid_signals}, got '{signal_choice}'")
        self.signal_choice = signal_choice
        
        # Operating mode
        self.target_temperature = target_temperature
        self.self_loop = self_loop
        
        # Determine tuning mode
        manual_gains_provided = (Kp is not None) or (Ti is not None) or (Td is not None)
        self.manual_tuning = manual_gains_provided
        
        if self.manual_tuning:
            # Manual tuning mode
            self.Kp = Kp if Kp is not None else 1.0
            self.Ti = Ti if Ti is not None else 10.0
            self.Td = Td if Td is not None else 2.0
            self.auto_tune = False
            print(f"PID Controller: Manual tuning mode")
            print(f"  Kp = {self.Kp:.4f}")
            print(f"  Ti = {self.Ti:.4f} ps")
            print(f"  Td = {self.Td:.4f} ps")
        else:
            # Auto-tuning mode
            self.Kp = 1.0  # Initial guess
            self.Ti = 10.0
            self.Td = 2.0
            self.auto_tune = auto_tune
            self.auto_tune_step_size = auto_tune_step_size
            self.auto_tune_duration_ps = auto_tune_duration_ps
            self.auto_tune_complete = False
            self.auto_tune_start_time = None
            self.auto_tune_data_time = []
            self.auto_tune_data_signal = []
            self.auto_tune_initial_bath = None
            print(f"PID Controller: Auto-tuning mode")
            print(f"  Step size: {self.auto_tune_step_size} K")
            print(f"  Duration: {self.auto_tune_duration_ps} ps")
        
        # Timing
        self.turn_on_time_ps = turn_on_time_ps
        self.turn_off_time_ps = turn_off_time_ps
        self.update_interval_ps = update_interval_ps
        self.last_update_time = -1.0
        self.last_console_time = -1.0
        self.console_output_period_ps = console_output_period_ps
        
        # Control application
        self.apply_to = apply_to
        
        # Constraints
        self.T_min = T_min
        self.T_max = T_max
        self.rate_limit_K_per_ps = rate_limit_K_per_ps
        
        # Derivative filter
        self.derivative_filter_N = derivative_filter_N
        
        # Anti-windup
        self.enable_anti_windup = enable_anti_windup
        
        # Controller state
        self.integral_state = 0.0
        self.derivative_filtered = 0.0
        self.error_prev = 0.0
        self.current_bath_temperature = None
        self.bath_temperature_prev = None  # Previous bath temp for self-loop mode
        self.is_active = False
        self.steady_state_offset = 0.0  # Offset between bath and signal (T_bath - T_signal)
        
        # Output
        self.output_file = output_file
        self._initialize_output_file()
        
        print(f"PID Controller initialized:")
        print(f"  Signal: {self.signal_choice}")
        print(f"  Mode: {'Self-loop' if self.self_loop else f'Setpoint tracking (T* = {self.target_temperature} K)'}")
        print(f"  Apply to: {self.apply_to}")
        print(f"  Update interval: {self.update_interval_ps} ps")
        print(f"  Turn-on time: {self.turn_on_time_ps} ps")
        if self.turn_off_time_ps is not None:
            print(f"  Turn-off time: {self.turn_off_time_ps} ps")
        print(f"  Output constraints: T_min = {self.T_min} K" + (f", T_max = {self.T_max} K" if self.T_max else ""))
        if self.rate_limit_K_per_ps:
            print(f"  Rate limit: {self.rate_limit_K_per_ps} K/ps")
        print(f"  Anti-windup: {'enabled' if self.enable_anti_windup else 'disabled'}")
        print(f"  Output file: {self.output_file}")
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        with open(self.output_file, 'w') as f:
            f.write("# PID Temperature Controller Output\n")
            f.write("# time_ps: simulation time in picoseconds\n")
            f.write("# T_signal_K: measured signal temperature\n")
            f.write("# T_setpoint_K: setpoint temperature\n")
            f.write("# error_K: tracking error (setpoint - signal)\n")
            f.write("# P_term: proportional term\n")
            f.write("# I_term: integral term\n")
            f.write("# D_term: derivative term\n")
            f.write("# T_bath_K: commanded bath temperature\n")
            f.write("# active: controller active (1) or inactive (0)\n")
            f.write("# phase: current phase (pre_init, auto_tune, control)\n")
            f.write("time_ps,T_signal_K,T_setpoint_K,error_K,P_term,I_term,D_term,T_bath_K,active,phase\n")
    
    def _get_signal_temperature(self) -> Optional[float]:
        """Get current signal temperature based on signal_choice."""
        if self.signal_choice == 'lj_coulombic':
            return self.temperature_tracker.lj_coulombic_fictive_temperature
        elif self.signal_choice == 'harmonic_fictive':
            return self.temperature_tracker.harmonic_fictive_temperature
        elif self.signal_choice == 'kinetic':
            return self.temperature_tracker.kinetic_temperature
        return None
    
    def _get_current_bath_temperature(self) -> float:
        """Get current bath temperature from thermostats."""
        kT = None
        
        # Try molecular thermostat first
        if self.molecular_thermostat is not None:
            if hasattr(self.molecular_thermostat, 'kT'):
                kT_value = self.molecular_thermostat.kT
                if hasattr(kT_value, 'value'):
                    kT = kT_value.value
                elif hasattr(kT_value, '__call__'):
                    kT = kT_value(0)
                else:
                    kT = float(kT_value)
            elif hasattr(self.molecular_thermostat, 'thermostat'):
                nested = self.molecular_thermostat.thermostat
                if hasattr(nested, 'kT'):
                    kT_value = nested.kT
                    if hasattr(kT_value, 'value'):
                        kT = kT_value.value
                    elif hasattr(kT_value, '__call__'):
                        kT = kT_value(0)
                    else:
                        kT = float(kT_value)
        
        # Try cavity thermostat if molecular failed
        if kT is None and self.cavity_thermostat is not None:
            if hasattr(self.cavity_thermostat, 'kT'):
                kT_value = self.cavity_thermostat.kT
                if hasattr(kT_value, 'value'):
                    kT = kT_value.value
                elif hasattr(kT_value, '__call__'):
                    kT = kT_value(0)
                else:
                    kT = float(kT_value)
            elif hasattr(self.cavity_thermostat, 'thermostat'):
                nested = self.cavity_thermostat.thermostat
                if hasattr(nested, 'kT'):
                    kT_value = nested.kT
                    if hasattr(kT_value, 'value'):
                        kT = kT_value.value
                    elif hasattr(kT_value, '__call__'):
                        kT = kT_value(0)
                    else:
                        kT = float(kT_value)
        
        if kT is None:
            return 300.0  # Fallback
        
        return kT / PhysicalConstants.KB_HARTREE_PER_K
    
    def _apply_control(self, T_bath_cmd: float):
        """Apply bath temperature command to thermostats."""
        kT_au = T_bath_cmd * PhysicalConstants.KB_HARTREE_PER_K
        
        if self.apply_to in ['molecular', 'both']:
            if self.molecular_thermostat is not None:
                self.molecular_thermostat.kT = kT_au
        
        if self.apply_to in ['cavity', 'both']:
            if self.cavity_thermostat is not None:
                self.cavity_thermostat.kT = kT_au
    
    def _fopdt_response(self, t, K, theta, tau, T0):
        """First-order plus dead time (FOPDT) step response model."""
        y = np.zeros_like(t)
        for i, ti in enumerate(t):
            if ti <= theta:
                y[i] = T0
            else:
                y[i] = T0 + K * (1.0 - np.exp(-(ti - theta) / tau))
        return y
    
    def _fit_fopdt_model(self, time_data, signal_data, step_size):
        """
        Fit FOPDT model to step response data.
        
        Returns
        -------
        K : float
            Process gain (dimensionless)
        theta : float
            Dead time in ps
        tau : float
            Time constant in ps
        """
        time_array = np.array(time_data)
        signal_array = np.array(signal_data)
        
        # Normalize time to start at 0
        time_array = time_array - time_array[0]
        
        # Initial guess
        T0_guess = signal_array[0]
        K_guess = (signal_array[-1] - signal_array[0]) / step_size
        theta_guess = 1.0  # ps
        tau_guess = 10.0  # ps
        
        # Fit FOPDT model
        popt, pcov = curve_fit(
            self._fopdt_response,
            time_array,
            signal_array,
            p0=[K_guess, theta_guess, tau_guess, T0_guess],
            bounds=([0.0, 0.0, 0.1, -np.inf], [10.0, 50.0, 500.0, np.inf]),
            maxfev=10000
        )
        
        K, theta, tau, T0 = popt
        
        return K, theta, tau
    
    def _calculate_zn_gains(self, K, theta, tau):
        """
        Calculate PID gains using Ziegler-Nichols tuning rules for FOPDT model.
        
        Parameters
        ----------
        K : float
            Process gain
        theta : float
            Dead time in ps
        tau : float
            Time constant in ps
            
        Returns
        -------
        Kp : float
            Proportional gain
        Ti : float
            Integral time in ps
        Td : float
            Derivative time in ps
        """
        # Ziegler-Nichols PID tuning rules for FOPDT
        Kp = 1.2 * tau / (K * theta)
        Ti = 2.0 * theta
        Td = 0.5 * theta
        
        return Kp, Ti, Td
    
    def _perform_auto_tuning(self):
        """Perform auto-tuning using collected step response data."""
        print("\n" + "="*80)
        print("PID AUTO-TUNING")
        print("="*80)
        
        if len(self.auto_tune_data_time) < 10:
            print(f"ERROR: Insufficient data for auto-tuning ({len(self.auto_tune_data_time)} points)")
            print("Using default gains: Kp=1.0, Ti=10.0, Td=2.0")
            self.Kp = 1.0
            self.Ti = 10.0
            self.Td = 2.0
            return
        
        # Fit FOPDT model
        K, theta, tau = self._fit_fopdt_model(
            self.auto_tune_data_time,
            self.auto_tune_data_signal,
            self.auto_tune_step_size
        )
        
        print(f"FOPDT Model Identified:")
        print(f"  Process gain K = {K:.4f} (dimensionless)")
        print(f"  Dead time θ = {theta:.4f} ps")
        print(f"  Time constant τ = {tau:.4f} ps")
        
        # Calculate PID gains
        self.Kp, self.Ti, self.Td = self._calculate_zn_gains(K, theta, tau)
        
        print(f"\nZiegler-Nichols PID Gains:")
        print(f"  Kp = {self.Kp:.4f}")
        print(f"  Ti = {self.Ti:.4f} ps")
        print(f"  Td = {self.Td:.4f} ps")
        
        # Estimate steady-state offset between bath and signal
        # Use the final values from the step response
        # NOTE: This is only valid for setpoint-tracking mode, not self-loop
        final_signal = np.mean(self.auto_tune_data_signal[-10:])  # Average last 10 points
        final_bath = self.auto_tune_initial_bath + self.auto_tune_step_size
        if self.self_loop:
            self.steady_state_offset = 0.0  # Not applicable in self-loop mode
            print(f"\nSteady-State Analysis:")
            print(f"  Final signal temperature: {final_signal:.2f} K")
            print(f"  Final bath temperature: {final_bath:.2f} K")
            print(f"  Steady-state offset: Not used in self-loop mode")
        else:
            self.steady_state_offset = final_bath - final_signal
            print(f"\nSteady-State Analysis:")
            print(f"  Final signal temperature: {final_signal:.2f} K")
            print(f"  Final bath temperature: {final_bath:.2f} K")
            print(f"  Steady-state offset (T_bath - T_signal): {self.steady_state_offset:.2f} K")
        print("="*80 + "\n")
    
    def act(self, timestep):
        """Execute PID control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if it's time to update
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        self.last_update_time = current_time_ps
        dt = self.update_interval_ps
        
        # Get signal temperature
        T_signal = self._get_signal_temperature()
        if T_signal is None:
            return
        
        # Get current bath temperature
        if self.current_bath_temperature is None:
            self.current_bath_temperature = self._get_current_bath_temperature()
        
        # Phase management
        if current_time_ps < self.turn_on_time_ps:
            # Pre-initialization phase
            phase = 'pre_init'
            self.is_active = False
            self._log_control_data(current_time_ps, T_signal, 0.0, 0.0, 0.0, 0.0, 0.0,
                                   self.current_bath_temperature, 0, phase)
            return
        
        # Check for turn-off
        if self.turn_off_time_ps is not None and current_time_ps >= self.turn_off_time_ps:
            phase = 'turned_off'
            self.is_active = False
            self._log_control_data(current_time_ps, T_signal, 0.0, 0.0, 0.0, 0.0, 0.0,
                                   self.current_bath_temperature, 0, phase)
            return
        
        # Auto-tuning phase
        if self.auto_tune and not self.auto_tune_complete:
            phase = 'auto_tune'
            
            if self.auto_tune_start_time is None:
                # Start auto-tuning
                self.auto_tune_start_time = current_time_ps
                self.auto_tune_initial_bath = self.current_bath_temperature
                
                # Apply step change
                step_bath = self.auto_tune_initial_bath + self.auto_tune_step_size
                self._apply_control(step_bath)
                self.current_bath_temperature = step_bath
                
                print(f"\n[{current_time_ps:.2f} ps] Starting auto-tuning")
                print(f"  Initial bath: {self.auto_tune_initial_bath:.2f} K")
                print(f"  Step to: {step_bath:.2f} K")
                print(f"  Duration: {self.auto_tune_duration_ps} ps")
            
            # Collect data
            self.auto_tune_data_time.append(current_time_ps)
            self.auto_tune_data_signal.append(T_signal)
            
            # Check if auto-tuning is complete
            if current_time_ps - self.auto_tune_start_time >= self.auto_tune_duration_ps:
                self._perform_auto_tuning()
                self.auto_tune_complete = True
                
                # Reset to initial bath temperature
                self._apply_control(self.auto_tune_initial_bath)
                self.current_bath_temperature = self.auto_tune_initial_bath
                
                # Reset controller state
                self.integral_state = 0.0
                self.derivative_filtered = 0.0
                self.error_prev = 0.0
            
            self._log_control_data(current_time_ps, T_signal, 0.0, 0.0, 0.0, 0.0, 0.0,
                                   self.current_bath_temperature, 0, phase)
            return
        
        # Control phase
        phase = 'control'
        self.is_active = True
        
        # In self-loop mode, use PREVIOUS bath temperature as setpoint
        # This prevents runaway positive feedback
        if self.self_loop:
            if self.bath_temperature_prev is None:
                self.bath_temperature_prev = self.current_bath_temperature
            T_setpoint = self.bath_temperature_prev  # Use previous bath as reference
        else:
            T_setpoint = self.target_temperature
        
        # Calculate error
        error = T_setpoint - T_signal
        
        # Proportional term
        P = self.Kp * error
        
        # Integral term (trapezoidal integration)
        if self.Ti > 0:
            self.integral_state += dt * (error + self.error_prev) / 2.0
            I = self.Kp * (1.0 / self.Ti) * self.integral_state
        else:
            I = 0.0
        
        # Derivative term (filtered)
        if self.Td > 0:
            tau_f = self.Td / self.derivative_filter_N
            alpha = tau_f / (tau_f + dt)
            derivative_raw = self.Kp * self.Td * (error - self.error_prev) / dt
            self.derivative_filtered = alpha * self.derivative_filtered + (1.0 - alpha) * derivative_raw
            D = self.derivative_filtered
        else:
            D = 0.0
        
        # PID control law:
        # In self-loop mode: Start from previous bath, add PID correction
        # In setpoint mode: Use feedforward to handle steady-state offset
        if self.self_loop:
            # Self-loop: incremental control from previous bath temperature
            u_unsat = self.bath_temperature_prev + (P + I + D)
        else:
            # Setpoint tracking: feedforward + PID correction
            feedforward = T_setpoint + self.steady_state_offset
            pid_correction = P + I + D
            u_unsat = feedforward + pid_correction
        
        # Apply constraints
        if self.T_max is not None:
            u_sat = np.clip(u_unsat, self.T_min, self.T_max)
        else:
            u_sat = np.maximum(u_unsat, self.T_min)
        
        # Apply rate limit
        if self.rate_limit_K_per_ps is not None:
            max_change = self.rate_limit_K_per_ps * dt
            delta = u_sat - self.current_bath_temperature
            if abs(delta) > max_change:
                u_sat = self.current_bath_temperature + np.sign(delta) * max_change
        
        # Anti-windup: back-calculate integral correction
        if self.enable_anti_windup and self.Ti > 0:
            windup_error = u_sat - u_unsat
            self.integral_state += (dt / self.Ti) * windup_error
        
        # Apply control
        self._apply_control(u_sat)
        
        # Update state for next iteration
        if self.self_loop:
            self.bath_temperature_prev = u_sat  # Store for next step
        self.current_bath_temperature = u_sat
        
        # Update previous error
        self.error_prev = error
        
        # Console output
        if current_time_ps - self.last_console_time >= self.console_output_period_ps:
            self.last_console_time = current_time_ps
            if self.self_loop:
                print(f"[{current_time_ps:.2f} ps] PID: T_signal={T_signal:.2f} K, "
                      f"T_setpoint={T_setpoint:.2f} K (bath), error={error:+.2f} K, "
                      f"T_bath={u_sat:.2f} K (P={P:.2f}, I={I:.2f}, D={D:.2f})")
            else:
                print(f"[{current_time_ps:.2f} ps] PID: T_signal={T_signal:.2f} K, "
                      f"T_setpoint={T_setpoint:.2f} K, error={error:+.2f} K, "
                      f"T_bath={u_sat:.2f} K (FF={feedforward:.2f}, P={P:.2f}, I={I:.2f}, D={D:.2f})")
        
        # Log data
        self._log_control_data(current_time_ps, T_signal, T_setpoint, error, P, I, D, u_sat, 1, phase)
    
    def _log_control_data(self, time_ps, T_signal, T_setpoint, error, P, I, D, T_bath, active, phase):
        """Write control data to CSV file."""
        with open(self.output_file, 'a') as f:
            f.write(f"{time_ps:.6f},{T_signal:.6f},{T_setpoint:.6f},{error:.6f},"
                   f"{P:.6f},{I:.6f},{D:.6f},{T_bath:.6f},{active},{phase}\n")

