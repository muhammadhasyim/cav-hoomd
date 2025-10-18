"""
Advanced control systems for cavity molecular dynamics simulations.

This module provides optimal and robust control algorithms for temperature regulation
in coupled thermal systems, including:
- LQR (Linear-Quadratic Regulator) with integral action
- Kalman filtering for state estimation
- In-situ system identification
"""

from typing import Optional, Dict, Any, Tuple
import hoomd
import numpy as np
import json
from datetime import datetime
from pathlib import Path

try:
    from scipy.linalg import solve_discrete_are
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("Warning: scipy not available. LQR controller will not be functional.")

from ..utils import PhysicalConstants


class DiffEqController(hoomd.custom.Action):
    """
    Differential equation temperature controller for cavity MD simulations.
    
    This controller implements a first-order differential equation to control
    the bath temperature based on a signal temperature:
    
    dT_bath(t)/dt = -(T_bath(t) - T_signal(t)) / τ
    
    where:
    - T_bath(t) is the bath temperature
    - T_signal(t) is the measured signal temperature (kinetic, LJ+Coulombic, etc.)
    - τ is the time constant controlling response speed
    
    The discretized update rule is:
    T_bath(t+dt) = T_bath(t) + dt * (-(T_bath(t) - T_signal(t)) / τ)
    T_bath(t+dt) = T_bath(t) - (dt/τ) * (T_bath(t) - T_signal(t))
    
    Parameters
    ----------
    temperature_method : str
        Temperature calculation method for signal ('kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_bath')
    time_constant_ps : float
        Time constant τ in picoseconds (controls response speed)
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    energy_tracker : EnergyTracker
        Energy tracker for temperature calculations
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Molecular thermostat to control
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Cavity thermostat to control
    apply_to : str, optional
        Which thermostats to control ('molecular', 'cavity', 'both'). Default: 'both'
    turn_on_time_ps : float, optional
        Time to start control in picoseconds. Default: 0.0
    turn_off_time_ps : float, optional
        Time to stop control in picoseconds. Default: None (never turn off)
    update_interval_ps : float, optional
        Interval between control updates in picoseconds. Default: 0.1
    T_min : float, optional
        Minimum allowed bath temperature in Kelvin. Default: 0.0
    T_max : float, optional
        Maximum allowed bath temperature in Kelvin. Default: None (no upper limit)
    rate_limit_K_per_ps : float, optional
        Maximum rate of temperature change in K/ps. Default: None (no rate limit)
    output_file : str, optional
        Output file for control data. Default: 'diffeq_feedback.csv'
    empirical_data_file : str, optional
        Path to empirical data file (required for certain methods)
    console_output_period_ps : float, optional
        Console output period in picoseconds. Default: 1.0
        
    Examples
    --------
    **Basic differential equation control following kinetic temperature:**
    
    >>> from hoomd.cavitymd.analysis import DiffEqController
    >>> 
    >>> diffeq_controller = DiffEqController(
    ...     temperature_method='kinetic',
    ...     time_constant_ps=5.0,  # 5 ps response time
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     apply_to='molecular'
    ... )
    
    **Fast response to LJ+Coulombic signal:**
    
    >>> diffeq_controller = DiffEqController(
    ...     temperature_method='lj_coulombic_bath',
    ...     time_constant_ps=1.0,  # Fast 1 ps response
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     empirical_data_file='lj_coulombic_calibration.txt',
    ...     apply_to='both'
    ... )
    
    Notes
    -----
    - Smaller time constants lead to faster response to signal changes
    - The bath temperature exponentially approaches the signal temperature
    - Unlike feedback controllers, this directly tracks the signal without a target
    - Compatible with all existing temperature calculation methods
    """
    
    def __init__(self,
                 temperature_method: str,
                 time_constant_ps: float,
                 time_tracker,
                 energy_tracker,
                 simulation=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 apply_to: str = 'both',
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 T_min: float = 0.0,
                 T_max: Optional[float] = None,
                 rate_limit_K_per_ps: Optional[float] = None,
                 output_file: str = 'diffeq_feedback.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0,
                 enable_bias_estimation: bool = True,
                 bias_process_noise: float = 1e-6,
                 bias_initial_covariance: float = 100.0):
        
        # Store core parameters
        self.temperature_method = temperature_method
        self.time_constant_ps = float(time_constant_ps)
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        
        # Control parameters
        self.apply_to = apply_to
        self.turn_on_time_ps = float(turn_on_time_ps)
        self.turn_off_time_ps = float(turn_off_time_ps) if turn_off_time_ps is not None else None
        self.update_interval_ps = float(update_interval_ps)
        self.T_min = float(T_min)
        self.T_max = float(T_max) if T_max is not None else None
        self.rate_limit_K_per_ps = float(rate_limit_K_per_ps) if rate_limit_K_per_ps is not None else None
        self.output_file = output_file
        self.console_output_period_ps = float(console_output_period_ps)
        
        # Kalman filter parameters for bias estimation
        self.enable_bias_estimation = bool(enable_bias_estimation)
        self.bias_process_noise = float(bias_process_noise)  # Q_bias: process noise for bias (should be small)
        self.bias_initial_covariance = float(bias_initial_covariance)  # P_bias initial uncertainty
        
        # Kalman filter state
        self.bias_estimate = 0.0  # Estimated measurement bias (K)
        self.bias_covariance = self.bias_initial_covariance  # Uncertainty in bias estimate
        
        # Load empirical data if needed
        self.empirical_data = None
        if empirical_data_file is not None and temperature_method in ['lj_coulombic', 'harmonic', 'lj_coulombic_bath']:
            try:
                energy_component = 'lj_coulombic' if 'lj_coulombic' in temperature_method else 'harmonic'
                self.empirical_data = EmpiricalTemperatureData(
                    data_file_path=empirical_data_file,
                    energy_component=energy_component
                )
            except Exception as e:
                print(f"Warning: Could not load empirical data file {empirical_data_file}: {e}")
        
        # Control state
        self.is_active = False
        self.current_bath_temperature = None
        self.last_update_time = 0.0
        self.last_console_output_time = 0.0
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("# Differential Equation Temperature Controller Output with Kalman Filter Bias Estimation\n")
                f.write("# time_ps: simulation time in picoseconds\n")
                f.write("# signal_temp_raw_K: raw measured signal temperature\n")
                f.write("# signal_temp_corrected_K: bias-corrected signal temperature\n")
                f.write("# bias_estimate_K: estimated measurement bias\n")
                f.write("# bias_covariance: uncertainty in bias estimate\n")
                f.write("# bath_temp_K: bath temperature setting\n")
                f.write("# dt_ps: timestep used for integration\n")
                f.write("# derivative_K_per_ps: dT_bath/dt calculated\n")
                f.write("# active: whether controller is active (1=active, 0=inactive)\n")
                f.write("time_ps,signal_temp_raw_K,signal_temp_corrected_K,bias_estimate_K,bias_covariance,"
                       f"bath_temp_K,dt_ps,derivative_K_per_ps,active\n")
        except Exception as e:
            print(f"Warning: Failed to initialize differential equation controller output file: {e}")
    
    def _calculate_system_temperature(self, current_time_ps: float) -> Optional[float]:
        """Calculate system temperature using the specified method."""
        # Use the same implementation as OffsetTemperatureController
        if self.temperature_method == 'kinetic':
            # Calculate kinetic temperature from particle velocities
            import numpy as np
            from ..utils import PhysicalConstants
            
            try:
                kinetic_energy = 0.0
                N_molecules = 0
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    for i, ptype in enumerate(snap.particles.typeid):
                        # Exclude cavity particle (assume it's the last particle)
                        if i < len(snap.particles.mass) - 1:
                            mass = snap.particles.mass[i]
                            velocity = np.array(snap.particles.velocity[i])
                            ke_particle = 0.5 * mass * np.sum(velocity**2)
                            kinetic_energy += ke_particle
                            N_molecules += 1
                
                if N_molecules > 0:
                    kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                    temperature_K = (2.0/3.0) * kinetic_energy / (N_molecules * kB_hartree)
                    return max(temperature_K, 0.0)
                else:
                    return None
                    
            except Exception as e:
                print(f"Warning: Could not calculate kinetic temperature: {e}")
                return None
        
        elif self.temperature_method == 'harmonic_equipartition':
            # Calculate harmonic temperature using equipartition theorem
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                harmonic_energy_hartree = energy_data.get('harmonic', 0.0)
                
                # Use equipartition theorem: E_harmonic = (3/2) * N * kB * T
                # Assume N=500 molecular particles (excluding cavity)
                N_molecules = 500
                from ..utils import PhysicalConstants
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                
                if harmonic_energy_hartree > 0:
                    temperature_K = (2.0/3.0) * harmonic_energy_hartree / (N_molecules * kB_hartree)
                    return max(temperature_K, 0.0)
                else:
                    return None
                    
            except Exception as e:
                print(f"Warning: Could not calculate harmonic equipartition temperature: {e}")
                return None
        
        elif self.temperature_method in ['lj_coulombic', 'lj_coulombic_bath']:
            # Calculate LJ+Coulombic fictive temperature using empirical data
            if self.empirical_data is None:
                print(f"Warning: No empirical data available for method '{self.temperature_method}'")
                return None
            
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                lj_energy = energy_data.get('lj', 0.0)
                coulomb_energy = energy_data.get('coulombic', 0.0)
                total_lj_coulomb = lj_energy + coulomb_energy
                
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    instantaneous_energy_hartree=total_lj_coulomb,
                    num_particles=500
                )
                return max(temperature_K, 0.0)
                
            except Exception as e:
                print(f"Warning: Could not calculate LJ+Coulombic temperature: {e}")
                return None
        
        elif self.temperature_method == 'harmonic':
            # Calculate harmonic fictive temperature using empirical data
            if self.empirical_data is None:
                print(f"Warning: No empirical data available for method '{self.temperature_method}'")
                return None
            
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                harmonic_energy = energy_data.get('harmonic', 0.0)
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    instantaneous_energy_hartree=harmonic_energy,
                    num_particles=500
                )
                return max(temperature_K, 0.0)
                
            except Exception as e:
                print(f"Warning: Could not calculate harmonic temperature: {e}")
                return None
        
        else:
            print(f"Warning: Unknown temperature method '{self.temperature_method}'")
            return None
    
    def _get_current_thermostat_temperature(self) -> float:
        """Get current thermostat temperature for initialization."""
        try:
            from ..utils import PhysicalConstants
            
            # Try molecular thermostat first
            if self.molecular_thermostat is not None:
                try:
                    if hasattr(self.molecular_thermostat, 'kT'):
                        kT = self.molecular_thermostat.kT
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.molecular_thermostat, 'thermostat'):
                        nested = self.molecular_thermostat.thermostat
                        if hasattr(nested, 'kT'):
                            kT = nested.kT
                            if hasattr(kT, 'value'):
                                kT = kT.value
                            elif callable(kT):
                                kT = kT()
                            return kT / PhysicalConstants.KB_HARTREE_PER_K
                except:
                    pass
            
            # Try cavity thermostat as backup
            if self.cavity_thermostat is not None:
                try:
                    if hasattr(self.cavity_thermostat, 'kT'):
                        kT = self.cavity_thermostat.kT
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.cavity_thermostat, 'thermostat'):
                        nested = self.cavity_thermostat.thermostat
                        if hasattr(nested, 'kT'):
                            kT = nested.kT
                            if hasattr(kT, 'value'):
                                kT = kT.value
                            elif callable(kT):
                                kT = kT()
                            return kT / PhysicalConstants.KB_HARTREE_PER_K
                except:
                    pass
            
            # Fallback to 300K
            return 300.0
            
        except Exception as e:
            print(f"Warning: Could not get current thermostat temperature: {e}")
            return 300.0
    
    def _update_bias_estimate(self, measured_signal: float, dt_ps: float):
        """
        Update bias estimate using Kalman filter.
        
        Model:
        ------
        State: b (bias in K)
        Dynamics: b[k+1] = b[k] + w  (random walk with tiny process noise)
        Measurement: y = T_signal_true + b
        
        At steady state (after controller converges):
            T_bath ≈ T_signal_true  (no bias case)
            Innovation: (measured_signal - T_bath) ≈ b
        
        Parameters
        ----------
        measured_signal : float
            Raw measured signal temperature (K)
        dt_ps : float
            Time step (ps)
        """
        if not self.enable_bias_estimation or self.current_bath_temperature is None:
            return
        
        # Kalman filter prediction step
        # State prediction: b_predict = b_estimate (constant bias assumption)
        b_predict = self.bias_estimate
        
        # Covariance prediction: P_predict = P + Q*dt
        P_predict = self.bias_covariance + self.bias_process_noise * dt_ps
        
        # Kalman filter update step
        # Innovation: difference between measured and expected signal
        # Expected signal = bath temperature (what the system should track without bias)
        # Measured signal = true signal + bias
        # Innovation = measured_signal - current_bath_temperature = (true_signal + bias) - bath
        # If controller is working well: true_signal ≈ bath, so innovation ≈ bias
        innovation = measured_signal - self.current_bath_temperature
        
        # Measurement noise estimate (assume ~1K RMS for temperature measurements)
        R = 1.0  # Measurement noise variance (K^2)
        
        # Kalman gain: K = P_predict / (P_predict + R)
        kalman_gain = P_predict / (P_predict + R)
        
        # State update: b_new = b_predict + K * innovation
        self.bias_estimate = b_predict + kalman_gain * innovation
        
        # Covariance update: P_new = (1 - K) * P_predict
        self.bias_covariance = (1.0 - kalman_gain) * P_predict

    def _apply_rate_limit(self, new_output: float, dt_ps: float) -> float:
        """Apply rate limiting to controller output."""
        if self.rate_limit_K_per_ps is None or self.current_bath_temperature is None:
            return new_output
        
        max_change = self.rate_limit_K_per_ps * dt_ps
        change = new_output - self.current_bath_temperature
        
        if abs(change) <= max_change:
            return new_output
        else:
            # Limit the change
            limited_change = max_change if change > 0 else -max_change
            return self.current_bath_temperature + limited_change
    
    def _update_thermostats(self, bath_temperature: float):
        """Update thermostat temperatures."""
        try:
            import hoomd
            from ..utils import PhysicalConstants
            kT_hartree = bath_temperature * PhysicalConstants.KB_HARTREE_PER_K
            
            if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
                if hasattr(self.molecular_thermostat, 'kT'):
                    self.molecular_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                elif hasattr(self.molecular_thermostat, 'thermostat'):
                    nested = self.molecular_thermostat.thermostat
                    if hasattr(nested, 'kT'):
                        nested.kT = hoomd.variant.Constant(kT_hartree)
            
            if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
                if hasattr(self.cavity_thermostat, 'kT'):
                    self.cavity_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                elif hasattr(self.cavity_thermostat, 'thermostat'):
                    nested = self.cavity_thermostat.thermostat
                    if hasattr(nested, 'kT'):
                        nested.kT = hoomd.variant.Constant(kT_hartree)
                        
        except Exception as e:
            print(f"Warning: Failed to update thermostats: {e}")
            import traceback
            traceback.print_exc()
    
    def act(self, timestep):
        """Execute differential equation control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should be active
        should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                          (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
        
        if should_be_active and not self.is_active:
            self.is_active = True
            # Initialize current_bath_temperature from thermostats when first activating
            if self.current_bath_temperature is None:
                self.current_bath_temperature = self._get_current_thermostat_temperature()
                print(f"Differential equation controller turned ON at t = {current_time_ps:.2f} ps")
                print(f"Initial bath temperature: {self.current_bath_temperature:.2f} K")
                print(f"Time constant: {self.time_constant_ps:.2f} ps")
            else:
                print(f"Differential equation controller turned ON at t = {current_time_ps:.2f} ps")
        elif not should_be_active and self.is_active:
            self.is_active = False
            print(f"Differential equation controller turned OFF at t = {current_time_ps:.2f} ps")
        
        # Skip if not active or not initialized
        if not self.is_active or self.current_bath_temperature is None:
            return
        
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        # Calculate dt_ps BEFORE updating last_update_time
        dt_ps = current_time_ps - self.last_update_time if self.last_update_time > 0 else self.update_interval_ps
        self.last_update_time = current_time_ps
        
        # Calculate raw signal temperature (may contain bias)
        signal_temperature_raw = self._calculate_system_temperature(current_time_ps)
        
        if signal_temperature_raw is None:
            return
        
        # Update bias estimate using Kalman filter
        self._update_bias_estimate(signal_temperature_raw, dt_ps)
        
        # Apply bias correction to get true signal temperature
        if self.enable_bias_estimation:
            signal_temperature_corrected = signal_temperature_raw - self.bias_estimate
        else:
            signal_temperature_corrected = signal_temperature_raw
        
        # Implement differential equation: dT_bath/dt = -(T_bath - T_signal_corrected) / τ
        # Discretized: T_bath(t+dt) = T_bath(t) + dt * (-(T_bath(t) - T_signal(t)) / τ)
        # Use CORRECTED signal to prevent runaway drift from bias
        derivative_K_per_ps = -(self.current_bath_temperature - signal_temperature_corrected) / self.time_constant_ps
        raw_output = self.current_bath_temperature + dt_ps * derivative_K_per_ps
        
        # Apply saturation limits
        if self.T_max is not None:
            saturated_output = max(self.T_min, min(self.T_max, raw_output))
        else:
            saturated_output = max(self.T_min, raw_output)
        
        # Apply rate limiting
        rate_limited_output = self._apply_rate_limit(saturated_output, dt_ps)
        
        # Update current bath temperature and thermostats
        self.current_bath_temperature = rate_limited_output
        self._update_thermostats(self.current_bath_temperature)
        
        # Console output
        if current_time_ps - self.last_console_output_time >= self.console_output_period_ps:
            if self.enable_bias_estimation:
                print(f"DiffEq Controller | Raw: {signal_temperature_raw:.1f}K | "
                      f"Bias: {self.bias_estimate:+.2f}K | Corrected: {signal_temperature_corrected:.1f}K | "
                      f"Bath: {self.current_bath_temperature:.1f}K | dT/dt: {derivative_K_per_ps:+.2f}K/ps")
            else:
                print(f"DiffEq Controller | {self.temperature_method}→{signal_temperature_raw:.1f}K | "
                      f"Bath: {self.current_bath_temperature:.1f}K | dT/dt: {derivative_K_per_ps:+.2f}K/ps")
            self.last_console_output_time = current_time_ps
        
        # Log data
        try:
            with open(self.output_file, 'a', encoding='utf-8') as f:
                f.write(f"{current_time_ps:.6f},{signal_temperature_raw:.6f},{signal_temperature_corrected:.6f},"
                       f"{self.bias_estimate:.6f},{self.bias_covariance:.6f},"
                       f"{self.current_bath_temperature:.6f},{dt_ps:.6f},{derivative_K_per_ps:.6f},{int(self.is_active)}\n")
        except Exception as e:
            print(f"Warning: Failed to write differential equation controller output: {e}")

