
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Analysis and tracking tools for cavity molecular dynamics simulations."""

from typing import Optional, List, Dict, Union, Any, Tuple
import hoomd
import numpy as np
import logging
import sys
import os
import time
from pathlib import Path
import enum
import scipy.fft

# CuPy import with fallback for CPU/GPU agnostic code
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

from ..utils import PhysicalConstants, unwrap_positions
from .timing import ElapsedTimeTracker
from .empirical import EmpiricalTemperatureData

class EmpiricalTemperatureFeedback(hoomd.custom.Action):
    """
    Empirical temperature feedback controller for cavity MD simulations.
    
    This class measures instantaneous systemic temperature from potential energy
    using empirical calibration data, averages over time windows, and updates
    both molecular and cavity thermostat temperatures to maintain equilibrium.
    
    The feedback system:
    1. Measures instantaneous LJ+Coulomb energy at each timestep
    2. Converts to systemic temperature using Rosenfeld-fitted relationship
    3. Averages systemic temperature over configurable time windows
    4. Updates thermostat target temperatures at regular intervals
    5. Optionally turns off at specified time to allow free evolution
        
        Parameters
        ----------
    empirical_data : EmpiricalTemperatureData
        Calibrated energy-temperature relationship
    energy_tracker : EnergyTracker
        Energy tracker for reading instantaneous energies
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Molecular thermostat to control (None if not applicable)
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Cavity thermostat to control (None if not applicable)
    apply_to : str
        Which thermostats to update ('molecular', 'cavity', 'both', 'none'). 
        'none' means measurement only with no feedback control. Default: 'both'
    output_period_ps : float
        Period for CSV output in picoseconds. Default: 0.1
    averaging_window_ps : float
        Time window for temperature averaging in picoseconds. Default: 5.0
    update_interval_ps : float
        Interval between thermostat updates in picoseconds. Default: 5.0
    T_min : float
        Minimum allowed temperature in Kelvin. Default: 0.0
    T_max : float, optional
        Maximum allowed temperature in Kelvin. Default: None (no upper limit)
    turn_off_time_ps : float, optional
        Time to turn off feedback (None = never turn off). Default: None
    switch_time_ps : float, optional
        Time when coupling switches on (feedback starts after this). Default: None
    output_file : str
        Output file for feedback CSV data. Default: 'empirical_feedback.csv'
    """
    
    def __init__(self, 
                 empirical_data: EmpiricalTemperatureData,
                 energy_tracker,  # EnergyTracker type hint
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 apply_to: str = 'both',
                 output_period_ps: float = 0.1,
                 averaging_window_ps: float = 5.0,
                 update_interval_ps: float = 5.0,
                 T_min: float = 0.0,
                 T_max: Optional[float] = None,
                 turn_off_time_ps: Optional[float] = None,
                 switch_time_ps: Optional[float] = None,
                 output_file: str = 'empirical_feedback.csv',
                 time_tracker=None,
                 initial_temperature: float = 100.0):
        
        super().__init__()
        
        self.empirical_data = empirical_data
        self.energy_tracker = energy_tracker
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.apply_to = apply_to
        self.output_period_ps = output_period_ps
        self.averaging_window_ps = averaging_window_ps
        self.update_interval_ps = update_interval_ps
        self.T_min = T_min
        self.T_max = T_max
        self.turn_off_time_ps = turn_off_time_ps
        self.switch_time_ps = switch_time_ps
        self.output_file = output_file
        self.time_tracker = time_tracker
        self.initial_temperature = initial_temperature
        
        # Create separate harmonic empirical data for harmonic fictive temperature
        # This follows the approach from plot_temperature_feedback.py but uses available parameters
        try:
            self.harmonic_empirical_data = EmpiricalTemperatureData(
                data_file_path=empirical_data.data_file_path,
                energy_component='harmonic',
                use_direct_harmonic=False,  # Use fitted function (extended model)
                create_plots=True  # Create diagnostic plots for fitting validation
            )
            print(f"Created harmonic empirical data for fictive temperature calculation")
        except Exception as e:
            print(f"Warning: Could not create harmonic empirical data: {e}")
            print(f"   Will use direct equipartition formula for harmonic fictive temperature")
            self.harmonic_empirical_data = None
        
        # Internal state
        self.instantaneous_temperatures = []
        self.measurement_times = []
        self.last_output_time = 0.0
        self.last_update_time = 0.0
        self.feedback_active = True
        self.last_applied_temperature = None
        
        # Validation
        if self.apply_to not in ['molecular', 'cavity', 'both', 'none']:
            raise ValueError("apply_to must be 'molecular', 'cavity', 'both', or 'none'")
        
        if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is None:
            raise ValueError("molecular_thermostat cannot be None when apply_to includes 'molecular'")
        
        if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is None:
            raise ValueError("cavity_thermostat cannot be None when apply_to includes 'cavity'")
        
        # Initialize CSV output file
        self._initialize_output_file()
        
        print("Empirical temperature feedback controller initialized:")
        print(f"   Energy component: {self.empirical_data.energy_component}")
        print(f"   Apply to: {self.apply_to}")
        print(f"   Averaging window: {self.averaging_window_ps:.1f} ps")
        print(f"   Update interval: {self.update_interval_ps:.1f} ps")
        print(f"   Temperature range: [{self.T_min:.1f}, {self.T_max:.1f}] K")
        if self.switch_time_ps:
            print(f"   Switch time: {self.switch_time_ps:.1f} ps (feedback starts after)")
        if self.turn_off_time_ps:
            print(f"   Turn-off time: {self.turn_off_time_ps:.1f} ps")
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        with open(self.output_file, 'w') as f:
            f.write("# Empirical Temperature Feedback Data\n")
            f.write("# time_ps: Simulation time in picoseconds\n")
            f.write("# real_energy_hartree: Instantaneous potential energy in hartree\n")
            f.write("# harmonic_energy_hartree: Expected harmonic oscillator energy in hartree\n")
            f.write("# instantaneous_T_fictive_K: Instantaneous systemic temperature from empirical fit in Kelvin\n")
            f.write("# harmonic_T_fictive_K: Fictive temperature if energy were purely harmonic in Kelvin\n")
            f.write("# averaged_T_fictive_K: Time-averaged systemic temperature in Kelvin\n")
            f.write("# applied_T_K: Target temperature applied to thermostats in Kelvin\n")
            f.write("# feedback_active: Whether feedback is currently active (0/1)\n")
            f.write("time_ps,real_energy_hartree,harmonic_energy_hartree,instantaneous_T_fictive_K,harmonic_T_fictive_K,averaged_T_fictive_K,applied_T_K,feedback_active\n")

    def act(self, timestep):
        """Execute feedback control at each timestep."""
        from ..utils import PhysicalConstants
        
        # Get current simulation time using time tracker (FIX: was using wrong calculation)
        if hasattr(self, 'time_tracker') and self.time_tracker is not None:
            current_time_ps = self.time_tracker.elapsed_time
        else:
            # Fallback: this is still wrong but better than before
            # timestep is the step NUMBER, not the timestep size!
            dt_ps = PhysicalConstants.atomic_units_to_ps(1.0)  # Assume 1 au timestep
            current_time_ps = timestep * dt_ps
            print(f"Warning: Using fallback time calculation - time may be inaccurate")
        
        # Always track and output data, but only apply feedback control after switch time
        feedback_should_be_active = True
        if self.switch_time_ps is not None and current_time_ps < self.switch_time_ps:
            feedback_should_be_active = False  # Track data but don't apply feedback yet
        
        # Check if we should turn off feedback
        if self.turn_off_time_ps is not None and current_time_ps >= self.turn_off_time_ps:
            if self.feedback_active:
                self.feedback_active = False
                print(f"🔴 Empirical feedback turned OFF at t = {current_time_ps:.2f} ps")
        
        # Get instantaneous energy from energy tracker
        try:
            if hasattr(self.energy_tracker, 'get_instantaneous_energy'):
                energy_data = self.energy_tracker.get_instantaneous_energy()
                
                if self.empirical_data.energy_component == 'lj_coulombic':
                    if 'lj' in energy_data and 'coulombic' in energy_data:
                        instantaneous_energy = energy_data['lj'] + energy_data['coulombic']
                    else:
                        print("Warning: LJ+Coulombic energy not available")
                        return
                elif self.empirical_data.energy_component == 'total_PE':
                    if 'total_potential' in energy_data:
                        instantaneous_energy = energy_data['total_potential']
                    else:
                        print("Warning: Total potential energy not available")
                        return
                else:
                    print(f"Warning: Energy component {self.empirical_data.energy_component} not supported")
                    return
            else:
                # Fallback: try to access energy directly (this may not work)
                print("Warning: Energy tracker doesn't have get_instantaneous_energy method")
                return
                
        except Exception as e:
            print(f"Warning: Failed to get energy data: {e}")
            return
        
        # Calculate instantaneous systemic temperature
        instantaneous_T = self.empirical_data.calculate_systemic_temperature(instantaneous_energy)
        
        # Calculate expected harmonic oscillator energy using equipartition theorem
        # E_harmonic = 0.25 * N * k_B * T for a harmonic solid with N particles
        N_particles = 500  # Number of molecular particles (excluding cavity particle)
        harmonic_energy = 0.25 * N_particles * PhysicalConstants.KB_HARTREE_PER_K * instantaneous_T
        
        # Calculate harmonic fictive temperature using T^(3/2) fitted function approach
        # This follows the exact approach from plot_temperature_feedback.py calculate_empirical_harmonic_temperature
        try:
            # Try to get actual harmonic energy from the energy tracker
            if hasattr(self, 'energy_tracker') and self.energy_tracker is not None:
                harmonic_energy_dict = self.energy_tracker.get_instantaneous_energy()
                actual_harmonic_energy = harmonic_energy_dict.get('harmonic', None)
                
                if actual_harmonic_energy is not None and actual_harmonic_energy > 0:
                    # Use T^(3/2) fitted function for harmonic energy if we have empirical data
                    if hasattr(self, 'harmonic_empirical_data') and self.harmonic_empirical_data is not None:
                        harmonic_fictive_T = self.harmonic_empirical_data.calculate_systemic_temperature(
                            instantaneous_energy_hartree=actual_harmonic_energy,
                            num_particles=N_particles
                        )
                    else:
                        # Fallback: use direct equipartition formula from plot_temperature_feedback.py
                        kb_hartree_per_k = 3.16681e-6  # Hartree/K (from plot_temperature_feedback.py)
                        harmonic_fictive_T = 4.0 * actual_harmonic_energy / (N_particles * kb_hartree_per_k)
                else:
                    # Fallback: use equipartition temperature
                    harmonic_fictive_T = instantaneous_T
            else:
                # No energy tracker: use equipartition temperature
                harmonic_fictive_T = instantaneous_T
                
        except Exception as e:
            # Any error: use equipartition temperature
            harmonic_fictive_T = instantaneous_T
        
        # Store measurement
        self.instantaneous_temperatures.append(instantaneous_T)
        self.measurement_times.append(current_time_ps)
        
        # Remove old data outside averaging window
        cutoff_time = current_time_ps - self.averaging_window_ps
        while len(self.measurement_times) > 0 and self.measurement_times[0] < cutoff_time:
            self.measurement_times.pop(0)
            self.instantaneous_temperatures.pop(0)
        
        # Calculate averaged temperature
        if len(self.instantaneous_temperatures) > 0:
            averaged_T = np.mean(self.instantaneous_temperatures)
        else:
            averaged_T = instantaneous_T
        
        # Update thermostat temperatures if it's time
        # Default to initial temperature when no feedback has been applied
        applied_T = self.last_applied_temperature if self.last_applied_temperature is not None else self.initial_temperature
        
        if feedback_should_be_active and self.feedback_active and (current_time_ps - self.last_update_time) >= self.update_interval_ps:
            # Clamp averaged temperature to allowed range
            target_T = np.clip(averaged_T, self.T_min, self.T_max)
                
            # Update thermostats (this may do nothing if apply_to='none')
            self._update_thermostat_temperatures(target_T)
            
            self.last_update_time = current_time_ps
            
            # Only update applied_T if we're actually applying feedback
            if self.apply_to != 'none':
                self.last_applied_temperature = target_T
                applied_T = target_T
        
        # Output to CSV file
        if (current_time_ps - self.last_output_time) >= self.output_period_ps:
            with open(self.output_file, 'a') as f:
                # Show feedback as active only when it should actually be controlling
                feedback_status = int(feedback_should_be_active and self.feedback_active)
                f.write(f"{current_time_ps:.6f},{instantaneous_energy:.8f},{harmonic_energy:.8f},"
                       f"{instantaneous_T:.6f},{harmonic_fictive_T:.6f},{averaged_T:.6f},{applied_T:.6f},{feedback_status}\n")
            
            self.last_output_time = current_time_ps

    def _update_thermostat_temperatures(self, target_temperature_K: float):
        """Update thermostat target temperatures."""
        from ..utils import PhysicalConstants
        
        # Convert temperature to HOOMD energy units (kT in atomic units)
        target_kT = target_temperature_K * PhysicalConstants.KB_HARTREE_PER_K
        
        # Only update thermostats if apply_to is not 'none'
        if self.apply_to == 'none':
            return  # Measurement only, no feedback
            
        try:
            updated_any = False
            
            if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
                if hasattr(self.molecular_thermostat, 'kT'):
                    # Direct kT assignment (may not work if it's a Variant)
                    self.molecular_thermostat.kT = target_kT
                    updated_any = True
                elif hasattr(self.molecular_thermostat, '_kT'):
                    # Try private attribute
                    self.molecular_thermostat._kT = target_kT
                    updated_any = True
                else:
                    print("Warning: Cannot update molecular thermostat temperature")
            
            if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
                if hasattr(self.cavity_thermostat, 'kT'):
                    self.cavity_thermostat.kT = target_kT
                    updated_any = True
                elif hasattr(self.cavity_thermostat, '_kT'):
                    self.cavity_thermostat._kT = target_kT
                    updated_any = True
                else:
                    print("Warning: Cannot update cavity thermostat temperature")
            
            # Only print if we actually updated something
            if updated_any:
                print(f"Updated thermostat temperatures to {target_temperature_K:.2f} K")
            
        except Exception as e:
            print(f"Warning: Failed to update thermostat temperatures: {e}")


class GradientDescentTemperatureFeedback(hoomd.custom.Action):
    """
    Gradient descent temperature feedback controller for cavity MD simulations.
    
    This class implements a discretized gradient descent algorithm to minimize
    the objective function J = 1/2 * (T_eff(t) - T_target)^2, where:
    - T_eff(t) = (T_measured(t) + T_bath(t)) / 2 is the effective system temperature
    - T_target is the target temperature
    
    The gradient descent update rule is:
    T_bath(t+1) = T_bath(t) - α * ∇J/∇T_bath
    T_bath(t+1) = T_bath(t) - α * 0.5 * (T_eff(t) - T_target)
    
    where α is the learning rate determined by the time constant.
    
    Parameters
    ----------
    temperature_method : str
        Temperature calculation method ('kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition')
    time_constant_ps : float
        Time constant for gradient descent in picoseconds (controls convergence speed)
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    energy_tracker : EnergyTracker  
        Energy tracker for temperature calculations
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Molecular thermostat to control
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Cavity thermostat to control
    target_temperature : float, optional
        Target temperature in Kelvin. Default: 100.0
    apply_to : str, optional
        Which thermostats to control ('molecular', 'cavity', 'both'). Default: 'both'
    turn_on_time_ps : float, optional
        Time to start gradient descent control in picoseconds. Default: 0.0
    turn_off_time_ps : float, optional
        Time to stop gradient descent control in picoseconds. Default: None (never turn off)
    update_interval_ps : float, optional
        Interval between control updates in picoseconds. Default: 0.1 (every MD timestep)
    T_min : float, optional
        Minimum allowed bath temperature in Kelvin. Default: 0.0
    T_max : float, optional
        Maximum allowed bath temperature in Kelvin. Default: None (no upper limit)
    rate_limit_K_per_ps : float, optional
        Maximum rate of temperature change in K/ps. Default: None (no rate limit)
    output_file : str, optional
        Output file for gradient descent control data. Default: 'gd_feedback.csv'
    empirical_data_file : str, optional
        Path to empirical data file (required for 'lj_coulombic' and 'harmonic' methods)
        
    Attributes
    ----------
    alpha : float
        Learning rate (computed from time constant and MD timestep)
    current_bath_temperature : float
        Current bath temperature setting
    is_active : bool
        Whether gradient descent control is currently active
        
    Examples
    --------
    **Basic kinetic temperature gradient descent control:**
    
    >>> from hoomd.cavitymd.analysis import GradientDescentTemperatureFeedback
    >>> 
    >>> gd_controller = GradientDescentTemperatureFeedback(
    ...     temperature_method='kinetic',
    ...     time_constant_ps=10.0,  # 10 ps time constant
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     target_temperature=100.0
    ... )
    
    **Harmonic fictive temperature control with fast convergence:**
    
    >>> gd_controller = GradientDescentTemperatureFeedback(
    ...     temperature_method='harmonic',
    ...     time_constant_ps=5.0,  # Fast convergence
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     target_temperature=80.0,
    ...     empirical_data_file='harmonic_calibration.txt',
    ...     apply_to='both'
    ... )
    
    **Harmonic equipartition temperature control (no empirical data needed):**
    
    >>> gd_controller = GradientDescentTemperatureFeedback(
    ...     temperature_method='harmonic_equipartition',
    ...     time_constant_ps=3.0,  # Fast response
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     target_temperature=100.0,
    ...     apply_to='molecular'
    ... )
            
        Notes
        -----
    - Models system-bath coupling as T_eff = (T_measured + T_bath) / 2
    - Smaller time constants lead to faster but potentially less stable convergence
    - Learning rate α = dt / τ where dt is MD timestep and τ is time constant
    - Gradient descent is simpler and more intuitive than PI control
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
                 target_temperature: float = 100.0,
                 apply_to: str = 'both',
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 T_min: float = 0.0,
                 T_max: Optional[float] = None,
                 rate_limit_K_per_ps: Optional[float] = None,
                 output_file: str = 'gd_feedback.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0):
        
        super().__init__()
        
        # Store core parameters
        self.temperature_method = temperature_method
        self.time_constant_ps = float(time_constant_ps)
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        
        # Control parameters
        self.target_temperature = float(target_temperature)
        self.apply_to = apply_to
        self.turn_on_time_ps = float(turn_on_time_ps)
        self.turn_off_time_ps = float(turn_off_time_ps) if turn_off_time_ps is not None else None
        self.update_interval_ps = float(update_interval_ps)
        self.T_min = float(T_min)
        self.T_max = float(T_max) if T_max is not None else None
        self.rate_limit_K_per_ps = float(rate_limit_K_per_ps) if rate_limit_K_per_ps is not None else None
        self.output_file = output_file
        self.console_output_period_ps = float(console_output_period_ps)
        
        # Compute learning rate: α = dt / τ 
        # Assuming MD timestep is typically 0.1 fs = 0.0001 ps
        md_timestep_ps = 0.0001  # This could be made configurable
        self.alpha = md_timestep_ps / self.time_constant_ps
        
        # Initialize state
        self.current_bath_temperature = None  # Will be set from thermostats when active
        self.last_update_time = 0.0
        self.last_console_output_time = 0.0
        self.is_active = False
        
        # Set up empirical data if needed
        self.empirical_data = None
        if empirical_data_file and self.temperature_method in ['lj_coulombic', 'harmonic', 'lj_coulombic_kinetic', 'lj_coulombic_bath']:
            try:
                if self.temperature_method == 'harmonic':
                    self.empirical_data = EmpiricalTemperatureData(
                        empirical_data_file, energy_component='harmonic', create_plots=True
                    )
                else:  # lj_coulombic, lj_coulombic_kinetic, lj_coulombic_bath
                    self.empirical_data = EmpiricalTemperatureData(
                        empirical_data_file, energy_component='lj_coulombic', create_plots=True
                    )
                print(f"Loaded empirical data for {self.temperature_method} temperature calculation")
            except Exception as e:
                print(f"Warning: Failed to load empirical data: {e}")
        
        # Validate configuration
        self._validate_configuration()
        
        # Initialize output file
        self._initialize_output_file()
        
        # Print configuration
        print(f"Gradient Descent Temperature Controller initialized:")
        print(f"   Temperature method: {self.temperature_method}")
        print(f"   Time constant: {self.time_constant_ps:.1f} ps")
        print(f"   Learning rate α: {self.alpha:.6f}")
        print(f"   Target temperature: {self.target_temperature:.1f} K")
        print(f"   Apply to: {self.apply_to}")
        print(f"   Update interval: {self.update_interval_ps:.1f} ps")
        T_max_str = f"{self.T_max:.1f}" if self.T_max is not None else "∞"
        print(f"   Temperature limits: [{self.T_min:.1f}, {T_max_str}] K")
        if self.rate_limit_K_per_ps is not None:
            print(f"   Rate limit: {self.rate_limit_K_per_ps:.3f} K/ps")
        if self.turn_on_time_ps > 0:
            print(f"   Turn-on time: {self.turn_on_time_ps:.1f} ps")
        if self.turn_off_time_ps is not None:
            print(f"   Turn-off time: {self.turn_off_time_ps:.1f} ps")
    
    def _validate_configuration(self):
        """Validate controller configuration."""
        if self.temperature_method not in ['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic', 'lj_coulombic_bath']:
            raise ValueError(f"Invalid temperature_method: {self.temperature_method}")
        
        if self.apply_to not in ['molecular', 'cavity', 'both']:
            raise ValueError(f"apply_to must be 'molecular', 'cavity', or 'both', got {self.apply_to}")
        
        if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is None:
            raise ValueError("molecular_thermostat cannot be None when apply_to includes 'molecular'")
        
        if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is None:
            raise ValueError("cavity_thermostat cannot be None when apply_to includes 'cavity'")
        
        if self.time_constant_ps <= 0:
            raise ValueError("time_constant_ps must be positive")
        
        if self.target_temperature <= 0:
            raise ValueError("target_temperature must be positive")
        
        # Validate reasonable parameter ranges
        if self.time_constant_ps > 10000:
            print(f"Warning: Very large time constant ({self.time_constant_ps:.1f} ps) will result in extremely slow convergence")
            print(f"  Learning rate α = {self.alpha:.2e} - consider using a smaller time constant")
        
        if self.update_interval_ps < 0.01:
            print(f"Warning: Very small update interval ({self.update_interval_ps:.3f} ps) may cause excessive computational overhead")
        
        if self.alpha < 1e-8:
            print(f"Warning: Learning rate α = {self.alpha:.2e} is extremely small - controller may not respond")
            print(f"  Consider reducing time constant from {self.time_constant_ps:.1f} ps")
        
        if self.alpha > 0.1:
            print(f"Warning: Learning rate α = {self.alpha:.2e} is very large - controller may be unstable")
            print(f"  Consider increasing time constant from {self.time_constant_ps:.1f} ps")
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("# Gradient Descent Temperature Feedback Control Data\n")
                f.write("# time_ps: Simulation time in picoseconds\n")
                f.write("# measured_T_K: Measured system temperature in Kelvin\n")
                f.write("# effective_T_K: Effective temperature = (measured + bath) / 2 in Kelvin\n")
                f.write("# target_T_K: Target (setpoint) temperature in Kelvin\n")
                f.write("# error_K: Temperature error (effective - target) in Kelvin\n")
                f.write("# gradient: Computed gradient dJ/dT_bath = 0.5 * error\n")
                f.write("# raw_output_K: Raw gradient descent output in Kelvin\n")
                f.write("# saturated_output_K: Saturated output (within limits) in Kelvin\n")
                f.write("# rate_limited_output_K: Final output after rate limiting in Kelvin\n")
                f.write("# is_active: Whether gradient descent control is active (0/1)\n")
                f.write("time_ps,measured_T_K,effective_T_K,target_T_K,error_K,gradient,raw_output_K,saturated_output_K,rate_limited_output_K,is_active\n")
        except Exception as e:
            print(f"Warning: Failed to initialize gradient descent output file: {e}")
    
    def _get_hoomd_simulation(self):
        """Get reference to HOOMD simulation object."""
        return self.simulation
    
    def _calculate_system_temperature(self, current_time_ps: float) -> Optional[float]:
        """Calculate system temperature using the specified method."""
        try:
            if self.temperature_method == 'kinetic':
                # Calculate kinetic temperature directly from particle velocities
                # (same method as TemperatureTracker._calculate_kinetic_temperature)
                import numpy as np
                from ..utils import PhysicalConstants
                
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    velocities = np.array(snap.particles.velocity)
                    masses = np.array(snap.particles.mass)
                    N_particles = len(masses)
                
                if N_particles == 0:
                    return 0.0
                
                # Calculate kinetic energy
                kinetic_energy = 0.5 * np.sum(masses[:, np.newaxis] * velocities**2)
                
                # Temperature from equipartition: (3/2)*N*kB*T = KE
                # T = (2*KE)/(3*N*kB)
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (2.0 * kinetic_energy) / (3.0 * N_particles * kB_hartree)
                return temperature_K
                
            elif self.temperature_method == 'lj_coulombic':
                if self.empirical_data is None:
                    return None
                
                # Get LJ + Coulombic energy
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                lj_energy = energy_dict.get('lj', 0.0)
                coulomb_energy = energy_dict.get('coulombic', 0.0)
                total_energy = lj_energy + coulomb_energy
                
                # Convert to temperature using empirical data
                temperature = self.empirical_data.calculate_systemic_temperature(total_energy)
                return max(temperature, 0.0)  # Ensure non-negative temperature
                
            elif self.temperature_method == 'lj_coulombic_kinetic':
                # Dual temperature control: use maximum error signal from either LJ+Coulombic or kinetic
                # Error signal = sign * max(|T_lj_coulombic - T_target|, |T_kinetic - T_target|)
                
                # Calculate kinetic temperature
                import numpy as np
                from ..utils import PhysicalConstants
                
                try:
                    # Get kinetic energy from molecular particles
                    kinetic_energy = 0.0
                    N_molecules = 0
                    with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                        for i, ptype in enumerate(snap.particles.typeid):
                            # Exclude cavity particle (type index depends on system setup)
                            if i < len(snap.particles.mass) - 1:  # Assume cavity is last particle
                                mass = snap.particles.mass[i]
                                velocity = snap.particles.velocity[i]
                                # KE = (1/2) * m * v^2
                                ke_particle = 0.5 * mass * np.sum(velocity**2)
                                kinetic_energy += ke_particle
                                N_molecules += 1
                    
                    # Convert to temperature: T = (2/3) * KE / (N * kB)
                    if N_molecules > 0:
                        kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                        T_kinetic_K = (2.0/3.0) * kinetic_energy / (N_molecules * kB_hartree)
                    else:
                        return None
                        
                except Exception as e:
                    print(f"Warning: Failed to calculate kinetic temperature: {e}")
                    return None
                
                # Calculate LJ+Coulombic fictive temperature from energy data
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if self.empirical_data is None:
                    print(f"Warning: No empirical data available for LJ+Coulombic temperature")
                    return None
                
                if 'lj' not in energy_data or 'coulombic' not in energy_data:
                    print(f"Warning: LJ+Coulombic energy not available")
                    return None
                
                lj_coulombic_energy = energy_data['lj'] + energy_data['coulombic']
                T_lj_coulombic_K = self.empirical_data.calculate_systemic_temperature(lj_coulombic_energy)
                if T_lj_coulombic_K is None:
                    return None
                T_lj_coulombic_K = max(T_lj_coulombic_K, 0.0)
                
                # Calculate signed errors for directional control
                error_kinetic_signed = T_kinetic_K - self.target_temperature
                error_lj_coulombic_signed = T_lj_coulombic_K - self.target_temperature
                
                # Use directionally-aware selection: choose temperature requiring larger control action
                error_kinetic_abs = abs(error_kinetic_signed)
                error_lj_coulombic_abs = abs(error_lj_coulombic_signed)
                
                # Selection logic: choose temperature with larger absolute error
                # In case of tie (within 0.1 K), prioritize LJ+Coulombic (more fundamental)
                tolerance = 0.1  # K
                if abs(error_kinetic_abs - error_lj_coulombic_abs) <= tolerance:
                    selected_temperature = T_lj_coulombic_K
                    selected_method = "lj_coulombic"
                    reason = f"tie-break (Δerr={error_kinetic_abs-error_lj_coulombic_abs:.2f}K)"
                elif error_kinetic_abs > error_lj_coulombic_abs:
                    selected_temperature = T_kinetic_K
                    selected_method = "kinetic"
                    reason = f"larger error ({error_kinetic_abs:.1f}K > {error_lj_coulombic_abs:.1f}K)"
                else:
                    selected_temperature = T_lj_coulombic_K
                    selected_method = "lj_coulombic" 
                    reason = f"larger error ({error_lj_coulombic_abs:.1f}K > {error_kinetic_abs:.1f}K)"
                
                # Log the selection (occasionally)
                current_time_ps = self.time_tracker.elapsed_time
                if hasattr(self, '_last_dual_log_time'):
                    if current_time_ps - self._last_dual_log_time > 10.0:  # Log every 10 ps
                        print(f"Dual control | T_kinetic={T_kinetic_K:.1f}K (err={error_kinetic_signed:+.1f}K) | "
                              f"T_lj_coul={T_lj_coulombic_K:.1f}K (err={error_lj_coulombic_signed:+.1f}K) | "
                              f"Using: {selected_method} ({reason})")
                        self._last_dual_log_time = current_time_ps
                else:
                    self._last_dual_log_time = current_time_ps
                
                return selected_temperature
            
            elif self.temperature_method == 'lj_coulombic_bath':
                # Dual temperature control: use maximum error signal from either LJ+Coulombic or bath
                # Error signal = sign * max(|T_lj_coulombic - T_target|, |T_bath - T_target|)
                
                # Get bath temperature from current thermostat
                T_bath_K = self._get_current_thermostat_temperature()
                if T_bath_K is None:
                    print("Warning: Bath temperature not available")
                    return None
                
                # Calculate LJ+Coulombic fictive temperature from energy data (same as lj_coulombic_kinetic)
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if self.empirical_data is None:
                    print(f"Warning: No empirical data available for LJ+Coulombic temperature")
                    return None
                
                if 'lj' not in energy_data or 'coulombic' not in energy_data:
                    print(f"Warning: LJ+Coulombic energy not available")
                    return None
                
                lj_coulombic_energy = energy_data['lj'] + energy_data['coulombic']
                T_lj_coulombic_K = self.empirical_data.calculate_systemic_temperature(lj_coulombic_energy)
                if T_lj_coulombic_K is None:
                    return None
                T_lj_coulombic_K = max(T_lj_coulombic_K, 0.0)
                
                # Calculate signed errors for directional control
                error_bath_signed = T_bath_K - self.target_temperature
                error_lj_coulombic_signed = T_lj_coulombic_K - self.target_temperature
                
                # Use directionally-aware selection: choose temperature requiring larger control action
                error_bath_abs = abs(error_bath_signed)
                error_lj_coulombic_abs = abs(error_lj_coulombic_signed)
                
                # Selection logic: choose temperature with larger absolute error
                # In case of tie (within 0.1 K), prioritize LJ+Coulombic (more fundamental)
                tolerance = 0.1  # K
                if abs(error_bath_abs - error_lj_coulombic_abs) <= tolerance:
                    selected_temp = T_lj_coulombic_K
                    selected_type = "lj_coulombic"
                    reason = f"tie-break (Δerr={error_bath_abs-error_lj_coulombic_abs:.2f}K)"
                elif error_bath_abs > error_lj_coulombic_abs:
                    selected_temp = T_bath_K
                    selected_type = "bath"
                    reason = f"larger error ({error_bath_abs:.1f}K > {error_lj_coulombic_abs:.1f}K)"
                else:
                    selected_temp = T_lj_coulombic_K
                    selected_type = "lj_coulombic" 
                    reason = f"larger error ({error_lj_coulombic_abs:.1f}K > {error_bath_abs:.1f}K)"
                
                # Log the selection (occasionally)
                current_time_ps = self.time_tracker.elapsed_time
                if hasattr(self, '_last_bath_dual_log_time'):
                    if current_time_ps - self._last_bath_dual_log_time > 10.0:  # Log every 10 ps
                        print(f"Bath dual control | T_bath={T_bath_K:.1f}K (err={error_bath_signed:+.1f}K) | "
                              f"T_lj_coul={T_lj_coulombic_K:.1f}K (err={error_lj_coulombic_signed:+.1f}K) | "
                              f"Using: {selected_type} ({reason})")
                        self._last_bath_dual_log_time = current_time_ps
                else:
                    self._last_bath_dual_log_time = current_time_ps
                
                return selected_temp
                
            elif self.temperature_method == 'harmonic':
                if self.empirical_data is None:
                    return None
                
                # Get harmonic energy
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                
                # Convert to temperature using empirical data
                temperature = self.empirical_data.calculate_systemic_temperature(harmonic_energy)
                return max(temperature, 0.0)  # Ensure non-negative temperature
                
            elif self.temperature_method == 'harmonic_equipartition':
                # Calculate harmonic equipartition temperature directly from harmonic energy
                # Using the same method as TemperatureTracker._calculate_harmonic_equipartition_temperature
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                
                if harmonic_energy <= 0:
                    return 0.0
                
                # Get number of particles
                from ..utils import PhysicalConstants
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    N_particles = len(snap.particles.mass)
                
                if N_particles == 0:
                    return 0.0
                
                # Direct equipartition calculation: T = 4*E/(N*kB) for 3D harmonic oscillator
                # Based on empirical observations from TemperatureTracker implementation
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (4.0 * harmonic_energy) / (N_particles * kB_hartree)
                return temperature_K
            
            else:
                print(f"Warning: Unknown temperature method '{self.temperature_method}'")
                return None
                
        except Exception as e:
            print(f"Warning: Failed to calculate {self.temperature_method} temperature: {e}")
            return None
    
    def _get_current_thermostat_temperature(self) -> float:
        """Get current bath temperature from thermostats."""
        try:
            from ..utils import PhysicalConstants
            
            # Try to get temperature from molecular thermostat first
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
            
            # Fallback to target temperature
            return self.target_temperature
            
        except Exception as e:
            print(f"Warning: Could not get current thermostat temperature: {e}")
            return self.target_temperature

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
    
    def _update_thermostats(self, bath_temperature: float, measured_temperature: float = None, effective_temperature: float = None):
        """Update thermostat temperatures using robust API detection."""
        from ..utils import PhysicalConstants
        target_kT = bath_temperature * PhysicalConstants.KB_HARTREE_PER_K
        
        updated_any = False
        
        # Update molecular thermostat
        if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
            try:
                if hasattr(self.molecular_thermostat, 'kT'):
                    # Direct kT attribute (Langevin thermostats)
                    self.molecular_thermostat.kT = target_kT
                    updated_any = True
                elif hasattr(self.molecular_thermostat, 'thermostat'):
                    # Nested thermostat (ConstantVolume with Bussi/MTTK)
                    nested_thermo = self.molecular_thermostat.thermostat
                    if hasattr(nested_thermo, 'kT'):
                        nested_thermo.kT = target_kT
                        updated_any = True
                    else:
                        print(f"Warning: Cannot update nested molecular thermostat kT")
                else:
                    print(f"Warning: Cannot find kT attribute in molecular thermostat")
            except Exception as e:
                print(f"Warning: Failed to update molecular thermostat: {e}")
        
        # Update cavity thermostat  
        if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
            try:
                if hasattr(self.cavity_thermostat, 'kT'):
                    # Direct kT attribute (Langevin thermostats)
                    self.cavity_thermostat.kT = target_kT
                    updated_any = True
                elif hasattr(self.cavity_thermostat, 'thermostat'):
                    # Nested thermostat (ConstantVolume with Bussi/MTTK)
                    nested_thermo = self.cavity_thermostat.thermostat
                    if hasattr(nested_thermo, 'kT'):
                        nested_thermo.kT = target_kT
                        updated_any = True
                    else:
                        print(f"Warning: Cannot update nested cavity thermostat kT")
                else:
                    print(f"Warning: Cannot find kT attribute in cavity thermostat")
            except Exception as e:
                print(f"Warning: Failed to update cavity thermostat: {e}")
        
        # Check if we should print detailed console output (sync with console output period)
        current_time_ps = self.time_tracker.elapsed_time
        should_print_console = (current_time_ps - self.last_console_output_time) >= self.console_output_period_ps
        
        if updated_any and should_print_console:
            # Print detailed information
            if measured_temperature is not None and effective_temperature is not None:
                print(f"GD Controller | Target: {self.target_temperature:.1f}K | Measured: {measured_temperature:.1f}K | "
                      f"Effective: {effective_temperature:.1f}K | Bath→{bath_temperature:.1f}K")
            else:
                print(f"GD Controller | Target: {self.target_temperature:.1f}K | Bath→{bath_temperature:.1f}K")
            self.last_console_output_time = current_time_ps

    def act(self, timestep):
        """Execute gradient descent control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should be active
        should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                          (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
        
        if should_be_active and not self.is_active:
            self.is_active = True
            # Initialize current_bath_temperature from thermostats when first activating
            if self.current_bath_temperature is None:
                self.current_bath_temperature = self._get_current_thermostat_temperature()
                print(f"Gradient descent feedback turned ON at t = {current_time_ps:.2f} ps")
                print(f"Initial bath temperature: {self.current_bath_temperature:.2f} K")
            else:
                print(f"Gradient descent feedback turned ON at t = {current_time_ps:.2f} ps")
        elif not should_be_active and self.is_active:
            self.is_active = False
            print(f"Gradient descent feedback turned OFF at t = {current_time_ps:.2f} ps")
        
        # Skip if not active or not initialized
        if not self.is_active or self.current_bath_temperature is None:
            return
        
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        # Calculate dt_ps BEFORE updating last_update_time
        dt_ps = current_time_ps - self.last_update_time if self.last_update_time > 0 else self.update_interval_ps
        self.last_update_time = current_time_ps
        
        # Calculate system temperature
        measured_temperature = self._calculate_system_temperature(current_time_ps)
        
        if measured_temperature is None:
            return
        
        # Calculate effective system temperature using bath-system coupling model
        # Model: T_effective = (T_measured + T_bath) / 2
        effective_temperature = (measured_temperature + self.current_bath_temperature) / 2.0
        
        # Calculate error and gradient w.r.t. T_bath
        error = effective_temperature - self.target_temperature  # e = T_effective - T_target
        gradient = 0.5 * error  # ∂J/∂T_bath = 0.5 * (T_effective - T_target)
        
        # Apply gradient descent update: T_bath(t+1) = T_bath(t) - α * gradient
        raw_output = self.current_bath_temperature - self.alpha * gradient
        
        # Apply saturation limits
        if self.T_max is not None:
            saturated_output = max(self.T_min, min(self.T_max, raw_output))
        else:
            saturated_output = max(self.T_min, raw_output)
        
        # Apply rate limiting
        rate_limited_output = self._apply_rate_limit(saturated_output, dt_ps)
        
        # Update current bath temperature and thermostats
        self.current_bath_temperature = rate_limited_output
        self._update_thermostats(self.current_bath_temperature, measured_temperature, effective_temperature)
        
        # Log data
        try:
            with open(self.output_file, 'a', encoding='utf-8') as f:
                f.write(f"{current_time_ps:.6f},{measured_temperature:.6f},{effective_temperature:.6f},"
                       f"{self.target_temperature:.6f},{error:.6f},{gradient:.6f},{raw_output:.6f},"
                       f"{saturated_output:.6f},{rate_limited_output:.6f},{int(self.is_active)}\n")
        except Exception as e:
            print(f"Warning: Failed to write gradient descent output: {e}")


