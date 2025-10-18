
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

class DualIndependentTemperatureFeedback(hoomd.custom.Action):
    """
    Dual independent temperature feedback controller for cavity MD simulations.
    
    This class implements two independent gradient descent loops - one for the cavity
    bath temperature and one for the molecular bath temperature. Each bath responds
    to a different temperature signal, allowing for precise control of different
    aspects of the system dynamics.
    
    The gradient descent update rule for each bath is:
    T_bath(t+1) = T_bath(t) - α * ∇J/∇T_bath
    T_bath(t+1) = T_bath(t) - α * 0.5 * (T_eff(t) - T_target)
    
    where T_eff(t) = (T_measured(t) + T_bath(t)) / 2 for each signal.
    
    Parameters
    ----------
    cavity_method : str
        Temperature calculation method for cavity bath control
        ('kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic')
    molecular_method : str
        Temperature calculation method for molecular bath control
        ('kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic')
    cavity_time_constant_ps : float
        Time constant for cavity bath gradient descent in picoseconds
    molecular_time_constant_ps : float
        Time constant for molecular bath gradient descent in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    energy_tracker : EnergyTracker
        Energy tracker for temperature calculations
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat
        Molecular thermostat to control
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat
        Cavity thermostat to control
    cavity_target_temperature : float, optional
        Target temperature for cavity bath in Kelvin. Default: 100.0
    molecular_target_temperature : float, optional
        Target temperature for molecular bath in Kelvin. Default: 100.0
    turn_on_time_ps : float, optional
        Time to start control in picoseconds. Default: 0.0
    turn_off_time_ps : float, optional
        Time to stop control in picoseconds. Default: None (never turn off)
    update_interval_ps : float, optional
        Interval between control updates in picoseconds. Default: 0.1
    cavity_T_min : float, optional
        Minimum allowed cavity bath temperature in Kelvin. Default: 0.0
    cavity_T_max : float, optional
        Maximum allowed cavity bath temperature in Kelvin. Default: None
    molecular_T_min : float, optional
        Minimum allowed molecular bath temperature in Kelvin. Default: 0.0
    molecular_T_max : float, optional
        Maximum allowed molecular bath temperature in Kelvin. Default: None
    output_file : str, optional
        Output file for control data. Default: 'dual_feedback.csv'
    empirical_data_file : str, optional
        Path to empirical data file (required for certain methods)
    console_output_period_ps : float, optional
        Console output period in picoseconds. Default: 1.0
        
    Examples
    --------
    **Basic dual independent control:**
    
    >>> from hoomd.cavitymd.analysis import DualIndependentTemperatureFeedback
    >>> 
    >>> dual_controller = DualIndependentTemperatureFeedback(
    ...     cavity_method='harmonic_equipartition',
    ...     molecular_method='lj_coulombic_kinetic',
    ...     cavity_time_constant_ps=5.0,
    ...     molecular_time_constant_ps=10.0,
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     cavity_target_temperature=80.0,
    ...     molecular_target_temperature=120.0
    ... )
    
    Notes
    -----
    - Each bath has independent control parameters and targets
    - Different temperature signals can be used for different physical insights
    - Maintains compatibility with all existing temperature calculation methods
    """
    
    def __init__(self, 
                 cavity_method: str,
                 molecular_method: str,
                 cavity_time_constant_ps: float,
                 molecular_time_constant_ps: float,
                 time_tracker,
                 energy_tracker,
                 simulation=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 cavity_target_temperature: float = 100.0,
                 molecular_target_temperature: float = 100.0,
                 cavity_dynamic_target: bool = False,
                 molecular_dynamic_target: bool = False,
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 cavity_T_min: float = 0.0,
                 cavity_T_max: Optional[float] = None,
                 molecular_T_min: float = 0.0,
                 molecular_T_max: Optional[float] = None,
                 output_file: str = 'dual_feedback.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0,
                 cavity_integral_time_constant_ps: Optional[float] = None,
                 molecular_integral_time_constant_ps: Optional[float] = None,
                 cavity_heating_gain_factor: float = 1.0,
                 cavity_cooling_gain_factor: float = 1.0,
                 molecular_heating_gain_factor: float = 1.0,
                 molecular_cooling_gain_factor: float = 1.0,
                 cavity_integral_heating_time_constant_ps: Optional[float] = None,
                 cavity_integral_cooling_time_constant_ps: Optional[float] = None,
                 molecular_integral_heating_time_constant_ps: Optional[float] = None,
                 molecular_integral_cooling_time_constant_ps: Optional[float] = None):
        
        # Store configuration
        self.cavity_method = cavity_method
        self.molecular_method = molecular_method
        self.cavity_time_constant_ps = cavity_time_constant_ps
        self.molecular_time_constant_ps = molecular_time_constant_ps
        self.cavity_integral_time_constant_ps = cavity_integral_time_constant_ps
        self.molecular_integral_time_constant_ps = molecular_integral_time_constant_ps
        
        # Asymmetric gain factors
        self.cavity_heating_gain_factor = cavity_heating_gain_factor
        self.cavity_cooling_gain_factor = cavity_cooling_gain_factor
        self.molecular_heating_gain_factor = molecular_heating_gain_factor
        self.molecular_cooling_gain_factor = molecular_cooling_gain_factor
        
        # Asymmetric integral time constants
        self.cavity_integral_heating_time_constant_ps = cavity_integral_heating_time_constant_ps
        self.cavity_integral_cooling_time_constant_ps = cavity_integral_cooling_time_constant_ps
        self.molecular_integral_heating_time_constant_ps = molecular_integral_heating_time_constant_ps
        self.molecular_integral_cooling_time_constant_ps = molecular_integral_cooling_time_constant_ps
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.cavity_target_temperature = cavity_target_temperature
        self.molecular_target_temperature = molecular_target_temperature
        self.cavity_dynamic_target = cavity_dynamic_target
        self.molecular_dynamic_target = molecular_dynamic_target
        self.turn_on_time_ps = turn_on_time_ps
        self.turn_off_time_ps = turn_off_time_ps
        self.update_interval_ps = update_interval_ps
        self.cavity_T_min = cavity_T_min
        self.cavity_T_max = cavity_T_max
        self.molecular_T_min = molecular_T_min
        self.molecular_T_max = molecular_T_max
        self.output_file = output_file
        self.console_output_period_ps = console_output_period_ps
        
        # Initialize state variables
        self.last_update_time = 0.0
        self.last_console_output_time = 0.0
        self.is_active = False
        self.current_cavity_temperature = None
        self.current_molecular_temperature = None
        self.dynamic_targets_set = False
        
        # Calculate learning rates (alpha = dt / time_constant)
        # We'll update these with actual dt during first timestep
        self.cavity_alpha = None
        self.molecular_alpha = None
        
        # Integral control parameters
        self.cavity_alpha_i = None  # Integral gain for cavity
        self.molecular_alpha_i = None  # Integral gain for molecular
        
        # Error accumulators for integral terms
        self.cavity_error_integral = 0.0
        self.molecular_error_integral = 0.0
        
        # Set up empirical data if needed
        self.cavity_empirical_data = None
        self.molecular_empirical_data = None
        
        if empirical_data_file:
            # Load empirical data for cavity method if needed
            if self.cavity_method in ['lj_coulombic', 'harmonic', 'lj_coulombic_kinetic', 'lj_coulombic_bath']:
                try:
                    if self.cavity_method == 'harmonic':
                        self.cavity_empirical_data = EmpiricalTemperatureData(
                            empirical_data_file, energy_component='harmonic', create_plots=True
                        )
                    else:
                        self.cavity_empirical_data = EmpiricalTemperatureData(
                            empirical_data_file, energy_component='lj_coulombic', create_plots=True
                        )
                    print(f"Loaded cavity empirical data for {self.cavity_method} method")
                except Exception as e:
                    print(f"Warning: Failed to load cavity empirical data: {e}")
            
            # Load empirical data for molecular method if needed
            if self.molecular_method in ['lj_coulombic', 'harmonic', 'lj_coulombic_kinetic', 'lj_coulombic_bath']:
                try:
                    if self.molecular_method == 'harmonic':
                        self.molecular_empirical_data = EmpiricalTemperatureData(
                            empirical_data_file, energy_component='harmonic', create_plots=True
                        )
                    else:
                        self.molecular_empirical_data = EmpiricalTemperatureData(
                            empirical_data_file, energy_component='lj_coulombic', create_plots=True
                        )
                    print(f"Loaded molecular empirical data for {self.molecular_method} method")
                except Exception as e:
                    print(f"Warning: Failed to load molecular empirical data: {e}")
        
        # Validate configuration
        self._validate_configuration()
        
        # Initialize output file
        self._initialize_output_file()
        
        # Print configuration
        print(f"DualIndependent Controller initialized:")
        print(f"   Cavity: {self.cavity_method} → target {self.cavity_target_temperature:.1f}K (τ={self.cavity_time_constant_ps:.1f}ps)")
        print(f"   Molecular: {self.molecular_method} → target {self.molecular_target_temperature:.1f}K (τ={self.molecular_time_constant_ps:.1f}ps)")
        print(f"   Update interval: {self.update_interval_ps:.3f} ps")
        print(f"   Turn-on time: {self.turn_on_time_ps:.1f} ps")
        if self.turn_off_time_ps is not None:
            print(f"   Turn-off time: {self.turn_off_time_ps:.1f} ps")
    
    def _validate_configuration(self):
        """Validate controller configuration."""
        valid_methods = ['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic', 'lj_coulombic_bath']
        
        if self.cavity_method not in valid_methods:
            raise ValueError(f"Invalid cavity_method: {self.cavity_method}")
        
        if self.molecular_method not in valid_methods:
            raise ValueError(f"Invalid molecular_method: {self.molecular_method}")
        
        if self.molecular_thermostat is None:
            raise ValueError("molecular_thermostat cannot be None")
        
        if self.cavity_thermostat is None:
            raise ValueError("cavity_thermostat cannot be None")
        
        if self.energy_tracker is None:
            raise ValueError("energy_tracker is required for temperature calculations")
        
        if self.cavity_time_constant_ps <= 0:
            raise ValueError("cavity_time_constant_ps must be positive")
        
        if self.molecular_time_constant_ps <= 0:
            raise ValueError("molecular_time_constant_ps must be positive")
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("time_ps,cavity_temp,molecular_temp,cavity_effective,molecular_effective,")
                f.write("cavity_target,molecular_target,cavity_error,molecular_error,")
                f.write("cavity_gradient,molecular_gradient,cavity_integral,molecular_integral,")
                f.write("cavity_integral_term,molecular_integral_term,")
                f.write("cavity_p_gain_factor,molecular_p_gain_factor,")
                f.write("cavity_tau_i,molecular_tau_i,cavity_output,molecular_output,active\n")
        except Exception as e:
            print(f"Warning: Failed to initialize dual feedback output file: {e}")
    
    def _calculate_system_temperature(self, method: str, empirical_data, current_time_ps: float) -> Optional[float]:
        """Calculate system temperature using specified method."""
        try:
            if method == 'kinetic':
                # Calculate kinetic temperature from particle velocities
                import numpy as np
                from ..utils import PhysicalConstants
                
                # Get velocities from simulation state
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    velocities = np.array(snap.particles.velocity)
                    masses = np.array(snap.particles.mass)
                    types = np.array(snap.particles.typeid)
                
                # Filter out cavity particles (assuming they have type 1, molecular particles type 0)
                molecular_mask = (types == 0)
                molecular_velocities = velocities[molecular_mask]
                molecular_masses = masses[molecular_mask]
                
                if len(molecular_velocities) == 0:
                    return 0.0
                
                # Calculate kinetic energy per particle
                ke_per_particle = 0.5 * np.sum(molecular_masses[:, np.newaxis] * molecular_velocities**2, axis=1)
                mean_ke_per_particle = np.mean(ke_per_particle)
                
                # Convert to temperature: <KE>/particle = (3/2) * kB * T
                # Therefore: T = (2/3) * <KE>/particle / kB
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (2.0/3.0) * mean_ke_per_particle / kB_hartree
                
                return temperature_K
                
            elif method == 'lj_coulombic':
                if empirical_data is None:
                    return None
                
                # Get LJ+Coulombic energy
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                if 'lj' not in energy_dict or 'coulombic' not in energy_dict:
                    return None
                
                total_energy = energy_dict['lj'] + energy_dict['coulombic']
                
                # Convert to temperature using empirical data
                temperature = empirical_data.calculate_systemic_temperature(total_energy)
                return max(temperature, 0.0)
                
            elif method == 'lj_coulombic_kinetic':
                # Dual signal: use maximum error from either LJ+Coulombic or kinetic
                
                # Calculate kinetic temperature
                kinetic_temp = self._calculate_system_temperature('kinetic', None, current_time_ps)
                if kinetic_temp is None:
                    return None
                
                # Calculate LJ+Coulombic temperature
                lj_coul_temp = self._calculate_system_temperature('lj_coulombic', empirical_data, current_time_ps)
                if lj_coul_temp is None:
                    return kinetic_temp
                
                # Return the temperature that gives the maximum error magnitude
                # This requires knowing the target, so we'll use a simplified approach
                # and return the average of the two signals
                return (kinetic_temp + lj_coul_temp) / 2.0
                
            elif method == 'harmonic':
                if empirical_data is None:
                    return None
                
                # Get harmonic energy
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                
                # Convert to temperature using empirical data
                temperature = empirical_data.calculate_systemic_temperature(harmonic_energy)
                return max(temperature, 0.0)
                
            elif method == 'harmonic_equipartition':
                # Calculate harmonic equipartition temperature directly
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
                
                # Direct equipartition calculation
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (4.0 * harmonic_energy) / (N_particles * kB_hartree)
                return temperature_K
                
            elif method == 'lj_coulombic_bath':
                # Calculate temperature from molecular bath thermostat using LJ+Coulombic coupling
                # This method combines the bath temperature with LJ+Coulombic energy weighting
                
                # Get current molecular bath temperature
                molecular_bath_temp = self._get_thermostat_temperature(self.molecular_thermostat)
                if molecular_bath_temp is None:
                    return None
                
                # Get LJ+Coulombic energy for weighting
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                if 'lj' not in energy_dict or 'coulombic' not in energy_dict:
                    # Fallback to bath temperature if energy components unavailable
                    return molecular_bath_temp
                
                lj_coulombic_energy = energy_dict['lj'] + energy_dict['coulombic']
                
                # Use empirical data if available for energy-based correction
                if empirical_data is not None:
                    try:
                        energy_temp = empirical_data.calculate_systemic_temperature(lj_coulombic_energy)
                        # Weighted average: 70% bath temperature, 30% energy-derived temperature
                        # This provides a stable signal that incorporates both thermostat control
                        # and instantaneous energy fluctuations
                        return 0.7 * molecular_bath_temp + 0.3 * max(energy_temp, 0.0)
                    except Exception:
                        # Fallback to bath temperature if empirical calculation fails
                        return molecular_bath_temp
                else:
                    # Without empirical data, return the bath temperature directly
                    # This ensures the method always provides a meaningful signal
                    return molecular_bath_temp
            
            else:
                print(f"Warning: Unknown temperature method '{method}'")
                return None
                
        except Exception as e:
            print(f"Warning: Could not calculate temperature using method '{method}': {e}")
            return None
    
    def _get_thermostat_temperature(self, thermostat) -> Optional[float]:
        """Get current temperature from a thermostat."""
        try:
            from ..utils import PhysicalConstants
            
            if hasattr(thermostat, 'kT'):
                # Direct kT attribute (Langevin thermostats)
                kT_value = thermostat.kT
                # Handle HOOMD variants (e.g., Constant) 
                if hasattr(kT_value, '__call__'):
                    # It's a variant, call it with timestep 0 to get current value
                    kT_hartree = kT_value(0)
                else:
                    # It's a plain number
                    kT_hartree = kT_value
                return kT_hartree / PhysicalConstants.KB_HARTREE_PER_K
            elif hasattr(thermostat, 'thermostat'):
                # Nested thermostat (ConstantVolume with Bussi/MTTK)
                nested_thermo = thermostat.thermostat
                if hasattr(nested_thermo, 'kT'):
                    kT_value = nested_thermo.kT
                    # Handle HOOMD variants
                    if hasattr(kT_value, '__call__'):
                        kT_hartree = kT_value(0)
                    else:
                        kT_hartree = kT_value
                    return kT_hartree / PhysicalConstants.KB_HARTREE_PER_K
            
            return None
            
        except Exception as e:
            print(f"Warning: Failed to get thermostat temperature: {e}")
            return None
    
    def _update_thermostat(self, thermostat, temperature_K: float) -> bool:
        """Update a thermostat temperature."""
        try:
            import hoomd
            from ..utils import PhysicalConstants
            target_kT = temperature_K * PhysicalConstants.KB_HARTREE_PER_K
            
            if hasattr(thermostat, 'kT'):
                # Direct kT attribute (Langevin thermostats)
                # Always use Constant variant for proper HOOMD compatibility
                thermostat.kT = hoomd.variant.Constant(target_kT)
                return True
            elif hasattr(thermostat, 'thermostat'):
                # Nested thermostat (ConstantVolume with Bussi/MTTK)
                nested_thermo = thermostat.thermostat
                if hasattr(nested_thermo, 'kT'):
                    nested_thermo.kT = hoomd.variant.Constant(target_kT)
                    return True
            
            return False
            
        except Exception as e:
            print(f"Warning: Failed to update thermostat: {e}")
            return False
    
    def _set_dynamic_targets(self, current_time_ps: float):
        """Set dynamic targets based on current bath temperatures when controller is activated."""
        try:
            # Get current bath temperatures from thermostats (not signal temperatures)
            if self.cavity_dynamic_target:
                cavity_bath_temp = self._get_thermostat_temperature(self.cavity_thermostat)
                if cavity_bath_temp is not None:
                    old_target = self.cavity_target_temperature
                    self.cavity_target_temperature = cavity_bath_temp
                    print(f"Dynamic cavity target set: {old_target:.1f} K → {cavity_bath_temp:.1f} K (bath temperature)")
                else:
                    print(f"Warning: Could not get cavity bath temperature for dynamic target, keeping {self.cavity_target_temperature:.1f} K")
            
            if self.molecular_dynamic_target:
                molecular_bath_temp = self._get_thermostat_temperature(self.molecular_thermostat)
                if molecular_bath_temp is not None:
                    old_target = self.molecular_target_temperature
                    self.molecular_target_temperature = molecular_bath_temp
                    print(f"Dynamic molecular target set: {old_target:.1f} K → {molecular_bath_temp:.1f} K (bath temperature)")
                else:
                    print(f"Warning: Could not get molecular bath temperature for dynamic target, keeping {self.molecular_target_temperature:.1f} K")
            
            self.dynamic_targets_set = True
            
        except Exception as e:
            print(f"Warning: Error setting dynamic targets: {e}")
            self.dynamic_targets_set = True  # Don't try again
        
    def act(self, timestep):
        """Execute dual independent control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should be active
        should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                          (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
        
        if should_be_active and not self.is_active:
            self.is_active = True
            # Initialize current temperatures from thermostats when first activating
            self.current_cavity_temperature = self._get_thermostat_temperature(self.cavity_thermostat)
            self.current_molecular_temperature = self._get_thermostat_temperature(self.molecular_thermostat)
            
            # Initialize learning rates based on actual timestep
            if self.cavity_alpha is None:
                # Estimate timestep from update interval (rough approximation)
                estimated_dt_ps = self.update_interval_ps
                
                # Proportional gains (existing)
                self.cavity_alpha = estimated_dt_ps / self.cavity_time_constant_ps
                self.molecular_alpha = estimated_dt_ps / self.molecular_time_constant_ps
                
                # Integral gains (new PI controller)
                if self.cavity_integral_time_constant_ps is not None:
                    self.cavity_alpha_i = estimated_dt_ps / self.cavity_integral_time_constant_ps
                    print(f"Cavity PI controller: τₚ={self.cavity_time_constant_ps:.1f}ps, τᵢ={self.cavity_integral_time_constant_ps:.1f}ps")
        else:
                self.cavity_alpha_i = 0.0
                print(f"Cavity P controller: τₚ={self.cavity_time_constant_ps:.1f}ps (no integral)")
                
                if self.molecular_integral_time_constant_ps is not None:
                    self.molecular_alpha_i = estimated_dt_ps / self.molecular_integral_time_constant_ps
                    print(f"Molecular PI controller: τₚ={self.molecular_time_constant_ps:.1f}ps, τᵢ={self.molecular_integral_time_constant_ps:.1f}ps")
                else:
                    self.molecular_alpha_i = 0.0
                    print(f"Molecular P controller: τₚ={self.molecular_time_constant_ps:.1f}ps (no integral)")
            
        print(f"Dual independent feedback turned ON at t = {current_time_ps:.2f} ps")
        cavity_temp_str = f"{self.current_cavity_temperature:.2f}" if self.current_cavity_temperature is not None else "None"
        molecular_temp_str = f"{self.current_molecular_temperature:.2f}" if self.current_molecular_temperature is not None else "None"
        print(f"Initial cavity temperature: {cavity_temp_str} K")
        print(f"Initial molecular temperature: {molecular_temp_str} K")
            
        # Set dynamic targets if enabled
        if not self.dynamic_targets_set:
            self._set_dynamic_targets(current_time_ps)
            
        elif not should_be_active and self.is_active:
            self.is_active = False
            print(f"Dual independent feedback turned OFF at t = {current_time_ps:.2f} ps")
        
        # Skip if not active or not initialized
        if not self.is_active or self.current_cavity_temperature is None or self.current_molecular_temperature is None:
            return
        
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        # Calculate dt_ps BEFORE updating last_update_time
        dt_ps = current_time_ps - self.last_update_time if self.last_update_time > 0 else self.update_interval_ps
        self.last_update_time = current_time_ps
        
        # Calculate system temperatures for each method
        cavity_measured_temp = self._calculate_system_temperature(
            self.cavity_method, self.cavity_empirical_data, current_time_ps)
        molecular_measured_temp = self._calculate_system_temperature(
            self.molecular_method, self.molecular_empirical_data, current_time_ps)
        
        if cavity_measured_temp is None or molecular_measured_temp is None:
            return
        
        # Calculate effective temperatures and errors for each bath
        cavity_effective_temp = (cavity_measured_temp + self.current_cavity_temperature) / 2.0
        molecular_effective_temp = (molecular_measured_temp + self.current_molecular_temperature) / 2.0
        
        cavity_error = cavity_effective_temp - self.cavity_target_temperature
        molecular_error = molecular_effective_temp - self.molecular_target_temperature
        
        # Update integral error accumulators
        self.cavity_error_integral += cavity_error * dt_ps
        self.molecular_error_integral += molecular_error * dt_ps
        
        # Apply asymmetric gains based on error direction
        # Error < 0: System too cold, need to HEAT (raise bath temperature) 
        # Error > 0: System too hot, need to COOL (lower bath temperature)
        
        # Cavity asymmetric gains
        if cavity_error < 0:  # Need to heat system
            cavity_p_gain_factor = self.cavity_heating_gain_factor
            cavity_tau_i = (self.cavity_integral_heating_time_constant_ps 
                           if self.cavity_integral_heating_time_constant_ps is not None 
                           else self.cavity_integral_time_constant_ps)
        else:  # Need to cool system
            cavity_p_gain_factor = self.cavity_cooling_gain_factor
            cavity_tau_i = (self.cavity_integral_cooling_time_constant_ps 
                           if self.cavity_integral_cooling_time_constant_ps is not None 
                           else self.cavity_integral_time_constant_ps)
        
        # Molecular asymmetric gains
        if molecular_error < 0:  # Need to heat system
            molecular_p_gain_factor = self.molecular_heating_gain_factor
            molecular_tau_i = (self.molecular_integral_heating_time_constant_ps 
                              if self.molecular_integral_heating_time_constant_ps is not None 
                              else self.molecular_integral_time_constant_ps)
        else:  # Need to cool system
            molecular_p_gain_factor = self.molecular_cooling_gain_factor
            molecular_tau_i = (self.molecular_integral_cooling_time_constant_ps 
                              if self.molecular_integral_cooling_time_constant_ps is not None 
                              else self.molecular_integral_time_constant_ps)
        
        # Calculate asymmetric gradients (proportional terms)
        cavity_gradient = 0.5 * cavity_error * cavity_p_gain_factor
        molecular_gradient = 0.5 * molecular_error * molecular_p_gain_factor
        
        # Calculate asymmetric integral terms
        if cavity_tau_i is not None:
            cavity_alpha_i_effective = dt_ps / cavity_tau_i
            cavity_integral_term = cavity_alpha_i_effective * self.cavity_error_integral
        else:
            cavity_integral_term = 0.0
            
        if molecular_tau_i is not None:
            molecular_alpha_i_effective = dt_ps / molecular_tau_i
            molecular_integral_term = molecular_alpha_i_effective * self.molecular_error_integral
        else:
            molecular_integral_term = 0.0
        
        # Apply asymmetric PI control updates
        cavity_raw_output = (self.current_cavity_temperature - 
                            self.cavity_alpha * cavity_gradient - 
                            cavity_integral_term)
        molecular_raw_output = (self.current_molecular_temperature - 
                               self.molecular_alpha * molecular_gradient - 
                               molecular_integral_term)
        
        # Apply saturation limits
        if self.cavity_T_max is not None:
            cavity_output = max(self.cavity_T_min, min(self.cavity_T_max, cavity_raw_output))
        else:
            cavity_output = max(self.cavity_T_min, cavity_raw_output)
        
        if self.molecular_T_max is not None:
            molecular_output = max(self.molecular_T_min, min(self.molecular_T_max, molecular_raw_output))
        else:
            molecular_output = max(self.molecular_T_min, molecular_raw_output)
        
        # Update current temperatures and thermostats
        self.current_cavity_temperature = cavity_output
        self.current_molecular_temperature = molecular_output
        
        cavity_updated = self._update_thermostat(self.cavity_thermostat, cavity_output)
        molecular_updated = self._update_thermostat(self.molecular_thermostat, molecular_output)
        
        # Console output (periodic)
        should_print_console = (current_time_ps - self.last_console_output_time) >= self.console_output_period_ps
        
        if (cavity_updated or molecular_updated) and should_print_console:
            print(f"Dual Controller | Cavity: {self.cavity_method}→{cavity_output:.1f}K | "
                  f"Molecular: {self.molecular_method}→{molecular_output:.1f}K")
            self.last_console_output_time = current_time_ps
        
        # Log data
        try:
            with open(self.output_file, 'a', encoding='utf-8') as f:
                f.write(f"{current_time_ps:.6f},{cavity_measured_temp:.6f},{molecular_measured_temp:.6f},")
                f.write(f"{cavity_effective_temp:.6f},{molecular_effective_temp:.6f},")
                f.write(f"{self.cavity_target_temperature:.6f},{self.molecular_target_temperature:.6f},")
                f.write(f"{cavity_error:.6f},{molecular_error:.6f},")
                f.write(f"{cavity_gradient:.6f},{molecular_gradient:.6f},")
                f.write(f"{self.cavity_error_integral:.6f},{self.molecular_error_integral:.6f},")
                f.write(f"{cavity_integral_term:.6f},{molecular_integral_term:.6f},")
                f.write(f"{cavity_p_gain_factor:.6f},{molecular_p_gain_factor:.6f},")
                f.write(f"{cavity_tau_i if cavity_tau_i is not None else -1:.6f},{molecular_tau_i if molecular_tau_i is not None else -1:.6f},")
                f.write(f"{cavity_output:.6f},{molecular_output:.6f},{int(self.is_active)}\n")
        except Exception as e:
            print(f"Warning: Failed to write dual feedback output: {e}")


