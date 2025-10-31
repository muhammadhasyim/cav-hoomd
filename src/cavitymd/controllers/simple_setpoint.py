# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""
Simple Setpoint Temperature Controller for cavity MD simulations.

This controller captures a signal temperature as a fixed setpoint at activation time,
then uses a simple first-order control law to drive kinetic temperature to that setpoint.
"""

from typing import Optional
import hoomd
import numpy as np
import csv
from pathlib import Path
from ..analysis.timing import ElapsedTimeTracker
from .empirical import EmpiricalTemperatureData
from ..utils import PhysicalConstants


class SimpleSetpointController(hoomd.custom.Action):
    """
    Simple setpoint-based temperature controller for cavity MD simulations.
    
    This controller implements a simple first-order control law:
    
    dT_bath/dt = -(T_kinetic - T_setpoint) / τ
    
    where:
    - T_kinetic is the measured kinetic temperature
    - T_setpoint is captured once from a signal temperature at turn-on time
    - τ is the time constant controlling response speed
    
    The setpoint is captured at turn-on time from the specified signal method
    (kinetic, lj_coulombic, or harmonic) and remains fixed throughout the simulation.
    
    Parameters
    ----------
    signal_method : str
        Temperature method for setpoint capture ('kinetic', 'lj_coulombic', 'harmonic')
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
    output_file : str, optional
        Output file for control data. Default: 'simple_setpoint_control.csv'
    empirical_data_file : str, optional
        Path to empirical data file (required for lj_coulombic method)
    console_output_period_ps : float, optional
        Console output period in picoseconds. Default: 1.0
    
    Examples
    --------
    Basic control using LJ+Coulombic setpoint:
    
    >>> simple_controller = SimpleSetpointController(
    ...     signal_method='lj_coulombic',
    ...     time_constant_ps=5.0,
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     apply_to='molecular',
    ...     empirical_data_file='lj_coulombic_calibration.txt'
    ... )
    
    Control with kinetic temperature setpoint:
    
    >>> simple_controller = SimpleSetpointController(
    ...     signal_method='kinetic',
    ...     time_constant_ps=2.0,
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     apply_to='both'
    ... )
    """
    
    def __init__(self,
                 signal_method: str,
                 time_constant_ps: float,
                 time_tracker: ElapsedTimeTracker,
                 energy_tracker,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 apply_to: str = 'both',
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 T_min: float = 0.0,
                 T_max: Optional[float] = None,
                 output_file: str = 'simple_setpoint_control.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0,
                 enable_csv_output: bool = False):
        
        # Validate signal method
        valid_methods = ['kinetic', 'lj_coulombic', 'harmonic']
        if signal_method not in valid_methods:
            raise ValueError(f"signal_method must be one of {valid_methods}, got '{signal_method}'")
        
        # Validate apply_to
        valid_apply_to = ['molecular', 'cavity', 'both']
        if apply_to not in valid_apply_to:
            raise ValueError(f"apply_to must be one of {valid_apply_to}, got '{apply_to}'")
        
        # Store core parameters
        self.signal_method = signal_method
        self.time_constant_ps = time_constant_ps
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.apply_to = apply_to
        self.turn_on_time_ps = turn_on_time_ps
        self.turn_off_time_ps = turn_off_time_ps
        self.update_interval_ps = update_interval_ps
        self.T_min = T_min
        self.T_max = T_max
        self.output_file = output_file
        self.console_output_period_ps = console_output_period_ps
        self.enable_csv_output = enable_csv_output
        
        # Control history for HDF5 output
        self.control_history = []
        
        # Load empirical data if needed
        self.empirical_data = None
        if signal_method == 'lj_coulombic' and empirical_data_file is not None:
            try:
                self.empirical_data = EmpiricalTemperatureData(empirical_data_file)
                print(f"SimpleSetpointController: Loaded empirical data from {empirical_data_file}")
            except Exception as e:
                print(f"Warning: Could not load empirical data from {empirical_data_file}: {e}")
        
        # Controller state
        self.is_active = False
        self.setpoint_captured = False
        self.setpoint_temperature = None
        self.current_molecular_temperature = None
        self.current_cavity_temperature = None
        self.last_update_time = 0.0
        self.last_console_output_time = 0.0
        
        # Initialize output file (only if CSV enabled)
        if self.enable_csv_output:
            self._initialize_output_file()
        else:
            self.output_file = None
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([
                    'time_ps',
                    'T_kinetic_K', 
                    'T_setpoint_K',
                    'T_bath_molecular_K',
                    'T_bath_cavity_K',
                    'error_K',
                    'dT_bath_dt_K_per_ps',
                    'controller_active'
                ])
            print(f"SimpleSetpointController: Initialized output file {self.output_file}")
        except Exception as e:
            print(f"Warning: Could not initialize output file {self.output_file}: {e}")
    
    def _get_hoomd_simulation(self):
        """Get the HOOMD simulation object."""
        try:
            # Navigate through the object hierarchy to find simulation
            # This matches the pattern used in other controllers
            if hasattr(self, '_simulation'):
                return self._simulation
            
            # Try to find it through thermostats
            if self.molecular_thermostat is not None:
                return self.molecular_thermostat._simulation
            elif self.cavity_thermostat is not None:
                return self.cavity_thermostat._simulation
            
            # Try to find it through hoomd.context
            import hoomd
            if hasattr(hoomd, 'context') and hasattr(hoomd.context, 'current'):
                return hoomd.context.current.system_definition
                
        except:
            return None
    
    def _calculate_kinetic_temperature(self) -> Optional[float]:
        """Calculate kinetic temperature from particle velocities."""
        try:
            kinetic_energy = 0.0
            N_molecules = 0
            
            sim = self._get_hoomd_simulation()
            if sim is None:
                return None
                
            with sim.state.cpu_local_snapshot as snap:
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
    
    def _calculate_lj_coulombic_temperature(self) -> Optional[float]:
        """Calculate LJ+Coulombic fictive temperature using empirical data."""
        if self.empirical_data is None:
            print("Warning: No empirical data available for LJ+Coulombic temperature")
            return None
            
        try:
            energy_data = self.energy_tracker.get_instantaneous_energy()
            
            # Get LJ energy
            lj_energy = energy_data.get('lj', 0.0)
            if 'lj' not in energy_data:
                print(f"Warning: LJ energy not available. Available keys: {list(energy_data.keys())}")
                return None
            
            # Get Coulombic energy (could be 'coulomb', 'ewald_short'+'ewald_long', or other variations)
            coulombic_energy = 0.0
            if 'coulomb' in energy_data:
                coulombic_energy = energy_data['coulomb']
            elif 'ewald_short' in energy_data and 'ewald_long' in energy_data:
                coulombic_energy = energy_data['ewald_short'] + energy_data['ewald_long']
            elif 'ewald' in energy_data:
                coulombic_energy = energy_data['ewald']
            else:
                print(f"Warning: Coulombic energy not available. Available keys: {list(energy_data.keys())}")
                return None
            
            lj_coulombic_energy = lj_energy + coulombic_energy
            temperature_K = self.empirical_data.calculate_systemic_temperature(lj_coulombic_energy)
            
            if temperature_K is not None:
                return max(temperature_K, 0.0)
            else:
                return None
                
        except Exception as e:
            print(f"Warning: Could not calculate LJ+Coulombic temperature: {e}")
            return None
    
    def _calculate_harmonic_temperature(self) -> Optional[float]:
        """Calculate harmonic temperature using empirical data."""
        if self.empirical_data is None:
            print("Warning: No empirical data available for harmonic temperature")
            return None
            
        try:
            energy_data = self.energy_tracker.get_instantaneous_energy()
            if 'harmonic' not in energy_data:
                print("Warning: Harmonic energy not available")
                return None
            
            harmonic_energy = energy_data['harmonic']
            temperature_K = self.empirical_data.calculate_systemic_temperature(harmonic_energy)
            
            if temperature_K is not None:
                return max(temperature_K, 0.0)
            else:
                return None
                
        except Exception as e:
            print(f"Warning: Could not calculate harmonic temperature: {e}")
            return None
    
    def _measure_signal_temperature(self) -> Optional[float]:
        """Measure temperature using the specified signal method."""
        if self.signal_method == 'kinetic':
            return self._calculate_kinetic_temperature()
        elif self.signal_method == 'lj_coulombic':
            return self._calculate_lj_coulombic_temperature()
        elif self.signal_method == 'harmonic':
            return self._calculate_harmonic_temperature()
        else:
            print(f"Warning: Unknown signal method: {self.signal_method}")
            return None
    
    def _update_thermostats(self, T_molecular: Optional[float], T_cavity: Optional[float]):
        """Update thermostat temperatures."""
        try:
            if self.apply_to in ['molecular', 'both'] and T_molecular is not None:
                if self.molecular_thermostat is not None:
                    self.molecular_thermostat.kT = T_molecular * PhysicalConstants.KB_HARTREE_PER_K
                    self.current_molecular_temperature = T_molecular
            
            if self.apply_to in ['cavity', 'both'] and T_cavity is not None:
                if self.cavity_thermostat is not None:
                    self.cavity_thermostat.kT = T_cavity * PhysicalConstants.KB_HARTREE_PER_K
                    self.current_cavity_temperature = T_cavity
                    
        except Exception as e:
            print(f"Warning: Failed to update thermostat temperatures: {e}")
    
    def _log_to_csv(self, current_time_ps: float, T_kinetic: Optional[float], 
                    error: float, dT_dt: float):
        """Log current state to CSV file."""
        try:
            with open(self.output_file, 'a', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([
                    f"{current_time_ps:.3f}",
                    f"{T_kinetic:.3f}" if T_kinetic is not None else "N/A",
                    f"{self.setpoint_temperature:.3f}" if self.setpoint_temperature is not None else "N/A",
                    f"{self.current_molecular_temperature:.3f}" if self.current_molecular_temperature is not None else "N/A",
                    f"{self.current_cavity_temperature:.3f}" if self.current_cavity_temperature is not None else "N/A",
                    f"{error:.3f}",
                    f"{dT_dt:.6f}",
                    str(self.is_active)
                ])
        except Exception as e:
            print(f"Warning: Could not log to CSV: {e}")
    
    def act(self, timestep):
        """Execute simple setpoint control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should be active
        should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                           (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
        
        if not should_be_active:
            # Controller is dormant
            if self.is_active:
                self.is_active = False
                print(f"SimpleSetpointController turned OFF at t = {current_time_ps:.2f} ps")
            return
        
        # Controller should be active
        if not self.is_active:
            self.is_active = True
            print(f"\n{'='*70}")
            print(f"SimpleSetpointController activated at t = {current_time_ps:.2f} ps")
            print(f"Signal method: {self.signal_method}")
            print(f"Time constant: {self.time_constant_ps:.2f} ps")
            print(f"Apply to: {self.apply_to}")
            print(f"{'='*70}")
        
        # Capture setpoint at first activation
        if not self.setpoint_captured:
            setpoint_temp = self._measure_signal_temperature()
            if setpoint_temp is not None:
                self.setpoint_temperature = setpoint_temp
                self.setpoint_captured = True
                print(f"SETPOINT CAPTURED: T_setpoint = {self.setpoint_temperature:.2f} K (from {self.signal_method})")
            else:
                print(f"Warning: Could not capture setpoint from {self.signal_method}, retrying...")
                return
        
        # Check update interval
        dt_since_last = current_time_ps - self.last_update_time
        if dt_since_last < self.update_interval_ps:
            return
        
        dt_ps = dt_since_last if self.last_update_time > 0 else self.update_interval_ps
        self.last_update_time = current_time_ps
        
        # Measure current kinetic temperature
        T_kinetic = self._calculate_kinetic_temperature()
        if T_kinetic is None:
            print("Warning: Could not measure kinetic temperature, skipping control update")
            return
        
        # Calculate error and control output
        error = T_kinetic - self.setpoint_temperature
        dT_bath_dt = -error / self.time_constant_ps
        
        # Calculate new bath temperatures
        T_molecular_new = None
        T_cavity_new = None
        
        if self.apply_to in ['molecular', 'both']:
            current_T_mol = self.current_molecular_temperature if self.current_molecular_temperature is not None else T_kinetic
            T_molecular_new = current_T_mol + dT_bath_dt * dt_ps
            # Apply bounds
            if self.T_max is not None:
                T_molecular_new = max(self.T_min, min(self.T_max, T_molecular_new))
            else:
                T_molecular_new = max(self.T_min, T_molecular_new)
        
        if self.apply_to in ['cavity', 'both']:
            current_T_cav = self.current_cavity_temperature if self.current_cavity_temperature is not None else T_kinetic
            T_cavity_new = current_T_cav + dT_bath_dt * dt_ps
            # Apply bounds
            if self.T_max is not None:
                T_cavity_new = max(self.T_min, min(self.T_max, T_cavity_new))
            else:
                T_cavity_new = max(self.T_min, T_cavity_new)
        
        # Update thermostats
        self._update_thermostats(T_molecular_new, T_cavity_new)
        
        # Store control history for HDF5
        self.control_history.append({
            'time_ps': current_time_ps,
            'signal_temperature': T_kinetic,
            'setpoint_temperature': self.setpoint_temperature,
            'bath_temperature': T_molecular_new if self.apply_to in ['molecular', 'both'] else T_cavity_new,
            'control_active': True
        })
        
        # Console output
        if current_time_ps - self.last_console_output_time >= self.console_output_period_ps:
            print(f"SimpleSetpoint | t={current_time_ps:.1f}ps | T_kinetic={T_kinetic:.1f}K | T_setpoint={self.setpoint_temperature:.1f}K | error={error:+.1f}K | dT/dt={dT_bath_dt:+.3f}K/ps")
            if T_molecular_new is not None:
                print(f"  Molecular bath: {T_molecular_new:.1f}K")
            if T_cavity_new is not None:
                print(f"  Cavity bath: {T_cavity_new:.1f}K")
            self.last_console_output_time = current_time_ps
        
        # Log to CSV (only if enabled)
        if self.enable_csv_output:
            self._log_to_csv(current_time_ps, T_kinetic, error, dT_bath_dt)
