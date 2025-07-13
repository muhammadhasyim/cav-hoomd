# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Main simulation framework for cavity molecular dynamics."""

import hoomd
import numpy as np
from .utils import PhysicalConstants


class AdaptiveTimestepUpdater(hoomd.custom.Action):
    """Update timestep adaptively based on energy conservation and stability."""
    
    def __init__(self, state, integrator, error_tolerance, time_constant_ps=50.0, 
                 initial_fraction=0.01, adaptiveerror=True, cavity_damping_factor=1.0, 
                 molecular_thermostat_tau=5.0, cavity_thermostat_tau=5.0, time_tracker=None,
                 switch_time_ps=None, timestep_change_threshold=0.1, max_timestep_change_factor=1.5,
                 shock_dampening_factor=1e-3, shock_dampening_enabled=True):
        super().__init__()
        
        self.state = state
        self.integrator = integrator
        self.target_error_tolerance = error_tolerance
        self.initial_error_tolerance = error_tolerance * initial_fraction
        self.time_constant_ps = time_constant_ps
        self.accumulated_time_ps = 0.0  # Fallback for backward compatibility
        self.last_timestep = 0  # Will be set correctly in first act() call
        self.adaptiveerror = adaptiveerror
        self.cavity_damping_factor = cavity_damping_factor
        self.molecular_thermostat_tau = molecular_thermostat_tau
        self.cavity_thermostat_tau = cavity_thermostat_tau
        self.time_tracker = time_tracker  # Reference to ElapsedTimeTracker for accurate timing
        self.switch_time_ps = switch_time_ps  # Time when ramping should start
        
        # NEW: Shock dampening parameters
        self.shock_dampening_factor = shock_dampening_factor
        self.shock_dampening_enabled = shock_dampening_enabled
        
        # Calculate shock dampening error tolerance
        if shock_dampening_enabled and switch_time_ps is not None:
            self.shock_error_tolerance = error_tolerance * shock_dampening_factor
        else:
            self.shock_error_tolerance = self.initial_error_tolerance
        
        # Initialize current_error_tolerance based on switch time behavior
        if switch_time_ps is not None:
            # Start with target tolerance before switch
            self.current_error_tolerance = self.target_error_tolerance
        else:
            # Original behavior: start with low tolerance for immediate ramping
            self.current_error_tolerance = self.initial_error_tolerance
        
        # Conservative timestep update parameters
        self.timestep_change_threshold = timestep_change_threshold  # Minimum relative change to trigger update
        self.max_timestep_change_factor = max_timestep_change_factor  # Maximum factor by which timestep can change
        self.last_timestep_update = 0  # Track when we last updated the timestep
        self.min_update_interval = 1000  # Minimum steps between timestep updates
        
        # NEW: Switch detection variables
        self.last_elapsed_time_ps = 0.0  # Track previous elapsed time to detect switch crossing
        self.switch_detected = False  # Flag to track if we've detected the switch
        self.switch_detection_tolerance = 0.001  # ps - tolerance for detecting switch crossing
        
        # Log the shock dampening behavior
        if shock_dampening_enabled and switch_time_ps is not None:
            print(f"Shock dampening enabled with time constant {time_constant_ps} ps")
            print(f"Before switch: error_tolerance = {self.target_error_tolerance:.2e} (normal tolerance)")
            print(f"At switch: error_tolerance drops to {self.shock_error_tolerance:.2e} (shock dampening factor: {shock_dampening_factor:.2e})")
            print(f"After switch: error_tolerance ramps from {self.shock_error_tolerance:.2e} to {self.target_error_tolerance:.2e} with τ = {time_constant_ps} ps")
            print(f"Switch detection: immediate timestep adjustment within {self.switch_detection_tolerance} ps of switch")
        elif switch_time_ps is not None:
            print(f"Error tolerance ramping will start at t = {self.switch_time_ps} ps")
            print(f"Before switch: error_tolerance = {self.target_error_tolerance:.2e} (final tolerance for efficiency)")
            print(f"At switch: error_tolerance drops to {self.initial_error_tolerance:.2e} (high precision)")
            print(f"After switch: error_tolerance ramps from {self.initial_error_tolerance:.2e} to {self.target_error_tolerance:.2e} with τ = {time_constant_ps} ps")
            print(f"Switch detection: immediate timestep adjustment within {self.switch_detection_tolerance} ps of switch")
        else:
            print("Error tolerance ramping starts immediately from t = 0 ps")
        
        # Log conservative timestep parameters
        print(f"Conservative timestep parameters:")
        print(f"  Change threshold: {self.timestep_change_threshold:.1%}")
        print(f"  Max change factor: {self.max_timestep_change_factor:.1f}")
        print(f"  Min update interval: {self.min_update_interval} steps")
        print(f"  Switch detection tolerance: {self.switch_detection_tolerance} ps")

    def act(self, timestep):
        """
        Custom action to update the timestep size.

        Parameters:
        - timestep: Current simulation timestep.
        """
        # Initialize last_timestep on first call to handle resuming from checkpoints
        if self.last_timestep == 0:
            self.last_timestep = timestep
        
        # Update accumulated simulation time (fallback method for backward compatibility)
        if timestep > self.last_timestep:
            # Convert dt to picoseconds using correct conversion
            dt_ps = PhysicalConstants.atomic_units_to_ps(self.integrator.dt)
            self.accumulated_time_ps += (timestep - self.last_timestep) * dt_ps
        self.last_timestep = timestep
        
        # Get accurate elapsed time for error tolerance ramping
        if self.time_tracker is not None:
            current_elapsed_time_ps = self.time_tracker.elapsed_time
        else:
            current_elapsed_time_ps = self.accumulated_time_ps
        
        # NEW: Detect if we've just crossed the switch time
        force_immediate_update = False
        if (self.switch_time_ps is not None and not self.switch_detected and
            self.last_elapsed_time_ps < self.switch_time_ps <= current_elapsed_time_ps):
            
            # We've just crossed the switch time!
            self.switch_detected = True
            force_immediate_update = True
            print(f"SWITCH DETECTED at timestep {timestep}: t = {current_elapsed_time_ps:.6f} ps")
            print(f"  Forcing immediate timestep adjustment for shock dampening")
        
        # Update last elapsed time for next iteration
        self.last_elapsed_time_ps = current_elapsed_time_ps
        
        # Update error tolerance based on shock dampening and exponential ramping
        if self.adaptiveerror:
            if self.switch_time_ps is not None:
                if current_elapsed_time_ps < self.switch_time_ps:
                    # Before switch: use target error tolerance
                    self.current_error_tolerance = self.target_error_tolerance
                else:
                    # After switch: implement shock dampening and exponential recovery
                    ramping_time = current_elapsed_time_ps - self.switch_time_ps
                    
                    if self.shock_dampening_enabled:
                        # NEW: Shock dampening mode - exponential recovery from shock tolerance
                        exp_factor = np.exp(-ramping_time / self.time_constant_ps)
                        self.current_error_tolerance = self.target_error_tolerance - \
                                                      (self.target_error_tolerance - self.shock_error_tolerance) * exp_factor
                    else:
                        # Original mode - exponential recovery from initial tolerance
                        exp_factor = np.exp(-ramping_time / self.time_constant_ps)
                        self.current_error_tolerance = self.target_error_tolerance - \
                                                      (self.target_error_tolerance - self.initial_error_tolerance) * exp_factor
            else:
                # Original behavior: start ramping immediately from t=0
                exp_factor = np.exp(-current_elapsed_time_ps / self.time_constant_ps)
                self.current_error_tolerance = self.target_error_tolerance - \
                                              (self.target_error_tolerance - self.initial_error_tolerance) * exp_factor
        else:
            self.current_error_tolerance = self.target_error_tolerance
        
        # Check if we should update the timestep (normal interval OR forced immediate update)
        steps_since_last_update = timestep - self.last_timestep_update
        should_update_timestep = (force_immediate_update or 
                                 steps_since_last_update >= self.min_update_interval)
        
        if not should_update_timestep:
            # Too soon to update timestep again (and no forced update)
            return
        
        # Collect forces and masses
        particle_data = self.state.get_snapshot().particles
        masses = np.array(particle_data.mass)
        n_particles = len(masses)
        
        # Initialize total forces array
        total_forces = np.zeros((n_particles, 3))
        
        # Sum forces from all force objects
        for force in self.integrator.forces:
            try:
                particle_forces = force.forces
                if particle_forces is not None and len(particle_forces) == n_particles:
                    total_forces += np.asarray(particle_forces)
            except (AttributeError, RuntimeError) as e:
                # Skip forces that can't be accessed (e.g., not yet computed)
                continue
        
        # Calculate sum |f_i| / m_i
        force_norm = np.array([np.linalg.norm(f) for f in total_forces])
        force_mass_sum = np.sum(force_norm / masses)
        
        # Compute optimal timestep using current error tolerance
        if force_mass_sum > 0:
            optimal_dt = np.sqrt(self.current_error_tolerance / force_mass_sum)
            current_dt = self.integrator.dt
            
            # Apply conservative timestep update logic
            dt_ratio = optimal_dt / current_dt
            
            # For forced updates (switch detected), be more aggressive with changes
            if force_immediate_update:
                # Allow larger changes when switch is detected for immediate stabilization
                effective_change_threshold = self.timestep_change_threshold * 0.1  # 10x more sensitive
                effective_max_change_factor = self.max_timestep_change_factor * 2.0  # Allow 2x larger changes
            else:
                # Normal conservative parameters
                effective_change_threshold = self.timestep_change_threshold
                effective_max_change_factor = self.max_timestep_change_factor
            
            # Only update if the change is significant (or forced)
            if force_immediate_update or abs(dt_ratio - 1.0) > effective_change_threshold:
                # Limit the maximum change to prevent dramatic jumps
                if dt_ratio > effective_max_change_factor:
                    new_dt = current_dt * effective_max_change_factor
                    clamped = True
                elif dt_ratio < (1.0 / effective_max_change_factor):
                    new_dt = current_dt / effective_max_change_factor
                    clamped = True
                else:
                    new_dt = optimal_dt
                    clamped = False
                
                # Apply the timestep change
                self.integrator.dt = new_dt
                self.last_timestep_update = timestep
                
                # Log the timestep change (always log for forced updates, otherwise reduced frequency)
                if force_immediate_update or timestep % 5000 == 0 or clamped:
                    update_type = "FORCED (switch)" if force_immediate_update else "regular"
                    print(f"Timestep updated at step {timestep} ({update_type}): {PhysicalConstants.atomic_units_to_ps(current_dt)*1e15:.1f} → {PhysicalConstants.atomic_units_to_ps(new_dt)*1e15:.1f} fs (error_tol: {self.current_error_tolerance:.2e})")
                    if clamped:
                        print(f"  CLAMPED by max change factor")
                    if force_immediate_update:
                        print(f"  Switch-triggered update: error_tolerance dropped to {self.current_error_tolerance:.2e}")
                
                # Update thermostat parameters only when timestep actually changes
                self._update_thermostat_parameters()
            else:
                # Timestep change is too small, don't update
                if timestep % 5000 == 0:  # Less frequent logging
                    print(f"Timestep change too small at step {timestep}: {dt_ratio:.3f} (threshold: {effective_change_threshold:.3f})")
        else:
            if timestep % 5000 == 0:
                print(f"WARNING: Zero force detected at step {timestep} - keeping current timestep")

    def _update_thermostat_parameters(self):
        """Update thermostat parameters when timestep changes."""
        # Convert thermostat time constants from ps to atomic units
        molecular_tau_au = PhysicalConstants.ps_to_atomic_units(self.molecular_thermostat_tau)
        cavity_tau_au = PhysicalConstants.ps_to_atomic_units(self.cavity_thermostat_tau)
        
        # Update thermostat parameters for both molecular and cavity methods
        # methods[0] is always molecular method, methods[1] is cavity method (if present)
        
        # Update molecular method thermostat parameters
        molecular_method = self.integrator.methods[0]
        if hasattr(molecular_method, 'default_gamma'):  # Langevin thermostat
            molecular_method.default_gamma = PhysicalConstants.gamma_from_tau_ps(self.molecular_thermostat_tau)
        elif hasattr(molecular_method, 'thermostat'):  # Bussi or MTTK thermostat
            if hasattr(molecular_method.thermostat, 'tau'):
                if hasattr(molecular_method.thermostat, 'translational_dof'):  # MTTK thermostat
                    molecular_method.thermostat.tau = molecular_tau_au  # Use configurable tau
                else:  # Bussi thermostat
                    molecular_method.thermostat.tau = molecular_tau_au  # Use configurable tau
        
        # Update cavity method thermostat parameters if present
        if len(self.integrator.methods) > 1:
            cavity_method = self.integrator.methods[1]
            if hasattr(cavity_method, 'default_gamma'):  # Langevin or Brownian thermostat
                if hasattr(cavity_method, 'gamma'):  # Brownian dynamics
                    # Brownian dynamics has per-type gamma parameters - apply damping factor
                    base_gamma = PhysicalConstants.gamma_from_tau_ps(self.cavity_thermostat_tau)
                    cavity_method.default_gamma = self.cavity_damping_factor * base_gamma
                else:  # Langevin thermostat
                    # Apply damping factor to base gamma calculated from tau
                    base_gamma = PhysicalConstants.gamma_from_tau_ps(self.cavity_thermostat_tau)
                    cavity_method.default_gamma = self.cavity_damping_factor * base_gamma
            elif hasattr(cavity_method, 'thermostat'):  # Bussi or MTTK thermostat
                if hasattr(cavity_method.thermostat, 'tau'):
                    if hasattr(cavity_method.thermostat, 'translational_dof'):  # MTTK thermostat
                        cavity_method.thermostat.tau = cavity_tau_au  # Use configurable tau
                    else:  # Bussi thermostat
                        cavity_method.thermostat.tau = cavity_tau_au  # Use configurable tau

    @property
    def error_tolerance(self):
        """Access the current effective error tolerance."""
        return self.current_error_tolerance

    @hoomd.logging.log
    def elapsed_time_ps(self):
        """Log the elapsed simulation time in picoseconds."""
        # Use accurate time from time_tracker if available, otherwise fallback
        if self.time_tracker is not None:
            return self.time_tracker.elapsed_time
        else:
            return self.accumulated_time_ps