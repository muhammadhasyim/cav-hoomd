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
                 switch_time_ps=None):
        super().__init__()
        print("Performing error tolerance ramping with time constant", time_constant_ps, "ps")
        self.state = state
        self.integrator = integrator
        self.target_error_tolerance = error_tolerance
        self.initial_error_tolerance = error_tolerance * initial_fraction
        # Initialize current_error_tolerance based on switch time behavior
        if switch_time_ps is not None:
            # Inverted behavior: start with high tolerance before switch
            self.current_error_tolerance = self.target_error_tolerance
        else:
            # Original behavior: start with low tolerance for immediate ramping
            self.current_error_tolerance = self.initial_error_tolerance
        self.time_constant_ps = time_constant_ps
        self.accumulated_time_ps = 0.0  # Fallback for backward compatibility
        self.last_timestep = 0  # Will be set correctly in first act() call
        self.adaptiveerror = adaptiveerror
        self.cavity_damping_factor = cavity_damping_factor
        self.molecular_thermostat_tau = molecular_thermostat_tau
        self.cavity_thermostat_tau = cavity_thermostat_tau
        self.time_tracker = time_tracker  # Reference to ElapsedTimeTracker for accurate timing
        self.switch_time_ps = switch_time_ps  # Time when ramping should start
        
        # Log the ramping behavior
        if self.switch_time_ps is not None:
            print(f"Error tolerance ramping will start at t = {self.switch_time_ps} ps")
            print(f"Before switch: error_tolerance = {self.target_error_tolerance:.2e} (final tolerance for efficiency)")
            print(f"At switch: error_tolerance drops to {self.initial_error_tolerance:.2e} (high precision)")
            print(f"After switch: error_tolerance ramps from {self.initial_error_tolerance:.2e} to {self.target_error_tolerance:.2e} with τ = {time_constant_ps} ps")
        else:
            print("Error tolerance ramping starts immediately from t = 0 ps")

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
        
        # DEBUG: Print debug information every 1000 steps
        if timestep % 1000 == 0:
            print(f"DEBUG AdaptiveTimestepUpdater at timestep {timestep}:")
            print(f"  current_elapsed_time_ps: {current_elapsed_time_ps:.6f}")
            print(f"  self.switch_time_ps: {self.switch_time_ps}")
            print(f"  self.adaptiveerror: {self.adaptiveerror}")
            print(f"  switch_time_ps is not None: {self.switch_time_ps is not None}")
            if self.switch_time_ps is not None:
                print(f"  current_elapsed_time_ps < switch_time_ps: {current_elapsed_time_ps < self.switch_time_ps}")
        
        # Update error tolerance based on exponential approach
        if self.adaptiveerror:
            if self.switch_time_ps is not None:
                # INVERTED BEHAVIOR: Use final tolerance before switch, then drop and ramp back up
                if current_elapsed_time_ps < self.switch_time_ps:
                    # Before switch: use final (target) error tolerance for efficiency
                    self.current_error_tolerance = self.target_error_tolerance
                    if timestep % 1000 == 0:
                        print(f"  Using BEFORE switch logic: error_tolerance = {self.current_error_tolerance:.6f} (final tolerance)")
                else:
                    # After switch: start ramping FROM reduced tolerance TO final tolerance
                    ramping_time = current_elapsed_time_ps - self.switch_time_ps
                    exp_factor = np.exp(-ramping_time / self.time_constant_ps)
                    self.current_error_tolerance = self.target_error_tolerance - \
                                                  (self.target_error_tolerance - self.initial_error_tolerance) * exp_factor
                    if timestep % 1000 == 0:
                        print(f"  Using AFTER switch logic: ramping_time = {ramping_time:.6f}, error_tolerance = {self.current_error_tolerance:.6f}")
            else:
                # Original behavior: start ramping immediately from t=0
                exp_factor = np.exp(-current_elapsed_time_ps / self.time_constant_ps)
                self.current_error_tolerance = self.target_error_tolerance - \
                                              (self.target_error_tolerance - self.initial_error_tolerance) * exp_factor
                if timestep % 1000 == 0:
                    print(f"  Using IMMEDIATE ramping logic: error_tolerance = {self.current_error_tolerance:.6f}")
        else:
            self.current_error_tolerance = self.target_error_tolerance
        
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
        
        # Update timestep using current error tolerance
        if force_mass_sum > 0:
            new_dt = np.sqrt(self.current_error_tolerance / force_mass_sum)
            
            self.integrator.dt = new_dt
            
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