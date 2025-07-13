# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Main simulation framework for cavity molecular dynamics."""

import hoomd
import numpy as np
import logging
import sys
import os
import time
from pathlib import Path
import gsd.hoomd
from hoomd.bussi_reservoir.thermostats import BussiReservoir as Bussi

from .utils import PhysicalConstants, unwrap_positions
from .forces import CavityForce
from .analysis import (
    Status, ElapsedTimeTracker, TimestepFormatter, FieldAutocorrelationTracker,
    CavityModeTracker, EnergyTracker, PerformanceTracker
)
from .updaters import CavityParticleDisplacer
from .variants import StepVariant


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
        """Return the elapsed time in picoseconds."""
        return self.accumulated_time_ps


class CavityMDSimulation:
    """
    A class to encapsulate cavity MD simulation setup and execution.
    
    Features smart cavity particle handling:
    - Automatically detects if cavity particles already exist in the GSD file
    - Only adds new cavity particles when none exist and add_cavity_particle=True
    - Validates existing cavity particles when found
    - Provides clear error messages for misconfigurations
    """
    
    def __init__(self, job_dir, replica, freq, couplstr, incavity, runtime_ps=500.0, 
                 input_gsd='molecular-0.gsd', frame=-1, name='prod', error_tolerance=0.01,
                 temperature=100.0, molecular_thermostat='bussi', cavity_thermostat='langevin',
                 cavity_damping_factor=1.0, use_brownian_overdamped=True, add_cavity_particle=True,
                 finite_q=False, molecular_thermostat_tau=5.0, cavity_thermostat_tau=5.0,
                 log_level='INFO', custom_log_file=None,
                 enable_fkt=True, fkt_kmag=1.0, fkt_num_wavevectors=50, fkt_reference_interval_ps=1.0, fkt_max_references=10,
                 max_energy_output_time_ps=None, enable_energy_tracking=False, dt_fs=None, device='CPU', gpu_id=0,
                 energy_output_period_ps=0.1, fkt_output_period_ps=1.0, gsd_output_period_ps=50.0, console_output_period_ps=1.0,
                 enable_text_output=False, text_output_file=None, truncate_gsd=False, seed=None, restart_velocities=True,
                 switch_time_ps=None, dissipation=0.0):
        """Initialize the CavityMDSimulation with simulation parameters."""
        self.job_dir = job_dir
        self.replica = replica
        self.freq = freq
        self.couplstr = couplstr
        self.incavity = incavity
        self.runtime_ps = runtime_ps
        self.input_gsd = input_gsd
        self.frame = frame
        self.name = name
        self.error_tolerance = error_tolerance
        self.temperature = temperature
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.cavity_damping_factor = cavity_damping_factor
        self.use_brownian_overdamped = use_brownian_overdamped
        self.add_cavity_particle = add_cavity_particle
        self.finite_q = finite_q
        self.molecular_thermostat_tau = molecular_thermostat_tau
        self.cavity_thermostat_tau = cavity_thermostat_tau
        
        # Time-varying parameters
        self.switch_time_ps = switch_time_ps
        self.dissipation = dissipation
        
        # Logging parameters - simplified to console only
        self.log_level = log_level
        self.custom_log_file = custom_log_file
        
        # F(k,t) calculation parameters
        self.enable_fkt = enable_fkt
        self.fkt_kmag = fkt_kmag
        self.fkt_num_wavevectors = fkt_num_wavevectors
        self.fkt_reference_interval_ps = fkt_reference_interval_ps
        self.fkt_max_references = fkt_max_references
        
        # Energy output limit parameter
        self.max_energy_output_time_ps = max_energy_output_time_ps
        
        # Energy tracking parameter
        self.enable_energy_tracking = enable_energy_tracking
        
        # Manual timestep parameter (in femtoseconds)
        self.dt_fs = dt_fs
        
        # Device configuration
        self.device = device.upper()
        self.gpu_id = gpu_id
        
        # Seed configuration
        self.seed = seed
        
        # Velocity restart configuration
        self.restart_velocities = restart_velocities
        
        # Physical constants
        self.kB = PhysicalConstants.KB_HARTREE_PER_K
        
        # Output control parameters - separate periods for different observables
        self.energy_output_period_ps = energy_output_period_ps
        self.fkt_output_period_ps = fkt_output_period_ps
        self.gsd_output_period_ps = gsd_output_period_ps
        self.console_output_period_ps = console_output_period_ps
        self.enable_text_output = enable_text_output
        self.text_output_file = text_output_file
        self.truncate_gsd = truncate_gsd
        
        # Initialize simulation components (will be set during setup)
        self.sim = None
        self.logger = None

    def setup_logging(self):
        """Set up simplified logging configuration for the simulation."""
        # Create a custom logger for this simulation
        logger_name = f"CavityMD_{self.name}_{self.replica}"
        self.logger = logging.getLogger(logger_name)
        self.logger.setLevel(getattr(logging, self.log_level.upper()))
        
        # Clear any existing handlers to avoid duplication
        self.logger.handlers.clear()
        
        # Create formatter for log messages
        formatter = logging.Formatter(
            '%(asctime)s | %(levelname)s | %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # Always set up console logging (simplified)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, self.log_level.upper()))
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        
        # Log initial setup information
        self.log_info("="*60)
        self.log_info("CAVITY MD SIMULATION STARTED")
        self.log_info("="*60)
        self.log_info(f"Simulation: {self.name}-{self.replica}")
        self.log_info(f"Device: {self.device}" + (f" (GPU {self.gpu_id})" if self.device == 'GPU' else ""))
        self.log_info(f"Runtime: {self.runtime_ps} ps")
        self.log_info(f"Temperature: {self.temperature} K")
        self.log_info(f"Random seed: {self.seed if self.seed is not None else 'replica-based'}")
        self.log_info(f"Cavity coupling: {'Enabled' if self.incavity else 'Disabled'}")
        if self.incavity:
            self.log_info(f"  Frequency: {self.freq} cm^-1")
            self.log_info(f"  Coupling strength: {self.couplstr}")
            self.log_info(f"  Finite-q mode: {self.finite_q}")
        self.log_info(f"Molecular thermostat: {self.molecular_thermostat} (tau={self.molecular_thermostat_tau} ps)")
        if self.incavity:
            self.log_info(f"Cavity thermostat: {self.cavity_thermostat} (tau={self.cavity_thermostat_tau} ps)")
        self.log_info("="*60)
    
    def log_info(self, message):
        """Log an info message."""
        if self.logger:
            self.logger.info(message)
        else:
            print(message)
    
    def log_warning(self, message):
        """Log a warning message."""
        if self.logger:
            self.logger.warning(message)
        else:
            print(f"WARNING: {message}")
    
    def log_error(self, message):
        """Log an error message."""
        if self.logger:
            self.logger.error(message)
        else:
            print(f"ERROR: {message}")

    def run(self):
        """Main orchestrator method that runs the complete simulation workflow."""
        try:
            # Phase 0: Setup logging
            self.setup_logging()
            
            # Phase 1: Setup simulation and state
            self.log_info("=== Phase 1: Setting up simulation ===")
            self.calculate_physical_parameters()
            snapshot = self.setup_simulation()

            # Phase 1.5: Setup time tracker (must be after sim setup, before force setup)
            self.log_info("=== Phase 1.5: Setting up time tracker ===")
            self.time_tracker = ElapsedTimeTracker(self.sim, self.runtime_ps)
            self.sim.operations.updaters.append(hoomd.update.CustomUpdater(
                action=self.time_tracker, trigger=hoomd.trigger.Periodic(1)
            ))
            
            # Phase 2: Configure forces and thermostats
            self.log_info("=== Phase 2: Configuring forces and thermostats ===")
            forces = self.setup_force_parameters(self.dt)
            molecular_method, cavity_method, thermostat_refs = self.setup_thermostat_parameters(self.dt)
            
            # Phase 3: Setup integrator and thermalization
            self.log_info("=== Phase 3: Setting up integrator and thermalization ===")
            methods = [molecular_method]
            if cavity_method is not None:
                methods.append(cavity_method)
            self.setup_integrator(forces, methods)
            self.thermalize_system()
            
            # Phase 3.5: Compute and set optimal timestep (after thermalization, before logging)
            self.log_info("=== Phase 3.5: Computing optimal timestep ===")
            self.compute_and_set_optimal_timestep()
            
            # Phase 4: Setup trackers and loggers
            self.log_info("=== Phase 4: Setting up trackers and loggers ===")
            self.setup_trackers_and_loggers()
            
            # Phase 5: Setup output writers
            self.log_info("=== Phase 5: Setting up output writers ===")
            self.setup_output_writers()
            
            # Phase 6: Run simulation
            self.log_info("=== Phase 6: Running simulation ===")
            self.run_simulation()
            
            # Phase 7: Cleanup
            self.log_info("=== Phase 7: Cleanup ===")
            self.cleanup()
            
            self.log_info("=== SIMULATION COMPLETED SUCCESSFULLY ===")
            return 0  # Success
            
        except Exception as e:
            self.log_error(f"CRITICAL ERROR in simulation: {str(e)}")
            import traceback
            self.log_error("Full traceback:")
            for line in traceback.format_exc().split('\n'):
                if line.strip():
                    self.log_error(line)
            self.cleanup()  # Try to cleanup even on error
            return 1  # Failure

    def cleanup(self):
        """Handle post-simulation cleanup and restore original directory."""
        # Note: Trackers from analysis.py write directly to files, no buffering needed
        self.log_info("Cleanup initiated...")
        
        # Restore original directory
        if hasattr(self, 'original_cwd'):
            os.chdir(self.original_cwd)
            self.log_info(f"Returned to original directory: {self.original_cwd}")
        
        self.log_info("Cleanup completed")

    def calculate_physical_parameters(self):
        """Calculate physical parameters and unit conversions."""
        # Determine the actual timestep to use for calculations
        if self.error_tolerance <= 0 and self.dt_fs is not None:
            # Fixed timestep mode with user-specified timestep
            dt_ps = self.dt_fs / 1000.0  # Convert fs to ps
            self.log_info(f"Using user-specified timestep for calculations: {self.dt_fs} fs = {dt_ps} ps")
        else:
            # Adaptive timestep mode or no user-specified timestep - use default
            dt_ps = 0.0001  # timestep in ps (0.1 fs - reasonable default for adaptive mode)
            self.log_info(f"Using default timestep for calculations: {dt_ps} ps")
        
        runtime_real = self.runtime_ps  # runtime in ps
        
        # Calculate different output periods in timesteps using the correct dt_ps
        energy_period = max(1, int(self.energy_output_period_ps / dt_ps))
        fkt_period = max(1, int(self.fkt_output_period_ps / dt_ps))
        gsd_period = max(1, int(self.gsd_output_period_ps / dt_ps))
        console_period = max(1, int(self.console_output_period_ps / dt_ps))
        
        # Convert timestep from ps to atomic units using helper method
        dt_au = PhysicalConstants.ps_to_atomic_units(dt_ps)
        
        # Calculate total steps needed for the specified runtime
        # Note: This is only used for fixed timestep mode
        # For adaptive timestep, ElapsedTimeTracker handles runtime termination
        total_steps_needed = int(runtime_real / dt_ps)
        
        # Store converted values as instance variables for later use
        self.dt = dt_au  # timestep in atomic units (for HOOMD)
        self.dt_ps = dt_ps  # timestep in ps (for calculations)
        self.runtime = total_steps_needed  # total steps needed
        self.energy_period = energy_period
        self.fkt_period = fkt_period
        self.gsd_period = gsd_period
        self.console_period = console_period
        
        self.log_info(f"Time conversions:")
        self.log_info(f"  Timestep: {dt_ps} ps = {dt_au:.6f} a.u.")
        self.log_info(f"  Runtime: {self.runtime_ps:.1f} ps = {total_steps_needed} steps")
        self.log_info(f"  Steps per ps: {1.0/dt_ps:.1f}")
        self.log_info(f"Output periods (steps):")
        self.log_info(f"  Energy: {energy_period} ({self.energy_output_period_ps:.3f} ps)")
        self.log_info(f"  F(k,t): {fkt_period} ({self.fkt_output_period_ps:.3f} ps)")
        self.log_info(f"  GSD: {gsd_period} ({self.gsd_output_period_ps:.3f} ps)")
        self.log_info(f"  Console: {console_period} ({self.console_output_period_ps:.3f} ps)")
        
        return dt_au, total_steps_needed, energy_period, fkt_period, gsd_period, console_period

    def setup_device(self):
        """Setup the HOOMD device (CPU or GPU)."""
        if self.device == 'GPU':
            try:
                # Try different GPU initialization methods for different HOOMD versions
                try:
                    # First try with gpu_ids parameter (newer HOOMD versions)
                    device = hoomd.device.GPU(gpu_ids=[self.gpu_id])
                except TypeError:
                    # Fall back to older parameter name
                    device = hoomd.device.GPU(gpu_id=self.gpu_id)
                self.log_info(f"Initializing simulation on GPU {self.gpu_id}")
            except Exception as e:
                self.log_warning(f"Failed to initialize GPU {self.gpu_id}: {str(e)}")
                self.log_warning("Falling back to CPU")
                device = hoomd.device.CPU()
                self.device = 'CPU'  # Update device setting
        elif self.device == 'CPU':
            device = hoomd.device.CPU()
            self.log_info("Initializing simulation on CPU")
        else:
            raise ValueError(f"Invalid device '{self.device}'. Must be 'CPU' or 'GPU'")
        
        return device

    def check_cavity_particle_exists(self, snapshot):
        """Check if a cavity particle already exists in the snapshot."""
        # Check if 'L' particle type exists
        if 'L' not in snapshot.particles.types:
            return False, 0
        
        # Check if any particles have typeid == 2 (cavity particles)
        cavity_count = np.sum(snapshot.particles.typeid == 2)
        return cavity_count > 0, cavity_count

    def create_cavity_particle(self, snapshot):
        """
        Add a new cavity particle to the simulation snapshot.
        
        Note: This method assumes no cavity particle already exists in the snapshot.
        Use check_cavity_particle_exists() to verify this before calling.
        """
        self.log_info("Creating new cavity particle and adding to system...")
        
        # Set numpy seed for reproducible cavity particle positioning if seed is provided
        if self.seed is not None:
            np.random.seed(self.seed + 2)  # Use seed+2 to differentiate from other random operations
            self.log_info(f"Using seed {self.seed + 2} for cavity particle positioning")
        
        # Calculate dipole moment and photon position
        positions = unwrap_positions(snapshot.particles.position, snapshot.particles.image, 
                                   snapshot.configuration.box[:3])
        dipmom = np.einsum('i,ij->j', snapshot.particles.charge, positions)

        omegac = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1
        
        # Determine coupling strength for initial placement
        initial_couplstr = 0.0 if self.switch_time_ps is not None else self.couplstr

        if self.finite_q:
            # Allow finite-q photon displacement based on dipole-coupling interaction
            if self.switch_time_ps is not None:
                # For instantaneous switching: start at zero coupling position
                # Displacement will happen at switch time via CavityParticleDisplacer
                newpos = np.array([0.0, 0.0, 0.0])
                self.log_info(f"Finite-q + instantaneous switching: Photon starts at origin")
                self.log_info(f"  CavityParticleDisplacer will handle displacement at switch time")
            else:
                # Original finite-q behavior for constant coupling
                newpos = -dipmom * initial_couplstr / omegac**2
                newpos[-1] = 0.0
                # Only add thermal fluctuations if coupling is non-zero
                if initial_couplstr != 0.0:
                    sigma = np.sqrt(self.kB * self.temperature / omegac**2)
                    newpos = np.random.normal(loc=newpos, scale=sigma, size=3)
                    self.log_info(f"Finite-q mode: Photon displaced by dipole interaction to {newpos} (with thermal fluctuations)")
                else:
                    self.log_info(f"Finite-q mode: Photon at equilibrium position {newpos} (no thermal fluctuations due to zero coupling)")
        else:
            # Start photon at origin (q=0 limit)
            newpos = np.array([0.0, 0.0, 0.0])
            # Only add thermal fluctuations if coupling is non-zero
            if initial_couplstr != 0.0:
                sigma = np.sqrt(self.kB * self.temperature / omegac**2)
                newpos = np.random.normal(loc=newpos, scale=sigma, size=3)
                self.log_info("q=0 mode: Photon positioned at origin + thermal fluctuations")
            else:
                self.log_info("q=0 mode: Photon positioned exactly at origin (no thermal fluctuations due to zero coupling)")
        
        # Wrap position and get image flags
        def wrap_position(x, L):
            # Compute the image flags (how many box lengths away from the primary box)
            image_flags = np.floor((x + L/2) / L)
            # Compute the wrapped position inside the primary box
            wrapped_position = x - image_flags * L
            return wrapped_position, image_flags.astype(int)
            
        newpos, image_flags = wrap_position(newpos, np.array(snapshot.configuration.box[:3]))
        
        # Add photon particle
        if 'L' not in snapshot.particles.types:
            snapshot.particles.types.append('L')
        snapshot.particles.N += 1
        snapshot.particles.typeid = np.append(snapshot.particles.typeid, [2])
        snapshot.particles.position = np.append(snapshot.particles.position, [newpos], axis=0)
        snapshot.particles.charge = np.append(snapshot.particles.charge, [0.0])
        snapshot.particles.mass = np.append(snapshot.particles.mass, [1.0])
        snapshot.particles.diameter = np.append(snapshot.particles.diameter, [1.0])
        snapshot.particles.image = np.vstack([snapshot.particles.image, image_flags])

        # Set additional particle properties for the photon
        if hasattr(snapshot.particles, "body"):
            snapshot.particles.body = np.append(snapshot.particles.body, [-1], axis=0)

        if hasattr(snapshot.particles, "orientation"):
            snapshot.particles.orientation = np.append(
                snapshot.particles.orientation,
                [[1.0, 0.0, 0.0, 0.0]],
                axis=0
            )

        if hasattr(snapshot.particles, "moment_inertia"):
            snapshot.particles.moment_inertia = np.vstack([
                snapshot.particles.moment_inertia,
                np.zeros((1, 3))
            ])

        if hasattr(snapshot.particles, "velocity"):
            snapshot.particles.velocity = np.vstack([
                snapshot.particles.velocity,
                np.zeros((1, 3))
            ])

        if hasattr(snapshot.particles, "angmom"):
            snapshot.particles.angmom = np.vstack([
                snapshot.particles.angmom,
                np.zeros((1, 4))
            ])
            
        self.log_info(f"Cavity particle added at position {newpos}")
        return snapshot

    def validate_cavity_particle(self):
        """Validate that cavity particle exists when required."""
        # Get particle types from simulation state, not local snapshot
        particle_types = self.sim.state.particle_types
        
        if 'L' not in particle_types:
            raise ValueError("ERROR: Cavity simulation requested but no cavity particle type 'L' found in GSD file.")
        
        with self.sim.state.cpu_local_snapshot as snap:
            if 2 not in snap.particles.typeid:
                raise ValueError("ERROR: Cavity simulation requested but no cavity particles found in GSD file.")
            
            cavity_count = np.sum(snap.particles.typeid == 2)
            if cavity_count != 1:
                raise ValueError(f"ERROR: Expected exactly 1 cavity particle but found {cavity_count} in GSD file.")
            
            cavity_index = np.where(snap.particles.typeid == 2)[0][0]
            cavity_position = snap.particles.position[cavity_index]
            self.log_info(f"Cavity particle validated at index {cavity_index}, position {cavity_position}")

    def run_simulation(self):
        """Execute the main simulation loop."""
        # For adaptive timestep, use a very large number of steps
        # The ElapsedTimeTracker will exit when runtime is reached
        if self.error_tolerance > 0:
            # Adaptive timestep mode - use effectively infinite steps
            total_steps = 999999999  # Very large number - will be stopped by ElapsedTimeTracker
            self.log_info(f"Starting adaptive timestep simulation for {self.runtime_ps:.1f} ps")
            self.log_info(f"Using max steps: {total_steps} (will exit when runtime reached)")
        else:
            # Fixed timestep mode - use calculated steps
            total_steps = self.runtime
            self.log_info(f"Starting fixed timestep simulation for {self.runtime_ps:.1f} ps ({total_steps} steps)")
        
        # Get the actual timestep being used (may have been updated by adaptive timestep computation)
        actual_dt = self.sim.operations.integrator.dt
        actual_dt_ps = PhysicalConstants.atomic_units_to_ps(actual_dt)
        
        self.log_info(f"Initial timestep: {actual_dt:.6f} a.u. ({actual_dt_ps * 1000:.3f} fs)")
        if self.error_tolerance > 0:
            self.log_info(f"Timestep will adapt dynamically (error_tolerance = {self.error_tolerance})")
        else:
            self.log_info(f"Fixed timestep mode - steps per ps: {1.0/actual_dt_ps:.1f}")
        
        # Final check for Bussi thermostat compatibility - ensure non-zero velocities
        if self.molecular_thermostat.lower() == 'bussi' or self.cavity_thermostat.lower() == 'bussi':
            self.log_info("Performing final velocity check for Bussi thermostat compatibility...")
            kT = self.kB * self.temperature
            # self.ensure_nonzero_velocities(kT)  # No longer needed - proper thermalization already done
        
        # Run the simulation
        self.sim.run(total_steps, write_at_start=True)
        
        self.log_info(f"Simulation completed successfully")

    def setup_simulation(self):
        """Create HOOMD simulation object and initialize state from GSD file."""
        # Save current directory and change to job directory
        self.original_cwd = os.getcwd()
        os.chdir(self.job_dir)
        
        # Setup device
        device = self.setup_device()
        
        # Create simulation object with seed control
        if self.seed is not None:
            simulation_seed = self.seed
            self.log_info(f"Using user-specified seed: {simulation_seed}")
        else:
            # Generate deterministic seed based on replica for reproducibility
            simulation_seed = hash(str(self.replica)) % (2**31)  # Use replica-based seed
            self.log_info(f"Using replica-based seed: {simulation_seed} (replica {self.replica})")
        
        self.sim = hoomd.Simulation(device=device, seed=simulation_seed)
        
        # Load GSD file and handle cavity particle
        with gsd.hoomd.open(self.input_gsd, 'r') as f:
            if self.frame < 0:
                self.frame = len(f) + self.frame  # Convert negative index to positive
                if self.frame < 0:  # Handle case where abs(frame) > len(f)
                    self.frame = 0
            snapshot = f[self.frame]
            
            if self.incavity:
                # Check if cavity particle already exists
                cavity_exists, cavity_count = self.check_cavity_particle_exists(snapshot)
                
                if cavity_exists:
                    # Cavity particle already exists - use original snapshot
                    self.log_info(f"Cavity particle already exists in GSD file ({cavity_count} found)")
                    if cavity_count > 1:
                        self.log_warning(f"Multiple cavity particles found ({cavity_count}). Expected exactly 1.")
                    self.sim.create_state_from_gsd(filename=self.input_gsd, frame=self.frame)
                    self.log_info(f"Simulation state created from original GSD file with existing cavity particle")
                    # Validate the existing cavity particle
                    self.validate_cavity_particle()
                    
                elif self.add_cavity_particle:
                    # No cavity particle exists, but we want to add one
                    self.log_info("No cavity particle found in GSD file - adding new cavity particle")
                    snapshot = self.create_cavity_particle(snapshot)
                    self.sim.create_state_from_snapshot(snapshot)
                    self.log_info(f"Simulation state created from modified snapshot with new cavity particle")
                    
                else:
                    # No cavity particle exists and we don't want to add one - this is an error
                    raise ValueError(
                        "ERROR: Cavity simulation requested but no cavity particle found in GSD file "
                        "and add_cavity_particle=False. Either:\n"
                        "  1. Set add_cavity_particle=True to automatically add a cavity particle, or\n"
                        "  2. Use a GSD file that already contains a cavity particle (type 'L' with typeid=2)"
                    )
            else:
                # No cavity simulation - use original GSD file
                self.sim.create_state_from_gsd(filename=self.input_gsd, frame=self.frame)
                self.log_info(f"Simulation state created from original GSD file (no cavity)")
        
        return snapshot