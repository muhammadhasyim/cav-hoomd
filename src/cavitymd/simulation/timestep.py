# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Adaptive timestep updater for cavity MD simulations."""

import hoomd
import numpy as np

from ..utils import PhysicalConstants


class AdaptiveTimestepUpdater(hoomd.custom.Action):
    """Update timestep adaptively based on energy conservation and stability."""
    
    def __init__(self, state, integrator, error_tolerance, time_constant_ps=50.0, 
                 initial_fraction=0.01, adaptiveerror=True, cavity_damping_factor=1.0, 
                 molecular_thermostat_tau=5.0, cavity_thermostat_tau=5.0, time_tracker=None,
                 switch_time_ps=None, timestep_change_threshold=0.1, max_timestep_change_factor=1.5,
                 shock_dampening_factor=1e-3, shock_dampening_enabled=True, shock_dampening_time_constant_ps=50.0,
                 # New dynamic coupling change detection parameters
                 dynamic_coupling_detection=True, coupling_change_threshold=1e-5):
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
        self.switch_time_ps = switch_time_ps #Time when ramping should start
        
        # NEW: Shock dampening parameters
        self.shock_dampening_factor = shock_dampening_factor
        self.shock_dampening_enabled = shock_dampening_enabled
        
        # Calculate shock dampening error tolerance
        if shock_dampening_enabled and switch_time_ps is not None:
            self.shock_error_tolerance = error_tolerance * shock_dampening_factor
        else:
            self.shock_error_tolerance = self.initial_error_tolerance
        
        # Initialize current_error_tolerance based on shock dampening mode
        if dynamic_coupling_detection:
            # Dynamic mode: start with target tolerance, only drop on actual shocks
            self.current_error_tolerance = self.target_error_tolerance
        elif switch_time_ps is not None:
            # Legacy switch time mode: start with target tolerance before switch
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
        
        # NEW: Dynamic coupling change detection
        self.dynamic_coupling_detection = dynamic_coupling_detection
        self.coupling_change_threshold = coupling_change_threshold  # a.u. - threshold for large coupling changes
        self.last_coupling_value = None  # Track previous coupling strength
        self.coupling_shock_detected_time = None  # When we last detected a coupling shock
        self.coupling_recovery_active = False  # Whether we're currently in recovery mode
        self.first_coupling_check = True  # Flag to initialize without triggering shock
        
        # Log the shock dampening behavior
        if dynamic_coupling_detection:
            print(f"Dynamic coupling change detection enabled", flush=True)
            print(f"Coupling change threshold: ±{coupling_change_threshold:.2e} a.u.", flush=True)
            print(f"Shock dampening: error_tolerance drops to {self.shock_error_tolerance:.2e} (factor: {shock_dampening_factor:.2e})", flush=True)
            print(f"Recovery: exponential ramp to {self.target_error_tolerance:.2e} with τ = {time_constant_ps} ps", flush=True)
        elif shock_dampening_enabled and switch_time_ps is not None:
            print(f"Shock dampening enabled with time constant {time_constant_ps} ps", flush=True)
            print(f"Before switch: error_tolerance = {self.target_error_tolerance:.2e} (normal tolerance)", flush=True)
            print(f"At switch: error_tolerance drops to {self.shock_error_tolerance:.2e} (shock dampening factor: {shock_dampening_factor:.2e})", flush=True)
            print(f"After switch: error_tolerance ramps from {self.shock_error_tolerance:.2e} to {self.target_error_tolerance:.2e} with τ = {time_constant_ps} ps", flush=True)
            print(f"Switch detection: immediate timestep adjustment within {self.switch_detection_tolerance} ps of switch", flush=True)
        elif switch_time_ps is not None:
            print(f"Error tolerance ramping will start at t = {self.switch_time_ps} ps", flush=True)
            print(f"Before switch: error_tolerance = {self.target_error_tolerance:.2e} (final tolerance for efficiency)", flush=True)
            if self.shock_dampening_factor > 0:
                print(f"At switch: error_tolerance drops to {self.initial_error_tolerance:.2e} (high precision)", flush=True)
                print(f"After switch: error_tolerance ramps from {self.initial_error_tolerance:.2e} to {self.target_error_tolerance:.2e} with τ = {time_constant_ps} ps", flush=True)
            else:
                print(f"At switch: error_tolerance remains at {self.target_error_tolerance:.2e} (no ramping, initial_fraction = 0.0)", flush=True)
            print(f"Switch detection: immediate timestep adjustment within {self.switch_detection_tolerance} ps of switch", flush=True)
        else:
            print("Error tolerance ramping starts immediately from t = 0 ps", flush=True)
        
        # Log conservative timestep parameters
        print(f"Conservative timestep parameters:", flush=True)
        print(f"  Change threshold: {self.timestep_change_threshold:.1%}", flush=True)
        print(f"  Max change factor: {self.max_timestep_change_factor:.1f}", flush=True)
        print(f"  Min update interval: {self.min_update_interval} steps", flush=True)
        print(f"  Switch detection tolerance: {self.switch_detection_tolerance} ps", flush=True)
    
    def _get_current_coupling_strength(self):
        """Get the current coupling strength from cavity force."""
        try:
            # Look for cavity force in the integrator's forces
            for force in self.integrator.forces:
                if hasattr(force, 'couplstr_variant'):
                    # This is the cavity force - get the variant
                    coupling_variant = force.couplstr_variant
                    if hasattr(coupling_variant, '__call__'):
                        # It's a variant - call it with current timestep
                        # Access timestep through the integrator's simulation
                        if hasattr(self.integrator, 'simulation'):
                            current_timestep = self.integrator.simulation.timestep
                        else:
                            # Fallback - use a reasonable default
                            current_timestep = 0
                        return coupling_variant(current_timestep)
                    else:
                        # It's a constant value
                        return float(coupling_variant.value) if hasattr(coupling_variant, 'value') else float(coupling_variant)
            return None
        except Exception as e:
            # If we can't get coupling strength, return None
            return None

        # DEBUG: Show key variables for shock dampening
        #print(f"DEBUG: shock_dampening_factor = {shock_dampening_factor:.2e}", flush=True)
        #print(f"DEBUG: shock_dampening_enabled = {shock_dampening_enabled}", flush=True)
        #print(f"DEBUG: initial_error_tolerance = {self.initial_error_tolerance:.2e}", flush=True)
        #print(f"DEBUG: shock_error_tolerance = {self.shock_error_tolerance:.2e}", flush=True)

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
        
        # NEW: Detect sudden coupling changes OR switch time crossing
        force_immediate_update = False
        coupling_shock_detected = False
        
        # Dynamic coupling change detection
        if self.dynamic_coupling_detection:
            # Get current coupling value from cavity force
            current_coupling = self._get_current_coupling_strength()
            
            if self.first_coupling_check:
                # First time - just store the value without triggering shock
                self.last_coupling_value = current_coupling
                self.first_coupling_check = False
            elif self.last_coupling_value is not None and current_coupling is not None:
                coupling_change = abs(current_coupling - self.last_coupling_value)
                
                
                if coupling_change >= self.coupling_change_threshold:
                    # Large coupling change detected!
                    coupling_shock_detected = True
                    force_immediate_update = True
                    self.coupling_shock_detected_time = current_elapsed_time_ps
                    self.coupling_recovery_active = True
                    
                    print(f"COUPLING SHOCK DETECTED at timestep {timestep}: t = {current_elapsed_time_ps:.6f} ps", flush=True)
                    print(f"  Coupling change: {self.last_coupling_value:.6e} → {current_coupling:.6e} a.u. (Δ = {coupling_change:.6e})", flush=True)
                    print(f"  Triggering shock dampening (threshold: ±{self.coupling_change_threshold:.2e} a.u.)", flush=True)
                
                # Update last coupling value for next iteration
                self.last_coupling_value = current_coupling
        
        # Legacy switch time detection (still supported)
        if (self.switch_time_ps is not None and not self.switch_detected and
            self.last_elapsed_time_ps < self.switch_time_ps <= current_elapsed_time_ps):
            
            # We've just crossed the switch time!
            self.switch_detected = True
            force_immediate_update = True
            print(f"SWITCH DETECTED at timestep {timestep}: t = {current_elapsed_time_ps:.6f} ps", flush=True)
            print(f"  Forcing immediate timestep adjustment for shock dampening", flush=True)
        
        # Update last elapsed time for next iteration
        self.last_elapsed_time_ps = current_elapsed_time_ps
        
        # Update error tolerance based on shock dampening and exponential ramping
        if self.adaptiveerror:
            # NEW: Dynamic coupling shock recovery
            if self.dynamic_coupling_detection:
                if coupling_shock_detected:
                    # Just detected a coupling shock - immediately apply shock tolerance
                    self.current_error_tolerance = self.shock_error_tolerance
                elif self.coupling_recovery_active:
                    # In recovery mode - exponential recovery from shock tolerance
                    recovery_time = current_elapsed_time_ps - self.coupling_shock_detected_time
                    if recovery_time >= 0:
                        exp_factor = np.exp(-recovery_time / self.time_constant_ps)
                        self.current_error_tolerance = self.target_error_tolerance - \
                                                      (self.target_error_tolerance - self.shock_error_tolerance) * exp_factor
                        
                        # Stop recovery when we're close to target tolerance
                        if exp_factor < 0.01:  # 99% recovered
                            self.coupling_recovery_active = False
                    else:
                        self.current_error_tolerance = self.shock_error_tolerance
                else:
                    # No shock detected and not in recovery - use target tolerance
                    self.current_error_tolerance = self.target_error_tolerance
            
            # Legacy switch time logic (still supported)
            elif self.switch_time_ps is not None:
                if current_elapsed_time_ps < self.switch_time_ps:
                    # Before switch: use target error tolerance
                    self.current_error_tolerance = self.target_error_tolerance
                else:
                    # After switch: implement shock dampening and exponential recovery
                    ramping_time = current_elapsed_time_ps - self.switch_time_ps
                    
                    if self.shock_dampening_enabled:
                        if force_immediate_update:
                            # At the exact moment of switch detection, use shock tolerance
                            self.current_error_tolerance = self.shock_error_tolerance
                        else:
                            # Normal exponential recovery after switch
                            exp_factor = np.exp(-ramping_time / self.time_constant_ps)
                            self.current_error_tolerance = self.target_error_tolerance - \
                                                            (self.target_error_tolerance - self.shock_error_tolerance) * exp_factor
                        #else:
                        #    # No ramping when initial_fraction = 0.0 - just use target tolerance
                        #    self.current_error_tolerance = self.target_error_tolerance
                    #if timestep % 1000 == 0:  # Debug output
                    #    print(f"DEBUG: shock_dampening mode, ramping_time={ramping_time:.6f}, current_error_tolerance={self.current_error_tolerance:.2e}", flush=True)
                    #elif self.shock_dampening_factor > 0:
                        # Original mode - exponential recovery from initial tolerance (only if initial_fraction > 0)
                    #    exp_factor = np.exp(-ramping_time / self.time_constant_ps)
                    #    self.current_error_tolerance = self.target_error_tolerance - \
                    #                                  (self.target_error_tolerance - self.initial_error_tolerance) * exp_factor
                        #if timestep % 1000 == 0:  # Debug output
                        #    print(f"DEBUG: ramping mode, ramping_time={ramping_time:.6f}, exp_factor={exp_factor:.6f}, current_error_tolerance={self.current_error_tolerance:.2e}", flush=True)
                        #if timestep % 1000 == 0:  # Debug output
                        #    print(f"DEBUG: no ramping mode, current_error_tolerance={self.current_error_tolerance:.2e}", flush=True)
            else:
                # Original behavior: start ramping immediately from t=0 (only if initial_fraction > 0)
                if self.shock_dampening_factor > 0:
                    exp_factor = np.exp(-current_elapsed_time_ps / self.time_constant_ps)
                    self.current_error_tolerance = self.target_error_tolerance - \
                                                  (self.target_error_tolerance - self.initial_error_tolerance) * exp_factor
                else:
                    # No ramping when initial_fraction = 0.0 - just use target tolerance
                    self.current_error_tolerance = self.target_error_tolerance
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
        
        # Sum forces from all force objects - let errors propagate clearly
        for force in self.integrator.forces:
            particle_forces = force.forces
            if particle_forces is not None and len(particle_forces) == n_particles:
                total_forces += np.asarray(particle_forces)
        
        # Calculate sum |f_i| / m_i
        force_norm = np.array([np.linalg.norm(f) for f in total_forces])
        force_max_norm = np.max(force_norm / masses) * n_particles
        
        # DEBUG: Show force information when timestep becomes very large
        # if timestep % 1000 == 0 or force_immediate_update:
        #     max_force = np.max(force_norm) if len(force_norm) > 0 else 0.0
        #     mean_force = np.mean(force_norm) if len(force_norm) > 0 else 0.0
        #     min_mass = np.min(masses) if len(masses) > 0 else 0.0
        #     max_mass = np.max(masses) if len(masses) > 0 else 0.0
            #print(f"DEBUG FORCES: max_force={max_force:.2e} a.u., mean_force={mean_force:.2e} a.u., force_mass_sum={force_mass_sum:.2e} a.u.^-1", flush=True)
            #print(f"DEBUG MASSES: min_mass={min_mass:.2e} a.u., max_mass={max_mass:.2e} a.u., n_particles={len(masses)}", flush=True)
            
        # Compute optimal timestep using current error tolerance
        if force_max_norm > 0:
            # Add minimum timestep protection to prevent numerical instability
            # Even during shock dampening, we need a non-zero timestep for the integrator
            min_dt = PhysicalConstants.fs_to_atomic_units(0.001)  # 0.001 fs minimum
            effective_error_tolerance = max(self.current_error_tolerance, 
                                          self.target_error_tolerance * 1e-10)  # Never below 1e-10 * target
            #print(f"DEBUG: effective_error_tolerance = {effective_error_tolerance}", flush=True)

            optimal_dt = np.sqrt(effective_error_tolerance / force_max_norm)
            current_dt = self.integrator.dt
            
            # Apply minimum timestep protection
            if optimal_dt < min_dt:
                optimal_dt = min_dt
                if force_immediate_update:
                    print(f"  WARNING: Computed timestep {np.sqrt(self.current_error_tolerance / force_max_norm):.2e} a.u. below minimum", flush=True)
                    print(f"  Clamping to minimum timestep: {min_dt:.2e} a.u. ({min_dt * PhysicalConstants.TIME_PS_CONVERSION * 1000:.3f} fs)", flush=True)
            
            # Apply maximum timestep protection  
            max_dt = PhysicalConstants.fs_to_atomic_units(10.0)  # 10 fs maximum
            if optimal_dt > max_dt:
                optimal_dt = max_dt
                if force_immediate_update:
                    print(f"  WARNING: Computed timestep above maximum, clamping to {max_dt * PhysicalConstants.TIME_PS_CONVERSION * 1000:.1f} fs", flush=True)
            
            # Apply the timestep change
            self.integrator.dt = optimal_dt
            self.last_timestep_update = timestep
                
        else:
            if timestep % 5000 == 0:
                print(f"WARNING: Zero force detected at step {timestep} - keeping current timestep", flush=True)

    @property
    def error_tolerance(self):
        """Access the current effective error tolerance."""
        return self.current_error_tolerance

    @hoomd.logging.log
    def elapsed_time_ps(self):
        """Return the elapsed time in picoseconds."""
        return self.accumulated_time_ps

