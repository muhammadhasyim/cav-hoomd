#!/usr/bin/env python3
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""
Advanced Cavity Molecular Dynamics Experiment Runner

This script provides a comprehensive framework for running cavity MD simulations
with various thermostat combinations and advanced analysis features.

FEATURES:
- Multiple thermostat combinations (Bussi, Langevin)
- Advanced analysis trackers (energy, F(k,t), cavity modes)
- Adaptive or fixed timestep control
- GPU and CPU support
- Comprehensive logging and output control
- SLURM array job support and local multi-replica execution

BASIC USAGE:
   # Run a single experiment with cavity coupling
   python 05_advanced_run.py --molecular-bath bussi --cavity-bath langevin --coupling 1e-3 --runtime 1000
   
   # Run without cavity (molecular-only simulation)
   python 05_advanced_run.py --molecular-bath bussi --no-cavity --runtime 1000
   
   # Run multiple replicas locally
   python 05_advanced_run.py --molecular-bath bussi --cavity-bath langevin --coupling 1e-3 --replicas 1-5 --runtime 500
   
   # Run with finite-q cavity mode
   python 05_advanced_run.py --molecular-bath bussi --cavity-bath langevin --finite-q --coupling 1e-3 --runtime 1000

THERMOSTAT OPTIONS:
- --molecular-bath: bussi, langevin, none (thermostat for molecular particles)
- --cavity-bath: bussi, langevin, none (thermostat for cavity particle)
- --finite-q: Enable finite-q cavity mode (default: q=0 mode)

REPLICA EXECUTION:
- Local execution: --replicas 1-5 or --replicas 1,2,3,4,5
- SLURM array jobs: Automatically detects SLURM_ARRAY_TASK_ID

ADVANCED FEATURES:
- --enable-energy-tracker: Detailed energy component tracking
- --enable-fkt: F(k,t) density correlation functions
- --fixed-timestep: Use fixed timestep instead of adaptive
- --device GPU: Run on GPU instead of CPU

OUTPUT CONTROL:
- Separate output periods for different observables (energy, F(k,t), trajectories, console)
- Organized directory structure
- Comprehensive logging with timestamps

See --help for all available options.
"""

import hoomd
import numpy as np
from numba import njit
from hoomd.bussi_reservoir.thermostats import BussiReservoir as Bussi
import logging
import sys
import os
import argparse
import time
from pathlib import Path
import datetime

# Import the CavityForce from the installed plugin
from hoomd.cavitymd import CavityForce

# Import analysis and tracking classes from the plugin
from hoomd.cavitymd import (
    CavityForce, PhysicalConstants, Status, ElapsedTimeTracker,
    TimestepFormatter, AdaptiveTimestepUpdater, FieldAutocorrelationTracker,
    CavityModeTracker, EnergyTracker
)

def unwrap_positions(positions, images, box_lengths):
    """Unwrap particle positions across periodic boundaries."""
    pos = np.asarray(positions)
    img = np.asarray(images)
    box = np.asarray(box_lengths)
    return pos + img * box[None, :]

# =============================================================================
# CUSTOM PERFORMANCE TRACKER FOR CONSOLE OUTPUT
# =============================================================================

class PerformanceTracker(hoomd.custom.Action):
    """Custom performance tracker to display ns/day and other metrics."""
    
    def __init__(self, simulation, runtime_ps, time_tracker=None):
        super().__init__()
        self.sim = simulation
        self.runtime_ps = runtime_ps
        self.time_tracker = time_tracker
        self.start_time = time.time()
        self.last_timestep = 0
        self.current_ns_per_day = 0.0
        self.current_eta = ""
        
    def act(self, timestep):
        if timestep <= 1:
            return
            
        # Calculate simulation progress
        if self.time_tracker is not None:
            simulation_time_ps = self.time_tracker.elapsed_time
        else:
            dt = float(self.sim.operations.integrator.dt)
            simulation_time_ps = PhysicalConstants.atomic_units_to_ps(dt * timestep)
        
        # Calculate wall time elapsed
        wall_time_elapsed = time.time() - self.start_time
        
        if wall_time_elapsed > 0:
            # Calculate ns/day
            ps_per_second = simulation_time_ps / wall_time_elapsed
            ns_per_second = ps_per_second / 1000.0
            self.current_ns_per_day = ns_per_second * 86400
            
            # Calculate ETA
            if simulation_time_ps > 0:
                total_wall_time_needed = (self.runtime_ps / simulation_time_ps) * wall_time_elapsed
                seconds_remaining = total_wall_time_needed - wall_time_elapsed
                if seconds_remaining > 0:
                    eta_td = datetime.timedelta(seconds=int(seconds_remaining))
                    self.current_eta = str(eta_td)
                else:
                    self.current_eta = "00:00:00"
            else:
                self.current_eta = "calculating..."
        
    @hoomd.logging.log
    def ns_per_day(self):
        return f"{self.current_ns_per_day:.2f}"
    
    @hoomd.logging.log
    def eta_remaining(self):
        return self.current_eta

# =============================================================================
# CAVITYMD SIMULATION CLASS (The only class we define locally)
# =============================================================================

class CavityMDSimulation:
    """
    A class to encapsulate cavity MD simulation setup and execution.
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
                 enable_text_output=False, text_output_file=None, truncate_gsd=False):
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

    def setup_simulation(self):
        """Create HOOMD simulation object and initialize state from GSD file."""
        import os
        import gsd.hoomd
        
        # Save current directory and change to job directory
        self.original_cwd = os.getcwd()
        os.chdir(self.job_dir)
        
        # Setup device
        device = self.setup_device()
        
        # Create simulation object
        self.sim = hoomd.Simulation(device=device, seed=np.random.randint(int(10**4)))
        
        # Load GSD file and handle cavity particle
        with gsd.hoomd.open(self.input_gsd, 'r') as f:
            if self.frame < 0:
                self.frame = len(f) + self.frame  # Convert negative index to positive
                if self.frame < 0:  # Handle case where abs(frame) > len(f)
                    self.frame = 0
            snapshot = f[self.frame]
            
            if self.incavity and self.add_cavity_particle:
                # Add new cavity particle
                self.log_info("Adding cavity particle to system...")
                snapshot = self.create_cavity_particle(snapshot)
                self.sim.create_state_from_snapshot(snapshot)
                self.log_info(f"Simulation state created from modified snapshot with cavity particle")
            else:
                # Use original GSD file
                self.sim.create_state_from_gsd(filename=self.input_gsd, frame=self.frame)
                self.log_info(f"Simulation state created from original GSD file frame {self.frame}")
                
                # Validate cavity particle if needed
                if self.incavity:
                    self.validate_cavity_particle()
        
        return snapshot

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

    def create_cavity_particle(self, snapshot):
        """Add a cavity particle to the simulation snapshot."""
        self.log_info("Adding cavity particle to system...")
        
        # Calculate dipole moment and photon position
        positions = unwrap_positions(snapshot.particles.position, snapshot.particles.image, 
                                   snapshot.configuration.box[:3])
        dipmom = np.einsum('i,ij->j', snapshot.particles.charge, positions)

        omegac = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1
        
        if self.finite_q:
            # Allow finite-q photon displacement based on dipole-coupling interaction
            newpos = -dipmom * self.couplstr / omegac**2
            newpos[-1] = 0.0
            # Only add thermal fluctuations if coupling is non-zero
            if self.couplstr != 0.0:
                sigma = np.sqrt(self.kB * self.temperature / omegac**2)
                newpos = np.random.normal(loc=newpos, scale=sigma, size=3)
                self.log_info(f"Finite-q mode: Photon displaced by dipole interaction to {newpos} (with thermal fluctuations)")
            else:
                self.log_info(f"Finite-q mode: Photon at equilibrium position {newpos} (no thermal fluctuations due to zero coupling)")
        else:
            # Start photon at origin (q=0 limit)
            newpos = np.array([0.0, 0.0, 0.0])
            # Only add thermal fluctuations if coupling is non-zero
            if self.couplstr != 0.0:
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
        with self.sim.state.cpu_local_snapshot as snap:
            if 'L' not in snap.particles.types:
                raise ValueError("ERROR: Cavity simulation requested but no cavity particle type 'L' found in GSD file.")
            
            if 2 not in snap.particles.typeid:
                raise ValueError("ERROR: Cavity simulation requested but no cavity particles found in GSD file.")
            
            cavity_count = np.sum(snap.particles.typeid == 2)
            if cavity_count != 1:
                raise ValueError(f"ERROR: Expected exactly 1 cavity particle but found {cavity_count} in GSD file.")
            
            cavity_index = np.where(snap.particles.typeid == 2)[0][0]
            cavity_position = snap.particles.position[cavity_index]
            self.log_info(f"Cavity particle validated at index {cavity_index}, position {cavity_position}")

    def setup_force_parameters(self, dt, rcut=15):
        """Set up force parameters for the simulation."""
        forces = []
        
        # Setup cavity force if requested
        if self.incavity:
            omegac = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1
            cavityforce = CavityForce(kvector=np.array([0,0,1]), couplstr=self.couplstr, omegac=omegac)
            forces.append(cavityforce)
        
        # Setup harmonic bonds
        harmonic = hoomd.md.bond.Harmonic()
        harmonic.params['O-O'] = dict(k=2*0.36602, r0=2.281655158)
        harmonic.params['N-N'] = dict(k=2*0.71625, r0=2.0743522177)
        forces.append(harmonic)
        
        # Setup neighbor list
        cell = hoomd.md.nlist.Cell(buffer=1.0, exclusions=('bond',))
        
        # Setup Lennard-Jones interactions
        lj = hoomd.md.pair.LJ(nlist=cell, mode='shift')
        lj.params[('O', 'O')] = dict(epsilon=0.00016685201, sigma=6.230426584)
        lj.r_cut[('O', 'O')] = rcut
        lj.params[('N', 'N')] = dict(epsilon=0.000083426, sigma=5.48277488)
        lj.r_cut[('N', 'N')] = rcut
        lj.params[('N', 'O')] = dict(epsilon=0.00025027802, sigma=4.9832074319)
        lj.r_cut[('N', 'O')] = rcut

        # Disable pair interaction with 'L' particle (photon)
        if self.incavity:
            lj.params[('L', 'N')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('L', 'N')] = 0.0
            lj.params[('N', 'L')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('N', 'L')] = 0.0
            lj.params[('O', 'L')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('O', 'L')] = 0.0
            lj.params[('L', 'O')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('L', 'O')] = 0.0
            lj.params[('L', 'L')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('L', 'L')] = 0.0
        forces.append(lj)

        # Setup long-range Coulomb interactions using PPPM method
        numpoints = 32
        order = 6
        short, long = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
            nlist=cell, resolution=[numpoints, numpoints, numpoints], 
            order=order, r_cut=rcut, alpha=0.0
        )
        forces.append(short)
        forces.append(long)
        
        return forces

    def setup_thermostat_parameters(self, dt):
        """Set up thermostat parameters for molecular and cavity systems."""
        kT = self.kB * self.temperature
        molecular_filter = hoomd.filter.Type(['O', 'N'])  # Molecular particles only
        
        # Convert thermostat time constants from ps to atomic units
        molecular_tau_au = PhysicalConstants.ps_to_atomic_units(self.molecular_thermostat_tau)
        cavity_tau_au = PhysicalConstants.ps_to_atomic_units(self.cavity_thermostat_tau)
        
        self.log_info(f"Thermostat time constant conversions:")
        self.log_info(f"  Molecular: {self.molecular_thermostat_tau:.3f} ps = {molecular_tau_au:.6f} a.u.")
        self.log_info(f"  Cavity: {self.cavity_thermostat_tau:.3f} ps = {cavity_tau_au:.6f} a.u.")
        
        # Validate tau=0.0 with Langevin thermostats
        if self.molecular_thermostat.lower() == 'langevin' and self.molecular_thermostat_tau <= 0.0:
            raise ValueError(
                f"ERROR: Cannot use Langevin thermostat with molecular_thermostat_tau={self.molecular_thermostat_tau} ps.\n"
                f"Langevin dynamics requires tau > 0 since gamma = 1/tau.\n"
                f"For overdamped dynamics (tau → 0), use Brownian dynamics instead."
            )
        
        if self.incavity and self.cavity_thermostat.lower() == 'langevin' and self.cavity_thermostat_tau <= 0.0:
            raise ValueError(
                f"ERROR: Cannot use Langevin thermostat with cavity_thermostat_tau={self.cavity_thermostat_tau} ps.\n"
                f"Langevin dynamics requires tau > 0 since gamma = 1/tau.\n"
                f"For overdamped dynamics (tau → 0), use Brownian dynamics instead."
            )
        
        # Store references to thermostats for energy tracking
        thermostat_refs = {
            'molecular_langevin': None,
            'cavity_langevin': None,
            'molecular_mttk': None,
            'cavity_mttk': None,
            'molecular_bussi': None,
            'cavity_bussi': None
        }
        
        # Configure molecular thermostat
        if self.molecular_thermostat.lower() == 'bussi':
            self.log_info("Running molecular system with Bussi thermostat (NVT ensemble)")
            molecular_bussi = Bussi(kT=kT, tau=molecular_tau_au)
            molecular_method = hoomd.md.methods.ConstantVolume(filter=molecular_filter, thermostat=molecular_bussi)
            thermostat_refs['molecular_bussi'] = molecular_bussi
            self.log_info(f"Molecular Bussi thermostat configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u., tau = {self.molecular_thermostat_tau:.3f} ps")
        elif self.molecular_thermostat.lower() == 'langevin':
            self.log_info("Running molecular system with Langevin thermostat (NVT ensemble)")
            molecular_gamma = PhysicalConstants.gamma_from_tau_ps(self.molecular_thermostat_tau)
            molecular_method = hoomd.md.methods.Langevin(filter=molecular_filter, kT=kT, default_gamma=molecular_gamma, tally_reservoir_energy=True)
            thermostat_refs['molecular_langevin'] = molecular_method
            self.log_info(f"Molecular Langevin thermostat configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u.")
            self.log_info(f"  gamma = {molecular_gamma:.6f} a.u.^-1 (tau = {self.molecular_thermostat_tau:.3f} ps)")
        elif self.molecular_thermostat.lower() == 'none':
            self.log_info("Running molecular system without thermostat (NVE ensemble)")
            molecular_method = hoomd.md.methods.ConstantVolume(filter=molecular_filter)
        else:
            raise ValueError(f"Invalid molecular_thermostat option: {self.molecular_thermostat}")
        
        # Set up thermostat for cavity particle if present
        cavity_method = None
        if self.incavity:
            cavity_filter = hoomd.filter.Type(['L'])  # Cavity particle only
            
            if self.cavity_thermostat.lower() == 'langevin':
                self.log_info("Running cavity with Langevin thermostat")
                base_gamma = PhysicalConstants.gamma_from_tau_ps(self.cavity_thermostat_tau)
                cavity_gamma = self.cavity_damping_factor * base_gamma
                cavity_method = hoomd.md.methods.Langevin(filter=cavity_filter, 
                                                         kT=kT, default_gamma=cavity_gamma, tally_reservoir_energy=True)
                thermostat_refs['cavity_langevin'] = cavity_method
                self.log_info(f"Cavity Langevin thermostat configured: T = {self.temperature:.1f} K, kT = {kT:.6f} a.u.")
                self.log_info(f"  base_gamma = {base_gamma:.6f} a.u.^-1 (tau = {self.cavity_thermostat_tau:.3f} ps)")
                self.log_info(f"  effective_gamma = {cavity_gamma:.6f} a.u.^-1 (damping_factor = {self.cavity_damping_factor:.1f}x)")
            elif self.cavity_thermostat.lower() == 'bussi':
                self.log_info("Running cavity with Bussi thermostat")
                cavity_bussi = Bussi(kT=kT, tau=cavity_tau_au)
                cavity_method = hoomd.md.methods.ConstantVolume(filter=cavity_filter, thermostat=cavity_bussi)
                thermostat_refs['cavity_bussi'] = cavity_bussi
                self.log_info(f"Cavity Bussi thermostat configured: kT = {kT:.6f} a.u., tau = {self.cavity_thermostat_tau:.3f} ps")
            elif self.cavity_thermostat.lower() == 'none':
                self.log_info("Running cavity without thermostat (NVE ensemble)")
                cavity_method = hoomd.md.methods.ConstantVolume(filter=cavity_filter)
            else:
                raise ValueError(f"Invalid cavity_thermostat option: {self.cavity_thermostat}")
        
        return molecular_method, cavity_method, thermostat_refs

    def setup_integrator(self, forces, methods):
        """Configure the integrator with forces and integration methods."""
        # Setup integrator with initial dt
        integrator = hoomd.md.Integrator(dt=self.dt, forces=forces)
        self.sim.operations.integrator = integrator
        
        # Set integration methods (filter out None methods)
        valid_methods = [method for method in methods if method is not None]
        self.sim.operations.integrator.methods = valid_methods
        
        self.log_info(f"Integrator configured with initial dt = {self.dt:.6f} a.u. ({self.dt_ps:.6f} ps)")
        self.log_info(f"Number of integration methods: {len(valid_methods)}")

    def thermalize_system(self):
        """Initialize particle velocities and thermostat degrees of freedom."""
        kT = self.kB * self.temperature
        
        # Thermalize particle momenta
        if self.incavity:
            # Only thermalize molecular particles, not cavity particle
            molecular_filter = hoomd.filter.Type(['O', 'N'])
            self.sim.state.thermalize_particle_momenta(kT=kT, filter=molecular_filter)
            self.log_info("Thermalized molecular particles only (cavity particle excluded)")
            
            # Initialize cavity particle velocity based on thermostat type
            with self.sim.state.cpu_local_snapshot as snap:
                cavity_indices = np.where(snap.particles.typeid == 2)[0]
                if len(cavity_indices) > 0:
                    cavity_idx = cavity_indices[0]  # Get first cavity particle
                    
                    # For 3D Maxwell-Boltzmann: each component has variance kT/m
                    # With mass = 1.0 a.u., std dev per component = sqrt(kT)
                    cavity_velocity = np.random.normal(0.0, np.sqrt(kT), size=3)
                    
                    # Calculate expected kinetic energy and temperature
                    expected_ke = 0.5 * 1.0 * np.sum(cavity_velocity**2)  # KE = (1/2) * m * v²
                    expected_temp = (2.0/3.0) * expected_ke / self.kB  # T = (2/3) * KE / kB for 3D

                    self.log_info(f"Cavity particle thermalization:")
                    self.log_info(f"  Target temperature: {self.temperature:.1f} K")
                    self.log_info(f"  kT = {kT:.6f} a.u.")
                    self.log_info(f"  Initial velocity: {cavity_velocity}")
                    self.log_info(f"  Initial KE: {expected_ke:.6f} a.u.")
                    self.log_info(f"  Expected temperature: {expected_temp:.1f} K")
                    self.log_info(f"  Thermostat: {self.cavity_thermostat}")
                    
                    snap.particles.velocity[cavity_idx] = cavity_velocity
                    
                else:
                    self.log_info("WARNING: No cavity particle found for thermalization!")
        else:
            # No cavity particle, thermalize all particles
            self.sim.state.thermalize_particle_momenta(kT=kT, filter=hoomd.filter.All())
            self.log_info("Thermalized all molecular particles")
        
        # Initialize reservoir energy logging (requires at least one simulation step)
        self.log_info("Initializing reservoir energy tracking...")
        self.sim.run(1)

    def compute_and_set_optimal_timestep(self):
        """Compute and set the optimal timestep after running one step to initialize forces."""
        if self.error_tolerance <= 0:
            # Fixed timestep mode
            if self.dt_fs is not None:
                # Use user-specified timestep
                dt_au = PhysicalConstants.ps_to_atomic_units(self.dt_fs / 1000.0)  # Convert fs to ps, then to a.u.
                self.sim.operations.integrator.dt = dt_au
                self.dt = dt_au
                self.log_info(f"Using user-specified fixed timestep: {dt_au:.6f} a.u. ({self.dt_fs:.3f} fs)")
            else:
                # Keep current dt (HOOMD default)
                self.log_info(f"Using default fixed timestep: {self.sim.operations.integrator.dt:.6f} a.u. ({PhysicalConstants.atomic_units_to_ps(self.sim.operations.integrator.dt) * 1000:.3f} fs)")
            return
        
        try:
            self.log_info("Computing optimal timestep...")
            
            # Run one step to initialize forces (required by HOOMD)
            self.sim.run(1)
            
            # Use initial error tolerance
            initial_error_tolerance = self.error_tolerance * 1e-3  # initial_fraction = 1e-3
            
            # Collect forces and masses
            particle_data = self.sim.state.get_snapshot().particles
            masses = np.array(particle_data.mass)
            n_particles = len(masses)
            
            # Initialize total forces array
            total_forces = np.zeros((n_particles, 3))
            
            # Sum forces from all force objects
            for force in self.sim.operations.integrator.forces:
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
                optimal_dt = np.sqrt(initial_error_tolerance / force_mass_sum)
                
                # Update integrator timestep
                self.sim.operations.integrator.dt = optimal_dt
                self.dt = optimal_dt  # Update stored dt
                
                self.log_info(f"Optimal timestep computed and set:")
                self.log_info(f"  Initial error tolerance: {initial_error_tolerance:.2e}")
                self.log_info(f"  Force/mass sum: {force_mass_sum:.6e}")
                self.log_info(f"  Optimal dt: {optimal_dt:.6f} a.u. ({PhysicalConstants.atomic_units_to_ps(optimal_dt) * 1000:.3f} fs)")
            else:
                self.log_info("WARNING: Force/mass sum is zero - keeping initial timestep")
                
        except Exception as e:
            self.log_info(f"WARNING: Failed to compute optimal timestep: {str(e)}")
            self.log_info(f"Using initial timestep: {self.sim.operations.integrator.dt:.6f} a.u.")

    def setup_trackers_and_loggers(self):
        """Set up comprehensive tracking and logging objects for the simulation."""
        # Create elapsed time tracker
        self.time_tracker = ElapsedTimeTracker(self.sim, self.runtime_ps)
        self.sim.operations.updaters.append(hoomd.update.CustomUpdater(
            action=self.time_tracker, trigger=hoomd.trigger.Periodic(1)
        ))
        
        # Create custom performance tracker
        self.performance_tracker = PerformanceTracker(self.sim, self.runtime_ps, self.time_tracker)
        self.sim.operations.updaters.append(hoomd.update.CustomUpdater(
            action=self.performance_tracker, trigger=hoomd.trigger.Periodic(100)
        ))
        
        # Set up adaptive timestep updater if error_tolerance is positive
        if self.error_tolerance > 0:
            self.log_info(f"Setting up adaptive timestep updater (error_tolerance = {self.error_tolerance})")
            self.adaptive_action = AdaptiveTimestepUpdater(
                state=self.sim.state,
                integrator=self.sim.operations.integrator,
                error_tolerance=self.error_tolerance,
                time_constant_ps=50.0,
                initial_fraction=1e-3,
                adaptiveerror=True,
                cavity_damping_factor=self.cavity_damping_factor,
                molecular_thermostat_tau=self.molecular_thermostat_tau,
                cavity_thermostat_tau=self.cavity_thermostat_tau,
                time_tracker=self.time_tracker
            )
            
            # Add adaptive updater - use energy period for adaptive timestep updates
            adaptive_updater = hoomd.update.CustomUpdater(
                action=self.adaptive_action,
                trigger=hoomd.trigger.Periodic(self.energy_period)
            )
            self.sim.operations.updaters.append(adaptive_updater)
            self.log_info("Adaptive timestep updater enabled")
        else:
            self.adaptive_action = None
            self.log_info("Fixed timestep mode (error_tolerance = 0)")
        
        # Create status tracker for performance metrics (kept for compatibility)
        self.status = Status(self.sim, self.runtime_ps, self.time_tracker)
        
        # Create timestep formatter for fs display
        self.timestep_formatter = TimestepFormatter(self.sim.operations.integrator)
        
        # Create comprehensive logger
        logger = hoomd.logging.Logger(categories=['scalar', 'string'])
        
        # Basic simulation quantities
        logger.add(self.sim, quantities=['timestep', 'tps'])
        
        # Time and performance information
        logger[('Time', 'elapsed_ps')] = (self.time_tracker, 'elapsed_time', 'scalar')
        logger[('Performance', 'ns_per_day')] = (self.performance_tracker, 'ns_per_day', 'string')
        logger[('Performance', 'eta')] = (self.performance_tracker, 'eta_remaining', 'string')
        logger[('Timestep', 'dt_fs')] = (self.timestep_formatter, 'dt_fs', 'scalar')
        
        # Add adaptive timestep logging if enabled
        if hasattr(self, 'adaptive_action') and self.adaptive_action is not None:
            logger[('Adaptive', 'error_tolerance')] = (self.adaptive_action, 'error_tolerance', 'scalar')
            logger[('Adaptive', 'elapsed_time_ps')] = (self.adaptive_action, 'elapsed_time_ps', 'scalar')
        
        # Add thermodynamic quantities
        try:
            # Add integrator energies - DISABLED: integrator doesn't expose these quantities
            # logger.add(self.sim.operations.integrator, quantities=['kinetic_energy', 'potential_energy'])
            
            # Create molecular thermodynamics computer
            molecular_filter = hoomd.filter.Type(['O', 'N'])
            self.molecular_thermo = hoomd.md.compute.ThermodynamicQuantities(filter=molecular_filter)
            self.sim.operations.computes.append(self.molecular_thermo)
            logger[('Molecular', 'temperature')] = (self.molecular_thermo, 'kinetic_temperature', 'scalar')
            logger[('Molecular', 'kinetic_energy')] = (self.molecular_thermo, 'kinetic_energy', 'scalar')
            
            if self.incavity:
                # Add cavity particle thermodynamics
                cavity_filter = hoomd.filter.Type(['L'])
                self.cavity_thermo = hoomd.md.compute.ThermodynamicQuantities(filter=cavity_filter)
                self.sim.operations.computes.append(self.cavity_thermo)
                logger[('Cavity', 'temperature')] = (self.cavity_thermo, 'kinetic_temperature', 'scalar')
                logger[('Cavity', 'kinetic_energy')] = (self.cavity_thermo, 'kinetic_energy', 'scalar')
            
        except Exception as e:
            self.log_warning(f"Could not add some thermodynamic quantities: {str(e)}")
        
        # Add thermostat energies if enabled and available
        if self.enable_energy_tracking:
            self.log_info("Setting up detailed energy tracking")
            try:
                # Track reservoir energies from integration methods
                method_count = 0
                for i, method in enumerate(self.sim.operations.integrator.methods):
                    if hasattr(method, 'reservoir_energy'):
                        try:
                            logger[('Method', f'reservoir_energy_{i}')] = (method, 'reservoir_energy', 'scalar')
                            method_count += 1
                            self.log_info(f"  Added reservoir energy tracking for method {i}")
                        except Exception as e:
                            self.log_warning(f"  Failed to add reservoir energy for method {i}: {e}")
                    
                    # Track thermostat energies
                    if hasattr(method, 'thermostat') and method.thermostat is not None:
                        if hasattr(method.thermostat, 'energy'):
                            try:
                                logger[('Thermostat', f'energy_{i}')] = (method.thermostat, 'energy', 'scalar')
                                self.log_info(f"  Added thermostat energy tracking for method {i}")
                            except Exception as e:
                                self.log_warning(f"  Failed to add thermostat energy for method {i}: {e}")
                        
                        # Track detailed BussiReservoir energies if available
                        if hasattr(method.thermostat, 'total_reservoir_energy'):
                            try:
                                logger[('Thermostat', f'total_reservoir_energy_{i}')] = (method.thermostat, 'total_reservoir_energy', 'scalar')
                                logger[('Thermostat', f'translational_energy_{i}')] = (method.thermostat, 'reservoir_energy_translational', 'scalar')
                                logger[('Thermostat', f'rotational_energy_{i}')] = (method.thermostat, 'reservoir_energy_rotational', 'scalar')
                                self.log_info(f"  Added detailed BussiReservoir energy tracking for method {i}")
                            except Exception as e:
                                self.log_warning(f"  Failed to add detailed thermostat energies for method {i}: {e}")
                
                # Track individual force energies
                force_count = 0
                for i, force in enumerate(self.sim.operations.integrator.forces):
                    force_name = type(force).__name__
                    try:
                        if hasattr(force, 'energy'):
                            logger[('Force', f'{force_name}_energy')] = (force, 'energy', 'scalar')
                            force_count += 1
                            self.log_info(f"  Added energy tracking for {force_name}")
                        
                        # Special handling for CavityForce with component energies
                        if hasattr(force, 'harmonic_energy'):
                            logger[('Cavity', 'harmonic_energy')] = (force, 'harmonic_energy', 'scalar')
                            logger[('Cavity', 'coupling_energy')] = (force, 'coupling_energy', 'scalar')
                            logger[('Cavity', 'dipole_self_energy')] = (force, 'dipole_self_energy', 'scalar')
                            self.log_info(f"  Added detailed cavity energy components")
                    except Exception as e:
                        self.log_warning(f"  Failed to add energy tracking for {force_name}: {e}")
                
                self.log_info(f"Energy tracking setup completed: {method_count} methods, {force_count} forces")
                
                # Set up cavity mode tracking if in cavity simulation
                if self.incavity:
                    # Find the cavity force object to pass to CavityModeTracker
                    cavityforce = None
                    for force in self.sim.operations.integrator.forces:
                        if 'cavity' in type(force).__name__.lower():
                            cavityforce = force
                            break
                    
                    if cavityforce is not None:
                        self.cavity_mode_tracker = CavityModeTracker(
                            simulation=self.sim,
                            cavityforce=cavityforce,
                            time_tracker=self.time_tracker,
                            output_prefix=f'{self.name}-{self.replica}',
                            output_period_steps=1
                        )
                        
                        # Add cavity mode tracker as updater
                        cavity_mode_updater = hoomd.update.CustomUpdater(
                            action=self.cavity_mode_tracker,
                            trigger=hoomd.trigger.Periodic(1)
                        )
                        self.sim.operations.updaters.append(cavity_mode_updater)
                        self.log_info("CavityModeTracker created and added to simulation")
                    else:
                        self.cavity_mode_tracker = None
                        self.log_warning("Cavity simulation requested but no cavity force found - cavity mode tracker disabled")
                else:
                    self.cavity_mode_tracker = None
                    self.log_info("No cavity simulation - cavity mode tracker not created")

                # Set up energy tracker from plugin for proper reservoir energy tracking
                try:
                    # Fix duplicate "energy" in filename - EnergyTracker adds "_energy_tracker.txt" automatically
                    output_prefix = f'{self.name}-{self.replica}'
                    actual_energy_filename = f'{output_prefix}_energy_tracker.txt'
                    self.log_info(f"Setting up energy tracker with output file: {actual_energy_filename}")
                    
                    # FIXED: Calculate reasonable step-based values but with smaller periods for better time accuracy
                    # Use small trigger period (1 step) and let tracker handle its internal timing
                    output_trigger_period = 1
                    
                    # Calculate max_timesteps from max_time_ps if specified
                    max_timesteps = None
                    if self.max_energy_output_time_ps:
                        # Use current timestep for calculation (better than initial estimate)
                        current_dt = self.sim.operations.integrator.dt
                        dt_ps = PhysicalConstants.atomic_units_to_ps(current_dt)
                        max_timesteps = int(self.max_energy_output_time_ps / dt_ps)
                        self.log_info(f"  Max timesteps calculated: {max_timesteps} (≈{self.max_energy_output_time_ps:.3f} ps)")
                    
                    self.log_info(f"  Output trigger: every {output_trigger_period} step")
                    self.log_info(f"  Target max time: {self.max_energy_output_time_ps} ps" if self.max_energy_output_time_ps else "  No max time limit")
                    
                    # Prepare individual force objects for EnergyTracker
                    force_objects = {}
                    thermostat_objects = {}
                    
                    for force in self.sim.operations.integrator.forces:
                        force_name = type(force).__name__.lower()
                        if 'cavity' in force_name:
                            force_objects['cavity'] = force
                        elif 'lj' in force_name or 'lennard' in force_name:
                            force_objects['lj'] = force
                        elif 'harmonic' in force_name or 'bond' in force_name:
                            force_objects['harmonic'] = force
                        elif 'ewald' in force_name:
                            force_objects['ewald_short'] = force  # Ewald is the short-range PPPM force
                        elif 'coulomb' in force_name:
                            force_objects['ewald_long'] = force   # Coulomb is the long-range PPPM force
                    
                    # Extract thermostat objects for EnergyTracker
                    for i, method in enumerate(self.sim.operations.integrator.methods):
                        method_name = type(method).__name__.lower()
                        
                        # Check for Langevin methods (have reservoir_energy)
                        if 'langevin' in method_name and hasattr(method, 'reservoir_energy'):
                            # Determine if this is molecular or cavity based on filter
                            if hasattr(method, 'filter'):
                                filter_types = getattr(method.filter, '_types', [])
                                if 'L' in filter_types:
                                    thermostat_objects['langevin_cavity'] = method
                                else:
                                    thermostat_objects['langevin_molecular'] = method
                        
                        # Check for Bussi thermostats
                        if hasattr(method, 'thermostat') and method.thermostat is not None:
                            thermostat_type = type(method.thermostat).__name__.lower()
                            if 'bussi' in thermostat_type:
                                # Determine if this is molecular or cavity based on filter
                                if hasattr(method, 'filter'):
                                    filter_types = getattr(method.filter, '_types', [])
                                    if 'L' in filter_types:
                                        thermostat_objects['bussi_cavity'] = method.thermostat
                                    else:
                                        thermostat_objects['bussi_molecular'] = method.thermostat
                    
                    # Get kinetic trackers 
                    kinetic_tracker = None  # Use internal kinetic computation
                    cavity_mode_tracker = getattr(self, 'cavity_mode_tracker', None) if self.incavity else None
                    
                    self.log_info(f"Energy tracker configuration:")
                    self.log_info(f"  Force objects: {list(force_objects.keys())}")
                    self.log_info(f"  Thermostat objects: {list(thermostat_objects.keys())}")
                    self.log_info(f"  Using internal kinetic computation (no external tracker needed)")
                    self.log_info(f"  Cavity mode tracker available: {cavity_mode_tracker is not None}")
                    
                    # Use time-based limit directly instead of converted timesteps
                    self.energy_tracker = EnergyTracker(
                        simulation=self.sim,
                        components=['kinetic', 'harmonic', 'lj', 'ewald_short', 'ewald_long', 'cavity'],
                        force_objects=force_objects,
                        thermostat_objects=thermostat_objects,
                        kinetic_tracker=kinetic_tracker,  # Use internal kinetic computation
                        cavity_mode_tracker=cavity_mode_tracker,
                        time_tracker=self.time_tracker,
                        output_prefix=output_prefix,
                        output_period_steps=1,  # Check every step for best timing accuracy
                        max_time_ps=self.max_energy_output_time_ps,  # Use time-based limit directly
                        compute_temperature=True,
                        track_reservoirs=True,
                        verbose='quiet'  # Reduce verbose output - use 'verbose' for full debug output
                    )
                    
                    # Add energy tracker to simulation with small trigger period
                    energy_updater = hoomd.update.CustomUpdater(
                        action=self.energy_tracker,
                        trigger=hoomd.trigger.Periodic(output_trigger_period)
                    )
                    self.sim.operations.updaters.append(energy_updater)
                    
                    self.log_info(f"Energy tracker setup completed successfully")
                    self.log_info(f"  Uses internal timing logic with max_time_ps for accurate time limits")
                    self.log_info(f"  Checking every step for precise timing control")
                    if self.max_energy_output_time_ps:
                        self.log_info(f"  Output limited to {self.max_energy_output_time_ps:.3f} ps (time-based)")
                    self.log_info(f"  Trigger period: {output_trigger_period} step")
                    
                except Exception as e:
                    self.log_error(f"Failed to setup energy tracker: {e}")
                    import traceback
                    self.log_error("Full traceback:")
                    for line in traceback.format_exc().split('\n'):
                        if line.strip():
                            self.log_error(line)
                    self.energy_tracker = None
                
            except Exception as e:
                self.log_warning(f"Could not complete energy tracking setup: {str(e)}")
                self.log_warning(f"Error details: {type(e).__name__}: {str(e)}")
                self.log_warning(f"Full traceback:")
                import traceback
                traceback.print_exc()
        else:
            self.log_info("Detailed energy tracking disabled")
            self.energy_tracker = None
        
        # Set up F(k,t) density correlation tracker if enabled
        if self.enable_fkt:
            try:
                self.log_info("Setting up F(k,t) density correlation tracker")
                self.log_info(f"  k magnitude: {self.fkt_kmag:.2f}")
                self.log_info(f"  Number of wavevectors: {self.fkt_num_wavevectors}")
                self.log_info(f"  Reference interval: {self.fkt_reference_interval_ps:.3f} ps")
                self.log_info(f"  Max references: {self.fkt_max_references}")
                self.log_info(f"  Output period: {self.fkt_output_period_ps:.3f} ps")
                
                # Use time-based intervals for adaptive timestep compatibility
                fkt_trigger_period = 1  # Check every step for best timing
                
                self.log_info(f"  Using time-based reference intervals for adaptive timestep compatibility")
                self.log_info(f"  Reference interval: {self.fkt_reference_interval_ps:.3f} ps")
                self.log_info(f"  Trigger period: {fkt_trigger_period} step")
                
                # Create density correlation tracker with time-based reference interval
                self.density_corr_tracker = FieldAutocorrelationTracker(
                    simulation=self.sim,
                    observable="density_correlation",
                    time_tracker=self.time_tracker,
                    output_period_steps=1,  # Check every step for best timing accuracy
                    output_prefix=f'{self.name}-{self.replica}',
                    reference_interval_ps=self.fkt_reference_interval_ps,  # Use time-based interval
                    max_references=self.fkt_max_references,
                    kmag=self.fkt_kmag,
                    num_wavevectors=self.fkt_num_wavevectors
                )
                
                # Add F(k,t) tracker to simulation with small trigger period
                fkt_updater = hoomd.update.CustomUpdater(
                    action=self.density_corr_tracker,
                    trigger=hoomd.trigger.Periodic(fkt_trigger_period)
                )
                self.sim.operations.updaters.append(fkt_updater)
                
                # Add F(k,t) data to logger
                logger[('F(k,t)', 'current_autocorr')] = (self.density_corr_tracker, 'current_autocorr', 'scalar')
                
                self.log_info("F(k,t) tracker successfully enabled")
                self.log_info(f"  Uses time-based intervals for accurate timing in adaptive mode")
                self.log_info(f"  Reference interval automatically adjusts to changing timesteps")
                
            except Exception as e:
                self.log_warning(f"Could not set up F(k,t) tracker: {str(e)}")
                self.density_corr_tracker = None
        else:
            self.density_corr_tracker = None
            self.log_info("F(k,t) tracking disabled")
        
        # Store logger for later use
        self.logger_hoomd = logger
        
        # Create console output with only performance and time metrics
        console_items = ["timestep", "tps", "elapsed_time", "ns_per_day", "eta", "dt(fs)"]
        
        # Add adaptive timestep info (performance related)
        if hasattr(self, 'adaptive_action') and self.adaptive_action is not None:
            console_items.append("adaptive_error_tolerance")
        
        self.log_info("CORRECTED tracking and logging setup completed")
        self.log_info(f"Console output includes: {', '.join(console_items)}")
        
        # Log detailed summary of what's enabled
        enabled_features = []
        if self.enable_energy_tracking:
            enabled_features.append("detailed energy tracking")
        if self.enable_fkt:
            enabled_features.append(f"F(k,t) density correlation (k={self.fkt_kmag})")
        if hasattr(self, 'adaptive_action') and self.adaptive_action is not None:
            enabled_features.append("adaptive timestep control")
        
        if enabled_features:
            self.log_info(f"Advanced features enabled: {', '.join(enabled_features)}")
        else:
            self.log_info("Running with basic tracking only")

    def setup_output_writers(self):
        """Configure GSD writer and console table for simulation output."""
        
        # Choose appropriate trigger for adaptive vs fixed timestep mode
        if self.error_tolerance > 0:
            # Adaptive timestep mode - use very small step intervals for better timing accuracy
            # Since timestep changes frequently, use small intervals to minimize timing error
            gsd_trigger_steps = max(1, int(self.gsd_output_period_ps / 0.001))  # Assume ~1 fs effective timestep
            console_trigger_steps = max(1, int(self.console_output_period_ps / 0.001))  # Assume ~1 fs effective timestep
            
            # Cap the intervals to reasonable values to avoid excessive overhead
            gsd_trigger_steps = min(gsd_trigger_steps, 10000)  # Max 10k steps
            console_trigger_steps = min(console_trigger_steps, 1000)  # Max 1k steps
            
            gsd_trigger = hoomd.trigger.Periodic(gsd_trigger_steps)
            console_trigger = hoomd.trigger.Periodic(period=console_trigger_steps)
            
            self.log_info("Using small step-based triggers for better timing accuracy in adaptive mode")
            self.log_info(f"  GSD trigger: every {gsd_trigger_steps} steps (target: {self.gsd_output_period_ps:.3f} ps)")
            self.log_info(f"  Console trigger: every {console_trigger_steps} steps (target: {self.console_output_period_ps:.3f} ps)")
            self.log_info("  Note: Actual timing may vary slightly due to adaptive timestep changes")
            
        else:
            # Fixed timestep mode - use calculated step-based triggers (most efficient)
            gsd_trigger = hoomd.trigger.Periodic(self.gsd_period)
            console_trigger = hoomd.trigger.Periodic(period=self.console_period)
            
            self.log_info("Using calculated step-based triggers for fixed timestep mode")
            self.log_info(f"  GSD trigger: every {self.gsd_period} steps ({self.gsd_output_period_ps:.3f} ps)")
            self.log_info(f"  Console trigger: every {self.console_period} steps ({self.console_output_period_ps:.3f} ps)")
        
        # Set up GSD writer with appropriate trigger
        gsd_writer = hoomd.write.GSD(
            filename=f'{self.name}-{self.replica}.gsd',
            trigger=gsd_trigger,
            dynamic=['property', 'momentum', 'particles/diameter', 'topology'],
            mode='wb',
            truncate=self.truncate_gsd,
            filter=hoomd.filter.All()
        )
        gsd_writer.logger = self.logger_hoomd
        
        # Write initial frame
        gsd_writer.write(self.sim.state, filename=f'{self.name}-{self.replica}.gsd',
                         mode='wb', filter=hoomd.filter.All(), logger=self.logger_hoomd)
        
        # Add GSD writer to simulation
        self.sim.operations.writers.append(gsd_writer)
        self.log_info(f"GSD writer added for file: {self.name}-{self.replica}.gsd")
        self.log_info(f"  GSD output period: {self.gsd_output_period_ps:.3f} ps")
        self.log_info(f"  GSD truncate mode: {self.truncate_gsd} ({'overwrite existing file' if self.truncate_gsd else 'append to existing file'})")
        
        # Create a separate logger for console output with only performance and time metrics
        console_logger = hoomd.logging.Logger(categories=['scalar', 'string'])
        
        # Basic simulation quantities
        console_logger.add(self.sim, quantities=['timestep', 'tps'])
        
        # Time and performance information
        console_logger[('Time', 'elapsed_ps')] = (self.time_tracker, 'elapsed_time', 'scalar')
        console_logger[('Performance', 'ns_per_day')] = (self.performance_tracker, 'ns_per_day', 'string')
        console_logger[('Performance', 'eta')] = (self.performance_tracker, 'eta_remaining', 'string')
        console_logger[('Timestep', 'dt_fs')] = (self.timestep_formatter, 'dt_fs', 'scalar')
        
        # Add adaptive timestep logging if enabled
        if hasattr(self, 'adaptive_action') and self.adaptive_action is not None:
            console_logger[('Adaptive', 'error_tolerance')] = (self.adaptive_action, 'error_tolerance', 'scalar')
        
        # Set up console output table with appropriate trigger
        table = hoomd.write.Table(
            trigger=console_trigger,
            logger=console_logger
        )
        self.sim.operations.writers.append(table)
        self.log_info(f"Console output period: {self.console_output_period_ps:.3f} ps")
        self.log_info("Console output restricted to performance and time metrics only")
        
        # Update/remove previous warnings since the issue is now much better
        if self.error_tolerance > 0:
            self.log_info("✅ IMPROVED: GSD and console output now use smaller step intervals in adaptive mode")
            self.log_info("  Output timing accuracy significantly improved (though may still vary slightly)")
            self.log_info("  F(k,t) and energy trackers use precise time-based logic and are fully accurate")
        else:
            self.log_info("Using step-based triggers for maximum efficiency in fixed timestep mode")

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
        
        # Run the simulation
        self.sim.run(total_steps, write_at_start=True)
        
        self.log_info(f"Simulation completed successfully")

    def cleanup(self):
        """Handle post-simulation cleanup and restore original directory."""
        # Note: Trackers from analysis.py write directly to files, no buffering needed
        self.log_info("Cleanup initiated...")
        
        # Restore original directory
        import os
        if hasattr(self, 'original_cwd'):
            os.chdir(self.original_cwd)
            self.log_info(f"Returned to original directory: {self.original_cwd}")
        
        self.log_info("Cleanup completed")

def get_slurm_info():
    """Get SLURM job information from environment variables."""
    task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    job_id = os.environ.get('SLURM_JOB_ID', 'unknown')
    
    if task_id is None:
        return None, job_id
    else:
        return int(task_id), job_id

def parse_replicas(replicas_str):
    """Parse replica specification string into list of integers."""
    if not replicas_str:
        return [1]  # Default single replica
    
    replicas = []
    parts = replicas_str.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            start, end = part.split('-', 1)
            start, end = int(start.strip()), int(end.strip())
            replicas.extend(range(start, end + 1))
        else:
            replicas.append(int(part))
    return sorted(list(set(replicas)))

def run_single_experiment(molecular_thermo, cavity_thermo, finite_q, 
                         coupling, temperature, frequency, replica, frame, 
                         runtime_ps, molecular_tau, cavity_tau, enable_fkt, fkt_kmag, fkt_wavevectors, 
                         fkt_ref_interval, fkt_max_refs, max_energy_output_time=None, 
                         device='CPU', gpu_id=0, incavity=True, fixed_timestep=False, 
                         timestep_fs=1.0, enable_energy_tracking=False, 
                         energy_output_period_ps=0.1, fkt_output_period_ps=1.0, 
                         gsd_output_period_ps=50.0, console_output_period_ps=1.0, truncate_gsd=False):
    """
    Run a single experiment using the CavityMDSimulation class.
    """
    
    try:
        # Create experiment directory with appropriate naming
        if incavity:
            # For cavity simulations, include coupling strength in directory name
            coupling_str = f"{coupling:.0e}".replace("-", "neg").replace("+", "pos")
            exp_dir = Path(f"cavity_coupling_{coupling_str}")
        else:
            # For non-cavity simulations
            exp_dir = Path("no_cavity")
        exp_dir.mkdir(exist_ok=True)
        
        print(f"Running experiment:")
        print(f"  Cavity coupling: {'Enabled' if incavity else 'Disabled'}")
        if incavity:
            print(f"  Coupling strength: {coupling}")
            print(f"  Molecular thermostat: {molecular_thermo}")
            print(f"  Cavity thermostat: {cavity_thermo}")
            print(f"  Finite-q mode: {finite_q}")
        else:
            print(f"  Molecular thermostat: {molecular_thermo}")
        print(f"  Replica: {replica}")
        print(f"  Frame: {frame}")
        print(f"  Output directory: {exp_dir}")
        
        # Set error tolerance based on timestepping mode
        error_tolerance = 0.0 if fixed_timestep else 1.0
        
        # Set timestep based on user preference (only used if fixed_timestep is True)
        dt_fs = timestep_fs if fixed_timestep else None
        
        # Create and run CavityMDSimulation
        sim = CavityMDSimulation(
            job_dir=str(exp_dir),
            replica=replica,
            freq=frequency,
            couplstr=coupling,
            incavity=incavity,
            runtime_ps=runtime_ps,
            input_gsd='../init-0.gsd',  # Use relative path from job directory back to parent
            frame=frame,
            name='prod',
            error_tolerance=error_tolerance,
            temperature=temperature,
            molecular_thermostat=molecular_thermo,
            cavity_thermostat=cavity_thermo,
            finite_q=finite_q,
            molecular_thermostat_tau=molecular_tau,
            cavity_thermostat_tau=cavity_tau,
            log_level='INFO',
            custom_log_file=None,
            enable_fkt=enable_fkt,
            fkt_kmag=fkt_kmag,
            fkt_num_wavevectors=fkt_wavevectors,
            fkt_reference_interval_ps=fkt_ref_interval,
            fkt_max_references=fkt_max_refs,
            max_energy_output_time_ps=max_energy_output_time,
            enable_energy_tracking=enable_energy_tracking,
            dt_fs=dt_fs,
            device=device,
            gpu_id=gpu_id,
            energy_output_period_ps=energy_output_period_ps,
            fkt_output_period_ps=fkt_output_period_ps,
            gsd_output_period_ps=gsd_output_period_ps,
            console_output_period_ps=console_output_period_ps,
            enable_text_output=False,
            text_output_file=None,
            truncate_gsd=truncate_gsd
        )
        
        # Run the simulation
        return sim.run() == 0  # Return True for success (exit code 0)
        
    except Exception as e:
        print(f"ERROR: Experiment failed: {e}")
        return False

def main():
    """Simplified main function for cavity MD experiments."""
    parser = argparse.ArgumentParser(
        description='Advanced Cavity MD Experiment Runner',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Basic simulation parameters
    parser.add_argument('--molecular-bath', type=str, default='bussi', choices=['bussi', 'langevin', 'none'], 
                       help='Molecular thermostat type (default: bussi)')
    parser.add_argument('--cavity-bath', type=str, default='langevin', choices=['bussi', 'langevin', 'none'], 
                       help='Cavity thermostat type (default: langevin)')
    parser.add_argument('--finite-q', action='store_true', 
                       help='Use finite-q cavity mode (default: q=0 mode)')
    parser.add_argument('--coupling', type=float, default=1e-3, 
                       help='Cavity coupling strength (default: 1e-3)')
    parser.add_argument('--temperature', type=float, default=100.0, 
                       help='Temperature in K (default: 100.0)')
    parser.add_argument('--frequency', type=float, default=2000.0, 
                       help='Cavity frequency in cm⁻¹ (default: 2000.0)')
    parser.add_argument('--runtime', type=float, default=500.0, 
                       help='Runtime in ps (default: 500.0)')
    parser.add_argument('--no-cavity', action='store_true', 
                       help='Disable cavity coupling (molecular-only simulation)')
    
    # Replica control
    parser.add_argument('--replicas', type=str, 
                       help='Replica specification (e.g., "1,2,3" or "1-5")')
    
    # Thermostat parameters
    parser.add_argument('--molecular-tau', type=float, default=5.0, 
                       help='Molecular thermostat tau in ps (default: 5.0)')
    parser.add_argument('--cavity-tau', type=float, default=5.0, 
                       help='Cavity thermostat tau in ps (default: 5.0)')
    
    # Timestep control
    parser.add_argument('--fixed-timestep', action='store_true', 
                       help='Use fixed timestep instead of adaptive')
    parser.add_argument('--timestep', type=float, default=1.0, 
                       help='Fixed timestep in fs (default: 1.0)')
    
    # Energy tracking
    parser.add_argument('--enable-energy-tracker', action='store_true', 
                       help='Enable comprehensive energy tracking')
    
    # Output control options - separate periods for different observables
    parser.add_argument('--energy-output-period-ps', type=float, default=0.1, 
                       help='Energy tracker output period in ps (default: 0.1)')
    parser.add_argument('--fkt-output-period-ps', type=float, default=1.0, 
                       help='F(k,t) tracker output period in ps (default: 1.0)')
    parser.add_argument('--gsd-output-period-ps', type=float, default=50.0, 
                       help='GSD trajectory output period in ps (default: 50.0)')
    parser.add_argument('--console-output-period-ps', type=float, default=1.0, 
                       help='Console output period in ps (default: 1.0)')
    
    # F(k,t) options
    parser.add_argument('--enable-fkt', action='store_true', 
                       help='Enable F(k,t) density correlation calculation')
    parser.add_argument('--fkt-kmag', type=float, default=1.0, 
                       help='F(k,t) k magnitude (default: 1.0)')
    parser.add_argument('--fkt-wavevectors', type=int, default=50, 
                       help='F(k,t) number of wavevectors (default: 50)')
    parser.add_argument('--fkt-ref-interval', type=float, default=1.0, 
                       help='F(k,t) reference interval in ps (default: 1.0)')
    parser.add_argument('--fkt-max-refs', type=int, default=10, 
                       help='F(k,t) maximum references (default: 10)')
    parser.add_argument('--max-energy-output-time', type=float, 
                       help='Maximum energy output time in ps (default: no limit)')
    
    # Device options
    parser.add_argument('--device', type=str, default='CPU', choices=['CPU', 'GPU'], 
                       help='Compute device (default: CPU)')
    parser.add_argument('--gpu-id', type=int, default=0, 
                       help='GPU ID when using GPU device (default: 0)')
    
    # GSD output control
    parser.add_argument('--truncate-gsd', action='store_true', 
                       help='Truncate GSD output file if it exists (default: append)')
    
    args = parser.parse_args()
    
    print("Advanced Cavity MD Experiment Runner")
    print("="*50)
    
    # Determine replica list
    task_id, job_id = get_slurm_info()
    
    if task_id is not None:
        # Running under SLURM array job
        replica_list = [task_id]
        print(f"SLURM array job detected: Task {task_id} (Job {job_id})")
    else:
        # Local execution - parse replicas
        replica_list = parse_replicas(args.replicas)
        print(f"Local execution: Replicas {replica_list}")
    
    # Set up simulation parameters
    incavity = not args.no_cavity
    molecular_thermo = args.molecular_bath
    cavity_thermo = args.cavity_bath if incavity else 'none'
    finite_q = args.finite_q
    
    print(f"\nSimulation Configuration:")
    print(f"  Cavity coupling: {'Enabled' if incavity else 'Disabled'}")
    if incavity:
        print(f"    Coupling strength: {args.coupling}")
        print(f"    Frequency: {args.frequency} cm⁻¹")
        print(f"    Finite-q mode: {finite_q}")
        print(f"    Cavity thermostat: {cavity_thermo}")
    print(f"  Molecular thermostat: {molecular_thermo}")
    print(f"  Temperature: {args.temperature} K")
    print(f"  Runtime: {args.runtime} ps")
    print(f"  Device: {args.device}")
    if args.device == 'GPU':
        print(f"    GPU ID: {args.gpu_id}")
    
    # Set up device configuration
    device = args.device.upper()
    
    # Performance tracking
    start_time = time.time()
    successful_experiments = 0
    failed_experiments = 0
    
    print(f"\nStarting execution for {len(replica_list)} replica(s)...")
    print("="*50)
    
    # Run replicas
    for replica in replica_list:
        frame = replica  # Use replica as frame number
        
        print(f"\nRunning replica {replica}...")
        
        # Run experiment
        success = run_single_experiment(
            molecular_thermo=molecular_thermo,
            cavity_thermo=cavity_thermo,
            finite_q=finite_q,
            coupling=args.coupling,
            temperature=args.temperature,
            frequency=args.frequency,
            replica=replica,
            frame=frame,
            runtime_ps=args.runtime,
            molecular_tau=args.molecular_tau,
            cavity_tau=args.cavity_tau,
            enable_fkt=args.enable_fkt,
            fkt_kmag=args.fkt_kmag,
            fkt_wavevectors=args.fkt_wavevectors,
            fkt_ref_interval=args.fkt_ref_interval,
            fkt_max_refs=args.fkt_max_refs,
            max_energy_output_time=args.max_energy_output_time,
            device=device,
            gpu_id=args.gpu_id,
            incavity=incavity,
            fixed_timestep=args.fixed_timestep,
            timestep_fs=args.timestep,
            enable_energy_tracking=args.enable_energy_tracker,
            energy_output_period_ps=args.energy_output_period_ps,
            fkt_output_period_ps=args.fkt_output_period_ps,
            gsd_output_period_ps=args.gsd_output_period_ps,
            console_output_period_ps=args.console_output_period_ps,
            truncate_gsd=args.truncate_gsd
        )
        
        if success:
            successful_experiments += 1
            print(f"SUCCESS: Replica {replica} completed successfully")
        else:
            failed_experiments += 1
            print(f"ERROR: Replica {replica} failed")
    
    # Final summary
    end_time = time.time()
    total_wall_time = end_time - start_time
    total_experiments = len(replica_list)
    
    print("\n" + "="*50)
    print("Execution Summary")
    print("="*50)
    print(f"Total replicas: {total_experiments}")
    print(f"Successful: {successful_experiments}")
    print(f"Failed: {failed_experiments}")
    print(f"Wall time: {total_wall_time:.2f} seconds")
    
    if failed_experiments > 0:
        print(f"\nWARNING: {failed_experiments} replicas failed - check individual logs for details")
        return 1
    else:
        print("\nAll replicas completed successfully!")
        return 0

if __name__ == '__main__':
    sys.exit(main())