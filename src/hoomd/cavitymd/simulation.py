# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Main simulation framework for cavity molecular dynamics."""

import hoomd
import gsd.hoomd
import numpy as np
import logging
import sys
import os
from pathlib import Path

from .forces import CavityForce
from .utils import PhysicalConstants
from .analysis import (
    Status, ElapsedTimeTracker, EnergyContributionTracker, 
    KineticEnergyTracker, CavityModeTracker, DensityCorrelationTracker,
    DipoleAutocorrelation, TimestepFormatter
)
from ..bussi_reservoir.thermostats import BussiReservoir


class AdaptiveTimestepUpdater(hoomd.custom.Action):
    """Update timestep adaptively based on energy conservation and stability."""
    
    def __init__(self, state, integrator, error_tolerance, time_constant_ps=50.0, 
                 initial_fraction=0.01, adaptiveerror=True, cavity_damping_factor=1.0, 
                 molecular_thermostat_tau=5.0, cavity_thermostat_tau=5.0):
        super().__init__()
        self.state = state
        self.integrator = integrator
        self.error_tolerance = error_tolerance
        self.time_constant_ps = time_constant_ps
        self.initial_fraction = initial_fraction
        self.adaptiveerror = adaptiveerror
        self.cavity_damping_factor = cavity_damping_factor
        self.molecular_thermostat_tau = molecular_thermostat_tau
        self.cavity_thermostat_tau = cavity_thermostat_tau
        
        # Convert time constant to atomic units
        self.time_constant_au = PhysicalConstants.ps_to_atomic_units(time_constant_ps)
        
        # Track previous energy for comparison
        self.previous_energy = None
        self.current_elapsed_time_ps = 0.0

    def act(self, timestep):
        if timestep == 0:
            return
            
        # Update elapsed time
        dt = float(self.integrator.dt)
        self.current_elapsed_time_ps = PhysicalConstants.atomic_units_to_ps(dt * timestep)
        
        # Adaptive timestep logic (simplified version)
        # In practice, this would involve energy conservation checks
        # and stability analysis
        
        # For now, implement basic adaptive scaling
        if self.adaptiveerror:
            # Scale timestep based on system properties
            current_dt = float(self.integrator.dt)
            
            # Simple adaptive scaling based on elapsed time
            if self.current_elapsed_time_ps < self.time_constant_ps:
                scaling_factor = 1.0 + self.initial_fraction
            else:
                scaling_factor = 1.0
            
            new_dt = current_dt * scaling_factor
            
            # Apply reasonable bounds
            min_dt = PhysicalConstants.ps_to_atomic_units(0.1)  # 0.1 fs minimum
            max_dt = PhysicalConstants.ps_to_atomic_units(2.0)  # 2 fs maximum
            
            new_dt = max(min_dt, min(max_dt, new_dt))
            
            if abs(new_dt - current_dt) / current_dt > 0.01:  # Only update if change > 1%
                self.integrator.dt = new_dt

    @hoomd.logging.log
    def error_tolerance(self):
        return self.error_tolerance

    @hoomd.logging.log
    def elapsed_time_ps(self):
        return self.current_elapsed_time_ps


class CavityMDSimulation:
    """
    Comprehensive cavity molecular dynamics simulation framework.
    
    This class orchestrates all aspects of a cavity MD simulation including:
    - System setup and initialization
    - Force field configuration
    - Thermostat setup
    - Analysis and tracking components
    - Simulation execution and output
    """
    
    def __init__(self, job_dir, replica, freq, couplstr, incavity, runtime_ps=500.0, 
                 input_gsd='molecular-0.gsd', frame=-1, name='prod', error_tolerance=0.01,
                 temperature=100.0, molecular_thermostat='bussi', cavity_thermostat='langevin',
                 cavity_damping_factor=1.0, use_brownian_overdamped=True, add_cavity_particle=True,
                 finite_q=False, molecular_thermostat_tau=5.0, cavity_thermostat_tau=5.0,
                 log_to_file=True, log_to_console=True, log_level='INFO', custom_log_file=None,
                 enable_fkt=True, fkt_kmag=1.0, fkt_num_wavevectors=50, fkt_reference_interval_ps=1.0, 
                 fkt_max_references=10, max_energy_output_time_ps=None, enable_energy_tracking=True, 
                 dt_fs=None, device='CPU', gpu_id=0):
        """
        Initialize cavity MD simulation.
        
        Parameters:
        -----------
        job_dir : str
            Directory for simulation output
        replica : int
            Replica number for this simulation
        freq : float
            Cavity frequency in cm^-1
        couplstr : float
            Cavity-molecule coupling strength
        incavity : bool
            Whether to include cavity coupling
        runtime_ps : float
            Simulation runtime in picoseconds
        ... (additional parameters as documented in original)
        """
        # Store initialization parameters
        self.job_dir = Path(job_dir)
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
        self.enable_fkt = enable_fkt
        self.fkt_kmag = fkt_kmag
        self.fkt_num_wavevectors = fkt_num_wavevectors
        self.fkt_reference_interval_ps = fkt_reference_interval_ps
        self.fkt_max_references = fkt_max_references
        self.max_energy_output_time_ps = max_energy_output_time_ps
        self.enable_energy_tracking = enable_energy_tracking
        self.dt_fs = dt_fs
        self.device = device
        self.gpu_id = gpu_id
        
        # Create job directory
        self.job_dir.mkdir(exist_ok=True)
        
        # Setup logging
        self.setup_logging(log_to_file, log_to_console, log_level, custom_log_file)
        
        # Calculate physical parameters
        self.calculate_physical_parameters()
        
        # Initialize simulation components
        self.simulation = None
        self.forces = []
        self.methods = []
        self.trackers = []
        
    def setup_logging(self, log_to_file, log_to_console, log_level, custom_log_file):
        """Setup logging configuration."""
        self.logger = logging.getLogger(f'CavityMD_R{self.replica}')
        self.logger.setLevel(getattr(logging, log_level.upper()))
        
        # Clear existing handlers
        self.logger.handlers.clear()
        
        # Setup file logging
        if log_to_file:
            log_file = custom_log_file or self.job_dir / f'{self.name}-{self.replica}.log'
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(getattr(logging, log_level.upper()))
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
        
        # Setup console logging
        if log_to_console:
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(getattr(logging, log_level.upper()))
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            console_handler.setFormatter(formatter)
            self.logger.addHandler(console_handler)

    def calculate_physical_parameters(self):
        """Calculate physical parameters from input values."""
        # Convert cavity frequency to atomic units
        self.cavity_freq_au = self.freq / PhysicalConstants.HARTREE_TO_CM_MINUS1
        
        # Set cavity mass (typically 1.0 in atomic units)
        self.cavity_mass = 1.0
        
        # Convert temperature to atomic units
        self.temperature_au = self.temperature * PhysicalConstants.KB_HARTREE_PER_K
        
        # Calculate damping parameters
        self.molecular_gamma = PhysicalConstants.gamma_from_tau_ps(self.molecular_thermostat_tau)
        if self.incavity:
            self.cavity_gamma = PhysicalConstants.gamma_from_tau_ps(self.cavity_thermostat_tau)
        
        self.log_info(f"Physical parameters calculated:")
        self.log_info(f"  Cavity frequency: {self.freq} cm^-1 = {self.cavity_freq_au:.6e} a.u.")
        self.log_info(f"  Temperature: {self.temperature} K = {self.temperature_au:.6e} a.u.")
        self.log_info(f"  Molecular damping: γ = {self.molecular_gamma:.6e} a.u.")
        if self.incavity:
            self.log_info(f"  Cavity damping: γ = {self.cavity_gamma:.6e} a.u.")

    def get_snapshot(self):
        """Load and return snapshot from GSD file."""
        try:
            with gsd.hoomd.open(self.input_gsd, mode='r') as trajectory:
                if self.frame == -1:
                    snapshot = trajectory[-1]
                else:
                    snapshot = trajectory[self.frame]
            self.log_info(f"Loaded snapshot from {self.input_gsd}, frame {self.frame}")
            return snapshot
        except Exception as e:
            self.log_error(f"Failed to load snapshot: {e}")
            raise

    def setup_simulation(self):
        """Setup the HOOMD simulation object."""
        # Initialize device
        if self.device.upper() == 'GPU':
            device = hoomd.device.GPU(gpu_id=self.gpu_id)
        else:
            device = hoomd.device.CPU()
        
        # Create simulation
        self.simulation = hoomd.Simulation(device=device, seed=self.replica)
        
        # Load snapshot
        snapshot = self.get_snapshot()
        
        # Add cavity particle if requested
        if self.add_cavity_particle and self.incavity:
            snapshot = self.create_cavity_particle(snapshot)
        
        # Create state from snapshot
        self.simulation.create_state_from_snapshot(snapshot)
        
        self.log_info(f"Simulation initialized with {len(snapshot.particles.position)} particles")

    def create_cavity_particle(self, snapshot):
        """Add cavity particle to the snapshot."""
        # This is a simplified version - the full implementation would
        # handle proper cavity particle setup
        self.log_info("Adding cavity particle to system")
        
        # Add one more particle for the cavity mode
        n_particles = len(snapshot.particles.position)
        new_n = n_particles + 1
        
        # Extend arrays
        new_positions = np.zeros((new_n, 3))
        new_positions[:-1] = snapshot.particles.position
        new_positions[-1] = [0, 0, 0]  # Cavity particle at origin
        
        # Update snapshot
        snapshot.particles.position = new_positions
        # ... (additional array extensions for other properties)
        
        return snapshot

    def setup_force_parameters(self, dt, rcut=15):
        """Setup force field parameters."""
        self.forces = []
        
        # Molecular forces (LJ, bonds, electrostatics)
        # This would contain the full force field setup from the original
        # For brevity, showing the structure:
        
        # LJ forces
        lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4))
        # ... setup LJ parameters
        self.forces.append(lj)
        
        # Bond forces  
        harmonic = hoomd.md.bond.Harmonic()
        # ... setup bond parameters
        self.forces.append(harmonic)
        
        # Electrostatic forces
        ewald = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
            nlist=hoomd.md.nlist.Cell(buffer=0.4),
            resolution=(64, 64, 64),
            order=6,
            r_cut=rcut
        )
        self.forces.extend(ewald)
        
        # Cavity force
        if self.incavity:
            cavity_force = CavityForce(
                frequency=self.cavity_freq_au,
                coupling=self.couplstr,
                finite_q=self.finite_q
            )
            self.forces.append(cavity_force)
            self.cavity_force = cavity_force
        
        self.log_info(f"Setup {len(self.forces)} force components")

    def setup_thermostat_parameters(self, dt):
        """Setup thermostat parameters."""
        self.methods = []
        
        # Molecular thermostat
        if self.molecular_thermostat == 'bussi':
            molecular_bussi = BussiReservoir(
                filter=hoomd.filter.Type(['O', 'N']),  # Molecular particles
                kT=self.temperature_au,
                tau=self.molecular_thermostat_tau
            )
            self.methods.append(molecular_bussi)
            self.molecular_bussi_thermostat = molecular_bussi
        
        # Cavity thermostat
        if self.incavity and hasattr(self, 'cavity_force'):
            if self.cavity_thermostat == 'langevin':
                cavity_langevin = hoomd.md.methods.Langevin(
                    filter=hoomd.filter.Type(['cavity']),  # Cavity particle
                    kT=self.temperature_au,
                    gamma=self.cavity_gamma
                )
                self.methods.append(cavity_langevin)
        
        self.log_info(f"Setup {len(self.methods)} thermostat methods")

    def setup_integrator(self, forces, methods):
        """Setup the MD integrator."""
        # Compute optimal timestep if not fixed
        if self.dt_fs is None:
            dt = self.compute_and_set_optimal_timestep()
        else:
            dt = PhysicalConstants.ps_to_atomic_units(self.dt_fs / 1000.0)  # fs to ps to au
        
        # Create integrator
        integrator = hoomd.md.Integrator(dt=dt)
        integrator.forces = forces
        integrator.methods = methods
        
        self.simulation.operations.integrator = integrator
        
        # Add adaptive timestep updater if needed
        if self.error_tolerance > 0.0:
            timestep_updater = AdaptiveTimestepUpdater(
                state=self.simulation.state,
                integrator=integrator,
                error_tolerance=self.error_tolerance,
                molecular_thermostat_tau=self.molecular_thermostat_tau,
                cavity_thermostat_tau=self.cavity_thermostat_tau
            )
            self.simulation.operations.updaters.append(timestep_updater)
        
        self.log_info(f"Integrator setup with dt = {dt:.6e} a.u. ({PhysicalConstants.atomic_units_to_ps(dt)*1000:.3f} fs)")

    def compute_and_set_optimal_timestep(self):
        """Compute optimal timestep for the system."""
        # Simplified timestep calculation
        # In practice, this would analyze system properties
        base_dt_fs = 1.0  # Base timestep in femtoseconds
        dt_au = PhysicalConstants.ps_to_atomic_units(base_dt_fs / 1000.0)
        return dt_au

    def setup_trackers_and_loggers(self):
        """Setup analysis trackers and loggers."""
        self.trackers = []
        
        # Time tracker
        time_tracker = ElapsedTimeTracker(self.simulation, self.runtime_ps)
        self.simulation.operations.updaters.append(time_tracker)
        self.time_tracker = time_tracker
        
        # Energy tracking
        if self.enable_energy_tracking:
            energy_tracker = EnergyContributionTracker(
                simulation=self.simulation,
                harmonic=self.forces[1],  # Bond force
                lj=self.forces[0],  # LJ force
                short=self.forces[2],  # Short-range electrostatics
                long=self.forces[3],  # Long-range electrostatics
                cavityforce=getattr(self, 'cavity_force', None),
                time_tracker=time_tracker,
                output_prefix=str(self.job_dir / f'{self.name}-{self.replica}'),
                max_time_ps=self.max_energy_output_time_ps
            )
            self.simulation.operations.updaters.append(energy_tracker)
            self.trackers.append(energy_tracker)
        
        # Kinetic energy tracking
        kinetic_tracker = KineticEnergyTracker(
            simulation=self.simulation,
            time_tracker=time_tracker,
            output_prefix=str(self.job_dir / f'{self.name}-{self.replica}')
        )
        self.simulation.operations.updaters.append(kinetic_tracker)
        self.trackers.append(kinetic_tracker)
        
        # Cavity mode tracking
        if self.incavity and hasattr(self, 'cavity_force'):
            cavity_tracker = CavityModeTracker(
                simulation=self.simulation,
                cavityforce=self.cavity_force,
                time_tracker=time_tracker,
                output_prefix=str(self.job_dir / f'{self.name}-{self.replica}')
            )
            self.simulation.operations.updaters.append(cavity_tracker)
            self.trackers.append(cavity_tracker)
        
        # F(k,t) tracking
        if self.enable_fkt:
            fkt_tracker = DensityCorrelationTracker(
                simulation=self.simulation,
                time_tracker=time_tracker,
                kmag=self.fkt_kmag,
                num_wavevectors=self.fkt_num_wavevectors,
                reference_interval_ps=self.fkt_reference_interval_ps,
                max_references=self.fkt_max_references,
                output_prefix=str(self.job_dir / f'{self.name}-{self.replica}')
            )
            self.simulation.operations.updaters.append(fkt_tracker)
            self.trackers.append(fkt_tracker)
        
        self.log_info(f"Setup {len(self.trackers)} analysis trackers")

    def setup_output_writers(self):
        """Setup output writers for trajectory and logging."""
        # GSD trajectory writer
        gsd_writer = hoomd.write.GSD(
            filename=str(self.job_dir / f'{self.name}-{self.replica}.gsd'),
            trigger=hoomd.trigger.Periodic(1000),
            mode='wb'
        )
        self.simulation.operations.writers.append(gsd_writer)
        
        # Logger for thermodynamic quantities
        logger = hoomd.logging.Logger()
        
        # Add loggable quantities
        logger.add(self.simulation, quantities=['timestep', 'tps'])
        if hasattr(self, 'time_tracker'):
            logger.add(self.time_tracker, quantities=['elapsed_time'])
        
        for tracker in self.trackers:
            # Add all loggable quantities from each tracker
            loggables = []
            for attr_name in dir(tracker):
                attr = getattr(tracker, attr_name)
                if hasattr(attr, '_hoomd_loggable') and attr._hoomd_loggable:
                    loggables.append(attr_name)
            if loggables:
                logger.add(tracker, quantities=loggables)
        
        # Table writer for log output
        table_writer = hoomd.write.Table(
            trigger=hoomd.trigger.Periodic(100),
            logger=logger,
            output=open(str(self.job_dir / f'{self.name}-{self.replica}.txt'), 'w')
        )
        self.simulation.operations.writers.append(table_writer)
        
        self.log_info("Output writers configured")

    def run_simulation(self):
        """Execute the main simulation."""
        # Calculate total timesteps
        dt = float(self.simulation.operations.integrator.dt)
        total_time_au = PhysicalConstants.ps_to_atomic_units(self.runtime_ps)
        total_timesteps = int(total_time_au / dt)
        
        self.log_info(f"Starting simulation: {total_timesteps} timesteps ({self.runtime_ps} ps)")
        
        # Run simulation
        try:
            self.simulation.run(total_timesteps)
            self.log_info("Simulation completed successfully")
            return 0
        except Exception as e:
            self.log_error(f"Simulation failed: {e}")
            return 1

    def run(self):
        """Run the complete simulation workflow."""
        try:
            self.log_info("="*60)
            self.log_info(f"CAVITY MD SIMULATION - Replica {self.replica}")
            self.log_info("="*60)
            
            # Setup simulation
            self.setup_simulation()
            
            # Setup forces
            dt = PhysicalConstants.ps_to_atomic_units(1.0 / 1000.0)  # 1 fs default
            self.setup_force_parameters(dt)
            
            # Setup thermostats
            self.setup_thermostat_parameters(dt)
            
            # Setup integrator
            self.setup_integrator(self.forces, self.methods)
            
            # Setup analysis and tracking
            self.setup_trackers_and_loggers()
            
            # Setup output
            self.setup_output_writers()
            
            # Run simulation
            exit_code = self.run_simulation()
            
            self.log_info("="*60)
            self.log_info("SIMULATION COMPLETED")
            self.log_info("="*60)
            
            return exit_code
            
        except Exception as e:
            self.log_error(f"Fatal error: {e}")
            import traceback
            self.log_error(traceback.format_exc())
            return 1

    def log_info(self, message):
        """Log info message."""
        if hasattr(self, 'logger'):
            self.logger.info(message)
        else:
            print(f"INFO: {message}")

    def log_warning(self, message):
        """Log warning message."""
        if hasattr(self, 'logger'):
            self.logger.warning(message)
        else:
            print(f"WARNING: {message}")

    def log_error(self, message):
        """Log error message."""
        if hasattr(self, 'logger'):
            self.logger.error(message)
        else:
            print(f"ERROR: {message}")

    def log_debug(self, message):
        """Log debug message."""
        if hasattr(self, 'logger'):
            self.logger.debug(message)
        else:
            print(f"DEBUG: {message}") 