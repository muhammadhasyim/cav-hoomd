# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Analysis and tracking components for cavity molecular dynamics simulations."""

import hoomd
import datetime
import numpy as np
import traceback
import sys

from .utils import PhysicalConstants, unwrap_positions

# =============================================================================
# OBSERVABLE LIBRARY
# =============================================================================

def compute_total_dipole_moment(snapshot):
    """Compute total dipole moment with proper position unwrapping."""
    box_lengths = np.array([
        snapshot.global_box.L[0],
        snapshot.global_box.L[1], 
        snapshot.global_box.L[2]
    ])
    unwrapped_positions = unwrap_positions(
        snapshot.particles.position,
        snapshot.particles.image,
        box_lengths
    )
    # Dipole = charges × positions
    return np.dot(snapshot.particles.charge, unwrapped_positions)


def compute_density_field(snapshot, wavevectors):
    """Compute density field ρ(k) for given wavevectors."""
    positions = snapshot.particles.position
    rhok_real = np.zeros(len(wavevectors))
    rhok_imag = np.zeros(len(wavevectors))
    
    for i, k_vec in enumerate(wavevectors):
        # Compute k·r for all particles
        kr = np.dot(positions, k_vec)
        # Compute ρ(k) = sum_j exp(i k·r_j)
        rhok_real[i] = np.sum(np.cos(kr))
        rhok_imag[i] = np.sum(np.sin(kr))
    
    return rhok_real + 1j * rhok_imag


def generate_fibonacci_sphere(samples=100):
    """Generate uniformly distributed points on a sphere using Fibonacci spiral."""
    points = np.zeros((samples, 3))
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * radius
        z = np.sin(theta) * radius

        points[i] = [x, y, z]

    return points


# Observable library definitions
SIMPLE_OBSERVABLES = {
    "dipole": compute_total_dipole_moment,
}

FIELD_OBSERVABLES = {
    "density_correlation": compute_density_field,
    # Add other field observables as needed
}

ENERGY_COMPONENTS = {
    "harmonic": lambda forces: forces.get("harmonic", None),
    "lj": lambda forces: forces.get("lj", None),
    "ewald_short": lambda forces: forces.get("ewald_short", None),
    "ewald_long": lambda forces: forces.get("ewald_long", None),
    "cavity": lambda forces: forces.get("cavity", None),
    # NOTE: Separate cavity components use the same cavity force object
    "cavity_harmonic": lambda forces: forces.get("cavity", None),
    "cavity_coupling": lambda forces: forces.get("cavity", None),
    "cavity_dipole_self": lambda forces: forces.get("cavity", None),
}

# Reservoir energy components for thermostat energy tracking
RESERVOIR_ENERGY_COMPONENTS = {
    "bussi_molecular_reservoir": lambda thermostats: thermostats.get("bussi_molecular", None),
    "bussi_cavity_reservoir": lambda thermostats: thermostats.get("bussi_cavity", None),
    "langevin_molecular_reservoir": lambda thermostats: thermostats.get("langevin_molecular", None),
    "langevin_cavity_reservoir": lambda thermostats: thermostats.get("langevin_cavity", None),
}


# =============================================================================
# BASE TRACKER CLASS
# =============================================================================

class BaseTracker(hoomd.custom.Action):
    """Base class for all tracking components with common infrastructure."""
    
    def __init__(self, simulation, time_tracker=None, output_prefix='tracker', output_period_steps=1000):
        """Initialize base tracker.
        
        Args:
            simulation: HOOMD simulation object
            time_tracker: Optional time tracker for accurate timing
            output_prefix: Prefix for output files
            output_period_steps: Output frequency in simulation steps (NOT time-based)
        """
        super().__init__()
        self.sim = simulation
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix
        self.output_period_steps = output_period_steps  # Always in steps, never time-based
        self.last_output_step = 0
        
        # Initialize current values for logging
        self._initialize_logging_values()
    
    def _initialize_logging_values(self):
        """Initialize values that will be logged. Override in subclasses."""
        pass
    
    def _get_current_time(self, timestep):
        """Get current simulation time."""
        if self.time_tracker is not None:
            return self.time_tracker.elapsed_time
        else:
            return timestep * self.sim.operations.integrator.dt * 0.02418884  # Convert to ps
    
    def _should_output(self, timestep):
        """Check if we should output at this timestep."""
        return timestep - self.last_output_step >= self.output_period_steps
    
    def _update_output_step(self, timestep):
        """Update the last output step."""
        self.last_output_step = timestep
    



# =============================================================================
# AUTOCORRELATION TRACKER (Simple Observables)
# =============================================================================

class AutocorrelationTracker(BaseTracker):
    """Generic autocorrelation tracker for simple observables."""
    
    def __init__(self, simulation, observable, time_tracker=None, output_prefix=None, output_period_steps=1000):
        """Initialize autocorrelation tracker.
        
        Args:
            simulation: HOOMD simulation object
            observable: Observable name from SIMPLE_OBSERVABLES
            time_tracker: Optional time tracker for accurate timing
            output_prefix: Prefix for output files
            output_period_steps: Output frequency in simulation steps
        """
        if observable not in SIMPLE_OBSERVABLES:
            raise ValueError(f"Unknown observable '{observable}'. Available: {list(SIMPLE_OBSERVABLES.keys())}")
        
        self.observable = observable
        self.observable_func = SIMPLE_OBSERVABLES[observable]
        
        if output_prefix is None:
            output_prefix = f'{observable}_autocorr'
        
        super().__init__(simulation, time_tracker, output_prefix, output_period_steps)
        
        # Initialize autocorrelation tracking
        self.output_file_number = 0
        self.output_file_path = f'{self.output_prefix}_{self.output_file_number}.txt'
        self.reference_time = 0.0
        self.reference_value = None
        self.current_autocorr_value = None
        
        # Initialize reference value and output file
        self._initialize_reference()
    
    def _initialize_reference(self):
        """Initialize reference value at t=0."""
        with self.sim.state.cpu_local_snapshot as snap:
            self.reference_value = self.observable_func(snap)
            self.current_autocorr_value = np.dot(self.reference_value, self.reference_value)
            
            # Write header and t=0 value
            with open(self.output_file_path, 'w') as f:
                f.write(f'# {self.observable.capitalize()} autocorrelation data\n')
                f.write(f'# Reference number: {self.output_file_number}\n')
                f.write(f'# Output period: {self.output_period_steps} steps\n')
                f.write('# timestep t(ps) C(t)\n')
                f.write(f'{0} {0.0:.6f} {self.current_autocorr_value:.6f}\n')
        
        print(f"{self.observable.capitalize()} autocorrelation tracker initialized. C(0) = {self.current_autocorr_value:.6e}")
        print(f"Output period: {self.output_period_steps} steps")
        print(f"Writing to file: {self.output_file_path}")
    
    def _initialize_logging_values(self):
        """Initialize logging values."""
        self.current_autocorr_value = 0.0
    
    def compute_autocorr(self, current_value):
        """Compute autocorrelation C(t) = observable(0)·observable(t)."""
        return np.dot(self.reference_value, current_value)
    
    def _start_new_reference(self, timestep):
        """Start a new reference file."""
        self.output_file_number += 1
        self.output_file_path = f'{self.output_prefix}_{self.output_file_number}.txt'
        with self.sim.state.cpu_local_snapshot as snap:
            self.reference_value = self.observable_func(snap)
            self.current_autocorr_value = np.dot(self.reference_value, self.reference_value)
            
            # Write header and t=0 value for new reference
            with open(self.output_file_path, 'w') as f:
                f.write(f'# {self.observable.capitalize()} autocorrelation data\n')
                f.write(f'# Reference number: {self.output_file_number}\n')
                f.write(f'# Output period: {self.output_period_steps} steps\n')
                f.write('# timestep t(ps) C(t)\n')
                f.write(f'{timestep} {self._get_current_time(timestep):.6f} {self.current_autocorr_value:.6f}\n')
    
    def act(self, timestep):
        if timestep == 0:
            return
            
        with self.sim.state.cpu_local_snapshot as snap:
            current_value = self.observable_func(snap)
            autocorr_value = self.compute_autocorr(current_value)
            self.current_autocorr_value = autocorr_value
            
            # Output periodically
            if self._should_output(timestep):
                current_time = self._get_current_time(timestep)
                with open(self.output_file_path, 'a') as f:
                    f.write(f'{timestep} {current_time:.6f} {autocorr_value:.6f}\n')
                self._update_output_step(timestep)
                
                # Start new reference file every 10000 steps
                if timestep % 10000 == 0:
                    self._start_new_reference(timestep)
    
    @hoomd.logging.log
    def current_autocorr(self):
        return self.current_autocorr_value if self.current_autocorr_value is not None else 0.0


# =============================================================================
# FIELD AUTOCORRELATION TRACKER (Field Observables)
# =============================================================================

class FieldAutocorrelationTracker(BaseTracker):
    """Generic autocorrelation tracker for field observables with k-space averaging."""
    
    def __init__(self, simulation, observable, time_tracker=None, output_prefix=None, 
                 output_period_steps=1000, reference_interval_steps=10000, max_references=10, **kwargs):
        """Initialize field autocorrelation tracker.
        
        Args:
            simulation: HOOMD simulation object
            observable: Observable name from FIELD_OBSERVABLES
            time_tracker: Optional time tracker for accurate timing
            output_prefix: Prefix for output files
            output_period_steps: Output frequency in simulation steps
            reference_interval_steps: Interval between new references in steps
            max_references: Maximum number of reference states to keep
            **kwargs: Observable-specific parameters
        """
        if observable not in FIELD_OBSERVABLES:
            raise ValueError(f"Unknown field observable '{observable}'. Available: {list(FIELD_OBSERVABLES.keys())}")
        
        self.observable = observable
        self.observable_func = FIELD_OBSERVABLES[observable]
        self.reference_interval_steps = reference_interval_steps
        self.max_references = max_references
        
        if output_prefix is None:
            output_prefix = f'{observable}_field_autocorr'
        
        super().__init__(simulation, time_tracker, output_prefix, output_period_steps)
        
        # Setup observable-specific parameters
        self._setup_parameters(kwargs)
        
        # Initialize field autocorrelation tracking
        self.references = []  # List of reference states
        self.last_reference_step = 0
        
        # Initialize first reference and output files
        self._initialize_new_reference_file(0)
    
    def _setup_parameters(self, kwargs):
        """Setup parameters specific to the observable."""
        if self.observable == 'density_correlation':
            # Default parameters for density correlation
            self.kmag = kwargs.get('kmag', 1.0)
            self.num_wavevectors = kwargs.get('num_wavevectors', 50)
            self.wavevectors = generate_fibonacci_sphere(self.num_wavevectors) * self.kmag
            print(f"Density correlation: k={self.kmag:.2f}, {self.num_wavevectors} wavevectors")
        else:
            # Add other observables as needed
            pass
    
    def _call_observable_func(self, snapshot):
        """Call the observable function with correct arguments."""
        if self.observable == 'density_correlation':
            return self.observable_func(snapshot, self.wavevectors)
        else:
            # For other observables, use the function directly
            return self.observable_func(snapshot)
    
    def _initialize_new_reference_file(self, ref_number):
        """Initialize a new reference file."""
        ref_filename = f'{self.output_prefix}_ref{ref_number}.txt'
        
        # Initialize reference state
        with self.sim.state.cpu_local_snapshot as snap:
            reference_field = self._call_observable_func(snap)
            current_time = self._get_current_time(self.sim.timestep)
            
            # Store reference
            self.references.append({
                'number': ref_number,
                'filename': ref_filename,
                'timestep': self.sim.timestep,
                'time': current_time,
                'field': reference_field
            })
            
            # Create output file with header
            with open(ref_filename, 'w') as f:
                f.write(f'# {self.observable.capitalize()} field autocorrelation\n')
                f.write(f'# Reference {ref_number} at t={current_time:.6f} ps\n')
                f.write(f'# Output period: {self.output_period_steps} steps\n')
                f.write('# timestep lag_time(ps) field_autocorr\n')
        
        print(f"Initialized {self.observable} field autocorr reference {ref_number}")
    
    def _initialize_logging_values(self):
        """Initialize logging values."""
        self.current_autocorr_value = 0.0
    
    def compute_field_autocorr(self, field0, field_t):
        """Compute field autocorrelation."""
        if isinstance(field0, np.ndarray) and isinstance(field_t, np.ndarray):
            return np.mean(np.real(field0 * np.conj(field_t)))
        else:
            return np.dot(field0, field_t)
    
    def act(self, timestep):
        # Get current time FIRST, outside snapshot context
        current_time = self._get_current_time(timestep)
        
        if timestep == 0:
            return
        
        # Compute current field value
        with self.sim.state.cpu_local_snapshot as snap:
            current_field = self._call_observable_func(snap)
        
        # Update autocorrelations for all active references
        for ref in self.references:
            lag_time = current_time - ref['time']
            autocorr_value = self.compute_field_autocorr(ref['field'], current_field)
            
            # Update current autocorr for first reference (for logging)
            if ref['number'] == 0:
                self.current_autocorr_value = autocorr_value
            
            # Output periodically
            if self._should_output(timestep):
                with open(ref['filename'], 'a') as f:
                    f.write(f'{timestep} {lag_time:.6f} {autocorr_value:.6f}\n')
        
        # Add new reference if interval has passed and we haven't hit max
        if (timestep - self.last_reference_step >= self.reference_interval_steps and 
            len(self.references) < self.max_references):
            ref_number = len(self.references)
            self._initialize_new_reference_file(ref_number)
            self.last_reference_step = timestep
        
        # Update output step
        if self._should_output(timestep):
            self._update_output_step(timestep)
    
    @hoomd.logging.log
    def current_autocorr(self):
        return self.current_autocorr_value if self.current_autocorr_value is not None else 0.0


# =============================================================================
# CORRECTED ENERGY TRACKER - Follows Working EnergyContributationTracker Pattern
# =============================================================================

class EnergyTracker(BaseTracker):
    """
    Energy tracking that exactly follows the working EnergyContributionTracker pattern.
    
    **KEY FIXES to match working code:**
    
    1. **Direct Force Energy Access**: Use direct `.energy` access like working code,
       no complex force computation initialization
    
    2. **Separate Kinetic Trackers**: Use separate kinetic energy trackers instead 
       of computing from snapshot (matches working code architecture)
    
    3. **Simple Error Handling**: Use simple AttributeError/DataAccessError handling
       like working code, not complex initialization
    
    4. **Exact Universe Total Calculation**: Follow working code's exact calculation:
       universe_total = system_total + total_thermostat_energy
    
    5. **Extensive Debug Output**: Add comprehensive debugging to identify issues
    """
    
    def __init__(self, simulation, components, force_objects=None, thermostat_objects=None, 
                 kinetic_tracker=None, cavity_mode_tracker=None,
                 time_tracker=None, output_prefix='energy', output_period_steps=1000, 
                 max_timesteps=None, max_time_ps=None, compute_temperature=True, track_reservoirs=True, verbose='normal'):
        """Initialize corrected energy tracker following working code pattern.
        
        Args:
            simulation: HOOMD simulation object
            components: List of energy components to track
            force_objects: Dictionary of force objects (harmonic, lj, ewald_short, ewald_long, cavity)
            thermostat_objects: Dictionary of thermostat/method objects for reservoir energy
            kinetic_tracker: KineticEnergyTracker object (like working code)
            cavity_mode_tracker: CavityModeTracker object (like working code)  
            time_tracker: ElapsedTimeTracker for accurate timing
            output_prefix: Prefix for output files
            output_period_steps: Output frequency in simulation steps
            max_timesteps: Maximum timesteps to track (ignored if max_time_ps is set)
            max_time_ps: Maximum simulation time in ps to track (more accurate than max_timesteps)
            compute_temperature: Whether to compute temperature
            track_reservoirs: Whether to track reservoir energies
            verbose: Verbosity level ('quiet', 'normal', 'verbose')
                    - 'quiet': No debug output, only essential messages
                    - 'normal': Minimal output, setup info and errors
                    - 'verbose': Full debug output (original behavior)
        """
        # Store configuration exactly like working code
        self.force_objects = force_objects or {}
        self.thermostat_objects = thermostat_objects or {}
        self.kinetic_tracker = kinetic_tracker  # Like working code
        self.cavity_mode_tracker = cavity_mode_tracker  # Like working code
        self.track_reservoirs = track_reservoirs
        self.max_timesteps = max_timesteps
        self.max_time_ps = max_time_ps
        self.compute_temperature = compute_temperature
        self.output_stopped = False
        
        # Set verbosity level
        self.verbose = verbose.lower() if isinstance(verbose, str) else 'normal'
        if self.verbose not in ['quiet', 'normal', 'verbose']:
            self.verbose = 'normal'
        
        # Validate components
        self.components = components
        
        super().__init__(simulation, time_tracker, output_prefix, output_period_steps)
        
        # Fix file naming to match working code
        self.output_file_path = f'{self.output_prefix}_energy_tracker.txt'
        
        # Only print setup info if not quiet
        if self.verbose != 'quiet':
            print(f"CORRECTED EnergyTracker (following working code pattern):")
            print(f"  Output file: {self.output_file_path}")
            print(f"  Components: {self.components}")
            print(f"  Force objects: {list(self.force_objects.keys())}")
            print(f"  Thermostat objects: {list(self.thermostat_objects.keys())}")
            print(f"  Kinetic tracker: {self.kinetic_tracker is not None}")
            print(f"  Cavity mode tracker: {self.cavity_mode_tracker is not None}")
            print(f"  Track reservoirs: {self.track_reservoirs}")
            print(f"  Output period: {self.output_period_steps} steps")
            print(f"  Verbosity: {self.verbose}")
            if self.max_time_ps:
                print(f"  Max time: {self.max_time_ps} ps (time-based limit)")
            elif self.max_timesteps:
                print(f"  Max timesteps: {self.max_timesteps} (step-based limit)")
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_logging_values(self):
        """Initialize energy values for logging."""
        # Current energy components
        self.current_harmonic_energy = 0.0
        self.current_lj_energy = 0.0
        self.current_ewald_short_energy = 0.0
        self.current_ewald_long_energy = 0.0
        self.current_cavity_harmonic_energy = 0.0
        self.current_cavity_coupling_energy = 0.0
        self.current_cavity_dipole_self_energy = 0.0
        self.current_cavity_total_potential_energy = 0.0
        
        # Kinetic energies
        self.current_molecular_kinetic_energy = 0.0
        self.current_cavity_kinetic_energy = 0.0
        self.current_total_kinetic_energy = 0.0
        
        # Reservoir energies
        self.current_molecular_reservoir_energy = 0.0
        self.current_cavity_reservoir_energy = 0.0
        self.current_total_reservoir_energy = 0.0
        
        # Totals
        self.current_total_potential_energy = 0.0
        self.current_system_total_energy = 0.0
        self.current_universe_total_energy = 0.0
        self.current_temperature = 0.0
    
    def _initialize_output_file(self):
        """Initialize output file with headers matching working code."""
        try:
            with open(self.output_file_path, 'w') as f:
                f.write('# CORRECTED Energy tracking following working EnergyContributionTracker pattern\n')
                f.write(f'# Output period: {self.output_period_steps} steps\n')  
                if self.max_time_ps:
                    f.write(f'# Max time: {self.max_time_ps} ps\n')
                elif self.max_timesteps:
                    f.write(f'# Max timesteps: {self.max_timesteps}\n')
                f.write('# All energies in Hartree (atomic units)\n')
                f.write('# Column definitions:\n')
                f.write('#   time(ps): simulation time in picoseconds\n')
                f.write('#   timestep: simulation timestep number\n')
                f.write('#   harmonic_energy: harmonic bond potential energy\n')
                f.write('#   lj_energy: Lennard-Jones potential energy\n')
                f.write('#   ewald_short_energy: short-range Coulomb energy\n')
                f.write('#   ewald_long_energy: long-range Coulomb energy\n')
                f.write('#   cavity_harmonic_energy: cavity harmonic potential energy\n')
                f.write('#   cavity_coupling_energy: cavity-molecule coupling energy\n')
                f.write('#   cavity_dipole_self_energy: dipole self-energy\n')
                f.write('#   cavity_total_potential_energy: total cavity potential energy\n')
                f.write('#   molecular_kinetic_energy: molecular kinetic energy\n')
                f.write('#   cavity_kinetic_energy: cavity kinetic energy\n')
                f.write('#   total_kinetic_energy: total kinetic energy\n')
                f.write('#   total_potential_energy: total potential energy\n')
                f.write('#   system_total_energy: total system energy (KE + PE)\n')
                f.write('#   molecular_reservoir_energy: molecular reservoir energy\n')
                f.write('#   cavity_reservoir_energy: cavity reservoir energy\n')
                f.write('#   total_reservoir_energy: total reservoir energy\n')
                f.write('#   universe_total_energy: universe total energy (system + reservoir) [CONSERVED]\n')
                if self.compute_temperature:
                    f.write('#   temperature: kinetic temperature (K)\n')
                
                # Create header line
                header = 'time(ps) timestep'
                header += ' harmonic_energy lj_energy ewald_short_energy ewald_long_energy'
                header += ' cavity_harmonic_energy cavity_coupling_energy cavity_dipole_self_energy cavity_total_potential_energy'
                header += ' molecular_kinetic_energy cavity_kinetic_energy total_kinetic_energy'
                header += ' total_potential_energy system_total_energy'
                header += ' molecular_reservoir_energy cavity_reservoir_energy total_reservoir_energy'
                header += ' universe_total_energy'
                if self.compute_temperature:
                    header += ' temperature'
                header += '\n'
                f.write(header)
                
            if self.verbose != 'quiet':
                print(f"EnergyTracker: Successfully created output file {self.output_file_path}")
        except Exception as e:
            # Always print errors regardless of verbosity
            print(f"EnergyTracker ERROR: Failed to create output file {self.output_file_path}: {e}")
    
    def act(self, timestep):
        """Main energy computation method following working code pattern exactly."""
        # Check if output has been stopped due to time or timestep limit
        if self.output_stopped:
            return
            
        if timestep == 0:
            return
        
        # Get current time for checking time-based limits
        if self.time_tracker is not None:
            current_time = self.time_tracker.elapsed_time
        else:
            current_time = timestep * self.sim.operations.integrator.dt * 0.02418884  # Convert to ps
            
        # Check time limit first (more accurate than timestep limit)
        if self.max_time_ps is not None:
            if current_time > self.max_time_ps:
                if not self.output_stopped:
                    self.output_stopped = True
                    if self.verbose != 'quiet':
                        print(f"Energy tracking stopped: reached time limit of {self.max_time_ps:.2f} ps at t={current_time:.4f} ps")
                return
        
        # Check timestep limit if specified and no time limit
        elif self.max_timesteps is not None:
            if timestep > self.max_timesteps:
                if not self.output_stopped:
                    self.output_stopped = True
                    if self.verbose != 'quiet':
                        print(f"Energy tracking stopped at timestep {timestep} (limit: {self.max_timesteps})")
                return
        
        # Only output periodically (matching working code pattern)
        if timestep - self.last_output_step < self.output_period_steps:
            return
            
        try:
            if self.verbose == 'verbose':
                print(f"\n=== ENERGY TRACKER DEBUG - Timestep {timestep} ===")
                print(f"Current time: {current_time:.6f} ps")
            
            # === 1. GET POTENTIAL ENERGY COMPONENTS (exactly like working code) ===
            if self.verbose == 'verbose':
                print("=== POTENTIAL ENERGY COMPONENTS ===")
            
            # Get individual potential energy contributions (direct access like working code)
            try:
                self.current_harmonic_energy = self.force_objects['harmonic'].energy if 'harmonic' in self.force_objects else 0.0
                if self.verbose == 'verbose':
                    print(f"Harmonic energy: {self.current_harmonic_energy:.6f} Hartree")
            except (AttributeError, KeyError) as e:
                self.current_harmonic_energy = 0.0
                if self.verbose in ['normal', 'verbose']:
                    print(f"Harmonic energy ERROR: {e}")
                
            try:
                self.current_lj_energy = self.force_objects['lj'].energy if 'lj' in self.force_objects else 0.0
                if self.verbose == 'verbose':
                    print(f"LJ energy: {self.current_lj_energy:.6f} Hartree")
            except (AttributeError, KeyError) as e:
                self.current_lj_energy = 0.0
                if self.verbose in ['normal', 'verbose']:
                    print(f"LJ energy ERROR: {e}")
                
            try:
                self.current_ewald_short_energy = self.force_objects['ewald_short'].energy if 'ewald_short' in self.force_objects else 0.0
                if self.verbose == 'verbose':
                    print(f"Ewald short energy: {self.current_ewald_short_energy:.6f} Hartree")
            except (AttributeError, KeyError) as e:
                self.current_ewald_short_energy = 0.0
                if self.verbose in ['normal', 'verbose']:
                    print(f"Ewald short energy ERROR: {e}")
                
            try:            
                self.current_ewald_long_energy = self.force_objects['ewald_long'].energy if 'ewald_long' in self.force_objects else 0.0
                if self.verbose == 'verbose':
                    print(f"Ewald long energy: {self.current_ewald_long_energy:.6f} Hartree")
            except (AttributeError, KeyError) as e:
                self.current_ewald_long_energy = 0.0
                if self.verbose in ['normal', 'verbose']:
                    print(f"Ewald long energy ERROR: {e}")
            
            # Calculate total potential energy (without cavity)
            molecular_potential_energy = (self.current_harmonic_energy + self.current_lj_energy + 
                                        self.current_ewald_short_energy + self.current_ewald_long_energy)
            if self.verbose == 'verbose':
                print(f"Molecular potential energy (harmonic + lj + ewald): {molecular_potential_energy:.6f} Hartree")
            
            # Get cavity potential energy components if present (exactly like working code)
            self.current_cavity_harmonic_energy = 0.0
            self.current_cavity_coupling_energy = 0.0
            self.current_cavity_dipole_self_energy = 0.0
            self.current_cavity_total_potential_energy = 0.0
            
            if 'cavity' in self.force_objects and self.force_objects['cavity'] is not None:
                cavityforce = self.force_objects['cavity']
                if self.verbose == 'verbose':
                    print("=== CAVITY ENERGY COMPONENTS ===")
                try:
                    # Use exact same pattern as working code
                    self.current_cavity_harmonic_energy = getattr(cavityforce, 'harmonic_energy', 0.0)
                    self.current_cavity_coupling_energy = getattr(cavityforce, 'coupling_energy', 0.0)
                    self.current_cavity_dipole_self_energy = getattr(cavityforce, 'dipole_self_energy', 0.0)
                    
                    if self.verbose == 'verbose':
                        print(f"Cavity harmonic energy: {self.current_cavity_harmonic_energy:.6f} Hartree")
                        print(f"Cavity coupling energy: {self.current_cavity_coupling_energy:.6f} Hartree")
                        print(f"Cavity dipole self energy: {self.current_cavity_dipole_self_energy:.6f} Hartree")
                    
                    # For total energy, try .energy property first, then sum components (exactly like working code)
                    if hasattr(cavityforce, 'energy'):
                        self.current_cavity_total_potential_energy = cavityforce.energy
                        if self.verbose == 'verbose':
                            print(f"Cavity total energy (from .energy): {self.current_cavity_total_potential_energy:.6f} Hartree")
                    else:
                        self.current_cavity_total_potential_energy = (self.current_cavity_harmonic_energy + 
                                                                   self.current_cavity_coupling_energy + 
                                                                   self.current_cavity_dipole_self_energy)
                        if self.verbose == 'verbose':
                            print(f"Cavity total energy (sum components): {self.current_cavity_total_potential_energy:.6f} Hartree")
                except Exception as e:
                    # If any error occurs, set all to zero and continue (exactly like working code)
                    self.current_cavity_harmonic_energy = 0.0
                    self.current_cavity_coupling_energy = 0.0
                    self.current_cavity_dipole_self_energy = 0.0
                    self.current_cavity_total_potential_energy = 0.0
                    if self.verbose in ['normal', 'verbose']:
                        print(f"ERROR accessing cavity energy components: {e}")
            else:
                if self.verbose == 'verbose':
                    print("No cavity force object - cavity energies set to zero")
            
            # Calculate total potential energy (exactly like working code)
            self.current_total_potential_energy = molecular_potential_energy + self.current_cavity_total_potential_energy
            if self.verbose == 'verbose':
                print(f"TOTAL POTENTIAL ENERGY: {self.current_total_potential_energy:.6f} Hartree")
            
            # === 2. GET KINETIC ENERGY COMPONENTS (exactly like working code) ===  
            if self.verbose == 'verbose':
                print("=== KINETIC ENERGY COMPONENTS ===")
            
            # Get molecular kinetic energy from kinetic tracker (exactly like working code)
            self.current_molecular_kinetic_energy = 0.0
            if self.kinetic_tracker is not None:
                try:
                    self.current_molecular_kinetic_energy = self.kinetic_tracker.kinetic_energy
                    if self.verbose == 'verbose':
                        print(f"Molecular kinetic energy (from tracker): {self.current_molecular_kinetic_energy:.6f} Hartree")
                except AttributeError as e:
                    self.current_molecular_kinetic_energy = 0.0
                    if self.verbose in ['normal', 'verbose']:
                        print(f"Molecular kinetic energy ERROR: {e}")
            else:
                if self.verbose == 'verbose':
                    print("No kinetic tracker - molecular kinetic energy set to zero")
            
            # Get cavity kinetic energy from cavity mode tracker (exactly like working code)
            self.current_cavity_kinetic_energy = 0.0
            if self.cavity_mode_tracker is not None:
                try:
                    self.current_cavity_kinetic_energy = self.cavity_mode_tracker.cavity_kinetic_energy
                    if self.verbose == 'verbose':
                        print(f"Cavity kinetic energy (from tracker): {self.current_cavity_kinetic_energy:.6f} Hartree")
                except AttributeError as e:
                    self.current_cavity_kinetic_energy = 0.0
                    if self.verbose in ['normal', 'verbose']:
                        print(f"Cavity kinetic energy ERROR: {e}")
            else:
                if self.verbose == 'verbose':
                    print("No cavity mode tracker - cavity kinetic energy set to zero")
            
            # Calculate total kinetic energy (exactly like working code)
            self.current_total_kinetic_energy = self.current_molecular_kinetic_energy + self.current_cavity_kinetic_energy
            if self.verbose == 'verbose':
                print(f"TOTAL KINETIC ENERGY: {self.current_total_kinetic_energy:.6f} Hartree")
            
            # === 3. GET RESERVOIR ENERGIES (exactly like working code) ===
            if self.verbose == 'verbose':
                print("=== RESERVOIR ENERGY COMPONENTS ===")
            
            # Get molecular reservoir energy if available (exactly like working code)
            self.current_molecular_reservoir_energy = 0.0
            
            # Check for molecular Langevin method
            if 'langevin_molecular' in self.thermostat_objects:
                try:
                    mol_langevin_reservoir = self.thermostat_objects['langevin_molecular'].reservoir_energy
                    self.current_molecular_reservoir_energy += mol_langevin_reservoir
                    if self.verbose == 'verbose':
                        print(f"Molecular Langevin reservoir energy: {mol_langevin_reservoir:.6f} Hartree")
                except AttributeError:
                    if self.verbose == 'verbose':
                        print("Molecular Langevin reservoir energy not available yet")
                    
            # Check for molecular Bussi thermostat
            if 'bussi_molecular' in self.thermostat_objects:
                try:
                    mol_bussi_reservoir = self.thermostat_objects['bussi_molecular'].total_reservoir_energy
                    self.current_molecular_reservoir_energy += mol_bussi_reservoir
                    if self.verbose == 'verbose':
                        print(f"Molecular Bussi reservoir energy: {mol_bussi_reservoir:.6f} Hartree")
                except (AttributeError, hoomd.error.DataAccessError):
                    if self.verbose == 'verbose':
                        print("Molecular Bussi reservoir energy not available yet")
            
            # Get cavity reservoir energy if available (exactly like working code)
            self.current_cavity_reservoir_energy = 0.0
            
            # Check for cavity Langevin method
            if 'langevin_cavity' in self.thermostat_objects:
                try:
                    cav_langevin_reservoir = self.thermostat_objects['langevin_cavity'].reservoir_energy
                    self.current_cavity_reservoir_energy += cav_langevin_reservoir
                    if self.verbose == 'verbose':
                        print(f"Cavity Langevin reservoir energy: {cav_langevin_reservoir:.6f} Hartree")
                except AttributeError:
                    if self.verbose == 'verbose':
                        print("Cavity Langevin reservoir energy not available yet")
                    
            # Check for cavity Bussi thermostat
            if 'bussi_cavity' in self.thermostat_objects:
                try:
                    cav_bussi_reservoir = self.thermostat_objects['bussi_cavity'].total_reservoir_energy
                    self.current_cavity_reservoir_energy += cav_bussi_reservoir
                    if self.verbose == 'verbose':
                        print(f"Cavity Bussi reservoir energy: {cav_bussi_reservoir:.6f} Hartree")
                except (AttributeError, hoomd.error.DataAccessError):
                    if self.verbose == 'verbose':
                        print("Cavity Bussi reservoir energy not available yet")
            
            # Calculate total reservoir energy (exactly like working code)
            self.current_total_reservoir_energy = self.current_molecular_reservoir_energy + self.current_cavity_reservoir_energy
            if self.verbose == 'verbose':
                print(f"TOTAL RESERVOIR ENERGY: {self.current_total_reservoir_energy:.6f} Hartree")
            
            # === 4. CALCULATE TOTAL ENERGIES (exactly like working code) ===
            if self.verbose == 'verbose':
                print("=== TOTAL ENERGY CALCULATIONS ===")
            
            # Calculate system total energy (exactly like working code)
            self.current_system_total_energy = self.current_total_potential_energy + self.current_total_kinetic_energy
            if self.verbose == 'verbose':
                print(f"SYSTEM TOTAL ENERGY (KE + PE): {self.current_system_total_energy:.6f} Hartree")
            
            # Calculate universe total energy (exactly like working code)
            # This is the conserved quantity: system energy + reservoir energy
            self.current_universe_total_energy = self.current_system_total_energy + self.current_total_reservoir_energy
            if self.verbose == 'verbose':
                print(f"UNIVERSE TOTAL ENERGY (system + reservoir): {self.current_universe_total_energy:.6f} Hartree")
                print(f"  (This should be conserved)")
            
            # Calculate temperature if requested
            if self.compute_temperature:
                # Use kinetic energy from trackers to compute temperature
                total_particles = 0
                if self.kinetic_tracker is not None:
                    try:
                        self.current_temperature = self.kinetic_tracker.temperature
                        if self.verbose == 'verbose':
                            print(f"Temperature (from kinetic tracker): {self.current_temperature:.2f} K")
                    except AttributeError:
                        self.current_temperature = 0.0
                        if self.verbose == 'verbose':
                            print("Temperature not available from kinetic tracker")
                else:
                    self.current_temperature = 0.0
                    if self.verbose == 'verbose':
                        print("No kinetic tracker - temperature set to zero")
            
            # === 5. WRITE OUTPUT DATA ===
            if self.verbose == 'verbose':
                print("=== WRITING OUTPUT DATA ===")
            self._write_energy_data(timestep, current_time)
            self.last_output_step = timestep
            
            if self.verbose == 'verbose':
                print(f"=== END ENERGY TRACKER DEBUG - Timestep {timestep} ===\n")
            
        except Exception as e:
            # Always print critical errors regardless of verbosity
            print(f"EnergyTracker CRITICAL ERROR at timestep {timestep}: {e}")
            import traceback
            traceback.print_exc()

    def _write_energy_data(self, timestep, current_time):
        """Write energy data to output file."""
        try:
            # Build output line exactly like working code format
            output_values = [
                current_time, timestep,
                # Potential energy components
                self.current_harmonic_energy,
                self.current_lj_energy, 
                self.current_ewald_short_energy,
                self.current_ewald_long_energy,
                self.current_cavity_harmonic_energy,
                self.current_cavity_coupling_energy,
                self.current_cavity_dipole_self_energy,
                self.current_cavity_total_potential_energy,
                # Kinetic energy components
                self.current_molecular_kinetic_energy,
                self.current_cavity_kinetic_energy,
                self.current_total_kinetic_energy,
                # Total energies
                self.current_total_potential_energy,
                self.current_system_total_energy,
                # Reservoir energies
                self.current_molecular_reservoir_energy,
                self.current_cavity_reservoir_energy,
                self.current_total_reservoir_energy,
                # Universe total (conserved quantity)
                self.current_universe_total_energy
            ]
            
            if self.compute_temperature:
                output_values.append(self.current_temperature)
            
            # Write to file
            with open(self.output_file_path, 'a') as f:
                formatted_values = [f'{val:.6f}' if isinstance(val, float) else str(val) for val in output_values]
                f.write(' '.join(formatted_values) + '\n')
                
            if self.verbose == 'verbose':
                print(f"Successfully wrote energy data to {self.output_file_path}")
                
        except Exception as e:
            # Always print errors regardless of verbosity
            print(f"EnergyTracker ERROR writing data at timestep {timestep}: {e}")
            import traceback
            traceback.print_exc()

    # Logging methods for HOOMD logger integration (matching working code)
    @hoomd.logging.log
    def total_energy(self):
        return self.current_system_total_energy
    
    @hoomd.logging.log
    def universe_total_energy(self):
        return self.current_universe_total_energy
    
    @hoomd.logging.log
    def total_potential_energy(self):
        return self.current_total_potential_energy
    
    @hoomd.logging.log
    def kinetic_energy(self):
        return self.current_total_kinetic_energy
    
    @hoomd.logging.log
    def total_reservoir_energy(self):
        return self.current_total_reservoir_energy
        
    @hoomd.logging.log
    def temperature(self):
        return self.current_temperature
    
    @hoomd.logging.log
    def harmonic_energy(self):
        return self.current_harmonic_energy
    
    @hoomd.logging.log
    def lj_energy(self):
        return self.current_lj_energy
    
    @hoomd.logging.log
    def ewald_short_energy(self):
        return self.current_ewald_short_energy
    
    @hoomd.logging.log
    def ewald_long_energy(self):
        return self.current_ewald_long_energy
    
    @hoomd.logging.log
    def cavity_harmonic_energy(self):
        return self.current_cavity_harmonic_energy
    
    @hoomd.logging.log
    def cavity_coupling_energy(self):
        return self.current_cavity_coupling_energy
    
    @hoomd.logging.log
    def cavity_dipole_self_energy(self):
        return self.current_cavity_dipole_self_energy
    
    @hoomd.logging.log
    def molecular_kinetic_energy(self):
        return self.current_molecular_kinetic_energy
    
    @hoomd.logging.log
    def cavity_kinetic_energy_separate(self):
        return self.current_cavity_kinetic_energy
    
    @hoomd.logging.log
    def molecular_reservoir_energy(self):
        return self.current_molecular_reservoir_energy
    
    @hoomd.logging.log
    def cavity_reservoir_energy_separate(self):
        return self.current_cavity_reservoir_energy


# =============================================================================
# UTILITY CLASSES (Keep as-is)
# =============================================================================

class Status:
    """Status monitoring for cavity MD simulations."""
    
    def __init__(self, simulation, chartime, time_tracker=None):
        self.simulation = simulation
        self.chartime = chartime
        self.starttime = datetime.datetime.now()
        self.time_tracker = time_tracker
        self.last_timestep = 0
        self.last_wall_time = datetime.datetime.now()

    @property
    def seconds_remaining(self):
        try:
            return (
                self.simulation.final_timestep - self.simulation.timestep
            ) / self.simulation.tps
        except ZeroDivisionError:
            return 0

    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining))
    
    def etr_string(self):
        """Get estimated time remaining as a string for logging."""
        return str(datetime.timedelta(seconds=self.seconds_remaining))
    
    @property
    def nsd(self):
        # Calculate nanoseconds per day based on actual simulation progress
        current_timestep = self.simulation.timestep
        if current_timestep <= 0:
            return "0.0"
        
        # Use time_tracker if available for more accurate timing
        if self.time_tracker is not None:
            simulation_time_ps = self.time_tracker.elapsed_time
        else:
            # Fallback to dt * timestep calculation
            dt = float(self.simulation.operations.integrator.dt)
            simulation_time_ps = PhysicalConstants.atomic_units_to_ps(dt * current_timestep)
        
        # Calculate wall time elapsed
        current_wall_time = datetime.datetime.now()
        wall_time_elapsed = (current_wall_time - self.starttime).total_seconds()
        
        if wall_time_elapsed <= 0:
            return "0.0"
        
        # Calculate simulation rate: ps per second of wall time
        ps_per_second = simulation_time_ps / wall_time_elapsed
        
        # Convert to nanoseconds per day
        ns_per_second = ps_per_second / 1000.0  # ps to ns conversion
        ns_per_day = ns_per_second * 86400  # seconds to day conversion
        
        return str(np.round(ns_per_day, 6))
    
    def ns_per_day(self):
        """Get nanoseconds per day performance metric for logging."""
        # Calculate nanoseconds per day based on actual simulation progress
        current_timestep = self.simulation.timestep
        if current_timestep <= 0:
            return "0.0"
        
        # Use time_tracker if available for more accurate timing
        if self.time_tracker is not None:
            simulation_time_ps = self.time_tracker.elapsed_time
        else:
            # Fallback to dt * timestep calculation
            dt = float(self.simulation.operations.integrator.dt)
            simulation_time_ps = PhysicalConstants.atomic_units_to_ps(dt * current_timestep)
        
        # Calculate wall time elapsed
        current_wall_time = datetime.datetime.now()
        wall_time_elapsed = (current_wall_time - self.starttime).total_seconds()
        
        if wall_time_elapsed <= 0:
            return "0.0"
        
        # Calculate simulation rate: ps per second of wall time
        ps_per_second = simulation_time_ps / wall_time_elapsed
        
        # Convert to nanoseconds per day
        ns_per_second = ps_per_second / 1000.0  # ps to ns conversion
        ns_per_day = ns_per_second * 86400  # seconds to day conversion
        
        return str(np.round(ns_per_day, 6))
    
    @property
    def Dt(self):
        return str(np.round(float(self.simulation.operations.integrator.dt*self.chartime*1000000),6))
    
    @property
    def elapsed(self):
        curtime = datetime.datetime.now()
        return str(curtime-self.starttime)


class ElapsedTimeTracker(hoomd.custom.Action):
    """Track elapsed simulation time in physical units and exit when runtime is reached."""
    
    def __init__(self, simulation, runtime):
        super().__init__()
        self.simulation = simulation
        self.runtime = runtime  # Target runtime in ps
        self.total_time = 0.0  # Store in atomic units like cavitymd.py
        self.last_timestep = 0  # Track previous timestep for proper accumulation
        self.initial_timestep = None  # Track the starting timestep to handle inheritance

    def act(self, timestep):
        """Update the total elapsed time by accumulating time increments."""
        # Get current timestep size
        dt = self.simulation.operations.integrator.dt
        
        # For the first call, handle initialization
        if self.last_timestep == 0:
            # Initialize - record the starting timestep but don't add its time
            self.initial_timestep = timestep
            self.last_timestep = timestep
            self.total_time = 0.0  # Always start elapsed time from 0, regardless of inherited timestep
            if timestep > 0:
                print(f"NOTICE: Starting from inherited timestep {timestep}")
                print(f"  Elapsed time will start from 0, not from inherited simulation time")
            return
        
        # Calculate time increment since last update
        if timestep > self.last_timestep:
            timestep_increment = timestep - self.last_timestep
            time_increment = timestep_increment * dt
            self.total_time += time_increment
        
        # Update last timestep for next iteration
        self.last_timestep = timestep
        
        # Check if we've reached the runtime and exit if so
        if PhysicalConstants.atomic_units_to_ps(self.total_time) >= self.runtime:
            print(f"Runtime {self.runtime} ps reached. Exiting simulation.")
            import sys
            sys.exit(0)

    @hoomd.logging.log
    def elapsed_time(self):
        """Expose the total elapsed time as a property in (ps)."""
        return PhysicalConstants.atomic_units_to_ps(self.total_time)


class TimestepFormatter(hoomd.custom.Action):
    """Format timestep information for logging."""
    
    def __init__(self, integrator):
        super().__init__()
        self.integrator = integrator

    def act(self, timestep):
        pass  # No action needed, just for logging

    @hoomd.logging.log
    def dt_fs(self):
        """Current timestep size in femtoseconds."""
        dt_au = self.integrator.dt
        dt_fs = PhysicalConstants.atomic_units_to_ps(dt_au) * 1000  # Convert ps to fs
        return dt_fs


class CavityModeTracker(hoomd.custom.Action):
    """Track cavity mode properties and energies (specialized tracker)."""
    
    def __init__(self, simulation, cavityforce, time_tracker=None, output_prefix='cavity_mode', output_period_steps=1000):
        """Initialize cavity mode tracker.
        
        Args:
            simulation: HOOMD simulation object
            cavityforce: Cavity force object to track
            time_tracker: Optional time tracker for accurate timing
            output_prefix: Prefix for output files
            output_period_steps: Output frequency in simulation steps
        """
        super().__init__()
        self.sim = simulation
        self.cavityforce = cavityforce
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix
        self.output_period_steps = output_period_steps
        self.output_file_path = f'{self.output_prefix}_cavity_mode.txt'
        
        # Track last output step
        self.last_output_step = 0
        
        # Initialize cavity values
        self.current_cavity_kinetic_energy = 0.0
        self.current_cavity_potential_energy = 0.0
        self.current_cavity_total_energy = 0.0
        self.current_cavity_temperature = 0.0
        
        # Initialize output file
        with open(self.output_file_path, 'w') as f:
            f.write('# Cavity mode tracking\n')
            f.write(f'# Output period: {self.output_period_steps} steps\n')
            f.write('# timestep time(ps) cavity_kinetic_energy cavity_potential_energy cavity_total_energy cavity_temperature\n')
        
        print(f"CavityModeTracker: Will write to {self.output_file_path}")
        print(f"CavityModeTracker: Output period = {self.output_period_steps} steps")

    def compute_cavity_properties(self):
        """Compute cavity mode kinetic and potential energies."""
        try:
            with self.sim.state.cpu_local_snapshot as snap:
                # Find photon particle (typeid == 2)
                photon_mask = snap.particles.typeid == 2
                
                if not np.any(photon_mask):
                    return 0.0, 0.0, 0.0, 0.0
                
                # Get photon properties
                photon_mass = snap.particles.mass[photon_mask][0]
                photon_velocity = snap.particles.velocity[photon_mask][0]
                
                # Get unwrapped position for the photon
                box_lengths = np.array([
                    snap.global_box.L[0],
                    snap.global_box.L[1],
                    snap.global_box.L[2]
                ])
                unwrapped_positions = unwrap_positions(
                    snap.particles.position,
                    snap.particles.image,
                    box_lengths
                )
                photon_position = unwrapped_positions[photon_mask][0]
                
                # Compute cavity kinetic energy: KE = (1/2) * m * v²
                vel_squared = np.sum(photon_velocity**2)
                kinetic_energy = 0.5 * photon_mass * vel_squared
                
                # Get harmonic potential energy from cavity force (without coupling and self-energy terms)
                if hasattr(self.cavityforce, 'harmonic_energy'):
                    potential_energy = self.cavityforce.harmonic_energy
                else:
                    potential_energy = 0.0
                
                # Total cavity oscillator energy (KE + harmonic PE only)
                total_energy = kinetic_energy + potential_energy
                
                # Temperature from kinetic energy (3 degrees of freedom)
                # T = (2/3) * KE / k_B for 3D
                temperature = (2.0/3.0) * kinetic_energy / PhysicalConstants.KB_HARTREE_PER_K
                
                return kinetic_energy, potential_energy, total_energy, temperature
                
        except Exception as e:
            # Fallback if cavity properties are not accessible
            print(f"CavityModeTracker error: {e}")
            return 0.0, 0.0, 0.0, 0.0

    def act(self, timestep):
        if timestep == 0:
            return
            
        # Compute cavity properties
        kinetic, potential, total, temperature = self.compute_cavity_properties()
        
        # Store current values
        self.current_cavity_kinetic_energy = kinetic
        self.current_cavity_potential_energy = potential
        self.current_cavity_total_energy = total
        self.current_cavity_temperature = temperature
        
        # Output periodically
        if timestep - self.last_output_step >= self.output_period_steps:
            # Get current time
            if self.time_tracker is not None:
                current_time = self.time_tracker.elapsed_time
            else:
                dt = float(self.sim.operations.integrator.dt)
                current_time = PhysicalConstants.atomic_units_to_ps(dt * timestep)
            
            with open(self.output_file_path, 'a') as f:
                f.write(f'{timestep} {current_time:.6f} {kinetic:.6f} {potential:.6f} {total:.6f} {temperature:.6f}\n')
            
            self.last_output_step = timestep

    @hoomd.logging.log
    def cavity_kinetic_energy(self):
        return self.current_cavity_kinetic_energy

    @hoomd.logging.log
    def cavity_potential_energy_harmonic(self):
        return self.current_cavity_potential_energy

    @hoomd.logging.log
    def cavity_total_energy(self):
        return self.current_cavity_total_energy

    @hoomd.logging.log
    def cavity_temperature(self):
        return self.current_cavity_temperature


# =============================================================================
# CONVENIENCE ALIASES AND BACKWARD COMPATIBILITY
# =============================================================================

class DipoleAutocorrelation(AutocorrelationTracker):
    """Dipole autocorrelation tracker - convenience wrapper around AutocorrelationTracker.
    
    This class provides backward compatibility for the DipoleAutocorrelation class
    that was previously defined in the monolithic cavitymd.py file.
    """
    
    def __init__(self, simulation, time_tracker=None, output_prefix='dipole_autocorr', output_period_steps=1000):
        """Initialize dipole autocorrelation tracker.
        
        Args:
            simulation: HOOMD simulation object
            time_tracker: Optional time tracker for accurate timing
            output_prefix: Prefix for output files
            output_period_steps: Output frequency in simulation steps
        """
        super().__init__(
            simulation=simulation, 
            observable="dipole", 
            time_tracker=time_tracker, 
            output_prefix=output_prefix, 
            output_period_steps=output_period_steps
        )