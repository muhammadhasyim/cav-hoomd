# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Analysis and tracking components for cavity molecular dynamics simulations."""

import hoomd
import datetime
import numpy as np
from numba import njit

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
    "density": {
        "compute": compute_density_field,
        "params": ["kmag", "num_wavevectors"],
        "averaging": lambda correlations, N: np.mean(correlations) / N if N > 0 else 0.0
    }
}

ENERGY_COMPONENTS = {
    "harmonic": lambda forces: forces.get("harmonic", None),
    "lj": lambda forces: forces.get("lj", None),
    "ewald_short": lambda forces: forces.get("ewald_short", None),
    "ewald_long": lambda forces: forces.get("ewald_long", None),
    "cavity": lambda forces: forces.get("cavity", None),
}


# =============================================================================
# BASE TRACKER CLASS
# =============================================================================

class BaseTracker(hoomd.custom.Action):
    """Base class for all tracking components with common infrastructure."""
    
    def __init__(self, simulation, time_tracker=None, output_prefix='tracker', output_period=1):
        super().__init__()
        self.sim = simulation
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix
        self.output_period = output_period
        self.last_output_step = 0
        
        # Initialize current values for logging
        self._initialize_logging_values()
    
    def _initialize_logging_values(self):
        """Initialize values that will be logged. Override in subclasses."""
        pass
    
    def _get_current_time(self, timestep):
        """Get current simulation time in picoseconds."""
        if self.time_tracker is not None:
            return self.time_tracker.elapsed_time
        else:
            dt = float(self.sim.operations.integrator.dt)
            return PhysicalConstants.atomic_units_to_ps(dt * timestep)
    
    def _should_output(self, timestep):
        """Check if this timestep should produce output."""
        return timestep - self.last_output_step >= self.output_period
    
    def _update_output_step(self, timestep):
        """Update the last output timestep."""
        self.last_output_step = timestep


# =============================================================================
# AUTOCORRELATION TRACKER (Simple Observables)
# =============================================================================

class AutocorrelationTracker(BaseTracker):
    """Generic autocorrelation tracker for simple observables."""
    
    def __init__(self, simulation, observable, time_tracker=None, output_prefix=None, output_period=1000):
        if observable not in SIMPLE_OBSERVABLES:
            raise ValueError(f"Unknown observable '{observable}'. Available: {list(SIMPLE_OBSERVABLES.keys())}")
        
        self.observable = observable
        self.observable_func = SIMPLE_OBSERVABLES[observable]
        
        if output_prefix is None:
            output_prefix = f'{observable}_autocorr'
        
        super().__init__(simulation, time_tracker, output_prefix, output_period)
        
        # Initialize autocorrelation tracking
        self.output_file_path = f'{self.output_prefix}.txt'
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
                f.write('# t0(ps) t(ps) C(t)\n')
                f.write(f'# Started at time 0.0 ps\n')
                f.write(f'{0.0:.6f} {0.0:.6f} {self.current_autocorr_value:.6f}\n')
        
        print(f"{self.observable.capitalize()} autocorrelation tracker initialized. C(0) = {self.current_autocorr_value:.6e}")
    
    def _initialize_logging_values(self):
        """Initialize logging values."""
        self.current_autocorr_value = 0.0
    
    def compute_autocorr(self, current_value):
        """Compute autocorrelation C(t) = observable(0)·observable(t)."""
        return np.dot(self.reference_value, current_value)
    
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
                    f.write(f'{self.reference_time:.6f} {current_time:.6f} {autocorr_value:.6f}\n')
                self._update_output_step(timestep)
    
    @hoomd.logging.log
    def current_autocorr(self):
        return self.current_autocorr_value if self.current_autocorr_value is not None else 0.0


# =============================================================================
# FIELD AUTOCORRELATION TRACKER (Field Observables)
# =============================================================================

class FieldAutocorrelationTracker(BaseTracker):
    """Generic autocorrelation tracker for field observables with k-space averaging."""
    
    def __init__(self, simulation, observable, time_tracker=None, output_prefix=None, 
                 output_period=1000, reference_interval_ps=1.0, max_references=10, **kwargs):
        if observable not in FIELD_OBSERVABLES:
            raise ValueError(f"Unknown field observable '{observable}'. Available: {list(FIELD_OBSERVABLES.keys())}")
        
        self.observable = observable
        self.observable_config = FIELD_OBSERVABLES[observable]
        self.observable_func = self.observable_config["compute"]
        self.averaging_func = self.observable_config["averaging"]
        
        if output_prefix is None:
            output_prefix = f'{observable}_autocorr'
            
        super().__init__(simulation, time_tracker, output_prefix, output_period)
        
        # Extract and validate parameters
        self._setup_parameters(kwargs)
        
        # Autocorrelation tracking
        self.output_file_path = f'{self.output_prefix}.txt'
        self.reference_interval_ps = reference_interval_ps
        self.max_references = max_references
        self.reference_data = []  # List of (time, field_values)
        self.last_reference_time = -self.reference_interval_ps  # Force first reference at t=0
        self.current_autocorr_value = 0.0
        
        # Initialize output file
        self._initialize_output_file()
    
    def _setup_parameters(self, kwargs):
        """Setup observable-specific parameters."""
        required_params = self.observable_config["params"]
        
        for param in required_params:
            if param not in kwargs:
                raise ValueError(f"Required parameter '{param}' not provided for observable '{self.observable}'")
            setattr(self, param, kwargs[param])
        
        # Setup wavevectors for density observable
        if self.observable == "density":
            self.wavevectors = generate_fibonacci_sphere(self.num_wavevectors) * self.kmag
    
    def _initialize_output_file(self):
        """Initialize output file with headers."""
        with open(self.output_file_path, 'w') as f:
            f.write(f'# {self.observable.capitalize()} autocorrelation tracking\n')
            if self.observable == "density":
                f.write(f'# kmag = {self.kmag}, num_wavevectors = {self.num_wavevectors}\n')
            f.write('# t0(ps) t(ps) F(k,t)\n')
    
    def _initialize_logging_values(self):
        """Initialize logging values."""
        self.current_autocorr_value = 0.0
    
    def compute_field_autocorr(self, field0, field_t):
        """Compute autocorrelation for field observables."""
        if self.observable == "density":
            # For complex density fields: ρ*(k,0) ρ(k,t)
            correlation = np.real(np.conj(field0) * field_t)
            return correlation
        else:
            # Generic real field autocorrelation
            return field0 * field_t
    
    def act(self, timestep):
        # Get current field values
        with self.sim.state.cpu_local_snapshot as snap:
            current_time = self._get_current_time(timestep)
            
            if self.observable == "density":
                current_field = self.observable_func(snap, self.wavevectors)
            else:
                current_field = self.observable_func(snap)
        
        # Check if we need to add a new reference
        if current_time - self.last_reference_time >= self.reference_interval_ps:
            # Add new reference
            self.reference_data.append((current_time, current_field.copy()))
            self.last_reference_time = current_time
            
            # Remove oldest reference if we have too many
            if len(self.reference_data) > self.max_references:
                self.reference_data.pop(0)
        
        # Compute autocorrelation for all references
        if len(self.reference_data) > 0:
            # For logging, use the first (oldest) reference
            ref_time, ref_field = self.reference_data[0]
            correlations = self.compute_field_autocorr(ref_field, current_field)
            
            # Apply averaging function
            with self.sim.state.cpu_local_snapshot as snap:
                num_particles = len(snap.particles.position)
            self.current_autocorr_value = self.averaging_func(correlations, num_particles)
        
        # Output periodically
        if self._should_output(timestep):
            # Output autocorrelation for all reference times
            for ref_time, ref_field in self.reference_data:
                correlations = self.compute_field_autocorr(ref_field, current_field)
                with self.sim.state.cpu_local_snapshot as snap:
                    num_particles = len(snap.particles.position)
                autocorr_value = self.averaging_func(correlations, num_particles)
                
                with open(self.output_file_path, 'a') as f:
                    f.write(f'{ref_time:.6f} {current_time:.6f} {autocorr_value:.6f}\n')
            
            self._update_output_step(timestep)
    
    @hoomd.logging.log
    def current_autocorr(self):
        return self.current_autocorr_value


# =============================================================================
# ENERGY TRACKER (Consolidated Energy Tracking)
# =============================================================================

class EnergyTracker(BaseTracker):
    """Consolidated energy tracking for all simulation components."""
    
    def __init__(self, simulation, components, force_objects=None, time_tracker=None, 
                 output_prefix='energy', output_period=1, max_time_ps=None, 
                 compute_temperature=True):
        super().__init__(simulation, time_tracker, output_prefix, output_period)
        
        self.components = components
        self.force_objects = force_objects or {}
        self.max_time_ps = max_time_ps
        self.compute_temperature = compute_temperature
        self.output_stopped = False
        
        self.output_file_path = f'{self.output_prefix}_energy.txt'
        
        # Initialize current energy values
        self._initialize_energy_values()
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_energy_values(self):
        """Initialize energy tracking values."""
        self.current_energies = {comp: 0.0 for comp in self.components}
        self.current_total_energy = 0.0
        self.current_temperature = 0.0
    
    def _initialize_logging_values(self):
        """Initialize logging values."""
        self._initialize_energy_values()
    
    def _initialize_output_file(self):
        """Initialize output file with headers."""
        with open(self.output_file_path, 'w') as f:
            f.write('# Energy tracking\n')
            header = '# time(ps)'
            for comp in self.components:
                header += f' {comp}_energy'
            if self.compute_temperature:
                header += ' temperature'
            header += ' total_energy\n'
            f.write(header)
    
    def _compute_kinetic_energy_and_temperature(self, snapshot):
        """Compute kinetic energy and temperature."""
        velocities = snapshot.particles.velocity
        masses = snapshot.particles.mass
        
        # Kinetic energy: KE = 0.5 * m * v^2
        v_squared = np.sum(velocities * velocities, axis=1)
        kinetic_energies = 0.5 * masses * v_squared
        total_kinetic = np.sum(kinetic_energies)
        
        # Temperature: T = 2*KE / (k_B * N_dof)
        num_particles = len(masses)
        if num_particles > 0:
            temperature = (2.0 * total_kinetic) / (3.0 * num_particles * PhysicalConstants.KB_HARTREE_PER_K)
        else:
            temperature = 0.0
            
        return total_kinetic, temperature
    
    def act(self, timestep):
        # Check if output has been stopped due to time limit
        if self.output_stopped:
            return
            
        if timestep == 0:
            return
            
        # Check time limit if specified
        if self.max_time_ps is not None:
            current_time = self._get_current_time(timestep)
            if current_time > self.max_time_ps:
                if not self.output_stopped:
                    self.output_stopped = True
                    print(f"Energy tracking stopped at {current_time:.2f} ps (limit: {self.max_time_ps:.2f} ps)")
                return
        
        # Compute energies
        total_energy = 0.0
        
        with self.sim.state.cpu_local_snapshot as snap:
            for component in self.components:
                if component == "kinetic":
                    kinetic_energy, temperature = self._compute_kinetic_energy_and_temperature(snap)
                    self.current_energies[component] = kinetic_energy
                    if self.compute_temperature:
                        self.current_temperature = temperature
                    total_energy += kinetic_energy
                elif component in ENERGY_COMPONENTS:
                    force_getter = ENERGY_COMPONENTS[component]
                    force_obj = force_getter(self.force_objects)
                    if force_obj is not None:
                        energy = force_obj.energy
                        self.current_energies[component] = energy
                        total_energy += energy
                elif component == "ewald":
                    # Special case: combine short and long range
                    ewald_energy = 0.0
                    for ewald_comp in ["ewald_short", "ewald_long"]:
                        force_obj = ENERGY_COMPONENTS[ewald_comp](self.force_objects)
                        if force_obj is not None:
                            ewald_energy += force_obj.energy
                    self.current_energies[component] = ewald_energy
                    total_energy += ewald_energy
        
        self.current_total_energy = total_energy
        
        # Output periodically
        if self._should_output(timestep):
            current_time = self._get_current_time(timestep)
            
            with open(self.output_file_path, 'a') as f:
                line = f'{current_time:.6f}'
                for comp in self.components:
                    line += f' {self.current_energies[comp]:.6f}'
                if self.compute_temperature:
                    line += f' {self.current_temperature:.6f}'
                line += f' {self.current_total_energy:.6f}\n'
                f.write(line)
            
            self._update_output_step(timestep)
    
    # Logging methods for each energy component
    def _create_energy_logging_method(self, component):
        """Create a logging method for a specific energy component."""
        @hoomd.logging.log
        def energy_method(self):
            return self.current_energies.get(component, 0.0)
        return energy_method
    
    def __init_subclass__(cls, **kwargs):
        """Dynamically create logging methods for energy components."""
        super().__init_subclass__(**kwargs)
        
    @hoomd.logging.log
    def total_energy(self):
        return self.current_total_energy
    
    @hoomd.logging.log
    def temperature(self):
        return self.current_temperature if self.compute_temperature else 0.0


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
    
    @property
    def Dt(self):
        return str(np.round(float(self.simulation.operations.integrator.dt*self.chartime*1000000),6))
    
    @property
    def elapsed(self):
        curtime = datetime.datetime.now()
        return str(curtime-self.starttime)


class ElapsedTimeTracker(hoomd.custom.Action):
    """Track elapsed simulation time in physical units."""
    
    def __init__(self, simulation, runtime):
        super().__init__()
        self.sim = simulation
        self.runtime = runtime
        self.current_elapsed_time = 0.0

    def act(self, timestep):
        if timestep == 0:
            self.current_elapsed_time = 0.0
            return
            
        # Calculate elapsed time in picoseconds
        dt = float(self.sim.operations.integrator.dt)
        self.current_elapsed_time = PhysicalConstants.atomic_units_to_ps(dt * timestep)

    @hoomd.logging.log
    def elapsed_time(self):
        return self.current_elapsed_time


class TimestepFormatter(hoomd.custom.Action):
    """Format timestep information for logging."""
    
    def __init__(self, integrator):
        super().__init__()
        self.integrator = integrator

    def act(self, timestep):
        pass

    @hoomd.logging.log
    def dt_fs(self):
        dt_au = float(self.integrator.dt)
        dt_fs = PhysicalConstants.atomic_units_to_ps(dt_au) * 1000  # Convert ps to fs
        return dt_fs


class CavityModeTracker(hoomd.custom.Action):
    """Track cavity mode properties and energies (specialized tracker)."""
    
    def __init__(self, simulation, cavityforce, time_tracker=None, output_prefix='cavity_mode', output_period=1):
        super().__init__()
        self.sim = simulation
        self.cavityforce = cavityforce
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix
        self.output_period = output_period
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
            f.write('# time(ps) cavity_kinetic_energy cavity_potential_energy cavity_total_energy cavity_temperature\n')

    def compute_cavity_properties(self):
        """Compute cavity mode kinetic and potential energies."""
        try:
            # Get cavity mode information from the force
            cavity_momentum = self.cavityforce.cavity_momentum
            cavity_position = self.cavityforce.cavity_position
            cavity_mass = self.cavityforce.cavity_mass
            cavity_frequency = self.cavityforce.cavity_frequency
            
            # Compute kinetic energy: KE = p^2 / (2*m)
            kinetic_energy = (cavity_momentum * cavity_momentum) / (2.0 * cavity_mass)
            
            # Compute potential energy: PE = 0.5 * k * x^2 = 0.5 * m * omega^2 * x^2
            k_spring = cavity_mass * cavity_frequency * cavity_frequency
            potential_energy = 0.5 * k_spring * cavity_position * cavity_position
            
            # Total energy
            total_energy = kinetic_energy + potential_energy
            
            # Temperature from kinetic energy (1 degree of freedom)
            # T = 2*KE / k_B for 1 DOF
            temperature = (2.0 * kinetic_energy) / PhysicalConstants.KB_HARTREE_PER_K
            
            return kinetic_energy, potential_energy, total_energy, temperature
            
        except AttributeError:
            # Fallback if cavity properties are not accessible
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
        if timestep - self.last_output_step >= self.output_period:
            # Get current time
            if self.time_tracker is not None:
                current_time = self.time_tracker.elapsed_time
            else:
                dt = float(self.sim.operations.integrator.dt)
                current_time = PhysicalConstants.atomic_units_to_ps(dt * timestep)
            
            with open(self.output_file_path, 'a') as f:
                f.write(f'{current_time:.6f} {kinetic:.6f} {potential:.6f} {total:.6f} {temperature:.6f}\n')
            
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