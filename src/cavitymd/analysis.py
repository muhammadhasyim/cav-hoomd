# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Analysis and tracking components for cavity molecular dynamics simulations."""

import hoomd
import datetime
import numpy as np
import traceback
import sys

# CuPy import with fallback for CPU/GPU agnostic code
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

from .utils import PhysicalConstants, unwrap_positions

# =============================================================================
# OBSERVABLE LIBRARY
# =============================================================================

def compute_total_dipole_moment(snapshot):
    """Compute total dipole moment with proper position unwrapping."""
    box_lengths = np.array(
        [snapshot.global_box.L[0], snapshot.global_box.L[1], snapshot.global_box.L[2]]
    )
    
    # Convert snapshot arrays to NumPy following CuPy guidelines
    if HAS_CUPY:
        # Use CuPy's asnumpy function for robust conversion
        positions_np = cp.asnumpy(snapshot.particles.position)
        images_np = cp.asnumpy(snapshot.particles.image)
        charges_np = cp.asnumpy(snapshot.particles.charge)
    else:
        # Fallback to NumPy arrays
        positions_np = np.asarray(snapshot.particles.position)
        images_np = np.asarray(snapshot.particles.image)
        charges_np = np.asarray(snapshot.particles.charge)
    
    unwrapped_positions = unwrap_positions(positions_np, images_np, box_lengths)
    # Dipole = charges × positions
    return np.dot(charges_np, unwrapped_positions)


def compute_density_field(snapshot, wavevectors):
    """Compute density field ρ(k) for given wavevectors."""
    # Convert snapshot arrays to NumPy following CuPy guidelines
    if HAS_CUPY:
        # Use CuPy's asnumpy function for robust conversion
        positions_np = cp.asnumpy(snapshot.particles.position)
    else:
        # Fallback to NumPy arrays
        positions_np = np.asarray(snapshot.particles.position)
    
    rhok_real = np.zeros(len(wavevectors))
    rhok_imag = np.zeros(len(wavevectors))

    for i, k_vec in enumerate(wavevectors):
        # Compute k·r for all particles
        kr = np.dot(positions_np, k_vec)
        # Compute ρ(k) = sum_j exp(i k·r_j)
        rhok_real[i] = np.sum(np.cos(kr))
        rhok_imag[i] = np.sum(np.sin(kr))

    return rhok_real + 1j * rhok_imag


def generate_fibonacci_sphere(samples=100):
    """Generate uniformly distributed points on a sphere using Fibonacci spiral."""
    points = np.zeros((samples, 3))
    phi = np.pi * (3.0 - np.sqrt(5.0))  # golden angle in radians

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
    "bussi_molecular_reservoir": lambda thermostats: thermostats.get(
        "bussi_molecular", None
    ),
    "bussi_cavity_reservoir": lambda thermostats: thermostats.get("bussi_cavity", None),
    "langevin_molecular_reservoir": lambda thermostats: thermostats.get(
        "langevin_molecular", None
    ),
    "langevin_cavity_reservoir": lambda thermostats: thermostats.get(
        "langevin_cavity", None
    ),
}


# =============================================================================
# BASE TRACKER CLASS
# =============================================================================


class BaseTracker(hoomd.custom.Action):
    """Base class for all tracking components with common infrastructure."""

    def __init__(
        self,
        simulation,
        time_tracker=None,
        output_prefix="tracker",
        output_period_steps=None,
        output_period_ps=None,
    ):
        """Initialize base tracker.

        Args:
            simulation: HOOMD simulation object
            time_tracker: Optional time tracker for accurate timing
            output_prefix: Prefix for output files
            output_period_steps: Output frequency in simulation steps (for fixed timestep mode)
            output_period_ps: Output frequency in ps (preferred for adaptive timestep mode)
        """
        super().__init__()
        self.sim = simulation
        self.time_tracker = time_tracker
        self.output_prefix = output_prefix

        # Support both step-based and time-based output periods
        self.output_period_steps = output_period_steps
        self.output_period_ps = output_period_ps

        # Validate that at least one period is specified
        if output_period_steps is None and output_period_ps is None:
            # Default to 1000 steps for backward compatibility
            self.output_period_steps = 1000
            self.output_period_ps = None

        # Prefer time-based periods when both are specified
        if self.output_period_ps is not None:
            self.use_time_based_output = True
        else:
            self.use_time_based_output = False

        # Track last output time/step
        self.last_output_step = 0
        self.last_output_time = 0.0

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
            return PhysicalConstants.atomic_units_to_ps(
                timestep * self.sim.operations.integrator.dt
            )

    def _should_output(self, timestep):
        """Check if we should output at this timestep (supports both time and step-based)."""
        if self.use_time_based_output:
            # Time-based output logic
            current_time = self._get_current_time(timestep)
            time_since_last = current_time - self.last_output_time
            return time_since_last >= self.output_period_ps
        else:
            # Step-based output logic (original behavior)
            return timestep - self.last_output_step >= self.output_period_steps

    def _update_output_step(self, timestep):
        """Update the last output step and time."""
        self.last_output_step = timestep
        if self.use_time_based_output:
            self.last_output_time = self._get_current_time(timestep)

    def get_local_snapshot(self):
        """
        Get the appropriate local snapshot context manager based on simulation device type.
        
        Returns
        -------
        Local snapshot context manager (CPU or GPU)
        """
        # Detect device type from simulation
        device = self.sim.device
        if hasattr(device, 'gpu_ids') or 'GPU' in str(type(device)):
            return self.sim.state.gpu_local_snapshot
        else:
            return self.sim.state.cpu_local_snapshot

    def get_particle_count(self, snap):
        """
        Get the number of particles from a local snapshot in a device-agnostic way.
        
        Parameters
        ----------
        snap : Local snapshot object
            Either CPU or GPU local snapshot
            
        Returns
        -------
        int
            Number of particles
        """
        # Detect device type from simulation
        device = self.sim.device
        if hasattr(device, 'gpu_ids') or 'GPU' in str(type(device)):
            # GPU snapshots don't have .N attribute, use length of an array
            return len(snap.particles.typeid)
        else:
            # CPU snapshots have .N attribute
            return snap.particles.N

    def get_array_module(self, *arrays):
        """
        Get the appropriate array module (numpy or cupy) based on the input arrays.
        
        This follows CuPy's guidelines for writing CPU/GPU agnostic code.
        If any input array is on GPU, returns cupy module, otherwise numpy.
        
        Parameters
        ----------
        *arrays : array-like
            Input arrays to check
            
        Returns
        -------
        module
            Either numpy or cupy module
        """
        if not HAS_CUPY:
            return np
            
        # Use CuPy's get_array_module for proper detection
        return cp.get_array_module(*arrays)

    def to_device_array(self, array_data, target_arrays=None):
        """
        Convert array to the appropriate device format (NumPy for CPU, CuPy for GPU).
        
        This follows CuPy's guidelines using cupy.asarray() for device conversion.
        
        Parameters
        ----------
        array_data : array-like
            Input array data
        target_arrays : array-like, optional
            Target arrays to match device type. If provided, uses their device.
            
        Returns
        -------
        array
            Array on the appropriate device (NumPy for CPU, CuPy for GPU)
        """
        if not HAS_CUPY:
            # Always use NumPy for CPU
            return np.asarray(array_data)
        
        if target_arrays is not None:
            # Use get_array_module to determine the appropriate module
            xp = self.get_array_module(target_arrays)
            return xp.asarray(array_data)
        else:
            # Fallback to cupy.asarray which handles both CPU and GPU cases
            return cp.asarray(array_data)

    def to_numpy(self, array_data):
        """
        Convert array to NumPy array, following CuPy guidelines using cupy.asnumpy().
        
        Parameters
        ----------
        array_data : array-like
            Input array data (NumPy or CuPy)
            
        Returns
        -------
        np.ndarray
            NumPy array
        """
        if HAS_CUPY:
            # Use cupy.asnumpy() for robust conversion
            return cp.asnumpy(array_data)
        else:
            # Fallback to numpy.asarray
            return np.asarray(array_data)


# =============================================================================
# AUTOCORRELATION TRACKER (Simple Observables)
# =============================================================================


class AutocorrelationTracker(BaseTracker):
    """Generic autocorrelation tracker for simple observables."""

    def __init__(
        self,
        simulation,
        observable,
        time_tracker=None,
        output_prefix=None,
        output_period_steps=1000,
        output_period_ps=None,
        reference_interval_steps=10000,
        max_references=10,
        reference_interval_ps=None,
    ):
        """Initialize autocorrelation tracker.

        Args:
            simulation: HOOMD simulation object
            observable: Observable name from SIMPLE_OBSERVABLES
            time_tracker: Optional time tracker for accurate timing
            output_prefix: Prefix for output files
            output_period_steps: Output frequency in simulation steps
            output_period_ps: Output frequency in picoseconds (preferred for adaptive timestep)
            reference_interval_steps: Interval between new references in steps (DEPRECATED in adaptive mode)
            max_references: Maximum number of reference states to keep
            reference_interval_ps: Interval between new references in ps (PREFERRED for adaptive mode)
        """
        if observable not in SIMPLE_OBSERVABLES:
            raise ValueError(
                f"Unknown observable '{observable}'. Available: {list(SIMPLE_OBSERVABLES.keys())}"
            )

        self.observable = observable
        self.observable_func = SIMPLE_OBSERVABLES[observable]
        self.reference_interval_steps = reference_interval_steps  # Fallback for fixed timestep
        self.reference_interval_ps = reference_interval_ps  # Preferred for adaptive timestep
        self.max_references = max_references

        if output_prefix is None:
            output_prefix = f"{observable}_autocorr"

        super().__init__(
            simulation, time_tracker, output_prefix, output_period_steps, output_period_ps
        )

        # Initialize autocorrelation tracking with multiple references
        self.references = []  # List of reference states
        self.last_reference_step = 0
        self.last_reference_time_ps = 0.0  # Track time of last reference for time-based intervals

        # Initialize first reference and output files
        self._initialize_new_reference_file(0)

    def _initialize_new_reference_file(self, ref_number):
        """Initialize a new reference file."""
        filename = f"{self.output_prefix}_{ref_number}.txt"
        with self.get_local_snapshot() as snap:
            reference_value = self.observable_func(snap)
            current_time = self._get_current_time(self.sim.timestep)
            
            # Store reference information
            ref_info = {
                "number": ref_number,
                "filename": filename,
                "value": reference_value,
                "time": current_time
            }
            self.references.append(ref_info)
            
            # Update last reference tracking
            self.last_reference_time_ps = current_time
            
            # Write header and t=0 value for this reference
            with open(filename, "w") as f:
                f.write(f"# {self.observable.capitalize()} autocorrelation data\n")
                f.write(f"# Reference number: {ref_number}\n")
                f.write(f"# Reference time: {current_time:.6f} ps\n")
                
                # Write the correct output period based on which mode is being used
                if self.use_time_based_output:
                    f.write(f"# Output period: {self.output_period_ps:.6f} ps (time-based)\n")
                else:
                    f.write(f"# Output period: {self.output_period_steps} steps (step-based)\n")
                    
                f.write("# timestep lag_time(ps) C(t)\n")
                
                # Compute autocorrelation at t=0 (should be C(0) = |reference|^2)
                autocorr_value = np.dot(reference_value, reference_value)
                f.write(f"{self.sim.timestep} 0.000000 {autocorr_value:.6f}\n")
                f.flush()
        
        print(f"Initialized {self.observable} autocorr reference {ref_number}", flush=True)

    def _should_create_new_reference(self, current_time_ps, timestep):
        """Determine if a new reference should be created based on time or step interval."""
        if len(self.references) >= self.max_references:
            return False

        # Prefer time-based intervals for accuracy in adaptive timestep mode
        if self.reference_interval_ps is not None:
            time_since_last = current_time_ps - self.last_reference_time_ps
            return time_since_last >= self.reference_interval_ps
        else:
            # Fallback to step-based for fixed timestep mode
            step_since_last = timestep - self.last_reference_step
            return step_since_last >= self.reference_interval_steps

    def _initialize_logging_values(self):
        """Initialize logging values."""
        self.current_autocorr_value = 0.0

    def compute_autocorr(self, reference_value, current_value):
        """Compute autocorrelation C(t) = observable(0)·observable(t)."""
        return np.dot(reference_value, current_value)

    def act(self, timestep):
        # Get current time FIRST, outside snapshot context
        current_time = self._get_current_time(timestep)

        if timestep == 0:
            return

        # Compute current observable value
        with self.get_local_snapshot() as snap:
            current_value = self.observable_func(snap)

        # Update autocorrelations for all active references
        for ref in self.references:
            lag_time = current_time - ref["time"]
            autocorr_value = self.compute_autocorr(ref["value"], current_value)

            # Update current autocorr for first reference (for logging)
            if ref["number"] == 0:
                self.current_autocorr_value = autocorr_value

            # Output periodically
            if self._should_output(timestep):
                with open(ref["filename"], "a") as f:
                    f.write(f"{timestep} {lag_time:.6f} {autocorr_value:.6f}\n")
                    f.flush()

        # Add new reference if interval has passed and we haven't hit max
        if self._should_create_new_reference(current_time, timestep):
            ref_number = len(self.references)
            self._initialize_new_reference_file(ref_number)
            self.last_reference_step = timestep

        # Update output step
        if self._should_output(timestep):
            self._update_output_step(timestep)

    @hoomd.logging.log
    def current_autocorr(self):
        return (
            self.current_autocorr_value
            if self.current_autocorr_value is not None
            else 0.0
        )


# =============================================================================
# FIELD AUTOCORRELATION TRACKER (Field Observables)
# =============================================================================


class FieldAutocorrelationTracker(BaseTracker):
    """Generic autocorrelation tracker for field observables with k-space averaging."""

    def __init__(
        self,
        simulation,
        observable,
        time_tracker=None,
        output_prefix=None,
        output_period_steps=1000,
        output_period_ps=None,
        reference_interval_steps=10000,
        max_references=10,
        reference_interval_ps=None,
        **kwargs,
    ):
        """Initialize field autocorrelation tracker.

        Args:
            simulation: HOOMD simulation object
            observable: Observable name from FIELD_OBSERVABLES
            time_tracker: Optional time tracker for accurate timing
            output_prefix: Prefix for output files
            output_period_steps: Output frequency in simulation steps
            reference_interval_steps: Interval between new references in steps (DEPRECATED in adaptive mode)
            max_references: Maximum number of reference states to keep
            reference_interval_ps: Interval between new references in ps (PREFERRED for adaptive mode)
            **kwargs: Observable-specific parameters
        """
        if observable not in FIELD_OBSERVABLES:
            raise ValueError(
                f"Unknown field observable '{observable}'. Available: {list(FIELD_OBSERVABLES.keys())}"
            )

        self.observable = observable
        self.observable_func = FIELD_OBSERVABLES[observable]
        self.reference_interval_steps = (
            reference_interval_steps  # Fallback for fixed timestep
        )
        self.reference_interval_ps = (
            reference_interval_ps  # Preferred for adaptive timestep
        )
        self.max_references = max_references

        if output_prefix is None:
            output_prefix = f"{observable}_field_autocorr"

        super().__init__(
            simulation,
            time_tracker,
            output_prefix,
            output_period_steps,
            output_period_ps,
        )

        # Setup observable-specific parameters
        self._setup_parameters(kwargs)

        # Initialize field autocorrelation tracking
        self.references = []  # List of reference states
        self.last_reference_step = 0
        self.last_reference_time_ps = (
            0.0  # Track time of last reference for time-based intervals
        )

        # Initialize first reference and output files
        self._initialize_new_reference_file(0)

    def _setup_parameters(self, kwargs):
        """Setup parameters specific to the observable."""
        if self.observable == "density_correlation":
            # Default parameters for density correlation
            self.kmag = kwargs.get("kmag", 1.0)
            self.num_wavevectors = kwargs.get("num_wavevectors", 50)
            self.wavevectors = (
                generate_fibonacci_sphere(self.num_wavevectors) * self.kmag
            )
            print(
                f"Density correlation: k={self.kmag:.2f}, {self.num_wavevectors} wavevectors"
            )
        else:
            # Add other observables as needed
            pass

    def _call_observable_func(self, snapshot):
        """Call the observable function with correct arguments."""
        if self.observable == "density_correlation":
            return self.observable_func(snapshot, self.wavevectors)
        else:
            # For other observables, use the function directly
            return self.observable_func(snapshot)

    def _initialize_new_reference_file(self, ref_number):
        """Initialize a new reference file."""
        ref_filename = f"{self.output_prefix}_ref{ref_number}.txt"

        # Initialize reference state
        with self.get_local_snapshot() as snap:
            reference_field = self._call_observable_func(snap)
            current_time = self._get_current_time(self.sim.timestep)

            # Store reference
            self.references.append(
                {
                    "number": ref_number,
                    "filename": ref_filename,
                    "timestep": self.sim.timestep,
                    "time": current_time,
                    "field": reference_field,
                }
            )

            # Update timing trackers
            self.last_reference_time_ps = current_time

            # Compute initial autocorrelation (field with itself at t=0)
            initial_autocorr = self.compute_field_autocorr(
                reference_field, reference_field
            )

            # Create output file with header and immediately write t=0 value
            with open(ref_filename, "w") as f:
                f.write(f"# {self.observable.capitalize()} field autocorrelation\n")
                f.write(f"# Reference {ref_number} at t={current_time:.6f} ps\n")
                if self.use_time_based_output:
                    f.write(f"# Output period: {self.output_period_ps:.3f} ps\n")
                else:
                    f.write(f"# Output period: {self.output_period_steps} steps\n")
                f.write("# timestep lag_time(ps) field_autocorr\n")
                # Write the t=0 correlation value immediately
                f.write(f"{self.sim.timestep} {0.0:.6f} {initial_autocorr:.6f}\n")
                f.flush()

        print(f"Initialized {self.observable} field autocorr reference {ref_number}", flush=True)

    def _initialize_logging_values(self):
        """Initialize logging values."""
        self.current_autocorr_value = 0.0

    def compute_field_autocorr(self, field0, field_t):
        """Compute field autocorrelation."""
        if isinstance(field0, np.ndarray) and isinstance(field_t, np.ndarray):
            return np.mean(np.real(field0 * np.conj(field_t)))
        else:
            return np.dot(field0, field_t)

    def _should_create_new_reference(self, current_time_ps, timestep):
        """Determine if a new reference should be created based on time or step interval."""
        if len(self.references) >= self.max_references:
            return False

        # Prefer time-based intervals for accuracy in adaptive timestep mode
        if self.reference_interval_ps is not None:
            time_since_last = current_time_ps - self.last_reference_time_ps
            return time_since_last >= self.reference_interval_ps
        else:
            # Fallback to step-based for fixed timestep mode
            step_since_last = timestep - self.last_reference_step
            return step_since_last >= self.reference_interval_steps

    def act(self, timestep):
        # Get current time FIRST, outside snapshot context
        current_time = self._get_current_time(timestep)

        if timestep == 0:
            return

        # Compute current field value
        with self.get_local_snapshot() as snap:
            current_field = self._call_observable_func(snap)

        # Update autocorrelations for all active references
        for ref in self.references:
            lag_time = current_time - ref["time"]
            autocorr_value = self.compute_field_autocorr(ref["field"], current_field)

            # Update current autocorr for first reference (for logging)
            if ref["number"] == 0:
                self.current_autocorr_value = autocorr_value

            # Output periodically
            if self._should_output(timestep):
                with open(ref["filename"], "a") as f:
                    f.write(f"{timestep} {lag_time:.6f} {autocorr_value:.6f}\n")
                    f.flush()

        # Add new reference if interval has passed and we haven't hit max
        if self._should_create_new_reference(current_time, timestep):
            ref_number = len(self.references)
            self._initialize_new_reference_file(ref_number)
            self.last_reference_step = timestep

        # Update output step
        if self._should_output(timestep):
            self._update_output_step(timestep)

    @hoomd.logging.log
    def current_autocorr(self):
        return (
            self.current_autocorr_value
            if self.current_autocorr_value is not None
            else 0.0
        )


# =============================================================================
# CORRECTED ENERGY TRACKER - Follows Working EnergyContributationTracker Pattern
# =============================================================================


class EnergyTracker(BaseTracker):
    r"""
    Comprehensive energy tracking for cavity molecular dynamics simulations.
    
    Monitors all energy components in the system to verify energy conservation
    and analyze energy redistribution during cavity coupling experiments. This
    is crucial for validating the physics and understanding energy flow in
    time-varying coupling scenarios.
    
    **Energy Components:**
    
    The total system energy in cavity MD includes:
    
    .. math::
        
        E_{\text{total}} = E_{\text{kinetic}} + E_{\text{potential}} + E_{\text{cavity}} + E_{\text{reservoir}}
    
    where:
    
    - :math:`E_{\text{kinetic}} = \frac{1}{2} \sum_i m_i \vec{v}_i^2` (molecular kinetic energy)
    - :math:`E_{\text{potential}}` includes harmonic bonds, LJ, Coulomb interactions
    - :math:`E_{\text{cavity}} = E_{\text{harmonic}} + E_{\text{coupling}} + E_{\text{dipole}}` (cavity contributions)
    - :math:`E_{\text{reservoir}}` is thermostat reservoir energy (if tracked)
    
    **Cavity Energy Components:**
    
    The cavity energy decomposes as:
    
    .. math::
        
        E_{\text{harmonic}} &= \frac{1}{2} K \tilde{q}_{0,\lambda}^2 \\
        E_{\text{coupling}} &= \tilde{\varepsilon}_{0,\lambda}(t) \tilde{q}_{0,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda} \\
        E_{\text{dipole}} &= \frac{\tilde{\varepsilon}_{0,\lambda}(t)^2}{2K} \left(\sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}\right)^2
    
    **Energy Conservation in Time-Varying Systems:**
    
    During coupling switching, individual components change but total energy is conserved:
    
    .. math::
        
        \frac{dE_{\text{total}}}{dt} = 0 \quad \text{(for conservative dynamics)}
    
    Energy redistribution occurs between molecular and cavity modes according to the new coupling.
    
    Parameters
    ----------
    simulation : hoomd.Simulation
        HOOMD simulation object
    components : List[str]
        Energy components to track. Options include:
        
        - 'kinetic': Molecular kinetic energy
        - 'harmonic': Harmonic bond energy  
        - 'lj': Lennard-Jones potential energy
        - 'ewald_short': Short-range electrostatic energy
        - 'ewald_long': Long-range electrostatic energy
        - 'cavity': Total cavity energy (all components)
        - Individual cavity components: 'cavity_harmonic', 'cavity_coupling', 'cavity_dipole'
    force_objects : Dict[str, hoomd.md.force.Force], optional
        Dictionary mapping component names to force objects for energy extraction
    thermostat_objects : Dict[str, object], optional  
        Dictionary of thermostat objects for reservoir energy tracking
    kinetic_tracker : object, optional
        External kinetic energy tracker (deprecated - use internal computation)
    time_tracker : ElapsedTimeTracker, optional
        Time tracker for accurate timing
    output_prefix : str, optional
        Prefix for output files. Default: 'energy'
    output_period_steps : int, optional
        Output frequency in timesteps (for fixed timestep mode)
    output_period_ps : float, optional
        Output frequency in picoseconds (preferred for adaptive timestep)
    max_timesteps : int, optional
        Maximum timesteps to track (deprecated)
    max_time_ps : float, optional
        Maximum simulation time in ps to track (more accurate than max_timesteps)
    compute_temperature : bool, optional
        Whether to compute temperature from kinetic energy. Default: True
    track_reservoirs : bool, optional  
        Whether to track thermostat reservoir energies. Default: True
    verbose : str, optional
        Verbosity level ('quiet', 'normal', 'verbose'). Default: 'normal'
        
    Attributes
    ----------
    output_file_path : str
        Path to the energy output file
    total_energy : float
        Current total system energy (logged)
    kinetic_energy : float
        Current kinetic energy (logged)
    potential_energy : float
        Current potential energy (logged)
    cavity_energy : float
        Current total cavity energy (logged)
    temperature : float
        Current temperature in Kelvin (logged, if computed)
        
    Examples
    --------
    **Basic energy tracking:**
    
    >>> from hoomd.cavitymd.analysis import EnergyTracker
    >>> 
    >>> # Track key energy components
    >>> energy_tracker = EnergyTracker(
    ...     simulation=sim,
    ...     components=['kinetic', 'harmonic', 'lj', 'cavity'],
    ...     force_objects={'cavity': cavity_force, 'harmonic': harmonic_force},
    ...     output_period_ps=0.1
    ... )
    >>> 
    >>> # Add to simulation
    >>> sim.operations.updaters.append(
    ...     hoomd.update.CustomUpdater(
    ...         action=energy_tracker,
    ...         trigger=hoomd.trigger.Periodic(energy_tracker.output_period_steps)
    ...     )
    ... )
    
    **Detailed cavity energy tracking:**
    
    >>> # Track individual cavity components
    >>> detailed_tracker = EnergyTracker(
    ...     simulation=sim,
    ...     components=['kinetic', 'cavity_harmonic', 'cavity_coupling', 'cavity_dipole'],
    ...     force_objects={'cavity': cavity_force},
    ...     verbose='verbose'
    ... )
    
    **Time-varying coupling validation:**
    
    >>> # Monitor energy conservation during switching
    >>> conservation_tracker = EnergyTracker(
    ...     simulation=sim,
    ...     components=['kinetic', 'harmonic', 'cavity'],
    ...     force_objects={'cavity': cavity_force, 'harmonic': harmonic_force},
    ...     time_tracker=time_tracker,
    ...     max_time_ps=switch_time_ps + 50.0,  # Track through switching
    ...     output_period_ps=0.01  # High frequency during switch
    ... )
    
    Notes
    -----
    - Essential for validating energy conservation in cavity MD simulations
    - Automatically computes kinetic energy internally when requested
    - Supports both detailed component tracking and total energy monitoring
    - Critical for debugging time-varying coupling experiments
    - Output files contain timestep, time, and all tracked energy components
    
    See Also
    --------
    hoomd.cavitymd.forces.CavityForce : For cavity energy components
    ElapsedTimeTracker : For accurate timing
    """

    def __init__(
        self,
        simulation,
        components,
        force_objects=None,
        thermostat_objects=None,
        kinetic_tracker=None,
        time_tracker=None,
        output_prefix="energy",
        output_period_steps=None,
        output_period_ps=None,
        max_timesteps=None,
        max_time_ps=None,
        compute_temperature=True,
        track_reservoirs=True,
        verbose="normal",
    ):
        """Initialize energy tracker with internal kinetic computation capability.

        Args:
            simulation: HOOMD simulation object
            components: List of energy components to track
                       ['kinetic', 'harmonic', 'lj', 'ewald_short', 'ewald_long', 'cavity']
                       NOTE: 'kinetic' is now a standard component!
            force_objects: Dictionary of force objects (harmonic, lj, ewald_short, ewald_long, cavity)
            thermostat_objects: Dictionary of thermostat/method objects for reservoir energy
            kinetic_tracker: [DEPRECATED] KineticEnergyTracker object (use internal computation instead)
                           If None and 'kinetic' in components, will compute kinetic energy internally
            time_tracker: ElapsedTimeTracker for accurate timing
            output_prefix: Prefix for output files
            output_period_steps: Output frequency in simulation steps
            max_timesteps: Maximum timesteps to track (ignored if max_time_ps is set)
            max_time_ps: Maximum simulation time in ps to track (more accurate than max_timesteps)
            compute_temperature: Whether to compute temperature
            track_reservoirs: Whether to track reservoir energies
            verbose: Verbosity level ('quiet', 'normal', 'verbose')
        """
        # Store configuration exactly like working code
        self.force_objects = force_objects or {}
        self.thermostat_objects = thermostat_objects or {}
        self.kinetic_tracker = kinetic_tracker  # Keep for backward compatibility
        self.track_reservoirs = track_reservoirs
        self.max_timesteps = max_timesteps
        self.max_time_ps = max_time_ps
        self.compute_temperature = compute_temperature
        self.output_stopped = False

        # Set verbosity level
        self.verbose = verbose.lower() if isinstance(verbose, str) else "normal"
        if self.verbose not in ["quiet", "normal", "verbose"]:
            self.verbose = "normal"

        # Validate components and check for internal kinetic computation
        self.components = components
        self.use_internal_kinetic = "kinetic" in components and kinetic_tracker is None

        # Warn about deprecated kinetic_tracker usage
        if kinetic_tracker is not None and "kinetic" in components:
            if self.verbose != "quiet":
                print(
                    "WARNING: Both kinetic_tracker and 'kinetic' component specified."
                )
                print("  Using external kinetic_tracker for backward compatibility.")
                print(
                    "  Consider removing kinetic_tracker and using internal computation."
                )

        super().__init__(
            simulation,
            time_tracker,
            output_prefix,
            output_period_steps,
            output_period_ps,
        )

        # Fix file naming to match working code
        self.output_file_path = f"{self.output_prefix}_energy_tracker.txt"

        # Print setup info
        if self.verbose != "quiet":
            print(f"REFACTORED EnergyTracker (Phase 1 - Backward Compatible):", flush=True)
            print(f"  Output file: {self.output_file_path}", flush=True)
            print(f"  Components: {self.components}", flush=True)
            print(f"  Internal kinetic computation: {self.use_internal_kinetic}", flush=True)
            if self.use_internal_kinetic:
                print(
                    "  → Kinetic energy will be computed internally (no external tracker needed)", flush=True
                )
            elif self.kinetic_tracker is not None:
                print("  → Using external kinetic_tracker (deprecated)", flush=True)
            print(f"  Force objects: {list(self.force_objects.keys())}", flush=True)
            print(f"  Thermostat objects: {list(self.thermostat_objects.keys())}", flush=True)
            print(f"  Track reservoirs: {self.track_reservoirs}", flush=True)
            if self.use_time_based_output:
                print(f"  Output period: {self.output_period_ps:.3f} ps", flush=True)
            else:
                print(f"  Output period: {self.output_period_steps} steps", flush=True)
            print(f"  Verbosity: {self.verbose}", flush=True)
            if self.max_time_ps:
                print(f"  Max time: {self.max_time_ps} ps (time-based limit)", flush=True)
            elif self.max_timesteps:
                print(f"  Max timesteps: {self.max_timesteps} (step-based limit)", flush=True)

        # Initialize output file
        self._initialize_output_file()

    def get_local_snapshot(self):
        """
        Get the appropriate local snapshot context manager based on simulation device type.
        
        Returns
        -------
        Local snapshot context manager (CPU or GPU)
        """
        # Detect device type from simulation
        device = self.sim.device
        if hasattr(device, 'gpu_ids') or 'GPU' in str(type(device)):
            return self.sim.state.gpu_local_snapshot
        else:
            return self.sim.state.cpu_local_snapshot

    def get_particle_count(self, snap):
        """
        Get the number of particles from a local snapshot in a device-agnostic way.
        
        Parameters
        ----------
        snap : Local snapshot object
            Either CPU or GPU local snapshot
            
        Returns
        -------
        int
            Number of particles
        """
        # Detect device type from simulation
        device = self.sim.device
        if hasattr(device, 'gpu_ids') or 'GPU' in str(type(device)):
            # GPU snapshots don't have .N attribute, use length of an array
            return len(snap.particles.typeid)
        else:
            # CPU snapshots have .N attribute
            return snap.particles.N

    def get_array_module(self, *arrays):
        """
        Get the appropriate array module (numpy or cupy) based on the input arrays.
        
        This follows CuPy's guidelines for writing CPU/GPU agnostic code.
        If any input array is on GPU, returns cupy module, otherwise numpy.
        
        Parameters
        ----------
        *arrays : array-like
            Input arrays to check
            
        Returns
        -------
        module
            Either numpy or cupy module
        """
        if not HAS_CUPY:
            return np
            
        # Use CuPy's get_array_module for proper detection
        return cp.get_array_module(*arrays)

    def to_device_array(self, array_data, target_arrays=None):
        """
        Convert array to the appropriate device format (NumPy for CPU, CuPy for GPU).
        
        This follows CuPy's guidelines using cupy.asarray() for device conversion.
        
        Parameters
        ----------
        array_data : array-like
            Input array data
        target_arrays : array-like, optional
            Target arrays to match device type. If provided, uses their device.
            
        Returns
        -------
        array
            Array on the appropriate device (NumPy for CPU, CuPy for GPU)
        """
        if not HAS_CUPY:
            # Always use NumPy for CPU
            return np.asarray(array_data)
        
        if target_arrays is not None:
            # Use get_array_module to determine the appropriate module
            xp = self.get_array_module(target_arrays)
            return xp.asarray(array_data)
        else:
            # Fallback to cupy.asarray which handles both CPU and GPU cases
            return cp.asarray(array_data)

    def to_numpy(self, array_data):
        """
        Convert array to NumPy array, following CuPy guidelines using cupy.asnumpy().
        
        Parameters
        ----------
        array_data : array-like
            Input array data (NumPy or CuPy)
            
        Returns
        -------
        np.ndarray
            NumPy array
        """
        if HAS_CUPY:
            # Use cupy.asnumpy() for robust conversion
            return cp.asnumpy(array_data)
        else:
            # Fallback to numpy.asarray
            return np.asarray(array_data)

    def _compute_molecular_kinetic_energy(self):
        """
        Compute molecular kinetic energy and temperature internally.

        Computes directly from simulation state, excluding cavity particles.

        Returns:
            tuple: (kinetic_energy, temperature) in atomic units
        """
        with self.get_local_snapshot() as snap:
            # Convert to numpy array for robust handling across CPU/GPU
            typeid_array = self.to_numpy(snap.particles.typeid)
            # Filter to molecular particles only (exclude cavity particle type 'L')
            molecular_mask = typeid_array != 2  # Type 2 is 'L' (cavity)

            if not np.any(molecular_mask):
                return 0.0, 0.0

            # Convert HOOMD arrays to NumPy first, then apply mask, then convert to device arrays
            velocities_np = self.to_numpy(snap.particles.velocity)[molecular_mask]
            masses_np = self.to_numpy(snap.particles.mass)[molecular_mask]
            
            # Get the appropriate array module for device-agnostic operations
            xp = self.get_array_module(snap.particles.velocity, snap.particles.mass)
            
            # Convert to device-appropriate arrays for operations
            velocities_device = self.to_device_array(velocities_np, snap.particles.velocity)
            masses_device = self.to_device_array(masses_np, snap.particles.mass)

            # Compute kinetic energy using device-appropriate arrays: KE = 0.5 * sum(m_i * v_i^2)
            kinetic_energy = 0.5 * xp.sum(masses_device[:, xp.newaxis] * velocities_device**2)

            # Convert result back to NumPy for final calculations
            kinetic_energy = self.to_numpy(kinetic_energy)

            # Compute temperature: T = (2/3) * KE / (N * k_B)
            N_dof = 3 * len(masses_device)  # 3 degrees of freedom per particle
            temperature = (
                (2.0)
                * kinetic_energy
                / (N_dof * PhysicalConstants.KB_HARTREE_PER_K)
            )

            return float(kinetic_energy), float(temperature)

    def _compute_cavity_kinetic_energy(self):
        """
        Compute cavity kinetic energy internally.

        Computes directly from simulation state.

        Returns:
            float: cavity kinetic energy in atomic units
        """
        # Compute directly from simulation state
        with self.get_local_snapshot() as snap:
            # Convert to numpy array for robust handling across CPU/GPU
            typeid_array = self.to_numpy(snap.particles.typeid)
            # Find cavity particle (type 'L', typeid 2)
            cavity_mask = typeid_array == 2

            if not np.any(cavity_mask):
                return 0.0

            # Convert HOOMD arrays to NumPy first, then apply mask, then convert to device arrays
            cavity_velocity_np = self.to_numpy(snap.particles.velocity)[cavity_mask][0]
            cavity_mass_np = self.to_numpy(snap.particles.mass)[cavity_mask][0]
            
            # Get the appropriate array module for device-agnostic operations
            xp = self.get_array_module(snap.particles.velocity, snap.particles.mass)
            
            # Convert to device-appropriate arrays for operations
            cavity_velocity_device = self.to_device_array(cavity_velocity_np, snap.particles.velocity)
            cavity_mass_device = self.to_device_array(cavity_mass_np, snap.particles.mass)

            # Compute cavity kinetic energy using device-appropriate arrays: KE = 0.5 * m * v^2
            kinetic_energy = 0.5 * cavity_mass_device * xp.sum(cavity_velocity_device**2)

            # Convert result back to NumPy for final calculations
            kinetic_energy = self.to_numpy(kinetic_energy)

            return float(kinetic_energy)

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
            with open(self.output_file_path, "w") as f:
                f.write(
                    "# CORRECTED Energy tracking following working EnergyContributionTracker pattern\n"
                )
                if self.use_time_based_output:
                    f.write(f"# Output period: {self.output_period_ps:.3f} ps\n")
                else:
                    f.write(f"# Output period: {self.output_period_steps} steps\n")
                if self.max_time_ps:
                    f.write(f"# Max time: {self.max_time_ps} ps\n")
                elif self.max_timesteps:
                    f.write(f"# Max timesteps: {self.max_timesteps}\n")
                f.write("# All energies in Hartree (atomic units)\n")
                f.write("# Column definitions:\n")
                f.write("#   time(ps): simulation time in picoseconds\n")
                f.write("#   timestep: simulation timestep number\n")
                f.write("#   harmonic_energy: harmonic bond potential energy\n")
                f.write("#   lj_energy: Lennard-Jones potential energy\n")
                f.write("#   ewald_short_energy: short-range Coulomb energy\n")
                f.write("#   ewald_long_energy: long-range Coulomb energy\n")
                f.write(
                    "#   cavity_harmonic_energy: cavity harmonic potential energy\n"
                )
                f.write("#   cavity_coupling_energy: cavity-molecule coupling energy\n")
                f.write("#   cavity_dipole_self_energy: dipole self-energy\n")
                f.write(
                    "#   cavity_total_potential_energy: total cavity potential energy\n"
                )
                f.write("#   molecular_kinetic_energy: molecular kinetic energy\n")
                f.write("#   cavity_kinetic_energy: cavity kinetic energy\n")
                f.write("#   total_kinetic_energy: total kinetic energy\n")
                f.write("#   total_potential_energy: total potential energy\n")
                f.write("#   system_total_energy: total system energy (KE + PE)\n")
                f.write("#   molecular_reservoir_energy: molecular reservoir energy\n")
                f.write("#   cavity_reservoir_energy: cavity reservoir energy\n")
                f.write("#   total_reservoir_energy: total reservoir energy\n")
                f.write(
                    "#   universe_total_energy: universe total energy (system + reservoir) [CONSERVED]\n"
                )
                if self.compute_temperature:
                    f.write("#   temperature: kinetic temperature (K)\n")

                # Create header line
                header = "time(ps) timestep"
                header += (
                    " harmonic_energy lj_energy ewald_short_energy ewald_long_energy"
                )
                header += " cavity_harmonic_energy cavity_coupling_energy cavity_dipole_self_energy cavity_total_potential_energy"
                header += " molecular_kinetic_energy cavity_kinetic_energy total_kinetic_energy"
                header += " total_potential_energy system_total_energy"
                header += " molecular_reservoir_energy cavity_reservoir_energy total_reservoir_energy"
                header += " universe_total_energy"
                if self.compute_temperature:
                    header += " temperature"
                header += "\n"
                f.write(header)
                f.flush()
            if self.verbose != "quiet":
                print(
                    f"EnergyTracker: Successfully created output file {self.output_file_path}", flush=True
                )
        except Exception as e:
            # Always print errors regardless of verbosity
            print(
                f"EnergyTracker ERROR: Failed to create output file {self.output_file_path}: {e}", flush=True
            )

    def act(self, timestep):
        """Main energy computation method with internal kinetic computation capability."""
        # Check if output has been stopped due to time or timestep limit
        if self.output_stopped:
            return

        if timestep == 0:
            return

        # Get current time for checking time-based limits
        if self.time_tracker is not None:
            current_time = self.time_tracker.elapsed_time
        else:
            current_time = PhysicalConstants.atomic_units_to_ps(
                timestep * self.sim.operations.integrator.dt
            )

        # Check time limit first (more accurate than timestep limit)
        if self.max_time_ps is not None:
            if current_time > self.max_time_ps:
                if not self.output_stopped:
                    self.output_stopped = True
                    if self.verbose != "quiet":
                        print(
                            f"Energy tracking stopped: reached time limit of {self.max_time_ps:.2f} ps at t={current_time:.4f} ps", flush=True
                        )
                return

        # Check timestep limit if specified and no time limit
        elif self.max_timesteps is not None:
            if timestep > self.max_timesteps:
                if not self.output_stopped:
                    self.output_stopped = True
                    if self.verbose != "quiet":
                        print(
                            f"Energy tracking stopped at timestep {timestep} (limit: {self.max_timesteps})", flush=True
                        )
                return

        # Only output periodically (use proper output logic from base class)
        if not self._should_output(timestep):
            return

        try:
            if self.verbose == "verbose":
                print(f"\n=== ENERGY TRACKER DEBUG - Timestep {timestep} ===", flush=True)
                print(f"Current time: {current_time:.6f} ps", flush=True)
                print(f"Internal kinetic computation: {self.use_internal_kinetic}", flush=True)

            # === 1. GET POTENTIAL ENERGY COMPONENTS ===
            if self.verbose == "verbose":
                print("=== POTENTIAL ENERGY COMPONENTS ===", flush=True)

            # Get individual potential energy contributions (direct access like working code)
            try:
                self.current_harmonic_energy = (
                    self.force_objects["harmonic"].energy
                    if "harmonic" in self.force_objects
                    else 0.0
                )
                if self.verbose == "verbose":
                    print(
                        f"Harmonic energy: {self.current_harmonic_energy:.6f} Hartree", flush=True
                    )
            except (AttributeError, KeyError) as e:
                self.current_harmonic_energy = 0.0
                if self.verbose in ["normal", "verbose"]:
                    print(f"Harmonic energy ERROR: {e}", flush=True)

            try:
                self.current_lj_energy = (
                    self.force_objects["lj"].energy
                    if "lj" in self.force_objects
                    else 0.0
                )
                if self.verbose == "verbose":
                    print(f"LJ energy: {self.current_lj_energy:.6f} Hartree", flush=True)
            except (AttributeError, KeyError) as e:
                self.current_lj_energy = 0.0
                if self.verbose in ["normal", "verbose"]:
                    print(f"LJ energy ERROR: {e}", flush=True)

            try:
                self.current_ewald_short_energy = (
                    self.force_objects["ewald_short"].energy
                    if "ewald_short" in self.force_objects
                    else 0.0
                )
                if self.verbose == "verbose":
                    print(
                        f"Ewald short energy: {self.current_ewald_short_energy:.6f} Hartree", flush=True
                    )
            except (AttributeError, KeyError) as e:
                self.current_ewald_short_energy = 0.0
                if self.verbose in ["normal", "verbose"]:
                    print(f"Ewald short energy ERROR: {e}", flush=True)

            try:
                self.current_ewald_long_energy = (
                    self.force_objects["ewald_long"].energy
                    if "ewald_long" in self.force_objects
                    else 0.0
                )
                if self.verbose == "verbose":
                    print(
                        f"Ewald long energy: {self.current_ewald_long_energy:.6f} Hartree", flush=True
                    )
            except (AttributeError, KeyError) as e:
                self.current_ewald_long_energy = 0.0
                if self.verbose in ["normal", "verbose"]:
                    print(f"Ewald long energy ERROR: {e}", flush=True)

            # Calculate total potential energy (without cavity)
            molecular_potential_energy = (
                self.current_harmonic_energy
                + self.current_lj_energy
                + self.current_ewald_short_energy
                + self.current_ewald_long_energy
            )
            if self.verbose == "verbose":
                print(
                    f"Molecular potential energy (harmonic + lj + ewald): {molecular_potential_energy:.6f} Hartree", flush=True
                )

            # Get cavity potential energy components if present
            self.current_cavity_harmonic_energy = 0.0
            self.current_cavity_coupling_energy = 0.0
            self.current_cavity_dipole_self_energy = 0.0
            self.current_cavity_total_potential_energy = 0.0

            if (
                "cavity" in self.force_objects
                and self.force_objects["cavity"] is not None
            ):
                cavityforce = self.force_objects["cavity"]
                if self.verbose == "verbose":
                    print("=== CAVITY ENERGY COMPONENTS ===", flush=True)
                try:
                    # Use the logged property methods directly instead of getattr
                    self.current_cavity_harmonic_energy = cavityforce.harmonic_energy
                    self.current_cavity_coupling_energy = cavityforce.coupling_energy
                    self.current_cavity_dipole_self_energy = cavityforce.dipole_self_energy

                    if self.verbose == "verbose":
                        print(
                            f"Cavity harmonic energy: {self.current_cavity_harmonic_energy:.6f} Hartree", flush=True
                        )
                        print(
                            f"Cavity coupling energy: {self.current_cavity_coupling_energy:.6f} Hartree", flush=True
                        )
                        print(
                            f"Cavity dipole self energy: {self.current_cavity_dipole_self_energy:.6f} Hartree", flush=True
                        )

                    # For total energy, try .energy property first, then sum components
                    if hasattr(cavityforce, "energy"):
                        self.current_cavity_total_potential_energy = cavityforce.energy
                        if self.verbose == "verbose":
                            print(
                                f"Cavity total energy (from .energy): {self.current_cavity_total_potential_energy:.6f} Hartree", flush=True
                            )
                    else:
                        self.current_cavity_total_potential_energy = (
                            self.current_cavity_harmonic_energy
                            + self.current_cavity_coupling_energy
                            + self.current_cavity_dipole_self_energy
                        )
                        if self.verbose == "verbose":
                            print(
                                f"Cavity total energy (sum components): {self.current_cavity_total_potential_energy:.6f} Hartree", flush=True
                            )
                except Exception as e:
                    self.current_cavity_harmonic_energy = 0.0
                    self.current_cavity_coupling_energy = 0.0
                    self.current_cavity_dipole_self_energy = 0.0
                    self.current_cavity_total_potential_energy = 0.0
                    if self.verbose in ["normal", "verbose"]:
                        print(f"ERROR accessing cavity energy components: {e}", flush=True)
                        print(f"  This might indicate a timing issue with the logged properties", flush=True)
                        print(f"  Cavity force implementation: {getattr(cavityforce, 'implementation', 'unknown')}", flush=True)
                        print(f"  Cavity force object: {type(cavityforce)}", flush=True)
            else:
                if self.verbose == "verbose":
                    print("No cavity force object - cavity energies set to zero", flush=True)

            # Calculate total potential energy
            self.current_total_potential_energy = (
                molecular_potential_energy + self.current_cavity_total_potential_energy
            )
            if self.verbose == "verbose":
                print(
                    f"TOTAL POTENTIAL ENERGY: {self.current_total_potential_energy:.6f} Hartree", flush=True
                )

            # === 2. GET KINETIC ENERGY COMPONENTS (NEW: Internal vs External) ===
            if self.verbose == "verbose":
                print("=== KINETIC ENERGY COMPONENTS ===", flush=True)

            # Get molecular kinetic energy - either internal or external
            self.current_molecular_kinetic_energy = 0.0

            if self.use_internal_kinetic:
                # NEW: Compute kinetic energy internally
                if self.verbose == "verbose":
                    print("Using INTERNAL kinetic energy computation", flush=True)

                if "kinetic" in self.components:
                    molecular_ke, molecular_temp = (
                        self._compute_molecular_kinetic_energy()
                    )
                    self.current_molecular_kinetic_energy = molecular_ke

                    if self.verbose == "verbose":
                        print(
                            f"Molecular kinetic energy (internal): {molecular_ke:.6f} Hartree", flush=True
                        )
                        print(
                            f"Molecular temperature (internal): {molecular_temp:.2f} K", flush=True
                        )

                    # Store temperature for later use
                    self._internal_molecular_temperature = molecular_temp

            else:
                # BACKWARD COMPATIBILITY: Use external kinetic tracker
                if self.verbose == "verbose":
                    print("Using EXTERNAL kinetic energy tracker (deprecated)", flush=True)

                if self.kinetic_tracker is not None:
                    try:
                        self.current_molecular_kinetic_energy = (
                            self.kinetic_tracker.kinetic_energy
                        )
                        if self.verbose == "verbose":
                            print(
                                f"Molecular kinetic energy (external tracker): {self.current_molecular_kinetic_energy:.6f} Hartree", flush=True
                            )
                    except AttributeError as e:
                        self.current_molecular_kinetic_energy = 0.0
                        if self.verbose in ["normal", "verbose"]:
                            print(f"Molecular kinetic energy ERROR: {e}", flush=True)
                else:
                    if self.verbose == "verbose":
                        print(
                            "No kinetic tracker - molecular kinetic energy set to zero", flush=True
                        )

            # Get cavity kinetic energy
            self.current_cavity_kinetic_energy = 0.0

            if self.use_internal_kinetic and "kinetic" in self.components:
                # NEW: Compute cavity kinetic energy internally
                cavity_ke = self._compute_cavity_kinetic_energy()
                self.current_cavity_kinetic_energy = cavity_ke
                if self.verbose == "verbose":
                    print(f"Cavity kinetic energy (internal): {cavity_ke:.6f} Hartree", flush=True)
            else:
                if self.verbose == "verbose":
                    print("No cavity kinetic energy computation", flush=True)

            # Calculate total kinetic energy
            self.current_total_kinetic_energy = (
                self.current_molecular_kinetic_energy
                + self.current_cavity_kinetic_energy
            )
            if self.verbose == "verbose":
                print(
                    f"TOTAL KINETIC ENERGY: {self.current_total_kinetic_energy:.6f} Hartree", flush=True
                )

            # === 3. GET RESERVOIR ENERGIES ===
            if self.verbose == "verbose":
                print("=== RESERVOIR ENERGY COMPONENTS ===", flush=True)

            # Get molecular reservoir energy if available
            molecular_reservoir_energy = 0.0

            # Check for molecular Langevin method
            if "langevin_molecular" in self.thermostat_objects:
                try:
                    mol_langevin_reservoir = self.thermostat_objects[
                        "langevin_molecular"
                    ].reservoir_energy
                    molecular_reservoir_energy += mol_langevin_reservoir
                    if self.verbose == "verbose":
                        print(
                            f"Molecular Langevin reservoir energy: {mol_langevin_reservoir:.6f} Hartree", flush=True
                        )
                except AttributeError:
                    if self.verbose == "verbose":
                        print("Molecular Langevin reservoir energy not available yet", flush=True)

            # Check for molecular Bussi thermostat
            if "bussi_molecular" in self.thermostat_objects:
                try:
                    mol_bussi_reservoir = self.thermostat_objects[
                        "bussi_molecular"
                    ].total_reservoir_energy
                    molecular_reservoir_energy += mol_bussi_reservoir
                    if self.verbose == "verbose":
                        print(
                            f"Molecular Bussi reservoir energy: {mol_bussi_reservoir:.6f} Hartree", flush=True
                        )
                except (AttributeError, hoomd.error.DataAccessError):
                    if self.verbose == "verbose":
                        print("Molecular Bussi reservoir energy not available yet", flush=True)

            self.current_molecular_reservoir_energy = molecular_reservoir_energy

            # Get cavity reservoir energy if available
            cavity_reservoir_energy = 0.0

            # Check for cavity Langevin method
            if "langevin_cavity" in self.thermostat_objects:
                try:
                    cav_langevin_reservoir = self.thermostat_objects[
                        "langevin_cavity"
                    ].reservoir_energy
                    cavity_reservoir_energy += cav_langevin_reservoir
                    if self.verbose == "verbose":
                        print(
                            f"Cavity Langevin reservoir energy: {cav_langevin_reservoir:.6f} Hartree", flush=True
                        )
                except AttributeError:
                    if self.verbose == "verbose":
                        print("Cavity Langevin reservoir energy not available yet", flush=True)

            # Check for cavity Bussi thermostat
            if "bussi_cavity" in self.thermostat_objects:
                try:
                    cav_bussi_reservoir = self.thermostat_objects[
                        "bussi_cavity"
                    ].total_reservoir_energy
                    cavity_reservoir_energy += cav_bussi_reservoir
                    if self.verbose == "verbose":
                        print(
                            f"Cavity Bussi reservoir energy: {cav_bussi_reservoir:.6f} Hartree", flush=True
                        )
                except (AttributeError, hoomd.error.DataAccessError):
                    if self.verbose == "verbose":
                        print("Cavity Bussi reservoir energy not available yet", flush=True)

            self.current_cavity_reservoir_energy = cavity_reservoir_energy

            # Calculate total reservoir energy
            self.current_total_reservoir_energy = (
                self.current_molecular_reservoir_energy
                + self.current_cavity_reservoir_energy
            )
            if self.verbose == "verbose":
                print(
                    f"TOTAL RESERVOIR ENERGY: {self.current_total_reservoir_energy:.6f} Hartree", flush=True
                )

            # === 4. CALCULATE TOTAL ENERGIES ===
            if self.verbose == "verbose":
                print("=== TOTAL ENERGY CALCULATIONS ===", flush=True)

            # Calculate system total energy
            self.current_system_total_energy = (
                self.current_total_potential_energy + self.current_total_kinetic_energy
            )
            if self.verbose == "verbose":
                print(
                    f"SYSTEM TOTAL ENERGY (KE + PE): {self.current_system_total_energy:.6f} Hartree", flush=True
                )

            # Calculate universe total energy
            self.current_universe_total_energy = (
                self.current_system_total_energy + self.current_total_reservoir_energy
            )
            if self.verbose == "verbose":
                print(
                    f"UNIVERSE TOTAL ENERGY (system + reservoir): {self.current_universe_total_energy:.6f} Hartree", flush=True
                )
                print(f"  (This should be conserved)", flush=True)

            # Calculate temperature if requested
            if self.compute_temperature:
                if self.use_internal_kinetic and hasattr(
                    self, "_internal_molecular_temperature"
                ):
                    # Use temperature from internal computation
                    self.current_temperature = self._internal_molecular_temperature
                    if self.verbose == "verbose":
                        print(
                            f"Temperature (from internal computation): {self.current_temperature:.2f} K", flush=True
                        )
                elif self.kinetic_tracker is not None:
                    # Use temperature from external tracker
                    try:
                        self.current_temperature = self.kinetic_tracker.temperature
                        if self.verbose == "verbose":
                            print(
                                f"Temperature (from external tracker): {self.current_temperature:.2f} K", flush=True
                            )
                    except AttributeError:
                        self.current_temperature = 0.0
                        if self.verbose == "verbose":
                            print("Temperature not available from external tracker", flush=True)
                else:
                    self.current_temperature = 0.0
                    if self.verbose == "verbose":
                        print("No kinetic computation - temperature set to zero", flush=True)

            # === 5. WRITE OUTPUT DATA ===
            if self.verbose == "verbose":
                print("=== WRITING OUTPUT DATA ===", flush=True)
            self._write_energy_data(timestep, current_time)
            self._update_output_step(timestep)

            if self.verbose == "verbose":
                print(f"=== END ENERGY TRACKER DEBUG - Timestep {timestep} ===\n", flush=True)

        except Exception as e:
            # Always print critical errors regardless of verbosity
            print(f"EnergyTracker CRITICAL ERROR at timestep {timestep}: {e}", flush=True)
            import traceback

            traceback.print_exc()

    def _write_energy_data(self, timestep, current_time):
        """Write energy data to output file."""
        try:
            # Build output line exactly like working code format
            output_values = [
                current_time,
                timestep,
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
                self.current_universe_total_energy,
            ]

            if self.compute_temperature:
                output_values.append(self.current_temperature)

            # Write to file
            with open(self.output_file_path, "a") as f:
                formatted_values = [
                    f"{val:.6f}" if isinstance(val, float) else str(val)
                    for val in output_values
                ]
                f.write(" ".join(formatted_values) + "\n")
                f.flush()

            if self.verbose == "verbose":
                print(f"Successfully wrote energy data to {self.output_file_path}", flush=True)

        except Exception as e:
            # Always print errors regardless of verbosity
            print(f"EnergyTracker ERROR writing data at timestep {timestep}: {e}", flush=True)
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
            simulation_time_ps = PhysicalConstants.atomic_units_to_ps(
                dt * current_timestep
            )

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
            simulation_time_ps = PhysicalConstants.atomic_units_to_ps(
                dt * current_timestep
            )

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
        return str(
            np.round(
                float(
                    self.simulation.operations.integrator.dt * self.chartime * 1000000
                ),
                6,
            )
        )

    @property
    def elapsed(self):
        curtime = datetime.datetime.now()
        return str(curtime - self.starttime)


class ElapsedTimeTracker(hoomd.custom.Action):
    r"""
    Track elapsed simulation time in physical units and exit when runtime is reached.
    
    This class provides accurate time tracking for adaptive timestep simulations,
    where the timestep size :math:`\Delta t` may vary during the simulation. It
    accumulates the actual elapsed time by integrating timestep increments:
    
    .. math::
        
        t_{\text{elapsed}} = \sum_{i=1}^{n} \Delta t_i
    
    where :math:`\Delta t_i` is the timestep size at step :math:`i`.
    
    **Adaptive Timestep Compatibility:**
    
    Unlike simple step counting, this tracker accounts for variable timesteps
    in adaptive integration schemes, ensuring accurate timing for:
    
    - Time-varying coupling experiments
    - Runtime termination conditions  
    - Performance metrics (ns/day calculations)
    - Output period controls
    
    **Integration with Time-Varying Parameters:**
    
    This tracker is essential for StepVariant and other time-dependent parameters
    that need to know the actual simulation time, not just the timestep number.
    
    Parameters
    ----------
    simulation : hoomd.Simulation
        HOOMD simulation object to track
    runtime : float
        Target runtime in picoseconds. Simulation exits when this time is reached.
        
    Attributes
    ----------
    runtime : float
        Target runtime in picoseconds
    total_time : float
        Current elapsed time in atomic units
    initial_timestep : int
        Starting timestep number (for inherited simulations)
    last_timestep : int
        Last processed timestep number
        
    Methods
    -------
    elapsed_time : float
        Current elapsed time in picoseconds (logged property)
    act(timestep) : None
        Update elapsed time and check for runtime completion
        
    Examples
    --------
    **Basic usage:**
    
    >>> from hoomd.cavitymd.analysis import ElapsedTimeTracker
    >>> 
    >>> # Track a 1000 ps simulation
    >>> time_tracker = ElapsedTimeTracker(sim, runtime=1000.0)
    >>> 
    >>> # Add to simulation updaters
    >>> sim.operations.updaters.append(
    ...     hoomd.update.CustomUpdater(
    ...         action=time_tracker, 
    ...         trigger=hoomd.trigger.Periodic(1)  # Update every step
    ...     )
    ... )
    
    **With time-varying coupling:**
    
    >>> # Use for accurate timing in step variants
    >>> time_tracker = ElapsedTimeTracker(sim, runtime=500.0)
    >>> coupling_variant = StepVariant(
    ...     target_value=0.001,
    ...     switch_time_ps=100.0,
    ...     time_tracker=time_tracker
    ... )
    
    Notes
    -----
    - Must be updated every timestep for accurate timing
    - Automatically exits simulation when runtime is reached
    - Handles inherited timesteps from restart simulations
    - Essential for time-varying parameter variants
    - Works correctly with both fixed and adaptive timestep schemes
    
    See Also
    --------
    hoomd.cavitymd.variants.StepVariant : For time-varying parameters
    PerformanceTracker : For simulation performance monitoring
    """
    
    def __init__(self, simulation: hoomd.Simulation, runtime: float) -> None:
        """
        Initialize the elapsed time tracker.
        
        Parameters
        ----------
        simulation : hoomd.Simulation
            HOOMD simulation object to track
        runtime : float
            Target runtime in picoseconds
        """
        super().__init__()
        self.simulation: hoomd.Simulation = simulation
        self.runtime: float = runtime
        
        # Initialize timing variables
        self.total_time: float = 0.0  # Total elapsed time in atomic units
        self.initial_timestep: int = 0  # Starting timestep (for inherited sims)
        self.last_timestep: int = 0  # Last processed timestep
        
        print(f"ElapsedTimeTracker initialized: target runtime = {runtime:.1f} ps", flush=True)

    def act(self, timestep: int) -> None:
        """
        Update the total elapsed time by accumulating time increments.
        
        This method is called by HOOMD at each timestep to update the elapsed time.
        For adaptive timestep simulations, it properly accounts for varying timestep sizes.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep number
        """
        # Get current timestep size
        dt = self.simulation.operations.integrator.dt

        # For the first call, handle initialization
        if self.last_timestep == 0:
            # Initialize - record the starting timestep but don't add its time
            self.initial_timestep = timestep
            self.last_timestep = timestep
            self.total_time = 0.0  # Always start elapsed time from 0, regardless of inherited timestep
            if timestep > 0:
                print(f"NOTICE: Starting from inherited timestep {timestep}", flush=True)
                print(
                    f"  Elapsed time will start from 0, not from inherited simulation time", flush=True
                )
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
            print(f"Runtime {self.runtime} ps reached. Exiting simulation.", flush=True)
            import sys
            sys.exit(0)

    @hoomd.logging.log
    def elapsed_time(self) -> float:
        """
        Current elapsed time in picoseconds.
        
        Returns
        -------
        float
            Elapsed simulation time in picoseconds
            
        Notes
        -----
        This property is logged and can be accessed by other components
        that need accurate timing information.
        """
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


# =============================================================================
# CONVENIENCE ALIASES AND BACKWARD COMPATIBILITY
# =============================================================================


class DipoleAutocorrelation(AutocorrelationTracker):
    """Dipole autocorrelation tracker - convenience wrapper around AutocorrelationTracker.

    This class provides backward compatibility for the DipoleAutocorrelation class
    that was previously defined in the monolithic cavitymd.py file.
    """

    def __init__(
        self,
        simulation,
        time_tracker=None,
        output_prefix="dipole_autocorr",
        output_period_steps=1000,
    ):
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
            output_period_steps=output_period_steps,
        )


# =============================================================================
# PERFORMANCE TRACKER
# =============================================================================


class PerformanceTracker(hoomd.custom.Action):
    """Custom performance tracker to display ns/day and other metrics.
    
    This tracker provides:
    - Nanoseconds per day (ns/day) performance metric
    - Estimated time to completion (ETA)
    - Dynamic precision for small ns/day values
    """
    
    def __init__(self, simulation, runtime_ps, time_tracker=None):
        """Initialize performance tracker.
        
        Args:
            simulation: HOOMD simulation object
            runtime_ps: Total runtime in picoseconds
            time_tracker: Optional time tracker for accurate timing
        """
        super().__init__()
        self.sim = simulation
        self.runtime_ps = runtime_ps
        self.time_tracker = time_tracker
        self.start_time = datetime.datetime.now().timestamp()
        self.last_timestep = 0
        self.current_ns_per_day = 0.0
        self.current_eta = ""
        
    def act(self, timestep):
        """Update performance metrics at each timestep."""
        # Calculate simulation progress
        if self.time_tracker is not None:
            simulation_time_ps = self.time_tracker.elapsed_time
        else:
            dt = float(self.sim.operations.integrator.dt)
            simulation_time_ps = PhysicalConstants.atomic_units_to_ps(dt * timestep)
        
        # Calculate wall time elapsed
        wall_time_elapsed = datetime.datetime.now().timestamp() - self.start_time
        
        # Calculate performance metrics (even if potentially unreliable early on)
        if wall_time_elapsed > 0 and simulation_time_ps > 0:
            # Calculate ns/day
            ps_per_second = simulation_time_ps / wall_time_elapsed
            ns_per_second = ps_per_second / 1000.0
            self.current_ns_per_day = ns_per_second * 86400
            
            # Calculate ETA
            total_wall_time_needed = (self.runtime_ps / simulation_time_ps) * wall_time_elapsed
            seconds_remaining = total_wall_time_needed - wall_time_elapsed
            
            if seconds_remaining > 0:
                eta_td = datetime.timedelta(seconds=int(seconds_remaining))
                self.current_eta = str(eta_td)
            else:
                self.current_eta = "0:00:00"

    @hoomd.logging.log
    def ns_per_day(self):
        """Return ns/day performance metric with appropriate precision."""
        # Use more decimal places for small values
        if self.current_ns_per_day < 0.01:
            return f"{self.current_ns_per_day:.4f}"
        else:
            return f"{self.current_ns_per_day:.2f}"
    
    @hoomd.logging.log
    def eta_remaining(self):
        """Return estimated time to completion."""
        return self.current_eta
