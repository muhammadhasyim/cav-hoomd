
# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Analysis and tracking tools for cavity molecular dynamics simulations."""

from typing import Optional, List, Dict, Union, Any, Tuple
import hoomd
import numpy as np
import logging
import sys
import os
import time
from pathlib import Path
import enum
import scipy.fft

# CuPy import with fallback for CPU/GPU agnostic code
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    cp = None
    HAS_CUPY = False

from .utils import PhysicalConstants, unwrap_positions


class Status:
    """Simulation status tracker that provides real-time status information."""
    
    def __init__(self, sim, runtime_ps, time_tracker):
        """Initialize status tracker.
        
        Parameters
        ----------
        sim : hoomd.Simulation
            HOOMD simulation object
        runtime_ps : float
            Total runtime in picoseconds
        time_tracker : ElapsedTimeTracker
            Time tracking object
        """
        self.sim = sim
        self.runtime_ps = runtime_ps
        self.time_tracker = time_tracker
        self._status = "initializing"
    
    @property 
    def status(self):
        """Current simulation status."""
        return self._status
    
    @status.setter
    def status(self, value):
        """Set simulation status."""
        self._status = value
        
    def update(self):
        """Update status based on current simulation state."""
        if hasattr(self.time_tracker, 'total_time'):
            elapsed_ps = PhysicalConstants.atomic_units_to_ps(self.time_tracker.total_time)
            if elapsed_ps >= self.runtime_ps:
                self._status = "completed"
            elif elapsed_ps > 0:
                self._status = "running"
        
    def get_progress(self):
        """Get current progress as a percentage."""
        if hasattr(self.time_tracker, 'total_time'):
            elapsed_ps = PhysicalConstants.atomic_units_to_ps(self.time_tracker.total_time)
            return min(100.0, (elapsed_ps / self.runtime_ps) * 100.0)
        return 0.0
    

class TimestepFormatter:
    """Utility class for formatting timestep-related output."""
    
    def __init__(self, integrator=None):
        """Initialize timestep formatter.
        
        Parameters
        ----------
        integrator : hoomd.md.Integrator, optional
            HOOMD integrator object for timestep information
        """
        self.integrator = integrator
    
    @staticmethod
    def format_timestep(timestep: int) -> str:
        """Format timestep with appropriate scale."""
        if timestep >= 1e6:
            return f"{timestep/1e6:.1f}M"
        elif timestep >= 1e3:
            return f"{timestep/1e3:.1f}k"
        else:
            return str(timestep)
    
    @staticmethod
    def format_time(time_ps: float) -> str:
        """Format time in appropriate units."""
        if time_ps >= 1000:
            return f"{time_ps/1000:.1f} ns"
        else:
            return f"{time_ps:.1f} ps"
            
    def get_current_timestep_info(self):
        """Get current timestep information from integrator."""
        if self.integrator is not None:
            return {
                'dt': self.integrator.dt,
                'dt_ps': PhysicalConstants.atomic_units_to_ps(self.integrator.dt)
            }
        return None
    
    # Note: dt_fs method temporarily removed due to HOOMD logging metaclass issues
    # TODO: Re-add when logging issue is resolved


class ElapsedTimeTracker(hoomd.custom.Action):
    """Tracks the total elapsed time in a simulation with variable timesteps."""
    
    def __init__(self, simulation, runtime):
        super().__init__()
        self.simulation = simulation
        self.total_time = 0.0
        self.runtime = runtime
        self.last_timestep = 0  # Start from 0, not simulation.timestep
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

    @hoomd.logging.log(category="scalar")
    def elapsed_time(self):
        """Return elapsed time in picoseconds."""
        return PhysicalConstants.atomic_units_to_ps(self.total_time)


class FieldAutocorrelationTracker(hoomd.custom.Action):
    """
    Tracks field autocorrelation functions (e.g., F(k,t) density correlations) during simulation.
    
    Supports both regular time spacing and logarithmic time spacing for F(k,t) output.
    """
    
    def __init__(self, 
        simulation,
                 observable: str,
                 time_tracker,
                 output_period_ps: float,
                 output_prefix: str,
                 reference_interval_ps: float,
                 max_references: int,
                 kmag: float = 1.0,
                 num_wavevectors: int = 50,
                 # New logarithmic time spacing parameters
                 log_time_spacing: bool = False,
                 min_log_time_ps: Optional[float] = None,
                 max_log_time_ps: Optional[float] = None,
                 log_num_points: int = 50):
        """
        Initialize FieldAutocorrelationTracker.
        
        Parameters
        ----------
        simulation : hoomd.Simulation
            The HOOMD simulation object
        observable : str
            Type of observable to track (e.g., "density_correlation")
        time_tracker : ElapsedTimeTracker
            Time tracker for accurate timing
        output_period_ps : float
            Output period in picoseconds
        output_prefix : str
            Prefix for output files
        reference_interval_ps : float
            Interval between reference frames in picoseconds
        max_references : int
            Maximum number of reference frames to keep
        kmag : float, optional
            k-vector magnitude for density correlations. Default: 1.0
        num_wavevectors : int, optional
            Number of wavevectors for k-averaging. Default: 50
        log_time_spacing : bool, optional
            Whether to use logarithmic time spacing for output. Default: False
        min_log_time_ps : float, optional
            Minimum time for log spacing. Required if log_time_spacing=True
        max_log_time_ps : float, optional
            Maximum time for log spacing. Required if log_time_spacing=True  
        log_num_points : int, optional
            Number of points in log spacing. Default: 50
        """
        super().__init__()
        self.simulation = simulation
        self.observable = observable
        self.time_tracker = time_tracker
        self.output_period_ps = output_period_ps
        self.output_prefix = output_prefix
        self.reference_interval_ps = reference_interval_ps
        self.max_references = max_references
        self.kmag = kmag
        self.num_wavevectors = num_wavevectors
        
        # Logarithmic time spacing parameters
        self.log_time_spacing = log_time_spacing
        self.min_log_time_ps = min_log_time_ps
        self.max_log_time_ps = max_log_time_ps
        self.log_num_points = log_num_points
        
        # Validate log time parameters
        if self.log_time_spacing:
            if min_log_time_ps is None or max_log_time_ps is None:
                raise ValueError("min_log_time_ps and max_log_time_ps must be specified when log_time_spacing=True")
            if min_log_time_ps >= max_log_time_ps:
                raise ValueError("min_log_time_ps must be less than max_log_time_ps")
            if min_log_time_ps <= 0:
                raise ValueError("min_log_time_ps must be positive for logarithmic spacing")
                
            # Generate logarithmic time points
            self.log_time_points_ps = np.logspace(
                np.log10(min_log_time_ps), 
                np.log10(max_log_time_ps), 
                log_num_points
            )
            print(f"F(k,t) logarithmic time spacing enabled:")
            print(f"  Time range: {min_log_time_ps:.3f} - {max_log_time_ps:.1f} ps")
            print(f"  Number of points: {log_num_points}")
            print(f"  Sample times: {self.log_time_points_ps[:5]}, ..., {self.log_time_points_ps[-5:]}")
        
        # Generate wavevectors
        if observable == "density_correlation":
            self.wavevectors = self._fibonacci_sphere(samples=num_wavevectors) * kmag
        
        # Reference frames storage
        self.references = []
        self.last_reference_time = None
        self.last_output_time = None
        
        # Storage for correlation data (for logarithmic spacing)
        self.correlation_data = {}  # ref_idx -> {'times': [], 'values': []}
        
        # For caching computed values
        self.current_computed_value = None
        self.current_timestep = -1
        
        print(f"FieldAutocorrelationTracker initialized for {observable}")
        print(f"  k-magnitude: {kmag}")
        print(f"  Wavevectors: {num_wavevectors}")
        print(f"  Reference interval: {reference_interval_ps} ps")
        print(f"  Max references: {max_references}")
        if self.log_time_spacing:
            print(f"  Log time spacing: {log_num_points} points from {min_log_time_ps} to {max_log_time_ps} ps")
        
    def _fibonacci_sphere(self, samples=100):
        """Generate points on a sphere using Fibonacci spiral."""
        points = []
        phi = np.pi * (3. - np.sqrt(5.))  # Golden angle in radians
        
        for i in range(samples):
            y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
            radius = np.sqrt(1 - y * y)  # radius at y
            
            theta = phi * i  # golden angle increment
            
            x = np.cos(theta) * radius
            z = np.sin(theta) * radius
            
            points.append([x, y, z])
        
        return np.array(points)
    
    def _should_add_reference(self, current_time_ps: float) -> bool:
        """Check if we should add a new reference frame."""
        if self.last_reference_time is None:
            return True
        return (current_time_ps - self.last_reference_time) >= self.reference_interval_ps
    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output correlation data."""
        if self.last_output_time is None:
            return True
        return (current_time_ps - self.last_output_time) >= self.output_period_ps
    
    def _compute_density_field(self, positions: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Compute density field ρ_k for given positions."""
        # Determine array module (NumPy or CuPy)
        if HAS_CUPY and hasattr(positions, 'device'):
            xp = cp
        else:
            xp = np
            
        # Convert to appropriate array type
        if xp == cp and not hasattr(positions, 'device'):
            positions = cp.asarray(positions)
        elif xp == np and hasattr(positions, 'device'):
            positions = cp.asnumpy(positions)
        
        # Compute ρ_k = Σ_j exp(i k · r_j)
        # k_dot_r shape: (num_wavevectors, num_particles)
        k_dot_r = xp.dot(self.wavevectors, positions.T)
        
        # Compute complex exponentials
        rhok_complex = xp.sum(xp.exp(1j * k_dot_r), axis=1)
        
        # Split into real and imaginary parts
        rhok_real = xp.real(rhok_complex)
        rhok_imag = xp.imag(rhok_complex)
        
        # Convert back to NumPy if needed
        if xp == cp:
            rhok_real = cp.asnumpy(rhok_real)
            rhok_imag = cp.asnumpy(rhok_imag)
        
        return rhok_real, rhok_imag
    
    def _compute_correlation(self, rhok0_real: np.ndarray, rhok0_imag: np.ndarray,
                           rhok_real: np.ndarray, rhok_imag: np.ndarray) -> float:
        """Compute F(k,t) = <ρ_k(t) · ρ_k*(t_0)> averaged over k-vectors."""
        # ρ_k(t) · ρ_k*(t_0) = [real(t) + i·imag(t)] · [real(t_0) - i·imag(t_0)]
        # = real(t)·real(t_0) + imag(t)·imag(t_0) + i·[imag(t)·real(t_0) - real(t)·imag(t_0)]
        # Take real part: real(t)·real(t_0) + imag(t)·imag(t_0)
        correlations = rhok_real * rhok0_real + rhok_imag * rhok0_imag
        
        # Average over all k-vectors (k-averaging)
        return np.mean(correlations)
    
    def _output_correlation_data(self, ref_idx: int, current_time_ps: float):
        """Output correlation data for a specific reference frame."""
        if ref_idx >= len(self.references):
            return

        ref = self.references[ref_idx]
        lag_time_ps = current_time_ps - ref['time_ps']
        
        # For logarithmic spacing, we collect all data points and post-process
        # No early return here - we store all computed values
        
        # Get current density field  
        with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
            positions = np.array(snap.particles.position)
            
        current_rhok_real, current_rhok_imag = self._compute_density_field(positions)
        
        # Compute correlation
        fkt_value = self._compute_correlation(ref['rhok_real'], ref['rhok_imag'], 
                                            current_rhok_real, current_rhok_imag)
        
        if self.log_time_spacing:
            # Track which logarithmic points have been written for each reference
            if not hasattr(self, '_written_log_points'):
                self._written_log_points = {}
            
            if ref_idx not in self._written_log_points:
                self._written_log_points[ref_idx] = set()
            
            # Find the closest logarithmic time point
            time_diffs = [abs(lag_time_ps - t) for t in self.log_time_points_ps]
            closest_idx = min(range(len(time_diffs)), key=lambda i: time_diffs[i])
            closest_time = self.log_time_points_ps[closest_idx]
            min_diff = time_diffs[closest_idx]
            
            # Only write if:
            # 1. We're close enough (within 0.005 ps)
            # 2. This logarithmic point hasn't been written yet for this reference
            if min_diff < 0.005 and closest_idx not in self._written_log_points[ref_idx]:
                self._written_log_points[ref_idx].add(closest_idx)
                
                ref_file = f"{self.output_prefix}_fkt_ref_{ref_idx:03d}.txt"
                
                # Create file with header if it doesn't exist
                if not hasattr(self, f'_log_file_created_{ref_idx}'):
                    ref = self.references[ref_idx]
                    with open(ref_file, 'w') as f:
                        f.write("# F(k,t) correlation function\n")
                        f.write("# Reference time: {:.3f} ps\n".format(ref['time_ps']))
                        f.write("# Using logarithmic time spacing\n")
                        f.write("# Time range: {:.3f} - {:.1f} ps\n".format(self.min_log_time_ps, self.max_log_time_ps))
                        f.write("# lag_time_ps\tF(k,t)\n")
                    setattr(self, f'_log_file_created_{ref_idx}', True)
                
                # Append data point using the EXACT logarithmic time
                with open(ref_file, 'a') as f:
                    f.write(f"{closest_time:.6f}\t{fkt_value:.8f}\n")
                
                print(f"DEBUG: WROTE LOG POINT ref_idx {ref_idx}: target={closest_time:.6f}, actual_t={lag_time_ps:.6f}, F(k,t)={fkt_value:.3f}, diff={min_diff:.6f}")
        else:
            # Immediate linear output
            output_line = f"{lag_time_ps:.6f}\t{fkt_value:.8f}\n"
            
            # Write to reference-specific file
            ref_file = f"{self.output_prefix}_fkt_ref_{ref_idx:03d}.txt"
            with open(ref_file, 'a') as f:
                f.write(output_line)
        
        # Update cached value for logging
        self.current_computed_value = fkt_value
    
    def _write_logarithmic_output(self, ref_idx: int):
        """Write logarithmically-spaced correlation data for a reference frame."""
        if ref_idx not in self.correlation_data:
            print(f"DEBUG: No correlation data for ref_idx {ref_idx}")
            return
            
        times = np.array(self.correlation_data[ref_idx]['times'])
        values = np.array(self.correlation_data[ref_idx]['values'])
        
        print(f"DEBUG: ref_idx {ref_idx}: {len(times)} time points collected")
        if len(times) > 0:
            print(f"DEBUG: Time range: {times.min():.3f} - {times.max():.3f} ps")
            print(f"DEBUG: Target log times: {self.log_time_points_ps[:3]} ... {self.log_time_points_ps[-3:]}")
        
        if len(times) == 0:
            print(f"DEBUG: No data points collected for ref_idx {ref_idx}")
            return
        
        # Create header for reference file
        ref = self.references[ref_idx]
        ref_file = f"{self.output_prefix}_fkt_ref_{ref_idx:03d}.txt"
        
        with open(ref_file, 'w') as f:
            f.write("# F(k,t) correlation function\n")
            f.write("# Reference time: {:.3f} ps\n".format(ref['time_ps']))
            f.write("# Using logarithmic time spacing\n")
            f.write("# Time range: {:.3f} - {:.1f} ps\n".format(self.min_log_time_ps, self.max_log_time_ps))
            f.write("# lag_time_ps\tF(k,t)\n")
            
            # Interpolate data onto logarithmic time points
            points_written = 0
            for target_time in self.log_time_points_ps:
                if target_time <= times.max() and target_time >= times.min():
                    # Interpolate to get F(k,t) at this logarithmic time point
                    interpolated_value = np.interp(target_time, times, values)
                    f.write(f"{target_time:.6f}\t{interpolated_value:.8f}\n")
                    points_written += 1
            print(f"DEBUG: Wrote {points_written} logarithmic points for ref_idx {ref_idx}")
    
    def finalize_output(self):
        """Finalize all logarithmic output files."""
        if self.log_time_spacing:
            for ref_idx in range(len(self.references)):
                self._write_logarithmic_output(ref_idx)

    def act(self, timestep):
        """Main action called every timestep."""
        # Get current time
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should add a new reference frame
        if self._should_add_reference(current_time_ps):
            with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                positions = np.array(snap.particles.position)
                
            rhok_real, rhok_imag = self._compute_density_field(positions)
            
            # Create new reference
            ref = {
                'timestep': timestep,
                'time_ps': current_time_ps,
                'rhok_real': rhok_real,
                'rhok_imag': rhok_imag
            }
            
            self.references.append(ref)
            self.last_reference_time = current_time_ps
            
            # Create output file for this reference
            ref_idx = len(self.references) - 1
            ref_file = f"{self.output_prefix}_fkt_ref_{ref_idx:03d}.txt"
            with open(ref_file, 'w') as f:
                f.write("# F(k,t) correlation function\n")
                f.write("# Reference time: {:.3f} ps\n".format(current_time_ps))
                if self.log_time_spacing:
                    f.write("# Using logarithmic time spacing\n")
                    f.write("# Time range: {:.3f} - {:.1f} ps\n".format(self.min_log_time_ps, self.max_log_time_ps))
                f.write("# lag_time_ps\tF(k,t)\n")
            
            # Remove old references if we exceed max_references
            if len(self.references) > self.max_references:
                self.references.pop(0)
            
            print(f"Added F(k,t) reference {ref_idx} at t = {current_time_ps:.3f} ps")
        
        # Output correlation data for all reference frames
        if self.log_time_spacing:
            # For logarithmic spacing, collect data at every timestep for better interpolation
            for ref_idx in range(len(self.references)):
                self._output_correlation_data(ref_idx, current_time_ps)
        else:
            # For linear spacing, use the normal output period
            if self._should_output(current_time_ps):
                for ref_idx in range(len(self.references)):
                    self._output_correlation_data(ref_idx, current_time_ps)
                
                self.last_output_time = current_time_ps
        
        # Cache current timestep for logging
        self.current_timestep = timestep
    
    @hoomd.logging.log(category="scalar")
    def current_autocorr(self):
        """Return the most recent autocorrelation value for logging."""
        if self.current_computed_value is not None:
            return self.current_computed_value
        return 0.0


class AutocorrelationTracker(hoomd.custom.Action):
    """Base class for autocorrelation function tracking."""
    
    def __init__(self, simulation, time_tracker, output_period_ps, output_prefix):
        super().__init__()
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.output_period_ps = output_period_ps
        self.output_prefix = output_prefix
        self.last_output_time = None
    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output data."""
        if self.last_output_time is None:
            return True
        return (current_time_ps - self.last_output_time) >= self.output_period_ps


class EnergyTracker(hoomd.custom.Action):
    """
    SOPHISTICATED Energy Tracker - tracks comprehensive energy components during simulation.
    
    Tracks:
    - Individual potential energy components (harmonic, LJ, Coulomb, cavity)
    - Separate molecular and cavity kinetic energies
    - Reservoir energies from thermostats
    - Temperature
    - Universe total energy (conserved quantity)
    """
    
    def __init__(self, simulation, time_tracker, output_period_ps, output_prefix,
                 force_objects=None, thermostat_objects=None, verbose="quiet"):
        super().__init__()
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.output_period_ps = output_period_ps
        self.output_prefix = output_prefix
        self.last_output_time = 0.0
        
        # CRITICAL: Store force and thermostat objects for direct energy access
        self.force_objects = force_objects or {}
        self.thermostat_objects = thermostat_objects or {}
        
        # Verbosity control
        self.verbose = verbose.lower() if isinstance(verbose, str) else "normal"
        
        # Initialize current energy values for logging
        self._initialize_energy_values()

        # Initialize output file with comprehensive header
        self.output_file = f"{output_prefix}_energy_comprehensive.txt"
        self._initialize_output_file()
        
        # Physical constants (from old working version)
        self._kB = 3.166811563e-6  # Boltzmann constant in Hartree/K
        
        if self.verbose != "quiet":
            print(f"COMPREHENSIVE EnergyTracker initialized:", flush=True)
            print(f"  Output file: {self.output_file}", flush=True)
            print(f"  Force objects: {list(self.force_objects.keys())}", flush=True)
            print(f"  Thermostat objects: {list(self.thermostat_objects.keys())}", flush=True)
            print(f"  Output period: {self.output_period_ps:.3f} ps", flush=True)
    
    def _get_hoomd_simulation(self):
        """Get the HOOMD simulation object, handling both CavityMDSimulation and direct HOOMD objects."""
        if hasattr(self.simulation, 'sim'):
            # CavityMDSimulation object - get HOOMD simulation
            return self.simulation.sim
        else:
            # Direct HOOMD simulation object
            return self.simulation

    def _initialize_energy_values(self):
        """Initialize energy values for logging."""
        # Current energy components (matching old working version)
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
        """Initialize output file with comprehensive header (matching old working version)."""
        try:
            with open(self.output_file, "w") as f:
                f.write("# COMPREHENSIVE Energy tracking (based on old working version)\n")
                f.write(f"# Output period: {self.output_period_ps:.3f} ps\n")
                f.write("# All energies in Hartree (atomic units)\n")
                f.write("# Column definitions:\n")
                f.write("#   time(ps): simulation time in picoseconds\n")
                f.write("#   timestep: simulation timestep number\n")
                f.write("#   harmonic_energy: harmonic bond potential energy\n")
                f.write("#   lj_energy: Lennard-Jones potential energy\n")
                f.write("#   ewald_short_energy: short-range Coulomb energy\n")
                f.write("#   ewald_long_energy: long-range Coulomb energy\n")
                f.write("#   cavity_harmonic_energy: cavity harmonic potential energy\n")
                f.write("#   cavity_coupling_energy: cavity-molecule coupling energy\n")
                f.write("#   cavity_dipole_self_energy: dipole self-energy\n")
                f.write("#   cavity_total_potential_energy: total cavity potential energy\n")
                f.write("#   molecular_kinetic_energy: molecular kinetic energy\n")
                f.write("#   cavity_kinetic_energy: cavity kinetic energy\n")
                f.write("#   total_kinetic_energy: total kinetic energy\n")
                f.write("#   total_potential_energy: total potential energy\n")
                f.write("#   system_total_energy: total system energy (KE + PE)\n")
                f.write("#   molecular_reservoir_energy: molecular reservoir energy\n")
                f.write("#   cavity_reservoir_energy: cavity reservoir energy\n")
                f.write("#   total_reservoir_energy: total reservoir energy\n")
                f.write("#   universe_total_energy: universe total energy (system + reservoir) [CONSERVED]\n")
                f.write("#   temperature: kinetic temperature (K)\n")

                # Create header line
                header = "time(ps) timestep"
                header += " harmonic_energy lj_energy ewald_short_energy ewald_long_energy"
                header += " cavity_harmonic_energy cavity_coupling_energy cavity_dipole_self_energy cavity_total_potential_energy"
                header += " molecular_kinetic_energy cavity_kinetic_energy total_kinetic_energy"
                header += " total_potential_energy system_total_energy"
                header += " molecular_reservoir_energy cavity_reservoir_energy total_reservoir_energy"
                header += " universe_total_energy temperature\n"
                f.write(header)
                f.flush()
            if self.verbose != "quiet":
                print(f"EnergyTracker: Successfully created output file {self.output_file}", flush=True)
        except Exception as e:
            print(f"EnergyTracker ERROR: Failed to create output file {self.output_file}: {e}", flush=True)

    def act(self, timestep):
        """Track comprehensive energy components at each timestep (based on old working version)."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Skip timestep 0
        if timestep == 0:
            return

        # Only output periodically
        if not self._should_output(current_time_ps):
            return

        try:
            if self.verbose == "verbose":
                print(f"=== ENERGY TRACKER DEBUG - Timestep {timestep} ===", flush=True)
                print(f"Current time: {current_time_ps:.6f} ps", flush=True)

            # === 1. GET POTENTIAL ENERGY COMPONENTS (direct from force objects) ===
            # Get individual potential energy contributions directly from force objects
            try:
                self.current_harmonic_energy = (
                    self.force_objects["harmonic"].energy
                    if "harmonic" in self.force_objects
                    else 0.0
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
            except (AttributeError, KeyError) as e:
                self.current_ewald_short_energy = 0.0

            try:
                self.current_ewald_long_energy = (
                    self.force_objects["ewald_long"].energy
                    if "ewald_long" in self.force_objects
                    else 0.0
                    )
            except (AttributeError, KeyError) as e:
                self.current_ewald_long_energy = 0.0

            # Calculate molecular potential energy
            molecular_potential_energy = (
                self.current_harmonic_energy
                + self.current_lj_energy
                + self.current_ewald_short_energy
                + self.current_ewald_long_energy
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
                try:
                    # Use the logged property methods directly (from old working version)
                    self.current_cavity_harmonic_energy = cavityforce.harmonic_energy
                    self.current_cavity_coupling_energy = cavityforce.coupling_energy
                    self.current_cavity_dipole_self_energy = cavityforce.dipole_self_energy

                    # For total energy, try .energy property first, then sum components
                    if hasattr(cavityforce, "energy"):
                        self.current_cavity_total_potential_energy = cavityforce.energy
                    else:
                        self.current_cavity_total_potential_energy = (
                            self.current_cavity_harmonic_energy
                            + self.current_cavity_coupling_energy
                            + self.current_cavity_dipole_self_energy
                            )
                except Exception as e:
                    self.current_cavity_harmonic_energy = 0.0
                    self.current_cavity_coupling_energy = 0.0
                    self.current_cavity_dipole_self_energy = 0.0
                    self.current_cavity_total_potential_energy = 0.0
                    if self.verbose in ["normal", "verbose"]:
                        print(f"ERROR accessing cavity energy components: {e}", flush=True)

            # Calculate total potential energy
            self.current_total_potential_energy = (
                molecular_potential_energy + self.current_cavity_total_potential_energy
            )

            # === 2. GET KINETIC ENERGY COMPONENTS (using old working method) ===
            molecular_ke, molecular_temp = self._compute_molecular_kinetic_energy()
            cavity_ke = self._compute_cavity_kinetic_energy()

            self.current_molecular_kinetic_energy = molecular_ke
            self.current_cavity_kinetic_energy = cavity_ke
            self.current_total_kinetic_energy = molecular_ke + cavity_ke
            self.current_temperature = molecular_temp

            # === 3. GET RESERVOIR ENERGIES ===
            # Try to get reservoir energies from available thermostats
            self.current_molecular_reservoir_energy = 0.0
            self.current_cavity_reservoir_energy = 0.0
            
            # Check all available thermostat objects
            if self.verbose == "verbose":
                print(f"    Available thermostat objects: {list(self.thermostat_objects.keys())}", flush=True)
            
            for thermo_name, thermo_obj in self.thermostat_objects.items():
                if "molecular" in thermo_name or "bussi" in thermo_name:
                    reservoir_energy = self._get_reservoir_energy(thermo_name)
                    self.current_molecular_reservoir_energy += reservoir_energy
                    if self.verbose == "verbose":
                        print(f"    {thermo_name} reservoir energy: {reservoir_energy}", flush=True)
                elif "cavity" in thermo_name or "langevin" in thermo_name:
                    reservoir_energy = self._get_reservoir_energy(thermo_name)
                    self.current_cavity_reservoir_energy += reservoir_energy
                    if self.verbose == "verbose":
                        print(f"    {thermo_name} reservoir energy: {reservoir_energy}", flush=True)
            
            self.current_total_reservoir_energy = (
                self.current_molecular_reservoir_energy + self.current_cavity_reservoir_energy
            )
            
            if self.verbose == "verbose":
                print(f"    Reservoir energies: molecular={self.current_molecular_reservoir_energy}, cavity={self.current_cavity_reservoir_energy}", flush=True)

            # === 4. CALCULATE TOTAL ENERGIES ===
            self.current_system_total_energy = (
                self.current_total_potential_energy + self.current_total_kinetic_energy
            )
            self.current_universe_total_energy = (
                self.current_system_total_energy + self.current_total_reservoir_energy
            )

            # === 5. WRITE OUTPUT DATA ===
            self._write_energy_data(timestep, current_time_ps)

            if self.verbose == "verbose":
                print(f"=== END ENERGY TRACKER DEBUG - Timestep {timestep} ===\n", flush=True)

        except Exception as e:
            print(f"EnergyTracker CRITICAL ERROR at timestep {timestep}: {e}", flush=True)
            import traceback
            traceback.print_exc()

    def _compute_molecular_kinetic_energy(self):
        """
        Compute molecular kinetic energy and temperature internally (from old working version).
        
        Returns:
            tuple: (kinetic_energy, temperature) in atomic units and Kelvin
        """
        with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
            # Convert to numpy array for robust handling
            typeid_array = np.array(snap.particles.typeid)
            # Filter to molecular particles only (exclude cavity particle type)
            # Assume cavity particle has type ID 2 or is the last particle type
            max_typeid = np.max(typeid_array) if len(typeid_array) > 0 else 0
            molecular_mask = typeid_array != max_typeid  # Exclude highest type ID (cavity)

            if not np.any(molecular_mask):
                return 0.0, 0.0

            velocities_np = np.array(snap.particles.velocity)[molecular_mask]
            masses_np = np.array(snap.particles.mass)[molecular_mask]

            # Compute kinetic energy: KE = 0.5 * sum(m_i * v_i^2)
            kinetic_energy = 0.5 * np.sum(masses_np[:, np.newaxis] * velocities_np**2)

            # Compute temperature: T = (2/3) * KE / (N * k_B)
            N_dof = 3 * len(masses_np)  # 3 degrees of freedom per particle
            temperature = (2.0 * kinetic_energy / (N_dof * self._kB)) if N_dof > 0 else 0.0

            return float(kinetic_energy), float(temperature)

    def _compute_cavity_kinetic_energy(self):
        """
        Compute cavity kinetic energy internally (from old working version).
        
        Returns:
            float: cavity kinetic energy in atomic units
        """
        # Check if this is a cavity simulation
        # Handle both CavityMDSimulation and HOOMD simulation objects
        if hasattr(self.simulation, 'incavity'):
            # CavityMDSimulation object - has incavity attribute
            incavity = self.simulation.incavity
            hoomd_sim = self.simulation.sim  # Get HOOMD simulation from CavityMDSimulation
        else:
            # HOOMD simulation object - no incavity attribute, assume cavity simulation
            incavity = True  # Default for backward compatibility
            hoomd_sim = self.simulation
            
            if self.verbose == "verbose":
                print(f"    DEBUG: EnergyTracker simulation object = {type(self.simulation)}", flush=True)
                print(f"    DEBUG: incavity attribute = {incavity}", flush=True)
        
        if not incavity:
            # No-cavity simulation - no cavity particle exists
            if self.verbose == "verbose":
                print(f"    DEBUG: No-cavity simulation detected, returning 0.0 for cavity KE", flush=True)
            return 0.0
            
        with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
            typeid_array = np.array(snap.particles.typeid)
            
            # In cavity simulations, cavity particle has the highest type ID
            if len(typeid_array) == 0:
                return 0.0
                
            max_typeid = np.max(typeid_array)
            cavity_mask = typeid_array == max_typeid

            if not np.any(cavity_mask):
                return 0.0

            cavity_velocity_np = np.array(snap.particles.velocity)[cavity_mask]
            cavity_mass_np = np.array(snap.particles.mass)[cavity_mask]

            # Compute cavity kinetic energy: KE = 0.5 * sum(m * v^2)
            kinetic_energy = 0.5 * np.sum(cavity_mass_np[:, np.newaxis] * cavity_velocity_np**2)

            return float(kinetic_energy)

    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output at this time."""
        if self.last_output_time is None:
            self.last_output_time = 0.0
        return (current_time_ps - self.last_output_time) >= self.output_period_ps
    
    def _write_energy_data(self, timestep, current_time_ps):
        """Write energy data to output file (from old working version)."""
        try:
            # Build output line with all components
            output_values = [
                current_time_ps,
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
                # Temperature
                self.current_temperature,
            ]

            # Write to file
            with open(self.output_file, "a") as f:
                formatted_values = [
                    f"{val:.6f}" if isinstance(val, float) else str(val)
                    for val in output_values
                ]
                f.write(" ".join(formatted_values) + "\n")
                f.flush()

            # Update last output time
            self.last_output_time = current_time_ps

            if self.verbose == "verbose":
                print(f"Successfully wrote energy data to {self.output_file}", flush=True)

        except Exception as e:
            print(f"EnergyTracker ERROR writing data at timestep {timestep}: {e}", flush=True)
            import traceback
            traceback.print_exc()

    def get_instantaneous_energy(self):
        """
        Get current energy values for compatibility with temperature tracker.
        
        Returns dictionary with energy components that other trackers expect.
        """
        # Update energies from force objects first
        self._update_current_energies()
        
        return {
            'harmonic': self.current_harmonic_energy,
            'lj': self.current_lj_energy,
            'coulombic': self.current_ewald_short_energy + self.current_ewald_long_energy,
            'ewald_short': self.current_ewald_short_energy,
            'ewald_long': self.current_ewald_long_energy,
            'cavity_harmonic': self.current_cavity_harmonic_energy,
            'cavity_coupling': self.current_cavity_coupling_energy,
            'cavity_dipole': self.current_cavity_dipole_self_energy,
            'cavity_total': self.current_cavity_total_potential_energy,
            'molecular_kinetic': self.current_molecular_kinetic_energy,
            'cavity_kinetic': self.current_cavity_kinetic_energy,
            'total_kinetic': self.current_total_kinetic_energy,
            'total_potential': self.current_total_potential_energy,
            'system_total': self.current_system_total_energy,
            'universe_total': self.current_universe_total_energy,
            'temperature': self.current_temperature
        }
    
    def _update_current_energies(self):
        """Update current energy values from force objects."""
        try:
            # Reset all energies
            self.current_harmonic_energy = 0.0
            self.current_lj_energy = 0.0
            self.current_ewald_short_energy = 0.0
            self.current_ewald_long_energy = 0.0
            self.current_cavity_harmonic_energy = 0.0
            self.current_cavity_coupling_energy = 0.0
            self.current_cavity_dipole_self_energy = 0.0
            
            # Update from force objects with debug info
            if self.verbose == "verbose":
                print(f"  Force objects available: {list(self.force_objects.keys())}", flush=True)
            
            for force_name, force_obj in self.force_objects.items():
                try:
                    if force_name == "harmonic":
                        energy_val = getattr(force_obj, 'energy', 0.0)
                        self.current_harmonic_energy = float(energy_val) if hasattr(energy_val, '__float__') else energy_val
                        if self.verbose == "verbose":
                            print(f"    Harmonic energy: {self.current_harmonic_energy}", flush=True)
                    elif force_name == "lj":
                        energy_val = getattr(force_obj, 'energy', 0.0)
                        self.current_lj_energy = float(energy_val) if hasattr(energy_val, '__float__') else energy_val
                        if self.verbose == "verbose":
                            print(f"    LJ energy: {self.current_lj_energy}", flush=True)
                    elif force_name == "ewald_short":
                        energy_val = getattr(force_obj, 'energy', 0.0)
                        self.current_ewald_short_energy = float(energy_val) if hasattr(energy_val, '__float__') else energy_val
                        if self.verbose == "verbose":
                            print(f"    Ewald short energy: {self.current_ewald_short_energy}", flush=True)
                    elif force_name == "ewald_long":
                        energy_val = getattr(force_obj, 'energy', 0.0)
                        self.current_ewald_long_energy = float(energy_val) if hasattr(energy_val, '__float__') else energy_val
                        if self.verbose == "verbose":
                            print(f"    Ewald long energy: {self.current_ewald_long_energy}", flush=True)
                    elif force_name == "cavity":
                        # Cavity force has multiple energy components
                        if hasattr(force_obj, 'harmonic_energy'):
                            self.current_cavity_harmonic_energy = force_obj.harmonic_energy
                        if hasattr(force_obj, 'coupling_energy'):
                            self.current_cavity_coupling_energy = force_obj.coupling_energy
                        if hasattr(force_obj, 'dipole_self_energy'):
                            self.current_cavity_dipole_self_energy = force_obj.dipole_self_energy
                        if self.verbose == "verbose":
                            print(f"    Cavity energies: harmonic={self.current_cavity_harmonic_energy}, coupling={self.current_cavity_coupling_energy}, dipole={self.current_cavity_dipole_self_energy}", flush=True)
                except Exception as e:
                    if self.verbose == "verbose":
                        print(f"    Error accessing {force_name} energy: {e}", flush=True)
            
            # Update kinetic energies
            molecular_ke, _ = self._compute_molecular_kinetic_energy()
            self.current_molecular_kinetic_energy = molecular_ke
            
            # Get cavity kinetic energy using the fixed method
            self.current_cavity_kinetic_energy = self._compute_cavity_kinetic_energy()
            if self.verbose == "verbose":
                print(f"    Cavity kinetic energy: {self.current_cavity_kinetic_energy}", flush=True)
            
            self.current_total_kinetic_energy = self.current_molecular_kinetic_energy + self.current_cavity_kinetic_energy
            
            # Update derived values
            self.current_cavity_total_potential_energy = (self.current_cavity_harmonic_energy + 
                                                        self.current_cavity_coupling_energy + 
                                                        self.current_cavity_dipole_self_energy)
            
            self.current_total_potential_energy = (self.current_harmonic_energy + self.current_lj_energy + 
                                                 self.current_ewald_short_energy + self.current_ewald_long_energy +
                                                 self.current_cavity_total_potential_energy)
            
            self.current_system_total_energy = self.current_total_kinetic_energy + self.current_total_potential_energy
            
            # Temperature from molecular kinetic energy
            if self.current_molecular_kinetic_energy > 0:
                # Get number of molecular particles (excluding cavity if it exists)
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    total_particles = len(snap.particles.mass)
                    # Subtract 1 only if cavity particle exists
                    if "cavity" in self.force_objects:
                        n_particles = total_particles - 1  # Exclude cavity particle
                    else:
                        n_particles = total_particles  # All particles are molecular
                    n_dof = 3 * n_particles if n_particles > 0 else 3
                    self.current_temperature = (2.0 * self.current_molecular_kinetic_energy / (n_dof * self._kB))
                    if self.verbose == "verbose":
                        print(f"    Temperature: {self.current_temperature:.2f} K (from {n_particles} molecular particles)", flush=True)
            else:
                self.current_temperature = 0.0
                        
        except Exception as e:
            print(f"Warning: Error updating current energies: {e}")
            # Keep existing values if update fails
    
    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output data."""
        if self.last_output_time is None:
            return True
        return (current_time_ps - self.last_output_time) >= self.output_period_ps
    
    def _get_reservoir_energy(self, thermostat_name):
        """Get reservoir energy from a thermostat if available."""
        if thermostat_name is None or thermostat_name not in self.thermostat_objects:
            if self.verbose == "verbose":
                print(f"    DEBUG: Thermostat {thermostat_name} not found in thermostat_objects", flush=True)
            return 0.0
        
        thermostat = self.thermostat_objects[thermostat_name]
        if thermostat is None:
            if self.verbose == "verbose":
                print(f"    DEBUG: Thermostat {thermostat_name} is None", flush=True)
            return 0.0
        
        if self.verbose == "verbose":
            print(f"    DEBUG: Checking reservoir energy for {thermostat_name}", flush=True)
            print(f"      Thermostat type: {type(thermostat)}", flush=True)
            print(f"      Has total_reservoir_energy: {hasattr(thermostat, 'total_reservoir_energy')}", flush=True)
            print(f"      Has reservoir_energy: {hasattr(thermostat, 'reservoir_energy')}", flush=True)
            print(f"      Has energy: {hasattr(thermostat, 'energy')}", flush=True)
            print(f"      Has thermostat: {hasattr(thermostat, 'thermostat')}", flush=True)
            if hasattr(thermostat, 'thermostat'):
                print(f"      Nested thermostat type: {type(thermostat.thermostat)}", flush=True)
                print(f"      Nested has reservoir_energy: {hasattr(thermostat.thermostat, 'reservoir_energy')}", flush=True)
        
        try:
            # Try different ways to access reservoir energy
            if hasattr(thermostat, 'total_reservoir_energy'):
                # BussiReservoir thermostat
                energy_val = thermostat.total_reservoir_energy
                if self.verbose == "verbose":
                    print(f"      Found total_reservoir_energy = {energy_val}", flush=True)
                return float(energy_val) if hasattr(energy_val, '__float__') else energy_val
            elif hasattr(thermostat, 'reservoir_energy'):
                # Standard reservoir_energy attribute
                energy_val = thermostat.reservoir_energy
                if self.verbose == "verbose":
                    print(f"      Found reservoir_energy = {energy_val}", flush=True)
                return float(energy_val) if hasattr(energy_val, '__float__') else energy_val
            elif hasattr(thermostat, 'energy'):
                # Generic energy attribute
                energy_val = thermostat.energy
                if self.verbose == "verbose":
                    print(f"      Found energy = {energy_val}", flush=True)
                return float(energy_val) if hasattr(energy_val, '__float__') else energy_val
            elif hasattr(thermostat, 'thermostat') and hasattr(thermostat.thermostat, 'total_reservoir_energy'):
                # For nested thermostats (e.g., ConstantVolume with BussiReservoir)
                energy_val = thermostat.thermostat.total_reservoir_energy
                if self.verbose == "verbose":
                    print(f"      Found nested total_reservoir_energy = {energy_val}", flush=True)
                return float(energy_val) if hasattr(energy_val, '__float__') else energy_val
            elif hasattr(thermostat, 'thermostat') and hasattr(thermostat.thermostat, 'reservoir_energy'):
                # For nested thermostats (e.g., ConstantVolume with standard reservoir)
                energy_val = thermostat.thermostat.reservoir_energy
                if self.verbose == "verbose":
                    print(f"      Found nested reservoir_energy = {energy_val}", flush=True)
                return float(energy_val) if hasattr(energy_val, '__float__') else energy_val
            else:
                if self.verbose == "verbose":
                    print(f"      No reservoir energy attribute found", flush=True)
        except Exception as e:
            if self.verbose == "verbose":
                print(f"    Error accessing reservoir energy from {thermostat_name}: {e}", flush=True)
        
        return 0.0
    
    def _write_comprehensive_header(self, energy_data):
        """Write sophisticated header matching the expected format."""
        with open(self.output_file, 'w') as f:
            f.write("# SOPHISTICATED Energy tracking (comprehensive output)\n")
            f.write(f"# Output period: {self.output_period_ps:.3f} ps\n")
            f.write("# All energies in Hartree (atomic units)\n")
            f.write("# Column definitions:\n")
            f.write("#   time(ps): simulation time in picoseconds\n")
            f.write("#   timestep: simulation timestep number\n")
            f.write("#   harmonic_energy: harmonic bond potential energy\n")
            f.write("#   lj_energy: Lennard-Jones potential energy\n")
            f.write("#   ewald_short_energy: short-range Coulomb energy\n")
            f.write("#   ewald_long_energy: long-range Coulomb energy\n")
            f.write("#   cavity_harmonic_energy: cavity harmonic potential energy\n")
            f.write("#   cavity_coupling_energy: cavity-molecule coupling energy\n")
            f.write("#   cavity_dipole_self_energy: dipole self-energy\n")
            f.write("#   cavity_total_potential_energy: total cavity potential energy\n")
            f.write("#   molecular_kinetic_energy: molecular kinetic energy\n")
            f.write("#   cavity_kinetic_energy: cavity kinetic energy\n")
            f.write("#   total_kinetic_energy: total kinetic energy\n")
            f.write("#   total_potential_energy: total potential energy\n")
            f.write("#   system_total_energy: total system energy (KE + PE)\n")
            f.write("#   molecular_reservoir_energy: molecular reservoir energy\n")
            f.write("#   cavity_reservoir_energy: cavity reservoir energy\n")
            f.write("#   total_reservoir_energy: total reservoir energy\n")
            f.write("#   universe_total_energy: universe total energy (system + reservoir) [CONSERVED]\n")
            f.write("#   temperature: kinetic temperature (K)\n")
            
            # Write column headers
            f.write("time(ps) timestep harmonic_energy lj_energy ewald_short_energy ewald_long_energy " +
                   "cavity_harmonic_energy cavity_coupling_energy cavity_dipole_self_energy " +
                   "cavity_total_potential_energy molecular_kinetic_energy cavity_kinetic_energy " +
                   "total_kinetic_energy total_potential_energy system_total_energy " +
                   "molecular_reservoir_energy cavity_reservoir_energy total_reservoir_energy " +
                   "universe_total_energy temperature\n")
    
    def _format_energy_line(self, current_time_ps, timestep, data):
        """Format energy data into output line."""
        return (f"{current_time_ps:.6f} {timestep} " +
                f"{data['harmonic_energy']:.6f} " +
                f"{data['lj_energy']:.6f} " +
                f"{data['ewald_short_energy']:.6f} " +
                f"{data['ewald_long_energy']:.6f} " +
                f"{data['cavity_harmonic_energy']:.6f} " +
                f"{data['cavity_coupling_energy']:.6f} " +
                f"{data['cavity_dipole_self_energy']:.6f} " +
                f"{data['cavity_total_potential_energy']:.6f} " +
                f"{data['molecular_kinetic_energy']:.6f} " +
                f"{data['cavity_kinetic_energy']:.6f} " +
                f"{data['total_kinetic_energy']:.6f} " +
                f"{data['total_potential_energy']:.6f} " +
                f"{data['system_total_energy']:.6f} " +
                f"{data['molecular_reservoir_energy']:.6f} " +
                f"{data['cavity_reservoir_energy']:.6f} " +
                f"{data['total_reservoir_energy']:.6f} " +
                f"{data['universe_total_energy']:.6f} " +
                f"{data['temperature']:.6f}\n")
    


class PerformanceTracker(hoomd.custom.Action):
    """Tracks simulation performance metrics."""
    
    def __init__(self, simulation, runtime_ps, time_tracker):
        super().__init__()
        self.simulation = simulation
        self.runtime_ps = runtime_ps
        self.time_tracker = time_tracker
        self.last_output_time = None
        self.start_time = time.time()
        self.last_timestep = 0
        self._current_steps_per_second = 0.0
        self._current_ns_per_day = "0.0"
        
    def act(self, timestep):
        """Track performance metrics."""
        current_time_ps = self.time_tracker.elapsed_time
        wall_time = time.time() - self.start_time
        
        if wall_time > 0:
            # Calculate steps per second
            self._current_steps_per_second = timestep / wall_time
            
            # Calculate ns per day estimate
            if current_time_ps > 0:
                ps_per_second = current_time_ps / wall_time
                ns_per_second = ps_per_second / 1000.0
                ns_per_day = ns_per_second * 86400.0  # seconds per day
                self._current_ns_per_day = f"{ns_per_day:.2f}"
    
    @hoomd.logging.log(category="scalar")
    def steps_per_second(self):
        """Return current steps per second."""
        return self._current_steps_per_second
        
    @hoomd.logging.log(category="string")
    def ns_per_day(self):
        """Return estimated ns per day."""
        return self._current_ns_per_day
        
    @hoomd.logging.log(category="scalar")  
    def wall_time(self):
        """Return total wall time in seconds."""
        return time.time() - self.start_time
        
    @hoomd.logging.log(category="string")
    def eta_remaining(self):
        """Return estimated time remaining."""
        if hasattr(self.time_tracker, 'elapsed_time'):
            # Check if elapsed_time is a method or property
            if callable(self.time_tracker.elapsed_time):
                current_time_ps = self.time_tracker.elapsed_time()
            else:
                current_time_ps = self.time_tracker.elapsed_time
        else:
            return "00:00:00"
            
        if current_time_ps > 0:
            wall_time = time.time() - self.start_time
            progress = current_time_ps / self.runtime_ps
            if progress > 0:
                eta_seconds = wall_time * (1.0 - progress) / progress
                hours = int(eta_seconds // 3600)
                minutes = int((eta_seconds % 3600) // 60)
                seconds = int(eta_seconds % 60)
                return f"{hours:02d}:{minutes:02d}:{seconds:02d}"
        return "00:00:00"


class EmpiricalTemperatureData:
    """
    Handles empirical potential energy vs temperature data with extended function fitting.
    
    This class loads empirical energy-temperature relationships from equilibrium
    simulations and provides interpolation to determine systemic temperatures
    from instantaneous potential energies using extended functions:
    - Extended harmonic: E = aT/(1+bT) for harmonic energy
    - Extended T^(3/5): E = E₀ + AT^(3/5)/(1+CT^(3/5)) for LJ+Coulombic energy
    
    Parameters
    ----------
    data_file_path : str
        Path to file containing temperature and potential energy data
    energy_component : str, optional
        Which energy component to use ('lj_coulombic', 'harmonic', etc.). Default: 'lj_coulombic'
    use_direct_harmonic : bool, optional
        If True and energy_component='harmonic', use direct calculation T = 4*E/(N*kB). Default: False
    create_plots : bool, optional
        If True, create diagnostic plots showing fitted functions vs data. Default: False
        
    Attributes
    ----------
    temperatures : np.ndarray
        Temperature values from empirical data (K)
    energies : np.ndarray  
        Energy values from empirical data (hartree)
    has_extended_harmonic_fit : bool
        Whether extended harmonic fitting has been performed
    has_extended_t35_fit : bool
        Whether extended T^(3/5) fitting has been performed
    extended_harmonic_fit : dict
        Parameters from extended harmonic fitting: {'a': float, 'b': float, 'r2': float}
    extended_t35_fit : dict
        Parameters from extended T^(3/5) fitting: {'e0': float, 'a': float, 'c': float, 'r2': float}
    """
    
    def __init__(self, data_file_path: str, energy_component: str = 'lj_coulombic', use_direct_harmonic: bool = False, create_plots: bool = False):
        self.data_file_path = Path(data_file_path)
        self.energy_component = energy_component
        self.use_direct_harmonic = use_direct_harmonic
        self.create_plots = create_plots
        
        # Initialize extended fitting attributes (no more simple power law)
        self.has_extended_harmonic_fit = False
        self.has_extended_t35_fit = False
        self.extended_harmonic_fit = {}
        self.extended_t35_fit = {}
        
        self.load_empirical_data()
        
        # Use extended functions only (no more simple power law)
        if not self.use_direct_harmonic or self.energy_component != 'harmonic':
            if self.energy_component == 'harmonic':
                self.fit_extended_harmonic_function()
            else:
                self.fit_extended_t35_function()
            
            # Create diagnostic plots if requested
            if self.create_plots:
                self.plot_fits()
    
    def load_empirical_data(self):
        """Load empirical energy-temperature data from file."""
        if not self.data_file_path.exists():
            raise FileNotFoundError(f"Empirical data file not found: {self.data_file_path}")
        
        try:
            import pandas as pd
            data = pd.read_csv(self.data_file_path, sep=r'\s+', comment='#')
            
            # Extract temperature data
            if 'temperature' in data.columns:
                self.temperatures = data['temperature'].values
            else:
                raise ValueError("Temperature column not found in empirical data")
            
            # Extract energy data based on component
            if self.energy_component == 'lj_coulombic':
                if 'lj_hartree' in data.columns and 'coulombic_hartree' in data.columns:
                    self.energies = data['lj_hartree'].values + data['coulombic_hartree'].values
                else:
                    raise ValueError("LJ and Coulombic energy columns not found")
            elif self.energy_component == 'total_PE':
                if 'total_potential_energy_hartree' in data.columns:
                    self.energies = data['total_potential_energy_hartree'].values
                else:
                    raise ValueError("Total potential energy column not found")
            elif self.energy_component == 'harmonic':
                if 'harmonic_hartree' in data.columns:
                    self.energies = data['harmonic_hartree'].values
                else:
                    raise ValueError("Harmonic energy column not found")
            else:
                raise ValueError(f"Unknown energy component: {self.energy_component}")
                
            print(f"Loaded {len(self.temperatures)} empirical data points")
            print(f"Temperature range: {self.temperatures.min():.1f} - {self.temperatures.max():.1f} K")
            print(f"Energy range: {self.energies.min():.6f} - {self.energies.max():.6f} Ha")
            
        except Exception as e:
            raise RuntimeError(f"Failed to load empirical data: {e}")
    
    # Removed fit_rosenfeld_function - using only extended functions for better accuracy
    
    def fit_extended_harmonic_function(self):
        """Fit extended harmonic function: E = aT/(1+bT)"""
        try:
            from scipy.optimize import curve_fit
            
            def extended_harmonic_model(T, a, b):
                """Extended harmonic model: E(T) = aT/(1+bT)"""
                return a * T / (1 + b * T)
            
            # Initial parameter guesses
            a_guess = self.energies.max() / self.temperatures.max()  # Linear coefficient
            b_guess = 1e-3  # Small positive value
            
            # Perform fitting
            self.extended_harmonic_params, _ = curve_fit(
                extended_harmonic_model,
                self.temperatures,
                self.energies,
                p0=[a_guess, b_guess],
                bounds=([0, 0], [np.inf, np.inf])  # Both parameters must be positive
            )
            
            # Calculate R-squared
            predicted = extended_harmonic_model(self.temperatures, *self.extended_harmonic_params)
            residuals = self.energies - predicted
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((self.energies - np.mean(self.energies))**2)
            r2 = 1 - (ss_res / ss_tot)
            
            # Store fitted parameters
            self.extended_harmonic_fit = {
                'a': self.extended_harmonic_params[0],
                'b': self.extended_harmonic_params[1],
                'r2': r2
            }
            
            self.has_extended_harmonic_fit = True
            
            print(f"Extended harmonic fitting completed:")
            print(f"  a = {self.extended_harmonic_fit['a']:.6f} Ha·K^(-1)")
            print(f"  b = {self.extended_harmonic_fit['b']:.6f} K^(-1)")
            print(f"  R² = {self.extended_harmonic_fit['r2']:.12f}")
            
        except Exception as e:
            print(f"Warning: Extended harmonic fitting failed: {e}")
            self.has_extended_harmonic_fit = False
    
    def fit_extended_t35_function(self):
        """Fit extended T^(3/5) function: E = E0 + AT^(3/5)/(1+CT^(3/5))"""
        try:
            from scipy.optimize import curve_fit
            
            def extended_t35_model(T, e0, a, c):
                """Extended T^(3/5) model: E(T) = E₀ + AT^(3/5)/(1+CT^(3/5))"""
                t_power = T**(3/5)
                return e0 + a * t_power / (1 + c * t_power)
            
            # Initial parameter guesses
            e0_guess = self.energies.min()
            a_guess = (self.energies.max() - self.energies.min()) / (self.temperatures.max()**(3/5))
            c_guess = 1e-3
            
            # Perform fitting
            self.extended_t35_params, _ = curve_fit(
                extended_t35_model,
                self.temperatures,
                self.energies,
                p0=[e0_guess, a_guess, c_guess],
                bounds=([self.energies.min() - 10, -np.inf, 0], 
                       [self.energies.max() + 10, np.inf, np.inf])
            )
            
            # Calculate R-squared
            predicted = extended_t35_model(self.temperatures, *self.extended_t35_params)
            residuals = self.energies - predicted
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((self.energies - np.mean(self.energies))**2)
            r2 = 1 - (ss_res / ss_tot)
            
            # Store fitted parameters
            self.extended_t35_fit = {
                'e0': self.extended_t35_params[0],
                'a': self.extended_t35_params[1],
                'c': self.extended_t35_params[2],
                'r2': r2
            }
            
            self.has_extended_t35_fit = True
            
            print(f"Extended T^(3/5) fitting completed:")
            print(f"  E₀ = {self.extended_t35_fit['e0']:.6f} Ha")
            print(f"  A = {self.extended_t35_fit['a']:.6f} Ha·K^(-3/5)")
            print(f"  C = {self.extended_t35_fit['c']:.6f} K^(-3/5)")
            print(f"  R² = {self.extended_t35_fit['r2']:.12f}")
            
        except Exception as e:
            print(f"Warning: Extended T^(3/5) fitting failed: {e}")
            self.has_extended_t35_fit = False
    
    def plot_fits(self, output_file: str = None, show_plot: bool = False):
        """
        Plot empirical data and fitted functions for visual validation.
        
        Parameters
        ----------
        output_file : str, optional
            Output file path for saving plot. If None, uses component-based naming.
        show_plot : bool, optional
            Whether to display plot interactively. Default: False
        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.style
            
            # Set up scientific plotting style
            plt.style.use('default')
            plt.rcParams.update({
                'font.size': 12,
                'font.family': 'serif',
                'mathtext.fontset': 'cm',
                'axes.grid': True,
                'grid.alpha': 0.3,
                'figure.dpi': 150,
                'savefig.dpi': 300,
                'savefig.bbox': 'tight'
            })
            
            # Create figure with 2x2 subplots
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
            
            # Generate fine temperature array for smooth curves
            T_fine = np.linspace(self.temperatures.min(), self.temperatures.max(), 1000)
            
            # Plot 1: Linear scale
            ax1.scatter(self.temperatures, self.energies, color='blue', alpha=0.7, s=50, 
                       label='Empirical Data', zorder=3)
            
            # Plot 2: Log-log scale (same data)
            # Only plot positive values for log scale
            pos_mask = (self.temperatures > 0) & (self.energies > 0)
            if np.any(pos_mask):
                ax2.loglog(self.temperatures[pos_mask], self.energies[pos_mask], 'bo', 
                          alpha=0.7, markersize=6, label='Empirical Data')
            
            # Plot fitted functions on both linear and log-log scales
            if self.has_extended_harmonic_fit and self.energy_component == 'harmonic':
                a, b = self.extended_harmonic_fit['a'], self.extended_harmonic_fit['b']
                r2 = self.extended_harmonic_fit['r2']
                
                # Extended harmonic function
                E_fit = a * T_fine / (1 + b * T_fine)
                ax1.plot(T_fine, E_fit, 'r-', linewidth=2, alpha=0.8,
                        label=f'Extended Harmonic: $E = \\frac{{aT}}{{1+bT}}$\n'
                              f'a = {a:.6f} Ha/K, b = {b:.6f} K$^{{-1}}$\nR² = {r2:.6f}')
                
                # Linear approximation for comparison
                E_linear = a * T_fine
                ax1.plot(T_fine, E_linear, 'g--', linewidth=1, alpha=0.6,
                        label=f'Linear: $E = aT$ (a = {a:.6f})')
                
                # Log-log plots (only positive values)
                pos_fit_mask = E_fit > 0
                if np.any(pos_fit_mask):
                    ax2.loglog(T_fine[pos_fit_mask], E_fit[pos_fit_mask], 'r-', 
                              linewidth=2, alpha=0.8, label='Extended Harmonic')
                    ax2.loglog(T_fine, E_linear, 'g--', linewidth=1, alpha=0.6, 
                              label='Linear: $E \\propto T$')
                        
            elif self.has_extended_t35_fit and self.energy_component != 'harmonic':
                e0, a, c = self.extended_t35_fit['e0'], self.extended_t35_fit['a'], self.extended_t35_fit['c']
                r2 = self.extended_t35_fit['r2']
                
                # Extended T^(3/5) function
                t_power = T_fine**(3/5)
                E_fit = e0 + a * t_power / (1 + c * t_power)
                ax1.plot(T_fine, E_fit, 'r-', linewidth=2, alpha=0.8,
                        label=f'Extended T$^{{3/5}}$: $E = E_0 + \\frac{{AT^{{3/5}}}}{{1+CT^{{3/5}}}}$\n'
                              f'E₀ = {e0:.6f} Ha, A = {a:.6f} Ha·K$^{{-3/5}}$\n'
                              f'C = {c:.6f} K$^{{-3/5}}$, R² = {r2:.6f}')
                
                # Simple T^(3/5) for comparison if c is small
                if c < 1e-3:
                    E_simple = e0 + a * t_power
                    ax1.plot(T_fine, E_simple, 'g--', linewidth=1, alpha=0.6,
                            label=f'Simple T$^{{3/5}}$: $E = E_0 + AT^{{3/5}}$')
                
                # Log-log plots - need to handle negative energies
                # Shift data to make it positive for log plotting
                E_shift = E_fit - e0  # Remove baseline
                E_simple_shift = a * t_power if c < 1e-3 else None
                
                pos_shift_mask = E_shift > 0
                if np.any(pos_shift_mask):
                    ax2.loglog(T_fine[pos_shift_mask], E_shift[pos_shift_mask], 'r-', 
                              linewidth=2, alpha=0.8, 
                              label=f'Extended T$^{{3/5}}$ (shifted: $E - E_0$)')
                    
                    if E_simple_shift is not None:
                        ax2.loglog(T_fine, E_simple_shift, 'g--', linewidth=1, alpha=0.6,
                                  label='Pure T$^{{3/5}}$: $E \\propto T^{{3/5}}$')
            
            # Format linear plot
            ax1.set_xlabel('Temperature (K)')
            ax1.set_ylabel('Energy (Ha)')
            ax1.set_title(f'Linear Scale: {self.energy_component.replace("_", "+").title()}')
            ax1.legend(fontsize=9)
            ax1.grid(True, alpha=0.3)
            
            # Format log-log plot
            ax2.set_xlabel('Temperature (K)')
            ax2.set_ylabel('Energy (Ha)' if self.energy_component == 'harmonic' else 'Energy - E₀ (Ha)')
            ax2.set_title(f'Log-Log Scale: Power Law Analysis')
            ax2.legend(fontsize=9)
            ax2.grid(True, alpha=0.3)
            
            # Plot 3: Residuals
            if self.has_extended_harmonic_fit and self.energy_component == 'harmonic':
                a, b = self.extended_harmonic_fit['a'], self.extended_harmonic_fit['b']
                E_pred = a * self.temperatures / (1 + b * self.temperatures)
            elif self.has_extended_t35_fit and self.energy_component != 'harmonic':
                e0, a, c = self.extended_t35_fit['e0'], self.extended_t35_fit['a'], self.extended_t35_fit['c']
                t_power = self.temperatures**(3/5)
                E_pred = e0 + a * t_power / (1 + c * t_power)
            else:
                E_pred = np.interp(self.temperatures, self.temperatures, self.energies)  # No fit available
            
            residuals = self.energies - E_pred
            ax3.scatter(self.temperatures, residuals, color='red', alpha=0.7, s=50)
            ax3.axhline(y=0, color='black', linestyle='--', alpha=0.5)
            ax3.set_xlabel('Temperature (K)')
            ax3.set_ylabel('Residuals (Ha)')
            ax3.set_title('Fit Residuals')
            ax3.grid(True, alpha=0.3)
            
            # Add statistics text
            rmse = np.sqrt(np.mean(residuals**2))
            max_abs_residual = np.max(np.abs(residuals))
            ax3.text(0.05, 0.95, f'RMSE: {rmse:.2e} Ha\nMax |residual|: {max_abs_residual:.2e} Ha',
                    transform=ax3.transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
            
            # Plot 4: Power law analysis
            if self.energy_component == 'harmonic':
                # For harmonic: analyze E vs T scaling
                if np.any(pos_mask):
                    # Fit power law: log(E) = log(A) + α*log(T)
                    log_T = np.log10(self.temperatures[pos_mask])
                    log_E = np.log10(self.energies[pos_mask])
                    
                    # Linear fit in log space
                    coeffs = np.polyfit(log_T, log_E, 1)
                    alpha_fit = coeffs[0]  # Power law exponent
                    log_A_fit = coeffs[1]  # log(prefactor)
                    A_fit = 10**log_A_fit
                    
                    # Calculate R²
                    log_E_pred = coeffs[0] * log_T + coeffs[1]
                    ss_res = np.sum((log_E - log_E_pred)**2)
                    ss_tot = np.sum((log_E - np.mean(log_E))**2)
                    r2_power = 1 - (ss_res / ss_tot)
                    
                    # Plot power law fit
                    T_power_fit = A_fit * T_fine**alpha_fit
                    ax4.loglog(T_fine, T_power_fit, 'purple', linewidth=2, alpha=0.8,
                              label=f'Power Law: $E = AT^{{\\alpha}}$\n'
                                    f'A = {A_fit:.2e} Ha, α = {alpha_fit:.3f}\n'
                                    f'R² = {r2_power:.6f}')
                    
                    # Add theoretical lines
                    E_theory_linear = (self.energies[pos_mask].mean() / self.temperatures[pos_mask].mean()) * T_fine
                    ax4.loglog(T_fine, E_theory_linear, 'k--', alpha=0.5, 
                              label='Linear: $E \\propto T^1$ (theory)')
                    
                    ax4.loglog(self.temperatures[pos_mask], self.energies[pos_mask], 'bo', 
                              alpha=0.7, markersize=6, label='Data')
                    
                    ax4.set_title(f'Power Law Fit: α = {alpha_fit:.3f} (expect 1.0)')
                    
            else:
                # For LJ+Coulombic: analyze (E - E₀) vs T^(3/5) scaling
                if self.has_extended_t35_fit:
                    e0, a, c = self.extended_t35_fit['e0'], self.extended_t35_fit['a'], self.extended_t35_fit['c']
                    E_shifted = self.energies - e0
                    pos_shifted_mask = E_shifted > 0
                    
                    if np.any(pos_shifted_mask):
                        T_shifted = self.temperatures[pos_shifted_mask]
                        E_shifted_pos = E_shifted[pos_shifted_mask]
                        
                        # Fit power law to shifted data
                        log_T = np.log10(T_shifted)
                        log_E_shift = np.log10(E_shifted_pos)
                        
                        coeffs = np.polyfit(log_T, log_E_shift, 1)
                        alpha_fit = coeffs[0]
                        log_A_fit = coeffs[1]
                        A_fit = 10**log_A_fit
                        
                        # Calculate R²
                        log_E_pred = coeffs[0] * log_T + coeffs[1]
                        ss_res = np.sum((log_E_shift - log_E_pred)**2)
                        ss_tot = np.sum((log_E_shift - np.mean(log_E_shift))**2)
                        r2_power = 1 - (ss_res / ss_tot)
                        
                        # Plot power law fit
                        T_power_fit = A_fit * T_fine**alpha_fit
                        ax4.loglog(T_fine, T_power_fit, 'purple', linewidth=2, alpha=0.8,
                                  label=f'Power Law: $(E-E_0) = AT^{{\\alpha}}$\n'
                                        f'A = {A_fit:.2e} Ha, α = {alpha_fit:.3f}\n'
                                        f'R² = {r2_power:.6f}')
                        
                        # Add theoretical T^(3/5) line
                        E_theory_t35 = a * T_fine**(3/5)
                        ax4.loglog(T_fine, E_theory_t35, 'k--', alpha=0.5,
                                  label='Theory: $E \\propto T^{3/5}$ (α = 0.6)')
                        
                        ax4.loglog(T_shifted, E_shifted_pos, 'bo', alpha=0.7, markersize=6, 
                                  label='Data (E - E₀)')
                        
                        ax4.set_title(f'Power Law Fit: α = {alpha_fit:.3f} (expect 0.6)')
            
            ax4.set_xlabel('Temperature (K)')
            ax4.set_ylabel('Energy (Ha)')
            ax4.legend(fontsize=9)
            ax4.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            # Save plot
            if output_file is None:
                output_file = f'empirical_fit_{self.energy_component}.png'
            
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"✓ Saved empirical fit plot to {output_file}")
            
            if show_plot:
                plt.show()
            else:
                plt.close()
                
        except ImportError:
            print("Warning: matplotlib not available, cannot create fit plots")
        except Exception as e:
            print(f"Warning: Failed to create fit plot: {e}")
    
    def calculate_systemic_temperature(self, instantaneous_energy_hartree: float, num_particles: int = None) -> float:
        """
        Calculate systemic temperature from instantaneous potential energy.
        
        Uses extended fitting functions when available for higher accuracy:
        - Extended harmonic: E = aT/(1+bT) -> T = E/(a-bE) 
        - Extended T^(3/5): E = E₀ + AT^(3/5)/(1+CT^(3/5)) -> solve numerically
    
    Parameters
    ----------
        instantaneous_energy_hartree : float
            Measured potential energy in hartree
        num_particles : int, optional
            Number of particles (required for harmonic direct calculation)
            
        Returns
    -------
        float
            Systemic temperature in Kelvin
        """
        if self.use_direct_harmonic and self.energy_component == 'harmonic':
            # Direct harmonic calculation: T = 4*E/(N*kB)
            if num_particles is None:
                raise ValueError("num_particles required for direct harmonic calculation")
            
            # Convert energy from hartree to Kelvin: E_hartree * hartree_to_K
            from .utils import PhysicalConstants
            energy_kelvin = instantaneous_energy_hartree * PhysicalConstants.hartree_to_kelvin()
            temperature = 4.0 * energy_kelvin / num_particles  # 4*E/(N*kB), kB=1 in these units
            return max(temperature, 0.0)  # Ensure positive temperature
        
        # Prioritize extended fits for better accuracy
        elif self.has_extended_harmonic_fit and self.energy_component == 'harmonic':
            # Extended harmonic: E = aT/(1+bT) -> T = E/(a-bE)
            a, b = self.extended_harmonic_fit['a'], self.extended_harmonic_fit['b']
            E = instantaneous_energy_hartree
            
            if E <= 0 or (a - b * E) <= 0:
                return 0.0
            
            try:
                temperature = E / (a - b * E)
                return max(temperature, 0.0)
            except (ValueError, ZeroDivisionError):
                return 0.0
        
        elif self.has_extended_t35_fit and self.energy_component != 'harmonic':
            # Extended T^(3/5): E = E₀ + AT^(3/5)/(1+CT^(3/5))
            # Solve numerically using scipy.optimize
            try:
                from scipy.optimize import fsolve
                
                e0, a, c = self.extended_t35_fit['e0'], self.extended_t35_fit['a'], self.extended_t35_fit['c']
                E = instantaneous_energy_hartree
                
                if E <= e0:
                    return 0.0
                
                def equation(T):
                    if T <= 0:
                        return float('inf')
                    t_power = T**(3/5)
                    return e0 + a * t_power / (1 + c * t_power) - E
                
                # Initial guess based on simple T^(3/5) scaling
                T_guess = max(((E - e0) / a) ** (5/3), 1.0) if a > 0 else 300.0
                
                solution = fsolve(equation, T_guess, full_output=True)
                temperature = solution[0][0]
                
                # Check if solution converged
                if solution[2] == 1 and temperature > 0:
                    return temperature
                else:
                    # Extended fit failed to converge, fall back to linear interpolation
                    pass
                    
            except Exception as e:
                print(f"Warning: Extended T^(3/5) temperature calculation failed: {e}")
                # Fall back to linear interpolation
                pass
        
        # Fallback to linear interpolation if no extended fits are available
        # This should rarely happen since we now always try to fit extended functions
        print(f"Warning: No extended fits available for {self.energy_component}, using linear interpolation")
        temperature = np.interp(instantaneous_energy_hartree, self.energies, self.temperatures)
        return max(temperature, 0.0)


class EmpiricalTemperatureFeedback(hoomd.custom.Action):
    """
    Empirical temperature feedback controller for cavity MD simulations.
    
    This class measures instantaneous systemic temperature from potential energy
    using empirical calibration data, averages over time windows, and updates
    both molecular and cavity thermostat temperatures to maintain equilibrium.
    
    The feedback system:
    1. Measures instantaneous LJ+Coulomb energy at each timestep
    2. Converts to systemic temperature using Rosenfeld-fitted relationship
    3. Averages systemic temperature over configurable time windows
    4. Updates thermostat target temperatures at regular intervals
    5. Optionally turns off at specified time to allow free evolution
        
        Parameters
        ----------
    empirical_data : EmpiricalTemperatureData
        Calibrated energy-temperature relationship
    energy_tracker : EnergyTracker
        Energy tracker for reading instantaneous energies
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Molecular thermostat to control (None if not applicable)
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Cavity thermostat to control (None if not applicable)
    apply_to : str
        Which thermostats to update ('molecular', 'cavity', 'both', 'none'). 
        'none' means measurement only with no feedback control. Default: 'both'
    output_period_ps : float
        Period for CSV output in picoseconds. Default: 0.1
    averaging_window_ps : float
        Time window for temperature averaging in picoseconds. Default: 5.0
    update_interval_ps : float
        Interval between thermostat updates in picoseconds. Default: 5.0
    T_min : float
        Minimum allowed temperature in Kelvin. Default: 0.0
    T_max : float, optional
        Maximum allowed temperature in Kelvin. Default: None (no upper limit)
    turn_off_time_ps : float, optional
        Time to turn off feedback (None = never turn off). Default: None
    switch_time_ps : float, optional
        Time when coupling switches on (feedback starts after this). Default: None
    output_file : str
        Output file for feedback CSV data. Default: 'empirical_feedback.csv'
    """
    
    def __init__(self, 
                 empirical_data: EmpiricalTemperatureData,
                 energy_tracker,  # EnergyTracker type hint
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 apply_to: str = 'both',
                 output_period_ps: float = 0.1,
                 averaging_window_ps: float = 5.0,
                 update_interval_ps: float = 5.0,
                 T_min: float = 0.0,
                 T_max: Optional[float] = None,
                 turn_off_time_ps: Optional[float] = None,
                 switch_time_ps: Optional[float] = None,
                 output_file: str = 'empirical_feedback.csv',
                 time_tracker=None,
                 initial_temperature: float = 100.0):
        
        super().__init__()
        
        self.empirical_data = empirical_data
        self.energy_tracker = energy_tracker
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.apply_to = apply_to
        self.output_period_ps = output_period_ps
        self.averaging_window_ps = averaging_window_ps
        self.update_interval_ps = update_interval_ps
        self.T_min = T_min
        self.T_max = T_max
        self.turn_off_time_ps = turn_off_time_ps
        self.switch_time_ps = switch_time_ps
        self.output_file = output_file
        self.time_tracker = time_tracker
        self.initial_temperature = initial_temperature
        
        # Create separate harmonic empirical data for harmonic fictive temperature
        # This follows the approach from plot_temperature_feedback.py but uses available parameters
        try:
            self.harmonic_empirical_data = EmpiricalTemperatureData(
                data_file_path=empirical_data.data_file_path,
                energy_component='harmonic',
                use_direct_harmonic=False,  # Use fitted function (extended model)
                create_plots=True  # Create diagnostic plots for fitting validation
            )
            print(f"Created harmonic empirical data for fictive temperature calculation")
        except Exception as e:
            print(f"Warning: Could not create harmonic empirical data: {e}")
            print(f"   Will use direct equipartition formula for harmonic fictive temperature")
            self.harmonic_empirical_data = None
        
        # Internal state
        self.instantaneous_temperatures = []
        self.measurement_times = []
        self.last_output_time = 0.0
        self.last_update_time = 0.0
        self.feedback_active = True
        self.last_applied_temperature = None
        
        # Validation
        if self.apply_to not in ['molecular', 'cavity', 'both', 'none']:
            raise ValueError("apply_to must be 'molecular', 'cavity', 'both', or 'none'")
        
        if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is None:
            raise ValueError("molecular_thermostat cannot be None when apply_to includes 'molecular'")
        
        if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is None:
            raise ValueError("cavity_thermostat cannot be None when apply_to includes 'cavity'")
        
        # Initialize CSV output file
        self._initialize_output_file()
        
        print("Empirical temperature feedback controller initialized:")
        print(f"   Energy component: {self.empirical_data.energy_component}")
        print(f"   Apply to: {self.apply_to}")
        print(f"   Averaging window: {self.averaging_window_ps:.1f} ps")
        print(f"   Update interval: {self.update_interval_ps:.1f} ps")
        print(f"   Temperature range: [{self.T_min:.1f}, {self.T_max:.1f}] K")
        if self.switch_time_ps:
            print(f"   Switch time: {self.switch_time_ps:.1f} ps (feedback starts after)")
        if self.turn_off_time_ps:
            print(f"   Turn-off time: {self.turn_off_time_ps:.1f} ps")
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        with open(self.output_file, 'w') as f:
            f.write("# Empirical Temperature Feedback Data\n")
            f.write("# time_ps: Simulation time in picoseconds\n")
            f.write("# real_energy_hartree: Instantaneous potential energy in hartree\n")
            f.write("# harmonic_energy_hartree: Expected harmonic oscillator energy in hartree\n")
            f.write("# instantaneous_T_fictive_K: Instantaneous systemic temperature from empirical fit in Kelvin\n")
            f.write("# harmonic_T_fictive_K: Fictive temperature if energy were purely harmonic in Kelvin\n")
            f.write("# averaged_T_fictive_K: Time-averaged systemic temperature in Kelvin\n")
            f.write("# applied_T_K: Target temperature applied to thermostats in Kelvin\n")
            f.write("# feedback_active: Whether feedback is currently active (0/1)\n")
            f.write("time_ps,real_energy_hartree,harmonic_energy_hartree,instantaneous_T_fictive_K,harmonic_T_fictive_K,averaged_T_fictive_K,applied_T_K,feedback_active\n")

    def act(self, timestep):
        """Execute feedback control at each timestep."""
        from .utils import PhysicalConstants
        
        # Get current simulation time using time tracker (FIX: was using wrong calculation)
        if hasattr(self, 'time_tracker') and self.time_tracker is not None:
            current_time_ps = self.time_tracker.elapsed_time
        else:
            # Fallback: this is still wrong but better than before
            # timestep is the step NUMBER, not the timestep size!
            dt_ps = PhysicalConstants.atomic_units_to_ps(1.0)  # Assume 1 au timestep
            current_time_ps = timestep * dt_ps
            print(f"Warning: Using fallback time calculation - time may be inaccurate")
        
        # Always track and output data, but only apply feedback control after switch time
        feedback_should_be_active = True
        if self.switch_time_ps is not None and current_time_ps < self.switch_time_ps:
            feedback_should_be_active = False  # Track data but don't apply feedback yet
        
        # Check if we should turn off feedback
        if self.turn_off_time_ps is not None and current_time_ps >= self.turn_off_time_ps:
            if self.feedback_active:
                self.feedback_active = False
                print(f"🔴 Empirical feedback turned OFF at t = {current_time_ps:.2f} ps")
        
        # Get instantaneous energy from energy tracker
        try:
            if hasattr(self.energy_tracker, 'get_instantaneous_energy'):
                energy_data = self.energy_tracker.get_instantaneous_energy()
                
                if self.empirical_data.energy_component == 'lj_coulombic':
                    if 'lj' in energy_data and 'coulombic' in energy_data:
                        instantaneous_energy = energy_data['lj'] + energy_data['coulombic']
                    else:
                        print("Warning: LJ+Coulombic energy not available")
                        return
                elif self.empirical_data.energy_component == 'total_PE':
                    if 'total_potential' in energy_data:
                        instantaneous_energy = energy_data['total_potential']
                    else:
                        print("Warning: Total potential energy not available")
                        return
                else:
                    print(f"Warning: Energy component {self.empirical_data.energy_component} not supported")
                    return
            else:
                # Fallback: try to access energy directly (this may not work)
                print("Warning: Energy tracker doesn't have get_instantaneous_energy method")
                return
                
        except Exception as e:
            print(f"Warning: Failed to get energy data: {e}")
            return
        
        # Calculate instantaneous systemic temperature
        instantaneous_T = self.empirical_data.calculate_systemic_temperature(instantaneous_energy)
        
        # Calculate expected harmonic oscillator energy using equipartition theorem
        # E_harmonic = 0.25 * N * k_B * T for a harmonic solid with N particles
        N_particles = 500  # Number of molecular particles (excluding cavity particle)
        harmonic_energy = 0.25 * N_particles * PhysicalConstants.KB_HARTREE_PER_K * instantaneous_T
        
        # Calculate harmonic fictive temperature using T^(3/2) fitted function approach
        # This follows the exact approach from plot_temperature_feedback.py calculate_empirical_harmonic_temperature
        try:
            # Try to get actual harmonic energy from the energy tracker
            if hasattr(self, 'energy_tracker') and self.energy_tracker is not None:
                harmonic_energy_dict = self.energy_tracker.get_instantaneous_energy()
                actual_harmonic_energy = harmonic_energy_dict.get('harmonic', None)
                
                if actual_harmonic_energy is not None and actual_harmonic_energy > 0:
                    # Use T^(3/2) fitted function for harmonic energy if we have empirical data
                    if hasattr(self, 'harmonic_empirical_data') and self.harmonic_empirical_data is not None:
                        harmonic_fictive_T = self.harmonic_empirical_data.calculate_systemic_temperature(
                            instantaneous_energy_hartree=actual_harmonic_energy,
                            num_particles=N_particles
                        )
                    else:
                        # Fallback: use direct equipartition formula from plot_temperature_feedback.py
                        kb_hartree_per_k = 3.16681e-6  # Hartree/K (from plot_temperature_feedback.py)
                        harmonic_fictive_T = 4.0 * actual_harmonic_energy / (N_particles * kb_hartree_per_k)
                else:
                    # Fallback: use equipartition temperature
                    harmonic_fictive_T = instantaneous_T
            else:
                # No energy tracker: use equipartition temperature
                harmonic_fictive_T = instantaneous_T
                
        except Exception as e:
            # Any error: use equipartition temperature
            harmonic_fictive_T = instantaneous_T
        
        # Store measurement
        self.instantaneous_temperatures.append(instantaneous_T)
        self.measurement_times.append(current_time_ps)
        
        # Remove old data outside averaging window
        cutoff_time = current_time_ps - self.averaging_window_ps
        while len(self.measurement_times) > 0 and self.measurement_times[0] < cutoff_time:
            self.measurement_times.pop(0)
            self.instantaneous_temperatures.pop(0)
        
        # Calculate averaged temperature
        if len(self.instantaneous_temperatures) > 0:
            averaged_T = np.mean(self.instantaneous_temperatures)
        else:
            averaged_T = instantaneous_T
        
        # Update thermostat temperatures if it's time
        # Default to initial temperature when no feedback has been applied
        applied_T = self.last_applied_temperature if self.last_applied_temperature is not None else self.initial_temperature
        
        if feedback_should_be_active and self.feedback_active and (current_time_ps - self.last_update_time) >= self.update_interval_ps:
            # Clamp averaged temperature to allowed range
            target_T = np.clip(averaged_T, self.T_min, self.T_max)
                
            # Update thermostats (this may do nothing if apply_to='none')
            self._update_thermostat_temperatures(target_T)
            
            self.last_update_time = current_time_ps
            
            # Only update applied_T if we're actually applying feedback
            if self.apply_to != 'none':
                self.last_applied_temperature = target_T
                applied_T = target_T
        
        # Output to CSV file
        if (current_time_ps - self.last_output_time) >= self.output_period_ps:
            with open(self.output_file, 'a') as f:
                # Show feedback as active only when it should actually be controlling
                feedback_status = int(feedback_should_be_active and self.feedback_active)
                f.write(f"{current_time_ps:.6f},{instantaneous_energy:.8f},{harmonic_energy:.8f},"
                       f"{instantaneous_T:.6f},{harmonic_fictive_T:.6f},{averaged_T:.6f},{applied_T:.6f},{feedback_status}\n")
            
            self.last_output_time = current_time_ps

    def _update_thermostat_temperatures(self, target_temperature_K: float):
        """Update thermostat target temperatures."""
        from .utils import PhysicalConstants
        
        # Convert temperature to HOOMD energy units (kT in atomic units)
        target_kT = target_temperature_K * PhysicalConstants.KB_HARTREE_PER_K
        
        # Only update thermostats if apply_to is not 'none'
        if self.apply_to == 'none':
            return  # Measurement only, no feedback
            
        try:
            updated_any = False
            
            if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
                if hasattr(self.molecular_thermostat, 'kT'):
                    # Direct kT assignment (may not work if it's a Variant)
                    self.molecular_thermostat.kT = target_kT
                    updated_any = True
                elif hasattr(self.molecular_thermostat, '_kT'):
                    # Try private attribute
                    self.molecular_thermostat._kT = target_kT
                    updated_any = True
                else:
                    print("Warning: Cannot update molecular thermostat temperature")
            
            if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
                if hasattr(self.cavity_thermostat, 'kT'):
                    self.cavity_thermostat.kT = target_kT
                    updated_any = True
                elif hasattr(self.cavity_thermostat, '_kT'):
                    self.cavity_thermostat._kT = target_kT
                    updated_any = True
                else:
                    print("Warning: Cannot update cavity thermostat temperature")
            
            # Only print if we actually updated something
            if updated_any:
                print(f"Updated thermostat temperatures to {target_temperature_K:.2f} K")
            
        except Exception as e:
            print(f"Warning: Failed to update thermostat temperatures: {e}")


class TemperatureTracker(hoomd.custom.Action):
    """
    Comprehensive temperature tracker for cavity MD simulations.
    
    Tracks and outputs all relevant temperature measurements:
    1. Kinetic temperature (from particle velocities)
    2. Harmonic fictive temperature (from harmonic energy via empirical data)
    3. LJ+Coulombic fictive temperature (from LJ+Coul energy via empirical data)
    4. Cavity bath temperature (thermostat kT)
    5. Molecular bath temperature (thermostat kT)
    6. Harmonic equipartition temperature (from harmonic energy via equipartition theorem)
        
        Parameters
        ----------
    simulation : hoomd.Simulation
        HOOMD simulation object
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    output_period_ps : float
        Output period in picoseconds
    output_file : str
        Output CSV file path
    energy_tracker : EnergyTracker, optional
        Energy tracker for energy-based temperatures
    molecular_thermostat : hoomd thermostat, optional
        Molecular thermostat object
    cavity_thermostat : hoomd thermostat, optional
        Cavity thermostat object
    empirical_data_file : str, optional
        Path to empirical data file for LJ+Coul fictive temperature
    """
    
    def __init__(self, 
                 simulation: hoomd.Simulation,
                 time_tracker: ElapsedTimeTracker,
                 output_period_ps: float,
                 output_file: str,
                 energy_tracker=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 empirical_data_file=None,
                 debug=False):
        
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.output_period_ps = output_period_ps
        self.output_file = output_file
        self.energy_tracker = energy_tracker
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.debug = debug
        
        # Load empirical data for fictive temperature calculations if provided
        self.empirical_data_harmonic = None
        self.empirical_data_lj_coul = None
        if empirical_data_file is not None:
            try:
                # Load empirical data for harmonic energy component
                self.empirical_data_harmonic = EmpiricalTemperatureData(
                    data_file_path=empirical_data_file,
                    energy_component='harmonic',
                    use_direct_harmonic=False,  # Use fitted function, not direct calculation
                    create_plots=True  # Create diagnostic plots
                )
                
                # Load empirical data for LJ+Coulombic energy component (uses extended T^(3/5) scaling)
                self.empirical_data_lj_coul = EmpiricalTemperatureData(
                    data_file_path=empirical_data_file,
                    energy_component='lj_coulombic',
                    create_plots=True  # Create diagnostic plots
                )
            except Exception as e:
                print(f"Warning: Could not load empirical data file {empirical_data_file}: {e}")
                print("Empirical fictive temperatures will not be available")
        
        # Tracking state
        self.last_output_time = None
        
        # Initialize temperature attributes for external access (e.g., AutoStopController)
        self.kinetic_temperature = None
        self.harmonic_fictive_temperature = None
        self.lj_coulombic_fictive_temperature = None
        self.cavity_bath_temperature = None
        self.molecular_bath_temperature = None
        self.harmonic_equipartition_temperature = None
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        with open(self.output_file, 'w') as f:
            f.write("# Comprehensive Temperature Tracking\n")
            f.write("# time_ps: simulation time in picoseconds\n")
            f.write("# kinetic_temp_K: kinetic temperature from particle velocities\n")
            f.write("# harmonic_fictive_K: fictive temperature from harmonic energy\n")
            f.write("# lj_coul_fictive_K: fictive temperature from LJ+Coulombic energy\n")
            f.write("# cavity_bath_K: cavity thermostat temperature\n")
            f.write("# molecular_bath_K: molecular thermostat temperature\n")
            f.write("# harmonic_equipartition_K: harmonic fictive temperature from equipartition theorem\n")
            f.write("time_ps,kinetic_temp_K,harmonic_fictive_K,lj_coul_fictive_K,cavity_bath_K,molecular_bath_K,harmonic_equipartition_K\n")
    
    def act(self, timestep):
        """Track all temperatures at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Always calculate temperatures for external access (regardless of output timing)
        # 1. Calculate kinetic temperature
        self.kinetic_temperature = self._calculate_kinetic_temperature()
        
        # 2. Calculate harmonic fictive temperature
        self.harmonic_fictive_temperature = self._calculate_harmonic_fictive_temperature()
        
        # 3. Calculate LJ+Coulombic fictive temperature
        self.lj_coulombic_fictive_temperature = self._calculate_lj_coul_fictive_temperature()
        
        # 4. Get cavity bath temperature
        self.cavity_bath_temperature = self._get_cavity_bath_temperature()
        
        # 5. Get molecular bath temperature
        self.molecular_bath_temperature = self._get_molecular_bath_temperature()
        
        # 6. Calculate harmonic equipartition temperature
        self.harmonic_equipartition_temperature = self._calculate_harmonic_equipartition_temperature()
        
        # Write to CSV only when needed
        if self._should_output(current_time_ps):
            with open(self.output_file, 'a') as f:
                f.write(f"{current_time_ps:.6f},{self.kinetic_temperature:.6f},{self.harmonic_fictive_temperature:.6f},"
                       f"{self.lj_coulombic_fictive_temperature:.6f},{self.cavity_bath_temperature:.6f},{self.molecular_bath_temperature:.6f},"
                       f"{self.harmonic_equipartition_temperature:.6f}\n")
            
            self.last_output_time = current_time_ps
    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output data."""
        if self.last_output_time is None:
            return True
        return (current_time_ps - self.last_output_time) >= self.output_period_ps
    
    def _get_hoomd_simulation(self):
        """Get the HOOMD simulation object, handling both CavityMDSimulation and direct HOOMD objects."""
        if hasattr(self.simulation, 'sim'):
            # CavityMDSimulation object - get HOOMD simulation
            return self.simulation.sim
        else:
            # Direct HOOMD simulation object
            return self.simulation
    
    def _calculate_kinetic_temperature(self) -> float:
        """Calculate kinetic temperature from particle velocities."""
        try:
            from .utils import PhysicalConstants
            
            with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                velocities = np.array(snap.particles.velocity)
                masses = np.array(snap.particles.mass)
                N_particles = len(masses)
            
            if N_particles == 0:
                if self.debug:
                    print("DEBUG: _calculate_kinetic_temperature: N_particles = 0")
                return 0.0
            
            # Calculate kinetic energy
            kinetic_energy = 0.5 * np.sum(masses[:, np.newaxis] * velocities**2)
            
            # Temperature from equipartition: (3/2)*N*kB*T = KE
            # T = (2*KE)/(3*N*kB)
            kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
            temperature_K = (2.0 * kinetic_energy) / (3.0 * N_particles * kB_hartree)
            
            if self.debug:
                print(f"DEBUG: _calculate_kinetic_temperature: N={N_particles}, KE={kinetic_energy:.6e}, T={temperature_K:.2f}K")
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not calculate kinetic temperature: {e}")
            import traceback
            traceback.print_exc()
            return 0.0
    
    def _calculate_harmonic_fictive_temperature(self) -> float:
        """Calculate harmonic fictive temperature using empirical data."""
        try:
            if self.energy_tracker is None:
                return 0.0
            
            # Get harmonic energy from energy tracker
            energy_dict = self.energy_tracker.get_instantaneous_energy()
            harmonic_energy = energy_dict.get('harmonic', 0.0)
            
            if harmonic_energy <= 0:
                return 0.0
            
            # Use empirical data if available (preferred method)
            if self.empirical_data_harmonic is not None:
                temperature_K = self.empirical_data_harmonic.calculate_systemic_temperature(harmonic_energy)
                return temperature_K
            
            # Fallback to direct calculation if no empirical data
            from .utils import PhysicalConstants
            with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                N_particles = len(snap.particles.mass)
            
            if N_particles == 0:
                return 0.0
            
            # Direct harmonic calculation: T = 4*E/(N*kB) for 3D harmonic oscillator
            # Use the exact same constant as in plot_temperature_feedback.py
            kB_hartree = 3.1668105e-6  # Hartree/K (Boltzmann constant)
            temperature_K = (4.0 * harmonic_energy) / (N_particles * kB_hartree)
            
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not calculate harmonic fictive temperature: {e}")
            return 0.0
    
    def _calculate_harmonic_equipartition_temperature(self) -> float:
        """Calculate harmonic fictive temperature using equipartition theorem only.
        
        This provides a direct equipartition-based estimate: T = 4*E_harmonic/(N*kB)
        for a 3D harmonic oscillator system. This becomes more accurate at low temperatures
        where the equipartition theorem is more valid.
        
        Returns
        -------
        float
            Harmonic equipartition temperature in Kelvin
        """
        try:
            if self.energy_tracker is None:
                if self.debug:
                    print("DEBUG: _calculate_harmonic_equipartition_temperature: energy_tracker is None")
                return 0.0
            
            # Get harmonic energy from energy tracker
            energy_dict = self.energy_tracker.get_instantaneous_energy()
            harmonic_energy = energy_dict.get('harmonic', 0.0)
            
            if harmonic_energy <= 0:
                if self.debug:
                    print(f"DEBUG: _calculate_harmonic_equipartition_temperature: harmonic_energy <= 0 ({harmonic_energy})")
                return 0.0
            
            # Get number of particles
            from .utils import PhysicalConstants
            with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                N_particles = len(snap.particles.mass)
            
            if N_particles == 0:
                if self.debug:
                    print("DEBUG: _calculate_harmonic_equipartition_temperature: N_particles = 0")
                return 0.0
            
            # Direct equipartition calculation: T = 4*E/(N*kB) for 3D harmonic oscillator
            # For a classical 3D harmonic oscillator: <E> = (3/2)*N*kB*T for kinetic + (3/2)*N*kB*T for potential
            # But here we only have the harmonic potential energy, so: <E_harmonic> = (3/2)*N*kB*T
            # Therefore: T = (2*E_harmonic)/(3*N*kB)
            # However, based on the empirical observations, the factor is 4 instead of 2/3
            kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
            temperature_K = (4.0 * harmonic_energy) / (N_particles * kB_hartree)
            
            if self.debug:
                print(f"DEBUG: _calculate_harmonic_equipartition_temperature: N={N_particles}, E_harm={harmonic_energy:.6e}, T={temperature_K:.2f}K")
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not calculate harmonic equipartition temperature: {e}")
            import traceback
            traceback.print_exc()
            return 0.0
    
    def _calculate_lj_coul_fictive_temperature(self) -> float:
        """Calculate LJ+Coulombic fictive temperature using empirical data."""
        try:
            if self.empirical_data_lj_coul is None or self.energy_tracker is None:
                return 0.0
            
            # Get LJ+Coulombic energy components
            energy_dict = self.energy_tracker.get_instantaneous_energy()
            lj_energy = energy_dict.get('lj', 0.0)
            coulombic_energy = energy_dict.get('coulombic', 0.0)
            total_lj_coul = lj_energy + coulombic_energy
            
            
            if total_lj_coul == 0:
                return 0.0
            
            # Use empirical data to convert LJ+Coul energy to systemic temperature
            temperature_K = self.empirical_data_lj_coul.calculate_systemic_temperature(total_lj_coul)
            
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not calculate LJ+Coul fictive temperature: {e}")
            return 0.0
    
    def _get_cavity_bath_temperature(self) -> float:
        """Get cavity thermostat temperature."""
        try:
            if self.cavity_thermostat is None:
                return 0.0
            
            from .utils import PhysicalConstants
            
            # Get kT from thermostat - handle different thermostat types
            kT = None
            if hasattr(self.cavity_thermostat, 'kT'):
                kT_value = self.cavity_thermostat.kT
                # Handle both Constant variants and plain floats
                if hasattr(kT_value, 'value'):
                    kT = kT_value.value
                elif hasattr(kT_value, '__call__'):
                    kT = kT_value(0)  # Evaluate at timestep 0
                else:
                    kT = float(kT_value)
            elif hasattr(self.cavity_thermostat, 'T'):
                # Some thermostats use T instead of kT
                T_value = self.cavity_thermostat.T
                if hasattr(T_value, 'value'):
                    temperature_K = T_value.value
                elif hasattr(T_value, '__call__'):
                    temperature_K = T_value(0)
                else:
                    temperature_K = float(T_value)
                return temperature_K
            elif hasattr(self.cavity_thermostat, 'thermostat'):
                # For ConstantVolume methods, check the thermostat attribute
                nested_thermo = self.cavity_thermostat.thermostat
                if hasattr(nested_thermo, 'kT'):
                    kT_value = nested_thermo.kT
                    if hasattr(kT_value, 'value'):
                        kT = kT_value.value
                    elif hasattr(kT_value, '__call__'):
                        kT = kT_value(0)
                    else:
                        kT = float(kT_value)
            
            if kT is None:
                return 0.0
                
            temperature_K = kT / PhysicalConstants.KB_HARTREE_PER_K
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not get cavity bath temperature: {e}")
            return 0.0
    
    def _get_molecular_bath_temperature(self) -> float:
        """Get molecular thermostat temperature."""
        try:
            if self.molecular_thermostat is None:
                return 0.0
            
            from .utils import PhysicalConstants
            
            # Get kT from thermostat - handle different thermostat types
            kT = None
            if hasattr(self.molecular_thermostat, 'kT'):
                kT_value = self.molecular_thermostat.kT
                # Handle both Constant variants and plain floats
                if hasattr(kT_value, 'value'):
                    kT = kT_value.value
                elif hasattr(kT_value, '__call__'):
                    kT = kT_value(0)  # Evaluate at timestep 0
                else:
                    kT = float(kT_value)
            elif hasattr(self.molecular_thermostat, 'T'):
                # Some thermostats use T instead of kT
                T_value = self.molecular_thermostat.T
                if hasattr(T_value, 'value'):
                    temperature_K = T_value.value
                elif hasattr(T_value, '__call__'):
                    temperature_K = T_value(0)
                else:
                    temperature_K = float(T_value)
                return temperature_K
            elif hasattr(self.molecular_thermostat, 'thermostat'):
                # For ConstantVolume methods, check the thermostat attribute
                nested_thermo = self.molecular_thermostat.thermostat
                if hasattr(nested_thermo, 'kT'):
                    kT_value = nested_thermo.kT
                    if hasattr(kT_value, 'value'):
                        kT = kT_value.value
                    elif hasattr(kT_value, '__call__'):
                        kT = kT_value(0)
                    else:
                        kT = float(kT_value)
            
            if kT is None:
                return 0.0
                
            temperature_K = kT / PhysicalConstants.KB_HARTREE_PER_K
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not get molecular bath temperature: {e}")
            return 0.0


class GradientDescentTemperatureFeedback(hoomd.custom.Action):
    """
    Gradient descent temperature feedback controller for cavity MD simulations.
    
    This class implements a discretized gradient descent algorithm to minimize
    the objective function J = 1/2 * (T_eff(t) - T_target)^2, where:
    - T_eff(t) = (T_measured(t) + T_bath(t)) / 2 is the effective system temperature
    - T_target is the target temperature
    
    The gradient descent update rule is:
    T_bath(t+1) = T_bath(t) - α * ∇J/∇T_bath
    T_bath(t+1) = T_bath(t) - α * 0.5 * (T_eff(t) - T_target)
    
    where α is the learning rate determined by the time constant.
    
    Parameters
    ----------
    temperature_method : str
        Temperature calculation method ('kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition')
    time_constant_ps : float
        Time constant for gradient descent in picoseconds (controls convergence speed)
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    energy_tracker : EnergyTracker  
        Energy tracker for temperature calculations
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Molecular thermostat to control
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Cavity thermostat to control
    target_temperature : float, optional
        Target temperature in Kelvin. Default: 100.0
    apply_to : str, optional
        Which thermostats to control ('molecular', 'cavity', 'both'). Default: 'both'
    turn_on_time_ps : float, optional
        Time to start gradient descent control in picoseconds. Default: 0.0
    turn_off_time_ps : float, optional
        Time to stop gradient descent control in picoseconds. Default: None (never turn off)
    update_interval_ps : float, optional
        Interval between control updates in picoseconds. Default: 0.1 (every MD timestep)
    T_min : float, optional
        Minimum allowed bath temperature in Kelvin. Default: 0.0
    T_max : float, optional
        Maximum allowed bath temperature in Kelvin. Default: None (no upper limit)
    rate_limit_K_per_ps : float, optional
        Maximum rate of temperature change in K/ps. Default: None (no rate limit)
    output_file : str, optional
        Output file for gradient descent control data. Default: 'gd_feedback.csv'
    empirical_data_file : str, optional
        Path to empirical data file (required for 'lj_coulombic' and 'harmonic' methods)
        
    Attributes
    ----------
    alpha : float
        Learning rate (computed from time constant and MD timestep)
    current_bath_temperature : float
        Current bath temperature setting
    is_active : bool
        Whether gradient descent control is currently active
        
    Examples
    --------
    **Basic kinetic temperature gradient descent control:**
    
    >>> from hoomd.cavitymd.analysis import GradientDescentTemperatureFeedback
    >>> 
    >>> gd_controller = GradientDescentTemperatureFeedback(
    ...     temperature_method='kinetic',
    ...     time_constant_ps=10.0,  # 10 ps time constant
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     target_temperature=100.0
    ... )
    
    **Harmonic fictive temperature control with fast convergence:**
    
    >>> gd_controller = GradientDescentTemperatureFeedback(
    ...     temperature_method='harmonic',
    ...     time_constant_ps=5.0,  # Fast convergence
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     target_temperature=80.0,
    ...     empirical_data_file='harmonic_calibration.txt',
    ...     apply_to='both'
    ... )
    
    **Harmonic equipartition temperature control (no empirical data needed):**
    
    >>> gd_controller = GradientDescentTemperatureFeedback(
    ...     temperature_method='harmonic_equipartition',
    ...     time_constant_ps=3.0,  # Fast response
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     target_temperature=100.0,
    ...     apply_to='molecular'
    ... )
            
        Notes
        -----
    - Models system-bath coupling as T_eff = (T_measured + T_bath) / 2
    - Smaller time constants lead to faster but potentially less stable convergence
    - Learning rate α = dt / τ where dt is MD timestep and τ is time constant
    - Gradient descent is simpler and more intuitive than PI control
    - Compatible with all existing temperature calculation methods
    """
    
    def __init__(self, 
                 temperature_method: str,
                 time_constant_ps: float,
                 time_tracker,
                 energy_tracker,
                 simulation=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 target_temperature: float = 100.0,
                 apply_to: str = 'both',
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 T_min: float = 0.0,
                 T_max: Optional[float] = None,
                 rate_limit_K_per_ps: Optional[float] = None,
                 output_file: str = 'gd_feedback.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0):
        
        super().__init__()
        
        # Store core parameters
        self.temperature_method = temperature_method
        self.time_constant_ps = float(time_constant_ps)
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        
        # Control parameters
        self.target_temperature = float(target_temperature)
        self.apply_to = apply_to
        self.turn_on_time_ps = float(turn_on_time_ps)
        self.turn_off_time_ps = float(turn_off_time_ps) if turn_off_time_ps is not None else None
        self.update_interval_ps = float(update_interval_ps)
        self.T_min = float(T_min)
        self.T_max = float(T_max) if T_max is not None else None
        self.rate_limit_K_per_ps = float(rate_limit_K_per_ps) if rate_limit_K_per_ps is not None else None
        self.output_file = output_file
        self.console_output_period_ps = float(console_output_period_ps)
        
        # Compute learning rate: α = dt / τ 
        # Assuming MD timestep is typically 0.1 fs = 0.0001 ps
        md_timestep_ps = 0.0001  # This could be made configurable
        self.alpha = md_timestep_ps / self.time_constant_ps
        
        # Initialize state
        self.current_bath_temperature = None  # Will be set from thermostats when active
        self.last_update_time = 0.0
        self.last_console_output_time = 0.0
        self.is_active = False
        
        # Set up empirical data if needed
        self.empirical_data = None
        if empirical_data_file and self.temperature_method in ['lj_coulombic', 'harmonic', 'lj_coulombic_kinetic', 'lj_coulombic_bath']:
            try:
                if self.temperature_method == 'harmonic':
                    self.empirical_data = EmpiricalTemperatureData(
                        empirical_data_file, energy_component='harmonic', create_plots=True
                    )
                else:  # lj_coulombic, lj_coulombic_kinetic, lj_coulombic_bath
                    self.empirical_data = EmpiricalTemperatureData(
                        empirical_data_file, energy_component='lj_coulombic', create_plots=True
                    )
                print(f"Loaded empirical data for {self.temperature_method} temperature calculation")
            except Exception as e:
                print(f"Warning: Failed to load empirical data: {e}")
        
        # Validate configuration
        self._validate_configuration()
        
        # Initialize output file
        self._initialize_output_file()
        
        # Print configuration
        print(f"Gradient Descent Temperature Controller initialized:")
        print(f"   Temperature method: {self.temperature_method}")
        print(f"   Time constant: {self.time_constant_ps:.1f} ps")
        print(f"   Learning rate α: {self.alpha:.6f}")
        print(f"   Target temperature: {self.target_temperature:.1f} K")
        print(f"   Apply to: {self.apply_to}")
        print(f"   Update interval: {self.update_interval_ps:.1f} ps")
        T_max_str = f"{self.T_max:.1f}" if self.T_max is not None else "∞"
        print(f"   Temperature limits: [{self.T_min:.1f}, {T_max_str}] K")
        if self.rate_limit_K_per_ps is not None:
            print(f"   Rate limit: {self.rate_limit_K_per_ps:.3f} K/ps")
        if self.turn_on_time_ps > 0:
            print(f"   Turn-on time: {self.turn_on_time_ps:.1f} ps")
        if self.turn_off_time_ps is not None:
            print(f"   Turn-off time: {self.turn_off_time_ps:.1f} ps")
    
    def _validate_configuration(self):
        """Validate controller configuration."""
        if self.temperature_method not in ['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic', 'lj_coulombic_bath']:
            raise ValueError(f"Invalid temperature_method: {self.temperature_method}")
        
        if self.apply_to not in ['molecular', 'cavity', 'both']:
            raise ValueError(f"apply_to must be 'molecular', 'cavity', or 'both', got {self.apply_to}")
        
        if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is None:
            raise ValueError("molecular_thermostat cannot be None when apply_to includes 'molecular'")
        
        if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is None:
            raise ValueError("cavity_thermostat cannot be None when apply_to includes 'cavity'")
        
        if self.time_constant_ps <= 0:
            raise ValueError("time_constant_ps must be positive")
        
        if self.target_temperature <= 0:
            raise ValueError("target_temperature must be positive")
        
        # Validate reasonable parameter ranges
        if self.time_constant_ps > 10000:
            print(f"Warning: Very large time constant ({self.time_constant_ps:.1f} ps) will result in extremely slow convergence")
            print(f"  Learning rate α = {self.alpha:.2e} - consider using a smaller time constant")
        
        if self.update_interval_ps < 0.01:
            print(f"Warning: Very small update interval ({self.update_interval_ps:.3f} ps) may cause excessive computational overhead")
        
        if self.alpha < 1e-8:
            print(f"Warning: Learning rate α = {self.alpha:.2e} is extremely small - controller may not respond")
            print(f"  Consider reducing time constant from {self.time_constant_ps:.1f} ps")
        
        if self.alpha > 0.1:
            print(f"Warning: Learning rate α = {self.alpha:.2e} is very large - controller may be unstable")
            print(f"  Consider increasing time constant from {self.time_constant_ps:.1f} ps")
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("# Gradient Descent Temperature Feedback Control Data\n")
                f.write("# time_ps: Simulation time in picoseconds\n")
                f.write("# measured_T_K: Measured system temperature in Kelvin\n")
                f.write("# effective_T_K: Effective temperature = (measured + bath) / 2 in Kelvin\n")
                f.write("# target_T_K: Target (setpoint) temperature in Kelvin\n")
                f.write("# error_K: Temperature error (effective - target) in Kelvin\n")
                f.write("# gradient: Computed gradient dJ/dT_bath = 0.5 * error\n")
                f.write("# raw_output_K: Raw gradient descent output in Kelvin\n")
                f.write("# saturated_output_K: Saturated output (within limits) in Kelvin\n")
                f.write("# rate_limited_output_K: Final output after rate limiting in Kelvin\n")
                f.write("# is_active: Whether gradient descent control is active (0/1)\n")
                f.write("time_ps,measured_T_K,effective_T_K,target_T_K,error_K,gradient,raw_output_K,saturated_output_K,rate_limited_output_K,is_active\n")
        except Exception as e:
            print(f"Warning: Failed to initialize gradient descent output file: {e}")
    
    def _calculate_system_temperature(self, current_time_ps: float) -> Optional[float]:
        """Calculate system temperature using the specified method."""
        try:
            if self.temperature_method == 'kinetic':
                # Calculate kinetic temperature directly from particle velocities
                # (same method as TemperatureTracker._calculate_kinetic_temperature)
                import numpy as np
                from .utils import PhysicalConstants
                
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    velocities = np.array(snap.particles.velocity)
                    masses = np.array(snap.particles.mass)
                    N_particles = len(masses)
                
                if N_particles == 0:
                    return 0.0
                
                # Calculate kinetic energy
                kinetic_energy = 0.5 * np.sum(masses[:, np.newaxis] * velocities**2)
                
                # Temperature from equipartition: (3/2)*N*kB*T = KE
                # T = (2*KE)/(3*N*kB)
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (2.0 * kinetic_energy) / (3.0 * N_particles * kB_hartree)
                return temperature_K
                
            elif self.temperature_method == 'lj_coulombic':
                if self.empirical_data is None:
                    return None
                
                # Get LJ + Coulombic energy
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                lj_energy = energy_dict.get('lj', 0.0)
                coulomb_energy = energy_dict.get('coulombic', 0.0)
                total_energy = lj_energy + coulomb_energy
                
                # Convert to temperature using empirical data
                temperature = self.empirical_data.calculate_systemic_temperature(total_energy)
                return max(temperature, 0.0)  # Ensure non-negative temperature
                
            elif self.temperature_method == 'lj_coulombic_kinetic':
                # Dual temperature control: use maximum error signal from either LJ+Coulombic or kinetic
                # Error signal = sign * max(|T_lj_coulombic - T_target|, |T_kinetic - T_target|)
                
                # Calculate kinetic temperature
                import numpy as np
                from .utils import PhysicalConstants
                
                try:
                    # Get kinetic energy from molecular particles
                    kinetic_energy = 0.0
                    N_molecules = 0
                    with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                        for i, ptype in enumerate(snap.particles.typeid):
                            # Exclude cavity particle (type index depends on system setup)
                            if i < len(snap.particles.mass) - 1:  # Assume cavity is last particle
                                mass = snap.particles.mass[i]
                                velocity = snap.particles.velocity[i]
                                # KE = (1/2) * m * v^2
                                ke_particle = 0.5 * mass * np.sum(velocity**2)
                                kinetic_energy += ke_particle
                                N_molecules += 1
                    
                    # Convert to temperature: T = (2/3) * KE / (N * kB)
                    if N_molecules > 0:
                        kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                        T_kinetic_K = (2.0/3.0) * kinetic_energy / (N_molecules * kB_hartree)
                    else:
                        return None
                        
                except Exception as e:
                    print(f"Warning: Failed to calculate kinetic temperature: {e}")
                    return None
                
                # Calculate LJ+Coulombic fictive temperature from energy data
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if self.empirical_data is None:
                    print(f"Warning: No empirical data available for LJ+Coulombic temperature")
                    return None
                
                if 'lj' not in energy_data or 'coulombic' not in energy_data:
                    print(f"Warning: LJ+Coulombic energy not available")
                    return None
                
                lj_coulombic_energy = energy_data['lj'] + energy_data['coulombic']
                T_lj_coulombic_K = self.empirical_data.calculate_systemic_temperature(lj_coulombic_energy)
                if T_lj_coulombic_K is None:
                    return None
                T_lj_coulombic_K = max(T_lj_coulombic_K, 0.0)
                
                # Calculate signed errors for directional control
                error_kinetic_signed = T_kinetic_K - self.target_temperature
                error_lj_coulombic_signed = T_lj_coulombic_K - self.target_temperature
                
                # Use directionally-aware selection: choose temperature requiring larger control action
                error_kinetic_abs = abs(error_kinetic_signed)
                error_lj_coulombic_abs = abs(error_lj_coulombic_signed)
                
                # Selection logic: choose temperature with larger absolute error
                # In case of tie (within 0.1 K), prioritize LJ+Coulombic (more fundamental)
                tolerance = 0.1  # K
                if abs(error_kinetic_abs - error_lj_coulombic_abs) <= tolerance:
                    selected_temperature = T_lj_coulombic_K
                    selected_method = "lj_coulombic"
                    reason = f"tie-break (Δerr={error_kinetic_abs-error_lj_coulombic_abs:.2f}K)"
                elif error_kinetic_abs > error_lj_coulombic_abs:
                    selected_temperature = T_kinetic_K
                    selected_method = "kinetic"
                    reason = f"larger error ({error_kinetic_abs:.1f}K > {error_lj_coulombic_abs:.1f}K)"
                else:
                    selected_temperature = T_lj_coulombic_K
                    selected_method = "lj_coulombic" 
                    reason = f"larger error ({error_lj_coulombic_abs:.1f}K > {error_kinetic_abs:.1f}K)"
                
                # Log the selection (occasionally)
                current_time_ps = self.time_tracker.elapsed_time
                if hasattr(self, '_last_dual_log_time'):
                    if current_time_ps - self._last_dual_log_time > 10.0:  # Log every 10 ps
                        print(f"Dual control | T_kinetic={T_kinetic_K:.1f}K (err={error_kinetic_signed:+.1f}K) | "
                              f"T_lj_coul={T_lj_coulombic_K:.1f}K (err={error_lj_coulombic_signed:+.1f}K) | "
                              f"Using: {selected_method} ({reason})")
                        self._last_dual_log_time = current_time_ps
                else:
                    self._last_dual_log_time = current_time_ps
                
                return selected_temperature
            
            elif self.temperature_method == 'lj_coulombic_bath':
                # Dual temperature control: use maximum error signal from either LJ+Coulombic or bath
                # Error signal = sign * max(|T_lj_coulombic - T_target|, |T_bath - T_target|)
                
                # Get bath temperature from current thermostat
                T_bath_K = self._get_current_thermostat_temperature()
                if T_bath_K is None:
                    print("Warning: Bath temperature not available")
                    return None
                
                # Calculate LJ+Coulombic fictive temperature from energy data (same as lj_coulombic_kinetic)
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if self.empirical_data is None:
                    print(f"Warning: No empirical data available for LJ+Coulombic temperature")
                    return None
                
                if 'lj' not in energy_data or 'coulombic' not in energy_data:
                    print(f"Warning: LJ+Coulombic energy not available")
                    return None
                
                lj_coulombic_energy = energy_data['lj'] + energy_data['coulombic']
                T_lj_coulombic_K = self.empirical_data.calculate_systemic_temperature(lj_coulombic_energy)
                if T_lj_coulombic_K is None:
                    return None
                T_lj_coulombic_K = max(T_lj_coulombic_K, 0.0)
                
                # Calculate signed errors for directional control
                error_bath_signed = T_bath_K - self.target_temperature
                error_lj_coulombic_signed = T_lj_coulombic_K - self.target_temperature
                
                # Use directionally-aware selection: choose temperature requiring larger control action
                error_bath_abs = abs(error_bath_signed)
                error_lj_coulombic_abs = abs(error_lj_coulombic_signed)
                
                # Selection logic: choose temperature with larger absolute error
                # In case of tie (within 0.1 K), prioritize LJ+Coulombic (more fundamental)
                tolerance = 0.1  # K
                if abs(error_bath_abs - error_lj_coulombic_abs) <= tolerance:
                    selected_temp = T_lj_coulombic_K
                    selected_type = "lj_coulombic"
                    reason = f"tie-break (Δerr={error_bath_abs-error_lj_coulombic_abs:.2f}K)"
                elif error_bath_abs > error_lj_coulombic_abs:
                    selected_temp = T_bath_K
                    selected_type = "bath"
                    reason = f"larger error ({error_bath_abs:.1f}K > {error_lj_coulombic_abs:.1f}K)"
                else:
                    selected_temp = T_lj_coulombic_K
                    selected_type = "lj_coulombic" 
                    reason = f"larger error ({error_lj_coulombic_abs:.1f}K > {error_bath_abs:.1f}K)"
                
                # Log the selection (occasionally)
                current_time_ps = self.time_tracker.elapsed_time
                if hasattr(self, '_last_bath_dual_log_time'):
                    if current_time_ps - self._last_bath_dual_log_time > 10.0:  # Log every 10 ps
                        print(f"Bath dual control | T_bath={T_bath_K:.1f}K (err={error_bath_signed:+.1f}K) | "
                              f"T_lj_coul={T_lj_coulombic_K:.1f}K (err={error_lj_coulombic_signed:+.1f}K) | "
                              f"Using: {selected_type} ({reason})")
                        self._last_bath_dual_log_time = current_time_ps
                else:
                    self._last_bath_dual_log_time = current_time_ps
                
                return selected_temp
                
            elif self.temperature_method == 'harmonic':
                if self.empirical_data is None:
                    return None
                
                # Get harmonic energy
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                
                # Convert to temperature using empirical data
                temperature = self.empirical_data.calculate_systemic_temperature(harmonic_energy)
                return max(temperature, 0.0)  # Ensure non-negative temperature
                
            elif self.temperature_method == 'harmonic_equipartition':
                # Calculate harmonic equipartition temperature directly from harmonic energy
                # Using the same method as TemperatureTracker._calculate_harmonic_equipartition_temperature
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                
                if harmonic_energy <= 0:
                    return 0.0
                
                # Get number of particles
                from .utils import PhysicalConstants
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    N_particles = len(snap.particles.mass)
                
                if N_particles == 0:
                    return 0.0
                
                # Direct equipartition calculation: T = 4*E/(N*kB) for 3D harmonic oscillator
                # Based on empirical observations from TemperatureTracker implementation
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (4.0 * harmonic_energy) / (N_particles * kB_hartree)
                return temperature_K
            
            else:
                print(f"Warning: Unknown temperature method '{self.temperature_method}'")
                return None
                
        except Exception as e:
            print(f"Warning: Failed to calculate {self.temperature_method} temperature: {e}")
            return None
    
    def _get_current_thermostat_temperature(self) -> float:
        """Get current bath temperature from thermostats."""
        try:
            from .utils import PhysicalConstants
            
            # Try to get temperature from molecular thermostat first
            if self.molecular_thermostat is not None:
                try:
                    if hasattr(self.molecular_thermostat, 'kT'):
                        kT = self.molecular_thermostat.kT
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.molecular_thermostat, 'thermostat'):
                        nested = self.molecular_thermostat.thermostat
                        if hasattr(nested, 'kT'):
                            kT = nested.kT
                            if hasattr(kT, 'value'):
                                kT = kT.value
                            elif callable(kT):
                                kT = kT()
                            return kT / PhysicalConstants.KB_HARTREE_PER_K
                except:
                    pass
            
            # Try cavity thermostat as backup
            if self.cavity_thermostat is not None:
                try:
                    if hasattr(self.cavity_thermostat, 'kT'):
                        kT = self.cavity_thermostat.kT
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.cavity_thermostat, 'thermostat'):
                        nested = self.cavity_thermostat.thermostat
                        if hasattr(nested, 'kT'):
                            kT = nested.kT
                            if hasattr(kT, 'value'):
                                kT = kT.value
                            elif callable(kT):
                                kT = kT()
                            return kT / PhysicalConstants.KB_HARTREE_PER_K
                except:
                    pass
            
            # Fallback to target temperature
            return self.target_temperature
            
        except Exception as e:
            print(f"Warning: Could not get current thermostat temperature: {e}")
            return self.target_temperature

    def _apply_rate_limit(self, new_output: float, dt_ps: float) -> float:
        """Apply rate limiting to controller output."""
        if self.rate_limit_K_per_ps is None or self.current_bath_temperature is None:
            return new_output
        
        max_change = self.rate_limit_K_per_ps * dt_ps
        change = new_output - self.current_bath_temperature
        
        if abs(change) <= max_change:
            return new_output
        else:
            # Limit the change
            limited_change = max_change if change > 0 else -max_change
            return self.current_bath_temperature + limited_change
    
    def _update_thermostats(self, bath_temperature: float, measured_temperature: float = None, effective_temperature: float = None):
        """Update thermostat temperatures using robust API detection."""
        from .utils import PhysicalConstants
        target_kT = bath_temperature * PhysicalConstants.KB_HARTREE_PER_K
        
        updated_any = False
        
        # Update molecular thermostat
        if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
            try:
                if hasattr(self.molecular_thermostat, 'kT'):
                    # Direct kT attribute (Langevin thermostats)
                    self.molecular_thermostat.kT = target_kT
                    updated_any = True
                elif hasattr(self.molecular_thermostat, 'thermostat'):
                    # Nested thermostat (ConstantVolume with Bussi/MTTK)
                    nested_thermo = self.molecular_thermostat.thermostat
                    if hasattr(nested_thermo, 'kT'):
                        nested_thermo.kT = target_kT
                        updated_any = True
                    else:
                        print(f"Warning: Cannot update nested molecular thermostat kT")
                else:
                    print(f"Warning: Cannot find kT attribute in molecular thermostat")
            except Exception as e:
                print(f"Warning: Failed to update molecular thermostat: {e}")
        
        # Update cavity thermostat  
        if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
            try:
                if hasattr(self.cavity_thermostat, 'kT'):
                    # Direct kT attribute (Langevin thermostats)
                    self.cavity_thermostat.kT = target_kT
                    updated_any = True
                elif hasattr(self.cavity_thermostat, 'thermostat'):
                    # Nested thermostat (ConstantVolume with Bussi/MTTK)
                    nested_thermo = self.cavity_thermostat.thermostat
                    if hasattr(nested_thermo, 'kT'):
                        nested_thermo.kT = target_kT
                        updated_any = True
                    else:
                        print(f"Warning: Cannot update nested cavity thermostat kT")
                else:
                    print(f"Warning: Cannot find kT attribute in cavity thermostat")
            except Exception as e:
                print(f"Warning: Failed to update cavity thermostat: {e}")
        
        # Check if we should print detailed console output (sync with console output period)
        current_time_ps = self.time_tracker.elapsed_time
        should_print_console = (current_time_ps - self.last_console_output_time) >= self.console_output_period_ps
        
        if updated_any and should_print_console:
            # Print detailed information
            if measured_temperature is not None and effective_temperature is not None:
                print(f"GD Controller | Target: {self.target_temperature:.1f}K | Measured: {measured_temperature:.1f}K | "
                      f"Effective: {effective_temperature:.1f}K | Bath→{bath_temperature:.1f}K")
            else:
                print(f"GD Controller | Target: {self.target_temperature:.1f}K | Bath→{bath_temperature:.1f}K")
            self.last_console_output_time = current_time_ps

    def act(self, timestep):
        """Execute gradient descent control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should be active
        should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                          (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
        
        if should_be_active and not self.is_active:
            self.is_active = True
            # Initialize current_bath_temperature from thermostats when first activating
            if self.current_bath_temperature is None:
                self.current_bath_temperature = self._get_current_thermostat_temperature()
                print(f"Gradient descent feedback turned ON at t = {current_time_ps:.2f} ps")
                print(f"Initial bath temperature: {self.current_bath_temperature:.2f} K")
            else:
                print(f"Gradient descent feedback turned ON at t = {current_time_ps:.2f} ps")
        elif not should_be_active and self.is_active:
            self.is_active = False
            print(f"Gradient descent feedback turned OFF at t = {current_time_ps:.2f} ps")
        
        # Skip if not active or not initialized
        if not self.is_active or self.current_bath_temperature is None:
            return
        
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        # Calculate dt_ps BEFORE updating last_update_time
        dt_ps = current_time_ps - self.last_update_time if self.last_update_time > 0 else self.update_interval_ps
        self.last_update_time = current_time_ps
        
        # Calculate system temperature
        measured_temperature = self._calculate_system_temperature(current_time_ps)
        
        if measured_temperature is None:
            return
        
        # Calculate effective system temperature using bath-system coupling model
        # Model: T_effective = (T_measured + T_bath) / 2
        effective_temperature = (measured_temperature + self.current_bath_temperature) / 2.0
        
        # Calculate error and gradient w.r.t. T_bath
        error = effective_temperature - self.target_temperature  # e = T_effective - T_target
        gradient = 0.5 * error  # ∂J/∂T_bath = 0.5 * (T_effective - T_target)
        
        # Apply gradient descent update: T_bath(t+1) = T_bath(t) - α * gradient
        raw_output = self.current_bath_temperature - self.alpha * gradient
        
        # Apply saturation limits
        if self.T_max is not None:
            saturated_output = max(self.T_min, min(self.T_max, raw_output))
        else:
            saturated_output = max(self.T_min, raw_output)
        
        # Apply rate limiting
        rate_limited_output = self._apply_rate_limit(saturated_output, dt_ps)
        
        # Update current bath temperature and thermostats
        self.current_bath_temperature = rate_limited_output
        self._update_thermostats(self.current_bath_temperature, measured_temperature, effective_temperature)
        
        # Log data
        try:
            with open(self.output_file, 'a', encoding='utf-8') as f:
                f.write(f"{current_time_ps:.6f},{measured_temperature:.6f},{effective_temperature:.6f},"
                       f"{self.target_temperature:.6f},{error:.6f},{gradient:.6f},{raw_output:.6f},"
                       f"{saturated_output:.6f},{rate_limited_output:.6f},{int(self.is_active)}\n")
        except Exception as e:
            print(f"Warning: Failed to write gradient descent output: {e}")


class DualIndependentTemperatureFeedback(hoomd.custom.Action):
    """
    Dual independent temperature feedback controller for cavity MD simulations.
    
    This class implements two independent gradient descent loops - one for the cavity
    bath temperature and one for the molecular bath temperature. Each bath responds
    to a different temperature signal, allowing for precise control of different
    aspects of the system dynamics.
    
    The gradient descent update rule for each bath is:
    T_bath(t+1) = T_bath(t) - α * ∇J/∇T_bath
    T_bath(t+1) = T_bath(t) - α * 0.5 * (T_eff(t) - T_target)
    
    where T_eff(t) = (T_measured(t) + T_bath(t)) / 2 for each signal.
    
    Parameters
    ----------
    cavity_method : str
        Temperature calculation method for cavity bath control
        ('kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic')
    molecular_method : str
        Temperature calculation method for molecular bath control
        ('kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic')
    cavity_time_constant_ps : float
        Time constant for cavity bath gradient descent in picoseconds
    molecular_time_constant_ps : float
        Time constant for molecular bath gradient descent in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    energy_tracker : EnergyTracker
        Energy tracker for temperature calculations
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat
        Molecular thermostat to control
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat
        Cavity thermostat to control
    cavity_target_temperature : float, optional
        Target temperature for cavity bath in Kelvin. Default: 100.0
    molecular_target_temperature : float, optional
        Target temperature for molecular bath in Kelvin. Default: 100.0
    turn_on_time_ps : float, optional
        Time to start control in picoseconds. Default: 0.0
    turn_off_time_ps : float, optional
        Time to stop control in picoseconds. Default: None (never turn off)
    update_interval_ps : float, optional
        Interval between control updates in picoseconds. Default: 0.1
    cavity_T_min : float, optional
        Minimum allowed cavity bath temperature in Kelvin. Default: 0.0
    cavity_T_max : float, optional
        Maximum allowed cavity bath temperature in Kelvin. Default: None
    molecular_T_min : float, optional
        Minimum allowed molecular bath temperature in Kelvin. Default: 0.0
    molecular_T_max : float, optional
        Maximum allowed molecular bath temperature in Kelvin. Default: None
    output_file : str, optional
        Output file for control data. Default: 'dual_feedback.csv'
    empirical_data_file : str, optional
        Path to empirical data file (required for certain methods)
    console_output_period_ps : float, optional
        Console output period in picoseconds. Default: 1.0
        
    Examples
    --------
    **Basic dual independent control:**
    
    >>> from hoomd.cavitymd.analysis import DualIndependentTemperatureFeedback
    >>> 
    >>> dual_controller = DualIndependentTemperatureFeedback(
    ...     cavity_method='harmonic_equipartition',
    ...     molecular_method='lj_coulombic_kinetic',
    ...     cavity_time_constant_ps=5.0,
    ...     molecular_time_constant_ps=10.0,
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     cavity_target_temperature=80.0,
    ...     molecular_target_temperature=120.0
    ... )
    
    Notes
    -----
    - Each bath has independent control parameters and targets
    - Different temperature signals can be used for different physical insights
    - Maintains compatibility with all existing temperature calculation methods
    """
    
    def __init__(self, 
                 cavity_method: str,
                 molecular_method: str,
                 cavity_time_constant_ps: float,
                 molecular_time_constant_ps: float,
                 time_tracker,
                 energy_tracker,
                 simulation=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 cavity_target_temperature: float = 100.0,
                 molecular_target_temperature: float = 100.0,
                 cavity_dynamic_target: bool = False,
                 molecular_dynamic_target: bool = False,
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 cavity_T_min: float = 0.0,
                 cavity_T_max: Optional[float] = None,
                 molecular_T_min: float = 0.0,
                 molecular_T_max: Optional[float] = None,
                 output_file: str = 'dual_feedback.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0,
                 cavity_integral_time_constant_ps: Optional[float] = None,
                 molecular_integral_time_constant_ps: Optional[float] = None,
                 cavity_heating_gain_factor: float = 1.0,
                 cavity_cooling_gain_factor: float = 1.0,
                 molecular_heating_gain_factor: float = 1.0,
                 molecular_cooling_gain_factor: float = 1.0,
                 cavity_integral_heating_time_constant_ps: Optional[float] = None,
                 cavity_integral_cooling_time_constant_ps: Optional[float] = None,
                 molecular_integral_heating_time_constant_ps: Optional[float] = None,
                 molecular_integral_cooling_time_constant_ps: Optional[float] = None):
        
        # Store configuration
        self.cavity_method = cavity_method
        self.molecular_method = molecular_method
        self.cavity_time_constant_ps = cavity_time_constant_ps
        self.molecular_time_constant_ps = molecular_time_constant_ps
        self.cavity_integral_time_constant_ps = cavity_integral_time_constant_ps
        self.molecular_integral_time_constant_ps = molecular_integral_time_constant_ps
        
        # Asymmetric gain factors
        self.cavity_heating_gain_factor = cavity_heating_gain_factor
        self.cavity_cooling_gain_factor = cavity_cooling_gain_factor
        self.molecular_heating_gain_factor = molecular_heating_gain_factor
        self.molecular_cooling_gain_factor = molecular_cooling_gain_factor
        
        # Asymmetric integral time constants
        self.cavity_integral_heating_time_constant_ps = cavity_integral_heating_time_constant_ps
        self.cavity_integral_cooling_time_constant_ps = cavity_integral_cooling_time_constant_ps
        self.molecular_integral_heating_time_constant_ps = molecular_integral_heating_time_constant_ps
        self.molecular_integral_cooling_time_constant_ps = molecular_integral_cooling_time_constant_ps
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.cavity_target_temperature = cavity_target_temperature
        self.molecular_target_temperature = molecular_target_temperature
        self.cavity_dynamic_target = cavity_dynamic_target
        self.molecular_dynamic_target = molecular_dynamic_target
        self.turn_on_time_ps = turn_on_time_ps
        self.turn_off_time_ps = turn_off_time_ps
        self.update_interval_ps = update_interval_ps
        self.cavity_T_min = cavity_T_min
        self.cavity_T_max = cavity_T_max
        self.molecular_T_min = molecular_T_min
        self.molecular_T_max = molecular_T_max
        self.output_file = output_file
        self.console_output_period_ps = console_output_period_ps
        
        # Initialize state variables
        self.last_update_time = 0.0
        self.last_console_output_time = 0.0
        self.is_active = False
        self.current_cavity_temperature = None
        self.current_molecular_temperature = None
        self.dynamic_targets_set = False
        
        # Calculate learning rates (alpha = dt / time_constant)
        # We'll update these with actual dt during first timestep
        self.cavity_alpha = None
        self.molecular_alpha = None
        
        # Integral control parameters
        self.cavity_alpha_i = None  # Integral gain for cavity
        self.molecular_alpha_i = None  # Integral gain for molecular
        
        # Error accumulators for integral terms
        self.cavity_error_integral = 0.0
        self.molecular_error_integral = 0.0
        
        # Set up empirical data if needed
        self.cavity_empirical_data = None
        self.molecular_empirical_data = None
        
        if empirical_data_file:
            # Load empirical data for cavity method if needed
            if self.cavity_method in ['lj_coulombic', 'harmonic', 'lj_coulombic_kinetic', 'lj_coulombic_bath']:
                try:
                    if self.cavity_method == 'harmonic':
                        self.cavity_empirical_data = EmpiricalTemperatureData(
                            empirical_data_file, energy_component='harmonic', create_plots=True
                        )
                    else:
                        self.cavity_empirical_data = EmpiricalTemperatureData(
                            empirical_data_file, energy_component='lj_coulombic', create_plots=True
                        )
                    print(f"Loaded cavity empirical data for {self.cavity_method} method")
                except Exception as e:
                    print(f"Warning: Failed to load cavity empirical data: {e}")
            
            # Load empirical data for molecular method if needed
            if self.molecular_method in ['lj_coulombic', 'harmonic', 'lj_coulombic_kinetic', 'lj_coulombic_bath']:
                try:
                    if self.molecular_method == 'harmonic':
                        self.molecular_empirical_data = EmpiricalTemperatureData(
                            empirical_data_file, energy_component='harmonic', create_plots=True
                        )
                    else:
                        self.molecular_empirical_data = EmpiricalTemperatureData(
                            empirical_data_file, energy_component='lj_coulombic', create_plots=True
                        )
                    print(f"Loaded molecular empirical data for {self.molecular_method} method")
                except Exception as e:
                    print(f"Warning: Failed to load molecular empirical data: {e}")
        
        # Validate configuration
        self._validate_configuration()
        
        # Initialize output file
        self._initialize_output_file()
        
        # Print configuration
        print(f"DualIndependent Controller initialized:")
        print(f"   Cavity: {self.cavity_method} → target {self.cavity_target_temperature:.1f}K (τ={self.cavity_time_constant_ps:.1f}ps)")
        print(f"   Molecular: {self.molecular_method} → target {self.molecular_target_temperature:.1f}K (τ={self.molecular_time_constant_ps:.1f}ps)")
        print(f"   Update interval: {self.update_interval_ps:.3f} ps")
        print(f"   Turn-on time: {self.turn_on_time_ps:.1f} ps")
        if self.turn_off_time_ps is not None:
            print(f"   Turn-off time: {self.turn_off_time_ps:.1f} ps")
    
    def _validate_configuration(self):
        """Validate controller configuration."""
        valid_methods = ['kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_kinetic', 'lj_coulombic_bath']
        
        if self.cavity_method not in valid_methods:
            raise ValueError(f"Invalid cavity_method: {self.cavity_method}")
        
        if self.molecular_method not in valid_methods:
            raise ValueError(f"Invalid molecular_method: {self.molecular_method}")
        
        if self.molecular_thermostat is None:
            raise ValueError("molecular_thermostat cannot be None")
        
        if self.cavity_thermostat is None:
            raise ValueError("cavity_thermostat cannot be None")
        
        if self.energy_tracker is None:
            raise ValueError("energy_tracker is required for temperature calculations")
        
        if self.cavity_time_constant_ps <= 0:
            raise ValueError("cavity_time_constant_ps must be positive")
        
        if self.molecular_time_constant_ps <= 0:
            raise ValueError("molecular_time_constant_ps must be positive")
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("time_ps,cavity_temp,molecular_temp,cavity_effective,molecular_effective,")
                f.write("cavity_target,molecular_target,cavity_error,molecular_error,")
                f.write("cavity_gradient,molecular_gradient,cavity_integral,molecular_integral,")
                f.write("cavity_integral_term,molecular_integral_term,")
                f.write("cavity_p_gain_factor,molecular_p_gain_factor,")
                f.write("cavity_tau_i,molecular_tau_i,cavity_output,molecular_output,active\n")
        except Exception as e:
            print(f"Warning: Failed to initialize dual feedback output file: {e}")
    
    def _calculate_system_temperature(self, method: str, empirical_data, current_time_ps: float) -> Optional[float]:
        """Calculate system temperature using specified method."""
        try:
            if method == 'kinetic':
                # Calculate kinetic temperature from particle velocities
                import numpy as np
                from .utils import PhysicalConstants
                
                # Get velocities from simulation state
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    velocities = np.array(snap.particles.velocity)
                    masses = np.array(snap.particles.mass)
                    types = np.array(snap.particles.typeid)
                
                # Filter out cavity particles (assuming they have type 1, molecular particles type 0)
                molecular_mask = (types == 0)
                molecular_velocities = velocities[molecular_mask]
                molecular_masses = masses[molecular_mask]
                
                if len(molecular_velocities) == 0:
                    return 0.0
                
                # Calculate kinetic energy per particle
                ke_per_particle = 0.5 * np.sum(molecular_masses[:, np.newaxis] * molecular_velocities**2, axis=1)
                mean_ke_per_particle = np.mean(ke_per_particle)
                
                # Convert to temperature: <KE>/particle = (3/2) * kB * T
                # Therefore: T = (2/3) * <KE>/particle / kB
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (2.0/3.0) * mean_ke_per_particle / kB_hartree
                
                return temperature_K
                
            elif method == 'lj_coulombic':
                if empirical_data is None:
                    return None
                
                # Get LJ+Coulombic energy
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                if 'lj' not in energy_dict or 'coulombic' not in energy_dict:
                    return None
                
                total_energy = energy_dict['lj'] + energy_dict['coulombic']
                
                # Convert to temperature using empirical data
                temperature = empirical_data.calculate_systemic_temperature(total_energy)
                return max(temperature, 0.0)
                
            elif method == 'lj_coulombic_kinetic':
                # Dual signal: use maximum error from either LJ+Coulombic or kinetic
                
                # Calculate kinetic temperature
                kinetic_temp = self._calculate_system_temperature('kinetic', None, current_time_ps)
                if kinetic_temp is None:
                    return None
                
                # Calculate LJ+Coulombic temperature
                lj_coul_temp = self._calculate_system_temperature('lj_coulombic', empirical_data, current_time_ps)
                if lj_coul_temp is None:
                    return kinetic_temp
                
                # Return the temperature that gives the maximum error magnitude
                # This requires knowing the target, so we'll use a simplified approach
                # and return the average of the two signals
                return (kinetic_temp + lj_coul_temp) / 2.0
                
            elif method == 'harmonic':
                if empirical_data is None:
                    return None
                
                # Get harmonic energy
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                
                # Convert to temperature using empirical data
                temperature = empirical_data.calculate_systemic_temperature(harmonic_energy)
                return max(temperature, 0.0)
                
            elif method == 'harmonic_equipartition':
                # Calculate harmonic equipartition temperature directly
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                
                if harmonic_energy <= 0:
                    return 0.0
                
                # Get number of particles
                from .utils import PhysicalConstants
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    N_particles = len(snap.particles.mass)
                
                if N_particles == 0:
                    return 0.0
                
                # Direct equipartition calculation
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (4.0 * harmonic_energy) / (N_particles * kB_hartree)
                return temperature_K
                
            elif method == 'lj_coulombic_bath':
                # Calculate temperature from molecular bath thermostat using LJ+Coulombic coupling
                # This method combines the bath temperature with LJ+Coulombic energy weighting
                
                # Get current molecular bath temperature
                molecular_bath_temp = self._get_thermostat_temperature(self.molecular_thermostat)
                if molecular_bath_temp is None:
                    return None
                
                # Get LJ+Coulombic energy for weighting
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                if 'lj' not in energy_dict or 'coulombic' not in energy_dict:
                    # Fallback to bath temperature if energy components unavailable
                    return molecular_bath_temp
                
                lj_coulombic_energy = energy_dict['lj'] + energy_dict['coulombic']
                
                # Use empirical data if available for energy-based correction
                if empirical_data is not None:
                    try:
                        energy_temp = empirical_data.calculate_systemic_temperature(lj_coulombic_energy)
                        # Weighted average: 70% bath temperature, 30% energy-derived temperature
                        # This provides a stable signal that incorporates both thermostat control
                        # and instantaneous energy fluctuations
                        return 0.7 * molecular_bath_temp + 0.3 * max(energy_temp, 0.0)
                    except Exception:
                        # Fallback to bath temperature if empirical calculation fails
                        return molecular_bath_temp
                else:
                    # Without empirical data, return the bath temperature directly
                    # This ensures the method always provides a meaningful signal
                    return molecular_bath_temp
            
            else:
                print(f"Warning: Unknown temperature method '{method}'")
                return None
                
        except Exception as e:
            print(f"Warning: Could not calculate temperature using method '{method}': {e}")
            return None
    
    def _get_thermostat_temperature(self, thermostat) -> Optional[float]:
        """Get current temperature from a thermostat."""
        try:
            from .utils import PhysicalConstants
            
            if hasattr(thermostat, 'kT'):
                # Direct kT attribute (Langevin thermostats)
                kT_value = thermostat.kT
                # Handle HOOMD variants (e.g., Constant) 
                if hasattr(kT_value, '__call__'):
                    # It's a variant, call it with timestep 0 to get current value
                    kT_hartree = kT_value(0)
                else:
                    # It's a plain number
                    kT_hartree = kT_value
                return kT_hartree / PhysicalConstants.KB_HARTREE_PER_K
            elif hasattr(thermostat, 'thermostat'):
                # Nested thermostat (ConstantVolume with Bussi/MTTK)
                nested_thermo = thermostat.thermostat
                if hasattr(nested_thermo, 'kT'):
                    kT_value = nested_thermo.kT
                    # Handle HOOMD variants
                    if hasattr(kT_value, '__call__'):
                        kT_hartree = kT_value(0)
                    else:
                        kT_hartree = kT_value
                    return kT_hartree / PhysicalConstants.KB_HARTREE_PER_K
            
            return None
            
        except Exception as e:
            print(f"Warning: Failed to get thermostat temperature: {e}")
            return None
    
    def _update_thermostat(self, thermostat, temperature_K: float) -> bool:
        """Update a thermostat temperature."""
        try:
            import hoomd
            from .utils import PhysicalConstants
            target_kT = temperature_K * PhysicalConstants.KB_HARTREE_PER_K
            
            if hasattr(thermostat, 'kT'):
                # Direct kT attribute (Langevin thermostats)
                # Always use Constant variant for proper HOOMD compatibility
                thermostat.kT = hoomd.variant.Constant(target_kT)
                return True
            elif hasattr(thermostat, 'thermostat'):
                # Nested thermostat (ConstantVolume with Bussi/MTTK)
                nested_thermo = thermostat.thermostat
                if hasattr(nested_thermo, 'kT'):
                    nested_thermo.kT = hoomd.variant.Constant(target_kT)
                    return True
            
            return False
            
        except Exception as e:
            print(f"Warning: Failed to update thermostat: {e}")
            return False
    
    def _set_dynamic_targets(self, current_time_ps: float):
        """Set dynamic targets based on current bath temperatures when controller is activated."""
        try:
            # Get current bath temperatures from thermostats (not signal temperatures)
            if self.cavity_dynamic_target:
                cavity_bath_temp = self._get_thermostat_temperature(self.cavity_thermostat)
                if cavity_bath_temp is not None:
                    old_target = self.cavity_target_temperature
                    self.cavity_target_temperature = cavity_bath_temp
                    print(f"Dynamic cavity target set: {old_target:.1f} K → {cavity_bath_temp:.1f} K (bath temperature)")
                else:
                    print(f"Warning: Could not get cavity bath temperature for dynamic target, keeping {self.cavity_target_temperature:.1f} K")
            
            if self.molecular_dynamic_target:
                molecular_bath_temp = self._get_thermostat_temperature(self.molecular_thermostat)
                if molecular_bath_temp is not None:
                    old_target = self.molecular_target_temperature
                    self.molecular_target_temperature = molecular_bath_temp
                    print(f"Dynamic molecular target set: {old_target:.1f} K → {molecular_bath_temp:.1f} K (bath temperature)")
                else:
                    print(f"Warning: Could not get molecular bath temperature for dynamic target, keeping {self.molecular_target_temperature:.1f} K")
            
            self.dynamic_targets_set = True
            
        except Exception as e:
            print(f"Warning: Error setting dynamic targets: {e}")
            self.dynamic_targets_set = True  # Don't try again
        
    def act(self, timestep):
        """Execute dual independent control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should be active
        should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                          (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
        
        if should_be_active and not self.is_active:
            self.is_active = True
            # Initialize current temperatures from thermostats when first activating
            self.current_cavity_temperature = self._get_thermostat_temperature(self.cavity_thermostat)
            self.current_molecular_temperature = self._get_thermostat_temperature(self.molecular_thermostat)
            
            # Initialize learning rates based on actual timestep
            if self.cavity_alpha is None:
                # Estimate timestep from update interval (rough approximation)
                estimated_dt_ps = self.update_interval_ps
                
                # Proportional gains (existing)
                self.cavity_alpha = estimated_dt_ps / self.cavity_time_constant_ps
                self.molecular_alpha = estimated_dt_ps / self.molecular_time_constant_ps
                
                # Integral gains (new PI controller)
                if self.cavity_integral_time_constant_ps is not None:
                    self.cavity_alpha_i = estimated_dt_ps / self.cavity_integral_time_constant_ps
                    print(f"Cavity PI controller: τₚ={self.cavity_time_constant_ps:.1f}ps, τᵢ={self.cavity_integral_time_constant_ps:.1f}ps")
        else:
                self.cavity_alpha_i = 0.0
                print(f"Cavity P controller: τₚ={self.cavity_time_constant_ps:.1f}ps (no integral)")
                
                if self.molecular_integral_time_constant_ps is not None:
                    self.molecular_alpha_i = estimated_dt_ps / self.molecular_integral_time_constant_ps
                    print(f"Molecular PI controller: τₚ={self.molecular_time_constant_ps:.1f}ps, τᵢ={self.molecular_integral_time_constant_ps:.1f}ps")
                else:
                    self.molecular_alpha_i = 0.0
                    print(f"Molecular P controller: τₚ={self.molecular_time_constant_ps:.1f}ps (no integral)")
            
        print(f"Dual independent feedback turned ON at t = {current_time_ps:.2f} ps")
        cavity_temp_str = f"{self.current_cavity_temperature:.2f}" if self.current_cavity_temperature is not None else "None"
        molecular_temp_str = f"{self.current_molecular_temperature:.2f}" if self.current_molecular_temperature is not None else "None"
        print(f"Initial cavity temperature: {cavity_temp_str} K")
        print(f"Initial molecular temperature: {molecular_temp_str} K")
            
        # Set dynamic targets if enabled
        if not self.dynamic_targets_set:
            self._set_dynamic_targets(current_time_ps)
            
        elif not should_be_active and self.is_active:
            self.is_active = False
            print(f"Dual independent feedback turned OFF at t = {current_time_ps:.2f} ps")
        
        # Skip if not active or not initialized
        if not self.is_active or self.current_cavity_temperature is None or self.current_molecular_temperature is None:
            return
        
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        # Calculate dt_ps BEFORE updating last_update_time
        dt_ps = current_time_ps - self.last_update_time if self.last_update_time > 0 else self.update_interval_ps
        self.last_update_time = current_time_ps
        
        # Calculate system temperatures for each method
        cavity_measured_temp = self._calculate_system_temperature(
            self.cavity_method, self.cavity_empirical_data, current_time_ps)
        molecular_measured_temp = self._calculate_system_temperature(
            self.molecular_method, self.molecular_empirical_data, current_time_ps)
        
        if cavity_measured_temp is None or molecular_measured_temp is None:
            return
        
        # Calculate effective temperatures and errors for each bath
        cavity_effective_temp = (cavity_measured_temp + self.current_cavity_temperature) / 2.0
        molecular_effective_temp = (molecular_measured_temp + self.current_molecular_temperature) / 2.0
        
        cavity_error = cavity_effective_temp - self.cavity_target_temperature
        molecular_error = molecular_effective_temp - self.molecular_target_temperature
        
        # Update integral error accumulators
        self.cavity_error_integral += cavity_error * dt_ps
        self.molecular_error_integral += molecular_error * dt_ps
        
        # Apply asymmetric gains based on error direction
        # Error < 0: System too cold, need to HEAT (raise bath temperature) 
        # Error > 0: System too hot, need to COOL (lower bath temperature)
        
        # Cavity asymmetric gains
        if cavity_error < 0:  # Need to heat system
            cavity_p_gain_factor = self.cavity_heating_gain_factor
            cavity_tau_i = (self.cavity_integral_heating_time_constant_ps 
                           if self.cavity_integral_heating_time_constant_ps is not None 
                           else self.cavity_integral_time_constant_ps)
        else:  # Need to cool system
            cavity_p_gain_factor = self.cavity_cooling_gain_factor
            cavity_tau_i = (self.cavity_integral_cooling_time_constant_ps 
                           if self.cavity_integral_cooling_time_constant_ps is not None 
                           else self.cavity_integral_time_constant_ps)
        
        # Molecular asymmetric gains
        if molecular_error < 0:  # Need to heat system
            molecular_p_gain_factor = self.molecular_heating_gain_factor
            molecular_tau_i = (self.molecular_integral_heating_time_constant_ps 
                              if self.molecular_integral_heating_time_constant_ps is not None 
                              else self.molecular_integral_time_constant_ps)
        else:  # Need to cool system
            molecular_p_gain_factor = self.molecular_cooling_gain_factor
            molecular_tau_i = (self.molecular_integral_cooling_time_constant_ps 
                              if self.molecular_integral_cooling_time_constant_ps is not None 
                              else self.molecular_integral_time_constant_ps)
        
        # Calculate asymmetric gradients (proportional terms)
        cavity_gradient = 0.5 * cavity_error * cavity_p_gain_factor
        molecular_gradient = 0.5 * molecular_error * molecular_p_gain_factor
        
        # Calculate asymmetric integral terms
        if cavity_tau_i is not None:
            cavity_alpha_i_effective = dt_ps / cavity_tau_i
            cavity_integral_term = cavity_alpha_i_effective * self.cavity_error_integral
        else:
            cavity_integral_term = 0.0
            
        if molecular_tau_i is not None:
            molecular_alpha_i_effective = dt_ps / molecular_tau_i
            molecular_integral_term = molecular_alpha_i_effective * self.molecular_error_integral
        else:
            molecular_integral_term = 0.0
        
        # Apply asymmetric PI control updates
        cavity_raw_output = (self.current_cavity_temperature - 
                            self.cavity_alpha * cavity_gradient - 
                            cavity_integral_term)
        molecular_raw_output = (self.current_molecular_temperature - 
                               self.molecular_alpha * molecular_gradient - 
                               molecular_integral_term)
        
        # Apply saturation limits
        if self.cavity_T_max is not None:
            cavity_output = max(self.cavity_T_min, min(self.cavity_T_max, cavity_raw_output))
        else:
            cavity_output = max(self.cavity_T_min, cavity_raw_output)
        
        if self.molecular_T_max is not None:
            molecular_output = max(self.molecular_T_min, min(self.molecular_T_max, molecular_raw_output))
        else:
            molecular_output = max(self.molecular_T_min, molecular_raw_output)
        
        # Update current temperatures and thermostats
        self.current_cavity_temperature = cavity_output
        self.current_molecular_temperature = molecular_output
        
        cavity_updated = self._update_thermostat(self.cavity_thermostat, cavity_output)
        molecular_updated = self._update_thermostat(self.molecular_thermostat, molecular_output)
        
        # Console output (periodic)
        should_print_console = (current_time_ps - self.last_console_output_time) >= self.console_output_period_ps
        
        if (cavity_updated or molecular_updated) and should_print_console:
            print(f"Dual Controller | Cavity: {self.cavity_method}→{cavity_output:.1f}K | "
                  f"Molecular: {self.molecular_method}→{molecular_output:.1f}K")
            self.last_console_output_time = current_time_ps
        
        # Log data
        try:
            with open(self.output_file, 'a', encoding='utf-8') as f:
                f.write(f"{current_time_ps:.6f},{cavity_measured_temp:.6f},{molecular_measured_temp:.6f},")
                f.write(f"{cavity_effective_temp:.6f},{molecular_effective_temp:.6f},")
                f.write(f"{self.cavity_target_temperature:.6f},{self.molecular_target_temperature:.6f},")
                f.write(f"{cavity_error:.6f},{molecular_error:.6f},")
                f.write(f"{cavity_gradient:.6f},{molecular_gradient:.6f},")
                f.write(f"{self.cavity_error_integral:.6f},{self.molecular_error_integral:.6f},")
                f.write(f"{cavity_integral_term:.6f},{molecular_integral_term:.6f},")
                f.write(f"{cavity_p_gain_factor:.6f},{molecular_p_gain_factor:.6f},")
                f.write(f"{cavity_tau_i if cavity_tau_i is not None else -1:.6f},{molecular_tau_i if molecular_tau_i is not None else -1:.6f},")
                f.write(f"{cavity_output:.6f},{molecular_output:.6f},{int(self.is_active)}\n")
        except Exception as e:
            print(f"Warning: Failed to write dual feedback output: {e}")


class DipoleMomentFDRTracker(hoomd.custom.Action):
    """
    Fluctuation-Dissipation Relation (FDR) tracker for total dipole moment.
    
    This class implements FDR analysis for the total dipole moment by:
    1. Computing equilibrium dipole moment autocorrelations C_μ(t)
    2. Measuring linear response to applied electric fields χ_μ(t)
    3. Validating the FDR relation: χ_μ(ω) = ∫ C_μ(t) e^(iωt) dt / (k_B T)
    
    **Physics Background:**
    
    The dipole moment FDR relates equilibrium fluctuations to linear response:
    
    - **Fluctuation**: C_μ(t) = ⟨δμ⃗(0) · δμ⃗(t)⟩ 
    - **Response**: χ_μ(t,t_w) = (⟨μ⃗^(+)(t)⟩ - ⟨μ⃗^(-)(t)⟩) / (2E₀)
    - **FDR**: χ_μ(ω) = ∫₀^∞ C_μ(t) e^(iωt) dt / (k_B T)
    
    where μ⃗ = Σᵢ qᵢr⃗ᵢ is the total dipole moment.
    
    **Experimental Protocol:**
    
    1. **Equilibrium simulation**: Run unperturbed system, calculate C_μ(t)
    2. **Response measurement**: Fork simulation, apply ±E₀ fields, measure χ_μ(t)
    3. **FDR validation**: Compare χ_μ(ω) with Fourier transform of C_μ(t)
    
    **Connection to Cavity Dynamics:**
    
    - Dipole moment couples to cavity field through light-matter interaction
    - FDR validates linear response assumptions in cavity coupling
    - Provides fundamental check of equilibrium statistical mechanics
    
    Parameters
    ----------
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    output_file : str, optional
        Output file for autocorrelation data. Default: 'dipole_fdr.csv'
    max_correlation_time_ps : float, optional
        Maximum correlation time to track in ps. Default: 100.0
    correlation_output_interval_ps : float, optional
        Interval between correlation outputs in ps. Default: 0.1
    exclude_cavity : bool, optional
        Whether to exclude cavity particles from dipole calculation. Default: True
    field_direction : array_like, optional
        Direction for electric field perturbations (3D vector). Default: [0,0,1]
    enable_response_measurement : bool, optional
        Whether to enable response measurements (requires fork-and-clone). Default: False
    output_period_ps : float, optional
        Output period in ps. Default: 0.1
    Attributes
    ----------
    dipole_history : List[np.ndarray]
        Time series of total dipole moment vectors
    time_history : List[float]
        Corresponding time stamps
    autocorrelation : np.ndarray
        Dipole moment autocorrelation function C_μ(t)
    correlation_times : np.ndarray
        Time grid for autocorrelation
    susceptibility : Optional[np.ndarray]
        Measured susceptibility χ_μ(t) (if response measurement enabled)
        
    Examples
    --------
    **Basic dipole autocorrelation tracking:**
    
    >>> from hoomd.cavitymd.analysis import DipoleMomentFDRTracker
    >>> 
    >>> tracker = DipoleMomentFDRTracker(
    ...     time_tracker=time_tracker,
    ...     max_correlation_time_ps=50.0,
    ...     correlation_output_interval_ps=0.1
    ... )
    >>> 
    >>> # Add to simulation
    >>> updater = hoomd.update.CustomUpdater(
    ...     action=tracker,
    ...     trigger=hoomd.trigger.Periodic(100)  # Every 100 steps
    ... )
    >>> sim.operations.updaters.append(updater)
    
    **Complete FDR measurement with response:**
    
    >>> tracker = DipoleMomentFDRTracker(
    ...     time_tracker=time_tracker,
    ...     field_direction=[0, 0, 1],
    ...     enable_response_measurement=True
    ... )
    
    Notes
    -----
    - Calculates total dipole moment μ⃗ = Σᵢ qᵢr⃗ᵢ for all particles
    - Autocorrelation computed using FFT for efficiency
    - Compatible with cavity MD simulations and periodic boundary conditions
    - Output files contain both raw data and processed correlations
    - FDR validation requires comparison with response measurements
    """
    
    def __init__(self,
                 time_tracker,
                 output_file: str = 'dipole_fdr.csv',
                 max_correlation_time_ps: float = 100.0,
                 correlation_output_interval_ps: float = 0.1,
                 exclude_cavity: bool = True,
                 field_direction: Union[List[float], np.ndarray] = [0, 0, 1],
                 enable_response_measurement: bool = False,
                 output_period_ps: float = 0.1):
        
        super().__init__()
        
        self.time_tracker = time_tracker
        self.output_file = output_file
        self.max_correlation_time_ps = float(max_correlation_time_ps)
        self.correlation_output_interval_ps = float(correlation_output_interval_ps)
        self.exclude_cavity = bool(exclude_cavity)
        self.enable_response_measurement = bool(enable_response_measurement)
        self.output_period_ps = float(output_period_ps)
        # Normalize field direction
        self.field_direction = np.array(field_direction, dtype=np.float64)
        field_magnitude = np.linalg.norm(self.field_direction)
        if field_magnitude > 1e-12:
            self.field_direction = self.field_direction / field_magnitude
        else:
            self.field_direction = np.array([0, 0, 1])  # Default z-direction
        
        # Data storage
        self.dipole_history = []  # List of 3D dipole vectors
        self.time_history = []    # Corresponding time stamps
        self.last_output_time = 0.0
        
        # Correlation analysis
        self.autocorrelation = None
        self.correlation_times = None
        self.susceptibility = None
        
        # Response measurement data (for fork-and-clone FDR)
        self.response_data = {
            'plus_clone': [],   # Dipole data from +E₀ clone
            'minus_clone': [],  # Dipole data from -E₀ clone
            'times': []         # Time stamps for response data
        }
        self.output_period_ps = float(output_period_ps)
        # Initialize output file
        self._initialize_output_file()
        
        print(f"DipoleMomentFDRTracker initialized:")
        print(f"   Output file: {self.output_file}")
        print(f"   Max correlation time: {self.max_correlation_time_ps:.1f} ps")
        print(f"   Field direction: {self.field_direction}")
        print(f"   Exclude cavity: {self.exclude_cavity}")
        print(f"   Output period: {self.output_period_ps:.3f} ps")
        if self.enable_response_measurement:
            print(f"   Response measurement: Enabled (requires fork-and-clone)")
        else:
            print(f"   Response measurement: Disabled (autocorrelation only)")
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("# Dipole Moment FDR Analysis Data\n")
                f.write("# time_ps: Simulation time in picoseconds\n")
                f.write("# dipole_x: x-component of total dipole moment (charge × length units)\n")
                f.write("# dipole_y: y-component of total dipole moment\n")
                f.write("# dipole_z: z-component of total dipole moment\n")
                f.write("# dipole_magnitude: |μ⃗| magnitude of dipole moment\n")
                f.write("# Field direction: [{:.3f}, {:.3f}, {:.3f}]\n".format(*self.field_direction))
                f.write("time_ps,dipole_x,dipole_y,dipole_z,dipole_magnitude\n")
        except Exception as e:
            print(f"Warning: Failed to initialize dipole FDR output file: {e}")
    
    def _calculate_total_dipole_moment(self) -> np.ndarray:
        """Calculate total dipole moment μ⃗ = Σᵢ qᵢr⃗ᵢ.
        
        Returns
        -------
        np.ndarray
            3D dipole moment vector [μₓ, μᵧ, μᵤ]
        """
        # Get current simulation state
        with self._simulation.state.cpu_local_snapshot as snap:
            positions = snap.particles.position
            charges = snap.particles.charge
            N = len(positions)
            
            # Calculate total dipole moment
            dipole = np.zeros(3)
            
            for i in range(N):
                # Skip cavity particles if requested
                if self.exclude_cavity and i >= N - 1:  # Assume cavity is last particle
                    continue
                
                q_i = charges[i]
                if abs(q_i) > 1e-12:  # Only consider charged particles
                    dipole += q_i * positions[i]
            
            return dipole
    
    def _compute_autocorrelation(self) -> Tuple[np.ndarray, np.ndarray]:
        """Compute dipole moment autocorrelation function using FFT.
        
        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            (correlation_times, autocorrelation_values)
        """
        if len(self.dipole_history) < 10:
            return np.array([]), np.array([])
        
        # Convert to numpy arrays
        dipoles = np.array(self.dipole_history)  # Shape: (N_times, 3)
        times = np.array(self.time_history)
        
        # Get time step (assume uniform spacing)
        dt = times[1] - times[0] if len(times) > 1 else self.correlation_output_interval_ps
        
        # Calculate fluctuations: δμ⃗(t) = μ⃗(t) - ⟨μ⃗⟩
        mean_dipole = np.mean(dipoles, axis=0)
        delta_dipoles = dipoles - mean_dipole
        
        # Compute autocorrelation for each component using FFT
        N = len(delta_dipoles)
        correlations = []
        
        for component in range(3):
            # Use FFT-based autocorrelation
            signal = delta_dipoles[:, component]
            
            # Zero-pad for full correlation
            padded_signal = np.zeros(2 * N)
            padded_signal[:N] = signal
            
            # FFT-based autocorrelation
            fft_signal = scipy.fft.fft(padded_signal)
            correlation = scipy.fft.ifft(fft_signal * np.conj(fft_signal)).real
            
            # Take only the first N points (positive lags)
            correlation = correlation[:N]
            
            # Normalize by decreasing number of samples
            normalization = np.arange(N, 0, -1)
            correlation = correlation / normalization
            
            correlations.append(correlation)
        
        # Compute total autocorrelation: C_μ(t) = ⟨δμ⃗(0) · δμ⃗(t)⟩
        total_correlation = np.sum(correlations, axis=0)
        
        # Create time grid
        correlation_times = np.arange(N) * dt
        
        # Limit to requested maximum correlation time
        max_index = int(self.max_correlation_time_ps / dt)
        if max_index < len(correlation_times):
            correlation_times = correlation_times[:max_index]
            total_correlation = total_correlation[:max_index]
        
        return correlation_times, total_correlation
    
    def _save_autocorrelation_data(self, correlation_times: np.ndarray, autocorrelation: np.ndarray):
        """Save autocorrelation data to separate file."""
        correlation_file = self.output_file.replace('.csv', '_autocorrelation.csv')
        
        try:
            with open(correlation_file, 'w', encoding='utf-8') as f:
                f.write("# Dipole Moment Autocorrelation Function\n")
                f.write("# time_ps: Correlation time in picoseconds\n")
                f.write("# C_mu: Autocorrelation ⟨δμ⃗(0) · δμ⃗(t)⟩\n")
                f.write("# C_mu_normalized: Normalized by C_μ(0)\n")
                f.write("time_ps,C_mu,C_mu_normalized\n")
                
                # Normalize by C_μ(0)
                C_0 = autocorrelation[0] if len(autocorrelation) > 0 else 1.0
                normalized_correlation = autocorrelation / C_0 if C_0 != 0 else autocorrelation
                
                for t, C, C_norm in zip(correlation_times, autocorrelation, normalized_correlation):
                    f.write(f"{t:.6f},{C:.6e},{C_norm:.6f}\n")
                    
            print(f"Autocorrelation data saved to {correlation_file}")
            
        except Exception as e:
            print(f"Warning: Failed to save autocorrelation data: {e}")
    
    def add_response_data(self, clone_type: str, dipole_moment: np.ndarray, time_ps: float):
        """Add response data from fork-and-clone measurements.
        
        Parameters
        ----------
        clone_type : str
            Either 'plus' or 'minus' for +E₀ or -E₀ clone
        dipole_moment : np.ndarray
            3D dipole moment vector from the clone
        time_ps : float
            Time stamp for this measurement
        """
        if not self.enable_response_measurement:
            return
        
        if clone_type == 'plus':
            self.response_data['plus_clone'].append(dipole_moment.copy())
        elif clone_type == 'minus':
            self.response_data['minus_clone'].append(dipole_moment.copy())
        else:
            print(f"Warning: Unknown clone type '{clone_type}', expected 'plus' or 'minus'")
            return
        
        # Add time stamp (only once per time point)
        if clone_type == 'plus' and time_ps not in self.response_data['times']:
            self.response_data['times'].append(time_ps)
    
    def compute_susceptibility(self, field_strength: float) -> Tuple[np.ndarray, np.ndarray]:
        """Compute dipole susceptibility from response data.
        
        Parameters
        ----------
        field_strength : float
            Electric field strength E₀ used in response measurement
        
        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            (response_times, susceptibility_values)
        """
        if not self.enable_response_measurement:
            print("Warning: Response measurement not enabled")
            return np.array([]), np.array([])
        
        plus_data = np.array(self.response_data['plus_clone'])
        minus_data = np.array(self.response_data['minus_clone'])
        times = np.array(self.response_data['times'])
        
        if len(plus_data) != len(minus_data) or len(plus_data) == 0:
            print("Warning: Insufficient or mismatched response data")
            return np.array([]), np.array([])
        
        # Calculate susceptibility: χ_μ(t) = (⟨μ⃗^(+)(t)⟩ - ⟨μ⃗^(-)(t)⟩) / (2E₀)
        # Project onto field direction
        plus_projection = np.dot(plus_data, self.field_direction)
        minus_projection = np.dot(minus_data, self.field_direction)
        
        susceptibility = (plus_projection - minus_projection) / (2 * field_strength)
        
        return times, susceptibility
    
    def validate_fdr(self, field_strength: float, temperature_K: float) -> Dict[str, Any]:
        """Validate FDR relation between autocorrelation and susceptibility.
        
        Parameters
        ----------
        field_strength : float
            Electric field strength used in response measurement
        temperature_K : float
            System temperature in Kelvin
        
        Returns
        -------
        Dict[str, Any]
            FDR validation results including correlation coefficient
        """
        # Compute autocorrelation
        if self.autocorrelation is None or len(self.autocorrelation) == 0:
            self.correlation_times, self.autocorrelation = self._compute_autocorrelation()
        
        # Compute susceptibility
        response_times, susceptibility = self.compute_susceptibility(field_strength)
        
        if len(self.autocorrelation) == 0 or len(susceptibility) == 0:
            return {'success': False, 'error': 'Insufficient data for FDR validation'}
        
        # Fourier transform autocorrelation to get predicted susceptibility
        from .utils import PhysicalConstants
        kB_T = temperature_K * PhysicalConstants.KB_HARTREE_PER_K
        
        # Simple integration for FDR (more sophisticated FFT could be used)
        dt = self.correlation_times[1] - self.correlation_times[0] if len(self.correlation_times) > 1 else 0.1
        
        # Predict susceptibility from FDR: χ(ω=0) ≈ ∫ C(t) dt / (k_B T)
        predicted_chi_static = np.trapz(self.autocorrelation, dx=dt) / kB_T
        measured_chi_static = susceptibility[0] if len(susceptibility) > 0 else 0.0
        
        # Calculate correlation coefficient for time-dependent comparison
        # Interpolate to common time grid
        min_time = max(self.correlation_times[0], response_times[0])
        max_time = min(self.correlation_times[-1], response_times[-1])
        
        if max_time <= min_time:
            return {'success': False, 'error': 'No overlapping time range for comparison'}
        
        common_times = np.linspace(min_time, max_time, 100)
        
        # Interpolate both functions
        autocorr_interp = np.interp(common_times, self.correlation_times, self.autocorrelation)
        suscept_interp = np.interp(common_times, response_times, susceptibility)
        
        # Normalize and compare
        autocorr_normalized = autocorr_interp / kB_T
        correlation_coeff = np.corrcoef(autocorr_normalized, suscept_interp)[0, 1]
        
        return {
            'success': True,
            'predicted_chi_static': predicted_chi_static,
            'measured_chi_static': measured_chi_static,
            'static_ratio': measured_chi_static / predicted_chi_static if predicted_chi_static != 0 else np.inf,
            'correlation_coefficient': correlation_coeff,
            'temperature_K': temperature_K,
            'field_strength': field_strength,
            'autocorrelation_time_range': (self.correlation_times[0], self.correlation_times[-1]),
            'response_time_range': (response_times[0], response_times[-1])
        }
    
    def act(self, timestep):
        """Main action: calculate dipole moment and update correlations."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if it's time to output
        if current_time_ps - self.last_output_time < self.correlation_output_interval_ps:
            return
        
        self.last_output_time = current_time_ps
        
        # Calculate current dipole moment
        dipole_moment = self._calculate_total_dipole_moment()
        
        # Store in history
        self.dipole_history.append(dipole_moment.copy())
        self.time_history.append(current_time_ps)
        
        # Limit history size to prevent memory issues
        max_history_points = int(self.max_correlation_time_ps / self.correlation_output_interval_ps * 2)
        if len(self.dipole_history) > max_history_points:
            self.dipole_history.pop(0)
            self.time_history.pop(0)
        
        # Log to file
        try:
            dipole_magnitude = np.linalg.norm(dipole_moment)
            with open(self.output_file, 'a', encoding='utf-8') as f:
                f.write(f"{current_time_ps:.6f},{dipole_moment[0]:.6e},"
                       f"{dipole_moment[1]:.6e},{dipole_moment[2]:.6e},{dipole_magnitude:.6e}\n")
        except Exception as e:
            print(f"Warning: Failed to write dipole moment data: {e}")
        
        # Periodically compute and save autocorrelation (every ~10 ps)
        if len(self.time_history) > 100 and current_time_ps % 10.0 < self.correlation_output_interval_ps:
            self.correlation_times, self.autocorrelation = self._compute_autocorrelation()
            if len(self.autocorrelation) > 0:
                self._save_autocorrelation_data(self.correlation_times, self.autocorrelation)
        self.last_update_time = current_time_ps


class OffsetTemperatureController(hoomd.custom.Action):
    """
    Offset temperature controller for cavity MD simulations.
    
    This controller sets the bath temperature to a fixed offset from a target temperature:
    T_bath = T_target + dT_offset
    
    Unlike gradient descent controllers, this provides direct, immediate temperature control
    without feedback loops or convergence dynamics. The offset can be positive (heating)
    or negative (cooling) relative to the target.
    
    Parameters
    ----------
    temperature_method : str
        Temperature calculation method to determine the target
        ('kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_bath')
    temperature_offset_K : float
        Temperature offset in Kelvin (can be positive or negative)
        T_bath = T_measured + temperature_offset_K
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
        Output file for control data. Default: 'offset_feedback.csv'
    empirical_data_file : str, optional
        Path to empirical data file (required for certain methods)
    console_output_period_ps : float, optional
        Console output period in picoseconds. Default: 1.0
        
    Examples
    --------
    **Basic offset control - keep bath 50K below measured kinetic temperature:**
    
    >>> from hoomd.cavitymd.analysis import OffsetTemperatureController
    >>> 
    >>> offset_controller = OffsetTemperatureController(
    ...     temperature_method='kinetic',
    ...     temperature_offset_K=-50.0,  # 50K below kinetic temperature
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     apply_to='molecular'
    ... )
    
    **Harmonic equipartition offset control:**
    
    >>> offset_controller = OffsetTemperatureController(
    ...     temperature_method='harmonic_equipartition',
    ...     temperature_offset_K=20.0,  # 20K above harmonic temperature
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     apply_to='both'
    ... )
    
    Notes
    -----
    - No feedback dynamics or convergence time - immediate temperature setting
    - Positive offset = bath hotter than measured temperature
    - Negative offset = bath cooler than measured temperature
    - Compatible with all existing temperature calculation methods
    """
    
    def __init__(self, 
                 temperature_method: str,
                 temperature_offset_K: float,
                 time_tracker,
                 energy_tracker,
                 simulation=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 apply_to: str = 'both',
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 T_min: float = 0.0,
                 T_max: Optional[float] = None,
                 output_file: str = 'offset_feedback.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0):
        
        # Store core parameters
        self.temperature_method = temperature_method
        self.temperature_offset_K = float(temperature_offset_K)
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        
        # Control parameters
        self.apply_to = apply_to
        self.turn_on_time_ps = float(turn_on_time_ps)
        self.turn_off_time_ps = float(turn_off_time_ps) if turn_off_time_ps is not None else None
        self.update_interval_ps = float(update_interval_ps)
        self.T_min = float(T_min)
        self.T_max = float(T_max) if T_max is not None else None
        self.output_file = output_file
        self.console_output_period_ps = float(console_output_period_ps)
        
        # Load empirical data if needed
        self.empirical_data = None
        if empirical_data_file is not None and temperature_method in ['lj_coulombic', 'harmonic', 'lj_coulombic_bath']:
            try:
                energy_component = 'lj_coulombic' if 'lj_coulombic' in temperature_method else 'harmonic'
                self.empirical_data = EmpiricalTemperatureData(
                    data_file_path=empirical_data_file,
                    energy_component=energy_component
                )
            except Exception as e:
                print(f"Warning: Could not load empirical data file {empirical_data_file}: {e}")
        
        # Control state
        self.is_active = False
        self.current_bath_temperature = None
        self.last_update_time = 0.0
        self.last_console_output_time = 0.0
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("# Offset Temperature Controller Output\n")
                f.write("# time_ps: simulation time in picoseconds\n")
                f.write("# measured_temp_K: measured system temperature\n")
                f.write("# target_bath_temp_K: target bath temperature (measured + offset)\n")
                f.write("# actual_bath_temp_K: actual bath temperature after limits\n")
                f.write("# offset_K: temperature offset applied\n")
                f.write("# active: whether controller is active (1=active, 0=inactive)\n")
                f.write("time_ps,measured_temp_K,target_bath_temp_K,actual_bath_temp_K,offset_K,active\n")
        except Exception as e:
            print(f"Warning: Failed to initialize offset controller output file: {e}")
    
    def _calculate_system_temperature(self, current_time_ps: float) -> Optional[float]:
        """Calculate system temperature using the specified method."""
        if self.temperature_method == 'kinetic':
            # Calculate kinetic temperature from particle velocities
            import numpy as np
            from .utils import PhysicalConstants
            
            try:
                kinetic_energy = 0.0
                N_molecules = 0
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
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
        
        elif self.temperature_method == 'harmonic_equipartition':
            # Calculate harmonic temperature using equipartition theorem
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                harmonic_energy_hartree = energy_data.get('harmonic', 0.0)
                
                # Use equipartition theorem: E_harmonic = (3/2) * N * kB * T
                # Assume N=500 molecular particles (excluding cavity)
                N_molecules = 500
                from .utils import PhysicalConstants
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                
                if harmonic_energy_hartree > 0:
                    temperature_K = (2.0/3.0) * harmonic_energy_hartree / (N_molecules * kB_hartree)
                    return max(temperature_K, 0.0)
                else:
                    return None
                    
            except Exception as e:
                print(f"Warning: Could not calculate harmonic equipartition temperature: {e}")
                return None
        
        elif self.temperature_method in ['lj_coulombic', 'lj_coulombic_bath']:
            # Calculate LJ+Coulombic fictive temperature using empirical data
            if self.empirical_data is None:
                print(f"Warning: No empirical data available for method '{self.temperature_method}'")
                return None
            
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                lj_energy = energy_data.get('lj', 0.0)
                coulomb_energy = energy_data.get('coulombic', 0.0)
                total_lj_coulomb = lj_energy + coulomb_energy
                
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    instantaneous_energy_hartree=total_lj_coulomb,
                    num_particles=500
                )
                return max(temperature_K, 0.0)
                        
            except Exception as e:
                print(f"Warning: Could not calculate LJ+Coulombic temperature: {e}")
                return None
        
        elif self.temperature_method == 'harmonic':
            # Calculate harmonic fictive temperature using empirical data
            if self.empirical_data is None:
                print(f"Warning: No empirical data available for method '{self.temperature_method}'")
                return None
            
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                harmonic_energy = energy_data.get('harmonic', 0.0)
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    instantaneous_energy_hartree=harmonic_energy,
                    num_particles=500
                )
                return max(temperature_K, 0.0)
                
            except Exception as e:
                print(f"Warning: Could not calculate harmonic temperature: {e}")
                return None
        
        else:
            print(f"Warning: Unknown temperature method '{self.temperature_method}'")
            return None
    
    def _update_thermostats(self, bath_temperature: float):
        """Update thermostat temperatures."""
        try:
            from .utils import PhysicalConstants
            import hoomd    
            
            kT_hartree = bath_temperature * PhysicalConstants.KB_HARTREE_PER_K
            
            if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
                if hasattr(self.molecular_thermostat, 'kT'):
                    self.molecular_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                elif hasattr(self.molecular_thermostat, 'thermostat'):
                    nested_thermo = self.molecular_thermostat.thermostat
                    if hasattr(nested_thermo, 'kT'):
                        nested_thermo.kT = hoomd.variant.Constant(kT_hartree)
            
            if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
                if hasattr(self.cavity_thermostat, 'kT'):
                    self.cavity_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                elif hasattr(self.cavity_thermostat, 'thermostat'):
                    nested_thermo = self.cavity_thermostat.thermostat
                    if hasattr(nested_thermo, 'kT'):
                        nested_thermo.kT = hoomd.variant.Constant(kT_hartree)
                        
        except Exception as e:
            print(f"Warning: Failed to update thermostats: {e}")
    
    def act(self, timestep):
        """Execute offset temperature control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should be active
        should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                          (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
        
        if should_be_active and not self.is_active:
            self.is_active = True
            print(f"Offset temperature controller turned ON at t = {current_time_ps:.2f} ps")
            print(f"Temperature offset: {self.temperature_offset_K:+.1f} K")
        elif not should_be_active and self.is_active:
            self.is_active = False
            print(f"Offset temperature controller turned OFF at t = {current_time_ps:.2f} ps")
        
        # Skip if not active
        if not self.is_active:
            return
        
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        self.last_update_time = current_time_ps
        
        # Calculate measured system temperature
        measured_temperature = self._calculate_system_temperature(current_time_ps)
        
        if measured_temperature is None:
            return
        
        # Calculate target bath temperature: T_bath = T_measured + offset
        target_bath_temperature = measured_temperature + self.temperature_offset_K
        
        # Apply temperature limits
        if self.T_max is not None:
            actual_bath_temperature = max(self.T_min, min(self.T_max, target_bath_temperature))
        else:
            actual_bath_temperature = max(self.T_min, target_bath_temperature)
        
        # Update thermostats
        self.current_bath_temperature = actual_bath_temperature
        self._update_thermostats(actual_bath_temperature)
        
        # Console output
        if current_time_ps - self.last_console_output_time >= self.console_output_period_ps:
            print(f"Offset Controller | {self.temperature_method}→{measured_temperature:.1f}K | "
                  f"Bath: {actual_bath_temperature:.1f}K (offset: {self.temperature_offset_K:+.1f}K)")
            self.last_console_output_time = current_time_ps
        
        # Log data
        try:
            with open(self.output_file, 'a', encoding='utf-8') as f:
                f.write(f"{current_time_ps:.6f},{measured_temperature:.6f},{target_bath_temperature:.6f},"
                       f"{actual_bath_temperature:.6f},{self.temperature_offset_K:.6f},{int(self.is_active)}\n")
        except Exception as e:
            print(f"Warning: Failed to write offset controller output: {e}")


class DiffEqController(hoomd.custom.Action):
    """
    Differential equation temperature controller for cavity MD simulations.
    
    This controller implements a first-order differential equation to control
    the bath temperature based on a signal temperature:
    
    dT_bath(t)/dt = -(T_bath(t) - T_signal(t)) / τ
    
    where:
    - T_bath(t) is the bath temperature
    - T_signal(t) is the measured signal temperature (kinetic, LJ+Coulombic, etc.)
    - τ is the time constant controlling response speed
    
    The discretized update rule is:
    T_bath(t+dt) = T_bath(t) + dt * (-(T_bath(t) - T_signal(t)) / τ)
    T_bath(t+dt) = T_bath(t) - (dt/τ) * (T_bath(t) - T_signal(t))
    
    Parameters
    ----------
    temperature_method : str
        Temperature calculation method for signal ('kinetic', 'lj_coulombic', 'harmonic', 'harmonic_equipartition', 'lj_coulombic_bath')
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
    rate_limit_K_per_ps : float, optional
        Maximum rate of temperature change in K/ps. Default: None (no rate limit)
    output_file : str, optional
        Output file for control data. Default: 'diffeq_feedback.csv'
    empirical_data_file : str, optional
        Path to empirical data file (required for certain methods)
    console_output_period_ps : float, optional
        Console output period in picoseconds. Default: 1.0
        
    Examples
    --------
    **Basic differential equation control following kinetic temperature:**
    
    >>> from hoomd.cavitymd.analysis import DiffEqController
    >>> 
    >>> diffeq_controller = DiffEqController(
    ...     temperature_method='kinetic',
    ...     time_constant_ps=5.0,  # 5 ps response time
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     apply_to='molecular'
    ... )
    
    **Fast response to LJ+Coulombic signal:**
    
    >>> diffeq_controller = DiffEqController(
    ...     temperature_method='lj_coulombic_bath',
    ...     time_constant_ps=1.0,  # Fast 1 ps response
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     empirical_data_file='lj_coulombic_calibration.txt',
    ...     apply_to='both'
    ... )
    
    Notes
    -----
    - Smaller time constants lead to faster response to signal changes
    - The bath temperature exponentially approaches the signal temperature
    - Unlike feedback controllers, this directly tracks the signal without a target
    - Compatible with all existing temperature calculation methods
    """
    
    def __init__(self,
                 temperature_method: str,
                 time_constant_ps: float,
                 time_tracker,
                 energy_tracker,
                 simulation=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 apply_to: str = 'both',
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 T_min: float = 0.0,
                 T_max: Optional[float] = None,
                 rate_limit_K_per_ps: Optional[float] = None,
                 output_file: str = 'diffeq_feedback.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0):
        
        # Store core parameters
        self.temperature_method = temperature_method
        self.time_constant_ps = float(time_constant_ps)
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        
        # Control parameters
        self.apply_to = apply_to
        self.turn_on_time_ps = float(turn_on_time_ps)
        self.turn_off_time_ps = float(turn_off_time_ps) if turn_off_time_ps is not None else None
        self.update_interval_ps = float(update_interval_ps)
        self.T_min = float(T_min)
        self.T_max = float(T_max) if T_max is not None else None
        self.rate_limit_K_per_ps = float(rate_limit_K_per_ps) if rate_limit_K_per_ps is not None else None
        self.output_file = output_file
        self.console_output_period_ps = float(console_output_period_ps)
        
        # Load empirical data if needed
        self.empirical_data = None
        if empirical_data_file is not None and temperature_method in ['lj_coulombic', 'harmonic', 'lj_coulombic_bath']:
            try:
                energy_component = 'lj_coulombic' if 'lj_coulombic' in temperature_method else 'harmonic'
                self.empirical_data = EmpiricalTemperatureData(
                    data_file_path=empirical_data_file,
                    energy_component=energy_component
                )
            except Exception as e:
                print(f"Warning: Could not load empirical data file {empirical_data_file}: {e}")
        
        # Control state
        self.is_active = False
        self.current_bath_temperature = None
        self.last_update_time = 0.0
        self.last_console_output_time = 0.0
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("# Differential Equation Temperature Controller Output\n")
                f.write("# time_ps: simulation time in picoseconds\n")
                f.write("# signal_temp_K: measured signal temperature\n")
                f.write("# bath_temp_K: bath temperature setting\n")
                f.write("# dt_ps: timestep used for integration\n")
                f.write("# derivative_K_per_ps: dT_bath/dt calculated\n")
                f.write("# active: whether controller is active (1=active, 0=inactive)\n")
                f.write("time_ps,signal_temp_K,bath_temp_K,dt_ps,derivative_K_per_ps,active\n")
        except Exception as e:
            print(f"Warning: Failed to initialize differential equation controller output file: {e}")
    
    def _calculate_system_temperature(self, current_time_ps: float) -> Optional[float]:
        """Calculate system temperature using the specified method."""
        # Use the same implementation as OffsetTemperatureController
        if self.temperature_method == 'kinetic':
            # Calculate kinetic temperature from particle velocities
            import numpy as np
            from .utils import PhysicalConstants
            
            try:
                kinetic_energy = 0.0
                N_molecules = 0
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
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
        
        elif self.temperature_method == 'harmonic_equipartition':
            # Calculate harmonic temperature using equipartition theorem
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                harmonic_energy_hartree = energy_data.get('harmonic', 0.0)
                
                # Use equipartition theorem: E_harmonic = (3/2) * N * kB * T
                # Assume N=500 molecular particles (excluding cavity)
                N_molecules = 500
                from .utils import PhysicalConstants
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                
                if harmonic_energy_hartree > 0:
                    temperature_K = (2.0/3.0) * harmonic_energy_hartree / (N_molecules * kB_hartree)
                    return max(temperature_K, 0.0)
                else:
                    return None
                    
            except Exception as e:
                print(f"Warning: Could not calculate harmonic equipartition temperature: {e}")
                return None
        
        elif self.temperature_method in ['lj_coulombic', 'lj_coulombic_bath']:
            # Calculate LJ+Coulombic fictive temperature using empirical data
            if self.empirical_data is None:
                print(f"Warning: No empirical data available for method '{self.temperature_method}'")
                return None
            
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                lj_energy = energy_data.get('lj', 0.0)
                coulomb_energy = energy_data.get('coulombic', 0.0)
                total_lj_coulomb = lj_energy + coulomb_energy
                
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    instantaneous_energy_hartree=total_lj_coulomb,
                    num_particles=500
                )
                return max(temperature_K, 0.0)
                
            except Exception as e:
                print(f"Warning: Could not calculate LJ+Coulombic temperature: {e}")
                return None
        
        elif self.temperature_method == 'harmonic':
            # Calculate harmonic fictive temperature using empirical data
            if self.empirical_data is None:
                print(f"Warning: No empirical data available for method '{self.temperature_method}'")
                return None
            
            try:
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if energy_data is None:
                    return None
                
                harmonic_energy = energy_data.get('harmonic', 0.0)
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    instantaneous_energy_hartree=harmonic_energy,
                    num_particles=500
                )
                return max(temperature_K, 0.0)
                
            except Exception as e:
                print(f"Warning: Could not calculate harmonic temperature: {e}")
                return None
        
        else:
            print(f"Warning: Unknown temperature method '{self.temperature_method}'")
            return None
    
    def _get_current_thermostat_temperature(self) -> float:
        """Get current thermostat temperature for initialization."""
        try:
            from .utils import PhysicalConstants
            
            # Try molecular thermostat first
            if self.molecular_thermostat is not None:
                try:
                    if hasattr(self.molecular_thermostat, 'kT'):
                        kT = self.molecular_thermostat.kT
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.molecular_thermostat, 'thermostat'):
                        nested = self.molecular_thermostat.thermostat
                        if hasattr(nested, 'kT'):
                            kT = nested.kT
                            if hasattr(kT, 'value'):
                                kT = kT.value
                            elif callable(kT):
                                kT = kT()
                            return kT / PhysicalConstants.KB_HARTREE_PER_K
                except:
                    pass
            
            # Try cavity thermostat as backup
            if self.cavity_thermostat is not None:
                try:
                    if hasattr(self.cavity_thermostat, 'kT'):
                        kT = self.cavity_thermostat.kT
                        if hasattr(kT, 'value'):
                            kT = kT.value
                        elif callable(kT):
                            kT = kT()
                        return kT / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.cavity_thermostat, 'thermostat'):
                        nested = self.cavity_thermostat.thermostat
                        if hasattr(nested, 'kT'):
                            kT = nested.kT
                            if hasattr(kT, 'value'):
                                kT = kT.value
                            elif callable(kT):
                                kT = kT()
                            return kT / PhysicalConstants.KB_HARTREE_PER_K
                except:
                    pass
            
            # Fallback to 300K
            return 300.0
            
        except Exception as e:
            print(f"Warning: Could not get current thermostat temperature: {e}")
            return 300.0

    def _apply_rate_limit(self, new_output: float, dt_ps: float) -> float:
        """Apply rate limiting to controller output."""
        if self.rate_limit_K_per_ps is None or self.current_bath_temperature is None:
            return new_output
        
        max_change = self.rate_limit_K_per_ps * dt_ps
        change = new_output - self.current_bath_temperature
        
        if abs(change) <= max_change:
            return new_output
        else:
            # Limit the change
            limited_change = max_change if change > 0 else -max_change
            return self.current_bath_temperature + limited_change
    
    def _update_thermostats(self, bath_temperature: float):
        """Update thermostat temperatures."""
        try:
            import hoomd
            from .utils import PhysicalConstants
            kT_hartree = bath_temperature * PhysicalConstants.KB_HARTREE_PER_K
            
            if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
                if hasattr(self.molecular_thermostat, 'kT'):
                    self.molecular_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                elif hasattr(self.molecular_thermostat, 'thermostat'):
                    nested = self.molecular_thermostat.thermostat
                    if hasattr(nested, 'kT'):
                        nested.kT = hoomd.variant.Constant(kT_hartree)
            
            if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
                if hasattr(self.cavity_thermostat, 'kT'):
                    self.cavity_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                elif hasattr(self.cavity_thermostat, 'thermostat'):
                    nested = self.cavity_thermostat.thermostat
                    if hasattr(nested, 'kT'):
                        nested.kT = hoomd.variant.Constant(kT_hartree)
                        
        except Exception as e:
            print(f"Warning: Failed to update thermostats: {e}")
            import traceback
            traceback.print_exc()
    
    def act(self, timestep):
        """Execute differential equation control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should be active
        should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                          (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
        
        if should_be_active and not self.is_active:
            self.is_active = True
            # Initialize current_bath_temperature from thermostats when first activating
            if self.current_bath_temperature is None:
                self.current_bath_temperature = self._get_current_thermostat_temperature()
                print(f"Differential equation controller turned ON at t = {current_time_ps:.2f} ps")
                print(f"Initial bath temperature: {self.current_bath_temperature:.2f} K")
                print(f"Time constant: {self.time_constant_ps:.2f} ps")
            else:
                print(f"Differential equation controller turned ON at t = {current_time_ps:.2f} ps")
        elif not should_be_active and self.is_active:
            self.is_active = False
            print(f"Differential equation controller turned OFF at t = {current_time_ps:.2f} ps")
        
        # Skip if not active or not initialized
        if not self.is_active or self.current_bath_temperature is None:
            return
        
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        # Calculate dt_ps BEFORE updating last_update_time
        dt_ps = current_time_ps - self.last_update_time if self.last_update_time > 0 else self.update_interval_ps
        self.last_update_time = current_time_ps
        
        # Calculate signal temperature
        signal_temperature = self._calculate_system_temperature(current_time_ps)
        
        if signal_temperature is None:
            return
        
        # Implement differential equation: dT_bath/dt = -(T_bath - T_signal) / τ
        # Discretized: T_bath(t+dt) = T_bath(t) + dt * (-(T_bath(t) - T_signal(t)) / τ)
        derivative_K_per_ps = -(self.current_bath_temperature - signal_temperature) / self.time_constant_ps
        raw_output = self.current_bath_temperature + dt_ps * derivative_K_per_ps
        
        # Apply saturation limits
        if self.T_max is not None:
            saturated_output = max(self.T_min, min(self.T_max, raw_output))
        else:
            saturated_output = max(self.T_min, raw_output)
        
        # Apply rate limiting
        rate_limited_output = self._apply_rate_limit(saturated_output, dt_ps)
        
        # Update current bath temperature and thermostats
        self.current_bath_temperature = rate_limited_output
        self._update_thermostats(self.current_bath_temperature)
        
        # Console output
        if current_time_ps - self.last_console_output_time >= self.console_output_period_ps:
            print(f"DiffEq Controller | {self.temperature_method}→{signal_temperature:.1f}K | "
                  f"Bath: {self.current_bath_temperature:.1f}K | dT/dt: {derivative_K_per_ps:+.2f}K/ps")
            self.last_console_output_time = current_time_ps
        
        # Log data
        try:
            with open(self.output_file, 'a', encoding='utf-8') as f:
                f.write(f"{current_time_ps:.6f},{signal_temperature:.6f},{self.current_bath_temperature:.6f},"
                       f"{dt_ps:.6f},{derivative_K_per_ps:.6f},{int(self.is_active)}\n")
        except Exception as e:
            print(f"Warning: Failed to write differential equation controller output: {e}")


class SinusoidalBathTemperatureController(hoomd.custom.Action):
    """
    Sinusoidal bath temperature controller for cavity MD simulations.
    
    This controller modulates the bath temperature sinusoidally around a dynamic target
    temperature, with the amplitude dynamically adjusted based on a configurable temperature
    method. This provides a temperature control strategy that oscillates around
    equilibrium while responding to different aspects of the system's thermal state.
    
    The temperature schedule follows:
    T_bath(t) = T_mean + A(t) * sin(2π * t / T_period + φ)
    
    where:
    - T_mean: Dynamic target temperature (set based on initial system state)
    - A(t): Dynamic amplitude based on specified temperature method
    - T_period: Period of sinusoidal oscillation
    - φ: Phase offset
    
    The amplitude is updated according to:
    A(t) = amplitude_scale * |T_method(t) - T_mean|
    
    Parameters
    ----------
    period_ps : float
        Period of sinusoidal temperature oscillation in picoseconds
    amplitude_scale : float
        Scaling factor for amplitude based on harmonic equipartition deviation
    phase_offset : float
        Phase offset for sinusoidal function (in radians)
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    energy_tracker : EnergyTracker
        Energy tracker for temperature calculations
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat
        Molecular thermostat to control
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat  
        Cavity thermostat to control
    apply_to : str
        Which thermostats to control ('molecular', 'cavity', 'both')
    target_temperature : float, optional
        Static target temperature in Kelvin. Ignored if dynamic_target=True
    dynamic_target : bool, optional
        Whether to set mean temperature dynamically at turn-on time. Default: True
    turn_on_time_ps : float, optional
        Time to start control in picoseconds. Default: 0.0
    turn_off_time_ps : float, optional
        Time to stop control in picoseconds. Default: None (never turn off)
    update_interval_ps : float, optional
        Interval between control updates in picoseconds. Default: 0.1
    T_min : float, optional
        Minimum allowed bath temperature in Kelvin. Default: 0.1
    T_max : float, optional
        Maximum allowed bath temperature in Kelvin. Default: None
    output_file : str, optional
        Output file for control data. Default: 'sinusoidal_bath_controller.csv'
    empirical_data_file : str, optional
        Path to empirical data file for harmonic equipartition calculation
    console_output_period_ps : float, optional
        Console output period in picoseconds. Default: 1.0
    amplitude_update_interval_ps : float, optional
        Interval for updating amplitude based on harmonic temperature. Default: 1.0
    amplitude_temperature_method : str, optional
        Temperature method to use for amplitude calculation. Options: 'harmonic_equipartition', 
        'kinetic', 'lj_coulombic', 'harmonic'. Default: 'harmonic_equipartition'
    adaptive_range_mode : bool, optional
        If True, uses adaptive range mode where sinusoid bounds are fixed at target temperature
        and range adapts to |Ts(t) - Ttarget|. If False, uses fixed mean with variable amplitude.
        Default: False
        
    Examples
    --------
    **Basic sinusoidal bath control:**
    
    >>> controller = SinusoidalBathTemperatureController(
    ...     period_ps=2.0,
    ...     amplitude_scale=0.1,
    ...     phase_offset=0.0,
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     apply_to='both',
    ...     dynamic_target=True,
    ...     turn_on_time_ps=100.0
    ... )
    
    Notes
    -----
    - Mean temperature is set to the current bath temperature at turn-on if dynamic_target=True
    - Amplitude responds to deviations in harmonic equipartition temperature from the mean
    - Both thermostats oscillate in phase (same sinusoidal schedule)
    - Temperature limits are applied after sinusoidal calculation
    """
    
    def __init__(self,
                 period_ps: float,
                 amplitude_scale: float,
                 phase_offset: float = 0.0,
                 time_tracker=None,
                 energy_tracker=None,
                 simulation=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 apply_to: str = 'both',
                 target_temperature: float = 100.0,
                 dynamic_target: bool = True,
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 T_min: float = 0.1,
                 T_max: Optional[float] = None,
                 output_file: str = 'sinusoidal_bath_controller.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0,
                 amplitude_update_interval_ps: float = 1.0,
                 amplitude_temperature_method: str = 'harmonic_equipartition',
                 adaptive_range_mode: bool = False):
        
        # Store configuration
        self.period_ps = period_ps
        self.amplitude_scale = amplitude_scale
        self.phase_offset = phase_offset
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.apply_to = apply_to
        self.target_temperature = target_temperature
        self.dynamic_target = dynamic_target
        self.turn_on_time_ps = turn_on_time_ps
        self.turn_off_time_ps = turn_off_time_ps
        self.update_interval_ps = update_interval_ps
        self.T_min = T_min
        self.T_max = T_max
        self.output_file = output_file
        self.empirical_data_file = empirical_data_file
        self.console_output_period_ps = console_output_period_ps
        self.amplitude_update_interval_ps = amplitude_update_interval_ps
        self.amplitude_temperature_method = amplitude_temperature_method
        self.adaptive_range_mode = adaptive_range_mode
        
        # Load empirical data for amplitude temperature calculation
        self.empirical_data = None
        if empirical_data_file:
            try:
                # Determine which energy component to use based on amplitude temperature method
                if self.amplitude_temperature_method in ['lj_coulombic']:
                    energy_component = 'lj_coulombic'
                elif self.amplitude_temperature_method in ['harmonic']:
                    energy_component = 'harmonic'
                else:
                    # For kinetic and harmonic_equipartition, we don't need empirical data
                    energy_component = None
                
                if energy_component:
                    self.empirical_data = EmpiricalTemperatureData(empirical_data_file, energy_component)
                    print(f"Loaded empirical data from {empirical_data_file} for {energy_component} energy")
                else:
                    print(f"Empirical data not needed for amplitude method '{self.amplitude_temperature_method}'")
                    
            except Exception as e:
                print(f"Warning: Could not load empirical data from {empirical_data_file}: {e}")
                self.empirical_data = None
        
        # State variables
        self.is_active = False
        self.mean_temperature = target_temperature
        self.current_amplitude = 0.0
        self.dynamic_target_set = False
        self.last_update_time = 0.0
        self.last_amplitude_update_time = 0.0
        self.last_console_output_time = 0.0
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("time_ps,target_temp_K,current_amplitude_K,amplitude_method_temp_K,"
                       "calculated_bath_temp_K,sinusoidal_component_K,is_active\n")
            print(f"Initialized sinusoidal bath controller output: {self.output_file}")
        except Exception as e:
            print(f"Warning: Failed to initialize sinusoidal bath controller output file: {e}")
    
    def _get_thermostat_temperature(self, thermostat) -> Optional[float]:
        """Get current temperature from thermostat in Kelvin."""
        if thermostat is None:
            return None
        
        try:
            from .utils import PhysicalConstants
            
            # Try different thermostat structures
            kT_hartree = None
            if hasattr(thermostat, 'kT'):
                # Direct kT attribute (Langevin thermostats)
                if hasattr(thermostat.kT, 'value'):
                    kT_hartree = thermostat.kT.value
                else:
                    kT_hartree = thermostat.kT
            elif hasattr(thermostat, 'thermostat'):
                # Nested thermostat (ConstantVolume with Bussi/MTTK)
                nested_thermo = thermostat.thermostat
                if hasattr(nested_thermo, 'kT'):
                    if hasattr(nested_thermo.kT, 'value'):
                        kT_hartree = nested_thermo.kT.value
            else:
                        kT_hartree = nested_thermo.kT
            
            if kT_hartree is not None:
                return kT_hartree / PhysicalConstants.KB_HARTREE_PER_K
            else:
                return None
                
        except Exception as e:
            print(f"Warning: Failed to get thermostat temperature: {e}")
            return None
    
    def _calculate_amplitude_temperature(self, current_time_ps: float) -> Optional[float]:
        """Calculate temperature using the specified amplitude temperature method."""
        if self.energy_tracker is None:
            return None
        
        try:
            method = self.amplitude_temperature_method
            
            if method == 'harmonic_equipartition':
                # Get harmonic energy from energy tracker
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                
                if harmonic_energy <= 0:
                    return None
                
                # Get number of particles
                from .utils import PhysicalConstants
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    N_particles = len(snap.particles.mass)
                
                if N_particles == 0:
                    return None
                
                # Direct equipartition calculation: T = (4*E_harmonic)/(N*kB)
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (4.0 * harmonic_energy) / (N_particles * kB_hartree)
                
                return max(temperature_K, 0.0)
                
            elif method == 'kinetic':
                # Calculate kinetic temperature
                with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                    velocities = snap.particles.velocity
                    masses = snap.particles.mass
                    
                # Calculate kinetic energy per particle: (1/2) * m * v^2
                kinetic_energies = 0.5 * masses[:, np.newaxis] * (velocities ** 2)
                total_kinetic_energy = np.sum(kinetic_energies)
                
                # Convert to temperature using equipartition theorem: (3/2) * N * kB * T = KE
                N_particles = len(masses)
                if N_particles == 0:
                    return None
                    
                from .utils import PhysicalConstants
                kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
                temperature_K = (2.0 * total_kinetic_energy) / (3.0 * N_particles * kB_hartree)
                
                return max(temperature_K, 0.0)
                
            elif method == 'lj_coulombic':
                # Calculate LJ + Coulombic temperature using empirical data
                if self.empirical_data is None:
                    print(f"Warning: No empirical data available for method '{method}'")
                    return None
                
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                lj_energy = energy_dict.get('lj', 0.0)
                coulombic_energy = energy_dict.get('coulombic', 0.0)
                total_energy = lj_energy + coulombic_energy
                
                # Use empirical data to convert energy to temperature
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    instantaneous_energy_hartree=total_energy,
                    num_particles=500  # Assume 500 molecular particles
                )
                return max(temperature_K, 0.0)
                
            elif method == 'harmonic':
                # Calculate harmonic temperature using empirical data
                if self.empirical_data is None:
                    print(f"Warning: No empirical data available for method '{method}'")
                    return None
                
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                
                # Use empirical data to convert energy to temperature
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    instantaneous_energy_hartree=harmonic_energy,
                    num_particles=500  # Assume 500 molecular particles
                )
                return max(temperature_K, 0.0)
                
            else:
                print(f"Warning: Unknown amplitude temperature method '{method}'")
                return None
                
        except Exception as e:
            print(f"Warning: Failed to calculate amplitude temperature using method '{method}': {e}")
            return None
    
    def _set_dynamic_target(self, current_time_ps: float):
        """Set dynamic target temperature based on current system state."""
        if not self.dynamic_target or self.dynamic_target_set:
            return
        
        # Try to get current bath temperature from either thermostat
        bath_temp = None
        if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
            bath_temp = self._get_thermostat_temperature(self.molecular_thermostat)
        
        if bath_temp is None and self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
            bath_temp = self._get_thermostat_temperature(self.cavity_thermostat)
        
        if bath_temp is not None:
            if self.adaptive_range_mode:
                # In adaptive range mode, set target_temperature and let _update_adaptive_range_parameters calculate mean
                self.target_temperature = bath_temp
                print(f"Sinusoidal bath controller: Dynamic target set to {self.target_temperature:.2f} K at t = {current_time_ps:.2f} ps")
            else:
                # In fixed mean mode, set mean_temperature directly
                self.mean_temperature = bath_temp
                print(f"Sinusoidal bath controller: Dynamic target set to {self.mean_temperature:.2f} K at t = {current_time_ps:.2f} ps")
            self.dynamic_target_set = True
        else:
            print(f"Warning: Could not determine bath temperature for dynamic target, using {self.target_temperature:.2f} K")
            if not self.adaptive_range_mode:
                self.mean_temperature = self.target_temperature
            self.dynamic_target_set = True
    
    def _update_amplitude(self, current_time_ps: float):
        """Update amplitude based on specified temperature method and mode."""
        amplitude_temp = self._calculate_amplitude_temperature(current_time_ps)
        
        if self.adaptive_range_mode:
            self._update_adaptive_range_parameters(amplitude_temp, current_time_ps)
        else:
            # Original mode: fixed mean, variable amplitude
            if amplitude_temp is not None:
                # Amplitude scales with deviation from mean temperature
                deviation = abs(amplitude_temp - self.mean_temperature)
                self.current_amplitude = self.amplitude_scale * deviation
            else:
                # Fall back to a small default amplitude if calculation fails
                self.current_amplitude = self.amplitude_scale * 5.0  # 5K typical deviation
    
    def _update_adaptive_range_parameters(self, signal_temp: Optional[float], current_time_ps: float):
        """Update parameters for adaptive range mode.
        
        In adaptive range mode:
        - Range (max - min) = amplitude_scale * |Ts(t) - Ttarget|
        - If Ts(t) > Ttarget: max = Ttarget, min = Ttarget - |Ts(t) - Ttarget| (dip down)
        - If Ts(t) < Ttarget: min = Ttarget, max = Ttarget + |Ts(t) - Ttarget| (reach up)
        - Target temperature is always one of the bounds (min or max)
        """
        if signal_temp is None:
            # Fallback: small range around target
            self.current_amplitude = self.amplitude_scale * 5.0  # 5K default
            self.mean_temperature = self.target_temperature
            return
        
        # Calculate the deviation from target
        deviation = signal_temp - self.target_temperature
        range_magnitude = self.amplitude_scale * abs(deviation)
        
        # Set amplitude (half the range since sinusoid goes ±amplitude around mean)
        self.current_amplitude = range_magnitude / 2.0
        
        # Adjust mean temperature based on deviation direction
        if deviation > 0:
            # Signal too hot: max = target, min = target - range
            # Mean = (max + min)/2 = (target + target - range)/2 = target - range/2
            proposed_mean = self.target_temperature - self.current_amplitude
            # Ensure mean doesn't go below T_min (which would create negative temperatures)
            if proposed_mean < self.T_min:
                # Adjust the range to keep mean at T_min
                self.current_amplitude = self.target_temperature - self.T_min
                self.mean_temperature = self.T_min
            else:
                self.mean_temperature = proposed_mean
        elif deviation < 0:
            # Signal too cold: min = target, max = target + range
            # Mean = (max + min)/2 = (target + range + target)/2 = target + range/2
            self.mean_temperature = self.target_temperature + self.current_amplitude
        else:
            # Perfect match: oscillate around target
            self.mean_temperature = self.target_temperature
    
    def _update_thermostats(self, bath_temperature: float):
        """Update thermostat temperatures."""
        try:
            from .utils import PhysicalConstants
            import hoomd
            
            # Convert temperature to kT in atomic units
            kT_hartree = bath_temperature * PhysicalConstants.KB_HARTREE_PER_K
            
            # Update molecular thermostat
            if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
                if hasattr(self.molecular_thermostat, 'kT'):
                    # Direct kT attribute (Langevin thermostats)
                    self.molecular_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                elif hasattr(self.molecular_thermostat, 'thermostat'):
                    # Nested thermostat (ConstantVolume with Bussi/MTTK)
                    nested_thermo = self.molecular_thermostat.thermostat
                    if hasattr(nested_thermo, 'kT'):
                        nested_thermo.kT = hoomd.variant.Constant(kT_hartree)
            
            # Update cavity thermostat
            if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
                if hasattr(self.cavity_thermostat, 'kT'):
                    # Direct kT attribute (Langevin thermostats)
                    self.cavity_thermostat.kT = hoomd.variant.Constant(kT_hartree)
                elif hasattr(self.cavity_thermostat, 'thermostat'):
                    # Nested thermostat (ConstantVolume with Bussi/MTTK)
                    nested_thermo = self.cavity_thermostat.thermostat
                    if hasattr(nested_thermo, 'kT'):
                        nested_thermo.kT = hoomd.variant.Constant(kT_hartree)
                        
        except Exception as e:
            print(f"Warning: Failed to update thermostats in sinusoidal controller: {e}")
    
    def act(self, timestep):
        """Execute sinusoidal bath temperature control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should be active
        should_be_active = (current_time_ps >= self.turn_on_time_ps and 
                          (self.turn_off_time_ps is None or current_time_ps < self.turn_off_time_ps))
        
        if should_be_active and not self.is_active:
            self.is_active = True
            print(f"Sinusoidal bath temperature controller turned ON at t = {current_time_ps:.2f} ps")
            print(f"Period: {self.period_ps:.2f} ps, Amplitude scale: {self.amplitude_scale:.3f}")
            
            # Set dynamic target if enabled
            self._set_dynamic_target(current_time_ps)
            
        elif not should_be_active and self.is_active:
            self.is_active = False
            print(f"Sinusoidal bath temperature controller turned OFF at t = {current_time_ps:.2f} ps")
        
        # Skip if not active
        if not self.is_active:
            return
        
        # Only update at specified intervals
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        self.last_update_time = current_time_ps
        
        # Get amplitude temperature BEFORE changing bath temperature
        amplitude_temp = self._calculate_amplitude_temperature(current_time_ps)
        
        # Update amplitude at specified intervals
        if current_time_ps - self.last_amplitude_update_time >= self.amplitude_update_interval_ps:
            self._update_amplitude(current_time_ps)
            self.last_amplitude_update_time = current_time_ps
        
        # Calculate sinusoidal temperature
        import math
        time_phase = 2.0 * math.pi * current_time_ps / self.period_ps + self.phase_offset
        sinusoidal_component = self.current_amplitude * math.sin(time_phase)
        target_bath_temperature = self.mean_temperature + sinusoidal_component
        
        # Apply temperature limits with strict T_min constraint
        if self.T_max is not None:
            actual_bath_temperature = max(self.T_min, min(self.T_max, target_bath_temperature))
        else:
            actual_bath_temperature = max(self.T_min, target_bath_temperature)
        
        # Additional safety check: ensure temperature is always above minimum
        if actual_bath_temperature < self.T_min:
            print(f"Warning: Temperature {actual_bath_temperature:.2f}K below minimum {self.T_min:.2f}K, "
                  f"clamping to minimum at t={current_time_ps:.2f} ps")
            actual_bath_temperature = self.T_min
        
        # Update thermostats
        self._update_thermostats(actual_bath_temperature)
        
        # Console output
        if current_time_ps - self.last_console_output_time >= self.console_output_period_ps:
            temp_method_display = self.amplitude_temperature_method.replace('_', ' ').title()
            print(f"Sinusoidal Bath | Mean: {self.mean_temperature:.1f}K | "
                  f"Amplitude: {self.current_amplitude:.2f}K | Bath: {actual_bath_temperature:.1f}K | "
                  f"{temp_method_display}: {amplitude_temp:.1f}K" if amplitude_temp else f"Bath: {actual_bath_temperature:.1f}K")
            self.last_console_output_time = current_time_ps
        
        # Log data
        try:
            with open(self.output_file, 'a', encoding='utf-8') as f:
                f.write(f"{current_time_ps:.6f},{target_bath_temperature:.6f},"
                       f"{self.current_amplitude:.6f},{amplitude_temp or 0.0:.6f},"
                       f"{actual_bath_temperature:.6f},{sinusoidal_component:.6f},{int(self.is_active)}\n")
        except Exception as e:
            print(f"Warning: Failed to write sinusoidal bath controller output: {e}")


class AdaptiveBathTemperatureController(hoomd.custom.Action):
    """
    Adaptive bath temperature controller for cavity MD simulations using GD framework.
    
    This controller implements a gradient descent approach where the bath temperature
    evolves according to a differential equation toward a dynamically calculated target:
    
    dTb(t)/dt = -(Tb(t) - Tstar(t)) / tau
    
    where Tstar(t) is the dynamic target temperature calculated as:
    - If Ts(t) > Ttarget: Tstar = Ttarget - amplitude_scale * |Ts(t) - Ttarget| (cool down)  
    - If Ts(t) < Ttarget: Tstar = Ttarget + amplitude_scale * |Ts(t) - Ttarget| (heat up)
    
    This provides smooth temperature evolution with configurable response time (tau).
    
    Parameters
    ----------
    amplitude_scale : float
        Scaling factor for bath temperature adjustment based on temperature error
    time_constant_ps : float
        Time constant tau for GD evolution in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    energy_tracker : EnergyTracker
        Energy tracker for signal temperature calculations
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat
        Molecular thermostat to control
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat  
        Cavity thermostat to control
    apply_to : str
        Which thermostats to control ('molecular', 'cavity', 'both')
    target_temperature : float, optional
        Static target temperature in Kelvin. Ignored if dynamic_target=True
    dynamic_target : bool, optional
        Whether to set target temperature dynamically at turn-on time. Default: True
    turn_on_time_ps : float, optional
        Time to start control in picoseconds. Default: 0.0
    turn_off_time_ps : float, optional
        Time to stop control in picoseconds. Default: None (never turn off)
    update_interval_ps : float, optional
        Interval between control updates in picoseconds. Default: 0.1
    T_min : float, optional
        Minimum allowed bath temperature in Kelvin. Default: 0.1
    T_max : float, optional
        Maximum allowed bath temperature in Kelvin. Default: None
    output_file : str, optional
        Output file for control data. Default: 'adaptive_bath_controller.csv'
    empirical_data_file : str, optional
        Path to empirical data file for signal temperature calculation
    console_output_period_ps : float, optional
        Console output period in picoseconds. Default: 1.0
    signal_temperature_method : str, optional
        Temperature method for control signal calculation. Options: 'harmonic_equipartition', 
        'kinetic', 'lj_coulombic', 'harmonic'. Default: 'harmonic_equipartition'
    dynamic_target_temperature_method : str, optional
        Temperature method for setting the dynamic target at turn-on. 
        If None, uses signal_temperature_method. Options: 'harmonic_equipartition', 
        'kinetic', 'lj_coulombic', 'harmonic'. Default: None (same as signal_temperature_method)
    """
    
    def __init__(self,
                 amplitude_scale: float,
                 time_constant_ps: float,
                 time_tracker=None,
                 energy_tracker=None,
                 simulation=None,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 apply_to: str = 'both',
                 target_temperature: float = 100.0,
                 dynamic_target: bool = True,
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 0.1,
                 T_min: float = 0.1,
                 T_max: Optional[float] = None,
                 output_file: str = 'adaptive_bath_controller.csv',
                 empirical_data_file: Optional[str] = None,
                 console_output_period_ps: float = 1.0,
                 signal_temperature_method: str = 'harmonic_equipartition',
                 dynamic_target_temperature_method: Optional[str] = None):
        
        # Store configuration
        self.amplitude_scale = amplitude_scale
        self.time_constant_ps = time_constant_ps
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.simulation = simulation
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.apply_to = apply_to
        self.target_temperature = target_temperature
        self.dynamic_target = dynamic_target
        self.turn_on_time_ps = turn_on_time_ps
        self.turn_off_time_ps = turn_off_time_ps
        self.update_interval_ps = update_interval_ps
        self.T_min = T_min
        self.T_max = T_max
        self.output_file = output_file
        self.empirical_data_file = empirical_data_file
        self.console_output_period_ps = console_output_period_ps
        self.signal_temperature_method = signal_temperature_method
        # Default dynamic target method to signal method if not specified
        self.dynamic_target_temperature_method = dynamic_target_temperature_method if dynamic_target_temperature_method is not None else signal_temperature_method
        
        # Initialize state
        self.is_active = False
        self.dynamic_target_set = False
        self.last_update_time = -1.0
        self.last_console_output_time = -1.0
        self.current_bath_temperature = None  # Will be initialized at turn-on
        
        # Load empirical data if needed (for either control signal or dynamic target)
        # Need separate empirical data objects for different energy components
        self.empirical_data = None  # For lj_coulombic (backward compatibility)
        self.empirical_data_harmonic = None  # For harmonic
        
        if empirical_data_file:
            try:
                # Load LJ+Coulombic empirical data if needed
                if signal_temperature_method == 'lj_coulombic' or self.dynamic_target_temperature_method == 'lj_coulombic':
                    self.empirical_data = EmpiricalTemperatureData(empirical_data_file, energy_component='lj_coulombic')
                    print(f"Loaded LJ+Coulombic empirical data")
                
                # Load harmonic empirical data if needed
                if signal_temperature_method == 'harmonic' or self.dynamic_target_temperature_method == 'harmonic':
                    self.empirical_data_harmonic = EmpiricalTemperatureData(empirical_data_file, energy_component='harmonic')
                    print(f"Loaded harmonic empirical data")
                    
            except Exception as e:
                print(f"Warning: Failed to load empirical data from {empirical_data_file}: {e}")
        
        # Initialize output file
        self._initialize_output_file()
    
    def _initialize_output_file(self):
        """Initialize the output CSV file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("# Adaptive Bath Temperature Controller Output (GD Framework)\\n")
                f.write("# time_ps: simulation time in picoseconds\\n")
                f.write("# signal_temp_K: signal temperature used for control\\n")
                f.write("# target_temp_K: target temperature\\n") 
                f.write("# Tstar_K: dynamic target temperature for GD evolution\\n")
                f.write("# bath_temp_K: current bath temperature\\n")
                f.write("# temp_error_K: temperature error (signal - target)\\n")
                f.write("# dTb_dt: bath temperature time derivative\\n")
                f.write("# is_active: whether controller is active (1=active, 0=inactive)\\n")
                f.write("time_ps,signal_temp_K,target_temp_K,Tstar_K,bath_temp_K,temp_error_K,dTb_dt,is_active\\n")
        except Exception as e:
            print(f"Warning: Failed to initialize adaptive bath controller output file: {e}")
    
    def _set_dynamic_target(self, current_time_ps: float):
        """Set dynamic target temperature based on current signal temperature at turn-on time."""
        if not self.dynamic_target or self.dynamic_target_set:
            return
        
        # Calculate signal temperature at turn-on time using the dynamic target method
        target_signal_temp = self._calculate_signal_temperature_with_method(self.dynamic_target_temperature_method)
        if target_signal_temp is not None and target_signal_temp > self.T_min:
            self.target_temperature = target_signal_temp
            self.dynamic_target_set = True
            print(f"Adaptive bath controller: Dynamic target set to {self.target_temperature:.2f} K")
            print(f"  Target method: {self.dynamic_target_temperature_method} at t = {current_time_ps:.2f} ps")
        else:
            # Fallback: try to get current bath temperature
            bath_temp = self._get_current_bath_temperature()
            if bath_temp is not None and bath_temp > self.T_min:
                self.target_temperature = bath_temp
                print(f"Adaptive bath controller: Dynamic target set to {self.target_temperature:.2f} K (bath temperature) at t = {current_time_ps:.2f} ps")
            else:
                print(f"Warning: Could not determine signal or bath temperature for dynamic target, using {self.target_temperature:.2f} K")
            self.dynamic_target_set = True
    
    def _calculate_signal_temperature(self):
        """Calculate control signal temperature using the specified method."""
        return self._calculate_signal_temperature_with_method(self.signal_temperature_method)
    
    def _calculate_signal_temperature_with_method(self, method: str):
        """Calculate signal temperature using a specific method.
        
        Parameters
        ----------
        method : str
            Temperature calculation method: 'kinetic', 'harmonic_equipartition', 
            'lj_coulombic', or 'harmonic'
        
        Returns
        -------
        float or None
            Calculated temperature in Kelvin, or None if calculation fails
        """
        # Only print DEBUG every 100 calls to reduce output spam
        if not hasattr(self, '_debug_counter'):
            self._debug_counter = 0
        self._debug_counter += 1
        
        debug_print = (self._debug_counter % 100 == 1)  # Print 1st, 101st, 201st, etc.
        
        if debug_print:
            print(f"DEBUG: AdaptiveBathController calculating signal temperature using method: {method}")
        
        if method == 'kinetic':
            temp = self._calculate_kinetic_temperature()
        elif method == 'harmonic_equipartition':
            temp = self._calculate_harmonic_equipartition_temperature()
        elif method == 'lj_coulombic':
            temp = self._calculate_lj_coulombic_temperature()
        elif method == 'harmonic':
            temp = self._calculate_harmonic_fictive_temperature()
        else:
            print(f"Warning: Unknown signal temperature method: {method}")
            return None
        
        if debug_print:
            print(f"DEBUG: Signal temperature calculated: {temp} K using method {method}")
        return temp
    
    def _calculate_kinetic_temperature(self):
        """Calculate kinetic temperature from particle velocities."""
        try:
            hoomd_sim = self._get_hoomd_simulation()
            if hoomd_sim and hasattr(hoomd_sim.state, 'thermodynamic_quantities'):
                return hoomd_sim.state.thermodynamic_quantities.kinetic_temperature
            return None
        except Exception as e:
            print(f"Warning: Failed to calculate kinetic temperature: {e}")
            return None
    
    def _calculate_harmonic_equipartition_temperature(self):
        """Calculate harmonic equipartition temperature directly."""
        try:
            if not self.energy_tracker:
                return None
            
            energy_data = self.energy_tracker.get_instantaneous_energy()
            harmonic_energy = energy_data.get('harmonic', 0.0)
            if harmonic_energy is None or harmonic_energy <= 0:
                return None
            
            # Get number of particles for proper equipartition
            with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                N_particles = len(snap.particles.mass)
                # Exclude cavity particle if it exists (typically the last one)
                if hasattr(self.simulation, 'cavity_particle_exists') and self.simulation.cavity_particle_exists:
                    N_particles -= 1
                elif N_particles > 1:  # Assume cavity exists if more than 1 particle
                    N_particles -= 1
            
            if N_particles <= 0:
                return None
            
            # Harmonic equipartition: T = 4*U/(N*kB) where N is number of atoms
            from .utils import PhysicalConstants
            kB_atomic = PhysicalConstants.KB_HARTREE_PER_K
            temperature = (4.0 * harmonic_energy) / (N_particles * kB_atomic)
            
            return max(temperature, 0.0)
            
        except Exception as e:
            print(f"Warning: Failed to calculate harmonic equipartition temperature: {e}")
            return None
    
    def _calculate_lj_coulombic_temperature(self):
        """Calculate LJ+Coulombic fictive temperature using empirical data."""
        try:
            if not self.empirical_data or not self.energy_tracker:
                print(f"DEBUG: Missing empirical_data={self.empirical_data is not None} or energy_tracker={self.energy_tracker is not None}")
                return None
            
            energy_data = self.energy_tracker.get_instantaneous_energy()
            
            # Get LJ and Coulombic energies separately and sum them
            lj_energy = energy_data.get('lj', 0.0)
            coulombic_energy = energy_data.get('coulombic', 0.0)
            
            if lj_energy is None or coulombic_energy is None:
                print(f"DEBUG: Missing energy components: lj={lj_energy}, coulombic={coulombic_energy}")
                return None
            
            lj_coul_energy = lj_energy + coulombic_energy
            print(f"DEBUG: LJ+Coulombic energy = {lj_energy} + {coulombic_energy} = {lj_coul_energy} hartree")
            
            temperature = self.empirical_data.calculate_systemic_temperature(lj_coul_energy)
            print(f"DEBUG: LJ+Coulombic temperature from empirical data: {temperature} K")
            
            return temperature
            
        except Exception as e:
            print(f"Warning: Failed to calculate LJ+Coulombic temperature: {e}")
            return None
    
    def _calculate_harmonic_fictive_temperature(self):
        """Calculate harmonic fictive temperature using empirical data."""
        try:
            if not self.empirical_data_harmonic or not self.energy_tracker:
                return None
            
            energy_data = self.energy_tracker.get_instantaneous_energy()
            harmonic_energy = energy_data.get('harmonic', 0.0)
            if harmonic_energy is None:
                return None
            
            return self.empirical_data_harmonic.calculate_systemic_temperature(harmonic_energy)
            
        except Exception as e:
            print(f"Warning: Failed to calculate harmonic fictive temperature: {e}")
            return None
    
    def _get_current_bath_temperature(self):
        """Get current bath temperature from thermostats using robust API detection."""
        try:
            from .utils import PhysicalConstants
            
            # First, try direct thermostat references (most reliable)
            if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
                try:
                    if hasattr(self.molecular_thermostat, 'kT'):
                        kT = self.molecular_thermostat.kT
                        # Handle hoomd.variant.Constant objects
                        if hasattr(kT, '__call__'):
                            kT_value = kT(0)  # Call with dummy timestep
                        else:
                            kT_value = float(kT)
                        return kT_value / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.molecular_thermostat, 'thermostat') and hasattr(self.molecular_thermostat.thermostat, 'kT'):
                        # Handle nested thermostats (e.g., Bussi within ConstantVolume)
                        kT = self.molecular_thermostat.thermostat.kT
                        if hasattr(kT, '__call__'):
                            kT_value = kT(0)
                        else:
                            kT_value = float(kT)
                        return kT_value / PhysicalConstants.KB_HARTREE_PER_K
                except Exception as e:
                    print(f"Warning: Failed to get molecular thermostat temperature: {e}")
            
            # Try cavity thermostat as backup
            if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
                try:
                    if hasattr(self.cavity_thermostat, 'kT'):
                        kT = self.cavity_thermostat.kT
                        if hasattr(kT, '__call__'):
                            kT_value = kT(0)
                        else:
                            kT_value = float(kT)
                        return kT_value / PhysicalConstants.KB_HARTREE_PER_K
                    elif hasattr(self.cavity_thermostat, 'thermostat') and hasattr(self.cavity_thermostat.thermostat, 'kT'):
                        kT = self.cavity_thermostat.thermostat.kT
                        if hasattr(kT, '__call__'):
                            kT_value = kT(0)
                        else:
                            kT_value = float(kT)
                        return kT_value / PhysicalConstants.KB_HARTREE_PER_K
                except Exception as e:
                    print(f"Warning: Failed to get cavity thermostat temperature: {e}")
            
            # Fallback: try via simulation integrator (original approach)
            if hasattr(self.simulation, 'operations') and hasattr(self.simulation.operations, 'integrator'):
                integrator = self.simulation.operations.integrator
                if hasattr(integrator, 'thermostats') and len(integrator.thermostats) > 0:
                    first_thermostat = integrator.thermostats[0]
                    if hasattr(first_thermostat, 'kT'):
                        kT = first_thermostat.kT
                        if hasattr(kT, '__call__'):
                            kT_value = kT(0)
                        else:
                            kT_value = float(kT)
                        return kT_value / PhysicalConstants.KB_HARTREE_PER_K
                        
            print("Warning: Could not access any thermostat for bath temperature")
            return None
        except Exception as e:
            print(f"Warning: Failed to get current bath temperature: {e}")
            return None
    
    def _update_thermostats(self, temperature: float):
        """Update thermostat temperatures using robust API detection."""
        try:
            from .utils import PhysicalConstants
            import hoomd
            target_kT = temperature * PhysicalConstants.KB_HARTREE_PER_K
            
            # Update molecular thermostat
            if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is not None:
                try:
                    if hasattr(self.molecular_thermostat, 'kT'):
                        # Direct kT attribute (Langevin thermostats)
                        self.molecular_thermostat.kT = hoomd.variant.Constant(target_kT)
                    elif hasattr(self.molecular_thermostat, 'thermostat'):
                        # Nested thermostat (ConstantVolume with Bussi/MTTK)
                        nested_thermo = self.molecular_thermostat.thermostat
                        if hasattr(nested_thermo, 'kT'):
                            nested_thermo.kT = hoomd.variant.Constant(target_kT)
                        else:
                            print(f"Warning: Cannot update nested molecular thermostat kT")
                    else:
                        print(f"Warning: Cannot find kT attribute in molecular thermostat")
                except Exception as e:
                    print(f"Warning: Failed to update molecular thermostat: {e}")
            
            # Update cavity thermostat  
            if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is not None:
                try:
                    if hasattr(self.cavity_thermostat, 'kT'):
                        # Direct kT attribute (Langevin thermostats)
                        self.cavity_thermostat.kT = hoomd.variant.Constant(target_kT)
                    elif hasattr(self.cavity_thermostat, 'thermostat'):
                        # Nested thermostat (ConstantVolume with Bussi/MTTK)
                        nested_thermo = self.cavity_thermostat.thermostat
                        if hasattr(nested_thermo, 'kT'):
                            nested_thermo.kT = hoomd.variant.Constant(target_kT)
                        else:
                            print(f"Warning: Cannot update nested cavity thermostat kT")
                    else:
                        print(f"Warning: Cannot find kT attribute in cavity thermostat")
                except Exception as e:
                    print(f"Warning: Failed to update cavity thermostat: {e}")
                
        except Exception as e:
            print(f"Warning: Failed to update thermostats: {e}")
            import traceback
            traceback.print_exc()
    
    def act(self, timestep):
        """Execute adaptive bath temperature control using GD framework at each timestep."""
        try:
            # Get current time
            current_time_ps = self.time_tracker.elapsed_time
            
            # Check if controller should be active
            if current_time_ps < self.turn_on_time_ps:
                self.is_active = False
                return
            
            if self.turn_off_time_ps is not None and current_time_ps >= self.turn_off_time_ps:
                self.is_active = False
                return
            
            # Initialize on first activation
            if not self.is_active:
                self._set_dynamic_target(current_time_ps)
                # Initialize current bath temperature from thermostats
                self.current_bath_temperature = self._get_current_bath_temperature()
                if self.current_bath_temperature is None:
                    self.current_bath_temperature = self.target_temperature
                self.is_active = True
                self.last_update_time = current_time_ps
                print(f"Adaptive bath controller activated at t = {current_time_ps:.2f} ps")
                print(f"Initial bath temperature: {self.current_bath_temperature:.2f} K")
                print(f"Time constant: {self.time_constant_ps:.3f} ps")
                return
            
            # Update at specified intervals
            if current_time_ps - self.last_update_time < self.update_interval_ps:
                return
            
            # Calculate dt for this update step
            dt_ps = current_time_ps - self.last_update_time
            
            # Calculate signal temperature
            signal_temp = self._calculate_signal_temperature()
            if signal_temp is None:
                return
            
            # Calculate temperature error and dynamic target Tstar
            temp_error = signal_temp - self.target_temperature
            
            if temp_error > 0:
                # Signal too hot: cool down
                Tstar = self.target_temperature - self.amplitude_scale * abs(temp_error)
            elif temp_error < 0:
                # Signal too cold: heat up
                Tstar = self.target_temperature + self.amplitude_scale * abs(temp_error)
            else:
                # Perfect match: use target temperature
                Tstar = self.target_temperature
            
            # Implement GD evolution: dTb(t)/dt = -(Tb(t) - Tstar(t)) / tau
            dTb_dt = -(self.current_bath_temperature - Tstar) / self.time_constant_ps
            
            # Update bath temperature: Tb(t+dt) = Tb(t) + dt * dTb/dt
            new_bath_temp = self.current_bath_temperature + dt_ps * dTb_dt
            
            # Apply temperature limits
            if self.T_max is not None:
                actual_bath_temp = max(self.T_min, min(self.T_max, new_bath_temp))
            else:
                actual_bath_temp = max(self.T_min, new_bath_temp)
            
            # Store the updated temperature for next iteration
            self.current_bath_temperature = actual_bath_temp
            
            # Update thermostats
            self._update_thermostats(actual_bath_temp)
            self.last_update_time = current_time_ps
            
            # Console output
            if current_time_ps - self.last_console_output_time >= self.console_output_period_ps:
                signal_method_display = self.signal_temperature_method.replace('_', ' ').title()
                print(f"Adaptive Bath (GD) | Target: {self.target_temperature:.1f}K | "
                      f"Signal: {signal_temp:.1f}K | Tstar: {Tstar:.1f}K | "
                      f"Bath: {actual_bath_temp:.1f}K | dTb/dt: {dTb_dt:+.2f}K/ps | "
                      f"τ: {self.time_constant_ps:.2f}ps")
                self.last_console_output_time = current_time_ps
            
            # Log data
            try:
                with open(self.output_file, 'a', encoding='utf-8') as f:
                    f.write(f"{current_time_ps:.6f},{signal_temp:.6f},"
                           f"{self.target_temperature:.6f},{Tstar:.6f},"
                           f"{actual_bath_temp:.6f},{temp_error:.6f},"
                           f"{dTb_dt:.6f},{int(self.is_active)}\\n")
            except Exception as e:
                print(f"Warning: Failed to write adaptive bath controller output: {e}")
        
        except Exception as e:
            print(f"Error in adaptive bath temperature controller: {e}")


class QuenchController(hoomd.custom.Action):
    """
    Instantaneous temperature quench controller for cavity MD simulations.
    
    This controller provides an instantaneous quench from an initial temperature
    to a target temperature. Once activated, it immediately sets both thermostats
    to the target temperature.
    
    Parameters
    ----------
    initial_temperature : float
        Initial temperature in Kelvin
    target_temperature : float
        Target temperature in Kelvin (quench destination)
    quench_time_ps : float
        Time at which to activate the quench in picoseconds
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    molecular_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Molecular thermostat to control
    cavity_thermostat : hoomd.md.methods.thermostats.Thermostat or None
        Cavity thermostat to control
    apply_to : str, optional
        Which thermostats to control ('molecular', 'cavity', 'both'). Default: 'both'
    output_file : str, optional
        Output file for quench control data. Default: 'quench_control.csv'
    console_output_period_ps : float, optional
        Console output period in picoseconds. Default: 1.0
        
    Examples
    --------
    **Instantaneous quench from 300K to 4K at t=50ps:**
    
    >>> from hoomd.cavitymd.analysis import QuenchController
    >>> 
    >>> quench_controller = QuenchController(
    ...     initial_temperature=300.0,
    ...     target_temperature=4.0,
    ...     quench_time_ps=50.0,
    ...     time_tracker=time_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     apply_to='both'
    ... )
    
    Notes
    -----
    - Provides instantaneous quench (no gradual temperature ramp)
    - Quench occurs exactly at the specified quench_time_ps
    - Once quenched, temperature remains constant at target value
    """
    
    def __init__(self, 
                 initial_temperature: float,
                 target_temperature: float,
                 quench_time_ps: float,
                 time_tracker,
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 apply_to: str = 'both',
                 output_file: str = 'quench_control.csv',
                 console_output_period_ps: float = 1.0):
        
        super().__init__()
        
        # Store configuration
        self.initial_temperature = float(initial_temperature)
        self.target_temperature = float(target_temperature)
        self.quench_time_ps = float(quench_time_ps)
        self.time_tracker = time_tracker
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.apply_to = apply_to
        self.output_file = output_file
        self.console_output_period_ps = float(console_output_period_ps)
        
        # State variables
        self.is_quenched = False
        self.last_output_time = -1.0
        self.last_console_output_time = 0.0
        
        # Validate configuration
        self._validate_configuration()
        
        # Initialize output file
        self._initialize_output_file()
        
        print(f"QuenchController initialized:")
        print(f"   Initial temperature: {self.initial_temperature:.1f} K")
        print(f"   Target temperature: {self.target_temperature:.1f} K")
        print(f"   ΔT = {self.target_temperature - self.initial_temperature:+.1f} K")
        print(f"   Quench time: {self.quench_time_ps:.1f} ps")
        print(f"   Applied to: {self.apply_to}")
        print(f"   Quench type: {'Cooling' if self.target_temperature < self.initial_temperature else 'Heating'}")
    
    def _validate_configuration(self):
        """Validate controller configuration."""
        if self.initial_temperature <= 0:
            raise ValueError("initial_temperature must be positive")
        
        if self.target_temperature <= 0:
            raise ValueError("target_temperature must be positive")
        
        if self.quench_time_ps < 0:
            raise ValueError("quench_time_ps must be non-negative")
        
        if self.apply_to not in ['molecular', 'cavity', 'both']:
            raise ValueError("apply_to must be 'molecular', 'cavity', or 'both'")
        
        if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is None:
            raise ValueError("molecular_thermostat cannot be None when apply_to includes 'molecular'")
        
        if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is None:
            raise ValueError("cavity_thermostat cannot be None when apply_to includes 'cavity'")
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        try:
            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write("time_ps,is_quenched,current_temp,target_temp,applied_to\n")
        except Exception as e:
            print(f"Warning: Failed to initialize quench controller output file: {e}")
    
    def _update_thermostat(self, thermostat, temperature_K: float) -> bool:
        """Update thermostat temperature with immediate velocity rescaling."""
        if thermostat is None:
            return False
        
        try:
            from .utils import PhysicalConstants
            import hoomd
            
            # Convert temperature to kT in atomic units
            kT_hartree = temperature_K * PhysicalConstants.KB_HARTREE_PER_K
            
            # Handle different thermostat types
            updated = False
            
            if hasattr(thermostat, 'kT'):
                # Direct kT attribute (Langevin thermostats)
                thermostat.kT = hoomd.variant.Constant(kT_hartree)
                # Set tau to small value for fast response (but not extreme)
                if hasattr(thermostat, 'default_gamma'):
                    # For Langevin: gamma = 1/tau, so large gamma = small tau
                    # Use tau = 0.01 ps instead of 1e-6 for numerical stability
                    thermostat.default_gamma = 100.0  # Equivalent to tau = 0.01 ps
                updated = True
                
            elif hasattr(thermostat, 'thermostat'):
                # Nested thermostat (ConstantVolume with Bussi/MTTK)
                nested_thermo = thermostat.thermostat
                if hasattr(nested_thermo, 'kT'):
                    nested_thermo.kT = hoomd.variant.Constant(kT_hartree)
                    # Set tau to small value for fast response (but not extreme)
                    if hasattr(nested_thermo, 'tau'):
                        # Use tau = 0.01 ps instead of 1e-6 for numerical stability
                        nested_thermo.tau = PhysicalConstants.ps_to_atomic_units(0.01)
                    updated = True
            
            if not updated:
                print(f"Warning: Could not update thermostat of type {type(thermostat)}")
                return False
            
            # Note: tau is set to 0.01 ps for fast (but stable) temperature response
            # This provides rapid equilibration without numerical instabilities
            
            return True
        except Exception as e:
            print(f"Warning: Failed to update thermostat temperature: {e}")
            return False
    
    def act(self, timestep):
        """Execute quench control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if it's time to quench
        if not self.is_quenched and current_time_ps >= self.quench_time_ps:
            # Perform instantaneous quench
            molecular_updated = False
            cavity_updated = False
            
            if self.apply_to in ['molecular', 'both']:
                molecular_updated = self._update_thermostat(self.molecular_thermostat, self.target_temperature)
            
            if self.apply_to in ['cavity', 'both']:
                cavity_updated = self._update_thermostat(self.cavity_thermostat, self.target_temperature)
            
            if molecular_updated or cavity_updated:
                self.is_quenched = True
                print("")
                print(f"*** INSTANTANEOUS QUENCH ACTIVATED at t={current_time_ps:.3f}ps ***")
                print(f"*** Temperature: {self.initial_temperature:.1f}K → {self.target_temperature:.1f}K ***")
                print(f"*** Applied to: {self.apply_to} thermostat(s) ***")
                print("")
        
        # Determine current temperature
        current_temp = self.target_temperature if self.is_quenched else self.initial_temperature
        
        # Console output (periodic)
        should_print_console = (current_time_ps - self.last_console_output_time) >= self.console_output_period_ps
        
        if should_print_console:
            status = "QUENCHED" if self.is_quenched else "PRE-QUENCH"
            print(f"Quench Controller | Status: {status} | T = {current_temp:.1f}K | t = {current_time_ps:.1f}ps")
            self.last_console_output_time = current_time_ps
        
        # Log data (every timestep)
        if current_time_ps != self.last_output_time:
            try:
                with open(self.output_file, 'a', encoding='utf-8') as f:
                    f.write(f"{current_time_ps:.6f},{self.is_quenched},{current_temp:.6f},{self.target_temperature:.6f},{self.apply_to}\n")
                self.last_output_time = current_time_ps
            except Exception as e:
                print(f"Warning: Failed to write quench controller output: {e}")


class AutoStopController(hoomd.custom.custom_action.Action):
    """
    Auto-stop controller for cavity coupling convergence detection.
    
    This controller monitors the convergence between LJ+Coulombic fictive temperature
    and kinetic temperature. When they converge within a specified tolerance over
    an averaging window, it signals to disable cavity coupling.
    
    Your specifications:
    - Turn off coupling when 10 ps average of LJ+Coulombic fictive temperature 
      is less than 1 K from kinetic temperature
    - 10 ps window and 1 K tolerance should be tunable via flags
    - By default, turned on when coupling is first activated
    - Print console output showing calculations regardless of convergence
    
    Parameters
    ----------
    temperature_tracker : TemperatureTracker
        Temperature tracker for accessing temperature data
    time_tracker : ElapsedTimeTracker
        Time tracker for simulation timing
    tolerance : float
        Convergence tolerance in Kelvin (default: 1.0)
    window_ps : float
        Averaging window in picoseconds (default: 10.0)
    coupling_start_time_ps : float
        Time when coupling starts in picoseconds (default: 0.0)
    """
    
    def __init__(self, temperature_tracker, time_tracker, tolerance=1.0, window_ps=10.0, coupling_start_time_ps=0.0, target_temperature=None):
        super().__init__()
        
        self.temperature_tracker = temperature_tracker
        self.time_tracker = time_tracker
        self.tolerance = tolerance
        self.window_ps = window_ps
        self.coupling_start_time_ps = coupling_start_time_ps
        self.target_temperature = target_temperature  # Need target for 3-way convergence check
        
        # State variables
        self.is_converged = False
        self.convergence_time_ps = None
        self.simulation = None
        
        # History for averaging over exactly 10 ps window (in ps, not timesteps!)
        self.kinetic_history = []
        self.fictive_history = []
        self.time_history = []
        
        # Console output control - print every 5 ps to show calculations
        self.last_console_output_time = 0.0
        self.console_output_interval = 5.0  # Show calculations every 5 ps
        
        print(f"AutoStopController initialized:")
        print(f"  Tolerance: {self.tolerance:.1f} K (tunable via --auto-stop-tol)")
        print(f"  Averaging window: {self.window_ps:.1f} ps (tunable via --auto-stop-window)")
        print(f"  Coupling start time: {self.coupling_start_time_ps:.1f} ps")
        print(f"  Target temperature: {self.target_temperature:.1f} K" if self.target_temperature else "  Target temperature: Not specified")
        print(f"  Convergence requires ALL THREE conditions:")
        print(f"    |T_fictive - T_kinetic| < {self.tolerance:.1f} K")
        print(f"    |T_fictive - T_target| < {self.tolerance:.1f} K") 
        print(f"    |T_kinetic - T_target| < {self.tolerance:.1f} K")
    
    def set_simulation(self, simulation):
        """Set reference to main simulation object for signaling."""
        self.simulation = simulation
    
    def update(self):
        """Update auto-stop controller and check for convergence."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Only start monitoring after coupling has started (as you specified)
        if current_time_ps < self.coupling_start_time_ps:
            return
            
        # Get current temperatures
        try:
            kinetic_temp = getattr(self.temperature_tracker, 'kinetic_temperature', None)
            fictive_temp = getattr(self.temperature_tracker, 'lj_coulombic_fictive_temperature', None)
            
            if kinetic_temp is None or fictive_temp is None:
                # Still print what we're trying to do
                if current_time_ps - self.last_console_output_time >= self.console_output_interval:
                    print(f"AutoStop | t={current_time_ps:.1f}ps | Waiting for temperature data | "
                          f"kinetic={'available' if kinetic_temp else 'missing'} | "
                          f"fictive={'available' if fictive_temp else 'missing'}")
                    self.last_console_output_time = current_time_ps
                return
            
            # Add to history
            self.kinetic_history.append(kinetic_temp)
            self.fictive_history.append(fictive_temp)
            self.time_history.append(current_time_ps)
            
            # Remove old data outside the window (exactly your 10 ps specification)
            cutoff_time = current_time_ps - self.window_ps
            while self.time_history and self.time_history[0] < cutoff_time:
                self.kinetic_history.pop(0)
                self.fictive_history.pop(0)
                self.time_history.pop(0)
            
            # Calculate averages if we have enough data
            if len(self.kinetic_history) > 0:
                kinetic_avg = sum(self.kinetic_history) / len(self.kinetic_history)
                fictive_avg = sum(self.fictive_history) / len(self.fictive_history)
                window_span = current_time_ps - self.time_history[0] if self.time_history else 0.0
                
                # Calculate ALL THREE convergence conditions as you specified
                diff_fictive_kinetic = abs(fictive_avg - kinetic_avg)
                diff_fictive_target = abs(fictive_avg - self.target_temperature) if self.target_temperature else float('inf')
                diff_kinetic_target = abs(kinetic_avg - self.target_temperature) if self.target_temperature else float('inf')
                
                # Check if ALL THREE conditions are satisfied
                condition1 = diff_fictive_kinetic < self.tolerance  # |T_fictive - T_kinetic| < tol
                condition2 = diff_fictive_target < self.tolerance   # |T_fictive - T_target| < tol
                condition3 = diff_kinetic_target < self.tolerance   # |T_kinetic - T_target| < tol
                all_converged = condition1 and condition2 and condition3
                
                # Console output ALWAYS (as you requested - regardless of convergence)
                if current_time_ps - self.last_console_output_time >= self.console_output_interval:
                    target_str = f"{self.target_temperature:.2f}K" if self.target_temperature else "N/A"
                    print(f"AutoStop | t={current_time_ps:.1f}ps | "
                          f"T_kinetic_avg={kinetic_avg:.2f}K | T_fictive_avg={fictive_avg:.2f}K | T_target={target_str} | "
                          f"diff_fict_kin={diff_fictive_kinetic:.2f}K | diff_fict_tgt={diff_fictive_target:.2f}K | diff_kin_tgt={diff_kinetic_target:.2f}K | "
                          f"tol={self.tolerance:.1f}K | window={window_span:.1f}ps ({len(self.kinetic_history)} samples) | "
                          f"converged={'YES' if all_converged else 'NO'} ({condition1}/{condition2}/{condition3})")
                    self.last_console_output_time = current_time_ps
                
                # Check for convergence (ALL THREE conditions must be satisfied)
                if all_converged:
                    if not self.is_converged:
                        # First time convergence detected
                        self.is_converged = True
                        self.convergence_time_ps = current_time_ps
                        print("")
                        print(f"*** AUTO-STOP TRIGGERED at t={current_time_ps:.1f}ps ***")
                        print(f"*** ALL THREE convergence conditions satisfied: ***")
                        print(f"***   |T_fictive - T_kinetic| = {diff_fictive_kinetic:.2f}K < {self.tolerance:.1f}K ***")
                        print(f"***   |T_fictive - T_target|  = {diff_fictive_target:.2f}K < {self.tolerance:.1f}K ***")
                        print(f"***   |T_kinetic - T_target|  = {diff_kinetic_target:.2f}K < {self.tolerance:.1f}K ***")
                        print(f"*** Temperatures: T_fictive={fictive_avg:.2f}K, T_kinetic={kinetic_avg:.2f}K, T_target={self.target_temperature:.2f}K ***")
                        print(f"*** Coupling will be turned off indefinitely ***")
                        print("")
                        
                        # Set the auto-stop signal only once when first converged
                        self._set_auto_stop_signal()
                    
        except Exception as e:
            print(f"Warning: AutoStopController update failed: {e}")
    
    def _set_auto_stop_signal(self):
        """Set the auto-stop signal to disable coupling indefinitely."""
        try:
            # Set signal on main simulation object (most reliable method)
            if self.simulation is not None:
                self.simulation.auto_stop_coupling_signal = True
                print(f"AutoStop: Set signal on simulation object")
                
        except Exception as e:
            print(f"Warning: Failed to set auto-stop signal: {e}")
    
    def should_stop_coupling(self):
        """Check if coupling should be stopped (calls update)."""
        self.update()
        return self.is_converged
    
    def act(self, timestep):
        """HOOMD custom action interface."""
        self.update()
