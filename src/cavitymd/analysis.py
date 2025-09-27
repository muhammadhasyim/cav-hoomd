
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
        with self.simulation.state.cpu_local_snapshot as snap:
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
            with self.simulation.state.cpu_local_snapshot as snap:
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
    """Tracks various energy components during simulation."""
    
    def __init__(self, simulation, time_tracker, output_period_ps, output_prefix):
        super().__init__()
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.output_period_ps = output_period_ps
        self.output_prefix = output_prefix
        self.last_output_time = None
        
        # Cache for force energy components (updated each timestep)
        self._current_energies = {}
        self._force_mapping = {}  # Maps force names to their objects

        # Initialize output file
        self.output_file = f"{output_prefix}_energy.txt"
        with open(self.output_file, 'w') as f:
            f.write("# Energy tracking data\n")
            f.write("# time_ps\tkinetic_energy\tpotential_energy\ttotal_energy\n")

    def act(self, timestep):
        """Track energy at each timestep and cache force energies."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Update cached energies every timestep for get_instantaneous_energy access
        self._update_energy_cache()
        
        if self._should_output(current_time_ps):
            # Get thermodynamic quantities
            with self.simulation.state.cpu_local_snapshot as snap:
                velocities = np.array(snap.particles.velocity)
                masses = np.array(snap.particles.mass)
            
            # Calculate kinetic energy
            kinetic_energy = 0.5 * np.sum(masses[:, np.newaxis] * velocities**2)
            
            # Calculate total potential energy from cached force energies
            potential_energy = sum(self._current_energies.values())
            total_energy = kinetic_energy + potential_energy
            
            # Output to file
            with open(self.output_file, 'a') as f:
                f.write(f"{current_time_ps:.6f}\t{kinetic_energy:.8f}\t{potential_energy:.8f}\t{total_energy:.8f}\n")
            
            self.last_output_time = current_time_ps
    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output data."""
        if self.last_output_time is None:
            return True
        return (current_time_ps - self.last_output_time) >= self.output_period_ps
    
    def _update_energy_cache(self):
        """Update cached force energy components from simulation."""
        self._current_energies.clear()
        
        try:
            # Access forces from the integrator
            if hasattr(self.simulation, 'operations') and hasattr(self.simulation.operations, 'integrator'):
                forces = self.simulation.operations.integrator.forces
                
                for force in forces:
                    force_type = type(force).__name__.lower()
                    
                    try:
                        # Get energy from force if available
                        if hasattr(force, 'energy'):
                            energy_value = force.energy
                            if isinstance(energy_value, (int, float)):
                                self._current_energies[force_type] = energy_value
                            elif hasattr(energy_value, '__float__'):
                                self._current_energies[force_type] = float(energy_value)
                                
                        # Special handling for cavity force components
                        if force_type == 'cavityforce':
                            if hasattr(force, 'harmonic_energy'):
                                self._current_energies['cavity_harmonic'] = force.harmonic_energy
                            if hasattr(force, 'coupling_energy'):
                                self._current_energies['cavity_coupling'] = force.coupling_energy
                            if hasattr(force, 'dipole_self_energy'):
                                self._current_energies['cavity_dipole'] = force.dipole_self_energy
                                
                    except Exception as e:
                        # Silently skip forces that can't provide energy
                        continue
                        
        except Exception as e:
            # If we can't access forces, keep previous energy values
            pass
    
    def get_instantaneous_energy(self) -> Dict[str, float]:
        """
        Get current instantaneous energy components for empirical feedback.
        
        Returns
        -------
        Dict[str, float]
            Dictionary containing energy components:
            - 'lj': Lennard-Jones energy in Hartree
            - 'coulombic': Coulombic energy in Hartree (includes ewald_short and ewald_long)
            - 'harmonic': Harmonic/bond energy in Hartree
            - 'cavity': Total cavity energy in Hartree
            - 'total_potential': Sum of all potential energy components in Hartree
            
        Notes
        -----
        This method is called by EmpiricalTemperatureFeedback to access
        real-time energy data for systemic temperature calculations.
        """
        energy_data = {}
        
        # Map HOOMD force names to standard energy components
        lj_energy = 0.0
        coulombic_energy = 0.0
        harmonic_energy = 0.0
        cavity_energy = 0.0
        
        # Process cached energies
        for force_name, energy in self._current_energies.items():
            if 'lj' in force_name or 'lennard' in force_name:
                lj_energy += energy
            elif 'ewald' in force_name or 'coulomb' in force_name or 'pppm' in force_name:
                coulombic_energy += energy
            elif 'harmonic' in force_name or 'bond' in force_name:
                if 'cavity' not in force_name:  # Exclude cavity harmonic from molecular harmonic
                    harmonic_energy += energy
            elif 'cavity' in force_name:
                cavity_energy += energy
        
        # Store individual components
        energy_data['lj'] = lj_energy
        energy_data['coulombic'] = coulombic_energy
        energy_data['harmonic'] = harmonic_energy
        energy_data['cavity'] = cavity_energy
        
        # Combined components for empirical feedback
        energy_data['lj_coulombic'] = lj_energy + coulombic_energy
        energy_data['total_potential'] = lj_energy + coulombic_energy + harmonic_energy + cavity_energy
        
        return energy_data


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
    Handles empirical potential energy vs temperature data with Rosenfeld function fitting.
    
    This class loads empirical energy-temperature relationships from equilibrium
    simulations and provides interpolation to determine systemic temperatures
    from instantaneous potential energies using the Rosenfeld-Tarazona scaling model.
    
    Parameters
    ----------
    data_file_path : str
        Path to file containing temperature and potential energy data
    energy_component : str, optional
        Which energy component to use ('lj_coulombic', 'total_PE', etc.). Default: 'lj_coulombic'
    use_direct_harmonic : bool, optional
        If True and energy_component='harmonic', use direct calculation T = 4*E/(N*kB). Default: False
        
    Attributes
    ----------
    temperatures : np.ndarray
        Temperature values from empirical data (K)
    energies : np.ndarray  
        Energy values from empirical data (hartree)
    fit_params : dict
        Parameters from Rosenfeld fitting: {'E0': float, 'A': float, 'alpha': float}
    is_fitted : bool
        Whether Rosenfeld fitting has been performed
    """
    
    def __init__(self, data_file_path: str, energy_component: str = 'lj_coulombic', use_direct_harmonic: bool = False):
        self.data_file_path = Path(data_file_path)
        self.energy_component = energy_component
        self.use_direct_harmonic = use_direct_harmonic
        self.fit_params = {}
        self.is_fitted = False
        
        self.load_empirical_data()
        if not self.use_direct_harmonic or self.energy_component != 'harmonic':
                self.fit_rosenfeld_function()
    
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
    
    def fit_rosenfeld_function(self):
        """Fit generalized Rosenfeld-Tarazona scaling function to the data."""
        try:
            from scipy.optimize import curve_fit
            
            def generalized_rt_model(T, E0, A, alpha):
                """Generalized RT scaling: E(T) = E₀ + A*T^α"""
                return E0 + A * T**alpha
            
            # Initial parameter guesses
            E0_guess = self.energies.min()  # Energy at 0K
            A_guess = (self.energies.max() - self.energies.min()) / (self.temperatures.max()**0.6)
            alpha_guess = 0.6  # Starting point between 3/5 and 3/2
            
            # Perform fitting
            self.fit_params, _ = curve_fit(
                generalized_rt_model,
                self.temperatures,
                self.energies,
                p0=[E0_guess, A_guess, alpha_guess],
                bounds=([self.energies.min() - 0.1, 0, 0.1], 
                       [self.energies.max(), np.inf, 2.0])
            )
            
            # Store fitted parameters
            self.fit_params = {
                'E0': self.fit_params[0],
                'A': self.fit_params[1], 
                'alpha': self.fit_params[2]
            }
            
            self.is_fitted = True
            
            print(f"Rosenfeld fitting completed:")
            print(f"  E₀ = {self.fit_params['E0']:.6f} Ha")
            print(f"  A = {self.fit_params['A']:.6f} Ha·K^(-α)")  
            print(f"  α = {self.fit_params['alpha']:.3f}")
            
        except Exception as e:
            print(f"Warning: Rosenfeld fitting failed: {e}")
            # Fallback to linear interpolation
            self.is_fitted = False
    
    def calculate_systemic_temperature(self, instantaneous_energy_hartree: float, num_particles: int = None) -> float:
        """
        Calculate systemic temperature from instantaneous potential energy.
        
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
        
        elif self.is_fitted:
                # Use fitted Rosenfeld function: E = E₀ + A*T^α
                # Solve for T: T = ((E - E₀)/A)^(1/α)
                E0, A, alpha = self.fit_params['E0'], self.fit_params['A'], self.fit_params['alpha']
                
                if instantaneous_energy_hartree <= E0:
                    return 0.0  # At or below 0K energy
                
                try:
                    temperature = ((instantaneous_energy_hartree - E0) / A) ** (1.0 / alpha)
                    return max(temperature, 0.0)
                except (ValueError, ZeroDivisionError):
                    return 0.0
                
        else:
            # Fallback to linear interpolation
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
                use_direct_harmonic=False  # Use fitted function (Rosenfeld model)
            )
            print(f"📊 Created harmonic empirical data for fictive temperature calculation")
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
        
        print("🌡️  Empirical temperature feedback controller initialized:")
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
                print(f"🌡️  Updated thermostat temperatures to {target_temperature_K:.2f} K")
            
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
                 empirical_data_file=None):
        
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.output_period_ps = output_period_ps
        self.output_file = output_file
        self.energy_tracker = energy_tracker
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        
        # Load empirical data for fictive temperature calculations if provided
        self.empirical_data_harmonic = None
        self.empirical_data_lj_coul = None
        if empirical_data_file is not None:
            try:
                # Load empirical data for harmonic energy component
                self.empirical_data_harmonic = EmpiricalTemperatureData(
                    data_file_path=empirical_data_file,
                    energy_component='harmonic',
                    use_direct_harmonic=False  # Use fitted function, not direct calculation
                )
                
                # Load empirical data for LJ+Coulombic energy component (uses Rosenfeld scaling)
                self.empirical_data_lj_coul = EmpiricalTemperatureData(
                    data_file_path=empirical_data_file,
                    energy_component='lj_coulombic'
                    # Rosenfeld scaling is the default - no additional parameters needed
                )
            except Exception as e:
                print(f"Warning: Could not load empirical data file {empirical_data_file}: {e}")
                print("Empirical fictive temperatures will not be available")
        
        # Tracking state
        self.last_output_time = None
        
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
        
        if self._should_output(current_time_ps):
            # 1. Calculate kinetic temperature
            kinetic_temp_K = self._calculate_kinetic_temperature()
            
            # 2. Calculate harmonic fictive temperature
            harmonic_fictive_K = self._calculate_harmonic_fictive_temperature()
            
            # 3. Calculate LJ+Coulombic fictive temperature
            lj_coul_fictive_K = self._calculate_lj_coul_fictive_temperature()
            
            # 4. Get cavity bath temperature
            cavity_bath_K = self._get_cavity_bath_temperature()
            
            # 5. Get molecular bath temperature
            molecular_bath_K = self._get_molecular_bath_temperature()
            
            # 6. Calculate harmonic equipartition temperature
            harmonic_equipartition_K = self._calculate_harmonic_equipartition_temperature()
            
            # Write to CSV
            with open(self.output_file, 'a') as f:
                f.write(f"{current_time_ps:.6f},{kinetic_temp_K:.6f},{harmonic_fictive_K:.6f},"
                       f"{lj_coul_fictive_K:.6f},{cavity_bath_K:.6f},{molecular_bath_K:.6f},"
                       f"{harmonic_equipartition_K:.6f}\n")
            
            self.last_output_time = current_time_ps
    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output data."""
        if self.last_output_time is None:
            return True
        return (current_time_ps - self.last_output_time) >= self.output_period_ps
    
    def _calculate_kinetic_temperature(self) -> float:
        """Calculate kinetic temperature from particle velocities."""
        try:
            from .utils import PhysicalConstants
            
            with self.simulation.state.cpu_local_snapshot as snap:
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
            
        except Exception as e:
            print(f"Warning: Could not calculate kinetic temperature: {e}")
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
            with self.simulation.state.cpu_local_snapshot as snap:
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
                return 0.0
            
            # Get harmonic energy from energy tracker
            energy_dict = self.energy_tracker.get_instantaneous_energy()
            harmonic_energy = energy_dict.get('harmonic', 0.0)
            
            if harmonic_energy <= 0:
                return 0.0
            
            # Get number of particles
            from .utils import PhysicalConstants
            with self.simulation.state.cpu_local_snapshot as snap:
                N_particles = len(snap.particles.mass)
            
            if N_particles == 0:
                return 0.0
            
            # Direct equipartition calculation: T = 4*E/(N*kB) for 3D harmonic oscillator
            # For a classical 3D harmonic oscillator: <E> = (3/2)*N*kB*T for kinetic + (3/2)*N*kB*T for potential
            # But here we only have the harmonic potential energy, so: <E_harmonic> = (3/2)*N*kB*T
            # Therefore: T = (2*E_harmonic)/(3*N*kB)
            # However, based on the empirical observations, the factor is 4 instead of 2/3
            kB_hartree = PhysicalConstants.KB_HARTREE_PER_K
            temperature_K = (4.0 * harmonic_energy) / (N_particles * kB_hartree)
            
            return temperature_K
            
        except Exception as e:
            print(f"Warning: Could not calculate harmonic equipartition temperature: {e}")
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
        Temperature calculation method ('kinetic', 'lj_coulombic', 'harmonic')
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
        if empirical_data_file and self.temperature_method in ['lj_coulombic', 'harmonic']:
            try:
                if self.temperature_method == 'harmonic':
                    self.empirical_data = EmpiricalTemperatureData(
                        empirical_data_file, energy_component='harmonic'
                    )
                else:  # lj_coulombic
                    self.empirical_data = EmpiricalTemperatureData(
                        empirical_data_file, energy_component='lj_coulombic'
                    )
                print(f"Loaded empirical data for {self.temperature_method} temperature calculation")
            except Exception as e:
                print(f"Warning: Failed to load empirical data: {e}")
        
        # Validate configuration
        self._validate_configuration()
        
        # Initialize output file
        self._initialize_output_file()
        
        # Print configuration
        print(f"🎯 Gradient Descent Temperature Controller initialized:")
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
        if self.temperature_method not in ['kinetic', 'lj_coulombic', 'harmonic']:
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
                
                with self.simulation.state.cpu_local_snapshot as snap:
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
                return self.empirical_data.calculate_systemic_temperature(total_energy)
                
            elif self.temperature_method == 'harmonic':
                if self.empirical_data is None:
                    return None
                
                # Get harmonic energy
                energy_dict = self.energy_tracker.get_instantaneous_energy()
                harmonic_energy = energy_dict.get('harmonic', 0.0)
                
                # Convert to temperature using empirical data
                return self.empirical_data.calculate_systemic_temperature(harmonic_energy)
            
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
                print(f"🎯 GD Controller | Target: {self.target_temperature:.1f}K | Measured: {measured_temperature:.1f}K | "
                      f"Effective: {effective_temperature:.1f}K | Bath→{bath_temperature:.1f}K")
            else:
                print(f"🎯 GD Controller | Target: {self.target_temperature:.1f}K | Bath→{bath_temperature:.1f}K")
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
                print(f"🎯 Gradient descent feedback turned ON at t = {current_time_ps:.2f} ps")
                print(f"🎯 Initial bath temperature: {self.current_bath_temperature:.2f} K")
            else:
                print(f"🎯 Gradient descent feedback turned ON at t = {current_time_ps:.2f} ps")
        elif not should_be_active and self.is_active:
            self.is_active = False
            print(f"🎯 Gradient descent feedback turned OFF at t = {current_time_ps:.2f} ps")
        
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


class PITemperatureFeedback(hoomd.custom.Action):
    """
    PI feedback controller for cavity MD simulations with IMC/λ auto-tuning.
    
    This class implements a true Proportional-Integral (PI) controller for managing
    thermostat temperatures based on system temperature measurements. It includes:
    
    - IMC/λ tuning methodology with molecular time constant
    - Anti-windup protection with back-calculation
    - Rate limiting for smooth temperature changes
    - Support for kinetic, LJ+Coulombic, and harmonic temperature methods
    - Configurable setpoint weighting (β parameter)
    
    The PI controller follows standard process control practices:
    
    .. math::
        
        u(t) = K_c [β r(t) - y(t)] + I(t)
        
        \\frac{dI}{dt} = \\frac{K_c}{T_i} e(t) + \\frac{u_{sat}(t) - u_{raw}(t)}{T_{aw}}
    
    where:
    - u(t) is the controller output (bath temperature)
    - r(t) is the setpoint (target temperature)  
    - y(t) is the measured system temperature
    - e(t) = r(t) - y(t) is the error
    - β is the setpoint weighting factor (default: 0.7)
    - T_aw is the anti-windup time constant
    
    **IMC/λ Auto-Tuning:**
    
    Based on the molecular thermal time constant τ, the controller parameters are:
    
    .. math::
        
        K_c = \\frac{τ}{λ + τ}, \\quad T_i = 2τ, \\quad λ = 2τ \\text{ to } 4τ
    
    This provides robust, well-damped response for cavity MD systems.
    
    Parameters
    ----------
    temperature_method : str
        Temperature calculation method ('kinetic', 'lj_coulombic', 'harmonic')
    molecular_tau_ps : float
        Molecular thermostat time constant in picoseconds (for auto-tuning)
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
    lambda_factor : float, optional
        Closed-loop time scale factor (λ = lambda_factor * τ). Default: 3.0
    beta : float, optional
        Setpoint weighting factor (0.0 to 1.0). Default: 0.7
    apply_to : str, optional
        Which thermostats to control ('molecular', 'cavity', 'both'). Default: 'both'
    turn_on_time_ps : float, optional
        Time to start PI control in picoseconds. Default: 0.0
    turn_off_time_ps : float, optional
        Time to stop PI control in picoseconds. Default: None (never turn off)
    update_interval_ps : float, optional
        Interval between control updates in picoseconds. Default: 5.0
    T_min : float, optional
        Minimum allowed bath temperature in Kelvin. Default: 0.0
    T_max : float, optional
        Maximum allowed bath temperature in Kelvin. Default: None (no upper limit)
    rate_limit_K_per_ps : float, optional
        Maximum rate of temperature change in K/ps. Default: None (no rate limit)
    output_file : str, optional
        Output file for PI control data. Default: 'pi_feedback.csv'
    empirical_data_file : str, optional
        Path to empirical data file (required for 'lj_coulombic' and 'harmonic' methods)
        
    Attributes
    ----------
    Kc : float
        Proportional gain (auto-tuned)
    Ti : float  
        Integral time constant in picoseconds (auto-tuned)
    integral_state : float
        Current integral state
    last_output : float
        Last controller output (for rate limiting)
    is_active : bool
        Whether PI control is currently active
        
    Examples
    --------
    **Basic kinetic temperature PI control:**
    
    >>> from hoomd.cavitymd.analysis import PITemperatureFeedback
    >>> 
    >>> pi_controller = PITemperatureFeedback(
    ...     temperature_method='kinetic',
    ...     molecular_tau_ps=5.0,  # 5 ps molecular tau
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     target_temperature=100.0
    ... )
    
    **LJ+Coulombic fictive temperature control:**
    
    >>> pi_controller = PITemperatureFeedback(
    ...     temperature_method='lj_coulombic',
    ...     molecular_tau_ps=5.0,
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     cavity_thermostat=cavity_thermo,
    ...     target_temperature=100.0,
    ...     empirical_data_file='equilibrium_data.txt',
    ...     apply_to='both',
    ...     lambda_factor=3.0  # Conservative tuning
    ... )
    
    **Harmonic fictive temperature with custom tuning:**
    
    >>> pi_controller = PITemperatureFeedback(
    ...     temperature_method='harmonic',
    ...     molecular_tau_ps=2.0,  # Fast molecular response
    ...     time_tracker=time_tracker,
    ...     energy_tracker=energy_tracker,
    ...     molecular_thermostat=molecular_thermo,
    ...     target_temperature=80.0,
    ...     empirical_data_file='harmonic_calibration.txt',
    ...     lambda_factor=2.0,  # More aggressive tuning
    ...     beta=0.5,  # Reduced setpoint weighting
    ...     rate_limit_K_per_ps=0.1  # Smooth temperature changes
    ... )
    
    Notes
    -----
    - Requires accurate molecular τ for proper auto-tuning
    - λ factor controls aggressiveness: 2τ (fast), 3τ (balanced), 4τ (conservative)
    - Anti-windup prevents integral saturation during temperature limits
    - Rate limiting smooths temperature changes to prevent simulation instability
    - Compatible with all existing temperature calculation methods
    """
    
    def __init__(self,
                 temperature_method: str,
                 molecular_tau_ps: float,
                 time_tracker,  # ElapsedTimeTracker
                 energy_tracker,  # EnergyTracker
                 molecular_thermostat=None,
                 cavity_thermostat=None,
                 target_temperature: float = 100.0,
                 lambda_factor: float = 3.0,
                 beta: float = 0.7,
                 apply_to: str = 'both',
                 turn_on_time_ps: float = 0.0,
                 turn_off_time_ps: Optional[float] = None,
                 update_interval_ps: float = 5.0,
                 T_min: float = 0.0,
                 T_max: Optional[float] = None,
                 rate_limit_K_per_ps: Optional[float] = None,
                 output_file: str = 'pi_feedback.csv',
                 empirical_data_file: Optional[str] = None):
        
        super().__init__()
        
        # Validate inputs
        if temperature_method not in ['kinetic', 'lj_coulombic', 'harmonic']:
            raise ValueError(f"temperature_method must be 'kinetic', 'lj_coulombic', or 'harmonic', got {temperature_method}")
        
        if temperature_method in ['lj_coulombic', 'harmonic'] and empirical_data_file is None:
            raise ValueError(f"empirical_data_file is required for temperature_method '{temperature_method}'")
        
        if apply_to not in ['molecular', 'cavity', 'both']:
            raise ValueError(f"apply_to must be 'molecular', 'cavity', or 'both', got {apply_to}")
        
        if not 0.0 <= beta <= 1.0:
            raise ValueError(f"beta must be between 0.0 and 1.0, got {beta}")
        
        # Store core parameters
        self.temperature_method = temperature_method
        self.molecular_tau_ps = float(molecular_tau_ps)
        self.time_tracker = time_tracker
        self.energy_tracker = energy_tracker
        self.molecular_thermostat = molecular_thermostat
        self.cavity_thermostat = cavity_thermostat
        self.target_temperature = float(target_temperature)
        self.lambda_factor = float(lambda_factor)
        self.beta = float(beta)
        self.apply_to = apply_to
        self.turn_on_time_ps = float(turn_on_time_ps)
        self.turn_off_time_ps = float(turn_off_time_ps) if turn_off_time_ps is not None else None
        self.update_interval_ps = float(update_interval_ps)
        self.T_min = float(T_min)
        self.T_max = float(T_max) if T_max is not None else None
        self.rate_limit_K_per_ps = float(rate_limit_K_per_ps) if rate_limit_K_per_ps is not None else None
        self.output_file = output_file
        
        # Set up empirical data if needed
        self.empirical_data = None
        if empirical_data_file is not None:
            self.empirical_data = EmpiricalTemperatureData(
                data_file_path=empirical_data_file,
                energy_component=temperature_method
            )
        
        # IMC/λ auto-tuning calculation
        self.lambda_ps = self.lambda_factor * self.molecular_tau_ps
        self.Kc = self.molecular_tau_ps / (self.lambda_ps + self.molecular_tau_ps)
        self.Ti = 2.0 * self.molecular_tau_ps  # Integral time constant
        self.Taw = self.Ti  # Anti-windup time constant (typically Ti)
        
        # PI controller state
        self.integral_state: float = 0.0
        self.last_output: float = self.target_temperature
        self.last_update_time: float = 0.0
        self.is_active: bool = False
        
        # Output tracking
        self.measurement_history = []
        self.output_history = []
        self.time_history = []
        
        # Validation
        if self.apply_to in ['molecular', 'both'] and self.molecular_thermostat is None:
            raise ValueError("molecular_thermostat cannot be None when apply_to includes 'molecular'")
        
        if self.apply_to in ['cavity', 'both'] and self.cavity_thermostat is None:
            raise ValueError("cavity_thermostat cannot be None when apply_to includes 'cavity'")
        
        # Initialize output file
        self._initialize_output_file()
        
        print("🎛️  PI Temperature Feedback Controller initialized:")
        print(f"   Temperature method: {self.temperature_method}")
        print(f"   Molecular τ: {self.molecular_tau_ps:.1f} ps")
        print(f"   IMC/λ Parameters:")
        print(f"     λ factor: {self.lambda_factor:.1f}")
        print(f"     λ time: {self.lambda_ps:.1f} ps")
        print(f"     Kc (proportional gain): {self.Kc:.3f}")
        print(f"     Ti (integral time): {self.Ti:.1f} ps")
        print(f"   Setpoint weighting β: {self.beta:.2f}")
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
    
    def _initialize_output_file(self):
        """Initialize CSV output file with headers."""
        with open(self.output_file, 'w') as f:
            f.write("# PI Temperature Feedback Control Data\n")
            f.write("# time_ps: Simulation time in picoseconds\n")
            f.write("# measured_T_K: Measured system temperature in Kelvin\n")
            f.write("# target_T_K: Target (setpoint) temperature in Kelvin\n")
            f.write("# error_K: Temperature error (target - measured) in Kelvin\n")
            f.write("# integral_state: Current integral state\n")
            f.write("# controller_output_K: Raw controller output in Kelvin\n")
            f.write("# saturated_output_K: Saturated output (within limits) in Kelvin\n")
            f.write("# rate_limited_output_K: Final output after rate limiting in Kelvin\n")
            f.write("# is_active: Whether PI control is active (0/1)\n")
            f.write("time_ps,measured_T_K,target_T_K,error_K,integral_state,controller_output_K,saturated_output_K,rate_limited_output_K,is_active\n")
    
    def _calculate_system_temperature(self, current_time_ps: float) -> Optional[float]:
        """Calculate system temperature using the specified method."""
        try:
            if self.temperature_method == 'kinetic':
                # Calculate kinetic temperature from energy tracker
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if 'molecular_kinetic' in energy_data:
                    # T = 2 * KE / (3 * N * kB) for 3D system
                    kinetic_energy = energy_data['molecular_kinetic']
                    # Assuming N particles - need to get from energy tracker or estimate
                    # For now, use simplified approach
                    from .utils import PhysicalConstants
                    kB = PhysicalConstants.KB_HARTREE_PER_K
                    # Estimate: typical 500 particles, 3 DOF each
                    estimated_dof = 1500  # This should ideally come from system info
                    temperature = 2.0 * kinetic_energy / (estimated_dof * kB)
                    return temperature
                else:
                    print(f"Warning: Kinetic energy not available for temperature calculation")
                    return None
                    
            elif self.temperature_method == 'lj_coulombic':
                # Calculate fictive temperature using empirical LJ+Coulombic data
                if self.empirical_data is None:
                    print(f"Warning: No empirical data available for LJ+Coulombic temperature")
                    return None
                
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if 'lj' in energy_data and 'coulombic' in energy_data:
                    lj_coulombic_energy = energy_data['lj'] + energy_data['coulombic']
                    temperature = self.empirical_data.calculate_systemic_temperature(lj_coulombic_energy)
                    return temperature
                else:
                    print(f"Warning: LJ+Coulombic energy not available for temperature calculation")
                    return None
                    
            elif self.temperature_method == 'harmonic':
                # Calculate fictive temperature using harmonic energy with T^(3/2) relation
                if self.empirical_data is None:
                    print(f"Warning: No empirical data available for harmonic temperature")
                    return None
                    
                energy_data = self.energy_tracker.get_instantaneous_energy()
                if 'harmonic' in energy_data:
                    harmonic_energy = energy_data['harmonic']
                    temperature = self.empirical_data.calculate_systemic_temperature(harmonic_energy)
                    return temperature
                else:
                    print(f"Warning: Harmonic energy not available for temperature calculation")
                    return None
                    
            else:
                print(f"Warning: Unknown temperature method '{self.temperature_method}'")
                return None
                
        except Exception as e:
            print(f"Warning: Failed to calculate {self.temperature_method} temperature: {e}")
            return None
    
    def _apply_rate_limit(self, new_output: float, dt_ps: float) -> float:
        """Apply rate limiting to controller output."""
        if self.rate_limit_K_per_ps is None:
            return new_output
        
        max_change = self.rate_limit_K_per_ps * dt_ps
        change = new_output - self.last_output
        
        if abs(change) <= max_change:
            return new_output
        else:
            # Limit the change
            limited_change = max_change if change > 0 else -max_change
            return self.last_output + limited_change
    
    def _update_thermostats(self, bath_temperature: float):
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
        
        if updated_any:
            print(f"🌡️  Updated thermostat temperatures to {bath_temperature:.2f} K (kT = {target_kT:.6f} a.u.)")
    
    def act(self, timestep):
        """Execute PI control at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Check if we should start
        if current_time_ps < self.turn_on_time_ps:
            return
        
        # Check if we should stop
        if self.turn_off_time_ps is not None and current_time_ps >= self.turn_off_time_ps:
            if self.is_active:
                print(f"🔴 PI feedback turned OFF at t = {current_time_ps:.2f} ps")
                self.is_active = False
            return
        
        # Check if it's time to update
        if current_time_ps - self.last_update_time < self.update_interval_ps:
            return
        
        if not self.is_active:
            print(f"🟢 PI feedback turned ON at t = {current_time_ps:.2f} ps")
            self.is_active = True
        
        # Calculate system temperature
        measured_temperature = self._calculate_system_temperature(current_time_ps)
        if measured_temperature is None:
            print(f"Warning: Could not measure temperature, skipping PI update")
            return
        
        # PI control calculations
        dt_ps = current_time_ps - self.last_update_time if self.last_update_time > 0 else self.update_interval_ps
        
        # Error calculation
        error = self.target_temperature - measured_temperature
        
        # Proportional term with setpoint weighting
        proportional_term = self.Kc * (self.beta * self.target_temperature - measured_temperature)
        
        # Integral term with anti-windup
        controller_output_raw = proportional_term + self.integral_state
        
        # Apply saturation limits
        if self.T_max is not None:
            controller_output_sat = max(self.T_min, min(self.T_max, controller_output_raw))
        else:
            controller_output_sat = max(self.T_min, controller_output_raw)
        
        # Anti-windup: back-calculate integral correction
        saturation_error = controller_output_sat - controller_output_raw
        integral_correction = saturation_error / self.Taw
        
        # Update integral state
        integral_update = (self.Kc / self.Ti) * error + integral_correction
        self.integral_state += integral_update * dt_ps
        
        # Apply rate limiting
        final_output = self._apply_rate_limit(controller_output_sat, dt_ps)
        
        # Update thermostats
        self._update_thermostats(final_output)
        
        # Store data for output
        self.measurement_history.append(measured_temperature)
        self.output_history.append(final_output)
        self.time_history.append(current_time_ps)
        
        # Write to CSV file
        with open(self.output_file, 'a') as f:
            f.write(f"{current_time_ps:.6f},{measured_temperature:.3f},{self.target_temperature:.3f},"
                   f"{error:.3f},{self.integral_state:.6f},{controller_output_raw:.3f},"
                   f"{controller_output_sat:.3f},{final_output:.3f},{int(self.is_active)}\n")
        
        # Update state
        self.last_output = final_output


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
                 enable_response_measurement: bool = False):
        
        super().__init__()
        
        self.time_tracker = time_tracker
        self.output_file = output_file
        self.max_correlation_time_ps = float(max_correlation_time_ps)
        self.correlation_output_interval_ps = float(correlation_output_interval_ps)
        self.exclude_cavity = bool(exclude_cavity)
        self.enable_response_measurement = bool(enable_response_measurement)
        
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
        
        # Initialize output file
        self._initialize_output_file()
        
        print(f"📊 DipoleMomentFDRTracker initialized:")
        print(f"   Output file: {self.output_file}")
        print(f"   Max correlation time: {self.max_correlation_time_ps:.1f} ps")
        print(f"   Field direction: {self.field_direction}")
        print(f"   Exclude cavity: {self.exclude_cavity}")
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
                    
            print(f"📊 Autocorrelation data saved to {correlation_file}")
            
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
