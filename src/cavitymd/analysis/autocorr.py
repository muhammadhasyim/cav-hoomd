
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

from ..utils import PhysicalConstants, unwrap_positions
from .timing import ElapsedTimeTracker

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
        
    def _get_hoomd_simulation(self):
        """Return the HOOMD simulation object."""
        return self.simulation
    
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
    
    def _get_hoomd_simulation(self):
        """Return the HOOMD simulation object."""
        return self.simulation
    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output data."""
        if self.last_output_time is None:
            return True
        return (current_time_ps - self.last_output_time) >= self.output_period_ps


