
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
                 output_period_ps: float = 0.1,
                 enable_csv_output: bool = False):
        
        super().__init__()
        
        self.time_tracker = time_tracker
        self.output_file = output_file
        self.max_correlation_time_ps = float(max_correlation_time_ps)
        self.correlation_output_interval_ps = float(correlation_output_interval_ps)
        self.exclude_cavity = bool(exclude_cavity)
        self.enable_response_measurement = bool(enable_response_measurement)
        self.output_period_ps = float(output_period_ps)
        self.enable_csv_output = enable_csv_output
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
        # Initialize output file (only if CSV enabled)
        if self.enable_csv_output:
            self._initialize_output_file()
        else:
            self.output_file = None
        
        print(f"DipoleMomentFDRTracker initialized:")
        if self.output_file:
            print(f"   Output file: {self.output_file}")
        else:
            print(f"   CSV output: DISABLED (HDF5 only)")
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
        from ..utils import PhysicalConstants
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
        
        # Log to file (only if CSV output enabled)
        if self.enable_csv_output:
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
            if len(self.autocorrelation) > 0 and self.enable_csv_output:
                self._save_autocorrelation_data(self.correlation_times, self.autocorrelation)
        self.last_update_time = current_time_ps




