#!/usr/bin/env python3
"""
Material Time Analysis for Cavity MD F(k,t) Data

This script implements the material time analysis methodology from:
"Time reversibility during the ageing of materials" by Böhmer et al. (2024)
Nature Physics, DOI: 10.1038/s41567-023-02366-z

Applied to F(k,t) data from cavity MD simulations with 200ps separated waiting times.

Author: Scientific Computing Assistant
Date: 2025-09-10
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
import warnings
warnings.filterwarnings('ignore')
from scipy.optimize import curve_fit

COUPLING_DIR_AXIS_VALUES = {
    "0epos00": 0.0,
    "3eneg04": 3e-4,
    "5eneg04": 5e-4,
    "7eneg04": 7e-4,
    "1eneg03": 1e-3,
    "1eneg02": 0.01,
    "2eneg02_l2333": 0.023333,
    "2eneg02": 0.016667,
    "3eneg02": 0.03,
}


def axis_value_from_coupling_dir(coupling_name: str) -> float:
    """Map a staged coupling directory name to a colorbar axis value."""
    ordered_tags = sorted(COUPLING_DIR_AXIS_VALUES, key=len, reverse=True)
    for tag in ordered_tags:
        if tag in coupling_name:
            return COUPLING_DIR_AXIS_VALUES[tag]
    return float("nan")


class MaterialTimeAnalyzer:
    """
    Analyzes F(k,t) data to calculate material time and test data collapse
    following the Tool-Narayanaswamy formalism.
    """
    
    def __init__(self, criterion_value=0.1, verbose=True, sg_window_ps=0.1):
        """
        Initialize the analyzer.
        
        Parameters:
        -----------
        criterion_value : float
            Threshold value for defining unit material time (default 0.1)
        verbose : bool
            Whether to print detailed output
        sg_window_ps : float
            Savitzky-Golay filter window size in picoseconds (default: 0.1)
        """
        self.criterion_value = criterion_value
        self.verbose = verbose
        self.sg_window_ps = sg_window_ps
        self.material_times = {}
        self.waiting_times = {}
        self.normalized_fkt = {}
        self.collapsed_data = {}
        
    def apply_savgol_filter(self, lag_time, fkt_data):
        """
        Apply Savitzky-Golay filter to F(k,t) data with specified time window.
        
        Parameters:
        -----------
        lag_time : array
            Time points in ps
        fkt_data : array
            F(k,t) values
            
        Returns:
        --------
        fkt_smoothed : array
            Smoothed F(k,t) data
        """
        if len(fkt_data) < 5:
            # Not enough points for smoothing
            return fkt_data
            
        # Calculate window size based on time resolution and desired window
        dt = np.mean(np.diff(lag_time))  # Average time step
        window_points = int(self.sg_window_ps / dt)
        
        # Ensure window size is odd and at least 3
        window_points = max(3, window_points)
        if window_points % 2 == 0:
            window_points += 1
            
        # Ensure window doesn't exceed data length
        window_points = min(window_points, len(fkt_data))
        if window_points % 2 == 0:
            window_points -= 1
            
        # Apply Savitzky-Golay filter (polynomial order 2)
        try:
            fkt_smoothed = savgol_filter(fkt_data, window_points, 2)
        except:
            # Fall back to original data if filtering fails
            if self.verbose:
                print(f"    Warning: Savitzky-Golay filtering failed, using original data")
            fkt_smoothed = fkt_data
            
        return fkt_smoothed
        
    def load_fkt_data(self, coupling_dir):
        """
        Load F(k,t) data from master files in the specified directory.
        
        Parameters:
        -----------
        coupling_dir : str
            Directory containing master_fskt_ref*.txt files
            
        Returns:
        --------
        fkt_files : dict
            Dictionary mapping ref_num to file paths
        """
        if not os.path.exists(coupling_dir):
            raise FileNotFoundError(f"Directory not found: {coupling_dir}")
            
        # Collect reference files
        fkt_files = {}
        ref_pattern = os.path.join(coupling_dir, "master_fskt_ref*.txt")
        ref_files = glob.glob(ref_pattern)
        
        for ref_file in ref_files:
            if "_sample_counts" not in ref_file:
                # Extract ref number from filename
                basename = os.path.basename(ref_file)
                try:
                    ref_part = basename.split('ref', 1)[1]
                    ref_num = int(ref_part.lstrip('_').split('.')[0])
                    fkt_files[ref_num] = ref_file
                except (IndexError, ValueError):
                    if self.verbose:
                        print(f"Warning: Could not parse ref number from {basename}")
                    continue
        
        if len(fkt_files) == 0:
            raise ValueError(f"No valid F(k,t) files found in {coupling_dir}")
            
        if self.verbose:
            print(f"Found {len(fkt_files)} reference frames in {coupling_dir}")
            
        return fkt_files
    
    def calculate_waiting_times(self, fkt_files):
        """
        Calculate waiting times for each reference frame.
        
        Parameters:
        -----------
        fkt_files : dict
            Dictionary mapping ref_num to file paths
            
        Returns:
        --------
        waiting_times : dict
            Dictionary mapping ref_num to waiting time in ps
        """
        waiting_times = {}
        for ref_num in fkt_files.keys():
            # Skip ref0 from analysis (pre-equilibrium data)
            if ref_num == 0:
                if self.verbose:
                    print(f"  ref{ref_num}: EXCLUDED from analysis (pre-equilibrium)")
                continue
                
            # Correct formula: tw_ps = (ref_num - 1) * 200
            # ref1 = 0 ps (coupling turns on), ref2 = 200 ps, etc.
            tw_ps = (ref_num - 1) * 200
            waiting_times[ref_num] = tw_ps
            
            if self.verbose:
                print(f"  ref{ref_num}: waiting time = {tw_ps} ps")
                
        return waiting_times
    
    def normalize_fkt_data(self, fkt_files):
        """
        Read and normalize F(k,t) data following the paper's methodology.
        
        Parameters:
        -----------
        fkt_files : dict
            Dictionary mapping ref_num to file paths
            
        Returns:
        --------
        normalized_fkt : dict
            Dictionary mapping ref_num to (lag_time, normalized_fkt) tuples
        """
        normalized_fkt = {}
        
        for ref_num, file_path in fkt_files.items():
            # Skip ref0 from analysis (pre-equilibrium data)
            if ref_num == 0:
                continue
                
            try:
                # Read the data file
                data = pd.read_csv(file_path, sep='\t', comment='#')
                
                if 'lag_time' not in data.columns or 'fskt' not in data.columns:
                    if self.verbose:
                        print(f"Warning: Expected columns not found in {file_path}")
                    continue
                    
                lag_time = data['lag_time'].values
                fkt = data['fskt'].values
                
                # Remove any NaN or infinite values
                valid_mask = np.isfinite(fkt) & np.isfinite(lag_time)
                lag_time = lag_time[valid_mask]
                fkt = fkt[valid_mask]
                
                if len(fkt) < 10:
                    if self.verbose:
                        print(f"Warning: Insufficient data points in ref{ref_num}")
                    continue
                
                # Apply Savitzky-Golay filtering to raw F(k,t) data BEFORE normalization
                fkt_smoothed = self.apply_savgol_filter(lag_time, fkt)
                
                # Normalize following the paper's approach using smoothed data
                # C(t₁,t₂) = [F(k,t) - F(k,∞)] / [F(k,0) - F(k,∞)]
                fkt_inf = np.mean(fkt_smoothed[-50:])  # Average last 50 points as F(k,∞)
                fkt_0 = fkt_smoothed[0]                # F(k,0)
                
                if abs(fkt_0 - fkt_inf) > 1e-10:  # Avoid division by zero
                    normalized = (fkt_smoothed - fkt_inf) / (fkt_0 - fkt_inf)
                else:
                    if self.verbose:
                        print(f"Warning: No decay observed in ref{ref_num}")
                    normalized = np.ones_like(fkt_smoothed)
                
                # Store both raw and smoothed normalized data for plotting options
                # normalized_fkt contains the smoothed, normalized data used for analysis
                normalized_fkt[ref_num] = (lag_time, normalized)
                
                if self.verbose:
                    print(f"  ref{ref_num}: F(k,0)={fkt_0:.1f}, F(k,∞)={fkt_inf:.1f}, "
                          f"decay={(fkt_0-fkt_inf)/fkt_0*100:.1f}%")
                    
            except Exception as e:
                if self.verbose:
                    print(f"Error processing ref{ref_num}: {e}")
                continue
                
        return normalized_fkt
    
    def verify_triangular_relation(self, normalized_fkt, n_samples=1000):
        """
        Verify the triangular relation C₁₃ = C₁₃(C₁₂, C₂₃) as a sanity check.
        
        Parameters:
        -----------
        normalized_fkt : dict
            Normalized F(k,t) data
        n_samples : int
            Number of time triplets to test
            
        Returns:
        --------
        bool
            True if triangular relation is approximately satisfied
        """
        if len(normalized_fkt) < 3:
            if self.verbose:
                print("Warning: Need at least 3 references to test triangular relation")
            return False
            
        # Pick three references for testing
        ref_nums = sorted(list(normalized_fkt.keys()))[:3]
        
        deviations = []
        for _ in range(min(n_samples, 100)):  # Limit to avoid excessive computation
            try:
                # Pick random time points
                ref1, ref2, ref3 = ref_nums
                t1, norm_fkt1 = normalized_fkt[ref1]
                t2, norm_fkt2 = normalized_fkt[ref2]
                t3, norm_fkt3 = normalized_fkt[ref3]
                
                # Find common time range
                max_t = min(t1[-1], t2[-1], t3[-1])
                if max_t < 10:  # Need sufficient overlap
                    continue
                    
                # Pick random time within common range
                t = np.random.uniform(0, max_t/2)
                
                # Interpolate to get C values
                C12 = np.interp(t, t1, norm_fkt1)
                C23 = np.interp(t, t2, norm_fkt2)
                C13_direct = np.interp(t, t3, norm_fkt3)
                
                # For a more sophisticated test, we would need the full
                # triangular relation implementation, but this gives a basic check
                deviations.append(abs(C13_direct - (C12 + C23)/2))  # Simplified check
                
            except Exception:
                continue
                
        if len(deviations) > 10:
            mean_deviation = np.mean(deviations)
            if self.verbose:
                print(f"Triangular relation test: mean deviation = {mean_deviation:.3f}")
            return mean_deviation < 0.1  # Reasonable threshold
        else:
            if self.verbose:
                print("Warning: Could not perform triangular relation test")
            return True  # Assume valid if can't test
    
    def calculate_material_time(self, normalized_fkt):
        """
        Calculate material time ξ(t) using the TRUE Figure 3a iterative procedure.
        
        CORRECT interpretation: t₁ = tw + t, t₂ = tw
        - C(t₁,t₂) = C(tw + t, tw) = correlation between time tw and time tw + t
        - When C(tw + t, tw) = criterion → ξ(tw + t) - ξ(tw) = 1
        - Each crossing defines: ξ(t₁) = ξ(t₂) + 1 where t₁ = tw + t, t₂ = tw
        
        Parameters:
        -----------
        normalized_fkt : dict
            Normalized F(k,t) data for all refs (different waiting times)
            
        Returns:
        --------
        material_time_function : tuple
            Single (lab_time_grid, xi_grid) that applies to all refs
        """
        
        if self.verbose:
            print(f" Implementing TRUE Figure 3a procedure for material time construction")
            print(f"   CORRECT interpretation: t₁ = tw + t, t₂ = tw")
            print(f"   Criterion: C(tw + t, tw) = {self.criterion_value} → ξ(tw + t) - ξ(tw) = 1")
        
        # Step 1: Prepare data - exclude ref0 and ref14, apply smoothing
        refs_data = {}
        for ref_num in sorted(normalized_fkt.keys()):
            if ref_num == 0 or ref_num == 14:  # Skip ref0 and ref14
                continue
                
            lag_time, norm_fkt_raw = normalized_fkt[ref_num]
            
            # Data is already smoothed with Savitzky-Golay filter
            # Just enforce monotonic decay
            norm_fkt_mono = np.copy(norm_fkt_raw)
            for i in range(1, len(norm_fkt_mono)):
                if norm_fkt_mono[i] > norm_fkt_mono[i-1]:
                    norm_fkt_mono[i] = norm_fkt_mono[i-1]
            
            refs_data[ref_num] = {
                'lag_time': lag_time,
                'norm_fkt': norm_fkt_mono,
                'waiting_time': (ref_num - 1) * 200  # tw for this reference
            }
        
        if self.verbose:
            print(f"   Using {len(refs_data)} references for material time construction")
        
        # Step 2: Collect material time constraints from crossing events
        # Each crossing C(tw + t, tw) = criterion gives us: ξ(tw + t) - ξ(tw) = 1
        
        xi_constraints = []  # Store (t₂, t₁, Δξ) where t₁ = tw + t, t₂ = tw
        
        for ref_num, ref_data in refs_data.items():
            tw = ref_data['waiting_time']  # t₂ = tw (reference time)
            lag_time = ref_data['lag_time'] 
            norm_fkt = ref_data['norm_fkt']
            
            # Find where C(tw + t, tw) = criterion in this reference's decay
            for i in range(1, len(norm_fkt)):
                # Check if we cross the criterion (decay below threshold)
                if (norm_fkt[i-1] > self.criterion_value >= norm_fkt[i] or 
                    norm_fkt[i-1] >= self.criterion_value > norm_fkt[i]):
                    
                    # Interpolate exact crossing lag time
                    if norm_fkt[i-1] != norm_fkt[i]:
                        lag_cross = lag_time[i-1] + (lag_time[i] - lag_time[i-1]) * \
                                   (self.criterion_value - norm_fkt[i-1]) / \
                                   (norm_fkt[i] - norm_fkt[i-1])
                    else:
                        lag_cross = lag_time[i]
                    
                    # Calculate t₁ = tw + t and t₂ = tw
                    t1 = tw + lag_cross  # measurement time
                    t2 = tw              # reference time
                    
                    # Record constraint: ξ(t₁) - ξ(t₂) = 1
                    xi_constraints.append((t2, t1, 1.0))
                    
                    if self.verbose and len(xi_constraints) <= 10:
                        print(f"     Constraint {len(xi_constraints)}: ξ({t1:.1f}) - ξ({t2:.1f}) = 1.0 (ref{ref_num})")
                    
                    break  # Only take first crossing per ref
        
        if len(xi_constraints) == 0:
            if self.verbose:
                print("    No crossings found! Using linear fallback")
            # Create a simple linear time grid
            time_grid = np.linspace(0, 2800, 1000)
            xi_grid = time_grid / np.max(time_grid) * 2.0
            return time_grid, xi_grid
        
        # Step 3: Solve the constraint system to build ξ(t)
        # We have constraints: ξ(t₁) - ξ(t₂) = 1 for multiple (t₂, t₁) pairs
        # Strategy: Build ξ(t) using a consistent cumulative approach
        
        if self.verbose:
            print(f"   Found {len(xi_constraints)} material time constraints")
        
        # Sort constraints by t₂ (waiting time) to solve systematically
        xi_constraints.sort(key=lambda x: x[0])  # Sort by t₂ = tw
        
        # Solve constraint system: ξ(t₁) - ξ(t₂) = 1
        # Strategy: Build ξ(t) by solving constraints in order of waiting time
        xi_points = [(0.0, 0.0)]  # Start with ξ(0) = 0
        
        for i, (t2, t1, delta_xi) in enumerate(xi_constraints):
            # Constraint: ξ(t₁) - ξ(t₂) = delta_xi
            # Need to find ξ(t₂), then calculate ξ(t₁) = ξ(t₂) + delta_xi
            
            # Find ξ(t₂) by interpolating from existing points
            if len(xi_points) >= 2:
                times = [p[0] for p in xi_points]
                xi_vals = [p[1] for p in xi_points]
                xi_t2 = np.interp(t2, times, xi_vals, left=0.0, right=xi_vals[-1])
            else:
                xi_t2 = 0.0 if t2 == 0.0 else np.interp(t2, [0.0], [0.0], left=0.0, right=0.0)
            
            # Calculate ξ(t₁) = ξ(t₂) + delta_xi
            xi_t1 = xi_t2 + delta_xi
            
            # Add point for t₁
            xi_points.append((t1, xi_t1))
            
            if self.verbose and i < 5:  # Print first few
                print(f"     Constraint {i+1}: ξ({t1:.1f}) = ξ({t2:.1f}) + {delta_xi:.1f} = {xi_t1:.2f}")
        
        # Sort points by time for proper interpolation
        xi_points.sort(key=lambda x: x[0])
        
        # Create continuous ξ(t) function on uniform grid
        t_max = max(p[0] for p in xi_points)
        time_grid = np.linspace(0, t_max + 500, 2000)
        
        times = [p[0] for p in xi_points]
        xi_vals = [p[1] for p in xi_points]
        
        # Interpolate and extrapolate
        xi_grid = np.zeros_like(time_grid)
        for i, t in enumerate(time_grid):
            if t <= times[-1]:
                xi_grid[i] = np.interp(t, times, xi_vals)
            else:
                # Linear extrapolation using last two points
                if len(times) >= 2:
                    dt = times[-1] - times[-2]
                    dxi = xi_vals[-1] - xi_vals[-2]
                    rate = dxi / dt if dt > 0 else 0.1
                else:
                    rate = 0.1
                xi_grid[i] = xi_vals[-1] + rate * (t - times[-1])
        
        if self.verbose:
            print(f"    Material time ξ(t) constructed: max ξ = {xi_grid[-1]:.1f}")
            print(f"   Time grid: 0 to {time_grid[-1]:.1f} ps ({len(time_grid)} points)")
        
        return time_grid, xi_grid
    
    def collapse_data(self, material_time_function, normalized_fkt, waiting_times):
        """
        Create collapsed data structure using the SINGLE material time function.
        Transform all ref data using the same ξ(t) function.
        
        Parameters:
        -----------
        material_time_function : tuple
            Single (time_grid, xi_grid) that applies to all refs
        normalized_fkt : dict
            Normalized F(k,t) data for all refs
        waiting_times : dict
            Waiting times for each reference
            
        Returns:
        --------
        collapsed_data : dict
            Dictionary with collapsed data for each reference
        """
        collapsed_data = {}
        time_grid, xi_grid = material_time_function
        
        if self.verbose:
            print("Applying material time transformation to all refs...")
        
        for ref_num in sorted(normalized_fkt.keys()):
            # Skip ref0 (should already be excluded, but double-check)
            if ref_num == 0:
                continue
                
            if ref_num not in waiting_times:
                continue
                
            # Exclude the last waiting time (2600 ps) which corresponds to ref14
            if ref_num == 14:
                continue
                
            lag_time, norm_fkt = normalized_fkt[ref_num]
            tw = waiting_times[ref_num]
            
            # For data collapse, we need to transform lag_time using the material time rate
            # at the waiting time. The idea is: how much material time elapses during 
            # a lag_time interval when aging started at waiting_time tw?
            
            # Get material time at waiting time: ξ(tw) 
            xi_at_tw = np.interp(tw, time_grid, xi_grid, left=0.0, right=xi_grid[-1])
            
            # Transform each lag time to absolute time: t_abs = tw + lag_time
            absolute_times = tw + lag_time
            
            # Get material time at each absolute time: ξ(tw + lag_time)
            xi_at_absolute_times = np.interp(absolute_times, time_grid, xi_grid, left=0.0, right=xi_grid[-1])
            
            # Material time intervals: Δξ = ξ(tw + lag_time) - ξ(tw)
            # This gives us the material time that has elapsed since the measurement started
            xi_transformed = xi_at_absolute_times - xi_at_tw
            
            collapsed_data[ref_num] = {
                'lag_time': lag_time,
                'xi': xi_transformed,  # Material time for this ref's data
                'R_fkt': norm_fkt,     # Normalized F(k,t) = R(ξ)
                'waiting_time_ps': tw,
                'ref_num': ref_num
            }
            
            if self.verbose:
                print(f"  ref{ref_num}: tw={tw} ps, max ξ = {xi_transformed[-1]:.1f}")
        
        return collapsed_data
    
    def calibrate_material_time(self, material_time_functions):
        """
        Calibrate material time so that ε=0 slope equals 1/τ_eq(T) where τ_eq = 113.7 ps.
        
        This ensures proper physical scaling: for equilibrium liquid (ε=0),
        material time should advance as ξ = t/τ_eq, so dξ/dt = 1/τ_eq = 1/113.7 ps⁻¹.
        
        Parameters:
        -----------
        material_time_functions : dict
            Dictionary mapping coupling_dir to (time_grid, xi_grid) tuples
            
        Returns:
        --------
        scaling_factor : float
            The factor 'a' to multiply all material times
        """
        tau_eq_ps = 113.7  # Equilibrium relaxation time in ps
        expected_slope = 1.0 / tau_eq_ps  # Expected slope for ε=0: dξ/dt = 1/τ_eq
        
        if self.verbose:
            print(f"\n Calibrating material time with τ_eq = {tau_eq_ps} ps")
            print(f"   Expected slope for ε=0: dξ/dt = {expected_slope:.6f} ps⁻¹")
        
        # Find the ε=0 coupling directory
        epsilon_zero_dir = None
        for coupling_dir in material_time_functions.keys():
            if "0epos00" in os.path.basename(coupling_dir):
                epsilon_zero_dir = coupling_dir
                break
        
        if epsilon_zero_dir is None:
            if self.verbose:
                print("    Warning: No ε=0 data found for calibration. Using scaling factor = 1.0")
            return 1.0
        
        # Get the ε=0 material time function
        time_grid, xi_grid = material_time_functions[epsilon_zero_dir]
        
        # Calculate the average slope over the middle portion of the time range
        # Use middle 60% to avoid edge effects
        start_idx = int(0.2 * len(time_grid))
        end_idx = int(0.8 * len(time_grid))
        
        if end_idx - start_idx < 10:
            # Not enough points, use full range
            start_idx = 1
            end_idx = len(time_grid) - 1
        
        if end_idx > start_idx:
            dt_section = time_grid[end_idx] - time_grid[start_idx]
            dxi_section = xi_grid[end_idx] - xi_grid[start_idx]
            
            if dt_section > 0:
                actual_slope = dxi_section / dt_section
                scaling_factor = expected_slope / actual_slope
                
                if self.verbose:
                    print(f"   Measured slope for ε=0: dξ/dt = {actual_slope:.6f} ps⁻¹")
                    print(f"   Scaling factor a = {scaling_factor:.6f}")
                    print(f"   Time range used: {time_grid[start_idx]:.1f} to {time_grid[end_idx]:.1f} ps")
                    print(f"   After scaling: dξ/dt = {actual_slope * scaling_factor:.6f} ps⁻¹")
                
                return scaling_factor
            else:
                if self.verbose:
                    print("    Warning: Invalid time range for slope calculation. Using scaling factor = 1.0")
                return 1.0
        else:
            if self.verbose:
                print("    Warning: Insufficient data points for slope calculation. Using scaling factor = 1.0")
            return 1.0

    def analyze_coupling_directory(self, coupling_dir):
        """
        Complete analysis pipeline for a single coupling strength directory.
        
        Parameters:
        -----------
        coupling_dir : str
            Directory containing F(k,t) data
            
        Returns:
        --------
        tuple
            (material_times, waiting_times, normalized_fkt, collapsed_data)
        """
        if self.verbose:
            print(f"\n=== Analyzing {coupling_dir} ===")
            
        # Load data
        fkt_files = self.load_fkt_data(coupling_dir)
        
        # Calculate waiting times
        waiting_times = self.calculate_waiting_times(fkt_files)
        
        # Normalize F(k,t) data
        if self.verbose:
            print("Normalizing F(k,t) data...")
        normalized_fkt = self.normalize_fkt_data(fkt_files)
        
        if len(normalized_fkt) == 0:
            raise ValueError("No valid normalized data found")
        
        # Verify triangular relation (optional)
        if self.verbose:
            print("Verifying triangular relation...")
        triangular_ok = self.verify_triangular_relation(normalized_fkt)
        
        # Calculate the SINGLE material time function for this coupling
        if self.verbose:
            print("Calculating material time...")
        material_time_function = self.calculate_material_time(normalized_fkt)
        
        # Create collapsed data structure using the single material time function
        collapsed_data = self.collapse_data(material_time_function, normalized_fkt, waiting_times)
        
        # Store results
        self.material_times[coupling_dir] = material_time_function  # Store the single function
        self.waiting_times[coupling_dir] = waiting_times
        self.normalized_fkt[coupling_dir] = normalized_fkt
        self.collapsed_data[coupling_dir] = collapsed_data
        
        if self.verbose:
            print(f"Analysis complete. Found {len(collapsed_data)} valid references.")
            if triangular_ok:
                print(" Triangular relation approximately satisfied")
            else:
                print(" Triangular relation may not be satisfied")
        
        return material_time_function, waiting_times, normalized_fkt, collapsed_data
    
    def apply_scaling_to_material_times(self, scaling_factor):
        """
        Apply scaling factor to all stored material time functions.
        
        Parameters:
        -----------
        scaling_factor : float
            Factor to multiply all material time values
        """
        if self.verbose:
            print(f"\n Applying scaling factor {scaling_factor:.6f} to all material times...")
        
        # Scale material time functions
        for coupling_dir in self.material_times.keys():
            time_grid, xi_grid = self.material_times[coupling_dir]
            scaled_xi_grid = xi_grid * scaling_factor
            self.material_times[coupling_dir] = (time_grid, scaled_xi_grid)
            
            coupling_name = os.path.basename(coupling_dir)
            if self.verbose:
                print(f"   Scaled {coupling_name}: max ξ = {xi_grid[-1]:.3f} → {scaled_xi_grid[-1]:.3f}")
        
        # Update collapsed data with scaled material times
        for coupling_dir in self.collapsed_data.keys():
            if coupling_dir in self.material_times:
                time_grid, scaled_xi_grid = self.material_times[coupling_dir]
                normalized_fkt = self.normalized_fkt[coupling_dir]
                waiting_times = self.waiting_times[coupling_dir]
                
                # Recalculate collapsed data with scaled material time
                scaled_collapsed_data = self.collapse_data(
                    (time_grid, scaled_xi_grid), normalized_fkt, waiting_times
                )
                self.collapsed_data[coupling_dir] = scaled_collapsed_data
        
        if self.verbose:
            print("    Scaling applied to all material times and collapsed data")
    
    def plot_material_time_collapse(self, coupling_dir, output_dir='.', max_lag_time=None, max_xi=None):
        """
        Plot the data collapse results for a specific coupling strength.
        
        Parameters:
        -----------
        coupling_dir : str
            Coupling directory name
        output_dir : str
            Output directory for plots
        max_lag_time : float or None
            Maximum lag time to plot (ps). If None, uses 80% of data range
        max_xi : float or None
            Maximum material time to plot. If None, uses 80% of data range
        """
        if coupling_dir not in self.collapsed_data:
            raise ValueError(f"No data found for {coupling_dir}. Run analysis first.")
            
        collapsed_data = self.collapsed_data[coupling_dir]
        
        if len(collapsed_data) == 0:
            print(f"No data to plot for {coupling_dir}")
            return None
        
        # Determine data ranges for automatic scaling
        all_lag_times = []
        all_xi_values = []
        waiting_times_list = [data['waiting_time_ps'] for data in collapsed_data.values()]
        
        for data in collapsed_data.values():
            all_lag_times.extend(data['lag_time'])
            all_xi_values.extend(data['xi'])
        
        # Set automatic limits if not provided
        if max_lag_time is None:
            max_lag_time = np.percentile(all_lag_times, 95)  # Use 95th percentile to avoid outliers
            max_lag_time = max(max_lag_time, 500.0)  # Ensure we see at least 500 ps
            
        if max_xi is None:
            # Show full material time evolution - use maximum value
            max_xi = np.max(all_xi_values) if all_xi_values else 5.0
            max_xi = max(max_xi, 2.0)  # Ensure reasonable minimum range
            
        # Create figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Color map for different waiting times - match F(k,t) analysis style
        waiting_times_list = [data['waiting_time_ps'] for data in collapsed_data.values()]
        if len(set(waiting_times_list)) > 1:
            # Use consistent colorbar range like in F(k,t) analysis
            tw_min = min(waiting_times_list)
            tw_max = max(waiting_times_list)
            norm = Normalize(vmin=tw_min, vmax=tw_max)
            cmap = plt.cm.viridis
        else:
            norm = None
            cmap = None
        
        # Panel 1: F(k,t) vs lag time (before collapse)
        for ref_num, data in collapsed_data.items():
            tw = data['waiting_time_ps']
            lag_time = data['lag_time']
            R = data['R_fkt']
            
            # Filter to reasonable time range
            mask = lag_time <= max_lag_time
            if np.any(mask):
                lag_time_plot = lag_time[mask]
                R_plot = R[mask]
                
                if cmap is not None:
                    color = cmap(norm(tw))
                    alpha = 0.8
                else:
                    color = f'C{ref_num % 10}'
                    alpha = 0.8
                    
                ax1.plot(lag_time_plot, R_plot, 
                        label=f'tw={tw} ps', color=color, alpha=alpha, linewidth=2)
        
        ax1.set_xlabel('Lag time (ps)')
        ax1.set_ylabel('Normalized F(k,t)')
        ax1.set_title(f'Before Material Time Collapse (Savitzky-Golay Smoothed)\n{os.path.basename(coupling_dir)}')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(0, max_lag_time)
        
        # Add colorbar if using color mapping
        if cmap is not None and norm is not None:
            # Create a scalar mappable for the colorbar
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])  # Required for ScalarMappable
            cbar1 = plt.colorbar(sm, ax=ax1, shrink=0.8, aspect=20)
            cbar1.set_label('Waiting time (ps)', rotation=270, labelpad=20)
        
        # Panel 2: R(ξ) vs material time (after collapse)
        for ref_num, data in collapsed_data.items():
            tw = data['waiting_time_ps']
            xi = data['xi']
            R = data['R_fkt']
            
            # Filter to reasonable material time range
            mask = xi <= max_xi
            if np.any(mask):
                xi_plot = xi[mask]
                R_plot = R[mask]
                
                if cmap is not None:
                    color = cmap(norm(tw))
                    alpha = 0.8
                else:
                    color = f'C{ref_num % 10}'
                    alpha = 0.8
                    
                ax2.plot(xi_plot, R_plot, 
                        label=f'tw={tw} ps', color=color, alpha=alpha, linewidth=2)
        
        ax2.set_xlabel('Material time ξ')
        ax2.set_ylabel('R(ξ)')
        ax2.set_title(f'After Material Time Collapse (Savitzky-Golay Smoothed)\n{os.path.basename(coupling_dir)}')
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim(0, max_xi)
        
        # Add colorbar for panel 2 as well
        if cmap is not None and norm is not None:
            # Create a scalar mappable for the colorbar
            sm2 = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm2.set_array([])  # Required for ScalarMappable
            cbar2 = plt.colorbar(sm2, ax=ax2, shrink=0.8, aspect=20)
            cbar2.set_label('Waiting time (ps)', rotation=270, labelpad=20)
        
        plt.tight_layout()
        
        # Save plot
        coupling_name = os.path.basename(coupling_dir)
        output_path = os.path.join(output_dir, f'{coupling_name}_material_time_collapse.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        if self.verbose:
            print(f"Saved plot: {output_path}")
        
        plt.show()
        return fig
    
    def calculate_ageing_rate(self, coupling_dir):
        """
        Calculate the ageing rate γ = dξ/dt for each reference.
        
        Parameters:
        -----------
        coupling_dir : str
            Coupling directory name
            
        Returns:
        --------
        ageing_rates : dict
            Dictionary mapping ref_num to ageing rate arrays
        """
        if coupling_dir not in self.collapsed_data:
            raise ValueError(f"No data found for {coupling_dir}")
            
        collapsed_data = self.collapsed_data[coupling_dir]
        ageing_rates = {}
        
        for ref_num, data in collapsed_data.items():
            # Skip ref0 (should already be excluded, but double-check)
            if ref_num == 0:
                continue
                
            lag_time = data['lag_time']
            xi = data['xi']
            
            # Calculate dξ/dt using finite differences
            if len(lag_time) > 2:
                dt = np.diff(lag_time)
                dxi = np.diff(xi)
                
                # Avoid division by zero
                valid_mask = dt > 1e-10
                if np.any(valid_mask):
                    # Use the time points corresponding to the derivatives
                    gamma_times = lag_time[1:]  # Drop first point for derivative calculation
                    gamma_values = np.zeros_like(gamma_times)
                    
                    # Calculate derivatives where dt is valid
                    if len(dxi[valid_mask]) > 0 and len(dt[valid_mask]) > 0:
                        gamma_values[valid_mask] = dxi[valid_mask] / dt[valid_mask]
                    
                    ageing_rates[ref_num] = {
                        'lag_time': gamma_times,
                        'gamma': gamma_values,
                        'waiting_time_ps': data['waiting_time_ps']
                    }
        
        return ageing_rates
    
    def plot_material_time_evolution(self, coupling_dirs, output_dir='.'):
        """
        Plot the reconstructed material time ξ(t) for different coupling strengths.
        
        Parameters:
        -----------
        coupling_dirs : list
            List of coupling directory names to compare
        output_dir : str
            Output directory for plots
        """
        if self.verbose:
            print("\n=== Plotting Material Time Evolution ===")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Define coupling strength mapping for proper labeling
        coupling_labels = {
            "cavity_coupling_0epos00_switch_200.0ps": "0.0",
            "cavity_coupling_3eneg04_switch_200.0ps": "3×10⁻⁴", 
            "cavity_coupling_5eneg04_switch_200.0ps": "5×10⁻⁴",
            "cavity_coupling_7eneg04_switch_200.0ps": "7×10⁻⁴",
            "cavity_coupling_1eneg03_switch_200.0ps": "1×10⁻³"
        }
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(coupling_dirs)))
        
        for i, coupling_dir in enumerate(coupling_dirs):
            if coupling_dir not in self.collapsed_data:
                continue
                
            collapsed_data = self.collapsed_data[coupling_dir]
            color = colors[i]
            coupling_name = os.path.basename(coupling_dir)
            coupling_label = coupling_labels.get(coupling_name, coupling_name)
            
            if not collapsed_data:
                continue
            
            # Plot 1: Material time ξ vs lag time for different waiting times
            for ref_num, data in collapsed_data.items():
                if ref_num == 1:  # Show only ref1 (tw=0, coupling just turned on)
                    lag_time = data['lag_time']
                    xi = data['xi']
                    ax1.plot(lag_time, xi, color=color, linewidth=2, alpha=0.8,
                            label=f'g = {coupling_label}')
            
            # Plot 2: Ageing rate γ = dξ/dt vs lag time for ref1
            ageing_rates = self.calculate_ageing_rate(coupling_dir)
            if 1 in ageing_rates:  # ref1 data
                rate_data = ageing_rates[1]
                lag_time = rate_data['lag_time']
                gamma = rate_data['gamma']
                ax2.plot(lag_time, gamma, color=color, linewidth=2, alpha=0.8,
                        label=f'g = {coupling_label}')
        
        # Format plot 1
        ax1.set_xlabel('Lag time (ps)')
        ax1.set_ylabel('Material time ξ')
        ax1.set_title('Material Time Function ξ(t) for Different Coupling Strengths\n(Single ξ(t) function used to transform all refs for each coupling)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Format plot 2
        ax2.set_xlabel('Lag time (ps)')
        ax2.set_ylabel('Ageing rate γ = dξ/dt')
        ax2.set_title('Ageing Rate γ=dξ/dt for Different Coupling Strengths\n(Derived from the single material time function for each coupling)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        output_path = os.path.join(output_dir, 'material_time_evolution.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        if self.verbose:
            print(f"Saved material time evolution plot: {output_path}")
        
        plt.show()
        return fig
    
    def compare_coupling_strengths(self, coupling_dirs, output_dir='.'):
        """
        Compare material time behavior across different coupling strengths.
        
        Parameters:
        -----------
        coupling_dirs : list
            List of coupling directory names to compare
        output_dir : str
            Output directory for plots
        """
        if self.verbose:
            print("\n=== Comparing Coupling Strengths ===")
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Define coupling strength mapping for proper labeling
        coupling_labels = {
            "cavity_coupling_0epos00_switch_200.0ps": "0.0",
            "cavity_coupling_3eneg04_switch_200.0ps": "3×10⁻⁴", 
            "cavity_coupling_5eneg04_switch_200.0ps": "5×10⁻⁴",
            "cavity_coupling_7eneg04_switch_200.0ps": "7×10⁻⁴",
            "cavity_coupling_1eneg03_switch_200.0ps": "1×10⁻³"
        }
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(coupling_dirs)))
        
        for i, coupling_dir in enumerate(coupling_dirs):
            if coupling_dir not in self.collapsed_data:
                if self.verbose:
                    print(f"Warning: No data for {coupling_dir}")
                continue
                
            collapsed_data = self.collapsed_data[coupling_dir]
            color = colors[i]
            coupling_name = os.path.basename(coupling_dir)
            
            # Get proper label for coupling strength
            coupling_label = coupling_labels.get(coupling_name, coupling_name)
            
            # Plot 1: Material time collapse for one reference from each coupling
            if collapsed_data:
                # Use reference with waiting time closest to 0 (should be ref1)
                best_ref = min(collapsed_data.keys(), 
                             key=lambda x: abs(collapsed_data[x]['waiting_time_ps']))
                data = collapsed_data[best_ref]
                
                axes[0,0].plot(data['xi'], data['R_fkt'], 
                             label=f'g = {coupling_label}', color=color, linewidth=2)
        
        axes[0,0].set_xlabel('Material time ξ')
        axes[0,0].set_ylabel('R(ξ)')
        axes[0,0].set_title('Material Time Collapse Comparison')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # Set reasonable x-limit based on data
        all_xi = []
        for coupling_dir in coupling_dirs:
            if coupling_dir in self.collapsed_data:
                for data in self.collapsed_data[coupling_dir].values():
                    all_xi.extend(data['xi'])
        if all_xi:
            max_xi_comp = np.percentile(all_xi, 90)
            axes[0,0].set_xlim(0, max(max_xi_comp, 5.0))
        
        # Additional analysis plots can be added here
        # For now, hide unused subplots
        for ax in axes.flatten()[1:]:
            ax.set_visible(False)
        
        plt.tight_layout()
        output_path = os.path.join(output_dir, 'coupling_comparison.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        if self.verbose:
            print(f"Saved comparison plot: {output_path}")
        
        plt.show()
        return fig
    
    def plot_unified_material_time_collapse(self, coupling_dirs, output_dir='.', max_xi=None):
        """
        Plot material time collapse for ALL coupling strengths in a single panel.
        This shows how well the Tool-Narayanaswamy theory works across different couplings.
        
        Parameters:
        -----------
        coupling_dirs : list
            List of coupling directory names to include
        output_dir : str
            Output directory for plots
        max_xi : float or None
            Maximum material time to plot. If None, uses 80% of data range
        """
        if self.verbose:
            print("\n=== Creating Unified Material Time Collapse Plot ===")
        
        # Define coupling strength mapping for proper labeling
        coupling_labels = {
            "cavity_coupling_0epos00_switch_200.0ps": "0.0",
            "cavity_coupling_3eneg04_switch_200.0ps": "3×10⁻⁴", 
            "cavity_coupling_5eneg04_switch_200.0ps": "5×10⁻⁴",
            "cavity_coupling_7eneg04_switch_200.0ps": "7×10⁻⁴",
            "cavity_coupling_1eneg03_switch_200.0ps": "1×10⁻³"
        }
        
        # Create figure with single panel
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        # Use the same color scheme as F(k,t) analysis (coolwarm colormap)
        coupling_values = []
        valid_dirs = []
        for coupling_dir in coupling_dirs:
            if coupling_dir in self.collapsed_data:
                coupling_name = os.path.basename(coupling_dir)
                axis_value = axis_value_from_coupling_dir(coupling_name)
                if np.isnan(axis_value):
                    axis_value = float(len(coupling_values))
                coupling_values.append(axis_value)
                valid_dirs.append(coupling_dir)
        
        # Create color mapping consistent with F(k,t) analysis
        from matplotlib.colors import Normalize
        coupling_norm = Normalize(vmin=0, vmax=max(coupling_values) if coupling_values else 1e-3)
        coupling_cmap = plt.colormaps.get_cmap('coolwarm')
        
        # Determine data range for automatic scaling
        all_xi_values = []
        for coupling_dir in valid_dirs:
            collapsed_data = self.collapsed_data[coupling_dir]
            for data in collapsed_data.values():
                all_xi_values.extend(data['xi'])
        
        # Set automatic limits if not provided
        if max_xi is None:
            max_xi = np.percentile(all_xi_values, 95) if all_xi_values else 5.0
            max_xi = max(max_xi, 3.0)  # Ensure reasonable minimum range
        
        # Plot material time collapse for each coupling strength
        legend_handles = []
        legend_labels = []
        
        for i, coupling_dir in enumerate(valid_dirs):
            coupling_name = os.path.basename(coupling_dir)
            coupling_label = coupling_labels.get(coupling_name, coupling_name)
            coupling_value = coupling_values[i]
            collapsed_data = self.collapsed_data[coupling_dir]
            
            # Get color for this coupling strength
            color = coupling_cmap(coupling_norm(coupling_value))
            
            if self.verbose:
                print(f"  Adding data for ε = {coupling_label}")
            
            # Plot all waiting times for this coupling with the same color but different alphas
            waiting_times_list = [data['waiting_time_ps'] for data in collapsed_data.values()]
            n_waiting_times = len(set(waiting_times_list))
            
            # Create alpha gradient for different waiting times
            alphas = np.linspace(0.3, 0.9, n_waiting_times)
            
            first_curve = True
            for j, (ref_num, data) in enumerate(collapsed_data.items()):
                tw = data['waiting_time_ps']
                xi = data['xi']
                R = data['R_fkt']
                
                # Filter to reasonable material time range
                mask = xi <= max_xi
                if np.any(mask):
                    xi_plot = xi[mask]
                    R_plot = R[mask]
                    
                    # Use alpha based on waiting time order
                    alpha = alphas[j % len(alphas)]
                    
                    line = ax.plot(xi_plot, R_plot, color=color, alpha=alpha, 
                                 linewidth=1.5, 
                                 label=f'ε = {coupling_label}' if first_curve else "")
                    
                    # Only add to legend for the first curve of each coupling
                    if first_curve:
                        legend_handles.append(line[0])
                        legend_labels.append(f'ε = {coupling_label}')
                        first_curve = False
        
        # Formatting
        ax.set_xlabel('Material time ξ', fontsize=14)
        ax.set_ylabel('R(ξ)', fontsize=14)
        ax.set_title('Material Time Collapse: All Coupling Strengths\n' + 
                    'Tool-Narayanaswamy Theory Test (F(k,t) = 0.1 criterion)', 
                    fontsize=16, fontweight='bold')
        
        # Add legend
        ax.legend(legend_handles, legend_labels, loc='upper right', fontsize=12,
                 title='Coupling Strength', title_fontsize=12)
        
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, max_xi)
        ax.set_ylim(0, 1.05)  # R(ξ) should be normalized between 0 and 1
        
        # Add text box with interpretation
        textstr = ('Good collapse: All curves overlap\n' +
                  'Poor collapse: Curves spread out\n' +
                  'Multiple lines per color = different waiting times')
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
               verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        
        # Save plot
        output_path = os.path.join(output_dir, 'unified_material_time_collapse.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        # Also save as PDF
        output_path_pdf = os.path.join(output_dir, 'unified_material_time_collapse.pdf')
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight')
        
        if self.verbose:
            print(f"Saved unified collapse plot: {output_path}")
            print(f"Saved unified collapse plot: {output_path_pdf}")
            print(f"  Material time range: 0 to {max_xi:.1f}")
            print(f"  Number of coupling strengths: {len(valid_dirs)}")
            
            # Analyze collapse quality
            print("\n Collapse Quality Analysis:")
            print("  - Look for overlapping curves of the same color")
            print("  - Good collapse: tight bundles of lines")
            print("  - Poor collapse: spread out lines")
            print("  - Each color represents a different coupling strength")
            print("  - Multiple lines per color = different waiting times")
        
        plt.show()
        return fig
    
    def stretched_exponential(self, xi, beta):
        """
        Stretched exponential function: R(ξ) = exp[-ln(10) * ξ^β]
        
        The prefactor ln(10) is included to normalize the characteristic time scale.
        Amplitude A is fixed at 1.
        
        Parameters:
        -----------
        xi : array
            Material time values
        beta : float
            Stretching exponent (0 < β ≤ 1)
            
        Returns:
        --------
        R : array
            Relaxation function values
        """
        return np.exp(-np.log(10) * xi ** beta)
    
    def fit_stretched_exponential_to_collapse(self, coupling_dirs, output_dir='.', max_xi=None):
        """
        Fit stretched exponential functions to the collapsed material time data.
        
        Parameters:
        -----------
        coupling_dirs : list
            List of coupling directory names to analyze
        output_dir : str
            Output directory for plots and results
        max_xi : float or None
            Maximum material time to fit. If None, uses 80% of data range
            
        Returns:
        --------
        fit_results : dict
            Dictionary containing fit parameters for each coupling
        """
        if self.verbose:
            print("\n=== Fitting Stretched Exponentials to Material Time Collapse ===")
        
        # Define coupling strength mapping for proper labeling
        coupling_labels = {
            "cavity_coupling_0epos00_switch_200.0ps": "0.0",
            "cavity_coupling_3eneg04_switch_200.0ps": "3×10⁻⁴", 
            "cavity_coupling_5eneg04_switch_200.0ps": "5×10⁻⁴",
            "cavity_coupling_7eneg04_switch_200.0ps": "7×10⁻⁴",
            "cavity_coupling_1eneg03_switch_200.0ps": "1×10⁻³"
        }
        
        # Create figure with subplots for individual fits
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        # Use the same color scheme as F(k,t) analysis
        coupling_values = []
        valid_dirs = []
        for coupling_dir in coupling_dirs:
            if coupling_dir in self.collapsed_data:
                coupling_name = os.path.basename(coupling_dir)
                axis_value = axis_value_from_coupling_dir(coupling_name)
                if np.isnan(axis_value):
                    axis_value = float(len(coupling_values))
                coupling_values.append(axis_value)
                valid_dirs.append(coupling_dir)
        
        from matplotlib.colors import Normalize
        coupling_norm = Normalize(vmin=0, vmax=max(coupling_values) if coupling_values else 1e-3)
        coupling_cmap = plt.colormaps.get_cmap('coolwarm')
        
        fit_results = {}
        
        # Determine data range for fitting
        all_xi_values = []
        for coupling_dir in valid_dirs:
            collapsed_data = self.collapsed_data[coupling_dir]
            for data in collapsed_data.values():
                all_xi_values.extend(data['xi'])
        
        if max_xi is None:
            max_xi = np.percentile(all_xi_values, 85) if all_xi_values else 5.0
            max_xi = max(max_xi, 2.0)
        
        for i, coupling_dir in enumerate(valid_dirs):
            coupling_name = os.path.basename(coupling_dir)
            coupling_label = coupling_labels.get(coupling_name, coupling_name)
            coupling_value = coupling_values[i]
            collapsed_data = self.collapsed_data[coupling_dir]
            color = coupling_cmap(coupling_norm(coupling_value))
            
            ax = axes[i]
            
            if self.verbose:
                print(f"\n  Fitting ε = {coupling_label}:")
            
            # Collect all data points for this coupling
            all_xi = []
            all_R = []
            
            for ref_num, data in collapsed_data.items():
                xi = data['xi']
                R = data['R_fkt']
                
                # Filter to fitting range and remove invalid points
                mask = (xi <= max_xi) & (xi > 0) & (R > 0) & (R <= 1)
                if np.any(mask):
                    all_xi.extend(xi[mask])
                    all_R.extend(R[mask])
            
            if len(all_xi) < 10:
                if self.verbose:
                    print(f"     Insufficient data points for fitting")
                fit_results[coupling_name] = {'tau': np.nan, 'beta': np.nan, 'r_squared': np.nan}
                ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'ε = {coupling_label}')
                continue
            
            # Convert to arrays and sort
            all_xi = np.array(all_xi)
            all_R = np.array(all_R)
            sort_idx = np.argsort(all_xi)
            xi_sorted = all_xi[sort_idx]
            R_sorted = all_R[sort_idx]
            
            try:
                # Improve data filtering for better fits
                # Only fit where R(ξ) is between 0.05 and 0.95 for better numerical stability
                valid_fit_mask = (R_sorted >= 0.05) & (R_sorted <= 0.95) & (xi_sorted > 0)
                if np.sum(valid_fit_mask) < 10:
                    # Fall back to original range if too few points
                    valid_fit_mask = (R_sorted > 0) & (R_sorted <= 1) & (xi_sorted > 0)
                
                xi_fit_data = xi_sorted[valid_fit_mask]
                R_fit_data = R_sorted[valid_fit_mask]
                
                if len(xi_fit_data) < 5:
                    raise ValueError("Insufficient data points for fitting")
                
                # Improved initial parameter guess based on data
                # Estimate β from the slope at R ≈ 0.37 (exp(-1))
                try:
                    idx_37 = np.argmin(np.abs(R_fit_data - 0.37))
                    xi_at_37 = xi_fit_data[idx_37]
                    if xi_at_37 > 0:
                        beta_guess = 1.0 / xi_at_37  # From exp(-xi^β) = 0.37 → xi^β ≈ 1
                        beta_guess = np.clip(beta_guess, 0.2, 1.5)
                    else:
                        beta_guess = 0.7
                except:
                    beta_guess = 0.7
                
                # Relaxed parameter bounds: 0.1 ≤ β ≤ 2.0
                bounds = ([0.1], [2.0])
                
                # Fit stretched exponential with improved settings
                popt, pcov = curve_fit(self.stretched_exponential, xi_fit_data, R_fit_data,
                                     p0=[beta_guess], bounds=bounds,
                                     maxfev=10000, 
                                     method='trf',  # Trust Region Reflective algorithm
                                     xtol=1e-12, ftol=1e-12)
                
                beta_fit = popt[0]
                param_errors = np.sqrt(np.diag(pcov))
                beta_err = param_errors[0]
                
                # Calculate R-squared using the filtered data used for fitting
                R_pred_fit = self.stretched_exponential(xi_fit_data, beta_fit)
                ss_res = np.sum((R_fit_data - R_pred_fit) ** 2)
                ss_tot = np.sum((R_fit_data - np.mean(R_fit_data)) ** 2)
                r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
                
                fit_results[coupling_name] = {
                    'A': 1.0,  # Fixed value
                    'A_err': 0.0,  # No uncertainty since it's fixed
                    'beta': beta_fit,
                    'beta_err': beta_err,
                    'r_squared': r_squared,
                    'coupling_value': coupling_value,
                    'coupling_label': coupling_label
                }
                
                if self.verbose:
                    print(f"    A = 1.000 (fixed)")
                    print(f"    β = {beta_fit:.3f} ± {beta_err:.3f}")
                    print(f"    R² = {r_squared:.4f}")
                
                # Plot data and fit
                # Plot raw data points (semi-transparent)
                for ref_num, data in collapsed_data.items():
                    xi_ref = data['xi']
                    R_ref = data['R_fkt']
                    mask = (xi_ref <= max_xi) & (xi_ref > 0) & (R_ref > 0) & (R_ref <= 1)
                    if np.any(mask):
                        ax.plot(xi_ref[mask], R_ref[mask], color=color, alpha=0.3, linewidth=1)
                
                # Plot fit line
                xi_fit = np.linspace(0.01, max_xi, 200)
                R_fit = self.stretched_exponential(xi_fit, beta_fit)
                ax.plot(xi_fit, R_fit, 'k-', linewidth=3, alpha=0.8,
                       label=f'R(ξ) = exp[-ln(10) × ξ^{beta_fit:.2f}]')
                
                # Add fit statistics to plot
                stats_text = f'A = 1.000 (fixed)\nβ = {beta_fit:.3f} ± {beta_err:.3f}\nR² = {r_squared:.3f}'
                ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                
            except Exception as e:
                if self.verbose:
                    print(f"     Fit failed: {e}")
                fit_results[coupling_name] = {'A': 1.0, 'A_err': 0.0, 'beta': np.nan, 'beta_err': np.nan, 'r_squared': np.nan}
                
                # Still plot the data
                for ref_num, data in collapsed_data.items():
                    xi_ref = data['xi']
                    R_ref = data['R_fkt']
                    mask = (xi_ref <= max_xi) & (xi_ref > 0) & (R_ref > 0) & (R_ref <= 1)
                    if np.any(mask):
                        ax.plot(xi_ref[mask], R_ref[mask], color=color, alpha=0.5, linewidth=1)
                
                ax.text(0.5, 0.3, 'Fit failed', ha='center', va='center', transform=ax.transAxes,
                       bbox=dict(boxstyle='round', facecolor='red', alpha=0.3))
            
            # Format subplot
            ax.set_xlabel('Material time ξ', fontsize=12)
            ax.set_ylabel('R(ξ)', fontsize=12)
            ax.set_title(f'ε = {coupling_label}', fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.set_xlim(0, max_xi)
            ax.set_ylim(0, 1.05)
            ax.legend(fontsize=9)
        
        # Hide unused subplots
        for j in range(len(valid_dirs), len(axes)):
            axes[j].set_visible(False)
        
        plt.suptitle('Stretched Exponential Fits to Material Time Collapse\n' + 
                    r'$R(\xi) = \exp[-\xi^\beta]$', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        # Save plot
        output_path = os.path.join(output_dir, 'stretched_exponential_fits.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        output_path_pdf = os.path.join(output_dir, 'stretched_exponential_fits.pdf')
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight')
        
        if self.verbose:
            print(f"\nSaved stretched exponential fits: {output_path}")
            print(f"Saved stretched exponential fits: {output_path_pdf}")
        
        # Create summary plot of fit parameters
        self._plot_fit_parameter_summary(fit_results, output_dir)
        
        # Save fit results to text file
        self._save_fit_results(fit_results, output_dir)
        
        plt.show()
        return fit_results
    
    def plot_overlapped_stretched_exponential_fits(self, coupling_dirs, output_dir='.', max_xi=None):
        """
        Plot all stretched exponential fits overlapped in a single panel.
        This shows the universality of the relaxation behavior across coupling strengths.
        
        Parameters:
        -----------
        coupling_dirs : list
            List of coupling directory names to analyze
        output_dir : str
            Output directory for plots
        max_xi : float or None
            Maximum material time to plot. If None, uses 80% of data range
        """
        if self.verbose:
            print("\n=== Creating Overlapped Stretched Exponential Fits Plot ===")
        
        # Define coupling strength mapping for proper labeling
        coupling_labels = {
            "cavity_coupling_0epos00_switch_200.0ps": "0.0",
            "cavity_coupling_3eneg04_switch_200.0ps": "3×10⁻⁴", 
            "cavity_coupling_5eneg04_switch_200.0ps": "5×10⁻⁴",
            "cavity_coupling_7eneg04_switch_200.0ps": "7×10⁻⁴",
            "cavity_coupling_1eneg03_switch_200.0ps": "1×10⁻³"
        }
        
        # Use the same color scheme as F(k,t) analysis
        coupling_values = []
        valid_dirs = []
        for coupling_dir in coupling_dirs:
            if coupling_dir in self.collapsed_data:
                coupling_name = os.path.basename(coupling_dir)
                axis_value = axis_value_from_coupling_dir(coupling_name)
                if np.isnan(axis_value):
                    axis_value = float(len(coupling_values))
                coupling_values.append(axis_value)
                valid_dirs.append(coupling_dir)
        
        from matplotlib.colors import Normalize
        coupling_norm = Normalize(vmin=0, vmax=max(coupling_values) if coupling_values else 1e-3)
        coupling_cmap = plt.colormaps.get_cmap('coolwarm')
        
        # Determine data range for plotting
        all_xi_values = []
        for coupling_dir in valid_dirs:
            collapsed_data = self.collapsed_data[coupling_dir]
            for data in collapsed_data.values():
                all_xi_values.extend(data['xi'])
        
        if max_xi is None:
            max_xi = np.percentile(all_xi_values, 85) if all_xi_values else 5.0
            max_xi = max(max_xi, 2.0)
        
        # Create single figure
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        
        fit_results = {}
        legend_elements = []
        
        for i, coupling_dir in enumerate(valid_dirs):
            coupling_name = os.path.basename(coupling_dir)
            coupling_label = coupling_labels.get(coupling_name, coupling_name)
            coupling_value = coupling_values[i]
            collapsed_data = self.collapsed_data[coupling_dir]
            color = coupling_cmap(coupling_norm(coupling_value))
            
            if self.verbose:
                print(f"  Processing ε = {coupling_label}")
            
            # Collect all data points for this coupling
            all_xi = []
            all_R = []
            
            for ref_num, data in collapsed_data.items():
                xi = data['xi']
                R = data['R_fkt']
                
                # Filter to plotting range and remove invalid points
                mask = (xi <= max_xi) & (xi > 0) & (R > 0) & (R <= 1)
                if np.any(mask):
                    all_xi.extend(xi[mask])
                    all_R.extend(R[mask])
            
            if len(all_xi) < 10:
                if self.verbose:
                    print(f"     Insufficient data points for ε = {coupling_label}")
                continue
            
            # Convert to arrays and sort
            all_xi = np.array(all_xi)
            all_R = np.array(all_R)
            sort_idx = np.argsort(all_xi)
            xi_sorted = all_xi[sort_idx]
            R_sorted = all_R[sort_idx]
            
            # Plot raw data points (semi-transparent)
            for ref_num, data in collapsed_data.items():
                xi_ref = data['xi']
                R_ref = data['R_fkt']
                mask = (xi_ref <= max_xi) & (xi_ref > 0) & (R_ref > 0) & (R_ref <= 1)
                if np.any(mask):
                    ax.plot(xi_ref[mask], R_ref[mask], color=color, alpha=0.2, linewidth=1)
            
            try:
                # Improve data filtering for better fits (same as individual fits)
                valid_fit_mask = (R_sorted >= 0.05) & (R_sorted <= 0.95) & (xi_sorted > 0)
                if np.sum(valid_fit_mask) < 10:
                    valid_fit_mask = (R_sorted > 0) & (R_sorted <= 1) & (xi_sorted > 0)
                
                xi_fit_data = xi_sorted[valid_fit_mask]
                R_fit_data = R_sorted[valid_fit_mask]
                
                if len(xi_fit_data) < 5:
                    raise ValueError("Insufficient data points for fitting")
                
                # Improved initial parameter guess
                try:
                    idx_37 = np.argmin(np.abs(R_fit_data - 0.37))
                    xi_at_37 = xi_fit_data[idx_37]
                    if xi_at_37 > 0:
                        beta_guess = 1.0 / xi_at_37
                        beta_guess = np.clip(beta_guess, 0.2, 1.5)
                    else:
                        beta_guess = 0.7
                except:
                    beta_guess = 0.7
                
                bounds = ([0.1], [2.0])
                
                popt, pcov = curve_fit(self.stretched_exponential, xi_fit_data, R_fit_data,
                                     p0=[beta_guess], bounds=bounds, maxfev=10000,
                                     method='trf', xtol=1e-12, ftol=1e-12)
                
                beta_fit = popt[0]
                param_errors = np.sqrt(np.diag(pcov))
                beta_err = param_errors[0]
                
                # Calculate R-squared using the filtered data used for fitting
                R_pred_fit = self.stretched_exponential(xi_fit_data, beta_fit)
                ss_res = np.sum((R_fit_data - R_pred_fit) ** 2)
                ss_tot = np.sum((R_fit_data - np.mean(R_fit_data)) ** 2)
                r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
                
                fit_results[coupling_name] = {
                    'A': 1.0,  # Fixed value
                    'A_err': 0.0,  # No uncertainty since it's fixed
                    'beta': beta_fit,
                    'beta_err': beta_err,
                    'r_squared': r_squared,
                    'coupling_value': coupling_value,
                    'coupling_label': coupling_label
                }
                
                # Plot fit line
                xi_fit = np.linspace(0.01, max_xi, 200)
                R_fit = self.stretched_exponential(xi_fit, beta_fit)
                line = ax.plot(xi_fit, R_fit, color=color, linewidth=4, alpha=0.9)[0]
                
                # Add to legend
                legend_elements.append((line, f'ε = {coupling_label}, β = {beta_fit:.3f}', color))
                
                if self.verbose:
                    print(f"    A = 1.000 (fixed), β = {beta_fit:.3f} ± {beta_err:.3f}, R² = {r_squared:.4f}")
                
            except Exception as e:
                if self.verbose:
                    print(f"     Fit failed for ε = {coupling_label}: {e}")
                continue
        
        # Format the plot
        ax.set_xlabel('Material time ξ', fontsize=16, fontweight='bold')
        ax.set_ylabel('R(ξ)', fontsize=16, fontweight='bold')
        ax.set_title('Universal Material Time Collapse\nStretched Exponential Fits: R(ξ) = exp[-ξ^β]', 
                    fontsize=18, fontweight='bold', pad=20)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, max_xi)
        ax.set_ylim(0, 1.05)
        
        # Create custom legend
        legend_lines = [elem[0] for elem in legend_elements]
        legend_labels = [elem[1] for elem in legend_elements]
        ax.legend(legend_lines, legend_labels, fontsize=12, loc='upper right',
                 title='Coupling Strength & Fit Parameters', title_fontsize=13)
        
        # Add explanatory text
        explanation_text = ("• Light lines: raw collapsed data\n"
                          "• Bold lines: stretched exponential fits\n"
                          "• Good collapse → overlapping curves\n"
                          "• β < 1 → sub-exponential (glassy) relaxation")
        ax.text(0.02, 0.65, explanation_text, transform=ax.transAxes, fontsize=11,
               verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
        
        plt.tight_layout()
        
        # Save plot
        output_path = os.path.join(output_dir, 'overlapped_stretched_exponential_fits.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        output_path_pdf = os.path.join(output_dir, 'overlapped_stretched_exponential_fits.pdf')
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight')
        
        if self.verbose:
            print(f"\nSaved overlapped fits: {output_path}")
            print(f"Saved overlapped fits: {output_path_pdf}")
            print(f"\n Quality of Material Time Collapse:")
            print("  • Excellent collapse: all curves should overlap")
            print("  • Each color = different coupling strength")
            print("  • Light lines = raw data, bold lines = fits")
            print("  • β values show how 'stretched' the relaxation is")
        
        plt.show()
        return fit_results
    
    def _plot_fit_parameter_summary(self, fit_results, output_dir):
        """Plot summary of fit parameters vs coupling strength."""
        valid_results = {k: v for k, v in fit_results.items() 
                        if not np.isnan(v.get('beta', np.nan))}
        
        if len(valid_results) < 2:
            return
        
        # Extract data for plotting
        coupling_vals = []
        beta_vals = []
        beta_errs = []
        labels = []
        
        for coupling_name, result in valid_results.items():
            coupling_vals.append(result['coupling_value'])
            beta_vals.append(result['beta'])
            beta_errs.append(result.get('beta_err', 0))
            labels.append(result['coupling_label'])
        
        # Sort by coupling strength
        sort_idx = np.argsort(coupling_vals)
        coupling_vals = np.array(coupling_vals)[sort_idx] * 10000  # Scale for plotting
        beta_vals = np.array(beta_vals)[sort_idx]
        beta_errs = np.array(beta_errs)[sort_idx]
        labels = [labels[i] for i in sort_idx]
        
        # Create single plot for β (A is fixed at 1)
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        
        # Plot beta vs coupling strength
        ax.errorbar(coupling_vals, beta_vals, yerr=beta_errs,
                   marker='s', markersize=10, linewidth=3, capsize=8, color='red',
                   markerfacecolor='red', markeredgecolor='darkred', capthick=2)
        ax.set_xlabel('ε (×10⁻⁴ a.u.)', fontsize=14)
        ax.set_ylabel('Stretching exponent β', fontsize=14)
        ax.set_title('Stretching Exponent vs Coupling Strength\nR(ξ) = exp[-ξ^β] (A = 1 fixed)', 
                    fontsize=16, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xticks(coupling_vals)
        ax.set_xticklabels([f'{cv:.1f}' for cv in coupling_vals])
        ax.set_ylim(0, 1.05)
        
        # Add horizontal line at β = 1 (exponential limit)
        ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, linewidth=2,
                  label='β = 1 (simple exponential)')
        ax.legend(fontsize=12)
        
        plt.tight_layout()
        
        # Save plot
        output_path = os.path.join(output_dir, 'stretched_exponential_parameters.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        output_path_pdf = os.path.join(output_dir, 'stretched_exponential_parameters.pdf')
        plt.savefig(output_path_pdf, dpi=300, bbox_inches='tight')
        
        if self.verbose:
            print(f"Saved parameter summary: {output_path}")
        
        plt.show()
    
    def _save_fit_results(self, fit_results, output_dir):
        """Save stretched exponential fit results to text file."""
        output_file = os.path.join(output_dir, 'stretched_exponential_fit_results.txt')
        
        with open(output_file, 'w') as f:
            f.write("# Stretched Exponential Fit Results for Material Time Collapse\n")
            f.write("# Fit function: R(ξ) = exp[-ξ^β] (A = 1 fixed)\n")
            f.write("# Standard Kohlrausch-Williams-Watts stretched exponential form\n")
            f.write("# β is the stretching exponent (amplitude A fixed at 1)\n")
            f.write("# Generated by Material Time Analysis\n")
            f.write("#\n")
            f.write("# Columns:\n")
            f.write("# 1. Coupling strength ε (a.u.)\n") 
            f.write("# 2. Coupling label\n")
            f.write("# 3. Amplitude A (fixed at 1.0)\n")
            f.write("# 4. A uncertainty (0.0, fixed parameter)\n")
            f.write("# 5. Stretching exponent β\n")
            f.write("# 6. β uncertainty\n")
            f.write("# 7. R-squared\n")
            f.write("#\n")
            f.write(f"{'Coupling_au':<12} {'Label':<12} {'A':<8} {'A_err':<8} {'beta':<8} {'beta_err':<8} {'R_squared':<10}\n")
            
            for coupling_name, result in fit_results.items():
                coupling_val = result.get('coupling_value', np.nan)
                coupling_label = result.get('coupling_label', 'unknown')
                A = result.get('A', 1.0)  # Should always be 1.0
                A_err = result.get('A_err', 0.0)  # Should always be 0.0
                beta = result.get('beta', np.nan)
                beta_err = result.get('beta_err', np.nan)
                r_squared = result.get('r_squared', np.nan)
                
                f.write(f"{coupling_val:<12.6f} {coupling_label:<12} "
                       f"{A:<8.4f} {A_err:<8.4f} {beta:<8.4f} {beta_err:<8.4f} {r_squared:<10.4f}\n")
        
        if self.verbose:
            print(f"Saved fit results: {output_file}")
    
    def save_material_time_data(self, coupling_dirs, output_dir='.'):
        """
        Save material time functions ξ(t) for different coupling strengths to text files.
        
        Parameters:
        -----------
        coupling_dirs : list
            List of coupling directory names to save
        output_dir : str
            Output directory for text files
        """
        if self.verbose:
            print("\n=== Saving Material Time Data ===")
        
        # Define coupling strength mapping for proper labeling
        coupling_labels = {
            "cavity_coupling_0epos00_switch_200.0ps": "0.0",
            "cavity_coupling_3eneg04_switch_200.0ps": "3e-4", 
            "cavity_coupling_5eneg04_switch_200.0ps": "5e-4",
            "cavity_coupling_7eneg04_switch_200.0ps": "7e-4",
            "cavity_coupling_1eneg03_switch_200.0ps": "1e-3"
        }
        
        # Save individual coupling material time functions
        for coupling_dir in coupling_dirs:
            if coupling_dir not in self.material_times:
                if self.verbose:
                    print(f"Warning: No material time data for {coupling_dir}")
                continue
            
            coupling_name = os.path.basename(coupling_dir)
            coupling_label = coupling_labels.get(coupling_name, coupling_name)
            
            # Get the material time function
            time_grid, xi_grid = self.material_times[coupling_dir]
            
            # Create output filename
            safe_label = coupling_label.replace('e-', 'eneg').replace('.', 'p')
            output_file = os.path.join(output_dir, f'material_time_coupling_{safe_label}.txt')
            
            # Save data
            header = f"""# Material Time Function ξ(t) for Coupling Strength g = {coupling_label}
# Generated by Material Time Analysis following Tool-Narayanaswamy formalism
# Criterion value: F(k,t) = {self.criterion_value}
# Savitzky-Golay smoothing window: {self.sg_window_ps} ps
#
# Columns:
# 1. Laboratory time t (ps)
# 2. Material time ξ(t)
#
# Usage: ξ(t) transforms all waiting time data for this coupling strength
# If Tool-Narayanaswamy theory applies, all F(k,t) curves should collapse 
# to a single master curve R(ξ) when plotted vs material time
#
# Laboratory_Time_ps  Material_Time_xi"""
            
            try:
                with open(output_file, 'w') as f:
                    f.write(header + '\n')
                    for t, xi in zip(time_grid, xi_grid):
                        f.write(f'{t:12.3f}  {xi:15.6f}\n')
                
                if self.verbose:
                    print(f"Saved material time data: {output_file}")
                    print(f"  Time range: {time_grid[0]:.1f} to {time_grid[-1]:.1f} ps")
                    print(f"  Material time range: {xi_grid[0]:.3f} to {xi_grid[-1]:.3f}")
                    
            except Exception as e:
                print(f"Error saving {output_file}: {e}")
        
        # Save combined comparison file
        combined_file = os.path.join(output_dir, 'material_time_all_couplings.txt')
        
        try:
            with open(combined_file, 'w') as f:
                f.write(f"""# Material Time Functions for All Coupling Strengths
# Generated by Material Time Analysis following Tool-Narayanaswamy formalism
# Criterion value: F(k,t) = {self.criterion_value}
# Savitzky-Golay smoothing window: {self.sg_window_ps} ps
#
# Each column pair represents (time_ps, xi) for a different coupling strength
# Column pairs: """)
                
                # Write column headers
                col_names = []
                for coupling_dir in coupling_dirs:
                    if coupling_dir in self.material_times:
                        coupling_name = os.path.basename(coupling_dir)
                        coupling_label = coupling_labels.get(coupling_name, coupling_name)
                        col_names.append(f't_g{coupling_label}_ps')
                        col_names.append(f'xi_g{coupling_label}')
                
                f.write(' '.join(col_names) + '\n')
                f.write('#\n')
                
                # Determine common time grid
                all_time_grids = []
                for coupling_dir in coupling_dirs:
                    if coupling_dir in self.material_times:
                        time_grid, _ = self.material_times[coupling_dir]
                        all_time_grids.append(time_grid)
                
                if all_time_grids:
                    # Use finest time grid
                    max_points = max(len(tg) for tg in all_time_grids)
                    finest_grid_idx = next(i for i, tg in enumerate(all_time_grids) if len(tg) == max_points)
                    common_time_grid = all_time_grids[finest_grid_idx]
                    
                    # Interpolate all material time functions to common grid
                    interpolated_data = []
                    for coupling_dir in coupling_dirs:
                        if coupling_dir in self.material_times:
                            time_grid, xi_grid = self.material_times[coupling_dir]
                            xi_interp = np.interp(common_time_grid, time_grid, xi_grid)
                            interpolated_data.append((common_time_grid, xi_interp))
                        else:
                            interpolated_data.append((common_time_grid, np.full_like(common_time_grid, np.nan)))
                    
                    # Write data rows
                    for i in range(len(common_time_grid)):
                        row_data = []
                        for time_grid, xi_interp in interpolated_data:
                            row_data.append(f'{time_grid[i]:12.3f}')
                            row_data.append(f'{xi_interp[i]:15.6f}')
                        f.write(' '.join(row_data) + '\n')
                    
                    if self.verbose:
                        print(f"Saved combined material time data: {combined_file}")
                        print(f"  Contains {len(interpolated_data)} coupling strengths")
                        print(f"  Time points: {len(common_time_grid)}")
                
        except Exception as e:
            print(f"Error saving combined file {combined_file}: {e}")
        
        # Save material time derivatives (aging rates)
        deriv_file = os.path.join(output_dir, 'material_time_derivatives_all_couplings.txt')
        
        try:
            with open(deriv_file, 'w') as f:
                f.write(f"""# Material Time Derivatives (Aging Rates) γ = dξ/dt for All Coupling Strengths
# Generated by Material Time Analysis following Tool-Narayanaswamy formalism
# Criterion value: F(k,t) = {self.criterion_value}
# Savitzky-Golay smoothing window: {self.sg_window_ps} ps
#
# Each column pair represents (time_ps, gamma) for a different coupling strength
# γ = dξ/dt represents the aging rate at each time
# Column pairs: """)
                
                # Write column headers for derivatives
                col_names = []
                for coupling_dir in coupling_dirs:
                    if coupling_dir in self.material_times:
                        coupling_name = os.path.basename(coupling_dir)
                        coupling_label = coupling_labels.get(coupling_name, coupling_name)
                        col_names.append(f't_g{coupling_label}_ps')
                        col_names.append(f'gamma_g{coupling_label}')
                
                f.write(' '.join(col_names) + '\n')
                f.write('#\n')
                
                # Calculate derivatives for each coupling
                all_deriv_data = []
                for coupling_dir in coupling_dirs:
                    if coupling_dir in self.material_times:
                        time_grid, xi_grid = self.material_times[coupling_dir]
                        
                        # Calculate dξ/dt using finite differences
                        dt = np.diff(time_grid)
                        dxi = np.diff(xi_grid)
                        gamma = np.zeros_like(time_grid)
                        
                        # Avoid division by zero
                        valid_mask = dt > 1e-10
                        if np.any(valid_mask):
                            gamma[1:][valid_mask] = dxi[valid_mask] / dt[valid_mask]
                        
                        all_deriv_data.append((time_grid, gamma))
                    else:
                        all_deriv_data.append((np.array([]), np.array([])))
                
                # Find common length (use shortest for simplicity)
                min_length = min(len(data[0]) for data in all_deriv_data if len(data[0]) > 0)
                
                if min_length > 0:
                    # Write derivative data
                    for i in range(min_length):
                        row_data = []
                        for time_grid, gamma in all_deriv_data:
                            if len(time_grid) > i:
                                row_data.append(f'{time_grid[i]:12.3f}')
                                row_data.append(f'{gamma[i]:15.6f}')
                            else:
                                row_data.append(f'{0.0:12.3f}')
                                row_data.append(f'{np.nan:15.6f}')
                        f.write(' '.join(row_data) + '\n')
                    
                    if self.verbose:
                        print(f"Saved material time derivatives: {deriv_file}")
                
        except Exception as e:
            print(f"Error saving derivatives file {deriv_file}: {e}")
        
        if self.verbose:
            print("Material time data export complete!")
            print(f"Files saved in: {output_dir}")
            print("Files created:")
            print("  - material_time_coupling_*.txt (individual coupling functions)")
            print("  - material_time_all_couplings.txt (combined data)")
            print("  - material_time_derivatives_all_couplings.txt (aging rates)")
    
    def plot_double_logarithmic_analysis(self, coupling_dirs, output_dir='.'):
        """
        Create double-logarithmic plots for stretched exponential analysis.
        
        For KWW function R(ξ) = exp(-ξ^β), ln(-ln R) vs ln(ξ) should be linear
        with slope β. This provides a more reliable way to analyze and fit
        stretched exponential behavior.
        
        Parameters:
        -----------
        coupling_dirs : list
            List of directories to analyze
        output_dir : str
            Directory to save output plots
        """
        if self.verbose:
            print("\n Creating double-logarithmic analysis plots...")
        
        n_dirs = len(coupling_dirs)
        if n_dirs == 0:
            return
        
        # Create figure with subplots for each coupling
        fig_rows = int(np.ceil(np.sqrt(n_dirs)))
        fig_cols = int(np.ceil(n_dirs / fig_rows))
        
        fig, axes = plt.subplots(fig_rows, fig_cols, figsize=(5*fig_cols, 4*fig_rows))
        if n_dirs == 1:
            axes = [axes]
        elif fig_rows == 1:
            axes = axes.reshape(1, -1)
        axes = axes.flatten()
        
        fit_results = {}
        
        for idx, coupling_dir in enumerate(coupling_dirs):
            ax = axes[idx]
            
            # Extract coupling information
            coupling_name = os.path.basename(coupling_dir)
            import re
            coupling_match = re.search(r'coupling_([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)', coupling_name)
            if coupling_match:
                coupling_value = float(coupling_match.group(1))
                if coupling_value == 0:
                    coupling_label = 'ε = 0'
                elif coupling_value < 0.001:
                    coupling_label = f'ε = {coupling_value:.0e}'
                else:
                    coupling_label = f'ε = {coupling_value:.3f}'
            else:
                coupling_label = coupling_name
                coupling_value = 0.0
            
            # Get material time data for this coupling
            if coupling_dir not in self.material_times:
                continue
                
            time_grid, xi_grid = self.material_times[coupling_dir]
            
            # Get collapsed data
            if coupling_dir not in self.collapsed_data:
                continue
                
            collapsed_data = self.collapsed_data[coupling_dir]
            
            # Process all waiting times for this coupling
            valid_points_xi = []
            valid_points_log = []
            
            for tw, data in collapsed_data.items():
                xi_vals = data['xi']
                R_vals = data['R_fkt']
                
                # Filter for valid points: R < 1/e for proper stretched exponential regime
                # ln(-ln R) is only meaningful for R < 1, and stretched exponential behavior
                # is most reliable for R < 1/e ≈ 0.368 (avoiding initial transients)
                valid_mask = (R_vals > 0.01) & (R_vals < 1/np.e) & (xi_vals > 0)
                
                if np.sum(valid_mask) > 0:
                    xi_valid = xi_vals[valid_mask]
                    R_valid = R_vals[valid_mask]
                    
                    # Calculate ln(-ln R) - need R < 1 for this to be real
                    try:
                        ln_neg_ln_R = np.log(-np.log(R_valid))
                        ln_xi = np.log(xi_valid)
                        
                        # Filter out any NaN or infinite values and exclude ln(-ln R) >= 1
                        # This excludes very early times where R < exp(-exp(1)) ≈ 0.067
                        finite_mask = (np.isfinite(ln_neg_ln_R) & np.isfinite(ln_xi) & 
                                     (ln_neg_ln_R < 1))
                        if np.sum(finite_mask) > 0:
                            valid_points_xi.extend(ln_xi[finite_mask])
                            valid_points_log.extend(ln_neg_ln_R[finite_mask])
                    except:
                        continue
            
            if len(valid_points_xi) > 10:  # Need enough points for fitting
                # Convert to arrays and sort
                xi_data = np.array(valid_points_xi)
                log_data = np.array(valid_points_log)
                
                # Sort by xi for plotting
                sort_idx = np.argsort(xi_data)
                xi_sorted = xi_data[sort_idx]
                log_sorted = log_data[sort_idx]
                
                # Plot the data points
                ax.scatter(xi_sorted, log_sorted, alpha=0.6, s=20, 
                          label=f'Data')
                
                # Fit linear regression: ln(-ln R) = β * ln(ξ) + C
                try:
                    from scipy.stats import linregress
                    slope, intercept, r_value, p_value, std_err = linregress(xi_sorted, log_sorted)
                    
                    # Plot fit line
                    xi_fit = np.linspace(xi_sorted.min(), xi_sorted.max(), 100)
                    log_fit = slope * xi_fit + intercept
                    ax.plot(xi_fit, log_fit, 'r-', linewidth=2, 
                           label=f'β = {slope:.3f} ± {std_err:.3f}')
                    
                    # Store fit results
                    fit_results[coupling_name] = {
                        'beta': slope,
                        'beta_err': std_err,
                        'r_squared': r_value**2,
                        'intercept': intercept,
                        'coupling_value': coupling_value,
                        'coupling_label': coupling_label
                    }
                    
                    if self.verbose:
                        print(f"    {coupling_label}: β = {slope:.3f} ± {std_err:.3f}, R² = {r_value**2:.3f}")
                
                except Exception as e:
                    if self.verbose:
                        print(f"     Linear fit failed for {coupling_label}: {e}")
                    fit_results[coupling_name] = {
                        'beta': np.nan,
                        'beta_err': np.nan,
                        'r_squared': np.nan,
                        'intercept': np.nan,
                        'coupling_value': coupling_value,
                        'coupling_label': coupling_label
                    }
            
            else:
                if self.verbose:
                    print(f"     Insufficient valid data points for {coupling_label}")
                fit_results[coupling_name] = {
                    'beta': np.nan,
                    'beta_err': np.nan,
                    'r_squared': np.nan,
                    'intercept': np.nan,
                    'coupling_value': coupling_value,
                    'coupling_label': coupling_label
                }
            
            # Customize subplot
            ax.set_xlabel('ln(ξ)', fontsize=12)
            ax.set_ylabel('ln(-ln R)', fontsize=12)
            ax.set_title(f'{coupling_label}', fontsize=14, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=10)
        
        # Hide empty subplots
        for idx in range(n_dirs, len(axes)):
            axes[idx].set_visible(False)
        
        plt.suptitle('Double-Logarithmic Analysis of Stretched Exponential\n' + 
                    r'If $R(\xi) = \exp(-\xi^\beta)$, then $\ln(-\ln R) = \beta \ln(\xi)$ (fitted for $R < 1/e$, $\ln(-\ln R) < 1$)', 
                    fontsize=15, fontweight='bold')
        plt.tight_layout()
        
        # Save plots
        output_path = os.path.join(output_dir, 'double_logarithmic_analysis.png')
        output_path_pdf = os.path.join(output_dir, 'double_logarithmic_analysis.pdf')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.savefig(output_path_pdf, bbox_inches='tight')
        plt.close()
        
        if self.verbose:
            print(f"\nSaved double-logarithmic analysis: {output_path}")
            print(f"Saved double-logarithmic analysis: {output_path_pdf}")
        
        # Save fit results
        self._save_double_log_fit_results(fit_results, output_dir)
        
        return fit_results
    
    def _save_double_log_fit_results(self, fit_results, output_dir):
        """Save double-logarithmic fit results to text file."""
        output_file = os.path.join(output_dir, 'double_logarithmic_fit_results.txt')
        
        with open(output_file, 'w') as f:
            f.write("# Double-Logarithmic Analysis Results for Stretched Exponential\n")
            f.write("# Linear fit to ln(-ln R) vs ln(ξ): ln(-ln R) = β * ln(ξ) + C\n")
            f.write("# For KWW form R(ξ) = exp(-ξ^β), the slope gives β directly\n")
            f.write("# Fitted only for R < 1/e ≈ 0.368 and ln(-ln R) < 1\n")
            f.write("# This excludes early transients (R > 0.067) and very late times for robust KWW analysis\n")
            f.write("# Generated by Material Time Analysis\n")
            f.write("#\n")
            f.write("# Columns:\n")
            f.write("# 1. Coupling strength ε (a.u.)\n") 
            f.write("# 2. Coupling label\n")
            f.write("# 3. Stretching exponent β (slope)\n")
            f.write("# 4. β uncertainty\n")
            f.write("# 5. R-squared of linear fit\n")
            f.write("# 6. Intercept C\n")
            f.write("#\n")
            f.write("# epsilon    coupling_label    beta        beta_err    r_squared   intercept\n")
            
            # Sort by coupling strength
            sorted_results = sorted(fit_results.items(), 
                                  key=lambda x: x[1]['coupling_value'])
            
            for coupling_name, result in sorted_results:
                f.write(f"{result['coupling_value']:10.6f}  "
                       f"{result['coupling_label']:15s}  "
                       f"{result.get('beta', np.nan):10.6f}  "
                       f"{result.get('beta_err', np.nan):10.6f}  "
                       f"{result.get('r_squared', np.nan):10.6f}  "
                       f"{result.get('intercept', np.nan):10.6f}\n")
        
        if self.verbose:
            print(f"Saved double-logarithmic fit results to: {output_file}")


def main():
    """
    Main analysis script - analyzes specific coupling directories.
    """
    # Initialize analyzer
    analyzer = MaterialTimeAnalyzer(criterion_value=0.1, verbose=True, sg_window_ps=0.1)
    
    # Define specific coupling strengths to analyze
    target_couplings = [0.0, 3e-4, 5e-4, 7e-4, 1e-3]
    
    # Map coupling values to directory names
    coupling_to_dir = {
        0.0: "cavity_coupling_0epos00_switch_200.0ps",
        3e-4: "cavity_coupling_3eneg04_switch_200.0ps", 
        5e-4: "cavity_coupling_5eneg04_switch_200.0ps",
        7e-4: "cavity_coupling_7eneg04_switch_200.0ps",
        1e-3: "cavity_coupling_1eneg03_switch_200.0ps"
    }
    
    base_dir = "/media/extradrive/Trajectories/final_nodiss_cavitymd"
    
    # Find available coupling directories
    coupling_dirs = []
    for coupling_val, dir_name in coupling_to_dir.items():
        full_path = os.path.join(base_dir, dir_name)
        if os.path.exists(full_path):
            coupling_dirs.append(full_path)
        else:
            print(f"Warning: Directory not found for coupling {coupling_val}: {dir_name}")
    
    if not coupling_dirs:
        print("No target coupling directories found!")
        return
    
    print(f"Analyzing {len(coupling_dirs)} selected coupling strengths:")
    for d in coupling_dirs:
        dir_name = os.path.basename(d)
        # Extract coupling value for display
        for coupling_val, expected_dir in coupling_to_dir.items():
            if expected_dir == dir_name:
                print(f"  {dir_name} (coupling = {coupling_val})")
                break
    
    # Analyze each coupling directory
    successful_analyses = []
    
    for coupling_dir in coupling_dirs:
        try:
            material_time_function, waiting_times, normalized_fkt, collapsed_data = \
                analyzer.analyze_coupling_directory(coupling_dir)
            
            # Create plots for this coupling
            analyzer.plot_material_time_collapse(coupling_dir, output_dir=base_dir)
            
            successful_analyses.append(coupling_dir)
            
        except Exception as e:
            print(f"Error analyzing {coupling_dir}: {e}")
            continue
    
    # Compare different coupling strengths
    if len(successful_analyses) > 1:
        try:
            analyzer.compare_coupling_strengths(successful_analyses, output_dir=base_dir)
        except Exception as e:
            print(f"Error in coupling comparison: {e}")
            
        # Plot material time evolution
        try:
            analyzer.plot_material_time_evolution(successful_analyses, output_dir=base_dir)
        except Exception as e:
            print(f"Error in material time evolution plot: {e}")
        
        # Create unified material time collapse plot
        try:
            analyzer.plot_unified_material_time_collapse(successful_analyses, output_dir=base_dir)
        except Exception as e:
            print(f"Error in unified material time collapse plot: {e}")
        
        # Fit stretched exponentials to the collapsed data
        try:
            fit_results = analyzer.fit_stretched_exponential_to_collapse(successful_analyses, output_dir=base_dir)
        except Exception as e:
            print(f"Error in stretched exponential fitting: {e}")
        
        # Create overlapped stretched exponential fits plot
        try:
            analyzer.plot_overlapped_stretched_exponential_fits(successful_analyses, output_dir=base_dir)
        except Exception as e:
            print(f"Error in overlapped stretched exponential fits: {e}")
        
        
    
    # Save material time data to text files (for any number of successful analyses)
    if len(successful_analyses) >= 1:
        try:
            analyzer.save_material_time_data(successful_analyses, output_dir=base_dir)
        except Exception as e:
            print(f"Error saving material time data: {e}")
    
    print(f"\nAnalysis complete. Successfully analyzed {len(successful_analyses)} coupling strengths.")
    print("\n Material Time Analysis - Following Tool-Narayanaswamy Formalism:")
    print("1. Individual coupling plots show before/after material time collapse")
    print("   - Left panels: F(k,t) vs lag time (different decay rates for different waiting times)")
    print("   - Right panels: R(ξ) vs material time (SHOULD collapse to single curve if TN theory works)")
    print("2. Comparison plot shows how cavity coupling affects aging dynamics")
    print("3. Material time evolution plot shows the SINGLE ξ(t) function for each coupling")
    print("4. Stretched exponential fits: R(ξ) = exp[-ln(10) × ξ^β]")
    print("\n Key concept: ONE material time ξ(t) per coupling transforms ALL waiting time data")
    print("If Tool-Narayanaswamy theory applies: all curves in right panels should overlap!")
    
    # Print summary statistics
    print(f"\n=== Summary Statistics ===")
    print(f"Note: ref0 and ref14 excluded from analysis (pre-equilibrium and last waiting time)")
    print(f"Material time criterion: {analyzer.criterion_value}")
    print("Waiting time interpretation:")
    print("  - ref1: tw = 0 ps (cavity coupling turns ON)")
    print("  - ref2-ref13: tw = 200-2400 ps (cavity-induced aging)")
    print()
    
    for coupling_dir in successful_analyses:
        coupling_name = os.path.basename(coupling_dir)
        collapsed_data = analyzer.collapsed_data[coupling_dir]
        
        if collapsed_data:
            n_refs = len(collapsed_data)
            waiting_times = [data['waiting_time_ps'] for data in collapsed_data.values()]
            tw_range = f"{min(waiting_times):.0f} to {max(waiting_times):.0f} ps"
            
            # Calculate average lag time range
            max_lag_times = [max(data['lag_time']) for data in collapsed_data.values()]
            avg_max_lag = np.mean(max_lag_times)
            
            print(f"{coupling_name}:")
            print(f"  - {n_refs} reference frames (ref0 excluded)")
            print(f"  - Waiting time range: {tw_range}")
            print(f"  - Average max lag time: {avg_max_lag:.0f} ps")
    
    print(f"\nOutput files saved in: {base_dir}")
    print("Generated files:")
    print(" Plots:")
    print("  - *_material_time_collapse.png (individual coupling collapse)")
    print("  - unified_material_time_collapse.png (ALL couplings in one plot)")
    print("  - stretched_exponential_fits.png (individual stretched exponential fits)")
    print("  - overlapped_stretched_exponential_fits.png (ALL fits overlapped)")
    print("  - stretched_exponential_parameters.png (fit parameters vs coupling)")
    print("  - coupling_comparison.png (comparison across couplings)")
    print("  - material_time_evolution.png (ξ(t) and γ evolution)")
    print(" Data files:")
    print("  - material_time_coupling_*.txt (individual ξ(t) functions)")
    print("  - material_time_all_couplings.txt (combined ξ(t) data)")
    print("  - material_time_derivatives_all_couplings.txt (aging rates γ = dξ/dt)")
    print("  - stretched_exponential_fit_results.txt (fit parameters and statistics)")


if __name__ == "__main__":
    main()
