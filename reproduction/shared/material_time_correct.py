#!/usr/bin/env python3
"""
Correct Material Time Reconstruction from F(k,t) Measurements

This script implements the CORRECT approach to reconstructing material time ξ(t)
from correlation function measurements F(k,t; tw) at regular waiting time intervals.

Physical Setup:
---------------
- Correlation measurements: F(k,t; tw) where tw is the waiting time
- Physical relation: F(k,t; tw) = Φ(ξ(t+tw) - ξ(tw))
- Criterion: When Φ = 0.1, we have ξ(t+tw) - ξ(tw) = 1
- Waiting times tw are at REGULAR LAB TIME intervals: 0, 200, 400, 600, ... ps

Key Insight:
------------
Cannot build ξ(t) sequentially with interpolation! Must solve all constraints
simultaneously using least squares with regularization.

Constraints:
------------
Each crossing where F(k,t; tw) = 0.1 gives: ξ(tcross) - ξ(tw) = 1

This is solved as:
  Minimize: ||A·ξ - b||² + α||D·ξ||²
  
Where:
  - A: constraint matrix (one row per crossing)
  - b: vector of 1's (target values)
  - D: finite difference operator for smoothness
  - α: regularization parameter

Author: Material Time Analysis
Date: November 2, 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.signal import savgol_filter
import os
import glob
import pandas as pd
from pathlib import Path


class MaterialTimeReconstructor:
    """
    Reconstructs material time ξ(t) from F(k,t) correlation measurements
    using the correct least squares approach.
    """
    
    def __init__(self, criterion_value=0.1, alpha_smooth=1.0, n_grid_points=500, verbose=True):
        """
        Initialize the material time reconstructor.
        
        Parameters:
        -----------
        criterion_value : float
            Threshold for Φ (typically 0.1, corresponding to unit material time)
        alpha_smooth : float
            Regularization parameter for smoothness penalty (larger = smoother)
        n_grid_points : int
            Number of regular grid points to solve for (more = smoother at cost of speed)
        verbose : bool
            Print diagnostic information
        """
        self.criterion_value = criterion_value
        self.alpha_smooth = alpha_smooth
        self.n_grid_points = n_grid_points
        self.verbose = verbose
        
    def load_fkt_data(self, coupling_dir):
        """
        Load F(k,t) data from coupling directory.
        
        Parameters:
        -----------
        coupling_dir : str
            Path to coupling directory containing master_fskt_ref*.txt files
            
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
        
        if self.verbose:
            print(f"Looking for F(k,t) files in: {coupling_dir}")
        
        for ref_file in ref_files:
            if "_sample_counts" not in ref_file:
                # Extract ref number from filename
                basename = os.path.basename(ref_file)
                try:
                    ref_num = int(basename.split('ref')[1].split('.')[0])
                    fkt_files[ref_num] = ref_file
                except (IndexError, ValueError):
                    if self.verbose:
                        print(f"  Warning: Could not parse ref number from {basename}")
                    continue
        
        if len(fkt_files) == 0:
            raise ValueError(f"No valid F(k,t) files found in {coupling_dir}")
            
        if self.verbose:
            print(f"  Found {len(fkt_files)} reference frames")
            
        return fkt_files
    
    def extract_normalized_fkt(self, fkt_files):
        """
        Extract and normalize F(k,t) data for material time analysis.
        
        Parameters:
        -----------
        fkt_files : dict
            Dictionary mapping ref_num to file paths
            
        Returns:
        --------
        normalized_fkt : dict
            Dictionary mapping ref_num -> (lag_time, normalized_fkt)
        waiting_times : dict
            Dictionary mapping ref_num -> waiting_time_ps
        """
        normalized_fkt = {}
        waiting_times = {}
        
        if self.verbose:
            print("\nExtracting and normalizing F(k,t) data...")
        
        for ref_num, file_path in fkt_files.items():
            # Skip ref0 (pre-equilibrium) and ref14 (last reference)
            if ref_num == 0:
                if self.verbose:
                    print(f"  ref{ref_num}: SKIPPED (pre-equilibrium)")
                continue
            
            if ref_num == 14:
                if self.verbose:
                    print(f"  ref{ref_num}: SKIPPED (last reference)")
                continue
                
            try:
                # Read the data file
                data = pd.read_csv(file_path, sep='\t', comment='#')
                
                if 'lag_time' not in data.columns or 'fskt' not in data.columns:
                    if self.verbose:
                        print(f"  Warning: Expected columns not found in {file_path}")
                    continue
                    
                lag_time = data['lag_time'].values
                fkt = data['fskt'].values
                
                # Remove any NaN or infinite values
                valid_mask = np.isfinite(fkt) & np.isfinite(lag_time)
                lag_time = lag_time[valid_mask]
                fkt = fkt[valid_mask]
                
                if len(fkt) < 10:
                    if self.verbose:
                        print(f"  Warning: Insufficient data points in ref{ref_num}")
                    continue
                
                # Apply Savitzky-Golay smoothing
                window = min(51, len(fkt) // 4 * 2 + 1)  # Ensure odd window
                if window >= 5:
                    fkt_smoothed = savgol_filter(fkt, window, 2)
                else:
                    fkt_smoothed = fkt
                
                # Normalize: C(t₁,t₂) = [F(k,t) - F(k,∞)] / [F(k,0) - F(k,∞)]
                fkt_inf = np.mean(fkt_smoothed[-50:])  # Average last 50 points as F(k,∞)
                fkt_0 = fkt_smoothed[0]                # F(k,0)
                
                if abs(fkt_0 - fkt_inf) > 1e-10:  # Avoid division by zero
                    normalized = (fkt_smoothed - fkt_inf) / (fkt_0 - fkt_inf)
                else:
                    if self.verbose:
                        print(f"  Warning: No decay observed in ref{ref_num}")
                    normalized = np.ones_like(fkt_smoothed)
                
                # Enforce monotonic decay
                fkt_mono = np.copy(normalized)
                for i in range(1, len(fkt_mono)):
                    if fkt_mono[i] > fkt_mono[i-1]:
                        fkt_mono[i] = fkt_mono[i-1]
                
                normalized_fkt[ref_num] = (lag_time, fkt_mono)
                
                # Calculate waiting time: tw_ps = (ref_num - 1) * 200
                # ref1 = 0 ps (coupling turns on), ref2 = 200 ps, etc.
                tw_ps = (ref_num - 1) * 200
                waiting_times[ref_num] = tw_ps
                
                if self.verbose:
                    print(f"  ref{ref_num}: tw = {tw_ps:.1f} ps, {len(lag_time)} points, "
                          f"F(k,0)={fkt_0:.1f}, F(k,∞)={fkt_inf:.1f}, decay={(fkt_0-fkt_inf)/fkt_0*100:.1f}%")
                    
            except Exception as e:
                if self.verbose:
                    print(f"  Error processing ref{ref_num}: {e}")
                continue
        
        return normalized_fkt, waiting_times
    
    def find_crossings(self, normalized_fkt, waiting_times):
        """
        Find where F(k,t; tw) crosses the criterion value.
        
        Parameters:
        -----------
        normalized_fkt : dict
            Normalized F(k,t) data for each reference
        waiting_times : dict
            Waiting times for each reference
            
        Returns:
        --------
        crossings : list of tuples
            List of (tw, tcross) where:
            - tw: waiting time (reference time)
            - tcross: lab time where F crosses criterion (tcross = tw + lag_cross)
        """
        crossings = []
        
        if self.verbose:
            print(f"\nFinding crossings at criterion = {self.criterion_value}")
        
        for ref_num in sorted(normalized_fkt.keys()):
            if ref_num == 14:  # Skip last reference
                continue
                
            lag_time, fkt = normalized_fkt[ref_num]
            tw = waiting_times[ref_num]
            
            # Find crossing point
            for i in range(1, len(fkt)):
                if (fkt[i-1] > self.criterion_value >= fkt[i] or 
                    fkt[i-1] >= self.criterion_value > fkt[i]):
                    
                    # Linear interpolation for exact crossing
                    if fkt[i-1] != fkt[i]:
                        lag_cross = lag_time[i-1] + (lag_time[i] - lag_time[i-1]) * \
                                   (self.criterion_value - fkt[i-1]) / (fkt[i] - fkt[i-1])
                    else:
                        lag_cross = lag_time[i]
                    
                    # tcross = tw + lag_cross (in laboratory time)
                    tcross = tw + lag_cross
                    
                    crossings.append((tw, tcross))
                    
                    if self.verbose:
                        print(f"  ref{ref_num}: tw = {tw:.1f}, lag = {lag_cross:.1f}, "
                              f"tcross = {tcross:.1f} ps")
                    
                    break
        
        if self.verbose:
            print(f"\nTotal crossings found: {len(crossings)}")
        
        return crossings
    
    def build_constraint_system(self, crossings):
        """
        Build the least squares constraint system from crossing data.
        
        Uses a DENSE REGULAR GRID for smooth reconstruction.
        Given constraints: ξ(tcross) - ξ(tw) = 1
        
        Parameters:
        -----------
        crossings : list of tuples
            List of (tw, tcross) pairs
            
        Returns:
        --------
        time_points : array
            Dense regular time grid
        A : array
            Constraint matrix (n_constraints x n_timepoints)
        b : array
            Target values (all 1's for unit material time differences)
        """
        # Create a DENSE REGULAR GRID instead of sparse constraint points
        # This is key for smooth reconstruction!
        t_max = max([max(tw, tcross) for tw, tcross in crossings])
        time_points = np.linspace(0, t_max, self.n_grid_points)
        n_times = len(time_points)
        
        if self.verbose:
            print(f"\nBuilding constraint system:")
            print(f"  Using DENSE REGULAR GRID for smooth reconstruction")
            print(f"  Grid points: {n_times}")
            print(f"  Time range: [{time_points[0]:.1f}, {time_points[-1]:.1f}] ps")
            print(f"  Grid spacing: {time_points[1] - time_points[0]:.2f} ps")
        
        # Build constraint matrix
        n_constraints = len(crossings) + 1  # +1 for ξ(0) = 0 boundary condition
        A = np.zeros((n_constraints, n_times))
        b = np.zeros(n_constraints)
        
        # Add crossing constraints: ξ(tcross) - ξ(tw) = 1
        # Use linear interpolation to map constraint times to grid indices
        for i, (tw, tcross) in enumerate(crossings):
            # For tw: find surrounding grid points and interpolation weights
            idx_tw_low = np.searchsorted(time_points, tw, side='right') - 1
            idx_tw_low = max(0, min(idx_tw_low, n_times - 2))
            idx_tw_high = idx_tw_low + 1
            
            if idx_tw_high < n_times:
                # Linear interpolation weights
                w_tw = (tw - time_points[idx_tw_low]) / (time_points[idx_tw_high] - time_points[idx_tw_low])
                A[i, idx_tw_low] -= (1 - w_tw)
                A[i, idx_tw_high] -= w_tw
            else:
                A[i, idx_tw_low] -= 1.0
            
            # For tcross: find surrounding grid points and interpolation weights
            idx_tc_low = np.searchsorted(time_points, tcross, side='right') - 1
            idx_tc_low = max(0, min(idx_tc_low, n_times - 2))
            idx_tc_high = idx_tc_low + 1
            
            if idx_tc_high < n_times:
                w_tc = (tcross - time_points[idx_tc_low]) / (time_points[idx_tc_high] - time_points[idx_tc_low])
                A[i, idx_tc_low] += (1 - w_tc)
                A[i, idx_tc_high] += w_tc
            else:
                A[i, idx_tc_low] += 1.0
            
            b[i] = 1.0
        
        # Add boundary condition: ξ(0) = 0
        A[-1, 0] = 1.0
        b[-1] = 0.0
        
        if self.verbose:
            print(f"  Constraint matrix shape: {A.shape}")
            print(f"  Matrix rank: {np.linalg.matrix_rank(A)}")
        
        return time_points, A, b
    
    def add_smoothness_regularization(self, n_times):
        """
        Create finite difference operator for smoothness regularization.
        
        Penalizes the second derivative: d²ξ/dt² ≈ ξ[i] - 2ξ[i+1] + ξ[i+2]
        
        Parameters:
        -----------
        n_times : int
            Number of time points
            
        Returns:
        --------
        D : array
            Finite difference operator matrix
        """
        D = np.zeros((n_times - 2, n_times))
        for i in range(n_times - 2):
            D[i, i] = 1.0
            D[i, i+1] = -2.0
            D[i, i+2] = 1.0
        
        return D
    
    def solve_least_squares(self, time_points, A, b):
        """
        Solve the regularized least squares problem.
        
        Minimize: ||A·ξ - b||² + α||D·ξ||²
        
        Parameters:
        -----------
        time_points : array
            Time points where ξ is evaluated
        A : array
            Constraint matrix
        b : array
            Target values
            
        Returns:
        --------
        xi_values : array
            Material time values at each time point
        """
        n_times = len(time_points)
        
        if self.alpha_smooth > 0:
            # Add smoothness regularization
            D = self.add_smoothness_regularization(n_times)
            
            # Combined system
            A_reg = np.vstack([A, self.alpha_smooth * D])
            b_reg = np.hstack([b, np.zeros(n_times - 2)])
            
            if self.verbose:
                print(f"\nSolving regularized least squares (α = {self.alpha_smooth})...")
            
            # Solve
            xi_values, residuals, rank, s = np.linalg.lstsq(A_reg, b_reg, rcond=None)
        else:
            # Pure least squares without regularization
            if self.verbose:
                print(f"\nSolving PURE least squares (NO regularization, α = 0)...")
            
            # Solve
            xi_values, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
        
        # Check solution quality
        constraint_errors = A @ xi_values - b
        
        if self.verbose:
            print(f"  Solution range: ξ ∈ [{xi_values.min():.3f}, {xi_values.max():.3f}]")
            print(f"  Max constraint error: {np.max(np.abs(constraint_errors)):.6f}")
            print(f"  RMS constraint error: {np.sqrt(np.mean(constraint_errors**2)):.6f}")
            print(f"  Residual norm: {np.linalg.norm(residuals) if len(residuals) > 0 else 'N/A'}")
            
            # Check for non-monotonicity
            dxi = np.diff(xi_values)
            n_decreasing = np.sum(dxi < 0)
            if n_decreasing > 0:
                print(f"  WARNING: Solution has {n_decreasing} non-monotonic points!")
        
        return xi_values
    
    def interpolate_to_fine_grid(self, time_points, xi_values, t_max_extend=1.05):
        """
        Interpolate material time to a fine regular grid.
        
        Parameters:
        -----------
        time_points : array
            Discrete time points from constraint system
        xi_values : array
            Material time values at those points
        t_max_extend : float
            Factor to extend time grid beyond last point (small extension only)
            
        Returns:
        --------
        time_grid : array
            Fine regular time grid
        xi_grid : array
            Interpolated material time values
        """
        if self.verbose:
            print(f"\nDEBUG: Interpolation diagnostics")
            print(f"  Input points: {len(time_points)}")
            print(f"  Time range: [{time_points[0]:.1f}, {time_points[-1]:.1f}]")
            print(f"  ξ range: [{xi_values.min():.3f}, {xi_values.max():.3f}]")
            print(f"  First 5 points: t={time_points[:5]}, ξ={xi_values[:5]}")
            print(f"  Last 5 points: t={time_points[-5:]}, ξ={xi_values[-5:]}")
        
        # Create fine grid with minimal extrapolation
        t_max = time_points[-1] * t_max_extend
        time_grid = np.linspace(0, t_max, 2000)
        
        # Use linear interpolation for stability (cubic splines can oscillate)
        # This is more robust for material time reconstruction
        interp_func = interp1d(time_points, xi_values, kind='linear', 
                               bounds_error=False, fill_value='extrapolate')
        xi_grid = interp_func(time_grid)
        
        if self.verbose:
            print(f"\nInterpolated to fine grid:")
            print(f"  Grid points: {len(time_grid)}")
            print(f"  Time range: [0, {t_max:.1f}] ps")
            print(f"  ξ range: [{xi_grid.min():.3f}, {xi_grid.max():.3f}]")
            
            # Check for monotonicity
            is_monotonic = np.all(np.diff(xi_grid) >= 0)
            print(f"  Monotonic: {is_monotonic}")
            if not is_monotonic:
                n_violations = np.sum(np.diff(xi_grid) < 0)
                print(f"  WARNING: {n_violations} monotonicity violations!")
        
        return time_grid, xi_grid
    
    def reconstruct_material_time(self, coupling_dir):
        """
        Complete pipeline: load data, find crossings, solve for ξ(t).
        
        Parameters:
        -----------
        coupling_dir : str
            Path to coupling directory
            
        Returns:
        --------
        results : dict
            Dictionary containing:
            - 'time_grid': fine time grid
            - 'xi_grid': material time on fine grid
            - 'time_points': discrete constraint points
            - 'xi_values': material time at constraint points
            - 'crossings': list of (tw, tcross) crossing data
            - 'normalized_fkt': normalized correlation data
            - 'waiting_times': waiting times for each ref
        """
        if self.verbose:
            print(f"\n{'='*70}")
            print(f"MATERIAL TIME RECONSTRUCTION: {os.path.basename(coupling_dir)}")
            print(f"{'='*70}")
        
        # Load and process data
        fkt_files = self.load_fkt_data(coupling_dir)
        normalized_fkt, waiting_times = self.extract_normalized_fkt(fkt_files)
        
        # Find crossings
        crossings = self.find_crossings(normalized_fkt, waiting_times)
        
        if len(crossings) == 0:
            raise ValueError("No crossings found! Cannot reconstruct material time.")
        
        # Build and solve constraint system
        time_points, A, b = self.build_constraint_system(crossings)
        xi_values = self.solve_least_squares(time_points, A, b)
        
        # Interpolate to fine grid
        time_grid, xi_grid = self.interpolate_to_fine_grid(time_points, xi_values)
        
        results = {
            'time_grid': time_grid,
            'xi_grid': xi_grid,
            'time_points': time_points,
            'xi_values': xi_values,
            'crossings': crossings,
            'normalized_fkt': normalized_fkt,
            'waiting_times': waiting_times,
            'coupling_dir': coupling_dir,
            'alpha_smooth': self.alpha_smooth
        }
        
        return results
    
    def calculate_ageing_rate(self, time_grid, xi_grid):
        """
        Calculate ageing rate γ = dξ/dt.
        
        Parameters:
        -----------
        time_grid : array
            Laboratory time points
        xi_grid : array
            Material time values
            
        Returns:
        --------
        gamma : array
            Ageing rate at each time point
        """
        # Use finite differences for derivative
        dt = time_grid[1] - time_grid[0]
        gamma = np.gradient(xi_grid, dt)
        
        return gamma


def plot_material_time_evolution(results_dict, output_dir='.'):
    """
    Plot material time evolution for multiple coupling strengths.
    
    Parameters:
    -----------
    results_dict : dict
        Dictionary mapping coupling_dir -> reconstruction results
    output_dir : str
        Output directory for plots
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Coupling labels for plotting
    coupling_labels = {
        "cavity_coupling_0epos00_switch_200.0ps": "0.0",
        "cavity_coupling_3eneg04_switch_200.0ps": "3×10⁻⁴",
        "cavity_coupling_5eneg04_switch_200.0ps": "5×10⁻⁴",
        "cavity_coupling_7eneg04_switch_200.0ps": "7×10⁻⁴",
        "cavity_coupling_1eneg03_switch_200.0ps": "1×10⁻³"
    }
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(results_dict)))
    
    for i, (coupling_dir, results) in enumerate(results_dict.items()):
        coupling_name = os.path.basename(coupling_dir)
        label = coupling_labels.get(coupling_name, coupling_name)
        color = colors[i]
        
        time_grid = results['time_grid']
        xi_grid = results['xi_grid']
        
        # Plot 1: Material time ξ(t)
        ax1.plot(time_grid, xi_grid, color=color, linewidth=2, alpha=0.8,
                label=f'g = {label}')
        
        # Plot constraint points
        time_points = results['time_points']
        xi_values = results['xi_values']
        ax1.scatter(time_points, xi_values, color=color, s=30, alpha=0.5, zorder=5)
        
        # Plot 2: Ageing rate γ = dξ/dt
        gamma = np.gradient(xi_grid, time_grid[1] - time_grid[0])
        ax2.plot(time_grid, gamma, color=color, linewidth=2, alpha=0.8,
                label=f'g = {label}')
    
    # Format plot 1
    ax1.set_xlabel('Laboratory time t (ps)', fontsize=12)
    ax1.set_ylabel('Material time ξ(t)', fontsize=12)
    ax1.set_title('Material Time Evolution (Corrected Reconstruction)\n' +
                  'Dots: constraint points from F(k,t) crossings', fontsize=13)
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    
    # Format plot 2
    ax2.set_xlabel('Laboratory time t (ps)', fontsize=12)
    ax2.set_ylabel('Ageing rate γ = dξ/dt', fontsize=12)
    ax2.set_title('Ageing Rate Evolution\n' +
                  'Derived from material time function', fontsize=13)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save with appropriate filename based on regularization
    alpha_val = list(results_dict.values())[0]['alpha_smooth'] if results_dict else 0.0
    if alpha_val == 0.0:
        filename = 'material_time_evolution_no_regularization.png'
    else:
        filename = f'material_time_evolution_alpha{alpha_val:.2f}.png'
    
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {output_path}")
    
    plt.show()
    return fig


# def compare_methods(coupling_dir, output_dir='.'):
#     """
#     Compare the correct method vs the old interpolation method.
    
#     Parameters:
#     -----------
#     coupling_dir : str
#         Path to coupling directory
#     output_dir : str
#         Output directory for comparison plot
#     """
#     reconstructor = MaterialTimeReconstructor(verbose=True, alpha_smooth=1e-8)
    
#     print("\n" + "="*70)
#     print("METHOD COMPARISON: Correct vs Old Interpolation")
#     print("="*70)
    
#     results = reconstructor.reconstruct_material_time(coupling_dir)
    
#     # Check constraint satisfaction
#     crossings = results['crossings']
#     time_grid = results['time_grid']
#     xi_grid = results['xi_grid']
    
#     # Interpolator for checking constraints
#     xi_interp = interp1d(time_grid, xi_grid, kind='cubic', fill_value='extrapolate')
    
#     print(f"\nConstraint Satisfaction Check:")
#     print(f"{'Constraint':<30} {'ξ(tcross) - ξ(tw)':<20} {'Target':<10} {'Error':<10}")
#     print("-" * 70)
    
#     errors = []
#     for tw, tcross in crossings[:10]:  # Show first 10
#         xi_tw = xi_interp(tw)
#         xi_tcross = xi_interp(tcross)
#         delta_xi = xi_tcross - xi_tw
#         error = abs(delta_xi - 1.0)
#         errors.append(error)
#         print(f"ξ({tcross:.1f}) - ξ({tw:.1f}): {delta_xi:18.4f} {1.0:<10.1f} {error:<10.6f}")
    
#     if len(crossings) > 10:
#         print(f"... and {len(crossings) - 10} more constraints")
    
#     # Calculate RMS error
#     all_errors = []
#     all_deltas = []
#     for tw, tcross in crossings:
#         xi_tw = xi_interp(tw)
#         xi_tcross = xi_interp(tcross)
#         delta_xi = xi_tcross - xi_tw
#         error = delta_xi - 1.0
#         all_errors.append(error**2)
#         all_deltas.append(delta_xi)
    
#     rms_error = np.sqrt(np.mean(all_errors))
#     max_error = np.max(np.abs([np.sqrt(e) for e in all_errors]))
    
#     print(f"\n" + "="*70)
#     print("FINAL CONSTRAINT SATISFACTION:")
#     print("="*70)
#     print(f"RMS constraint error: {rms_error:.6f}")
#     print(f"Max constraint error: {max_error:.6f}")
#     print(f"Mean Δξ: {np.mean(all_deltas):.6f} (target: 1.0)")
#     print(f"Std Δξ: {np.std(all_deltas):.6f}")
#     print(f"Min Δξ: {np.min(all_deltas):.6f}")
#     print(f"Max Δξ: {np.max(all_deltas):.6f}")
#     print("="*70)
    
#     return results


if __name__ == "__main__":
    # Example usage
    base_dir = Path("/home/mh7373/GitRepos/cav-hoomd/examples/final_nodiss_cavitymd")
    
    coupling_dirs = [
        base_dir / "cavity_coupling_0epos00_switch_200.0ps",
        base_dir / "cavity_coupling_3eneg04_switch_200.0ps",
        base_dir / "cavity_coupling_5eneg04_switch_200.0ps",
        base_dir / "cavity_coupling_7eneg04_switch_200.0ps",
        base_dir / "cavity_coupling_1eneg03_switch_200.0ps"
    ]
    
    # Reconstruct material time for all coupling strengths
    # Use small regularization to balance constraint satisfaction and smoothness
    reconstructor = MaterialTimeReconstructor(
        criterion_value=0.1, 
        alpha_smooth=7e-11,  # Regularization: larger = smoother (try 0.01, 0.1, 1.0, 10.0)
        n_grid_points=500,  # Dense grid for smooth interpolation
        verbose=True
    )
    
    results_dict = {}
    for coupling_dir in coupling_dirs:
        if not coupling_dir.exists():
            print(f"Warning: {coupling_dir} does not exist, skipping")
            continue
        
        try:
            results = reconstructor.reconstruct_material_time(str(coupling_dir))
            results_dict[str(coupling_dir)] = results
        except Exception as e:
            print(f"Error processing {coupling_dir.name}: {e}")
            continue
    
    # Plot evolution for all couplings
    if results_dict:
        plot_material_time_evolution(results_dict, output_dir=str(base_dir))
        
        # Also do detailed comparison for one coupling
        first_coupling = list(results_dict.keys())[0]
        print("\n" + "="*70)
        print("Detailed analysis for first coupling:")
        #compare_methods(first_coupling, output_dir=str(base_dir))

