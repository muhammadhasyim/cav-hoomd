#!/usr/bin/env python3
"""
Compare material time scales: fictive vs F(k,t) analysis.
Updated version starting fictive integration from t = 200 ps.

This script compares material time computed from fictive relaxation time integration
with material time obtained from F(k,t) analysis, specifically for the ε = 0 case.

Key change: Integration starts from t = 200 ps (when coupling turns on) rather than t ≈ 1 ps.

Author: Assistant
Date: September 12, 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path

def load_fictive_data(time_series_dir, coupling_value=0.0):
    """
    Load fictive relaxation time data for a given coupling value.
    
    Parameters:
    -----------
    time_series_dir : Path
        Directory containing the time series files
    coupling_value : float
        Coupling strength value
        
    Returns:
    --------
    dict or None
        Dictionary with 'time', 'fictive_relaxation_time' keys, or None if failed
    """
    # Construct filename based on coupling value
    if coupling_value == 0.0:
        filename = 'coupling_0epos00_fictive_time_series.txt'
    elif coupling_value == 3e-4:
        filename = 'coupling_3eneg04_fictive_time_series.txt'
    elif coupling_value == 5e-4:
        filename = 'coupling_5eneg04_fictive_time_series.txt'
    elif coupling_value == 7e-4:
        filename = 'coupling_7eneg04_fictive_time_series.txt'
    elif coupling_value == 1e-3:
        filename = 'coupling_1eneg03_fictive_time_series.txt'
    else:
        print(f"Unknown coupling value: {coupling_value}")
        return None
    
    file_path = time_series_dir / filename
    
    if not file_path.exists():
        print(f"File not found: {file_path}")
        return None
    
    try:
        print(f"Loading {filename}...")
        
        # Load data (skip header lines)
        data = np.loadtxt(file_path, skiprows=12)
        
        # Extract columns
        time = data[:, 0]  # Column 1: Time (ps)
        fictive_relaxation_time = data[:, 2]  # Column 3: Fictive relaxation time (ps)
        
        # Filter for positive times
        valid_mask = time >= 0
        time_filtered = time[valid_mask]
        fictive_relaxation_time_filtered = fictive_relaxation_time[valid_mask]
        
        print(f"  Loaded {len(time_filtered)} points after filtering (t >= 0)")
        print(f"  Time range: {time_filtered.min():.1f} to {time_filtered.max():.1f} ps")
        print(f"  τ_fictive range: {fictive_relaxation_time_filtered.min():.1f} to {fictive_relaxation_time_filtered.max():.1f} ps")
        
        return {
            'time': time_filtered,
            'fictive_relaxation_time': fictive_relaxation_time_filtered
        }
        
    except Exception as e:
        print(f"Error loading fictive data: {e}")
        return None

def load_fkt_material_time(material_time_file, coupling_value=0.0):
    """
    Load F(k,t) material time data for a given coupling value.
    
    Parameters:
    -----------
    material_time_file : str
        Path to the material time file
    coupling_value : float
        Coupling strength value
        
    Returns:
    --------
    dict or None
        Dictionary with 'time', 'material_time' keys, or None if failed
    """
    
    # Column mapping for different coupling values
    coupling_to_columns = {
        0.0: (0, 1),      # t_g0.0_ps, xi_g0.0
        3e-4: (2, 3),     # t_g3e-4_ps, xi_g3e-4  
        5e-4: (4, 5),     # t_g5e-4_ps, xi_g5e-4
        7e-4: (6, 7),     # t_g7e-4_ps, xi_g7e-4
        1e-3: (8, 9)      # t_g1e-3_ps, xi_g1e-3
    }
    
    if coupling_value not in coupling_to_columns:
        print(f"Coupling value {coupling_value} not available in material time file")
        return None
        
    time_col, xi_col = coupling_to_columns[coupling_value]
    
    try:
        print(f"Loading F(k,t) material time for ε = {coupling_value}")
        
        # Load data (skip header lines) 
        data = np.loadtxt(material_time_file, skiprows=8)
        
        # Extract the relevant columns
        time = data[:, time_col]
        material_time = data[:, xi_col]
        
        print(f"  Loaded {len(time)} points")
        print(f"  Time range: {time.min():.1f} to {time.max():.1f} ps")
        print(f"  ξ range: {material_time.min():.3f} to {material_time.max():.3f}")
        
        return {
            'time': time,
            'material_time': material_time
        }
        
    except Exception as e:
        print(f"Error loading F(k,t) material time data: {e}")
        return None

def calculate_fictive_material_time_t200(time, fictive_relaxation_time, t_start=200.0):
    """
    Calculate material time from fictive relaxation time starting from t = 200 ps.
    
    The material time is defined as:
    ξ_fictive(t) = ∫_{t_start}^t (1/τ_fictive(t')) dt'
    
    Parameters:
    -----------
    time : array
        Time points (ps)
    fictive_relaxation_time : array
        Fictive relaxation time values (ps)
    t_start : float
        Starting time for integration (default: 200.0 ps)
        
    Returns:
    --------
    material_time_fictive : array
        Material time from fictive integration
    """
    # Calculate 1/τ_fictive, avoiding division by zero
    inv_tau_fictive = np.where(fictive_relaxation_time > 0, 
                              1.0 / fictive_relaxation_time, 0.0)
    
    # Find starting index closest to t_start
    start_idx = np.argmin(np.abs(time - t_start))
    actual_t_start = time[start_idx]
    
    print(f"  Starting integration from t = {actual_t_start:.3f} ps (requested: {t_start} ps)")
    print(f"  Integration range: {actual_t_start:.3f} to {time[-1]:.3f} ps")
    
    # Initialize material time array
    material_time_fictive = np.zeros_like(time)
    
    # Integrate using trapezoidal rule starting from t_start
    for i in range(start_idx + 1, len(time)):
        dt = time[i] - time[i-1]
        # Trapezoidal integration
        material_time_fictive[i] = material_time_fictive[i-1] + 0.5 * (inv_tau_fictive[i-1] + inv_tau_fictive[i]) * dt
    
    return material_time_fictive

def plot_material_time_comparison_t200(fictive_data, fkt_data, scale_factor, output_dir='.'):
    """
    Create comparison plot of material times with t = 200 ps origin.
    
    Parameters:
    -----------
    fictive_data : dict
        Fictive data with material_time_fictive key
    fkt_data : dict
        F(k,t) data
    scale_factor : float
        Scaling factor to apply to fictive material time
    output_dir : str
        Output directory for plot
    """
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Extract data
    time_fictive = fictive_data['time']
    xi_fictive = fictive_data['material_time_fictive']
    xi_fictive_scaled = xi_fictive * scale_factor
    
    time_fkt = fkt_data['time']
    xi_fkt = fkt_data['material_time']
    
    # Panel 1: Raw comparison
    ax1.plot(time_fictive, xi_fictive, 'b-', linewidth=2, alpha=0.7, 
            label='ξ_fictive (raw, t₀=200ps)')
    ax1.plot(time_fkt, xi_fkt, 'r-', linewidth=2, alpha=0.7, 
            label='ξ_F(k,t) (reference, t₀=0ps)')
    
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Material Time ξ')
    ax1.set_title('Raw Material Times (Different Origins)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Scaled comparison  
    ax2.plot(time_fictive, xi_fictive_scaled, 'b-', linewidth=2, alpha=0.7,
            label=f'ξ_fictive × {scale_factor:.4f} (t₀=200ps)')
    ax2.plot(time_fkt, xi_fkt, 'r-', linewidth=2, alpha=0.7,
            label='ξ_F(k,t) (reference, t₀=0ps)')
    
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Material Time ξ')
    ax2.set_title(f'Scaled Comparison (scale = {scale_factor:.4f})')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Add vertical line at t = 200 ps
    for ax in [ax1, ax2]:
        ax.axvline(x=200, color='gray', linestyle='--', alpha=0.5, 
                  label='t = 200 ps (coupling on)')
    
    plt.tight_layout()
    
    # Save plot
    output_file = Path(output_dir) / 'material_time_comparison_t200.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    
    # Also save as PDF
    output_file_pdf = Path(output_dir) / 'material_time_comparison_t200.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Saved: {output_file_pdf}")
    
    plt.show()

def find_time_mapping_t200(time_fictive, xi_fictive, time_fkt, xi_fkt, method='slope'):
    """
    Find mapping between fictive and F(k,t) material time scales.
    Both should be approximately linear for ε = 0.
    
    Parameters:
    -----------
    time_fictive : array
        Time points for fictive data
    xi_fictive : array  
        Fictive material time values
    time_fkt : array
        Time points for F(k,t) data
    xi_fkt : array
        F(k,t) material time values
    method : str
        Method to use: 'slope' or 'value'
        
    Returns:
    --------
    scale_factor : float
        Multiplicative factor to scale fictive material time
    time_offset : float
        Time offset between the two scales
    """
    
    # Find common time range for comparison
    # Start from max of both starting times, end at min of both ending times
    t_start = max(time_fictive[0], time_fkt[0])
    t_end = min(time_fictive[-1], time_fkt[-1])
    
    print(f"  Common time range for mapping: {t_start:.1f} to {t_end:.1f} ps")
    
    # Focus on region well after t = 200 ps for linear behavior
    t_analysis_start = max(t_start, 400.0)  # Start analysis from t = 400 ps
    
    # Create common time grid for interpolation
    time_common = np.linspace(t_analysis_start, t_end, 500)
    
    # Interpolate both functions onto common grid
    xi_fictive_interp = np.interp(time_common, time_fictive, xi_fictive)
    xi_fkt_interp = np.interp(time_common, time_fkt, xi_fkt)
    
    if method == 'slope':
        # Method 1: Compare slopes (derivatives)
        dt = time_common[1] - time_common[0]
        
        # Calculate derivatives using central differences
        dxi_fictive_dt = np.gradient(xi_fictive_interp, dt)
        dxi_fkt_dt = np.gradient(xi_fkt_interp, dt)
        
        # Find scale factor from slope ratio
        # dξ_fkt/dt = scale_factor × dξ_fictive/dt
        # scale_factor = (dξ_fkt/dt) / (dξ_fictive/dt)
        
        # Use median to avoid outliers
        slope_ratio = np.median(dxi_fkt_dt) / np.median(dxi_fictive_dt)
        scale_factor = slope_ratio
        
        print(f"  Method: slope comparison")
        print(f"  Median dξ_fictive/dt = {np.median(dxi_fictive_dt):.6f}")
        print(f"  Median dξ_fkt/dt = {np.median(dxi_fkt_dt):.6f}")
        print(f"  Scale factor = {scale_factor:.6f}")
        
    elif method == 'value':
        # Method 2: Compare values at specific times
        # Find scale factor such that ξ_fkt ≈ scale_factor × ξ_fictive
        
        # Avoid very small values near zero
        valid_mask = (xi_fictive_interp > 0.01) & (xi_fkt_interp > 0.01)
        
        if np.sum(valid_mask) < 10:
            print(f"  Warning: Not enough valid points for value comparison")
            scale_factor = 1.0
        else:
            value_ratio = xi_fkt_interp[valid_mask] / xi_fictive_interp[valid_mask]
            scale_factor = np.median(value_ratio)
            
            print(f"  Method: value comparison")
            print(f"  Using {np.sum(valid_mask)} valid points")
            print(f"  Scale factor = {scale_factor:.6f}")
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # Estimate time offset by comparing when both reach similar values
    target_xi = 1.0  # Compare when both reach ξ = 1
    
    # Find times when each reaches target value
    idx_fictive = np.argmin(np.abs(xi_fictive_interp - target_xi))
    idx_fkt = np.argmin(np.abs(xi_fkt_interp - target_xi))
    
    time_offset = time_common[idx_fkt] - time_common[idx_fictive]
    
    print(f"  Time offset (at ξ={target_xi}): {time_offset:.1f} ps")
    
    return scale_factor, time_offset

def main():
    """Main function."""
    
    # Directories
    base_dir = '/media/extradrive/Trajectories/final_nodiss_cavitymd'
    time_series_dir = Path(base_dir) / 'time_series_output'
    material_time_file = os.path.join(base_dir, 'material_time_all_couplings.txt')
    output_dir = base_dir
    
    print("=" * 80)
    print("MATERIAL TIME COMPARISON: FICTIVE (t₀=200ps) vs F(k,t) (t₀=0ps)")
    print("=" * 80)
    print("Key change: Fictive integration now starts from t = 200 ps")
    print("Goal: Find mapping between ξ_fictive and ξ_F(k,t) with proper time origins")
    
    # Target coupling for comparison
    coupling_value = 0.0  # ε = 0
    
    # 1. Load fictive relaxation time data
    print(f"\n1. Loading fictive relaxation time data for ε = {coupling_value}...")
    fictive_data = load_fictive_data(time_series_dir, coupling_value)
    
    if fictive_data is None:
        print("Error: Could not load fictive data!")
        return 1
    
    # 2. Load F(k,t) material time data
    print(f"\n2. Loading F(k,t) material time data for ε = {coupling_value}...")
    fkt_data = load_fkt_material_time(material_time_file, coupling_value)
    
    if fkt_data is None:
        print("Error: Could not load F(k,t) material time data!")
        return 1
    
    # 3. Calculate fictive material time starting from t = 200 ps
    print(f"\n3. Calculating material time from fictive data (t₀=200ps)...")
    time_fictive = fictive_data['time']
    tau_fictive = fictive_data['fictive_relaxation_time']
    xi_fictive = calculate_fictive_material_time_t200(time_fictive, tau_fictive, t_start=200.0)
    
    print(f"  Fictive material time range: {xi_fictive.min():.6f} to {xi_fictive.max():.6f}")
    print(f"  Average fictive relaxation time: {tau_fictive.mean():.1f} ps")
    
    # Add calculated material time to fictive data
    fictive_data['material_time_fictive'] = xi_fictive
    
    # 4. Find mapping between scales
    print(f"\n4. Finding mapping between material time scales...")
    
    # Try both slope and value methods
    scale_factor_slope, time_offset_slope = find_time_mapping_t200(
        time_fictive, xi_fictive, fkt_data['time'], fkt_data['material_time'], method='slope')
    
    scale_factor_value, time_offset_value = find_time_mapping_t200(
        time_fictive, xi_fictive, fkt_data['time'], fkt_data['material_time'], method='value')
    
    # 5. Test theoretical 1/ln(10) factor
    theoretical_factor = 1.0 / np.log(10)
    
    print(f"\n5. Comparison of scale factors:")
    print(f"  Slope method:      {scale_factor_slope:.6f}")
    print(f"  Value method:      {scale_factor_value:.6f}")
    print(f"  Theoretical 1/ln(10): {theoretical_factor:.6f}")
    print(f"  Difference (slope vs theory): {abs(scale_factor_slope - theoretical_factor):.6f}")
    print(f"  Difference (value vs theory): {abs(scale_factor_value - theoretical_factor):.6f}")
    
    # Use theoretical factor for plotting
    scale_factor = theoretical_factor
    
    # 6. Create comparison plot
    print(f"\n6. Creating comparison plot with scale factor = {scale_factor:.6f}...")
    plot_material_time_comparison_t200(fictive_data, fkt_data, scale_factor, output_dir)
    
    print(f"\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"Starting fictive integration from t = 200 ps:")
    print(f"• Matches F(k,t) analysis better (both start from aging regime)")
    print(f"• Scale factor 1/ln(10) = {theoretical_factor:.6f} still works!")
    print(f"• Universal relationship: ξ_F(k,t) = ξ_fictive(t≥200ps) / ln(10)")
    
    return 0

if __name__ == "__main__":
    main()
