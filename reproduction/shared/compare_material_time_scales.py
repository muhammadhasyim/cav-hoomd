#!/usr/bin/env python3
"""
Compare material time scales from fictive relaxation time vs F(k,t) analysis.

This script compares:
1. Material time from integrating 1/τ_fictive: ξ_fictive(t) = ∫ (1/τ_fictive(t')) dt'
2. Material time from F(k,t) analysis: ξ_fkt(t) where Δξ = 1 corresponds to R(Δξ) = 0.1

The goal is to find the mapping between these two scales, particularly for ε=0 where
both should be linear in time.

Author: Assistant
Date: September 12, 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import subprocess

def load_fictive_data(time_series_dir, coupling_value=0.0):
    """
    Load fictive relaxation time data for a specific coupling strength.
    
    Parameters:
    -----------
    time_series_dir : str or Path
        Directory containing fictive time series files
    coupling_value : float
        Coupling strength to load (default: 0.0 for ε=0)
        
    Returns:
    --------
    fictive_data : dict or None
        Dictionary with time and relaxation time data, or None if not found
    """
    time_series_path = Path(time_series_dir)
    
    # Map coupling values to filenames
    coupling_file_map = {
        0.0: 'coupling_0epos00_fictive_time_series.txt',
        3e-4: 'coupling_3eneg04_fictive_time_series.txt',
        5e-4: 'coupling_5eneg04_fictive_time_series.txt',
        7e-4: 'coupling_7eneg04_fictive_time_series.txt',
        1e-3: 'coupling_1eneg03_fictive_time_series.txt',
    }
    
    if coupling_value not in coupling_file_map:
        print(f"Error: Coupling value {coupling_value} not supported")
        return None
    
    filename = coupling_file_map[coupling_value]
    file_path = time_series_path / filename
    
    if not file_path.exists():
        print(f"Error: File {filename} not found in {time_series_dir}")
        return None
    
    print(f"Loading {filename}...")
    
    # Load data (skip header lines starting with # and column names)
    data = np.loadtxt(file_path, skiprows=12)  # Skip the header and column names
    
    time_ps = data[:, 0]  # Column 1: Time (ps)
    fictive_temp = data[:, 1]  # Column 2: Fictive temperature (K)
    fictive_relax_time = data[:, 2]  # Column 3: Fictive relaxation time (ps)
    
    # Apply transformations as in plot_fictive_vs_waiting_relaxation.py:
    # 1. Shift time origin by -400 ps (cavity coupling starts at t=200 ps, analysis starts at t=400 ps)
    time_shifted = time_ps - 400.0
    
    # 2. Only keep data after t = 0 (after shifting)
    valid_mask = time_shifted >= 0.0
    time_final = time_shifted[valid_mask]
    fictive_relax_final = fictive_relax_time[valid_mask]
    
    print(f"  Loaded {len(time_final)} points after filtering (t >= 0)")
    print(f"  Time range: {time_final[0]:.1f} to {time_final[-1]:.1f} ps")
    print(f"  τ_fictive range: {fictive_relax_final.min():.1f} to {fictive_relax_final.max():.1f} ps")
    
    return {
        'time': time_final,
        'fictive_relaxation_time': fictive_relax_final,
        'coupling_value': coupling_value
    }

def load_fkt_material_time(material_time_file, coupling_value=0.0):
    """
    Load F(k,t) material time data for a specific coupling strength.
    
    Parameters:
    -----------
    material_time_file : str
        Path to the material time data file
    coupling_value : float
        Coupling strength to load (default: 0.0 for ε=0)
        
    Returns:
    --------
    fkt_material_time_data : dict or None
        Dictionary with time and material time data, or None if not found
    """
    if not os.path.exists(material_time_file):
        print(f"Error: F(k,t) material time file {material_time_file} not found")
        return None
    
    try:
        # Load the data, skipping comment lines
        data = np.loadtxt(material_time_file, comments='#')
        
        # Column structure from header:
        # t_g0.0_ps xi_g0.0 t_g3e-4_ps xi_g3e-4 t_g5e-4_ps xi_g5e-4 t_g7e-4_ps xi_g7e-4 t_g1e-3_ps xi_g1e-3
        
        coupling_column_map = {
            0.0: (0, 1),      # columns t_g0.0_ps, xi_g0.0
            3e-4: (2, 3),     # columns t_g3e-4_ps, xi_g3e-4
            5e-4: (4, 5),     # columns t_g5e-4_ps, xi_g5e-4
            7e-4: (6, 7),     # columns t_g7e-4_ps, xi_g7e-4
            1e-3: (8, 9)      # columns t_g1e-3_ps, xi_g1e-3
        }
        
        if coupling_value not in coupling_column_map:
            print(f"Error: Coupling value {coupling_value} not found in F(k,t) material time data")
            return None
        
        time_col, xi_col = coupling_column_map[coupling_value]
        
        if xi_col >= data.shape[1]:
            print(f"Error: Column {xi_col} not found in F(k,t) material time data")
            return None
        
        time_data = data[:, time_col]
        xi_data = data[:, xi_col]
        
        # Filter out any NaN or invalid data
        valid_mask = np.isfinite(time_data) & np.isfinite(xi_data)
        time_clean = time_data[valid_mask]
        xi_clean = xi_data[valid_mask]
        
        print(f"Loading F(k,t) material time for ε = {coupling_value}")
        print(f"  Loaded {len(time_clean)} points")
        print(f"  Time range: {time_clean[0]:.1f} to {time_clean[-1]:.1f} ps")
        print(f"  ξ range: {xi_clean.min():.3f} to {xi_clean.max():.3f}")
        
        return {
            'time': time_clean,
            'material_time': xi_clean,
            'coupling_value': coupling_value
        }
        
    except Exception as e:
        print(f"Error loading F(k,t) material time data: {e}")
        return None

def calculate_fictive_material_time(time, fictive_relaxation_time):
    """
    Calculate material time from fictive relaxation time.
    
    The material time is defined as:
    ξ_fictive(t) = ∫₀ᵗ (1/τ_fictive(t')) dt'
    
    Parameters:
    -----------
    time : array
        Time points (ps)
    fictive_relaxation_time : array
        Fictive relaxation time values (ps)
        
    Returns:
    --------
    material_time_fictive : array
        Material time from fictive integration
    """
    # Calculate 1/τ_fictive, avoiding division by zero
    inv_tau_fictive = np.where(fictive_relaxation_time > 0, 
                              1.0 / fictive_relaxation_time, 0.0)
    
    # Integrate using trapezoidal rule from t = 0
    material_time_fictive = np.zeros_like(time)
    for i in range(1, len(time)):
        dt = time[i] - time[i-1]
        # Trapezoidal integration
        material_time_fictive[i] = material_time_fictive[i-1] + 0.5 * (inv_tau_fictive[i-1] + inv_tau_fictive[i]) * dt
    
    return material_time_fictive

def find_time_mapping(time_fictive, xi_fictive, time_fkt, xi_fkt, method='slope'):
    """
    Find the mapping between fictive and F(k,t) material time scales.
    
    Parameters:
    -----------
    time_fictive : array
        Time points for fictive data
    xi_fictive : array
        Fictive material time
    time_fkt : array
        Time points for F(k,t) data
    xi_fkt : array
        F(k,t) material time
    method : str
        Method for finding mapping ('slope' or 'value')
        
    Returns:
    --------
    scale_factor : float
        Scale factor to convert ξ_fictive to ξ_fkt scale
    time_offset : float
        Time offset between the two datasets (ps)
    """
    print(f"\nFinding time mapping between fictive and F(k,t) material times...")
    print(f"Method: {method}")
    
    if method == 'slope':
        # Method 1: Compare slopes (assuming both are linear for ε=0)
        
        # Calculate slope of fictive material time
        # Use middle portion to avoid edge effects
        mid_start = len(time_fictive) // 4
        mid_end = 3 * len(time_fictive) // 4
        if mid_end > mid_start + 10:
            time_fit_fictive = time_fictive[mid_start:mid_end]
            xi_fit_fictive = xi_fictive[mid_start:mid_end]
            slope_fictive = np.polyfit(time_fit_fictive, xi_fit_fictive, 1)[0]
        else:
            slope_fictive = (xi_fictive[-1] - xi_fictive[0]) / (time_fictive[-1] - time_fictive[0])
        
        # Calculate slope of F(k,t) material time
        mid_start = len(time_fkt) // 4
        mid_end = 3 * len(time_fkt) // 4
        if mid_end > mid_start + 10:
            time_fit_fkt = time_fkt[mid_start:mid_end]
            xi_fit_fkt = xi_fkt[mid_start:mid_end]
            slope_fkt = np.polyfit(time_fit_fkt, xi_fit_fkt, 1)[0]
        else:
            slope_fkt = (xi_fkt[-1] - xi_fkt[0]) / (time_fkt[-1] - time_fkt[0])
        
        # Scale factor to match slopes
        scale_factor = slope_fkt / slope_fictive if slope_fictive != 0 else 1.0
        
        print(f"  Fictive slope: dξ/dt = {slope_fictive:.6f} ps⁻¹")
        print(f"  F(k,t) slope: dξ/dt = {slope_fkt:.6f} ps⁻¹")
        print(f"  Scale factor (slope ratio): {scale_factor:.6f}")
        
        # Time offset (assume they start at the same time)
        time_offset = 0.0
        
    elif method == 'value':
        # Method 2: Match values at a specific time point
        
        # Find a common time point (e.g., 1000 ps after t=0)
        target_time = 1000.0
        
        # Interpolate to find ξ values at target time
        if target_time <= time_fictive.max() and target_time <= time_fkt.max():
            xi_fictive_at_target = np.interp(target_time, time_fictive, xi_fictive)
            xi_fkt_at_target = np.interp(target_time, time_fkt, xi_fkt)
            
            scale_factor = xi_fkt_at_target / xi_fictive_at_target if xi_fictive_at_target != 0 else 1.0
            
            print(f"  At t = {target_time} ps:")
            print(f"    ξ_fictive = {xi_fictive_at_target:.6f}")
            print(f"    ξ_fkt = {xi_fkt_at_target:.6f}")
            print(f"    Scale factor (value ratio): {scale_factor:.6f}")
        else:
            print(f"  Warning: Target time {target_time} ps not available in both datasets")
            scale_factor = 1.0
        
        time_offset = 0.0
    
    else:
        print(f"  Error: Unknown mapping method '{method}'")
        scale_factor = 1.0
        time_offset = 0.0
    
    return scale_factor, time_offset

def plot_material_time_comparison(fictive_data, fkt_data, scale_factor, output_dir='.'):
    """
    Create comparison plot of material time scales.
    
    Parameters:
    -----------
    fictive_data : dict
        Fictive relaxation time data
    fkt_data : dict
        F(k,t) material time data
    scale_factor : float
        Scale factor to convert fictive to F(k,t) scale
    output_dir : str
        Output directory for plots
    """
    print("\nCreating material time comparison plot...")
    
    # Set up LaTeX rendering
    plt.style.use('classic')
    use_latex = True
    try:
        subprocess.check_output(['latex', '--version'], stderr=subprocess.DEVNULL)
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Computer Modern Roman']
        plt.rcParams['mathtext.fontset'] = 'cm'
        plt.rcParams['text.latex.preamble'] = r'\usepackage{lmodern}'
        print("  Using LaTeX with Computer Modern fonts")
    except (subprocess.CalledProcessError, FileNotFoundError, ImportError):
        print("  Warning: LaTeX not available, using default matplotlib rendering")
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'DejaVu Sans'
        use_latex = False
    
    # Helper function for LaTeX formatting
    def latex_format(text, fallback_text=None):
        if use_latex:
            return text
        else:
            return fallback_text if fallback_text else text.replace('$', '').replace(r'\xi', 'ξ').replace(r'\tau', 'τ')
    
    # Create figure with three panels
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    
    # Extract data
    time_fictive = fictive_data['time']
    tau_fictive = fictive_data['fictive_relaxation_time']
    
    time_fkt = fkt_data['time']
    xi_fkt = fkt_data['material_time']
    
    # Calculate fictive material time
    xi_fictive = calculate_fictive_material_time(time_fictive, tau_fictive)
    
    # Panel 1: Fictive relaxation time vs time
    ax1.plot(time_fictive, tau_fictive, 'b-', linewidth=2, alpha=0.8, label=r'$\tau_{\mathrm{fictive}}(t)$')
    ax1.set_xlabel(latex_format(r'Time (ps)'))
    ax1.set_ylabel(latex_format(r'$\tau_{\mathrm{fictive}}$ (ps)', 'τ_fictive (ps)'))
    ax1.set_title(latex_format(r'Fictive Relaxation Time ($\varepsilon = 0$)', 'Fictive Relaxation Time (ε = 0)'))
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Panel 2: Both material times (raw and scaled)
    ax2.plot(time_fictive, xi_fictive, 'b-', linewidth=2, alpha=0.8, 
            label=latex_format(r'$\xi_{\mathrm{fictive}}(t)$ (raw)', 'ξ_fictive(t) (raw)'))
    ax2.plot(time_fictive, xi_fictive * scale_factor, 'r--', linewidth=2, alpha=0.8,
            label=latex_format(f'$\\xi_{{\\mathrm{{fictive}}}}(t) \\times {scale_factor:.3f}$', 
                             f'ξ_fictive(t) × {scale_factor:.3f}'))
    ax2.plot(time_fkt, xi_fkt, 'g-', linewidth=2, alpha=0.8,
            label=latex_format(r'$\xi_{\mathrm{F(k,t)}}(t)$', 'ξ_F(k,t)(t)'))
    
    ax2.set_xlabel(latex_format(r'Time (ps)'))
    ax2.set_ylabel(latex_format(r'Material Time $\xi$', 'Material Time ξ'))
    ax2.set_title(latex_format(r'Material Time Comparison ($\varepsilon = 0$)', 'Material Time Comparison (ε = 0)'))
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Panel 3: Difference between scaled fictive and F(k,t) material times
    # Interpolate to common time grid for comparison
    time_common = np.linspace(max(time_fictive[0], time_fkt[0]), 
                             min(time_fictive[-1], time_fkt[-1]), 500)
    
    xi_fictive_interp = np.interp(time_common, time_fictive, xi_fictive * scale_factor)
    xi_fkt_interp = np.interp(time_common, time_fkt, xi_fkt)
    
    difference = xi_fictive_interp - xi_fkt_interp
    relative_difference = difference / xi_fkt_interp * 100  # Percentage
    
    ax3.plot(time_common, difference, 'k-', linewidth=2, alpha=0.8, label='Absolute difference')
    ax3_twin = ax3.twinx()
    ax3_twin.plot(time_common, relative_difference, 'orange', linewidth=1, alpha=0.7, label='Relative difference (%)')
    
    ax3.set_xlabel(latex_format(r'Time (ps)'))
    ax3.set_ylabel(latex_format(r'$\xi_{\mathrm{fictive}} - \xi_{\mathrm{F(k,t)}}$', 'ξ_fictive - ξ_F(k,t)'), color='k')
    ax3_twin.set_ylabel('Relative Difference (%)', color='orange')
    ax3.set_title(latex_format(r'Difference (Scaled Fictive - F(k,t))', 'Difference (Scaled Fictive - F(k,t))'))
    ax3.grid(True, alpha=0.3)
    
    # Add zero line
    ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    # Calculate and display statistics
    mean_abs_diff = np.mean(np.abs(difference))
    max_abs_diff = np.max(np.abs(difference))
    mean_rel_diff = np.mean(np.abs(relative_difference))
    
    stats_text = f'Mean |diff| = {mean_abs_diff:.4f}\nMax |diff| = {max_abs_diff:.4f}\nMean rel diff = {mean_rel_diff:.2f}%'
    ax3.text(0.02, 0.98, stats_text, transform=ax3.transAxes, fontsize=9,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    output_file = Path(output_dir) / 'material_time_scale_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    
    # Also save as PDF
    output_file_pdf = Path(output_dir) / 'material_time_scale_comparison.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Saved: {output_file_pdf}")
    
    plt.show()
    
    # Reset matplotlib settings
    plt.rcParams.update(plt.rcParamsDefault)
    
    return {
        'time_common': time_common,
        'xi_fictive_scaled': xi_fictive_interp,
        'xi_fkt': xi_fkt_interp,
        'difference': difference,
        'relative_difference': relative_difference,
        'mean_abs_diff': mean_abs_diff,
        'max_abs_diff': max_abs_diff,
        'mean_rel_diff': mean_rel_diff
    }

def main():
    """Main function."""
    
    # Directories
    base_dir = '/media/extradrive/Trajectories/final_nodiss_cavitymd'
    time_series_dir = Path(base_dir) / 'time_series_output'
    material_time_file = os.path.join(base_dir, 'material_time_all_couplings.txt')
    output_dir = base_dir
    
    print("=" * 70)
    print("MATERIAL TIME SCALE COMPARISON: FICTIVE vs F(k,t) ANALYSIS")
    print("=" * 70)
    print("Goal: Find mapping between ξ_fictive = ∫(1/τ_fictive)dt and ξ_F(k,t)")
    print("Focus: ε = 0 case where both should be linear in time")
    
    # Target coupling for initial comparison
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
    
    # 3. Calculate fictive material time
    print(f"\n3. Calculating material time from fictive data...")
    time_fictive = fictive_data['time']
    tau_fictive = fictive_data['fictive_relaxation_time']
    xi_fictive = calculate_fictive_material_time(time_fictive, tau_fictive)
    
    print(f"  Fictive material time range: {xi_fictive.min():.6f} to {xi_fictive.max():.6f}")
    print(f"  Average fictive relaxation time: {tau_fictive.mean():.1f} ps")
    
    # Add calculated material time to fictive data
    fictive_data['material_time_fictive'] = xi_fictive
    
    # 4. Find mapping between scales
    print(f"\n4. Finding mapping between material time scales...")
    
    # Try both slope and value methods
    scale_factor_slope, time_offset_slope = find_time_mapping(
        time_fictive, xi_fictive, fkt_data['time'], fkt_data['material_time'], method='slope')
    
    scale_factor_value, time_offset_value = find_time_mapping(
        time_fictive, xi_fictive, fkt_data['time'], fkt_data['material_time'], method='value')
    
    # Use slope method as primary (more robust for linear trends)
    scale_factor = scale_factor_slope
    print(f"\nUsing slope-based scale factor: {scale_factor:.6f}")
    
    # 5. Create comparison plot
    print(f"\n5. Creating comparison plots...")
    comparison_results = plot_material_time_comparison(fictive_data, fkt_data, scale_factor, output_dir)
    
    # 6. Summary and interpretation
    print(f"\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"Coupling strength: ε = {coupling_value}")
    print(f"Scale factor (slope method): {scale_factor:.6f}")
    print(f"Scale factor (value method): {scale_factor_value:.6f}")
    print(f"\nMapping relationship:")
    print(f"  ξ_F(k,t) ≈ {scale_factor:.6f} × ξ_fictive")
    print(f"  or equivalently: ξ_fictive ≈ {1/scale_factor:.6f} × ξ_F(k,t)")
    
    print(f"\nComparison statistics (after scaling):")
    print(f"  Mean absolute difference: {comparison_results['mean_abs_diff']:.6f}")
    print(f"  Maximum absolute difference: {comparison_results['max_abs_diff']:.6f}")
    print(f"  Mean relative difference: {comparison_results['mean_rel_diff']:.2f}%")
    
    # Interpretation
    print(f"\nInterpretation:")
    if comparison_results['mean_rel_diff'] < 5:
        print("   Excellent agreement between fictive and F(k,t) material times!")
    elif comparison_results['mean_rel_diff'] < 15:
        print("   Good agreement between fictive and F(k,t) material times")
    else:
        print("    Significant differences between the two material time scales")
    
    print(f"\nPhysical meaning:")
    print(f"  Both methods calculate ξ(t) = ∫(1/τ(t'))dt', but with different τ sources:")
    print(f"  - Fictive: τ_fictive(t) from fictive temperature relaxation")
    print(f"  - F(k,t): Effective τ from material time construction (Δξ=1 ↔ R=0.1)")
    print(f"  The scale factor {scale_factor:.3f} relates these two time scales.")
    
    print(f"\nGenerated files:")
    print(f"  - material_time_scale_comparison.png")
    print(f"  - material_time_scale_comparison.pdf")
    
    return 0

if __name__ == "__main__":
    exit(main())
