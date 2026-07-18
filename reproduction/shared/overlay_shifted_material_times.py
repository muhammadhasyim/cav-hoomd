#!/usr/bin/env python3
"""
Overlay material times with proper time shifting.

Key idea: Shift fictive material time so that ξ_fictive(t=200ps) becomes ξ_fictive(t=0ps).
This creates a proper comparison where both material times start from their respective origins.

Author: Assistant
Date: September 12, 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from compare_material_time_t200 import load_fictive_data, load_fkt_material_time, calculate_fictive_material_time_t200

def calculate_time_shifted_fictive_material_time(time, fictive_relaxation_time, t_shift=200.0):
    """
    Calculate fictive material time and shift both time and material time origins.
    
    Time shift: t_new = t - t_shift (so t=200ps becomes t=0ps)
    Material time shift: ξ_new(t) = ξ(t) - ξ(t_shift) (so ξ(t=200ps) becomes ξ=0)
    
    Parameters:
    -----------
    time : array
        Original time points (ps)
    fictive_relaxation_time : array
        Fictive relaxation time values (ps)
    t_shift : float
        Time point to shift to new origin (default: 200.0 ps)
        
    Returns:
    --------
    time_shifted : array
        Time shifted to new origin: t_new = t - t_shift
    xi_shifted : array  
        Material time shifted to new origin: ξ_new = ξ - ξ(t_shift)
    """
    # Calculate material time starting from t_shift
    xi_raw = calculate_fictive_material_time_t200(time, fictive_relaxation_time, t_start=t_shift)
    
    # Find the shift point
    shift_idx = np.argmin(np.abs(time - t_shift))
    actual_t_shift = time[shift_idx]
    xi_at_shift = xi_raw[shift_idx]
    
    print(f"  Time shifting: t = {actual_t_shift:.3f} ps → t = 0 ps")
    print(f"  Material time shifting: ξ = {xi_at_shift:.6f} → ξ = 0")
    
    # Apply shifts
    # Only keep data for t >= t_shift (where we have valid material time data)
    valid_mask = time >= actual_t_shift
    
    time_shifted = time[valid_mask] - actual_t_shift  # t_new = t - t_shift
    xi_shifted = xi_raw[valid_mask] - xi_at_shift     # ξ_new = ξ - ξ(t_shift)
    
    print(f"  Shifted time range: {time_shifted.min():.3f} to {time_shifted.max():.3f} ps")
    print(f"  Shifted material time range: {xi_shifted.min():.6f} to {xi_shifted.max():.6f}")
    
    return time_shifted, xi_shifted

def create_shifted_overlay_plot():
    """Create overlay plot with proper time shifting."""
    
    # Directories
    base_dir = '/media/extradrive/Trajectories/final_nodiss_cavitymd'
    time_series_dir = Path(base_dir) / 'time_series_output'
    material_time_file = os.path.join(base_dir, 'material_time_all_couplings.txt')
    output_dir = base_dir
    
    # All available coupling values
    coupling_values = [0.0, 3e-4, 5e-4, 7e-4, 1e-3]
    
    # Theoretical scale factor
    ln10_factor = 1.0 / np.log(10)
    
    print("=" * 80)
    print("TIME-SHIFTED MATERIAL TIME OVERLAY")
    print("=" * 80)
    print(f"Fictive time shift: t = 200 ps → t = 0 ps")
    print(f"Fictive material time shift: ξ(t=200ps) → ξ = 0")
    print(f"F(k,t) reference: unchanged (already starts at t = 0 ps)")
    print(f"Scale factor: 1/ln(10) = {ln10_factor:.6f}")
    
    # Use simple matplotlib settings to avoid LaTeX issues
    plt.style.use('default')
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.alpha'] = 0.3
    
    # Create figure with two panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Color scheme
    from matplotlib.colors import Normalize
    coupling_norm = Normalize(vmin=0, vmax=max(coupling_values))
    coupling_cmap = plt.colormaps.get_cmap('coolwarm')
    
    def get_coupling_color(coupling_value):
        return coupling_cmap(coupling_norm(coupling_value))
    
    # Process each coupling strength
    all_data = {}
    
    for coupling_value in coupling_values:
        print(f"\nProcessing ε = {coupling_value}...")
        
        # Load data
        fictive_data = load_fictive_data(time_series_dir, coupling_value)
        fkt_data = load_fkt_material_time(material_time_file, coupling_value)
        
        if fictive_data is None or fkt_data is None:
            print(f"  Could not load data for ε = {coupling_value}")
            continue
        
        # Calculate time-shifted fictive material time
        time_fictive = fictive_data['time']
        tau_fictive = fictive_data['fictive_relaxation_time']
        
        time_fictive_shifted, xi_fictive_shifted = calculate_time_shifted_fictive_material_time(
            time_fictive, tau_fictive, t_shift=200.0)
        
        # Apply 1/ln(10) scaling
        xi_fictive_scaled = xi_fictive_shifted * ln10_factor
        
        # Store data
        all_data[coupling_value] = {
            'time_fictive_shifted': time_fictive_shifted,
            'xi_fictive_shifted': xi_fictive_shifted,
            'xi_fictive_scaled': xi_fictive_scaled,
            'time_fkt': fkt_data['time'],
            'xi_fkt': fkt_data['material_time']
        }
        
        # Calculate comparison statistics for t > 200 ps (in shifted time: t > 0)
        if len(time_fictive_shifted) > 50 and len(fkt_data['time']) > 50:
            t_analysis_start = 200.0  # In shifted fictive time, this is t = 0
            t_analysis_end = min(time_fictive_shifted[-1] + 200, fkt_data['time'][-1])
            
            # Common time grid (in lab time coordinates)
            time_common = np.linspace(t_analysis_start, t_analysis_end, 200)
            
            # Interpolate fictive data (convert back to lab time for interpolation)
            time_fictive_lab = time_fictive_shifted + 200.0  # Convert back to lab time
            xi_fictive_interp = np.interp(time_common, time_fictive_lab, xi_fictive_scaled)
            
            # Interpolate F(k,t) data
            xi_fkt_interp = np.interp(time_common, fkt_data['time'], fkt_data['material_time'])
            
            # Calculate relative difference
            valid_mask = xi_fkt_interp > 0.01  # Avoid division by very small numbers
            if np.sum(valid_mask) > 10:
                relative_diff = np.abs((xi_fictive_interp[valid_mask] - xi_fkt_interp[valid_mask]) / xi_fkt_interp[valid_mask] * 100)
                mean_rel_diff = np.mean(relative_diff)
                print(f"  Mean relative difference (t>200ps): {mean_rel_diff:.2f}%")
    
    # Panel 1: Detailed comparison for ε = 0
    print(f"\nCreating detailed comparison for ε = 0...")
    
    if 0.0 in all_data:
        data_eps0 = all_data[0.0]
        
        # Plot F(k,t) material time (reference)
        ax1.plot(data_eps0['time_fkt'], data_eps0['xi_fkt'], 
                'g-', linewidth=3, alpha=0.8, 
                label='ξ_F(k,t) (reference, t₀=0ps)')
        
        # Plot time-shifted and scaled fictive material time
        # Convert fictive time back to lab time for plotting
        time_fictive_lab = data_eps0['time_fictive_shifted'] + 200.0
        ax1.plot(time_fictive_lab, data_eps0['xi_fictive_scaled'], 
                'r--', linewidth=2, alpha=0.8, 
                label='ξ_fictive/ln(10) (shifted, t₀=200ps→0ps)')
        
        # Add vertical line at t = 200 ps
        ax1.axvline(x=200, color='gray', linestyle=':', alpha=0.7, linewidth=2,
                   label='Original t = 200 ps')
        
        ax1.set_xlabel('Lab Time (ps)')
        ax1.set_ylabel('Material Time ξ')
        ax1.set_title('ε = 0: Time-Shifted Comparison')
        ax1.legend(loc='upper left')
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(0, 3000)
        
        # Add statistics text box
        time_common = np.linspace(400, 2800, 200)
        xi_fictive_interp = np.interp(time_common, time_fictive_lab, data_eps0['xi_fictive_scaled'])
        xi_fkt_interp = np.interp(time_common, data_eps0['time_fkt'], data_eps0['xi_fkt'])
        
        valid_mask = xi_fkt_interp > 0.01
        if np.sum(valid_mask) > 10:
            mean_abs_diff = np.mean(np.abs(xi_fictive_interp[valid_mask] - xi_fkt_interp[valid_mask]))
            mean_rel_diff = np.mean(np.abs((xi_fictive_interp[valid_mask] - xi_fkt_interp[valid_mask]) / xi_fkt_interp[valid_mask] * 100))
            
            stats_text = f'1/ln(10) = {ln10_factor:.4f}\nMean |diff| = {mean_abs_diff:.4f}\nMean rel diff = {mean_rel_diff:.2f}%\n(t > 400 ps)'
            ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=10,
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Panel 2: Overlay with F(k,t) in original coordinates, fictive shifted
    print(f"\nCreating overlay: F(k,t) original + fictive shifted...")
    
    for coupling_value in coupling_values:
        if coupling_value not in all_data:
            continue
        
        data = all_data[coupling_value]
        color = get_coupling_color(coupling_value)
        
        # Format coupling label
        if coupling_value == 0.0:
            label_base = 'ε = 0'
        elif coupling_value >= 1e-3:
            label_base = f'ε = {coupling_value:.0e}'
        else:
            scaled_val = coupling_value * 10000
            label_base = f'ε = {scaled_val:.0f}×10⁻⁴'
        
        # Plot F(k,t) material time (original coordinates: no shift)
        ax2.plot(data['time_fkt'], data['xi_fkt'], 
                '-', color=color, linewidth=2, alpha=0.8, 
                label=f'{label_base} (F(k,t))')
        
        # Plot scaled fictive material time (keep shifted coordinates)
        # The fictive time is already shifted: time_fictive_shifted = t - 200
        # The fictive ξ is already shifted: xi_fictive_scaled with ξ(200ps) = 0
        # Plot in original lab time but with the shifted ξ values
        time_fictive_lab = data['time_fictive_shifted'] + 200.0  # Convert time back to lab coordinates
        ax2.plot(time_fictive_lab-200, data['xi_fictive_scaled'], 
                '--', color=color, linewidth=2, alpha=0.8, 
                label=f'{label_base} (fictive)')
    
    # Add vertical line at t = 200 ps to show the fictive starting point
    ax2.axvline(x=200, color='gray', linestyle=':', alpha=0.7, linewidth=2,
               label='t = 200 ps (fictive origin)')
    
    ax2.set_xlabel('Lab Time (ps)')
    ax2.set_ylabel('Material Time ξ')
    ax2.set_title('All Couplings: F(k,t) vs Shifted Fictive')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 3000)
    
    # Simplified legend for panel 2
    handles, labels = ax2.get_legend_handles_labels()
    # Show only F(k,t) entries to avoid clutter
    fkt_handles = []
    fkt_labels = []
    for handle, label in zip(handles, labels):
        if 'F(k,t)' in label:
            fkt_labels.append(label.replace(' (F(k,t))', ''))
            fkt_handles.append(handle)
    
    ax2.legend(fkt_handles, fkt_labels, loc='upper left', fontsize=10, title='Coupling Strengths')
    
    # Add text box explaining the transformation
    transform_text = 'Fictive: t\' = t - 200ps, ξ\'(200ps) = 0\nF(k,t): original coordinates\nScale: ξ_F(k,t) = ξ_fictive / ln(10)'
    ax2.text(0.02, 0.75, transform_text, transform=ax2.transAxes, fontsize=9,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    output_file = Path(output_dir) / 'material_time_shifted_overlay.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved: {output_file}")
    
    # Also save as PDF
    output_file_pdf = Path(output_dir) / 'material_time_shifted_overlay.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Saved: {output_file_pdf}")
    
    plt.show()
    
    # Print summary
    print(f"\n" + "=" * 80)
    print("TIME-SHIFTED OVERLAY ANALYSIS SUMMARY")
    print("=" * 80)
    print(f"Time transformation:")
    print(f"  • Fictive: t' = t - 200 ps, ξ'(t) = ξ(t) - ξ(200ps)")
    print(f"    (so fictive origin becomes t=200ps → ξ=0)")
    print(f"  • F(k,t): UNCHANGED (keeps original t=0ps → ξ=0 origin)")
    print(f"")
    print(f"Universal scaling relationship:")
    print(f"  ξ_F(k,t)(t) = ξ'_fictive(t-200ps) / ln(10)")
    print(f"  Scale factor: 1/ln(10) = {ln10_factor:.6f}")
    print(f"")
    print(f"Key advantages:")
    print(f"   Fictive starts aging from its natural t=200ps origin")
    print(f"   F(k,t) maintains its natural t=0ps origin")
    print(f"   Direct comparison shows both track same dynamics")
    print(f"   Physically meaningful: fictive aging offset matches coupling turn-on")
    
    return all_data


def load_actual_fkt_data(coupling_value=0.0, ref_num=1):
    """
    Load actual F(k,t) data for a specific coupling and reference time.
    
    Parameters:
    -----------
    coupling_value : float
        Coupling strength (default: 0.0 for ε=0)
    ref_num : int
        Reference time number (default: 1 for tw=200ps)
        
    Returns:
    --------
    time_lag : array
        Lag time values (ps)
    fkt_actual : array
        Actual F(k,t) values
    """
    # Find the appropriate coupling directory
    coupling_dirs = [d for d in Path('.').iterdir() if d.is_dir() and 'cavity_coupling' in d.name and 'switch_200.0ps' in d.name]
    
    target_dir = None
    for d in coupling_dirs:
        if coupling_value == 0.0 and '0epos00' in d.name:
            target_dir = d
            break
        elif coupling_value > 0:
            # Parse coupling from directory name
            if 'eneg' in d.name:
                # Extract the numeric part before 'eneg'
                coupling_part = d.name.replace('cavity_coupling_', '').replace('_switch_200.0ps', '')
                parts = coupling_part.split('eneg')
                if len(parts) == 2:
                    try:
                        mantissa = int(parts[0])
                        exponent_part = parts[1]
                        exponent = int(exponent_part)  # Full exponent number (e.g., "04" -> 4)
                        dir_coupling = mantissa * (10 ** (-exponent))
                        print(f"    Checking directory {d.name}: parsed coupling = {dir_coupling:.6f}")
                        if abs(dir_coupling - coupling_value) < 1e-6:
                            target_dir = d
                            break
                    except Exception as e:
                        print(f"    Failed to parse {d.name}: {e}")
                        continue
    
    if target_dir is None:
        raise FileNotFoundError(f"Could not find directory for coupling {coupling_value}")
    
    # Load the specific reference file
    ref_file = target_dir / f'master_fskt_ref{ref_num}.txt'
    if not ref_file.exists():
        raise FileNotFoundError(f"Could not find {ref_file}")
    
    print(f"Loading actual F(k,t) data from: {ref_file}")
    
    # Read the data, skipping header lines
    data = np.loadtxt(ref_file, skiprows=9)  # Skip the header lines
    time_lag = data[:, 0]  # First column is lag time
    fkt_actual = data[:, 1]  # Second column is F(k,t)
    
    return time_lag, fkt_actual


def generate_fkt_from_fictive(time_lab, xi_fictive_scaled, tw=200.0, beta=0.55):
    """
    Generate F(k,t; tw) from fictive material time using stretched exponential.
    
    The key insight is that F(k,t; tw) depends on the material time difference:
    F(k,t; tw) = exp(-ln(10) * [ξ(t+tw) - ξ(tw)]^β)
    
    Parameters:
    -----------
    time_lab : array
        Laboratory time points (ps)
    xi_fictive_scaled : array
        Scaled fictive material time ξ = ∫ (1/τ_fictive) dt / ln(10)
    tw : float
        Waiting time (ps), default 200.0 for ref1
    beta : float
        Stretched exponential exponent, default 0.55
        
    Returns:
    --------
    time_lag : array
        Lag time values (t for F(k,t; tw))
    fkt_fictive : array
        F(k,t; tw) calculated from fictive material time
    """
    print(f"Generating F(k,t; tw={tw}ps) from fictive material time with β={beta}")
    
    # Interpolate to get exact ξ(tw) at the requested waiting time
    from scipy.interpolate import interp1d
    
    # Create interpolation function for material time
    xi_interp = interp1d(time_lab, xi_fictive_scaled, kind='linear', 
                        bounds_error=False, fill_value='extrapolate')
    
    # Get exact ξ(tw) at the requested waiting time
    xi_tw = xi_interp(tw)
    
    print(f"  Waiting time: tw = {tw:.1f} ps, ξ(tw) = {xi_tw:.6f} (interpolated)")
    
    # For times t >= tw, calculate lag time and material time difference
    valid_mask = time_lab >= tw
    time_valid = time_lab[valid_mask]
    xi_valid = xi_fictive_scaled[valid_mask]
    
    # Calculate lag time: t_lag = t - tw
    time_lag = time_valid - tw
    
    # Calculate material time difference: Δξ = ξ(t+tw) - ξ(tw)
    # Note: xi_valid already represents ξ(t+tw) for times t >= tw
    delta_xi = xi_valid - xi_tw
    
    # Calculate F(k,t; tw) using stretched exponential
    # F(k,t; tw) = exp(-ln(10) * Δξ^β)
    fkt_fictive = np.exp(-np.log(10) * delta_xi**beta)
    
    print(f"  Generated {len(fkt_fictive)} F(k,t) points")
    print(f"  Lag time range: {time_lag[0]:.1f} - {time_lag[-1]:.1f} ps")
    print(f"  Material time difference range: {delta_xi[0]:.6f} - {delta_xi[-1]:.6f}")
    print(f"  F(k,t) range: {fkt_fictive[-1]:.6f} - {fkt_fictive[0]:.6f}")
    
    return time_lag, fkt_fictive


def compare_fkt_fictive_vs_actual(all_data, coupling_value=0.0, ref_num=1, beta=0.55, output_dir='.'):
    """
    Compare F(k,t; tw) generated from fictive material time with actual F(k,t) data.
    
    Parameters:
    -----------
    all_data : dict
        Data from create_shifted_overlay_plot()
    coupling_value : float
        Coupling strength to compare (default: 0.0 for ε=0)
    ref_num : int
        Reference time number (default: 1 for tw=200ps)
    beta : float
        Stretched exponential exponent (default: 0.55)
    output_dir : str
        Output directory for plots
    """
    print(f"\n" + "=" * 80)
    print(f"COMPARING F(k,t; tw) FROM FICTIVE VS ACTUAL DATA")
    print("=" * 80)
    print(f"Coupling strength: {coupling_value}")
    print(f"Reference time: ref{ref_num} (tw = {ref_num * 200} ps)")
    print(f"Stretched exponential exponent: β = {beta}")
    
    # Find the matching coupling in all_data
    coupling_key = None
    for key in all_data.keys():
        if abs(key - coupling_value) < 1e-6:
            coupling_key = key
            break
    
    if coupling_key is None:
        raise ValueError(f"Could not find coupling {coupling_value} in loaded data")
    
    coupling_data = all_data[coupling_key]
    
    # Generate F(k,t; tw) from fictive material time
    tw = ref_num * 200.0  # Convert ref number to actual time
    time_lag_fictive, fkt_fictive = generate_fkt_from_fictive(
        coupling_data['time_fictive_shifted'] + 200.0,  # Convert back to lab time
        coupling_data['xi_fictive_scaled'],
        tw=tw,
        beta=beta
    )
    
    # Load actual F(k,t) data
    print(f"\nLoading actual F(k,t) data...")
    time_lag_actual, fkt_actual_raw = load_actual_fkt_data(coupling_value, ref_num)
    
    # Normalize actual F(k,t) data to start at 1.0 (same as fictive)
    fkt_actual = fkt_actual_raw / fkt_actual_raw[0]
    print(f"  Loaded {len(fkt_actual)} F(k,t) points")
    print(f"  Raw F(k,t) range: {fkt_actual_raw[0]:.1f} - {fkt_actual_raw[-1]:.1f}")
    print(f"  Normalized F(k,t) range: {fkt_actual[0]:.6f} - {fkt_actual[-1]:.6f}")
    
    # Create comparison plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Panel 1: Direct comparison
    ax1.plot(time_lag_actual, fkt_actual, 'b-', linewidth=2, alpha=0.8, 
             label=f'Actual F(k,t; tw={tw:.0f}ps) [normalized]')
    ax1.plot(time_lag_fictive, fkt_fictive, 'r--', linewidth=2, alpha=0.8,
             label=f'Fictive-based F(k,t; tw={tw:.0f}ps)')
    
    ax1.set_xlabel('Lag Time t (ps)')
    ax1.set_ylabel('F(k,t) [normalized]')
    ax1.set_title(f'F(k,t) Comparison: Fictive vs Actual (Normalized)\nCoupling g = {coupling_value}, tw = {tw:.0f} ps, β = {beta}')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1.1)
    
    # Panel 2: Log-scale comparison
    ax2.semilogy(time_lag_actual, fkt_actual, 'b-', linewidth=2, alpha=0.8,
                 label=f'Actual F(k,t; tw={tw:.0f}ps) [normalized]')
    ax2.semilogy(time_lag_fictive, fkt_fictive, 'r--', linewidth=2, alpha=0.8,
                 label=f'Fictive-based F(k,t; tw={tw:.0f}ps)')
    
    ax2.set_xlabel('Lag Time t (ps)')
    ax2.set_ylabel('F(k,t) [normalized, log scale]')
    ax2.set_title(f'F(k,t) Comparison (Log Scale, Normalized)\nCoupling g = {coupling_value}, tw = {tw:.0f} ps, β = {beta}')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(1e-3, 1.1)
    
    plt.tight_layout()
    
    # Save plot
    coupling_str = f"eps{coupling_value:.4f}".replace(".", "p")
    output_file = Path(output_dir) / f'fkt_comparison_fictive_vs_actual_{coupling_str}_ref{ref_num}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved F(k,t) comparison: {output_file}")
    
    # Also save as PDF
    output_file_pdf = Path(output_dir) / f'fkt_comparison_fictive_vs_actual_{coupling_str}_ref{ref_num}.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Saved: {output_file_pdf}")
    
    plt.show()
    
    # Calculate and report statistics
    print(f"\n" + "-" * 60)
    print("COMPARISON STATISTICS")
    print("-" * 60)
    
    # Interpolate to common time grid for comparison
    time_common = np.linspace(max(time_lag_actual.min(), time_lag_fictive.min()),
                             min(time_lag_actual.max(), time_lag_fictive.max()), 100)
    
    fkt_actual_interp = np.interp(time_common, time_lag_actual, fkt_actual)
    fkt_fictive_interp = np.interp(time_common, time_lag_fictive, fkt_fictive)
    
    # Calculate correlation coefficient
    correlation = np.corrcoef(fkt_actual_interp, fkt_fictive_interp)[0, 1]
    
    # Calculate mean squared error
    mse = np.mean((fkt_actual_interp - fkt_fictive_interp)**2)
    rmse = np.sqrt(mse)
    
    # Calculate relative error
    rel_error = np.mean(np.abs(fkt_actual_interp - fkt_fictive_interp) / fkt_actual_interp) * 100
    
    print(f"Correlation coefficient: {correlation:.6f}")
    print(f"Root mean squared error: {rmse:.6f}")
    print(f"Mean relative error: {rel_error:.2f}%")
    
    if correlation > 0.99:
        print(" EXCELLENT agreement between fictive and actual F(k,t)")
    elif correlation > 0.95:
        print(" GOOD agreement between fictive and actual F(k,t)")
    elif correlation > 0.90:
        print("  MODERATE agreement between fictive and actual F(k,t)")
    else:
        print(" POOR agreement between fictive and actual F(k,t)")


def create_all_coupling_overlay_plots(all_data, beta=0.55, output_dir='.'):
    """
    Create overlay plots showing all coupling strengths together for each reference time.
    
    Parameters:
    -----------
    all_data : dict
        Data from create_shifted_overlay_plot()
    beta : float
        Stretched exponential exponent (default: 0.55)
    output_dir : str
        Output directory for plots
    """
    print(f"\n" + "=" * 80)
    print("CREATING ALL-COUPLING OVERLAY PLOTS FOR ALL REFERENCE TIMES")
    print("=" * 80)
    
    coupling_values = sorted(all_data.keys())
    print(f"Processing {len(coupling_values)} coupling strengths: {coupling_values}")
    
    # Define color scheme for coupling strengths
    from matplotlib.colors import Normalize
    import matplotlib.cm as cm
    coupling_norm = Normalize(vmin=0, vmax=max(coupling_values))
    coupling_cmap = plt.colormaps.get_cmap('coolwarm')
    
    def get_coupling_color(coupling_value):
        return coupling_cmap(coupling_norm(coupling_value))
    
    # Find all available reference times by checking the first coupling
    first_coupling = coupling_values[0]
    
    # Check what reference files are available
    available_refs = []
    coupling_dirs = [d for d in Path('.').iterdir() if d.is_dir() and 'cavity_coupling' in d.name and 'switch_200.0ps' in d.name]
    
    # Use the first directory to check available reference files
    if coupling_dirs:
        first_dir = coupling_dirs[0]
        ref_files = list(first_dir.glob('master_fskt_ref*.txt'))
        for ref_file in ref_files:
            # Extract reference number from filename
            stem = ref_file.stem
            if stem.startswith('master_fskt_ref') and not '_sample_counts' in stem:
                try:
                    ref_part = stem.replace('master_fskt_ref', '')
                    ref_num = int(ref_part)
                    if ref_num >= 1:  # Start from ref1
                        available_refs.append(ref_num)
                except ValueError:
                    continue
    
    available_refs = sorted(available_refs)
    print(f"Found reference times: {available_refs} (corresponding to tw = {[r*200 for r in available_refs]} ps)")
    
    generated_files = []
    
    for ref_num in available_refs:
        tw = ref_num * 200.0
        print(f"\n  -> Creating overlay for ref{ref_num} (tw = {tw:.0f} ps)...")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        
        plot_data = []
        
        for coupling_val in coupling_values:
            try:
                print(f"    Processing ε = {coupling_val}...")
                
                # Generate F(k,t) from fictive material time
                coupling_data = all_data[coupling_val]
                time_lag_fictive, fkt_fictive = generate_fkt_from_fictive(
                    coupling_data['time_fictive_shifted'] + 200.0,
                    coupling_data['xi_fictive_scaled'],
                    tw=tw,
                    beta=beta
                )
                
                # Load actual F(k,t) data
                time_lag_actual, fkt_actual_raw = load_actual_fkt_data(coupling_val, ref_num)
                fkt_actual = fkt_actual_raw / fkt_actual_raw[0]  # Normalize
                
                color = get_coupling_color(coupling_val)
                
                # Plot on both panels
                coupling_label = f'ε = {coupling_val:.4f}' if coupling_val > 0 else 'ε = 0.0'
                
                # Panel 1: Linear scale
                ax1.plot(time_lag_actual, fkt_actual, '-', color=color, linewidth=2, alpha=0.7,
                         label=f'{coupling_label} (actual)')
                ax1.plot(time_lag_fictive, fkt_fictive, '--', color=color, linewidth=2, alpha=0.8,
                         label=f'{coupling_label} (fictive)')
                
                # Panel 2: Log scale
                ax2.semilogy(time_lag_actual, fkt_actual, '-', color=color, linewidth=2, alpha=0.7)
                ax2.semilogy(time_lag_fictive, fkt_fictive, '--', color=color, linewidth=2, alpha=0.8)
                
                # Calculate correlation for summary
                time_common = np.linspace(max(time_lag_actual.min(), time_lag_fictive.min()),
                                         min(time_lag_actual.max(), time_lag_fictive.max()), 100)
                fkt_actual_interp = np.interp(time_common, time_lag_actual, fkt_actual)
                fkt_fictive_interp = np.interp(time_common, time_lag_fictive, fkt_fictive)
                correlation = np.corrcoef(fkt_actual_interp, fkt_fictive_interp)[0, 1]
                
                plot_data.append({
                    'coupling': coupling_val,
                    'correlation': correlation
                })
                
            except Exception as e:
                print(f"      Warning: Could not process ε = {coupling_val}: {e}")
                continue
        
        # Customize Panel 1
        ax1.set_xlabel('Lag Time t (ps)')
        ax1.set_ylabel('F(k,t) [normalized]')
        ax1.set_title(f'F(k,t) All Couplings: Fictive vs Actual\nReference tw = {tw:.0f} ps, β = {beta}')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, 1.1)
        
        # Create custom legend for Panel 1 (show only coupling strengths)
        coupling_legend_elements = []
        for coupling_val in coupling_values:
            color = get_coupling_color(coupling_val)
            coupling_label = f'ε = {coupling_val:.4f}' if coupling_val > 0 else 'ε = 0.0'
            coupling_legend_elements.append(plt.Line2D([0], [0], color=color, linewidth=2, label=coupling_label))
        
        ax1.legend(handles=coupling_legend_elements, loc='upper right', fontsize=9, title='Coupling Strengths')
        
        # Customize Panel 2
        ax2.set_xlabel('Lag Time t (ps)')
        ax2.set_ylabel('F(k,t) [normalized, log scale]')
        ax2.set_title(f'F(k,t) All Couplings (Log Scale)\nReference tw = {tw:.0f} ps, β = {beta}')
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(1e-3, 1.1)
        
        # Add legend for line styles on Panel 2
        style_legend_elements = [
            plt.Line2D([0], [0], color='gray', linestyle='-', linewidth=2, label='Actual F(k,t)'),
            plt.Line2D([0], [0], color='gray', linestyle='--', linewidth=2, label='Fictive-based F(k,t)')
        ]
        ax2.legend(handles=style_legend_elements, loc='upper right', fontsize=9, title='Data Types')
        
        plt.tight_layout()
        
        # Save plot
        output_file = Path(output_dir) / f'fkt_all_couplings_overlay_ref{ref_num}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"    Saved: {output_file}")
        generated_files.append(output_file.name)
        
        # Also save as PDF
        output_file_pdf = Path(output_dir) / f'fkt_all_couplings_overlay_ref{ref_num}.pdf'
        plt.savefig(output_file_pdf, bbox_inches='tight')
        
        plt.show()
        
        # Print correlation summary for this reference
        print(f"    Correlation summary for ref{ref_num}:")
        for data in plot_data:
            print(f"      ε = {data['coupling']:.4f}: r = {data['correlation']:.4f}")
    
    return generated_files


def find_relaxation_time_01_criterion(time, fkt):
    """
    Find relaxation time where F(k,t) = 0.1 using linear interpolation.
    
    Parameters:
    -----------
    time : array
        Time values (ps)
    fkt : array
        F(k,t) values (normalized to start at 1.0)
        
    Returns:
    --------
    tau : float
        Relaxation time where F(k,t) = 0.1, or np.nan if not found
    """
    target = 0.1
    
    # Find where F(k,t) crosses 0.1
    if fkt[0] < target or fkt[-1] > target:
        return np.nan
    
    # Find the crossing point
    for i in range(len(fkt) - 1):
        if fkt[i] >= target and fkt[i+1] <= target:
            # Linear interpolation between the two points
            t1, t2 = time[i], time[i+1]
            f1, f2 = fkt[i], fkt[i+1]
            tau = t1 + (target - f1) * (t2 - t1) / (f2 - f1)
            return tau
    
    return np.nan


def compare_relaxation_times_analysis(all_data, beta=0.55, output_dir='.'):
    """
    Compare relaxation times (F(k,t) = 0.1 criterion) from actual vs fictive-based F(k,t).
    Create plots similar to normalized_relaxation_analysis_filtered.png.
    
    Parameters:
    -----------
    all_data : dict
        Data from create_shifted_overlay_plot()
    beta : float
        Stretched exponential exponent (default: 0.55)
    output_dir : str
        Output directory for plots
    """
    print(f"\n" + "=" * 80)
    print("RELAXATION TIME COMPARISON: ACTUAL vs FICTIVE-BASED F(k,t)")
    print("=" * 80)
    print(f"Using F(k,t) = 0.1 criterion, β = {beta}")
    
    coupling_values = sorted(all_data.keys())
    print(f"Processing {len(coupling_values)} coupling strengths: {coupling_values}")
    
    # Find all available reference times
    coupling_dirs = [d for d in Path('.').iterdir() if d.is_dir() and 'cavity_coupling' in d.name and 'switch_200.0ps' in d.name]
    available_refs = []
    
    if coupling_dirs:
        first_dir = coupling_dirs[0]
        ref_files = list(first_dir.glob('master_fskt_ref*.txt'))
        for ref_file in ref_files:
            stem = ref_file.stem
            if stem.startswith('master_fskt_ref') and not '_sample_counts' in stem:
                try:
                    ref_part = stem.replace('master_fskt_ref', '')
                    ref_num = int(ref_part)
                    if ref_num >= 1:  # Start from ref1
                        available_refs.append(ref_num)
                except ValueError:
                    continue
    
    available_refs = sorted(available_refs)
    print(f"Found reference times: {available_refs} (corresponding to tw = {[r*200 for r in available_refs]} ps)")
    
    # Collect relaxation time data
    relaxation_data = {
        'actual': {},
        'fictive': {}
    }
    
    for coupling_val in coupling_values:
        print(f"\n  Processing ε = {coupling_val}...")
        relaxation_data['actual'][coupling_val] = {}
        relaxation_data['fictive'][coupling_val] = {}
        
        coupling_data = all_data[coupling_val]
        
        for ref_num in available_refs:
            tw = ref_num * 200.0
            
            try:
                # Generate fictive F(k,t)
                time_lag_fictive, fkt_fictive = generate_fkt_from_fictive(
                    coupling_data['time_fictive_shifted'] + 200.0,
                    coupling_data['xi_fictive_scaled'],
                    tw=tw,
                    beta=beta
                )
                
                # Load actual F(k,t)
                time_lag_actual, fkt_actual_raw = load_actual_fkt_data(coupling_val, ref_num)
                fkt_actual = fkt_actual_raw / fkt_actual_raw[0]  # Normalize
                
                # Find relaxation times
                tau_fictive = find_relaxation_time_01_criterion(time_lag_fictive, fkt_fictive)
                tau_actual = find_relaxation_time_01_criterion(time_lag_actual, fkt_actual)
                
                relaxation_data['fictive'][coupling_val][ref_num] = tau_fictive
                relaxation_data['actual'][coupling_val][ref_num] = tau_actual
                
                if not np.isnan(tau_actual) and not np.isnan(tau_fictive):
                    print(f"    ref{ref_num} (tw={tw:.0f}ps): τ_actual={tau_actual:.1f}ps, τ_fictive={tau_fictive:.1f}ps")
                
            except Exception as e:
                print(f"    Warning: Could not process ref{ref_num}: {e}")
                relaxation_data['fictive'][coupling_val][ref_num] = np.nan
                relaxation_data['actual'][coupling_val][ref_num] = np.nan
                continue
    
    # Create the comparison plot (similar to normalized_relaxation_analysis_filtered.png)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    # Color scheme for coupling strengths
    from matplotlib.colors import Normalize
    coupling_norm = Normalize(vmin=0, vmax=max(coupling_values))
    coupling_cmap = plt.colormaps.get_cmap('coolwarm')
    
    def get_coupling_color(coupling_value):
        return coupling_cmap(coupling_norm(coupling_value))
    
    # Panel 1: Relaxation time vs coupling strength for different references
    print("\n  Creating Panel 1: Relaxation time vs coupling strength...")
    
    for ref_num in available_refs[::2]:  # Use every other reference to avoid overcrowding
        ref_time_ps = ref_num * 200
        tau_actual_vals = []
        tau_fictive_vals = []
        coupling_vals_plot = []
        
        for coupling_val in coupling_values:
            tau_actual = relaxation_data['actual'][coupling_val].get(ref_num, np.nan)
            tau_fictive = relaxation_data['fictive'][coupling_val].get(ref_num, np.nan)
            
            if not np.isnan(tau_actual) and not np.isnan(tau_fictive):
                tau_actual_vals.append(tau_actual)
                tau_fictive_vals.append(tau_fictive)
                coupling_vals_plot.append(coupling_val)
        
        if len(tau_actual_vals) > 0:
            color = plt.colormaps.get_cmap('viridis')(ref_num / max(available_refs))
            
            ax1.plot(coupling_vals_plot, tau_actual_vals, 'o-', color=color, linewidth=2, 
                    markersize=6, alpha=0.8, label=f'Actual, tw={ref_time_ps}ps')
            ax1.plot(coupling_vals_plot, tau_fictive_vals, 's--', color=color, linewidth=2, 
                    markersize=6, alpha=0.8, label=f'Fictive, tw={ref_time_ps}ps')
    
    ax1.set_xlabel('Coupling Strength ε')
    ax1.set_ylabel('Relaxation Time τ (ps) [F(k,t) = 0.1]')
    ax1.set_title('Relaxation Time vs Coupling Strength\nActual vs Fictive-based F(k,t)')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8, ncol=2)
    
    # Panel 2: Relaxation time vs reference time for different couplings
    print("  Creating Panel 2: Relaxation time vs reference time...")
    
    for coupling_val in coupling_values:
        ref_times = []
        tau_actual_vals = []
        tau_fictive_vals = []
        
        for ref_num in available_refs:
            tau_actual = relaxation_data['actual'][coupling_val].get(ref_num, np.nan)
            tau_fictive = relaxation_data['fictive'][coupling_val].get(ref_num, np.nan)
            
            if not np.isnan(tau_actual) and not np.isnan(tau_fictive):
                ref_times.append(ref_num * 200)
                tau_actual_vals.append(tau_actual)
                tau_fictive_vals.append(tau_fictive)
        
        if len(tau_actual_vals) > 0:
            color = get_coupling_color(coupling_val)
            coupling_label = f'ε = {coupling_val:.4f}' if coupling_val > 0 else 'ε = 0.0'
            
            ax2.plot(ref_times, tau_actual_vals, 'o-', color=color, linewidth=2, 
                    markersize=6, alpha=0.8, label=f'{coupling_label} (actual)')
            ax2.plot(ref_times, tau_fictive_vals, 's--', color=color, linewidth=2, 
                    markersize=6, alpha=0.8, label=f'{coupling_label} (fictive)')
    
    # Add vertical line at coupling turn-on
    ax2.axvline(x=200, color='red', linestyle=':', linewidth=2, alpha=0.8,
                label='Coupling Turn-On (tw = 200 ps)')
    
    ax2.set_xlabel('Reference Time tw (ps)')
    ax2.set_ylabel('Relaxation Time τ (ps) [F(k,t) = 0.1]')
    ax2.set_title('Relaxation Time vs Reference Time\nActual vs Fictive-based F(k,t)')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8, ncol=2)
    ax2.set_xlim(left=0)
    
    plt.suptitle(f'Relaxation Time Analysis: F(k,t) = 0.1 Criterion\nActual vs Fictive-based F(k,t), β = {beta}', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save plot
    output_file = Path(output_dir) / 'relaxation_time_comparison_actual_vs_fictive.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved relaxation time comparison: {output_file}")
    
    # Also save as PDF
    output_file_pdf = Path(output_dir) / 'relaxation_time_comparison_actual_vs_fictive.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"Saved: {output_file_pdf}")
    
    plt.show()
    
    # Calculate and report statistics
    print(f"\n" + "-" * 60)
    print("RELAXATION TIME COMPARISON STATISTICS")
    print("-" * 60)
    
    all_tau_actual = []
    all_tau_fictive = []
    
    for coupling_val in coupling_values:
        for ref_num in available_refs:
            tau_actual = relaxation_data['actual'][coupling_val].get(ref_num, np.nan)
            tau_fictive = relaxation_data['fictive'][coupling_val].get(ref_num, np.nan)
            
            if not np.isnan(tau_actual) and not np.isnan(tau_fictive):
                all_tau_actual.append(tau_actual)
                all_tau_fictive.append(tau_fictive)
    
    if len(all_tau_actual) > 0:
        correlation = np.corrcoef(all_tau_actual, all_tau_fictive)[0, 1]
        rel_error = np.mean(np.abs(np.array(all_tau_actual) - np.array(all_tau_fictive)) / np.array(all_tau_actual)) * 100
        
        print(f"Overall correlation coefficient: {correlation:.6f}")
        print(f"Mean relative error: {rel_error:.2f}%")
        print(f"Total data points compared: {len(all_tau_actual)}")
        
        if correlation > 0.99:
            print(" EXCELLENT agreement between actual and fictive relaxation times")
        elif correlation > 0.95:
            print(" GOOD agreement between actual and fictive relaxation times")
        else:
            print("  MODERATE agreement between actual and fictive relaxation times")
    
    return output_file.name


def main():
    """Main function to run both overlay and F(k,t) comparison."""
    print("=" * 80)
    print("MATERIAL TIME OVERLAY AND F(k,t) COMPARISON ANALYSIS")
    print("=" * 80)
    
    # Step 1: Create material time overlay
    print("\nStep 1: Creating material time overlay...")
    all_data = create_shifted_overlay_plot()
    
    # Step 2: Compare F(k,t) from fictive vs actual data for all coupling strengths
    print("\nStep 2: Comparing F(k,t) from fictive material time vs actual data...")
    coupling_values_to_test = list(all_data.keys())
    print(f"Will generate individual F(k,t) comparisons for {len(coupling_values_to_test)} coupling strengths: {coupling_values_to_test}")
    
    individual_files = []
    for coupling_val in coupling_values_to_test:
        print(f"\n  -> Processing coupling ε = {coupling_val}...")
        compare_fkt_fictive_vs_actual(all_data, coupling_value=coupling_val, ref_num=1, beta=0.55)
        coupling_str = f"eps{coupling_val:.4f}".replace(".", "p")
        individual_files.append(f"fkt_comparison_fictive_vs_actual_{coupling_str}_ref1.png")
    
    # Step 3: Create all-coupling overlay plots for all reference times
    print("\nStep 3: Creating all-coupling overlay plots for all reference times...")
    overlay_files = create_all_coupling_overlay_plots(all_data, beta=0.55)
    
    # Step 4: Compare relaxation times (F(k,t) = 0.1 criterion)
    print("\nStep 4: Comparing relaxation times from actual vs fictive-based F(k,t)...")
    relaxation_comparison_file = compare_relaxation_times_analysis(all_data, beta=0.55)
    
    print(f"\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print("Generated files:")
    print("   material_time_shifted_overlay.png - Material time overlay")
    print(f"   {relaxation_comparison_file} - Relaxation time comparison (F(k,t) = 0.1)")
    print("\nIndividual coupling comparisons (ref1 only):")
    for filename in individual_files:
        print(f"   {filename}")
    print("\nAll-coupling overlay plots (all reference times):")
    for filename in overlay_files:
        print(f"   {filename}")
    print("   PDF versions of all plots")


if __name__ == "__main__":
    main()
