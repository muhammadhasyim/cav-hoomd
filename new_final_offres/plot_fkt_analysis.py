#!/usr/bin/env python3
"""
F(k,t) Analysis and Plotting Script

This script analyzes F(k,t) data from averaged replica files and creates three types of plots:
1. Multi-panel plot: F(k,t) vs time for different ref's (one panel per coupling strength)
2. Multi-panel plot: F(k,t) vs time for different coupling strengths (one panel per ref)
3. Two-panel plot: Relaxation times analysis

Usage: python plot_fkt_analysis.py [--base_dir ./]
"""

import argparse
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for remote servers
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import re
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import subprocess

def coupling_to_frequency(coupling):
    """Convert coupling strength to actual target frequency in cm⁻¹."""
    # Direct mapping from actual coupling values (from directories) to target frequencies
    # Based on constant ratio ε/ω² = 1e-3/(1560)² and available simulation folders
    coupling_to_freq_map = {
        # Actual coupling → Target frequency (8 frequencies total)
        6.79487e-04: 1060.00000000,        # coupling_679487e-09_switch_200.0ps
        7.70513e-04: 1202.85714286,        # coupling_770513e-09_switch_200.0ps
        8.62820e-04: 1345.71428571,        # coupling_862820e-09_switch_200.0ps
        9.54487e-04: 1488.57142857,        # coupling_954487e-09_switch_200.0ps
        1.04551e-03: 1631.42857143,        # coupling_104551e-08_switch_200.0ps
        1.13718e-03: 1774.28571429,        # coupling_113718e-08_switch_200.0ps
        1.22885e-03: 1917.14285714,        # coupling_122885e-08_switch_200.0ps
        1.32051e-03: 2060.00000000,        # coupling_132051e-08_switch_200.0ps
    }
    
    # Find the closest coupling value in our map
    if coupling <= 0:
        return 0.0
    
    # Find the closest match
    best_match = min(coupling_to_freq_map.keys(), key=lambda x: abs(x - coupling))
    
    # Only return the mapped frequency if the coupling is reasonably close
    if abs(best_match - coupling) / best_match < 0.5:  # Within 50%
        return coupling_to_freq_map[best_match]
    
    # Fallback to approximate calculation if no close match
    CM1_TO_AU = 0.000004556335
    REF_FREQ_CM1 = 1560.0
    REF_COUPLING = 0.001
    REF_FREQ_AU = REF_FREQ_CM1 * CM1_TO_AU
    RATIO_CONSTANT = REF_COUPLING / (REF_FREQ_AU * REF_FREQ_AU)
    
    freq_au_squared = coupling / RATIO_CONSTANT
    freq_au = (freq_au_squared) ** 0.5
    freq_cm1 = freq_au / CM1_TO_AU
    return freq_cm1

def parse_coupling_strength(folder_name, use_latex=False):
    """Parse coupling strength from folder name."""
    # Handle both cavity_coupling_* and coupling_* folder formats
    if 'cavity_coupling_' in folder_name:
        coupling_part = folder_name.replace('cavity_coupling_', '')
    elif 'coupling_' in folder_name:
        coupling_part = folder_name.replace('coupling_', '')
    else:
        return float('nan'), folder_name
    
    # Remove the switch part
    coupling_part = coupling_part.split('_switch_')[0]
    
    if '0epos00' in coupling_part:
        return 0.0, '0.0'
    elif 'eneg' in coupling_part:
        parts = coupling_part.split('eneg')
        if len(parts) == 2:
            mantissa = int(parts[0])
            exponent = int(parts[1][1])
            value = mantissa * (10 ** (-exponent))
            
            # Generate appropriate label format
            if use_latex:
                label = f'{mantissa}\\times10^{{-{exponent}}}'
            else:
                label = f'{mantissa}×10⁻{exponent}'
            
            return value, label
    elif 'e-' in coupling_part:
        # Handle scientific notation like 109310e-08
        try:
            mantissa_str, exponent_str = coupling_part.split('e-')
            mantissa_int = int(mantissa_str)
            exponent = int(exponent_str)
            
            # Convert format like 109310e-08 to 1.09310e-03
            mantissa_decimal = mantissa_int / (10**(len(mantissa_str)-1))
            actual_exponent = -(exponent - (len(mantissa_str)-1))
            value = mantissa_decimal * (10 ** actual_exponent)
            
            if use_latex:
                label = f'{mantissa_decimal:.5f}\\times10^{{-{abs(actual_exponent)}}}'
            else:
                label = f'{mantissa_decimal:.5f}×10⁻{abs(actual_exponent)}'
            
            return value, label
        except Exception:
            pass
    
    return float('nan'), coupling_part

def filter_coupling_strengths(coupling_info):
    """Filter to only include specific coupling strengths corresponding to target frequencies."""
    # Target couplings for the 8 frequencies (excluding switch_1.0ps folders)
    target_couplings = [
        6.79487e-04,  # 1060.00 cm⁻¹
        7.70513e-04,  # 1202.86 cm⁻¹ 
        8.62820e-04,  # 1345.71 cm⁻¹
        9.54487e-04,  # 1488.57 cm⁻¹
        1.04551e-03,  # 1631.43 cm⁻¹
        1.13718e-03,  # 1774.29 cm⁻¹
        1.22885e-03,  # 1917.14 cm⁻¹
        1.32051e-03,  # 2060.00 cm⁻¹
    ]
    
    filtered_couplings = {}
    for coupling_name, (coupling_value, coupling_label) in coupling_info.items():
        # Skip folders with switch_1.0ps
        if 'switch_1.0ps' in coupling_name:
            continue
            
        if not np.isnan(coupling_value):
            # Check if coupling is close to any target (within 15% tolerance)
            for target in target_couplings:
                if abs(coupling_value - target) / target < 0.15:
                    filtered_couplings[coupling_name] = (coupling_value, coupling_label)
                    break
    
    return filtered_couplings

def process_fkt_data(time, fkt, normalization_value=None):
    """
    Process F(k,t) data: normalize and remove values below 0.001.
    Also auto-detect maximum time for plotting.
    """
    time = np.array(time, dtype=np.float64)
    fkt = np.array(fkt, dtype=np.float64)

    # Remove NaN values and zero values
    valid_mask = ~(np.isnan(time) | np.isnan(fkt)) & (fkt != 0)
    time_clean = time[valid_mask]
    fkt_clean = fkt[valid_mask]
    
    if len(time_clean) < 2:
        return None, None, None
    
    # Sort by time to ensure monotonic behavior
    sort_indices = np.argsort(time_clean)
    time_sorted = time_clean[sort_indices]
    fkt_sorted = fkt_clean[sort_indices]
    
    # Normalize using provided normalization value
    if normalization_value is not None and normalization_value != 0:
        fkt_normalized = fkt_sorted / normalization_value
    else:
        # Fallback: normalize by first value
        if len(fkt_sorted) > 0:
            fkt_normalized = fkt_sorted / fkt_sorted[0]
        else:
            return None, None, None
    
    # Remove values below 0.001 after normalization
    above_threshold_mask = fkt_normalized >= 0.001
    if not np.any(above_threshold_mask):
        return None, None, None
    
    time_filtered = time_sorted[above_threshold_mask]
    fkt_filtered = fkt_normalized[above_threshold_mask]
    
    # Auto-detect maximum time for plotting
    # Use the last time point where F(k,t) >= 0.001
    max_time = time_filtered[-1] if len(time_filtered) > 0 else None
    
    return time_filtered, fkt_filtered, max_time

def read_fkt_data(file_path):
    """Read F(k,t) data from averaged replica file."""
    try:
        # Read with header=0 to properly handle the header line
        data = pd.read_csv(file_path, sep=r'\s+', comment='#', header=0)
        
        # Debug: Print first few values to check data types
        # print(f"Debug - file: {file_path}")
        # print(f"Debug - columns: {data.columns.tolist()}")
        # print(f"Debug - first 3 rows:\n{data.head(3)}")
        
        # Use column names or indices to access the data
        if 'lag_time' in data.columns and 'fskt' in data.columns:
            time = data['lag_time'].values
            fkt = data['fskt'].values
        else:
            # Fallback to positional indexing if column names are unexpected
            time = data.iloc[:, 0].values
            fkt = data.iloc[:, 1].values
        
        # Additional check: ensure data types are numeric
        if len(time) > 0:
            # Convert to numeric, forcing errors to NaN (which can be filtered out later)
            time = pd.to_numeric(time, errors='coerce')
            fkt = pd.to_numeric(fkt, errors='coerce')
            
            # Remove any NaN values that resulted from conversion errors
            valid_mask = ~(pd.isna(time) | pd.isna(fkt))
            time = time[valid_mask]
            fkt = fkt[valid_mask]
        
        return time, fkt
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None, None

def stretched_exponential(t, A, tau, beta):
    """
    Stretched exponential (Kohlrausch-Williams-Watts) function.
    
    Parameters:
    t : array-like, time
    A : float, amplitude (pre-exponential factor)
    tau : float, relaxation time
    beta : float, stretching exponent (0 < beta <= 1)
    
    Returns:
    F(t) = A * exp(-(t/tau)^beta)
    """
    return A * np.exp(-np.power(t / tau, beta))

def stretched_exponential_offset(t, A, tau, beta, t0):
    """
    Stretched exponential function with time offset for fitting decay region.
    
    Parameters:
    t : array-like, time
    A : float, amplitude at t=t0
    tau : float, relaxation time
    beta : float, stretching exponent (0 < beta <= 1)
    t0 : float, time offset
    
    Returns:
    F(t) = A * exp(-((t-t0)/tau)^beta) for t >= t0
    """
    # Only apply for t >= t0
    result = np.zeros_like(t)
    mask = t >= t0
    result[mask] = A * np.exp(-np.power((t[mask] - t0) / tau, beta))
    return result

def fit_stretched_exponential(time, fkt, normalization_value=None):
    """
    Fit F(k,t) data to a stretched exponential function.
    
    Parameters:
    time : array-like, time points
    fkt : array-like, F(k,t) values
    normalization_value : float, optional, value to normalize F(k,t) to at t=0
    
    Returns:
    dict with fitting results: {'tau': relaxation_time, 'beta': stretch_exponent, 
                               'A': amplitude, 'r_squared': goodness_of_fit,
                               'success': bool, 'fit_curve': fitted_values}
    """
    try:
        # Clean and sort data
        valid_mask = ~(np.isnan(time) | np.isnan(fkt)) & (fkt != 0) & (time >= 0)
        time_clean = time[valid_mask]
        fkt_clean = fkt[valid_mask]
        
        if len(time_clean) < 10:  # Need sufficient points for reliable fitting
            return {'success': False, 'tau': np.nan, 'beta': np.nan, 'A': np.nan, 
                   'r_squared': np.nan, 'fit_curve': None}
        
        # Sort by time
        sort_indices = np.argsort(time_clean)
        time_sorted = time_clean[sort_indices]
        fkt_sorted = fkt_clean[sort_indices]
        
        # Normalize if requested
        if normalization_value is not None and normalization_value != 0:
            fkt_sorted = fkt_sorted / normalization_value
        
        # Focus on the decay region from 0.5 to 0.0 for more robust fitting
        decay_mask = (fkt_sorted <= 0.5) & (fkt_sorted > 0.0)
        if np.sum(decay_mask) < 10:
            return {'success': False, 'tau': np.nan, 'beta': np.nan, 'A': np.nan, 
                   'r_squared': np.nan, 'fit_curve': None}
        
        time_fit = time_sorted[decay_mask]
        fkt_fit = fkt_sorted[decay_mask]
        
        # Initial parameter guesses
        # Since we're fitting the decay region (0.5 to 0.0), use offset fitting
        A_guess = fkt_fit[0] if len(fkt_fit) > 0 else 0.5  # Amplitude at first point
        t0_guess = time_fit[0] if len(time_fit) > 0 else 0.0  # Time offset
        
        # Estimate tau from the time span of decay
        time_span = time_fit[-1] - time_fit[0]
        tau_guess = time_span / 2  # Rough estimate
        
        beta_guess = 1.0  # Start with simple exponential
        
        # Parameter bounds: A>0, tau>0, 0<beta<=1, t0 around the start time
        bounds = ([0.05, 1e-6, 0.1, time_fit[0] - 100], 
                 [1.0, np.inf, 1.0, time_fit[0] + 100])
        
        # Fit the stretched exponential with offset
        popt, pcov = curve_fit(
            stretched_exponential_offset, 
            time_fit, 
            fkt_fit,
            p0=[A_guess, tau_guess, beta_guess, t0_guess],
            bounds=bounds,
            maxfev=5000
        )
        
        A_fit, tau_fit, beta_fit, t0_fit = popt
        
        # Calculate goodness of fit (R²)
        fkt_predicted = stretched_exponential_offset(time_fit, A_fit, tau_fit, beta_fit, t0_fit)
        ss_res = np.sum((fkt_fit - fkt_predicted) ** 2)
        ss_tot = np.sum((fkt_fit - np.mean(fkt_fit)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        # Generate fit curve for plotting
        time_extended = np.linspace(time_fit[0], time_fit[-1], 200)
        fit_curve = {
            'time': time_extended,
            'fkt': stretched_exponential_offset(time_extended, A_fit, tau_fit, beta_fit, t0_fit)
        }
        
        return {
            'success': True,
            'tau': tau_fit,
            'beta': beta_fit,
            'A': A_fit,
            'r_squared': r_squared,
            'fit_curve': fit_curve,
            'fitted_time': time_fit,
            'fitted_fkt': fkt_fit,
            'predicted_fkt': fkt_predicted
        }
        
    except Exception as e:
        print(f"Stretched exponential fitting failed: {e}")
        return {'success': False, 'tau': np.nan, 'beta': np.nan, 'A': np.nan, 
               'r_squared': np.nan, 'fit_curve': None}

def find_relaxation_time(time, fkt, target_value=0.1, normalization_value=None, use_fitting=False):
    """Find the relaxation time using direct threshold criterion (default: F(k,t) = 0.1)."""
    try:
        # Always use direct threshold method (disabled fitting)
        # Remove NaN values and zero values (zero values are not meaningful for relaxation analysis)
        valid_mask = ~(np.isnan(time) | np.isnan(fkt)) & (fkt != 0)
        time_clean = time[valid_mask]
        fkt_clean = fkt[valid_mask]
        
        if len(time_clean) < 2:
            return np.nan
        
        # Sort by time to ensure monotonic behavior
        sort_indices = np.argsort(time_clean)
        time_sorted = time_clean[sort_indices]
        fkt_sorted = fkt_clean[sort_indices]
        
        # Don't artificially truncate - use all available non-zero data
        # Only remove last few points if they seem like noise (less aggressive)
        if len(fkt_sorted) > 10:
            n_keep = int(0.995 * len(fkt_sorted))  # Only remove 0.5% instead of 1%
            fkt_sorted = fkt_sorted[:n_keep]
            time_sorted = time_sorted[:n_keep]
        
        # Use provided normalization value
        if normalization_value is not None and normalization_value != 0:
            fkt_sorted = fkt_sorted / normalization_value
        else:
            # Fallback: normalize by first value (should not happen now since we have ref0 norm)
            if len(fkt_sorted) > 0:
                fkt_sorted = fkt_sorted / fkt_sorted[0]
            else:
                return np.nan
        
        # Check if target value is within the range of fkt
        if target_value < fkt_sorted.min() or target_value > fkt_sorted.max():
            return np.nan
        
        # Remove duplicate F(k,t) values by keeping the first occurrence
        _, unique_indices = np.unique(fkt_sorted, return_index=True)
        unique_indices = np.sort(unique_indices)
        fkt_unique = fkt_sorted[unique_indices]
        time_unique = time_sorted[unique_indices]
        
        if len(fkt_unique) < 2:
            return np.nan
        
        # Check again if target is in range after removing duplicates
        if target_value < fkt_unique.min() or target_value > fkt_unique.max():
            return np.nan
        
        # For F(k,t) data, we typically expect monotonic decay
        # If F(k,t) is not monotonic, we need to handle this carefully
        if not np.all(np.diff(fkt_unique) <= 0):  # Not monotonically decreasing
            # Find the part where F(k,t) crosses the target value (first crossing)
            crossing_indices = np.where((fkt_sorted[:-1] >= target_value) & (fkt_sorted[1:] <= target_value))[0]
            if len(crossing_indices) > 0:
                # Use linear interpolation for the first crossing
                idx = crossing_indices[0]
                t1, t2 = time_sorted[idx], time_sorted[idx+1]
                f1, f2 = fkt_sorted[idx], fkt_sorted[idx+1]
                
                if f1 != f2:  # Avoid division by zero
                    relaxation_time = t1 + (target_value - f1) * (t2 - t1) / (f2 - f1)
                    return relaxation_time if relaxation_time >= 0 else np.nan
                else:
                    return t1 if t1 >= 0 else np.nan
            else:
                return np.nan
        
        # If monotonic, use interpolation
        if len(fkt_unique) > 3:
            interp_func = interp1d(fkt_unique, time_unique, kind='cubic', 
                                 bounds_error=False, fill_value=np.nan)
        else:
            interp_func = interp1d(fkt_unique, time_unique, kind='linear',
                                 bounds_error=False, fill_value=np.nan)
        
        relaxation_time = interp_func(target_value)
        return relaxation_time if not np.isnan(relaxation_time) and relaxation_time >= 0 else np.nan
        
    except Exception as e:
        print(f"Error finding relaxation time: {e}")
        return np.nan

def collect_data(base_dir):
    """Collect all F(k,t) data from the coupling strength folders."""
    base_path = Path(base_dir)
    
    # Look for both cavity_coupling_* and coupling_* folders with switch_200.0ps
    cavity_coupling_folders = list(base_path.glob('cavity_coupling_*200.0ps'))
    coupling_folders = list(base_path.glob('coupling_*200.0ps'))
    
    # Combine both types of folders
    all_coupling_folders = cavity_coupling_folders + coupling_folders
    
    if not all_coupling_folders:
        print(f"No coupling folders found in {base_dir}")
        return {}, {}
    
    data_dict = {}
    coupling_info = {}
    
    for folder in all_coupling_folders:
        coupling_name = folder.name
        coupling_value, coupling_label = parse_coupling_strength(coupling_name)
        coupling_info[coupling_name] = (coupling_value, coupling_label)
        
        print(f"Processing folder: {coupling_name} (coupling = {coupling_label})")
        
        ref_files = list(folder.glob('master_fskt_ref*.txt'))
        
        if not ref_files:
            print(f"  No master_fskt_ref*.txt files found in {folder}")
            continue
        
        folder_data = {}
        
        for ref_file in ref_files:
            ref_match = re.search(r'master_fskt_ref(\d+)\.txt', ref_file.name)
            if ref_match:
                ref_number = int(ref_match.group(1))
                
                time, fkt = read_fkt_data(ref_file)
                
                if time is not None and fkt is not None:
                    folder_data[ref_number] = (time, fkt)
                    print(f"  Read ref{ref_number}: {len(time)} data points")
                else:
                    print(f"  Failed to read ref{ref_number}")
        
        if folder_data:
            data_dict[coupling_name] = folder_data
        
        print(f"  Total refs found: {len(folder_data)}")
    
    return data_dict, coupling_info

def plot_fkt_by_coupling(data_dict, coupling_info, output_dir='.'):
    """Create multi-panel plot: F(k,t) vs time for different ref's."""
    print("\nCreating F(k,t) by coupling strength plot...")
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    n_couplings = len(sorted_couplings)
    if n_couplings == 0:
        print("No target coupling data to plot")
        return
    
    if n_couplings <= 2:
        rows, cols = 1, n_couplings
    elif n_couplings <= 4:
        rows, cols = 2, 2
    else:
        rows = int(np.ceil(n_couplings / 3))
        cols = 3
    
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows))
    if n_couplings == 1:
        axes = [axes]
    elif rows == 1 or cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    # Track maximum time across all data for consistent x-axis
    global_max_time = 0
    
    for i, coupling_name in enumerate(sorted_couplings):
        if i >= len(axes):
            break
            
        ax = axes[i]
        
        # Check if this coupling has data
        if coupling_name not in data_dict:
            coupling_value, coupling_label = filtered_coupling_info[coupling_name]
            ax.text(0.5, 0.5, 'No data files found', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Coupling: {coupling_label}')
            continue
            
        folder_data = data_dict[coupling_name]
        coupling_value, coupling_label = filtered_coupling_info[coupling_name]
        
        ref_numbers = sorted(folder_data.keys())
        
        if not ref_numbers:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Coupling: {coupling_label}')
            continue
        
        # Get normalization value from ref0
        normalization_value = None
        if 0 in folder_data:
            time_ref0, fkt_ref0 = folder_data[0]
            # Find first non-zero value from ref0 for normalization
            nonzero_mask = fkt_ref0 != 0
            if np.any(nonzero_mask):
                first_nonzero_idx = np.where(nonzero_mask)[0][0]
                normalization_value = fkt_ref0[first_nonzero_idx]
        
        norm = Normalize(vmin=0, vmax=max(ref_numbers))
        cmap = plt.colormaps.get_cmap('coolwarm')
        
        panel_max_time = 0
        
        for ref_num in ref_numbers:
            time, fkt = folder_data[ref_num]
            color = cmap(norm(ref_num))
            
            # Process F(k,t) data with filtering and normalization
            time_processed, fkt_processed, max_time = process_fkt_data(time, fkt, normalization_value)
            
            if time_processed is not None and fkt_processed is not None:
                ax.plot(time_processed, fkt_processed, color=color, linewidth=2, alpha=0.8, 
                       label=f'ref{ref_num}')
                
                if max_time is not None:
                    panel_max_time = max(panel_max_time, max_time)
                    global_max_time = max(global_max_time, max_time)
        
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('F(k,t) (normalized, ≥0.001)')
        ax.set_title(f'Coupling: {coupling_label}')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, ncol=2)
        ax.set_ylim(bottom=0.001)  # Start y-axis at threshold
        
        # Set x-axis limit based on data for this panel
        if panel_max_time > 0:
            ax.set_xlim(0, panel_max_time * 1.05)  # Add 5% padding
    
    for i in range(n_couplings, len(axes)):
        fig.delaxes(axes[i])
    
    plt.suptitle('F(k,t) vs Time: Different References by Coupling Strength\n(8 Frequencies: 1060-2060 cm⁻¹)', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    output_file = Path(output_dir) / 'fkt_by_coupling_filtered.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_fkt_by_ref(data_dict, coupling_info, output_dir='.'):
    """Create multi-panel plot: F(k,t) vs time for different coupling strengths."""
    print("\nCreating F(k,t) by reference plot...")
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    all_refs = set()
    for coupling_name in filtered_coupling_info.keys():
        if coupling_name in data_dict:
            all_refs.update(data_dict[coupling_name].keys())
    
    all_refs = sorted(all_refs)
    
    if not all_refs:
        print("No reference data to plot")
        return
    
    n_refs = len(all_refs)
    
    if n_refs <= 2:
        rows, cols = 1, n_refs
    elif n_refs <= 4:
        rows, cols = 2, 2
    else:
        rows = int(np.ceil(n_refs / 3))
        cols = 3
    
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4*rows))
    if n_refs == 1:
        axes = [axes]
    elif rows == 1 or cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    coupling_values = [filtered_coupling_info[c][0] for c in sorted_couplings]
    norm = Normalize(vmin=min(coupling_values), vmax=max(coupling_values))
    cmap = plt.colormaps.get_cmap('coolwarm')
    
    # Track maximum time across all data for consistent x-axis
    global_max_time = 0
    
    for i, ref_num in enumerate(all_refs):
        if i >= len(axes):
            break
            
        ax = axes[i]
        has_data = False
        panel_max_time = 0
        
        for coupling_name in sorted_couplings:
            if coupling_name in data_dict and ref_num in data_dict[coupling_name]:
                time, fkt = data_dict[coupling_name][ref_num]
                coupling_value, coupling_label = filtered_coupling_info[coupling_name]
                
                color = cmap(norm(coupling_value))
                
                # Get normalization value from ref0 of this coupling
                normalization_value = None
                if 0 in data_dict[coupling_name]:
                    time_ref0, fkt_ref0 = data_dict[coupling_name][0]
                    # Find first non-zero value from ref0 for normalization
                    nonzero_mask = fkt_ref0 != 0
                    if np.any(nonzero_mask):
                        first_nonzero_idx = np.where(nonzero_mask)[0][0]
                        normalization_value = fkt_ref0[first_nonzero_idx]
                
                # Process F(k,t) data with filtering and normalization
                time_processed, fkt_processed, max_time = process_fkt_data(time, fkt, normalization_value)
                
                if time_processed is not None and fkt_processed is not None:
                    ax.plot(time_processed, fkt_processed, color=color, linewidth=2, alpha=0.8, 
                       label=f'{coupling_label}')
                    has_data = True
                    
                    if max_time is not None:
                        panel_max_time = max(panel_max_time, max_time)
                        global_max_time = max(global_max_time, max_time)
        
        if not has_data:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('F(k,t) (normalized, ≥0.001)')
        ax.set_title(f'ref{ref_num}')
        ax.grid(True, alpha=0.3)
        
        if has_data:
            ax.legend(fontsize=8)
            ax.set_ylim(bottom=0.001)  # Start y-axis at threshold
            
            # Set x-axis limit based on data for this panel
            if panel_max_time > 0:
                ax.set_xlim(0, panel_max_time * 1.05)  # Add 5% padding
    
    for i in range(n_refs, len(axes)):
        fig.delaxes(axes[i])
    
    plt.suptitle('F(k,t) vs Time: Different Coupling Strengths by Reference\n(8 Frequencies: 1060-2060 cm⁻¹)', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    output_file = Path(output_dir) / 'fkt_by_ref_filtered.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_fkt_diagnostic(data_dict, coupling_info, output_dir='.'):
    """Create diagnostic plot to understand F(k,t) behavior."""
    print("\nCreating F(k,t) diagnostic plot...")
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    n_couplings = len(sorted_couplings)
    if n_couplings == 0:
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    target_line = 0.1  # 0.1 criterion line for reference
    global_max_time = 0
    
    for i, coupling_name in enumerate(sorted_couplings[:4]):  # Limit to 4 for 2x2 grid
        if i >= len(axes):
            break
            
        ax = axes[i]
        
        # Check if this coupling has data
        if coupling_name not in data_dict:
            coupling_value, coupling_label = filtered_coupling_info[coupling_name]
            ax.text(0.5, 0.5, 'No data files found', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Coupling: {coupling_label}')
            continue
            
        folder_data = data_dict[coupling_name]
        coupling_value, coupling_label = filtered_coupling_info[coupling_name]
        
        if 0 in folder_data:  # Only plot ref0 for diagnostic
            time, fkt = folder_data[0]
            
            # Get normalization value from ref0 (same as current data since this is ref0)
            normalization_value = None
            nonzero_mask = fkt != 0
            if np.any(nonzero_mask):
                first_nonzero_idx = np.where(nonzero_mask)[0][0]
                normalization_value = fkt[first_nonzero_idx]
            
            # Process F(k,t) data with filtering and normalization
            time_processed, fkt_processed, max_time = process_fkt_data(time, fkt, normalization_value)
            
            if time_processed is not None and fkt_processed is not None:
                # Plot F(k,t) vs time
                ax.plot(time_processed, fkt_processed, 'b-', linewidth=2, alpha=0.8, label='F(k,t)')
                
                # Add 0.1 reference line (only if within range)
                if target_line >= 0.001:
                    ax.axhline(y=target_line, color='red', linestyle='--', alpha=0.7, 
                          label=f'F(k,t) = {target_line:.1f}')
                
                # Calculate and show relaxation time
                rel_time = find_relaxation_time(time, fkt, normalization_value=normalization_value)
                
                if not np.isnan(rel_time) and rel_time > 0 and target_line >= 0.001:
                    ax.axvline(x=rel_time, color='green', linestyle=':', alpha=0.7,
                              label=f'τ = {rel_time:.1f} ps')
                    ax.plot(rel_time, target_line, 'go', markersize=8, alpha=0.8)
                
                if max_time is not None:
                    global_max_time = max(global_max_time, max_time)
                    ax.set_xlim(0, max_time * 1.05)  # Add 5% padding
            else:
                # If all zeros, show message
                ax.text(0.5, 0.5, 'No valid data\n(all zeros or below threshold)', 
                       ha='center', va='center', transform=ax.transAxes)
            
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('F(k,t) (normalized, ≥0.001)')
            ax.set_title(f'Coupling: {coupling_label}')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=9)
            ax.set_ylim(bottom=0.001)  # Start y-axis at threshold
            
            # Add text with range info
            if time_processed is not None and fkt_processed is not None:
                fkt_min, fkt_max = fkt_processed.min(), fkt_processed.max()
                info_text = f'Range: {fkt_min:.3f} - {fkt_max:.3f}'
            ax.text(0.02, 0.98, info_text, transform=ax.transAxes, 
                   verticalalignment='top', fontsize=8,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        else:
            ax.text(0.5, 0.5, 'No ref0 data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Coupling: {coupling_label}')
    
    # Remove unused subplots
    for i in range(len(sorted_couplings), len(axes)):
        fig.delaxes(axes[i])
    
    plt.suptitle('F(k,t) Diagnostic: Understanding Relaxation Behavior\n(8 Frequencies: 1060-2060 cm⁻¹)', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_file = Path(output_dir) / 'fkt_diagnostic_filtered.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_relaxation_analysis(data_dict, coupling_info, output_dir='.'):
    """Create two-panel plot: Relaxation time analysis."""
    print("\nCreating relaxation time analysis plot...")
    
    # Set matplotlib style to classic and try to enable LaTeX
    plt.style.use('classic')
    try:
        # Test if LaTeX is available by trying to compile a simple expression
        subprocess.check_output(['latex', '--version'], stderr=subprocess.DEVNULL)
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        use_latex = True
        print("  Using LaTeX for mathematical typesetting")
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
            return fallback_text if fallback_text else text.replace('$', '').replace(r'\tau', 'τ').replace(r'\times', '×')
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    # Regenerate coupling labels with LaTeX format if needed
    if use_latex:
        latex_coupling_info = {}
        for coupling_name, (coupling_value, _) in filtered_coupling_info.items():
            _, latex_label = parse_coupling_strength(coupling_name, use_latex=True)
            latex_coupling_info[coupling_name] = (coupling_value, latex_label)
        filtered_coupling_info = latex_coupling_info
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    # Determine all available references across all couplings
    all_refs = set()
    for coupling_name in sorted_couplings:
        if coupling_name in data_dict:
            all_refs.update(data_dict[coupling_name].keys())
    all_refs = sorted(all_refs)
    
    # First, compute ref1 (200 ps) relaxation times for each coupling for normalization
    print("  Computing reference relaxation times (ref1, t_w = 200 ps) for normalization...")
    ref1_relaxation_times = {}
    
    for coupling_name in sorted_couplings:
        if coupling_name in data_dict and 1 in data_dict[coupling_name]:
            time_ref1, fkt_ref1 = data_dict[coupling_name][1]
            
            # Get normalization value from ref0 (for F(k,t) normalization)
            normalization_value = None
            if 0 in data_dict[coupling_name]:
                time_ref0, fkt_ref0 = data_dict[coupling_name][0]
                nonzero_mask = fkt_ref0 != 0
                if np.any(nonzero_mask):
                    first_nonzero_idx = np.where(nonzero_mask)[0][0]
                    normalization_value = fkt_ref0[first_nonzero_idx]
            
            tau_ref1 = find_relaxation_time(time_ref1, fkt_ref1, normalization_value=normalization_value)
            
            if not np.isnan(tau_ref1):
                ref1_relaxation_times[coupling_name] = tau_ref1
                coupling_value, _ = filtered_coupling_info[coupling_name]
                freq_value = coupling_to_frequency(coupling_value)
                print(f"    Coupling {freq_value:.0f} cm⁻¹: τ₀ = {tau_ref1:.1f} ps")
    
    # Panel 1: Relaxation time vs frequency for different references
    print("  Calculating normalized relaxation times vs frequency for different references...")
    
    coupling_values = [filtered_coupling_info[c][0] for c in sorted_couplings if c in data_dict]
    coupling_labels = [filtered_coupling_info[c][1] for c in sorted_couplings if c in data_dict]
    
    # Convert coupling values to frequencies
    frequency_values = [coupling_to_frequency(cv) for cv in coupling_values]
    frequency_labels = [f'{fv:.0f}' for fv in frequency_values]
    
    valid_couplings = [c for c in sorted_couplings if c in data_dict]
    
    # Create color mapping for reference times (Panel 1)
    ref_times_ps = [ref_num * 200 for ref_num in all_refs]
    ref_time_norm = Normalize(vmin=0, vmax=max(ref_times_ps))
    ref_time_cmap = plt.colormaps.get_cmap('viridis')
    
    # Store data for Panel 1 colorbar
    panel1_scatter_points = []
    
    for i, ref_num in enumerate(all_refs):
        ref_relaxation_times = []
        ref_frequency_values = []
        ref_frequency_labels = []
        
        for coupling_name in valid_couplings:
            if ref_num in data_dict[coupling_name]:
                time, fkt = data_dict[coupling_name][ref_num]
                coupling_value, coupling_label = filtered_coupling_info[coupling_name]
                
                # Get normalization value from ref0 of this coupling
                normalization_value = None
                if 0 in data_dict[coupling_name]:
                    time_ref0, fkt_ref0 = data_dict[coupling_name][0]
                    nonzero_mask = fkt_ref0 != 0
                    if np.any(nonzero_mask):
                        first_nonzero_idx = np.where(nonzero_mask)[0][0]
                        normalization_value = fkt_ref0[first_nonzero_idx]
                
                rel_time = find_relaxation_time(time, fkt, normalization_value=normalization_value)
                
                if not np.isnan(rel_time):  # Only include valid relaxation times
                    frequency_value = coupling_to_frequency(coupling_value)
                    
                    # Normalize by ref1 (200 ps) relaxation time if available
                    if coupling_name in ref1_relaxation_times and ref1_relaxation_times[coupling_name] > 0:
                        rel_time_normalized = rel_time / ref1_relaxation_times[coupling_name]
                    else:
                        rel_time_normalized = rel_time  # Fallback to absolute value
                    
                    ref_relaxation_times.append(rel_time_normalized)
                    ref_frequency_values.append(frequency_value)
                    ref_frequency_labels.append(f'{frequency_value:.0f}')
        
        if ref_relaxation_times:  # Only plot if we have data
            ref_time_ps = ref_num * 200
            color = ref_time_cmap(ref_time_norm(ref_time_ps))
            
            scatter = ax1.scatter(ref_frequency_values, ref_relaxation_times, 
                               c=[ref_time_ps]*len(ref_relaxation_times), 
                               cmap='viridis', norm=ref_time_norm,
                               s=80, alpha=0.8, edgecolors='black', linewidth=0.5,
                               label=latex_format(f'$t = {ref_time_ps}$ ps', f't = {ref_time_ps} ps'))
            
            # Connect points with lines
            ax1.plot(ref_frequency_values, ref_relaxation_times, '-', 
                    color=color, linewidth=2, alpha=0.6)
            
            panel1_scatter_points.append(scatter)
            
            print(f"    ref{ref_num} (t={ref_time_ps} ps): {len(ref_relaxation_times)} valid frequency points")
    
    ax1.set_xlabel(latex_format(r'Frequency (cm$^{-1}$)', 'Frequency (cm⁻¹)'))
    ax1.set_ylabel(latex_format(r'Normalized Relaxation Time $\tau/\tau_0$', 'Normalized Relaxation Time τ/τ₀'))
    ax1.set_title(latex_format(r'Normalized Relaxation Time vs Frequency' + '\n' + r'(Direct Threshold: $F(k,t) = 0.1$)', 
                              'Normalized Relaxation Time vs Frequency\n(Direct Threshold: F(k,t) = 0.1)'))
    ax1.grid(True, alpha=0.3)
    
    if frequency_values:
        ax1.set_xticks(frequency_values)
        ax1.set_xticklabels([f'{fv:.0f}' for fv in frequency_values], rotation=45)
    
    # Set x-axis range for Panel 1 (frequencies)
    if frequency_values:
        min_freq = min(frequency_values)
        max_freq = max(frequency_values)
        freq_range = max_freq - min_freq
        ax1.set_xlim(left=min_freq - 0.05*freq_range, right=max_freq + 0.05*freq_range)
    
    # Add colorbar for Panel 1 (reference times)
    if panel1_scatter_points:
        cbar1 = plt.colorbar(panel1_scatter_points[0], ax=ax1, shrink=0.8)
        cbar1.set_label(latex_format(r'Reference Time $t$ (ps)', 'Reference Time t (ps)'), rotation=270, labelpad=20)
    
    # Panel 2: Relaxation time vs reference time for different coupling strengths
    print("  Calculating normalized relaxation times vs reference time for different couplings...")
    
    # Create color mapping for coupling strengths (Panel 2)
    coupling_values_sorted = sorted([filtered_coupling_info[c][0] for c in valid_couplings])
    coupling_norm = Normalize(vmin=0, vmax=max(coupling_values_sorted))
    coupling_cmap = plt.colormaps.get_cmap('coolwarm')
    
    # Store data for Panel 2 colorbar
    panel2_scatter_points = []
    
    for i, coupling_name in enumerate(valid_couplings):
        coupling_value, coupling_label = filtered_coupling_info[coupling_name]
        folder_data = data_dict[coupling_name]
        
        # Get normalization value from ref0 of this coupling
        normalization_value = None
        if 0 in folder_data:
            time_ref0, fkt_ref0 = folder_data[0]
            nonzero_mask = fkt_ref0 != 0
            if np.any(nonzero_mask):
                first_nonzero_idx = np.where(nonzero_mask)[0][0]
                normalization_value = fkt_ref0[first_nonzero_idx]
        
        ref_numbers = sorted(folder_data.keys())
        ref_times = []
        relaxation_times = []
        
        for ref_num in ref_numbers:
            time, fkt = folder_data[ref_num]
            ref_time_ps = ref_num * 200
            
            rel_time = find_relaxation_time(time, fkt, normalization_value=normalization_value)
            
            if not np.isnan(rel_time):  # Only include valid relaxation times
                # Normalize by ref1 (200 ps) relaxation time if available
                if coupling_name in ref1_relaxation_times and ref1_relaxation_times[coupling_name] > 0:
                    rel_time_normalized = rel_time / ref1_relaxation_times[coupling_name]
                else:
                    rel_time_normalized = rel_time  # Fallback to absolute value
                
                ref_times.append(ref_time_ps)
                relaxation_times.append(rel_time_normalized)
        
        if relaxation_times:  # Only plot if we have data
            color = coupling_cmap(coupling_norm(coupling_value))
            frequency_value = coupling_to_frequency(coupling_value)
            
            scatter = ax2.scatter(ref_times, relaxation_times,
                               c=[coupling_value]*len(relaxation_times),
                               cmap='coolwarm', norm=coupling_norm,
                               s=80, alpha=0.8, edgecolors='black', linewidth=0.5,
                               label=latex_format(f'$\\nu = {frequency_value:.0f}$ cm$^{{-1}}$', 
                                                 f'ν = {frequency_value:.0f} cm⁻¹'))
            
            # Connect points with lines
            ax2.plot(ref_times, relaxation_times, '-',
                    color=color, linewidth=2, alpha=0.6)
            
            panel2_scatter_points.append(scatter)
            
            print(f"    Coupling {frequency_value:.0f} cm⁻¹: {len(relaxation_times)} valid reference points")
    
    # Highlight coupling turn-on time at 200 ps
    ax2.axvline(x=200, color='red', linestyle='--', linewidth=2, alpha=0.8,
                label=latex_format(r'Coupling Turn-On ($t = 200$ ps)', 'Coupling Turn-On (t = 200 ps)'))
    
    # Add shaded region to emphasize before/after coupling
    if ax2.get_xlim()[0] < 200:
        ax2.axvspan(ax2.get_xlim()[0], 200, alpha=0.1, color='gray', 
                   label=latex_format(r'Pre-Coupling Phase', 'Pre-Coupling Phase'))
    
    ax2.set_xlabel(latex_format(r'Reference Time $t$ (ps)', 'Reference Time t (ps)'))
    ax2.set_ylabel(latex_format(r'Normalized Relaxation Time $\tau/\tau_0$', 'Normalized Relaxation Time τ/τ₀'))
    ax2.set_title(latex_format(r'Normalized Relaxation Time vs Reference Time' + '\n' + r'(Direct Threshold: $F(k,t) = 0.1$)',
                              'Normalized Relaxation Time vs Reference Time\n(Direct Threshold: F(k,t) = 0.1)'))
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9, loc='best')
    
    # Set x-axis to start at 0 for Panel 2
    ax2.set_xlim(left=0)
    
    # Add colorbar for Panel 2 (frequencies)
    if panel2_scatter_points:
        cbar2 = plt.colorbar(panel2_scatter_points[0], ax=ax2, shrink=0.8)
        cbar2.set_label(latex_format(r'Frequency $\nu$ (cm$^{-1}$)', 'Frequency ν (cm⁻¹)'), rotation=270, labelpad=20)
        
        # Format colorbar ticks for frequencies
        coupling_ticks = coupling_values_sorted
        frequency_ticks = [coupling_to_frequency(cv) for cv in coupling_ticks]
        frequency_tick_labels = [f'{fv:.0f}' for fv in frequency_ticks]
        cbar2.set_ticks(coupling_ticks)
        cbar2.set_ticklabels(frequency_tick_labels)
    
    plt.suptitle(latex_format(r'Normalized Relaxation Time Analysis: Direct Threshold $F(k,t) = 0.1$' + '\n' + 
                             r'Cavity MD with Coupling Turn-On at $t = 200$ ps (Normalized by $\tau_0$ at $t_w = 200$ ps)',
                             'Normalized Relaxation Time Analysis: Direct Threshold F(k,t) = 0.1\nCavity MD with Coupling Turn-On at t = 200 ps (Normalized by τ₀ at t_w = 200 ps)'),
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    output_file = Path(output_dir) / 'relaxation_analysis_filtered.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()
    
    # Reset matplotlib settings
    plt.rcParams.update(plt.rcParamsDefault)

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Analyze and plot F(k,t) data from averaged replica files')
    parser.add_argument('--base_dir', default='./', 
                       help='Base directory containing coupling folders (default: ./)')
    parser.add_argument('--output_dir', default='./', 
                       help='Output directory for plots (default: ./)')
    
    args = parser.parse_args()
    
    # Set matplotlib style to classic globally
    plt.style.use('classic')
    
    print("F(k,t) Analysis and Plotting Script")
    print("=" * 50)
    print(f"Base directory: {Path(args.base_dir).resolve()}")
    print(f"Output directory: {Path(args.output_dir).resolve()}")
    
    # Collect all data
    print("\nCollecting F(k,t) data...")
    data_dict, coupling_info = collect_data(args.base_dir)
    
    if not data_dict:
        print("No data found! Make sure you have cavity_coupling_* folders with master_fskt_ref*.txt files.")
        return 1
    
    print(f"\nSummary:")
    print(f"  Found {len(data_dict)} coupling strengths")
    
    total_refs = set()
    for folder_data in data_dict.values():
        total_refs.update(folder_data.keys())
    
    print(f"  Total unique references: {sorted(total_refs)}")
    
    # Create output directory if it doesn't exist
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Generate all plots
    plot_fkt_by_coupling(data_dict, coupling_info, args.output_dir)
    plot_fkt_by_ref(data_dict, coupling_info, args.output_dir)
    plot_fkt_diagnostic(data_dict, coupling_info, args.output_dir)
    plot_relaxation_analysis(data_dict, coupling_info, args.output_dir)
    
    print("\n" + "=" * 50)
    print("Analysis complete! Generated plots:")
    print(f"  1. fkt_by_coupling_filtered.png - F(k,t) for different refs (by coupling)")
    print(f"  2. fkt_by_ref_filtered.png - F(k,t) for different couplings (by ref)")
    print(f"  3. fkt_diagnostic_filtered.png - F(k,t) diagnostic with 1/e lines")
    print(f"  4. relaxation_analysis_filtered.png - Relaxation time analysis with color gradations")
    
    return 0

if __name__ == "__main__":
    exit(main())
