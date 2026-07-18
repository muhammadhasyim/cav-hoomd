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
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import subprocess

def parse_coupling_strength(folder_name):
    """Parse coupling strength from folder name."""
    coupling_part = folder_name.replace('cavity_coupling_', '')
    
    if '0epos00' in coupling_part:
        return 0.0, '0.0'
    elif 'eneg' in coupling_part:
        parts = coupling_part.split('eneg')
        if len(parts) == 2:
            mantissa = int(parts[0])
            exponent = int(parts[1][1])
            value = mantissa * (10 ** (-exponent))
            return value, f'{mantissa}×10⁻{exponent}'
    
    return float('nan'), coupling_part

def filter_coupling_strengths(coupling_info):
    """Filter to only include specific coupling strengths: 0.0, 3e-4, 5e-4, 1e-3."""
    target_couplings = [0.0, 3e-4, 5e-4, 7e-4, 1e-3, 2e-3, 4e-3, 1e-2]
    
    filtered_couplings = {}
    for coupling_name, (coupling_value, coupling_label) in coupling_info.items():
        if not np.isnan(coupling_value) and any(abs(coupling_value - target) < 1e-6 for target in target_couplings):
            filtered_couplings[coupling_name] = (coupling_value, coupling_label)
    
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
        data = pd.read_csv(file_path, sep=r'\s+', comment='#')#, header=None)
        
        time = data.iloc[:, 0].values
        fkt = data.iloc[:, 1].values
        
        return time, fkt
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None, None

def find_relaxation_time(time, fkt, target_value=1/np.e, normalization_value=None):
    """Find the relaxation time where F(k,t) = target_value."""
    try:
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
    coupling_folders = list(base_path.glob('cavity_coupling_*200.0ps'))
    
    if not coupling_folders:
        print(f"No coupling folders found in {base_dir}")
        return {}, {}
    
    data_dict = {}
    coupling_info = {}
    
    for folder in coupling_folders:
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
    
    plt.suptitle('F(k,t) vs Time: Different References by Coupling Strength\n(Filtered: 0.0, 3×10⁻⁴, 5×10⁻⁴, 1×10⁻³)', 
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
    
    plt.suptitle('F(k,t) vs Time: Different Coupling Strengths by Reference\n(Filtered: 0.0, 3×10⁻⁴, 5×10⁻⁴, 1×10⁻³)', 
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
    
    target_line = 1/np.e  # 1/e line for reference
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
                
                # Add 1/e reference line (only if within range)
                if target_line >= 0.001:
                    ax.axhline(y=target_line, color='red', linestyle='--', alpha=0.7, 
                          label=f'1/e ≈ {target_line:.3f}')
                
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
    
    plt.suptitle('F(k,t) Diagnostic: Understanding Relaxation Behavior\n(Filtered: 0.0, 3×10⁻⁴, 5×10⁻⁴, 1×10⁻³)', 
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
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    # Determine all available references across all couplings
    all_refs = set()
    for coupling_name in sorted_couplings:
        if coupling_name in data_dict:
            all_refs.update(data_dict[coupling_name].keys())
    all_refs = sorted(all_refs)
    
    # Panel 1: Relaxation time vs coupling strength for different references
    print("  Calculating relaxation times vs coupling strength for different references...")
    
    coupling_values = [filtered_coupling_info[c][0] for c in sorted_couplings if c in data_dict]
    coupling_labels = [filtered_coupling_info[c][1] for c in sorted_couplings if c in data_dict]
    
    valid_couplings = [c for c in sorted_couplings if c in data_dict]
    
    # Create color mapping for reference times (Panel 1)
    ref_times_ps = [ref_num * 200 for ref_num in all_refs]
    ref_time_norm = Normalize(vmin=0, vmax=max(ref_times_ps))
    ref_time_cmap = plt.colormaps.get_cmap('viridis')
    
    # Store data for Panel 1 colorbar
    panel1_scatter_points = []
    
    for i, ref_num in enumerate(all_refs):
        ref_relaxation_times = []
        ref_coupling_values = []
        ref_coupling_labels = []
        
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
                    ref_relaxation_times.append(rel_time)
                    ref_coupling_values.append(coupling_value)
                    ref_coupling_labels.append(coupling_label)
        
        if ref_relaxation_times:  # Only plot if we have data
            ref_time_ps = ref_num * 200
            color = ref_time_cmap(ref_time_norm(ref_time_ps))
            
            scatter = ax1.scatter(ref_coupling_values, ref_relaxation_times, 
                               c=[ref_time_ps]*len(ref_relaxation_times), 
                               cmap='viridis', norm=ref_time_norm,
                               s=80, alpha=0.8, edgecolors='black', linewidth=0.5,
                               label=latex_format(f'$t = {ref_time_ps}$ ps', f't = {ref_time_ps} ps'))
            
            # Connect points with lines
            ax1.plot(ref_coupling_values, ref_relaxation_times, '-', 
                    color=color, linewidth=2, alpha=0.6)
            
            panel1_scatter_points.append(scatter)
            
            print(f"    ref{ref_num} (t={ref_time_ps} ps): {len(ref_relaxation_times)} valid coupling points")
    
    ax1.set_xlabel(latex_format(r'Coupling Strength', 'Coupling Strength'))
    ax1.set_ylabel(latex_format(r'Relaxation Time $\tau$ (ps)', 'Relaxation Time τ (ps)'))
    ax1.set_title(latex_format(r'Relaxation Time vs Coupling Strength' + '\n' + r'($F(k,t) = 1/e$ criterion)', 
                              'Relaxation Time vs Coupling Strength\n(F(k,t) = 1/e criterion)'))
    ax1.grid(True, alpha=0.3)
    
    if coupling_values:
        ax1.set_xticks(coupling_values)
        if use_latex:
            ax1.set_xticklabels([r'$' + label.replace('×', r'\times') + r'$' for label in coupling_labels], rotation=45)
        else:
            ax1.set_xticklabels(coupling_labels, rotation=45)
    
    # Set x-axis to start at 0 for Panel 1
    if coupling_values:
        ax1.set_xlim(left=0)
    
    # Add colorbar for Panel 1 (reference times)
    if panel1_scatter_points:
        cbar1 = plt.colorbar(panel1_scatter_points[0], ax=ax1, shrink=0.8)
        cbar1.set_label(latex_format(r'Reference Time $t$ (ps)', 'Reference Time t (ps)'), rotation=270, labelpad=20)
    
    # Panel 2: Relaxation time vs reference time for different coupling strengths
    print("  Calculating relaxation times vs reference time for different couplings...")
    
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
                ref_times.append(ref_time_ps)
                relaxation_times.append(rel_time)
        
        if relaxation_times:  # Only plot if we have data
            color = coupling_cmap(coupling_norm(coupling_value))
            
            scatter = ax2.scatter(ref_times, relaxation_times,
                               c=[coupling_value]*len(relaxation_times),
                               cmap='coolwarm', norm=coupling_norm,
                               s=80, alpha=0.8, edgecolors='black', linewidth=0.5,
                               label=latex_format(f'$g = {coupling_label.replace("×", r"\\times")}$', 
                                                 f'g = {coupling_label}'))
            
            # Connect points with lines
            ax2.plot(ref_times, relaxation_times, '-',
                    color=color, linewidth=2, alpha=0.6)
            
            panel2_scatter_points.append(scatter)
            
            print(f"    Coupling {coupling_label}: {len(relaxation_times)} valid reference points")
    
    # Highlight coupling turn-on time at 200 ps
    ax2.axvline(x=200, color='red', linestyle='--', linewidth=2, alpha=0.8,
                label=latex_format(r'Coupling Turn-On ($t = 200$ ps)', 'Coupling Turn-On (t = 200 ps)'))
    
    # Add shaded region to emphasize before/after coupling
    if ax2.get_xlim()[0] < 200:
        ax2.axvspan(ax2.get_xlim()[0], 200, alpha=0.1, color='gray', 
                   label=latex_format(r'Pre-Coupling Phase', 'Pre-Coupling Phase'))
    
    ax2.set_xlabel(latex_format(r'Reference Time $t$ (ps)', 'Reference Time t (ps)'))
    ax2.set_ylabel(latex_format(r'Relaxation Time $\tau$ (ps)', 'Relaxation Time τ (ps)'))
    ax2.set_title(latex_format(r'Relaxation Time vs Reference Time' + '\n' + r'($F(k,t) = 1/e$ criterion)',
                              'Relaxation Time vs Reference Time\n(F(k,t) = 1/e criterion)'))
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9, loc='best')
    
    # Set x-axis to start at 0 for Panel 2
    ax2.set_xlim(left=0)
    
    # Add colorbar for Panel 2 (coupling strengths)
    if panel2_scatter_points:
        cbar2 = plt.colorbar(panel2_scatter_points[0], ax=ax2, shrink=0.8)
        cbar2.set_label(latex_format(r'Coupling Strength $g$', 'Coupling Strength g'), rotation=270, labelpad=20)
        
        # Format colorbar ticks for coupling strengths
        coupling_ticks = coupling_values_sorted
        if use_latex:
            coupling_tick_labels = [filtered_coupling_info[c][1].replace('×', r'$\times$') 
                                   for c in valid_couplings]
            cbar2.set_ticks(coupling_ticks)
            cbar2.set_ticklabels([f'${label}$' for label in coupling_tick_labels])
        else:
            coupling_tick_labels = [filtered_coupling_info[c][1] for c in valid_couplings]
            cbar2.set_ticks(coupling_ticks)
            cbar2.set_ticklabels(coupling_tick_labels)
    
    plt.suptitle(latex_format(r'Relaxation Time Analysis: $F(k,t) = 1/e$ Criterion' + '\n' + 
                             r'Cavity MD with Coupling Turn-On at $t = 200$ ps',
                             'Relaxation Time Analysis: F(k,t) = 1/e Criterion\nCavity MD with Coupling Turn-On at t = 200 ps'),
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
