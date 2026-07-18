#!/usr/bin/env python3
"""
Plot fictive relaxation time vs waiting time relaxation analysis.

This script compares:
1. Fictive relaxation time from time_series_output (normalized by first 200 ps, shifted by -400 ps)
2. Relaxation time analysis from F(k,t) data over waiting times

Author: Assistant
Date: September 10, 2025
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import subprocess

# Note: plot_fkt_analysis functions not used in current implementation

def load_fkt_data(base_dir):
    """
    Load F(k,t) data from cavity coupling directories.
    
    This is a simple replacement for the missing collect_data function.
    Looks for F(k,t) data files in cavity_coupling_* directories.
    
    Parameters:
    -----------
    base_dir : str
        Base directory containing cavity_coupling_* folders
        
    Returns:
    --------
    data_dict : dict
        Dictionary with F(k,t) data for each coupling
    coupling_info : dict
        Dictionary with coupling strength information
    """
    import glob
    import os
    
    data_dict = {}
    coupling_info = {}
    
    # Look for cavity coupling directories
    pattern = os.path.join(base_dir, "cavity_coupling_*_switch_*")
    coupling_dirs = glob.glob(pattern)
    
    for coupling_dir in coupling_dirs:
        try:
            # Extract coupling information from directory name
            dir_name = os.path.basename(coupling_dir)
            parts = dir_name.split('_')
            
            # Find coupling strength part (e.g., "1eneg03", "0epos00")
            coupling_part = None
            for part in parts:
                if 'eneg' in part or 'epos' in part:
                    coupling_part = part
                    break
            
            if coupling_part is None:
                print(f"Warning: Could not parse coupling from {dir_name}")
                continue
                
            # Convert coupling string to float
            if 'eneg' in coupling_part:
                # e.g., "1eneg03" -> 1e-03
                mantissa = coupling_part.split('eneg')[0]
                exponent = coupling_part.split('eneg')[1]
                coupling_value = float(mantissa) * (10 ** (-int(exponent)))
            elif 'epos' in coupling_part:
                # e.g., "0epos00" -> 0e+00 = 0
                mantissa = coupling_part.split('epos')[0]
                exponent = coupling_part.split('epos')[1]
                coupling_value = float(mantissa) * (10 ** int(exponent))
            else:
                continue
                
            # Look for ref1 file (F(k,t) with waiting time tw = 200 ps)
            ref1_file = os.path.join(coupling_dir, "master_fskt_ref1.txt")
            
            if os.path.exists(ref1_file):
                print(f"Loading F(k,t) data (tw=200ps) from {os.path.relpath(ref1_file, base_dir)}")
                
                # Load the F(k,t) data
                try:
                    # Read and find data start (skip comments and column headers)
                    with open(ref1_file, 'r') as f:
                        lines = f.readlines()
                    
                    # Find the first data line
                    data_start = 0
                    for i, line in enumerate(lines):
                        line = line.strip()
                        if line and not line.startswith('#') and not line.startswith('lag_time'):
                            # Check if this line contains numbers
                            try:
                                float(line.split()[0])
                                data_start = i
                                break
                            except (ValueError, IndexError):
                                continue
                    
                    data = np.loadtxt(ref1_file, skiprows=data_start)
                    lag_times = data[:, 0]  # First column: lag time
                    fkt_values = data[:, 1]  # Second column: F(k,t)
                    
                    # Store data
                    data_dict[dir_name] = {
                        'lag_times': lag_times,
                        'fkt_values': fkt_values,
                        'coupling_value': coupling_value
                    }
                    
                    coupling_info[dir_name] = {
                        'coupling_value': coupling_value,
                        'directory': coupling_dir
                    }
                    
                except Exception as e:
                    print(f"Warning: Could not load F(k,t) data from {ref1_file}: {e}")
                    continue
            else:
                print(f"Warning: master_fskt_ref1.txt not found in {coupling_dir}")
                
        except Exception as e:
            print(f"Warning: Error processing {coupling_dir}: {e}")
            continue
    
    print(f"Successfully loaded F(k,t) data from {len(data_dict)} coupling directories")
    return data_dict, coupling_info

def load_fictive_data(time_series_dir):
    """
    Load fictive relaxation time data from time_series_output directory.
    
    Parameters:
    -----------
    time_series_dir : str or Path
        Directory containing fictive time series files
        
    Returns:
    --------
    fictive_data : dict
        Dictionary with coupling names as keys and (time, relaxation_time) as values
    """
    fictive_data = {}
    time_series_path = Path(time_series_dir)
    
    # Map file patterns to coupling strengths
    coupling_map = {
        'coupling_0epos00_fictive_time_series.txt': ('cavity_coupling_0epos00_switch_200.0ps', 0.0),
        'coupling_3eneg04_fictive_time_series.txt': ('cavity_coupling_3eneg04_switch_200.0ps', 3e-4),
        'coupling_5eneg04_fictive_time_series.txt': ('cavity_coupling_5eneg04_switch_200.0ps', 5e-4),
        'coupling_7eneg04_fictive_time_series.txt': ('cavity_coupling_7eneg04_switch_200.0ps', 7e-4),
        'coupling_1eneg03_fictive_time_series.txt': ('cavity_coupling_1eneg03_switch_200.0ps', 1e-3),
    }
    
    for filename, (coupling_name, coupling_value) in coupling_map.items():
        file_path = time_series_path / filename
        
        if file_path.exists():
            print(f"Loading {filename}...")
            
            # Load data (skip header lines starting with # and column names)
            data = np.loadtxt(file_path, skiprows=12)  # Skip the header and column names
            
            time_ps = data[:, 0]  # Column 1: Time (ps)
            fictive_temp = data[:, 1]  # Column 2: Fictive temperature (K)
            fictive_relax_time = data[:, 2]  # Column 3: Fictive relaxation time (ps)
            
            # Apply user's transformations:
            # 1. For display purposes, calculate normalization by value in first 200 ps
            first_200ps_mask = time_ps <= 200.0
            if np.any(first_200ps_mask):
                normalization_value = np.mean(fictive_relax_time[first_200ps_mask])
                fictive_relax_time_normalized = fictive_relax_time / normalization_value
            else:
                print(f"Warning: No data in first 200 ps for {filename}")
                normalization_value = fictive_relax_time[0]
                fictive_relax_time_normalized = fictive_relax_time / fictive_relax_time[0]
            
            # 2. Shift time origin by -400 ps
            time_shifted = time_ps - 400.0
            
            # 3. Only keep data after t = 0 (after shifting)
            valid_mask = time_shifted >= 0.0
            time_final = time_shifted[valid_mask]
            fictive_relax_final_normalized = fictive_relax_time_normalized[valid_mask]  # For display
            fictive_relax_final_actual = fictive_relax_time[valid_mask]  # Actual values in ps for material time
            
            fictive_data[coupling_name] = {
                'time': time_final,
                'relaxation_time': fictive_relax_final_normalized,  # For display (normalized)
                'relaxation_time_actual_ps': fictive_relax_final_actual,  # For material time calculation (actual ps)
                'coupling_value': coupling_value,
                'normalization_value': normalization_value
            }
            
            print(f"  Loaded {len(time_final)} points after filtering (t >= 0)")
            print(f"  Time range: {time_final[0]:.1f} to {time_final[-1]:.1f} ps")
            print(f"  Normalization value: {normalization_value:.3f} ps")
    
    return fictive_data

def calculate_waiting_time_relaxation(data_dict, coupling_info):
    """
    Calculate relaxation times as a function of waiting time from F(k,t) data.
    This replicates the analysis from plot_fkt_analysis.py.
    
    Parameters:
    -----------
    data_dict : dict
        F(k,t) data from collect_data()
    coupling_info : dict
        Coupling information from collect_data()
        
    Returns:
    --------
    waiting_time_data : dict
        Dictionary with coupling names as keys and waiting time analysis
    """
    waiting_time_data = {}
    
    # Filter to specific coupling strengths (same as in normalized_relaxation_analysis_filtered.png)
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    # Get all available reference frames
    all_refs = set()
    for folder_data in data_dict.values():
        all_refs.update(folder_data.keys())
    all_refs = sorted(all_refs)
    
    print(f"\nCalculating waiting time relaxation analysis...")
    print(f"Available references: {all_refs}")
    
    # Calculate relaxation times for each coupling
    for coupling_name in filtered_coupling_info:
        if coupling_name not in data_dict:
            continue
            
        coupling_value, coupling_label = filtered_coupling_info[coupling_name]
        print(f"\nProcessing {coupling_name} (g = {coupling_value})")
        
        # Get normalization value from ref1 (tw = 0)
        normalization_value = None
        if 1 in data_dict[coupling_name]:
            time_ref1, fkt_ref1 = data_dict[coupling_name][1]
            nonzero_mask = fkt_ref1 != 0
            if np.any(nonzero_mask):
                first_nonzero_idx = np.where(nonzero_mask)[0][0]
                normalization_value = fkt_ref1[first_nonzero_idx]
        
        waiting_times = []
        relaxation_times = []
        
        for ref_num in all_refs:
            if ref_num == 0:  # Skip ref0 (pre-equilibrium)
                continue
                
            if ref_num in data_dict[coupling_name]:
                time, fkt = data_dict[coupling_name][ref_num]
                
                # Calculate waiting time: tw = (ref_num - 1) * 200 ps
                waiting_time = (ref_num - 1) * 200.0
                
                # Find relaxation time (where F(k,t) = 1/e of normalized value)
                rel_time = find_relaxation_time(time, fkt, target_value=1/np.exp(1), 
                                                normalization_value=normalization_value)
                    
                if not np.isnan(rel_time):
                    waiting_times.append(waiting_time)
                    relaxation_times.append(rel_time)
                    print(f"  ref{ref_num}: tw = {waiting_time} ps, τ = {rel_time:.2f} ps")
        
        if waiting_times:
            waiting_time_data[coupling_name] = {
                'waiting_times': np.array(waiting_times),
                'relaxation_times': np.array(relaxation_times),
                'coupling_value': coupling_value,
                'coupling_label': coupling_label
            }
    
    return waiting_time_data

def load_material_time_reference_data(file_path):
    """
    Load material time reference data from the provided text file.
    
    The file contains material time evolution data for different coupling strengths
    in column pairs: (time_ps, xi) for each coupling.
    
    Parameters:
    -----------
    file_path : str
        Path to the material time data file
        
    Returns:
    --------
    material_time_ref_data : dict
        Dictionary with material time reference data for each coupling
    """
    material_time_ref_data = {}
    
    try:
        # Load the data, skipping comment lines
        data = np.loadtxt(file_path, comments='#')
        
        # Column structure from header:
        # t_g0.0_ps xi_g0.0 t_g3e-4_ps xi_g3e-4 t_g5e-4_ps xi_g5e-4 t_g7e-4_ps xi_g7e-4 t_g1e-3_ps xi_g1e-3
        
        coupling_mapping = {
            0: (0.0, 'g0.0'),      # columns 0,1
            1: (3e-4, 'g3e-4'),    # columns 2,3
            2: (5e-4, 'g5e-4'),    # columns 4,5
            3: (7e-4, 'g7e-4'),    # columns 6,7
            4: (1e-3, 'g1e-3')     # columns 8,9
        }
        
        for idx, (coupling_value, coupling_label) in coupling_mapping.items():
            time_col = idx * 2
            xi_col = idx * 2 + 1
            
            if xi_col < data.shape[1]:
                time_data = data[:, time_col]
                xi_data = data[:, xi_col]
                
                # Store in a format similar to our material time data
                material_time_ref_data[coupling_label] = {
                    'time': time_data,
                    'material_time': xi_data,
                    'coupling_value': coupling_value,
                    'coupling_label': coupling_label
                }
        
        print(f"Loaded reference material time data for {len(material_time_ref_data)} coupling strengths")
        return material_time_ref_data
        
    except Exception as e:
        print(f"Error loading material time reference data: {e}")
        return {}

def calculate_waiting_time_relaxation(data_dict, coupling_info):
    """
    Calculate waiting time relaxation analysis from F(k,t) data.
    
    For ref1 data, the waiting time is tw = 200 ps.
    We need to normalize F(k,t) and find the decay to exp(-1) ≈ 0.368.
    
    Parameters:
    -----------
    data_dict : dict
        Dictionary with F(k,t) data for each coupling
    coupling_info : dict
        Dictionary with coupling strength information
        
    Returns:
    --------
    waiting_time_data : dict
        Dictionary with waiting time relaxation analysis results
    """
    waiting_time_data = {}
    
    for dir_name, data in data_dict.items():
        coupling_value = data['coupling_value']
        lag_times = data['lag_times']
        fkt_values = data['fkt_values']
        
        # Normalize F(k,t) by its initial value (t=0)
        if len(fkt_values) > 0 and fkt_values[0] > 0:
            normalized_fkt = fkt_values / fkt_values[0]
        else:
            print(f"Warning: Invalid F(k,t) data for {dir_name}")
            continue
        
        # Find relaxation time (when F(k,t) decays to exp(-1) ≈ 0.368)
        target_value = np.exp(-1)  # ≈ 0.368
        
        # Find where normalized F(k,t) crosses the target value
        below_target = normalized_fkt <= target_value
        if np.any(below_target):
            first_below_idx = np.where(below_target)[0][0]
            if first_below_idx > 0:
                # Linear interpolation to get precise crossing time
                idx_before = first_below_idx - 1
                idx_after = first_below_idx
                
                t1, t2 = lag_times[idx_before], lag_times[idx_after]
                y1, y2 = normalized_fkt[idx_before], normalized_fkt[idx_after]
                
                if y2 != y1:
                    relaxation_time = t1 + (target_value - y1) * (t2 - t1) / (y2 - y1)
                else:
                    relaxation_time = t1
            else:
                relaxation_time = lag_times[0]
        else:
            # If never reaches target, use NaN
            relaxation_time = np.nan
        
        # Store results in waiting time data format
        # For ref1 data, we have one waiting time (tw = 200 ps) and one relaxation time
        waiting_time_data[dir_name] = {
            'waiting_times': np.array([200.0]),  # tw = 200 ps for ref1
            'relaxation_times': np.array([relaxation_time]),
            'coupling_value': coupling_value,
            'lag_times': lag_times,
            'normalized_fkt': normalized_fkt
        }
    
    return waiting_time_data

def calculate_material_time(time, fictive_relaxation_time_actual_ps, scaling_factor=1.0):
    """
    Calculate material time as the integral of 1/fictive_relaxation_time starting from t=0 ps.
    This ensures consistency with the F(k,t) material time analysis.
    
    The material time is defined as:
    Ξ(t) = ∫₀ᵗ (1/τ_fictive(t')) dt' for t >= 0 ps
    ξ(t) = Ξ(t) × scaling_factor
    
    Then we calculate the material time difference ξ(t) - ξ(t₀) where t₀ = 0 ps.
    
    Parameters:
    -----------
    time : array
        Time points (ps) - original time starting from 0
    fictive_relaxation_time_actual_ps : array
        Fictive relaxation time values in picoseconds (actual values, not normalized)
    scaling_factor : float
        Multiplicative factor to rescale material time
        
    Returns:
    --------
    material_time_true : array
        True material time Ξ(t)
    material_time_scaled : array
        Scaled material time ξ(t) = Ξ(t) × scaling_factor
    material_time_diff : array
        Material time difference ξ(t) - ξ(0) starting from t=0
    exp_neg_material_time_diff : array
        exp(-(ξ(t) - ξ(0))^0.55) starting from t=0
    """
    # Calculate 1/τ_fictive, avoiding division by zero
    inv_tau_fictive = np.where(fictive_relaxation_time_actual_ps > 0, 
                              1.0 / fictive_relaxation_time_actual_ps, 0.0)
    
    # Start integration from t = 0 ps (beginning of simulation)
    t_start = 0.0  # ps
    start_idx = np.argmin(np.abs(time - t_start))
    
    # Integrate using trapezoidal rule starting from t = 0 ps
    # Ξ(t) = ∫₀ᵗ (1/τ_fictive(t')) dt' for t >= 0 ps
    integral_inv_tau_fictive = np.zeros_like(time)
    for i in range(start_idx + 1, len(time)):
        dt = time[i] - time[i-1]
        # Trapezoidal integration of 1/τ_fictive starting from t=0 ps
        integral_inv_tau_fictive[i] = integral_inv_tau_fictive[i-1] + 0.5 * (inv_tau_fictive[i-1] + inv_tau_fictive[i]) * dt
    
    # True material time: Ξ(t) = ∫₀ᵗ (1/τ_fictive(t')) dt'
    material_time_true = integral_inv_tau_fictive
    
    # Scaled material time: ξ(t) = Ξ(t) × scaling_factor
    material_time_scaled = material_time_true * scaling_factor
    
    # Calculate material time difference: ξ(t) - ξ(t₀) with t₀ = 0 ps
    # t₀ = 0 ps corresponds to start_idx in our array
    t0_idx = start_idx
    xi_t0 = material_time_scaled[t0_idx]  # ξ(t₀) where t₀ = 0 ps
    
    # Calculate ξ(t) - ξ(t₀) for all time points
    # This gives material time accumulated since t₀ = 0 ps
    material_time_diff = material_time_scaled - xi_t0
    
    # Calculate exp(-(ξ(t) - ξ(t₀))^0.55) = exp(-material_time_diff^0.55)
    material_time_diff_power = np.where(material_time_diff >= 0, 
                                       np.power(material_time_diff, 0.55), 0.0)
    exp_neg_material_time_diff = np.exp(-material_time_diff_power)
    
    return material_time_true, material_time_scaled, material_time_diff, exp_neg_material_time_diff

def find_decay_time(time, exp_neg_material_time_diff, target_value=0.1):
    """
    Find the time when exp(-(ξ(t+tw) - ξ(tw))^0.55) decays to a target value.
    
    Parameters:
    -----------
    time : array
        Time points (ps)
    exp_neg_material_time_diff : array
        exp(-(ξ(t+tw) - ξ(tw))^0.55) values
    target_value : float
        Target value to find (default 0.1)
        
    Returns:
    --------
    decay_time : float
        Time when the function reaches target_value (ps), or NaN if not found
    """
    # Find indices where the function crosses the target value
    # Look for the first time it goes below the target
    below_target = exp_neg_material_time_diff <= target_value
    
    if not np.any(below_target):
        # Never reaches target value
        return np.nan
    
    # Find first index where it goes below target
    first_below_idx = np.where(below_target)[0][0]
    
    if first_below_idx == 0:
        # Already below target at t=0
        return time[0]
    
    # Interpolate between the two points to get more precise crossing time
    idx_before = first_below_idx - 1
    idx_after = first_below_idx
    
    t1, t2 = time[idx_before], time[idx_after]
    y1, y2 = exp_neg_material_time_diff[idx_before], exp_neg_material_time_diff[idx_after]
    
    # Linear interpolation to find exact crossing point
    if y2 != y1:  # Avoid division by zero
        decay_time = t1 + (target_value - y1) * (t2 - t1) / (y2 - y1)
    else:
        decay_time = t1
    
    return decay_time

def plot_fictive_vs_waiting_comparison(fictive_data, waiting_time_data, output_dir='.'):
    """
    Create comparison plot of fictive relaxation time vs waiting time relaxation analysis,
    and material time analysis.
    
    Parameters:
    -----------
    fictive_data : dict
        Fictive relaxation time data
    waiting_time_data : dict
        Waiting time relaxation analysis data
    output_dir : str
        Output directory for plots
    """
    print("\nCreating fictive vs waiting time relaxation comparison plot...")
    
    # Set up LaTeX rendering (same as in plot_fkt_analysis.py)
    plt.style.use('classic')
    use_latex = True
    try:
        subprocess.check_output(['latex', '--version'], stderr=subprocess.DEVNULL)
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Computer Modern Roman']
        plt.rcParams['font.sans-serif'] = ['Computer Modern Sans Serif']
        plt.rcParams['font.monospace'] = ['Computer Modern Typewriter']
        plt.rcParams['mathtext.fontset'] = 'cm'
        plt.rcParams['mathtext.rm'] = 'serif'
        plt.rcParams['axes.unicode_minus'] = False
        plt.rcParams['text.latex.preamble'] = r'\usepackage{lmodern}'
        print("  Using LaTeX with Computer Modern fonts")
    except (subprocess.CalledProcessError, FileNotFoundError, ImportError):
        print("  Warning: LaTeX not available, using default matplotlib rendering")
        plt.rcParams['text.usetex'] = False
        plt.rcParams['font.family'] = 'DejaVu Sans'
        plt.rcParams['mathtext.fontset'] = 'dejavusans'
        use_latex = False
    
    # Helper function for LaTeX formatting
    def latex_format(text, fallback_text=None):
        if use_latex:
            return text
        else:
            return fallback_text if fallback_text else text.replace('$', '').replace(r'\tau', 'τ').replace(r'\times', '×').replace(r'\mathrm{', '').replace('}', '')
    
    # Create figure with two panels: relaxation comparison and material time analysis
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Use the same color mapping as F(k,t) analysis (coolwarm colormap)
    # Get all coupling values to set up normalization
    all_coupling_values = []
    for data in fictive_data.values():
        all_coupling_values.append(data['coupling_value'])
    for data in waiting_time_data.values():
        all_coupling_values.append(data['coupling_value'])
    
    # Remove duplicates and sort
    unique_coupling_values = sorted(set(all_coupling_values))
    
    # Create the same color mapping as in F(k,t) analysis
    from matplotlib.colors import Normalize
    coupling_norm = Normalize(vmin=0, vmax=max(unique_coupling_values))
    coupling_cmap = plt.colormaps.get_cmap('coolwarm')
    
    def get_coupling_color(coupling_value):
        """Get color for a given coupling value using the same scheme as F(k,t) analysis."""
        return coupling_cmap(coupling_norm(coupling_value))
    
    # First, normalize F(k,t) relaxation time data by g=0 values
    # Check if we have F(k,t) data
    if not waiting_time_data:
        print("  No F(k,t) data available, showing only fictive and material time analysis...")
    else:
        print("  Normalizing F(k,t) data by g=0 values...")
        
        # Find g=0 data for normalization
        g0_data = None
        for coupling_name, data in waiting_time_data.items():
            if data['coupling_value'] == 0.0:
                g0_data = data
                break
        
        if g0_data is None:
            print("  Warning: No g=0 data found for normalization!")
            print("  Proceeding with fictive and material time analysis only...")
            waiting_time_data = {}  # Clear to skip F(k,t) plotting
    
    # Only extract g=0 data if we have waiting time data
    # For material time reference data, we'll skip the first panel F(k,t) comparison
    # and focus on the material time comparison in the second panel
    if waiting_time_data and 'waiting_times' in g0_data:
        g0_waiting_times = g0_data['waiting_times']
        g0_relaxation_times = g0_data['relaxation_times']
    else:
        # Material time reference data - skip F(k,t) plotting in first panel
        waiting_time_data_for_plotting = {}
    
    # Plot both datasets in the first panel for direct comparison
    print("  Plotting both fictive and F(k,t) relaxation data together...")
    
    # Calculate material time data for all fictive datasets
    material_time_data = {}
    
    # Material time calculation using actual fictive relaxation time (ps)
    print("    Using actual fictive relaxation time (column 3) for material time calculation...")
    
    # Step 1: Calculate scaling factor based on REFERENCE ε=0 material time to have slope 1
    print("    Calculating scaling factor to make reference ε=0 material time have slope 1...")
    
    # First, check if we have reference material time data
    if not waiting_time_data:
        print("    Error: No reference material time data available for scaling!")
        return
    
    # Find reference ε=0 data
    ref_epsilon_0_data = None
    for coupling_name, data in waiting_time_data.items():
        if data['coupling_value'] == 0.0:
            ref_epsilon_0_data = data
            break
    
    if ref_epsilon_0_data is None:
        print("    Error: Could not find reference ε=0 data for scaling!")
        return
    
    # Calculate the slope of the reference ε=0 material time
    ref_time = ref_epsilon_0_data['time']
    ref_material_time = ref_epsilon_0_data['material_time']
    
    # Calculate slope dξ/dt for reference data (using linear fit)
    # Use middle portion of data to avoid edge effects
    mid_start = len(ref_time) // 4
    mid_end = 3 * len(ref_time) // 4
    if mid_end > mid_start + 10:  # Need enough points
        fit_time = ref_time[mid_start:mid_end] - ref_time[mid_start]  # Start from 0
        fit_xi = ref_material_time[mid_start:mid_end] - ref_material_time[mid_start]  # Start from 0
        
        # Linear fit to get slope
        if len(fit_time) > 1 and np.max(fit_time) > 0:
            ref_slope = np.polyfit(fit_time, fit_xi, 1)[0]  # Slope coefficient
        else:
            ref_slope = (ref_material_time[-1] - ref_material_time[0]) / (ref_time[-1] - ref_time[0])
    else:
        ref_slope = (ref_material_time[-1] - ref_material_time[0]) / (ref_time[-1] - ref_time[0])
    
    print(f"    Reference ε=0 material time slope: {ref_slope:.6f}")
    print(f"    Target slope: 1.0")
    
    # Calculate fictive material time for ε=0 with unit scaling first
    epsilon_0_fictive_data = None
    for coupling_name, data in fictive_data.items():
        if data['coupling_value'] == 0.0:
            epsilon_0_fictive_data = data
            break
    
    if epsilon_0_fictive_data is None:
        print("    Error: Could not find ε=0 fictive data for scaling!")
        return
    
    # Calculate unscaled fictive material time for ε=0
    time_epsilon_0 = epsilon_0_fictive_data['time']
    relaxation_time_epsilon_0_actual = epsilon_0_fictive_data['relaxation_time_actual_ps']
    
    # Calculate with unit scaling first
    material_time_true_eps0, material_time_unit_eps0, _, _ = calculate_material_time(time_epsilon_0, relaxation_time_epsilon_0_actual, scaling_factor=1.0)
    
    # Calculate the slope of our fictive ε=0 material time
    if len(material_time_unit_eps0) > 1:
        fictive_slope = (material_time_unit_eps0[-1] - material_time_unit_eps0[0]) / (time_epsilon_0[-1] - time_epsilon_0[0])
    else:
        fictive_slope = 1.0  # Fallback
    
    print(f"    Fictive ε=0 material time slope (unscaled): {fictive_slope:.6f}")
    
    # Calculate scaling factor to make reference slope = 1.0
    # We want: ref_slope_normalized = 1.0
    # Since slope = dξ/dt, to make slope = 1, we need to multiply time by ref_slope
    time_normalization_factor = ref_slope if ref_slope > 0 else 1.0
    
    # For the fictive material time, we want it to have the same behavior as reference
    # We'll scale it so that when time is normalized, our fictive ε=0 also has slope ≈ 1
    scaling_factor = 1.0  # Keep natural material time units
    
    print(f"    Time normalization factor: {time_normalization_factor:.6f}")
    print(f"    This will make reference ε=0 slope = {ref_slope / time_normalization_factor:.6f} ≈ 1.0")
    print(f"    Fictive material time scaling factor: {scaling_factor:.6f}")
    
    # Store reference time scale for plotting
    reference_time_scale = 1.0 / time_normalization_factor  # This is the natural time scale of the reference
    
    # Step 2: Calculate scaled material time for all datasets using the proper scaling factor
    print("    Calculating material time for all coupling strengths...")
    for coupling_name, data in fictive_data.items():
        coupling_value = data['coupling_value']
        time = data['time']
        relaxation_time_actual = data['relaxation_time_actual_ps']
        
        # Calculate material time with proper scaling
        material_time_true, material_time_scaled, material_time_diff, exp_neg_material_time_diff = calculate_material_time(time, relaxation_time_actual, scaling_factor)
        material_time_data[coupling_name] = {
            'time': time,
            'material_time_true': material_time_true,
            'material_time_scaled': material_time_scaled,
            'material_time_diff': material_time_diff,
            'exp_neg_material_time_diff': exp_neg_material_time_diff,
            'coupling_value': coupling_value
        }
        
        # Report material time characteristics for this coupling
        if len(material_time_scaled) > 1:
            time_span = time[-1] - time[0]
            xi_span = material_time_scaled[-1] - material_time_scaled[0]
            avg_rate = xi_span / time_span if time_span > 0 else 0
            avg_tau = np.mean(relaxation_time_actual)
            
            coupling_label = f"ε = {coupling_value:.0e}" if coupling_value > 0 else "ε = 0"
            print(f"      {coupling_label}: <τ> = {avg_tau:.1f} ps, dξ/dt = {avg_rate:.6f} ps⁻¹")
    
    # First plot fictive relaxation time data in panel 1
    print("    Adding fictive relaxation time data...")
    for coupling_name, data in fictive_data.items():
        coupling_value = data['coupling_value']
        time = data['time']
        relaxation_time_normalized = data['relaxation_time'] * scaling_factor  # Apply the same scaling factor
        
        color = get_coupling_color(coupling_value)
        
        if coupling_value == 0.0:
            label = r'$\varepsilon = 0$ (fictive)'
        else:
            # Format coupling value in scientific notation to match F(k,t) labels
            if coupling_value >= 1e-3:
                label = f'$\\varepsilon = {coupling_value:.0e}$ (fictive)'
            else:
                # Convert to the same format as F(k,t) analysis (×10^-4 notation)
                scaled_val = coupling_value * 10000
                label = f'$\\varepsilon = {scaled_val:.0f}\\times10^{{-4}}$ (fictive)'
        
        ax1.plot(time-time[0], relaxation_time_normalized, color=color, linewidth=2, alpha=0.7, 
               linestyle='-', label=label)
    
    # Then plot normalized F(k,t) relaxation time (discrete points) if available
    # Only plot F(k,t) data if it has the waiting_times structure (not material time reference data)
    if waiting_time_data and any('waiting_times' in data for data in waiting_time_data.values()):
        print("    Adding normalized F(k,t) relaxation data...")
        for coupling_name, data in waiting_time_data.items():
            if 'waiting_times' not in data:  # Skip material time reference data
                continue
            coupling_value = data['coupling_value']
            waiting_times = data['waiting_times']
            relaxation_times = data['relaxation_times']
            
            color = get_coupling_color(coupling_value)
            
            if coupling_value == 0.0:
                label = r'$\varepsilon = 0$ (F(k,t))'
            else:
                # Format coupling value in scientific notation to match F(k,t) labels
                if coupling_value >= 1e-3:
                    label = f'$\\varepsilon = {coupling_value:.0e}$ (F(k,t))'
                else:
                    # Convert to the same format as F(k,t) analysis (×10^-4 notation)
                    scaled_val = coupling_value * 10000
                    label = f'$\\varepsilon = {scaled_val:.0f}\\times10^{{-4}}$ (F(k,t))'
            
            # Normalize by g=0 values at corresponding waiting times
            normalized_relaxation_times = []
            normalized_waiting_times = []
            
            for i, tw in enumerate(waiting_times):
                # Find corresponding g=0 value at same waiting time
                g0_idx = np.where(np.abs(g0_waiting_times - tw) < 1e-6)[0]
                if len(g0_idx) > 0:
                    g0_tau = g0_relaxation_times[g0_idx[0]]
                    normalized_tau = relaxation_times[i] / g0_tau
                    normalized_relaxation_times.append(normalized_tau)
                    normalized_waiting_times.append(tw)
                else:
                    print(f"      Warning: No g=0 reference found for tw={tw} ps")
            
            if normalized_relaxation_times:
                # Plot as points with different markers for distinction
                ax1.plot(normalized_waiting_times, normalized_relaxation_times, 'o', 
                       color=color, markersize=8, alpha=0.9, markerfacecolor='white',
                       markeredgewidth=2, markeredgecolor=color, label=label)
    else:
        print("    No F(k,t) data available for first panel. Showing only fictive relaxation time data.")
    
    # Set labels and formatting for first panel
    ax1.set_xlabel(latex_format(r'Time (ps)'))
    ax1.set_ylabel(latex_format(r'$\tau / \tau_0$', 'τ / τ_0'))
    ax1.set_title(latex_format(r'Relaxation Time Comparison: Fictive vs F(k,t) Analysis', 
                              'Relaxation Time Comparison: Fictive vs F(k,t) Analysis'))
    
    # Create a more organized legend for first panel
    handles, labels = ax1.get_legend_handles_labels()
    
    # Separate fictive and F(k,t) entries
    fictive_handles = []
    fictive_labels = []
    fkt_handles = []
    fkt_labels = []
    
    for handle, label in zip(handles, labels):
        if 'fictive' in label:
            fictive_handles.append(handle)
            fictive_labels.append(label.replace(' (fictive)', ''))
        else:
            fkt_handles.append(handle)
            fkt_labels.append(label.replace(' (F(k,t))', ''))
    
    # Create legend with two columns
    legend1 = ax1.legend(fictive_handles, fictive_labels, loc='upper right', 
                        title=latex_format(r'Fictive $\tau$ (lines)', 'Fictive τ (lines)'),
                        fontsize=8, title_fontsize=8)
    ax1.add_artist(legend1)
    
    legend2 = ax1.legend(fkt_handles, fkt_labels, loc='center right', 
                        title=latex_format(r'F(k,t) $\tau$ (points)', 'F(k,t) τ (points)'),
                        fontsize=8, title_fontsize=8)
    
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, None)
    ax1.set_ylim(0.5, None)  # Start y-axis from 0.5 to better show the data range
    
    # Second panel: Material time analysis
    print("    Adding material time analysis...")
    
    # Report material time scaling analysis
    print(f"\n    MATERIAL TIME SCALING ANALYSIS:")
    print(f"    Material time ξ(t) is properly scaled using reference time scale τ_ref = {reference_time_scale:.2f} ps")
    print(f"    Scaling factor = {scaling_factor:.6f} ps⁻¹ ensures ξ has correct time units")
    print(f"    For constant τ: ξ(t) = (t - t₀)/τ_ref, giving dξ/dt = 1/τ_ref")
    print(f"    Material time differences: ξ(t) - ξ(t₀) starting from t₀ = 0 ps")
    
    # Calculate and report characteristic material time spans for each coupling
    print(f"\n    Material time characteristics for each coupling:")
    print(f"    {'Coupling':<25} {'<τ> (ps)':<12} {'dξ/dt (ps⁻¹)':<15} {'Max ξ':<10}")
    print(f"    {'-'*25} {'-'*12} {'-'*15} {'-'*10}")
    
    for coupling_name, data in material_time_data.items():
        coupling_value = data['coupling_value']
        time = data['time']
        material_time_scaled = data['material_time_scaled']
        relaxation_time_actual = fictive_data[coupling_name]['relaxation_time_actual_ps']
        
        # Calculate characteristics
        avg_tau = np.mean(relaxation_time_actual)
        if len(material_time_scaled) > 1:
            time_span = time[-1] - time[0]
            xi_span = material_time_scaled[-1] - material_time_scaled[0]
            avg_rate = xi_span / time_span if time_span > 0 else 0
            max_xi = np.max(material_time_scaled) - np.min(material_time_scaled)
        else:
            avg_rate = 0
            max_xi = 0
        
        # Format coupling label for console output
        if coupling_value == 0.0:
            coupling_label = 'ε = 0'
        elif coupling_value >= 1e-3:
            coupling_label = f'ε = {coupling_value:.0e}'
        else:
            scaled_val = coupling_value * 10000
            coupling_label = f'ε = {scaled_val:.0f}×10⁻⁴'
        
        print(f"    {coupling_label:<25} {avg_tau:<12.1f} {avg_rate:<15.6f} {max_xi:<10.3f}")
    
    # Use the time normalization factor calculated from reference ε=0 slope
    print(f"    For plotting: using time normalization factor = {time_normalization_factor:.6f}")
    print(f"    This ensures reference ε=0 material time has slope 1.0")
    print(f"    Time axis: t' = t × {time_normalization_factor:.6f}")
    print(f"    Material time ξ keeps natural units from ∫(1/τ)dt")
    
    # Now plot all the material time curves with normalized time axis
    for coupling_name, data in material_time_data.items():
        coupling_value = data['coupling_value']
        time_data = data['time']
        exp_neg_material_time_diff_data = data['exp_neg_material_time_diff']
        
        # Apply time normalization and subtract time origin
        normalized_time_data = (time_data - time_data[0]) * time_normalization_factor
        
        color = get_coupling_color(coupling_value)
        
        if coupling_value == 0.0:
            label = r'$\varepsilon = 0$'
        else:
            # Format coupling value in scientific notation
            if coupling_value >= 1e-3:
                label = f'$\\varepsilon = {coupling_value:.0e}$'
            else:
                # Convert to the same format as F(k,t) analysis (×10^-4 notation)
                scaled_val = coupling_value * 10000
                label = f'$\\varepsilon = {scaled_val:.0f}\\times10^{{-4}}$'
        
        # Plot material time directly instead of stretched exponential
        material_time_diff_data = data['material_time_diff']
        
        # Subtract the origin to set the initial value to zero
        if len(material_time_diff_data) > 0:
            material_time_diff_zeroed = material_time_diff_data - material_time_diff_data[0]
        else:
            material_time_diff_zeroed = material_time_diff_data
        
        ax2.plot(normalized_time_data, material_time_diff_zeroed, color=color, linewidth=2, alpha=0.8, 
                linestyle='-', label=label)
    
    # Overlay reference material time data if available
    if waiting_time_data:
        print("    Overlaying reference material time evolution...")
        print(f"    Using time normalization factor from our ε=0 fictive data: a = {time_normalization_factor:.6f}")
        print(f"    All time axes: (t - t₀) × a to normalize our ε=0 slope to 1 and start at t=0")
        print(f"    Subtracting initial material time values to set ξ origin to zero")
        
        for coupling_name, data in waiting_time_data.items():
            coupling_value = data['coupling_value']
            ref_time = data['time']
            ref_material_time = data['material_time']
            
            color = get_coupling_color(coupling_value)
            
            # Apply the same time normalization that makes reference ε=0 have slope 1
            # Start from time origin (t=0 for reference data)
            normalized_ref_time = (ref_time - ref_time[0]) * time_normalization_factor
            
            # Filter to reasonable time range to avoid overcrowding
            max_time = 2000.0 * time_normalization_factor  # Adjust max time for normalized axis
            time_mask = normalized_ref_time <= max_time
            if np.any(time_mask):
                plot_times = normalized_ref_time[time_mask]
                plot_xi = ref_material_time[time_mask]
                
                # Subtract the origin to set the initial value to zero
                if len(plot_xi) > 0:
                    plot_xi_zeroed = plot_xi - plot_xi[0]
                else:
                    plot_xi_zeroed = plot_xi
                
                # Create reference material time label
                if coupling_value == 0.0:
                    ref_label = r'$\varepsilon = 0$ (Ref $\xi$)'
                else:
                    if coupling_value >= 1e-3:
                        ref_label = f'$\\varepsilon = {coupling_value:.0e}$ (Ref $\\xi$)'
                    else:
                        scaled_val = coupling_value * 10000
                        ref_label = f'$\\varepsilon = {scaled_val:.0f}\\times10^{{-4}}$ (Ref $\\xi$)'
                
                # Plot reference material time with dashed line to distinguish from our calculation
                ax2.plot(plot_times, plot_xi_zeroed, color=color, linewidth=1.5, alpha=0.6, 
                        linestyle='--', label=ref_label)
    
    # No horizontal reference line needed for material time comparison
    
    # Set labels and formatting for second panel
    ax2.set_xlabel(latex_format(r'Normalized Time (a×t)'))
    ax2.set_ylabel(latex_format(r'Material Time $\xi(t)$', 'Material Time ξ(t)'))
    
    if waiting_time_data:
        title_text = 'Material Time Analysis: Our Calculation vs Reference'
        title_latex = r'Material Time Analysis: Our Calculation vs Reference'
    else:
        title_text = 'Material Time Analysis: exp(-(ξ(t+tw) - ξ(tw))^0.55) vs Time'
        title_latex = r'Material Time Analysis: $\exp(-(\xi(t+t_w) - \xi(t_w))^{0.55})$ vs Time'
    
    ax2.set_title(latex_format(title_latex, title_text))
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, None)
    
    # Use linear y-scale with automatic limits
    ax2.set_ylim(bottom=0)  # Start from 0 but let y-axis expand as needed
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    output_file = Path(output_dir) / 'fictive_vs_waiting_relaxation_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    
    # Also save as PDF
    output_file_pdf = Path(output_dir) / 'fictive_vs_waiting_relaxation_comparison.pdf'
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file_pdf}")
    
    #plt.show()
    
    # Reset matplotlib settings
    plt.rcParams.update(plt.rcParamsDefault)

def main():
    """Main function."""
    
    # Directories
    base_dir = '/media/extradrive/Trajectories/final_nodiss_cavitymd'
    time_series_dir = Path(base_dir) / 'time_series_output'
    output_dir = base_dir
    
    print("=" * 60)
    print("FICTIVE vs WAITING TIME RELAXATION COMPARISON")
    print("=" * 60)
    
    # 1. Load fictive relaxation time data
    print("\n1. Loading fictive relaxation time data...")
    fictive_data = load_fictive_data(time_series_dir)
    
    if not fictive_data:
        print("Error: No fictive data found!")
        return 1
    
    print(f"Loaded fictive data for {len(fictive_data)} coupling strengths")
    
    # 2. Load material time reference data for comparison
    print("\n2. Loading material time reference data for comparison...")
    material_time_ref_file = os.path.join(base_dir, "material_time_all_couplings.txt")
    material_time_ref_data = load_material_time_reference_data(material_time_ref_file)
    
    if not material_time_ref_data:
        print("Warning: No material time reference data found! Proceeding with fictive analysis only...")
        waiting_time_data = {}
    else:
        print(f"Loaded material time reference data for {len(material_time_ref_data)} coupling strengths")
        waiting_time_data = material_time_ref_data  # Use this for comparison
    
    # Skip waiting time analysis validation since we're focusing on material time
    print("Proceeding to material time analysis plot...")
    
    # 3. Create material time analysis plot
    print("\n3. Creating comparison plot...")
    plot_fictive_vs_waiting_comparison(fictive_data, waiting_time_data, output_dir)
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print("=" * 60)
    print("\nGenerated plots:")
    print("  - fictive_vs_waiting_relaxation_comparison.png")
    print("  - fictive_vs_waiting_relaxation_comparison.pdf")
    print("\nPanel 1: Fictive relaxation time (continuous, normalized by first 200 ps) vs F(k,t) relaxation time analysis (discrete waiting times)")
    if waiting_time_data:
        print("Panel 2: Material time evolution comparison")
        print("         Solid lines: ξ(t+tw) - ξ(tw) from our fictive temperature analysis")
        print("         Dashed lines: Reference material time ξ(t) from Tool-Narayanaswamy analysis")
        print("         Our ξ(t): Material time difference with tw = 400 ps (200 ps after coupling)")
        print("         Reference ξ(t): Material time from F(k,t) analysis with criterion F(k,t) = 0.1")
        print("         Both calculated as ξ(t) = ∫ (1/τ) dt but using different τ sources")
    else:
        print("Panel 2: Material time analysis - ξ(t) - ξ(tw) vs time, where tw = 400 ps")
        print("         ξ(t) = Ξ(t) × scaling_factor, Ξ(t) = ∫₄₀₀ᵗ (1/τ_fictive(t')) dt'")
        print("         Material time difference starting from tw = 400 ps (200 ps after coupling)")
    
    return 0

if __name__ == "__main__":
    exit(main())
