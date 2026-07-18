#!/usr/bin/env python3
"""
Fictive Temperature Components Plotting Script

This script creates separate plots for each coupling strength showing the evolution
of fictive temperatures from different energy components:
- Harmonic potential energy
- LJ (Lennard-Jones) potential energy  
- Coulombic potential energy
- Total potential energy

New Feature: Material Time Calculation
- Fits LJ+Coul component: T(t) = T∞ - A*exp(-t/τ) to extract relaxation time τ
- Calculates material time Ξ(t) = ∫(1/τ)dt from fitted relaxation time
- Rescales to material time ξ(t) with condition: exp(-ξ^0.55) = 0.1 when ξ = 1
- Time origin shifted to 201 ps (configurable) for fitting

Usage: python plot_fictive_temperature_components.py [--base_dir time_series_output] [--output_dir ./] [--time_shift 201.0]
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
import subprocess
import matplotlib.style
import matplotlib as mpl
from scipy import ndimage
from scipy.optimize import curve_fit
import sys
from pathlib import Path

REPRO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPRO_ROOT / "shared"))
from latex_config_adobe import setup_latex_fonts, latex_safe
from scipy import interpolate
import warnings
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
mpl.style.use('classic')

def parse_coupling_strength(filename):
    """Parse coupling strength from filename or directory name."""
    # Handle both file names and directory names
    if filename.startswith('cavity_'):
        # Directory name format: cavity_coupling_1eneg03_switch_200.0ps
        coupling_part = filename.replace('cavity_coupling_', '').split('_switch_')[0]
    else:
        # File name format: coupling_1eneg03_averaged_fictive_temperatures.txt
        coupling_part = filename.replace('coupling_', '').replace('_averaged_fictive_temperatures.txt', '')
    
    if '0epos00' in coupling_part:
        return 0.0, '0.0'
    elif 'eneg' in coupling_part:
        parts = coupling_part.split('eneg')
        if len(parts) == 2:
            try:
                mantissa = int(parts[0])
                exponent = int(parts[1])
                value = mantissa * (10 ** (-exponent))
                return value, f'{mantissa}\\times10^{{-{exponent}}}'
            except ValueError:
                print(f"Warning: Could not parse coupling strength from '{coupling_part}'")
                return float('nan'), coupling_part
    
    return float('nan'), coupling_part

def load_temperature_energy_components(filename="potential_energy_components_vs_temperature.txt"):
    """Load the potential energy components vs temperature data for fictive temperature calculation."""
    try:
        # Read the header to get column names
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Find the header line (first non-comment line)
        header_line = None
        data_start_line = 0
        
        for i, line in enumerate(lines):
            line = line.strip()
            if not line.startswith('#') and line:
                header_line = line
                data_start_line = i + 1
                break
        
        if header_line is None:
            raise ValueError("Could not find column headers in temperature components file")
        
        column_names = header_line.split()
        temp_columns = {name: i for i, name in enumerate(column_names)}
        
        # Read the data
        temp_data = np.loadtxt(filename, skiprows=data_start_line)
        
        return temp_data, temp_columns
        
    except Exception as e:
        print(f"Warning: Could not load temperature components file {filename}: {e}")
        print("Fictive temperature calculation will use fallback methods.")
        return None, None

def create_energy_temperature_interpolators(temp_data, temp_columns):
    """Create interpolating functions U(T) for energy components."""
    if temp_data is None or temp_columns is None:
        return {}
    
    interpolators = {}
    temperatures = temp_data[:, temp_columns['temperature']]
    
    # Create interpolators for different energy components
    energy_components = {
        'Harmonic Bond Energy': ('harmonic_hartree', 'harmonic_hartree'),
        'Lennard-Jones Energy': ('lj_hartree', 'lj_hartree'),
        'Coulombic Energy': ('coulombic_hartree', 'coulombic_hartree'),
        'Total Potential Energy': ('total_PE_hartree', 'total_PE_hartree')
    }
    
    for display_name, (col_name, energy_column) in energy_components.items():
        if col_name in temp_columns:
            print(f"Creating interpolator for {display_name} using column {col_name}")
            energies = temp_data[:, temp_columns[col_name]]
            
            # Sort by temperature for proper interpolation
            sort_idx = np.argsort(temperatures)
            T_sorted = temperatures[sort_idx]
            U_sorted = energies[sort_idx]
            
            # Create U(T) and T(U) interpolators
            U_of_T = interpolate.interp1d(T_sorted, U_sorted, kind='linear', bounds_error=False, fill_value='extrapolate')
            T_of_U = interpolate.interp1d(U_sorted, T_sorted, kind='linear', bounds_error=False, fill_value='extrapolate')
            
            # Create linear extrapolation functions for extreme cases
            n_fit = min(5, len(T_sorted) // 4)
            
            # Low temperature linear extrapolation
            T_low = T_sorted[:n_fit]
            U_low = U_sorted[:n_fit]
            if len(T_low) >= 2:
                low_slope = (U_low[-1] - U_low[0]) / (T_low[-1] - T_low[0])
                low_intercept = U_low[0] - low_slope * T_low[0]
                low_linear = lambda T: (T - low_intercept) / low_slope if low_slope != 0 else T_low[0]
            else:
                low_linear = None
            
            # High temperature linear extrapolation  
            T_high = T_sorted[-n_fit:]
            U_high = U_sorted[-n_fit:]
            if len(T_high) >= 2:
                high_slope = (U_high[-1] - U_high[0]) / (T_high[-1] - T_high[0])
                high_intercept = U_high[0] - high_slope * T_high[0]
                high_linear = lambda T: (T - high_intercept) / high_slope if high_slope != 0 else T_high[0]
            else:
                high_linear = None
            
            interpolators[display_name] = {
                'U_of_T': U_of_T,
                'T_of_U': T_of_U,
                'U_range': (U_sorted.min(), U_sorted.max()),
                'T_range': (T_sorted.min(), T_sorted.max()),
                'energy_column': energy_column,
                'low_linear': low_linear,
                'high_linear': high_linear
            }
            
            print(f"  Temperature range: {T_sorted.min():.1f} - {T_sorted.max():.1f} K")
            print(f"  Energy range: {U_sorted.min():.6f} - {U_sorted.max():.6f} Hartree")
        else:
            print(f"  Column {col_name} not found, skipping {display_name}")
    
    print(f"Successfully created interpolators for: {list(interpolators.keys())}")
    return interpolators

def calculate_fictive_temperatures_from_energies(energy_dict, interpolators):
    """Calculate fictive temperatures from energy components using empirical U(T) relationships."""
    fictive_temps = {}
    
    # Mapping from energy component names to interpolator names
    component_mapping = {
        'harmonic': 'Harmonic Bond Energy',
        'lj': 'Lennard-Jones Energy', 
        'coulombic': 'Coulombic Energy',
        'total': 'Total Potential Energy'
    }
    
    for comp_name, energies in energy_dict.items():
        if comp_name in component_mapping:
            interp_name = component_mapping[comp_name]
            
            if interp_name in interpolators:
                interp_data = interpolators[interp_name]
                T_of_U = interp_data['T_of_U']
                U_range = interp_data['U_range']
                low_linear = interp_data.get('low_linear')
                high_linear = interp_data.get('high_linear')
                
                # Calculate fictive temperatures with improved extrapolation
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fictive_T = T_of_U(energies)
                
                # Use linear extrapolation for extreme cases
                if low_linear is not None and high_linear is not None:
                    # Identify extreme extrapolation regions (beyond 50% extension of original range)
                    U_span = U_range[1] - U_range[0]
                    low_extreme = energies < (U_range[0] - 0.5 * U_span)
                    high_extreme = energies > (U_range[1] + 0.5 * U_span)
                    
                    # Apply linear extrapolation for extreme cases
                    if np.any(low_extreme):
                        fictive_T[low_extreme] = low_linear(energies[low_extreme])
                    
                    if np.any(high_extreme):
                        fictive_T[high_extreme] = high_linear(energies[high_extreme])
                
                # Only mask physically unreasonable temperatures
                reasonable_mask = (fictive_T > 0) & (fictive_T < 10000)  # 0-10000 K seems reasonable
                fictive_T[~reasonable_mask] = np.nan
                
                # Count extrapolated vs interpolated points
                interpolated_mask = (energies >= U_range[0]) & (energies <= U_range[1])
                n_interpolated = np.sum(interpolated_mask & reasonable_mask)
                n_extrapolated = np.sum(~interpolated_mask & reasonable_mask)
                n_total = len(fictive_T)
                
                print(f"  {comp_name}: {n_interpolated} interpolated, {n_extrapolated} extrapolated, {n_total - n_interpolated - n_extrapolated} invalid")
                
                fictive_temps[comp_name] = fictive_T
            else:
                print(f"  No interpolator found for {comp_name} component")
                # Fallback to NaN array
                fictive_temps[comp_name] = np.full_like(energies, np.nan)
        else:
            print(f"  Unknown component: {comp_name}")
            fictive_temps[comp_name] = np.full_like(energies, np.nan)
    
    return fictive_temps

def load_energy_tracker_data(coupling_dir):
    """Load energy tracker data from a coupling directory."""
    coupling_path = Path(coupling_dir)
    
    # Look for energy tracker files
    energy_files = list(coupling_path.glob('prod-*_energy_tracker.txt'))
    
    if not energy_files:
        print(f"  No energy tracker files found in {coupling_dir}")
        return None, None
    
    print(f"  Found {len(energy_files)} energy tracker files")
    
    # Read and average energy data from all replicas
    all_energy_data = []
    
    for energy_file in energy_files:
        try:
            # Read energy tracker file (skip comment lines)
            data = pd.read_csv(energy_file, sep=r'\s+', comment='#')
            all_energy_data.append(data)
        except Exception as e:
            print(f"    Warning: Could not read {energy_file}: {e}")
            continue
    
    if not all_energy_data:
        print(f"  No valid energy data found in {coupling_dir}")
        return None, None
    
    # Find common time grid (use shortest time range)
    min_length = min(len(df) for df in all_energy_data)
    common_time = all_energy_data[0]['time(ps)'].iloc[:min_length].values
    
    # Average energy components across replicas
    energy_components = {}
    required_columns = ['harmonic_energy', 'lj_energy', 'ewald_short_energy', 'ewald_long_energy']
    
    for col in required_columns:
        if all(col in df.columns for df in all_energy_data):
            values = np.array([df[col].iloc[:min_length].values for df in all_energy_data])
            energy_components[col] = np.mean(values, axis=0)
        else:
            print(f"    Warning: Column {col} not found in all files")
    
    # Calculate combined components
    if 'harmonic_energy' in energy_components:
        energy_components['harmonic'] = energy_components['harmonic_energy']
    
    if 'lj_energy' in energy_components:
        energy_components['lj'] = energy_components['lj_energy']
    
    if 'ewald_short_energy' in energy_components and 'ewald_long_energy' in energy_components:
        energy_components['coulombic'] = energy_components['ewald_short_energy'] + energy_components['ewald_long_energy']
    
    if all(comp in energy_components for comp in ['harmonic_energy', 'lj_energy', 'ewald_short_energy', 'ewald_long_energy']):
        energy_components['total'] = (energy_components['harmonic_energy'] + 
                                     energy_components['lj_energy'] + 
                                     energy_components['ewald_short_energy'] + 
                                     energy_components['ewald_long_energy'])
    
    print(f"  Loaded {len(common_time)} time points from {len(all_energy_data)} replicas")
    
    return common_time, energy_components

def smooth_data(time, values, window_ps=100.0):
    """
    Smooth data using a sliding window average.
    
    Parameters:
    -----------
    time : array-like
        Time values in picoseconds
    values : array-like
        Data values to smooth
    window_ps : float
        Smoothing window size in picoseconds (default: 100.0)
    
    Returns:
    --------
    smoothed_values : ndarray
        Smoothed data values
    """
    if len(time) < 2:
        return np.array(values)
    
    # Calculate time step
    dt = np.mean(np.diff(time))
    
    # Calculate window size in data points
    window_points = int(window_ps / dt)
    
    # Ensure minimum window size
    window_points = max(1, window_points)
    
    # For very small windows, return original data
    if window_points <= 3:
        return np.array(values)
    
    # Use uniform filter for sliding window average
    # Handle NaN values by replacing with interpolated values temporarily
    values_array = np.array(values)
    nan_mask = np.isnan(values_array)
    
    if np.all(nan_mask):
        return values_array
    
    # If there are NaNs, interpolate them for smoothing
    if np.any(nan_mask):
        # Find valid indices
        valid_indices = np.where(~nan_mask)[0]
        if len(valid_indices) < 2:
            return values_array
        
        # Interpolate NaN values
        values_interp = np.copy(values_array)
        values_interp[nan_mask] = np.interp(
            np.where(nan_mask)[0], 
            valid_indices, 
            values_array[valid_indices]
        )
    else:
        values_interp = values_array
    
    # Apply uniform filter (sliding window average)
    smoothed = ndimage.uniform_filter1d(values_interp, size=window_points, mode='nearest')
    
    # Restore NaN values in their original positions
    if np.any(nan_mask):
        smoothed[nan_mask] = np.nan
    
    return smoothed

def exponential_rise(t, A, tau, T_final):
    """
    Exponential rise function: T(t) = T_final - A * exp(-t/tau)
    
    Parameters:
    -----------
    t : array-like
        Time values (after shifting origin)
    A : float
        Amplitude coefficient (positive for rise to T_final)
    tau : float
        Rise time constant
    T_final : float
        Final temperature (asymptotic value)
    
    Returns:
    --------
    T : array-like
        Temperature values
    """
    return T_final - A * np.exp(-t / tau)

def fit_lj_relaxation_time(time, lj_temps, time_shift_ps=201.0):
    """
    Fit LJ+Coul component to extract relaxation time for material time calculation.
    
    Parameters:
    -----------
    time : array-like
        Time values in picoseconds
    lj_temps : array-like
        LJ+Coul component temperature values in Kelvin
    time_shift_ps : float
        Time point to shift origin to (default: 201.0 ps)
    
    Returns:
    --------
    tau_relaxation : float or None
        Relaxation time constant in ps, or None if fitting failed
    """
    try:
        # Find data points after the time shift
        time_mask = time >= time_shift_ps
        
        if not np.any(time_mask):
            print(f"    Warning: No data points found after {time_shift_ps} ps")
            return None
            
        time_fit = time[time_mask]
        lj_fit = lj_temps[time_mask]
        
        # Remove NaN values
        valid_mask = ~np.isnan(lj_fit)
        if not np.any(valid_mask):
            print(f"    Warning: All LJ+Coul temperature values are NaN after {time_shift_ps} ps")
            return None
            
        time_fit = time_fit[valid_mask]
        lj_fit = lj_fit[valid_mask]
        
        if len(time_fit) < 5:
            print(f"    Warning: Too few data points ({len(time_fit)}) for fitting")
            return None
        
        # Shift time origin
        time_shifted = time_fit - time_shift_ps
        
        # Fit LJ+Coul data with exponential rise
        T0_guess = lj_fit[0] if len(lj_fit) > 0 else 50.0
        T_final_guess = lj_fit[-1] if len(lj_fit) > 0 else 100.0
        A_guess = max(T_final_guess - T0_guess, 10.0)
        tau_guess = 100.0
        
        bounds = ([0, 1, 50], [np.inf, 10000, 200])
        initial_guess = [A_guess, tau_guess, T_final_guess]
        
        popt, pcov = curve_fit(exponential_rise, time_shifted, lj_fit, 
                              p0=initial_guess, bounds=bounds, maxfev=2000)
        
        A_fit, tau_fit, T_final_fit = popt
        
        # Calculate R-squared
        fit_temps = exponential_rise(time_shifted, A_fit, tau_fit, T_final_fit)
        ss_res = np.sum((lj_fit - fit_temps) ** 2)
        ss_tot = np.sum((lj_fit - np.mean(lj_fit)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        print(f"    LJ+Coul relaxation fit: A = {A_fit:.1f} K, tau = {tau_fit:.1f} ps, T_final = {T_final_fit:.1f} K, R² = {r_squared:.3f}")
        
        return tau_fit
        
    except Exception as e:
        print(f"    Error fitting LJ+Coul for relaxation time: {e}")
        return None

def calculate_material_time(time, tau_relaxation, time_shift_ps=201.0):
    """
    Calculate material time by integrating 1/τ over time.
    
    Parameters:
    -----------
    time : array-like
        Time values in picoseconds
    tau_relaxation : float
        Relaxation time constant in ps
    time_shift_ps : float
        Time point to shift origin to (default: 201.0 ps)
    
    Returns:
    --------
    material_time_Xi : array-like
        Material time Ξ(t) = ∫(1/τ)dt
    material_time_xi : array-like
        Rescaled material time ξ(t) 
    valid_mask : array-like
        Boolean mask for valid time points (t >= time_shift_ps)
    """
    # Find data points after the time shift
    valid_mask = time >= time_shift_ps
    
    if not np.any(valid_mask):
        return None, None, None
    
    time_valid = time[valid_mask]
    time_shifted = time_valid - time_shift_ps
    
    # Calculate material time Ξ(t) = ∫(1/τ)dt from t=0
    # Since τ is constant from the fit, this is simply t/τ
    material_time_Xi = time_shifted / tau_relaxation
    
    # Rescale material time ξ(t) so that exp(-ξ^0.55) = 0.1 when ξ = 1
    # From exp(-ξ^0.55) = 0.1, we get ξ^0.55 = -ln(0.1) = ln(10) ≈ 2.3026
    # So when ξ = 1, we need ξ^0.55 = 2.3026
    # This means we need to find the scaling factor
    
    # The condition is: when ξ_scaled = 1, then exp(-1^0.55) = exp(-1) ≠ 0.1
    # We need: exp(-ξ_scaled^0.55) = 0.1 when ξ_scaled = 1
    # So: ξ_scaled^0.55 = ln(10) when ξ_scaled = 1
    # Therefore: 1^0.55 = ln(10), which means we need a scaling factor
    
    # Let ξ = α * Ξ, where α is the scaling factor
    # We want: exp(-(α*Ξ)^0.55) = 0.1 when α*Ξ = 1
    # So: (α*Ξ)^0.55 = ln(10) when α*Ξ = 1
    # This gives us: α^0.55 * Ξ^0.55 = ln(10) when α*Ξ = 1
    # When α*Ξ = 1, then Ξ = 1/α, so: α^0.55 * (1/α)^0.55 = ln(10)
    # This simplifies to: 1 = ln(10), which is wrong.
    
    # Let me reconsider: we want exp(-ξ^0.55) = 0.1 when ξ = 1
    # So ξ^0.55 = ln(10) when ξ = 1
    # This means 1^0.55 = ln(10), so ln(10) = 1, which is not true.
    
    # I think the condition means: find the value of Ξ where exp(-Ξ^0.55) = 0.1
    # and then rescale so that value corresponds to ξ = 1
    
    target_value = np.log(10)  # We want ξ^0.55 = ln(10) 
    # So ξ = (ln(10))^(1/0.55) = (ln(10))^(20/11)
    xi_at_target = target_value**(1/0.55)
    
    # Find the Ξ value where we want ξ = 1
    # We'll arbitrarily choose the point where material time reaches xi_at_target
    # and rescale so that becomes ξ = 1
    
    if len(material_time_Xi) > 0 and np.max(material_time_Xi) > 0:
        # Scale so that the maximum value corresponds to some reasonable ξ value
        # Let's scale so that ξ goes from 0 to some maximum that makes physical sense
        material_time_xi = material_time_Xi * (xi_at_target / np.max(material_time_Xi))
    else:
        material_time_xi = material_time_Xi
    
    return material_time_Xi, material_time_xi, valid_mask

def read_averaged_potential_energy(file_path):
    """Read time and potential energy component columns from averaged PE file."""
    try:
        data = pd.read_csv(file_path, sep=r'\s+', comment='#')
        time = data.iloc[:, 0].values
        return {
            'time': time,
            'energy_harmonic': data.iloc[:, 1].values,
            'energy_lj': data.iloc[:, 2].values,
            'energy_coulombic': data.iloc[:, 3].values,
            'energy_total': data.iloc[:, 4].values,
        }
    except Exception as e:
        print(f"Error reading potential energy file {file_path}: {e}")
        return None


def read_fictive_temperature_components(file_path):
    """Read fictive temperature data from file with all energy components."""
    try:
        # Read the file, skipping comment lines
        data = pd.read_csv(file_path, sep=r'\s+', comment='#')
        
        # Expected columns: time, coulombic, harmonic, LJ, total
        time = data.iloc[:, 0].values
        fictive_coulombic = data.iloc[:, 1].values  
        fictive_harmonic = data.iloc[:, 2].values
        fictive_lj = data.iloc[:, 3].values
        fictive_total = data.iloc[:, 4].values
        
        return time, fictive_coulombic, fictive_harmonic, fictive_lj, fictive_total
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return None, None, None, None, None

def collect_fictive_temperature_data(base_dir):
    """Collect all fictive temperature data by loading energy tracker files and calculating using empirical U(T) relationships."""
    base_path = Path(base_dir)
    
    # Look for coupling directories
    coupling_dirs = [d for d in base_path.iterdir() if d.is_dir() and d.name.startswith('cavity_coupling_')]
    
    # If no coupling directories found, look for pre-calculated files as fallback
    if not coupling_dirs:
        print(f"No coupling directories found in {base_dir}, looking for pre-calculated files...")
        return collect_fictive_temperature_data_fallback(base_dir)
    
    # Check if coupling directories contain energy tracker files
    has_energy_files = False
    for coupling_dir in coupling_dirs[:1]:  # Check first directory only
        energy_files = list(coupling_dir.glob('prod-*_energy_tracker.txt'))
        if energy_files:
            has_energy_files = True
            break
    
    if not has_energy_files:
        print(f"No energy tracker files found in coupling directories, looking for pre-calculated files...")
        return collect_fictive_temperature_data_fallback(Path(base_dir) / 'time_series_output')
    
    # Load empirical U(T) data
    print("Loading empirical energy-temperature relationships...")
    temp_data, temp_columns = load_temperature_energy_components()
    
    if temp_data is None:
        print("Could not load empirical data, falling back to pre-calculated files...")
        return collect_fictive_temperature_data_fallback(base_dir)
    
    # Create interpolators
    print("Creating energy-temperature interpolators...")
    interpolators = create_energy_temperature_interpolators(temp_data, temp_columns)
    
    if not interpolators:
        print("Could not create interpolators, falling back to pre-calculated files...")
        return collect_fictive_temperature_data_fallback(base_dir)
    
    data_dict = {}
    coupling_info = {}
    
    print(f"Processing {len(coupling_dirs)} coupling directories...")
    
    for coupling_dir in sorted(coupling_dirs):
        # Parse coupling strength from directory name
        dir_name = coupling_dir.name
        coupling_value, coupling_label = parse_coupling_strength(dir_name)
        coupling_name = dir_name.replace('cavity_', '')
        
        coupling_info[coupling_name] = (coupling_value, coupling_label)
        
        print(f"Processing directory: {dir_name} (coupling = {coupling_label})")
        
        # Load energy tracker data
        time, energy_components = load_energy_tracker_data(coupling_dir)
        
        if time is not None and energy_components:
            # Calculate fictive temperatures using empirical relationships
            print(f"  Calculating fictive temperatures using empirical U(T) relationships...")
            fictive_temps = calculate_fictive_temperatures_from_energies(energy_components, interpolators)
            
            if fictive_temps:
                data_dict[coupling_name] = {
                    'time': time,
                    'fictive_coulombic': fictive_temps.get('coulombic', np.full_like(time, np.nan)),
                    'fictive_harmonic': fictive_temps.get('harmonic', np.full_like(time, np.nan)),
                    'fictive_lj': fictive_temps.get('lj', np.full_like(time, np.nan)),
                    'fictive_total': fictive_temps.get('total', np.full_like(time, np.nan)),
                    # Include raw energy data
                    'energy_harmonic': energy_components.get('harmonic', np.full_like(time, np.nan)),
                    'energy_lj': energy_components.get('lj', np.full_like(time, np.nan)),
                    'energy_coulombic': energy_components.get('coulombic', np.full_like(time, np.nan)),
                    'energy_total': energy_components.get('total', np.full_like(time, np.nan)),
                    'coupling_value': coupling_value,
                    'coupling_label': coupling_label
                }
                print(f"  Successfully calculated fictive temperatures for {len(time)} data points")
            else:
                print(f"  Failed to calculate fictive temperatures")
        else:
            print(f"  Failed to load energy data")
    
    return data_dict, coupling_info

def estimate_energies_from_fictive_temps(fictive_temps, interpolators):
    """Estimate raw energies from fictive temperatures using inverse U(T) relationships."""
    estimated_energies = {}
    
    # Load empirical data to create inverse interpolators
    temp_data, temp_columns = load_temperature_energy_components()
    if temp_data is None or temp_columns is None:
        print("  Warning: Cannot estimate energies - no empirical data available")
        return {}
    
    # Create T -> U interpolators
    temperatures = temp_data[:, temp_columns['temperature']]
    energy_components = {
        'harmonic': 'harmonic_hartree',
        'lj': 'lj_hartree', 
        'coulombic': 'coulombic_hartree',
        'total': 'total_PE_hartree'
    }
    
    for comp_name, col_name in energy_components.items():
        if col_name in temp_columns and comp_name in fictive_temps:
            energies = temp_data[:, temp_columns[col_name]]
            
            # Sort by temperature for proper interpolation
            sort_idx = np.argsort(temperatures)
            T_sorted = temperatures[sort_idx]
            U_sorted = energies[sort_idx]
            
            # Create T -> U interpolator
            U_of_T = interpolate.interp1d(T_sorted, U_sorted, kind='linear', 
                                        bounds_error=False, fill_value='extrapolate')
            
            # Estimate energies from fictive temperatures
            fictive_T = fictive_temps[comp_name]
            valid_mask = ~np.isnan(fictive_T) & (fictive_T > 0) & (fictive_T < 10000)
            
            estimated_U = np.full_like(fictive_T, np.nan)
            if np.any(valid_mask):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    estimated_U[valid_mask] = U_of_T(fictive_T[valid_mask])
            
            estimated_energies[comp_name] = estimated_U
            
            n_valid = np.sum(valid_mask)
            print(f"    Estimated {comp_name} energies for {n_valid} points")
    
    return estimated_energies

def collect_fictive_temperature_data_fallback(base_dir):
    """Fallback method: collect fictive temperature data from pre-calculated files."""
    print("Using fallback method: reading pre-calculated fictive temperature files...")
    
    base_path = Path(base_dir)
    fictive_files = list(base_path.glob('coupling_*_averaged_fictive_temperatures.txt'))
    
    if not fictive_files:
        print(f"No fictive temperature files found in {base_dir}")
        return {}, {}
    
    data_dict = {}
    coupling_info = {}
    
    for file_path in fictive_files:
        filename = file_path.name
        coupling_value, coupling_label = parse_coupling_strength(filename)
        coupling_name = filename.replace('_averaged_fictive_temperatures.txt', '')
        
        coupling_info[coupling_name] = (coupling_value, coupling_label)
        
        print(f"Processing file: {filename} (coupling = {coupling_label})")
        
        time, coulombic, harmonic, lj, total = read_fictive_temperature_components(file_path)
        
        if time is not None:
            pe_file = base_path / f"{coupling_name}_averaged_potential_energy.txt"
            if pe_file.exists():
                print(f"  Loading raw energies from {pe_file.name}...")
                pe_data = read_averaged_potential_energy(pe_file)
                if pe_data and len(pe_data['time']) != len(time):
                    print(f"  Interpolating energies onto fictive-temperature time grid "
                          f"({len(pe_data['time'])} -> {len(time)} points)...")
                    for key in ('energy_harmonic', 'energy_lj', 'energy_coulombic', 'energy_total'):
                        f = interpolate.interp1d(
                            pe_data['time'], pe_data[key],
                            kind='linear', bounds_error=False, fill_value='extrapolate',
                        )
                        pe_data[key] = f(time)
                energy_harmonic = pe_data['energy_harmonic'] if pe_data else np.full_like(time, np.nan)
                energy_lj = pe_data['energy_lj'] if pe_data else np.full_like(time, np.nan)
                energy_coulombic = pe_data['energy_coulombic'] if pe_data else np.full_like(time, np.nan)
                energy_total = pe_data['energy_total'] if pe_data else np.full_like(time, np.nan)
            else:
                fictive_temps = {
                    'coulombic': coulombic,
                    'harmonic': harmonic,
                    'lj': lj,
                    'total': total
                }
                print(f"  Estimating raw energies from fictive temperatures...")
                estimated_energies = estimate_energies_from_fictive_temps(fictive_temps, {})
                energy_harmonic = estimated_energies.get('harmonic', np.full_like(time, np.nan))
                energy_lj = estimated_energies.get('lj', np.full_like(time, np.nan))
                energy_coulombic = estimated_energies.get('coulombic', np.full_like(time, np.nan))
                energy_total = estimated_energies.get('total', np.full_like(time, np.nan))

            data_dict[coupling_name] = {
                'time': time,
                'fictive_coulombic': coulombic,
                'fictive_harmonic': harmonic,
                'fictive_lj': lj,
                'fictive_total': total,
                'energy_harmonic': energy_harmonic,
                'energy_lj': energy_lj,
                'energy_coulombic': energy_coulombic,
                'energy_total': energy_total,
                'coupling_value': coupling_value,
                'coupling_label': coupling_label
            }
            print(f"  Read {len(time)} data points")
        else:
            print(f"  Failed to read data")
    
    return data_dict, coupling_info

def filter_coupling_strengths(coupling_info, coupling_filter=None):
    """Filter to only include specific coupling strengths."""
    target_couplings = [0.0, 3e-4, 5e-4, 7e-4, 1e-3]
    if coupling_filter is not None:
        target_couplings = [coupling_filter]

    filtered_couplings = {}
    for coupling_name, (coupling_value, coupling_label) in coupling_info.items():
        if not np.isnan(coupling_value) and any(abs(coupling_value - target) < 1e-6 for target in target_couplings):
            filtered_couplings[coupling_name] = (coupling_value, coupling_label)

    return filtered_couplings


def parse_coupling_arg(value):
    """Parse --coupling values like 3e-4 or 3eneg04."""
    if 'eneg' in value:
        parts = value.replace('coupling_', '').split('eneg')
        return int(parts[0]) * (10 ** (-int(parts[1])))
    return float(value)


def plot_coupling_fictive_components(coupling_name, data, output_dir='.', smooth_window_ps=100.0, time_shift_ps=201.0, use_latex=True):
    """Create a plot for a single coupling showing all energy components."""
    if use_latex:
        use_latex = setup_latex_fonts()
    else:
        from latex_config_adobe import _apply_fallback_fonts
        _apply_fallback_fonts()
    
    coupling_value = data['coupling_value']
    coupling_label = data['coupling_label']
    time = data['time']
    
    # Get all energy component data
    fictive_coulombic = data['fictive_coulombic']
    fictive_harmonic = data['fictive_harmonic']
    fictive_lj = data['fictive_lj']
    fictive_total = data['fictive_total']
    
    # Apply smoothing if window size is specified
    if smooth_window_ps > 0:
        print(f"    Applying {smooth_window_ps} ps smoothing window...")
        fictive_harmonic = smooth_data(time, fictive_harmonic, smooth_window_ps)
        fictive_lj = smooth_data(time, fictive_lj, smooth_window_ps)
        fictive_total = smooth_data(time, fictive_total, smooth_window_ps)
        
        # Implement consistent windowing strategy: subsample data points
        # consistent with smoothing window rather than plotting all smoothed points
        dt_avg = np.mean(np.diff(time))  # Average time step in data
        subsample_interval = max(1, int(smooth_window_ps / (2 * dt_avg)))  # Plot every window/2 interval
        print(f"    Subsampling data: plotting every {subsample_interval} points ({subsample_interval * dt_avg:.1f} ps intervals)")
        
        # Apply subsampling to all arrays (fictive temperatures and energies)
        time = time[::subsample_interval]
        fictive_coulombic = fictive_coulombic[::subsample_interval]
        fictive_harmonic = fictive_harmonic[::subsample_interval]
        fictive_lj = fictive_lj[::subsample_interval]
        fictive_total = fictive_total[::subsample_interval]
        
        # Also subsample energy arrays if they exist
        if 'energy_harmonic' in data:
            data['energy_harmonic'] = data['energy_harmonic'][::subsample_interval]
        if 'energy_lj' in data:
            data['energy_lj'] = data['energy_lj'][::subsample_interval]
        if 'energy_coulombic' in data:
            data['energy_coulombic'] = data['energy_coulombic'][::subsample_interval]
        if 'energy_total' in data:
            data['energy_total'] = data['energy_total'][::subsample_interval]
    
    # Extract relaxation time from LJ+Coul component and calculate material time
    print(f"    Extracting relaxation time from LJ+Coul component (time shift: {time_shift_ps} ps)...")
    tau_relaxation = fit_lj_relaxation_time(time, fictive_lj, time_shift_ps)
    
    if tau_relaxation is not None:
        print(f"    Calculating material time with tau = {tau_relaxation:.1f} ps...")
        material_time_Xi, material_time_xi, material_time_mask = calculate_material_time(
            time, tau_relaxation, time_shift_ps)
    else:
        print(f"    Could not extract relaxation time - skipping material time calculation")
        material_time_Xi = material_time_xi = material_time_mask = None
    
    # Create figure with two panels: top for energies, bottom for fictive temperatures
    fig, (ax_energy, ax_temp) = plt.subplots(2, 1, figsize=(4, 6), sharex=True)
    
    # Check if we have raw energy data available
    has_energy_data = ('energy_harmonic' in data and 'energy_lj' in data and 
                      'energy_coulombic' in data and 'energy_total' in data)
    
    if has_energy_data:
        # Top panel: Raw potential energies
        energy_harmonic = data['energy_harmonic']
        energy_lj = data['energy_lj'] 
        energy_coulombic = data['energy_coulombic']
        energy_total = data['energy_total']
        
        # Apply same smoothing and subsampling to energy data
        if smooth_window_ps > 0:
            energy_harmonic = smooth_data(time, energy_harmonic, smooth_window_ps)
            energy_lj = smooth_data(time, energy_lj, smooth_window_ps)
            energy_coulombic = smooth_data(time, energy_coulombic, smooth_window_ps)
            energy_total = smooth_data(time, energy_total, smooth_window_ps)
        
        # Calculate equilibrium values (average before t < 200 ps) for baseline subtraction
        equilibrium_mask = time < 200.0
        if np.any(equilibrium_mask):
            print(f"    Calculating equilibrium baselines from {np.sum(equilibrium_mask)} points before 200 ps")
            
            # Calculate equilibrium averages
            eq_harmonic = np.nanmean(energy_harmonic[equilibrium_mask])
            eq_lj = np.nanmean(energy_lj[equilibrium_mask])
            eq_coulombic = np.nanmean(energy_coulombic[equilibrium_mask])
            eq_total = np.nanmean(energy_total[equilibrium_mask])
            
            print(f"    Equilibrium values: Harmonic = {eq_harmonic:.6f} Hartree, LJ = {eq_lj:.6f} Hartree")
            print(f"                        Coulombic = {eq_coulombic:.6f} Hartree, Total = {eq_total:.6f} Hartree")
            
            # Subtract equilibrium values (keep in Hartree)
            energy_harmonic_delta = energy_harmonic - eq_harmonic
            energy_lj_coulomb_delta = (energy_lj - eq_lj) + (energy_coulombic - eq_coulombic)  # Combine LJ+Coulomb
            energy_total_delta = energy_total - eq_total
        else:
            print(f"    Warning: No data points before 200 ps for equilibrium baseline")
            # Use raw values if no equilibrium data (keep in Hartree)
            energy_harmonic_delta = energy_harmonic
            energy_lj_coulomb_delta = energy_lj + energy_coulombic
            energy_total_delta = energy_total
        
        # Define energy components for plotting (now as deviations from equilibrium)
        # Using matplotlib classic colors: blue, red, green
        energy_components = [
            ('Harmonic', energy_harmonic_delta, 'blue', '-'),
            ('LJ+Coulomb', energy_lj_coulomb_delta, 'red', '-'),
            ('Total', energy_total_delta, 'green', '-', 3)  # Thicker line for total
        ]
        
        # Plot energy components
        for i, (name, values, color, style, *width) in enumerate(energy_components):
            linewidth = width[0] if width else 2
            
            # Remove NaN values
            valid_mask = ~np.isnan(values)
            if np.any(valid_mask):
                time_clean = time[valid_mask]
                values_clean = values[valid_mask]
                
                ax_energy.plot(time_clean, values_clean, color=color, linestyle=style, 
                             linewidth=linewidth, alpha=0.8, label=name)
                
                # Print energy statistics
                mean_val = np.mean(values_clean)
                std_val = np.std(values_clean)
                print(f"    {name} Energy: Mean = {mean_val:.6f} Hartree, Std = {std_val:.6f} Hartree, Points = {len(values_clean)}")
        
        # Add horizontal reference line at zero for equilibrium
        ax_energy.axhline(y=0, color='black', linestyle='--', alpha=0.7, linewidth=1, label='Equilibrium')
        
        # Format energy panel
        ax_energy.set_ylabel(latex_safe(r'$\Delta V$ (a.u.)', 'ΔV (a.u.)', use_latex), fontsize=14)
        ax_energy.grid(True, alpha=0.3)
        ax_energy.legend(fontsize=10, loc='best')
        # Remove title as requested
        
        # Add vertical line at time shift point for energy panel
        ax_energy.axvline(x=time_shift_ps, color='red', linestyle=':', alpha=0.7, linewidth=2)
    else:
        # Hide energy panel if no raw energy data available
        ax_energy.text(0.5, 0.5, 'Raw energy data not available\n(using pre-calculated fictive temperatures)', 
                      transform=ax_energy.transAxes, ha='center', va='center', 
                      fontsize=12, style='italic', color='gray')
        ax_energy.set_xlim(0, np.max(time) if len(time) > 0 else 1000)
        ax_energy.set_ylabel(latex_safe(r'$\Delta V$ (a.u.)', 'ΔV (a.u.)', use_latex), fontsize=14)
        ax_energy.grid(True, alpha=0.3)
        # Remove title as requested
    
    # Bottom panel: Fictive temperatures (existing code)
    # Using matplotlib classic colors: blue, red, green
    components = [
        ('Harmonic', fictive_harmonic, 'blue', '-'),
        ('LJ+Coulomb', fictive_lj, 'red', '-'),
        ('Total', fictive_total, 'green', '-', 3)  # Thicker line for total
    ]
    
    # Plot fictive temperature components
    for i, (name, values, color, style, *width) in enumerate(components):
        linewidth = width[0] if width else 2
        
        # Remove NaN values
        valid_mask = ~np.isnan(values)
        if np.any(valid_mask):
            time_clean = time[valid_mask]
            values_clean = values[valid_mask]
            
            ax_temp.plot(time_clean, values_clean, color=color, linestyle=style, 
                        linewidth=linewidth, alpha=0.8, label=name)
            
            # Print statistics
            mean_val = np.mean(values_clean)
            std_val = np.std(values_clean)
            print(f"    {name}: Mean = {mean_val:.1f} K, Std = {std_val:.1f} K, Points = {len(values_clean)}")
        else:
            print(f"    {name}: All NaN values")
    
    # Add horizontal reference line at 100 K
    ax_temp.axhline(y=100, color='black', linestyle='--', alpha=0.7, linewidth=1, label='100 K')
    
    # Format fictive temperature panel
    ax_temp.set_xlabel(latex_safe(r'$t$ (ps)', 't (ps)', use_latex), fontsize=14)
    ax_temp.set_ylabel(latex_safe(r'$T_{\mathrm{S}}$ (K)', 'T_S (K)', use_latex), fontsize=14)
    # Remove logarithmic scale as requested
    ax_temp.grid(True, alpha=0.3)
    # Remove legend as requested
    # Remove title as requested
    
    # Remove overall title as requested
    
    # Add vertical line at time shift point to fictive temperature panel
    ax_temp.axvline(x=time_shift_ps, color='red', linestyle=':', alpha=0.7, linewidth=2, 
                   label=f'Fit origin ({time_shift_ps} ps)')
    
    # Set reasonable limits for fictive temperature panel
    ax_temp.set_xlim(left=0)
    ax_temp.set_ylim(bottom=0)  # Set minimum to 0 K as requested
    
    # Create inset to zoom in on LJ+Coulomb temperature drop (200-750 ps)
    inset_ax = inset_axes(ax_temp, width="50%", height="50%", loc='upper right')
    
    # Plot only LJ+Coulomb and Total in the inset for clarity
    inset_components = [
        ('LJ+Coulomb', fictive_lj, 'red', '-'),
        ('Total', fictive_total, 'green', '-', 3)
    ]
    
    for name, values, color, style, *width in inset_components:
        linewidth = width[0] if width else 2
        
        # Remove NaN values and apply zoom range
        valid_mask = ~np.isnan(values)
        if np.any(valid_mask):
            time_clean = time[valid_mask]
            values_clean = values[valid_mask]
            
            # Filter for zoom range (200-750 ps)
            zoom_mask = (time_clean >= 200) & (time_clean <= 750)
            if np.any(zoom_mask):
                time_zoom = time_clean[zoom_mask]
                values_zoom = values_clean[zoom_mask]
                
                inset_ax.plot(time_zoom, values_zoom, color=color, linestyle=style, 
                             linewidth=linewidth, alpha=0.8, label=name)
    
    # Add reference lines to inset
    inset_ax.axhline(y=100, color='black', linestyle='--', alpha=0.7, linewidth=1)
    inset_ax.axvline(x=time_shift_ps, color='red', linestyle=':', alpha=0.7, linewidth=1)
    
    # Format inset
    inset_ax.set_xlim(200, 750)
    inset_ax.set_ylim(bottom=70, top=100)  # Set y-axis range from 70 K to 100 K
    inset_ax.grid(True, alpha=0.3)
    inset_ax.tick_params(axis='both', which='major', labelsize=12)  # Increased font size
    inset_ax.set_xlabel(latex_safe(r'$t$ (ps)', 't (ps)', use_latex), fontsize=12)  # Increased font size
    inset_ax.set_ylabel(latex_safe(r'$T_{\mathrm{S}}$ (K)', 'T_S (K)', use_latex), fontsize=12)  # Increased font size
    
    # Add a subtle border to the inset
    for spine in inset_ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_edgecolor('gray')
    
    # Add connection lines between inset and main panel region
    mark_inset(ax_temp, inset_ax, loc1=2, loc2=4, fc="none", ec="gray", linestyle='--', alpha=0.7)
    
    # Text box removed as requested
    
    # Skip tight_layout due to inset compatibility
    
    # Save plot
    safe_coupling_label = coupling_label.replace('\\times', 'x').replace('^{-', 'neg').replace('}', '').replace('.', 'p')
    
    # Add smoothing info to filename if smoothing is applied
    if smooth_window_ps > 0:
        smooth_suffix = f'_smooth{int(smooth_window_ps)}ps'
    else:
        smooth_suffix = ''
    
    output_file = Path(output_dir) / f'fictive_components_coupling_{safe_coupling_label}{smooth_suffix}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    
    output_file_pdf = Path(output_dir) / f'fictive_components_coupling_{safe_coupling_label}{smooth_suffix}.pdf'
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file_pdf}")
    
    plt.close()

def plot_all_coupling_fictive_components(data_dict, coupling_info, output_dir='.', smooth_window_ps=100.0, time_shift_ps=201.0, use_latex=True, coupling_filter=None):
    """Create plots for all coupling strengths showing energy components."""
    print(f"\nCreating fictive temperature component plots for each coupling...")
    if smooth_window_ps > 0:
        print(f"Using {smooth_window_ps} ps smoothing window")
    print(f"Using {time_shift_ps} ps as time origin for exponential fitting")
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info, coupling_filter)
    
    if not filtered_coupling_info:
        print("No target coupling data to plot")
        return
    
    # Sort couplings by value
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    plots_created = 0
    
    for coupling_name in sorted_couplings:
        if coupling_name not in data_dict:
            continue
            
        data = data_dict[coupling_name]
        coupling_label = data['coupling_label']
        
        print(f"\nProcessing coupling {coupling_label}:")
        plot_coupling_fictive_components(
            coupling_name, data, output_dir, smooth_window_ps, time_shift_ps, use_latex=use_latex
        )
        plots_created += 1
    
    print(f"\nCreated {plots_created} component plots")

def create_summary_comparison(data_dict, coupling_info, output_dir='.', smooth_window_ps=100.0, use_latex=True):
    """Create a summary plot comparing one component across all couplings."""
    if use_latex:
        use_latex = setup_latex_fonts()
    else:
        from latex_config_adobe import _apply_fallback_fonts
        _apply_fallback_fonts()
    
    print("\nCreating summary comparison plot...")
    if smooth_window_ps > 0:
        print(f"Using {smooth_window_ps} ps smoothing window for summary")
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    if not filtered_coupling_info:
        print("No target coupling data to plot")
        return
    
    # Sort couplings by value
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    # Create figure with subplots for different components
    fig, axes = plt.subplots(2, 2, figsize=(6, 8))
    axes = axes.flatten()
    
    components = [
        ('Harmonic', 'fictive_harmonic', 'blue'),
        ('LJ', 'fictive_lj', 'red'), 
        ('Total PE', 'fictive_total', 'green'),
        ('All Components', 'all', 'various')  # Special case for overlay
    ]
    
    for comp_idx, (comp_name, comp_key, comp_color) in enumerate(components):
        ax = axes[comp_idx]
        
        if comp_key == 'all':
            # Special case: overlay all components for zero coupling
            zero_coupling_data = None
            for coupling_name in sorted_couplings:
                if coupling_name in data_dict:
                    coupling_value = data_dict[coupling_name]['coupling_value']
                    if abs(coupling_value) < 1e-10:  # Zero coupling
                        zero_coupling_data = data_dict[coupling_name]
                        break
            
            if zero_coupling_data:
                time = zero_coupling_data['time']
                
                overlay_components = [
                    ('Harmonic', zero_coupling_data['fictive_harmonic'], 'blue', '-'),
                    ('LJ', zero_coupling_data['fictive_lj'], 'red', '-'),
                    ('Total PE', zero_coupling_data['fictive_total'], 'green', '-')
                ]
                
                for name, values, color, style in overlay_components:
                    # Apply smoothing if specified
                    if smooth_window_ps > 0:
                        values = smooth_data(time, values, smooth_window_ps)
                    
                    valid_mask = ~np.isnan(values)
                    if np.any(valid_mask):
                        time_clean = time[valid_mask]
                        values_clean = values[valid_mask]
                        ax.plot(time_clean, values_clean, color=color, linestyle=style,
                               linewidth=2, alpha=0.8, label=name)
                
                ax.set_title(latex_safe('All Components ($\\varepsilon = 0$)', 'All Components (ε = 0)', use_latex), fontsize=14, fontweight='bold')
                ax.legend(fontsize=10)
        else:
            # Plot specific component for all couplings
            for coupling_name in sorted_couplings:
                if coupling_name not in data_dict:
                    continue
                    
                data = data_dict[coupling_name]
                coupling_value = data['coupling_value']
                coupling_label = data['coupling_label']
                time = data['time']
                values = data[comp_key]
                
                # Apply smoothing if specified
                if smooth_window_ps > 0:
                    values = smooth_data(time, values, smooth_window_ps)
                
                # Set color based on coupling strength
                if coupling_value == 0.0:
                    color = '#1f77b4'  # Blue
                elif coupling_value <= 5e-4:
                    color = '#ff7f0e'  # Orange
                else:
                    color = '#d62728'  # Red
                
                # Remove NaN values
                valid_mask = ~np.isnan(values)
                if np.any(valid_mask):
                    time_clean = time[valid_mask]
                    values_clean = values[valid_mask]
                    
                    if coupling_value == 0.0:
                        label = latex_safe('$\\varepsilon = 0$', 'ε = 0', use_latex)
                    else:
                        label = latex_safe(f'$\\varepsilon = {coupling_label}$', f'ε = {coupling_label}', use_latex)
                    
                    ax.plot(time_clean, values_clean, color=color, linewidth=2, 
                           alpha=0.7, label=label)
            
            ax.set_title(f'{comp_name} Component', fontsize=14, fontweight='bold')
            ax.legend(fontsize=10)
        
        # Format subplot
        ax.set_xlabel('t (ps)', fontsize=12)
        ax.set_ylabel('Fictive Temperature (K)', fontsize=12)
        # Remove logarithmic scale as requested
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=75)  # Set minimum to 75 K as requested
    
    plt.suptitle('Fictive Temperature Components Comparison', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    # Save plot
    if smooth_window_ps > 0:
        smooth_suffix = f'_smooth{int(smooth_window_ps)}ps'
    else:
        smooth_suffix = ''
    
    output_file = Path(output_dir) / f'fictive_components_summary{smooth_suffix}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved summary: {output_file}")
    
    output_file_pdf = Path(output_dir) / f'fictive_components_summary{smooth_suffix}.pdf'
    plt.savefig(output_file_pdf, dpi=300, bbox_inches='tight')
    print(f"Saved summary: {output_file_pdf}")
    
    plt.close()

def analyze_component_statistics(data_dict, coupling_info):
    """Analyze and print statistics for each component and coupling."""
    print("\n=== Fictive Temperature Component Statistics ===")
    
    # Filter to specific coupling strengths
    filtered_coupling_info = filter_coupling_strengths(coupling_info)
    
    # Sort couplings by value
    sorted_couplings = sorted(filtered_coupling_info.keys(), 
                            key=lambda x: filtered_coupling_info[x][0])
    
    components = ['fictive_harmonic', 'fictive_lj', 'fictive_total']
    component_names = ['Harmonic', 'LJ', 'Total PE']
    
    for comp_idx, (comp_key, comp_name) in enumerate(zip(components, component_names)):
        print(f"\n--- {comp_name} Component ---")
        print(f"{'Coupling':<15} {'Mean (K)':<10} {'Std (K)':<10} {'Min (K)':<10} {'Max (K)':<10} {'Points':<8}")
        print("-" * 70)
        
        for coupling_name in sorted_couplings:
            if coupling_name not in data_dict:
                continue
                
            data = data_dict[coupling_name]
            coupling_label = data['coupling_label']
            values = data[comp_key]
            
            # Handle NaN values
            valid_mask = ~np.isnan(values)
            if np.any(valid_mask):
                values_clean = values[valid_mask]
                mean_val = np.mean(values_clean)
                std_val = np.std(values_clean)
                min_val = np.min(values_clean)
                max_val = np.max(values_clean)
                n_points = len(values_clean)
            else:
                mean_val = std_val = min_val = max_val = np.nan
                n_points = 0
            
            print(f"{coupling_label:<15} {mean_val:<10.1f} {std_val:<10.1f} {min_val:<10.1f} {max_val:<10.1f} {n_points:<8}")

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Analyze and plot fictive temperature components for each coupling using empirical U(T) relationships')
    parser.add_argument('--base_dir', default='.', 
                       help='Base directory containing coupling directories or time_series_output (default: .)')
    parser.add_argument('--output_dir', default='./', 
                       help='Output directory for plots (default: ./)')
    parser.add_argument('--smooth_window', type=float, default=100.0,
                       help='Smoothing window size in picoseconds (default: 100.0, use 0 to disable smoothing)')
    parser.add_argument('--time_shift', type=float, default=201.0,
                       help='Time origin shift for exponential fitting in picoseconds (default: 201.0)')
    parser.add_argument('--use_empirical', action='store_true', default=True,
                       help='Use empirical U(T) relationships for fictive temperature calculation (default: True)')
    parser.add_argument('--empirical_data_file', default='potential_energy_components_vs_temperature.txt',
                       help='Path to empirical energy-temperature calibration data (default: potential_energy_components_vs_temperature.txt)')
    parser.add_argument('--no-latex', action='store_true',
                       help='Disable LaTeX rendering (use when TeX fonts are unavailable)')
    parser.add_argument('--coupling', type=str, default=None,
                       help='Plot only this coupling strength (e.g. 3e-4 or 3eneg04)')
    
    args = parser.parse_args()
    use_latex = not args.no_latex
    coupling_filter = parse_coupling_arg(args.coupling) if args.coupling else None
    
    print("Fictive Temperature Components Analysis and Plotting Script")
    print("=" * 70)
    print(f"Base directory: {Path(args.base_dir).resolve()}")
    print(f"Output directory: {Path(args.output_dir).resolve()}")
    if args.smooth_window > 0:
        print(f"Smoothing window: {args.smooth_window} ps")
    else:
        print("Smoothing: Disabled")
    print(f"Time origin for fitting: {args.time_shift} ps")
    print(f"Empirical U(T) relationships: {'Enabled' if args.use_empirical else 'Disabled'}")
    if args.use_empirical:
        print(f"Empirical data file: {args.empirical_data_file}")
    
    # Collect all data
    print("\nCollecting fictive temperature component data...")
    data_dict, coupling_info = collect_fictive_temperature_data(args.base_dir)
    
    if not data_dict:
        print("No data found! Make sure you have coupling_*_averaged_fictive_temperatures.txt files.")
        return 1
    
    print(f"\nSummary:")
    print(f"  Found {len(data_dict)} coupling strengths")
    
    # Create output directory if it doesn't exist
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # Analyze statistics (on original data, before smoothing)
    analyze_component_statistics(data_dict, coupling_info)
    
    # Generate plots for each coupling
    plot_all_coupling_fictive_components(
        data_dict, coupling_info, args.output_dir, args.smooth_window, args.time_shift,
        use_latex=use_latex, coupling_filter=coupling_filter,
    )
    
    # Create summary comparison (skip when plotting a single coupling)
    if coupling_filter is None:
        create_summary_comparison(
            data_dict, coupling_info, args.output_dir, args.smooth_window, use_latex=use_latex
        )
    
    print("\n" + "=" * 70)
    print("Analysis complete! Generated plots:")
    print(f"  - fictive_components_coupling_*.png/pdf - Individual coupling component plots with material time info")
    print(f"  - fictive_components_summary.png/pdf - Summary comparison plot")
    print(f"\nMaterial time calculation details:")
    print(f"  - LJ+Coul component fit: T(t) = T∞ - A*exp(-t/τ) (extracts relaxation time τ)")
    print(f"  - Material time Ξ(t) = ∫(1/τ)dt calculated from relaxation time")
    print(f"  - Rescaled material time ξ(t) with condition: exp(-ξ^0.55) = 0.1 when ξ = 1")
    print(f"  - Time origin shifted to {args.time_shift} ps for fitting")
    
    return 0

if __name__ == "__main__":
    exit(main())
