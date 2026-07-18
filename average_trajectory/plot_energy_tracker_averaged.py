#!/usr/bin/env python3
"""
Plot Averaged Energy Tracker Output from Multiple Replicas

This script finds all energy tracker files (prod-*_energy_tracker.txt) in a directory,
averages them using interpolation onto a uniform time grid, and creates a multi-panel plot:
1. Energy components vs time (molecular, cavity, reservoir, universe) with error bars
2. Temperature vs time with error bars  
3. Fictive temperatures vs time (optional, derived from energy components using U(T) relationships)

NEW FEATURES:
- Save averaged data to CSV or NPZ formats
- Load previously saved averaged data for quick replotting
- Comprehensive error propagation and statistics
- Fictive temperature calculation with extrapolation filtering
- Explicit grid specification for performance optimization

This script is based on the existing averaging procedure from average_energy_replicas.py
and the plotting functionality from plot_energy_tracker_simple.py.

Usage:
    python plot_energy_tracker_averaged.py <directory> [options]
    
Examples:
    # Basic usage (reads all files to determine time range)
    python plot_energy_tracker_averaged.py cavity_coupling_1eneg03_switch_1.0ps/
    
    # Fast preview with limited resolution
    python plot_energy_tracker_averaged.py cavity_coupling_1eneg03_switch_1.0ps/ --fast
    
    # EXPLICIT GRID (FASTEST - no file reading for grid determination)
    python plot_energy_tracker_averaged.py cavity_coupling_1eneg03_switch_1.0ps/ --time-min 0.0 --time-max 100.0 --grid-points 5000
    
    # Explicit grid with fixed time step (no file reading)
    python plot_energy_tracker_averaged.py cavity_coupling_1eneg03_switch_1.0ps/ --time-min 10.0 --time-max 50.0 --time-step 0.01
    
    # Custom time range with auto-detected resolution (some file reading)
    python plot_energy_tracker_averaged.py cavity_coupling_1eneg03_switch_1.0ps/ --time-min 10.0 --time-max 50.0
    
    # Save averaged data
    python plot_energy_tracker_averaged.py cavity_coupling_1eneg03_switch_1.0ps/ --save-data averaged_results --save-format both
    
    # Load and replot saved data
    python plot_energy_tracker_averaged.py . --load-data averaged_results.csv --output replot.png
    
    # Include fictive temperatures (interpolated values only)
    python plot_energy_tracker_averaged.py cavity_coupling_1eneg03_switch_1.0ps/ --fictive-temp
    
    # Include fictive temperatures with extrapolated values
    python plot_energy_tracker_averaged.py cavity_coupling_1eneg03_switch_1.0ps/ --fictive-temp --include-extrapolated
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
import glob
import pandas as pd
import time
import warnings
from pathlib import Path
from scipy.interpolate import interp1d
from tqdm import tqdm

def find_energy_replica_files(folder_path):
    """
    Find all energy tracker replica files in the specified folder.
    
    Parameters:
    -----------
    folder_path : str
        Path to the folder containing replica files
        
    Returns:
    --------
    replica_files : list
        List of file paths for energy tracker files
    """
    pattern = os.path.join(folder_path, "prod-*_energy_tracker.txt")
    return sorted(glob.glob(pattern))

def read_time_column_only(filepath):
    """
    Fast extraction of just the time column from an energy tracker file.
    Used for determining the uniform time grid without reading all data.
    
    Parameters:
    -----------
    filepath : str
        Path to the energy tracker file
        
    Returns:
    --------
    time : np.array or None
        Array of time values, or None if error
    """
    try:
        # Use numpy's optimized reading - much faster than manual parsing
        # Read only the first column (time), skip header comments
        time_values = np.loadtxt(filepath, usecols=0, comments='#')
        return time_values
            
    except Exception as e:
        # Fallback to manual reading if numpy fails
        try:
            time_values = []
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        try:
                            first_value = float(line.split()[0])
                            time_values.append(first_value)
                        except (ValueError, IndexError):
                            continue
            
            if time_values:
                return np.array(time_values)
            else:
                return None
                
        except Exception as e2:
            print(f"Error reading time from {filepath}: {e2}")
            return None

def read_energy_file(filepath):
    """
    Read an energy tracker file and extract time and all energy columns.
    
    Parameters:
    -----------
    filepath : str
        Path to the energy tracker file
        
    Returns:
    --------
    time : np.array
        Array of time values
    data : np.array
        2D array of all data columns (including time)
    column_names : list
        List of column names
    """
    try:
        # Read the file, skipping header lines that start with '#'
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Find the header line with column names
        header_line = None
        data_start = 0
        for i, line in enumerate(lines):
            if line.strip() and not line.startswith('#'):
                header_line = line.strip()
                data_start = i
                break
        
        if header_line is None:
            raise ValueError("Could not find header line")
        
        # Extract column names
        column_names = header_line.split()
        
        # Read the data
        data = []
        for line in lines[data_start + 1:]:
            if line.strip() and not line.startswith('#'):
                values = [float(x) for x in line.split()]
                data.append(values)
        
        data = np.array(data)
        
        # Extract time column (first column)
        time = data[:, 0]
        
        return time, data, column_names
        
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return None, None, None

def create_uniform_time_grid(time_arrays, num_points=None, time_min=None, time_max=None, 
                           time_step=None, fast_mode=False):
    """
    Create a uniform time grid based on the common range of time values or explicit specification.
    
    Parameters:
    -----------
    time_arrays : list of np.arrays
        List of time arrays from all replicas (can be empty if time_min/max specified)
    num_points : int, optional
        Number of grid points. If None, uses the median number of points from replicas
    time_min : float, optional
        Minimum time for grid. If None, uses intersection of all replica ranges
    time_max : float, optional
        Maximum time for grid. If None, uses intersection of all replica ranges
    time_step : float, optional
        Time step for grid. If provided, overrides num_points
    fast_mode : bool, optional
        If True, limit to maximum 1000 points for quick preview
        
    Returns:
    --------
    uniform_grid : np.array
        Uniform grid of time values
    """
    # Handle explicit grid specification
    if time_min is not None and time_max is not None:
        min_time = time_min
        max_time = time_max
        print(f"Using explicitly specified time range: {min_time} to {max_time} ps")
    elif time_arrays:
        # Find the natural common range (intersection of all ranges)
        natural_min_time = max([time_array.min() for time_array in time_arrays])
        natural_max_time = min([time_array.max() for time_array in time_arrays])
        
        # Use user-specified range or natural range
        min_time = time_min if time_min is not None else natural_min_time
        max_time = time_max if time_max is not None else natural_max_time
    else:
        raise ValueError("Must provide either time_arrays or explicit time_min/time_max")
    
    # Validate time range
    if min_time >= max_time:
        raise ValueError(f"Invalid time range: min_time ({min_time}) >= max_time ({max_time})")
    
    # Warn if user range extends beyond natural range (only if we have time_arrays)
    if time_arrays and time_min is not None and time_max is not None:
        natural_min_time = max([time_array.min() for time_array in time_arrays])
        natural_max_time = min([time_array.max() for time_array in time_arrays])
        if time_min < natural_min_time:
            print(f"Warning: Requested min_time ({time_min}) is before natural range ({natural_min_time})")
        if time_max > natural_max_time:
            print(f"Warning: Requested max_time ({time_max}) is after natural range ({natural_max_time})")
    
    # Determine number of points
    if time_step is not None:
        # Use time step to determine number of points
        num_points = int(np.ceil((max_time - min_time) / time_step)) + 1
        print(f"Using time step {time_step} ps → {num_points} grid points")
    elif num_points is None:
        if time_arrays:
            # Use median number of points from all replicas
            num_points = int(np.median([len(time_array) for time_array in time_arrays]))
            print(f"Auto-determined grid points: {num_points} (median of replica lengths)")
        else:
            # Default when no time arrays and no explicit grid points
            raise ValueError("Must specify either --grid-points or --time-step when using explicit time range")
    
    # Apply fast mode limit
    if fast_mode and num_points > 1000:
        num_points = 1000
        print(f"Fast mode: Limited to {num_points} grid points")
    
    # Create the grid
    uniform_grid = np.linspace(min_time, max_time, num_points)
    
    # Print grid information
    actual_time_step = (max_time - min_time) / (num_points - 1) if num_points > 1 else 0
    print(f"Created uniform time grid:")
    print(f"  Range: {min_time:.6f} to {max_time:.6f} ps")
    print(f"  Points: {num_points}")
    print(f"  Effective time step: {actual_time_step:.6f} ps")
    
    return uniform_grid

def interpolate_to_time_grid(time, data, uniform_grid):
    """
    Interpolate all data columns onto a uniform time grid.
    
    Parameters:
    -----------
    time : np.array
        Original time array
    data : np.array
        Original data array (2D: time_points x columns)
    uniform_grid : np.array
        Target uniform time grid
        
    Returns:
    --------
    interpolated_data : np.array
        Data interpolated onto the uniform grid (2D: grid_points x columns)
    valid_mask : np.array
        Boolean mask indicating which grid points have valid interpolated data
    """
    try:
        # Find the range of the current file
        min_time = time.min()
        max_time = time.max()
        
        # Create mask for grid points within the file's time range
        valid_mask = (uniform_grid >= min_time) & (uniform_grid <= max_time)
        
        # Initialize output array with NaN
        interpolated_data = np.full((len(uniform_grid), data.shape[1]), np.nan)
        
        if np.any(valid_mask):
            # Interpolate each column
            for col_idx in range(data.shape[1]):
                # Create interpolation function
                interp_func = interp1d(time, data[:, col_idx], kind='linear', 
                                      bounds_error=False, fill_value=np.nan, assume_sorted=True)
                
                # Interpolate only for valid grid points
                interpolated_data[valid_mask, col_idx] = interp_func(uniform_grid[valid_mask])
        
        return interpolated_data, valid_mask
        
    except Exception as e:
        print(f"Error during interpolation: {e}")
        return None, None

def update_online_stats_with_mask(mean, variance, count_per_point, new_data, valid_mask):
    """
    Update running mean and variance using Welford's online algorithm with variable sample counts.
    
    Parameters:
    -----------
    mean : np.array
        Current running mean (2D: time_points x columns)
    variance : np.array  
        Current running variance (2D: time_points x columns)
    count_per_point : np.array
        Current count of samples per time point
    new_data : np.array
        New data to incorporate (2D: time_points x columns, may contain NaN)
    valid_mask : np.array
        Boolean mask indicating which points have valid data
        
    Returns:
    --------
    updated_mean : np.array
        Updated running mean
    updated_variance : np.array
        Updated running variance
    updated_count : np.array
        Updated count per time point
    """
    for i in range(len(valid_mask)):
        if valid_mask[i]:
            # Only update count once per time point (not per column)
            count_per_point[i] += 1
            n = count_per_point[i]
            
            # Update each column independently
            for j in range(new_data.shape[1]):
                if not np.isnan(new_data[i, j]):
                    # Calculate delta
                    delta = new_data[i, j] - mean[i, j]
                    
                    # Update mean
                    mean[i, j] += delta / n
                    
                    # Update variance
                    if n > 1:
                        delta2 = new_data[i, j] - mean[i, j]
                        variance[i, j] += delta * delta2
    
    return mean, variance, count_per_point

def average_energy_replicas(folder_path, grid_points=None, time_min=None, time_max=None, 
                          time_step=None, fast_mode=False):
    """
    Average energy replica files using memory-efficient online averaging.
    
    Parameters:
    -----------
    folder_path : str
        Path to the folder containing replica files
    grid_points : int, optional
        Number of grid points for interpolation
    time_min : float, optional
        Minimum time for interpolation grid (ps)
    time_max : float, optional
        Maximum time for interpolation grid (ps)
    time_step : float, optional
        Time step for uniform grid (ps)
    fast_mode : bool, optional
        Fast mode with limited grid points
        
    Returns:
    --------
    uniform_grid : np.array
        Time grid
    averaged_data : np.array
        Averaged data
    std_data : np.array
        Standard deviation data
    column_names : list
        Column names
    count_per_point : np.array
        Sample count per time point
    """
    # Find all energy replica files
    replica_files = find_energy_replica_files(folder_path)
    
    if not replica_files:
        print(f"No files matching pattern 'prod-*_energy_tracker.txt' found in {folder_path}")
        return None, None, None, None, None
    
    print(f"\nProcessing energy replica files in {folder_path}:")
    print(f"Found {len(replica_files)} replica files:")
    for f in sorted(replica_files):
        print(f"  - {os.path.basename(f)}")
    
    # Determine uniform time grid
    if time_min is not None and time_max is not None and (grid_points is not None or time_step is not None):
        # Grid fully specified - no need to read files for grid determination
        print("Using explicitly specified time grid...")
        uniform_grid = create_uniform_time_grid([], grid_points, time_min, time_max, 
                                              time_step, fast_mode)
        
        # Read one file completely to get column names
        print("Reading column names from first file...")
        column_names = None
        for filepath in replica_files:
            _, _, cols = read_energy_file(filepath)
            if cols is not None:
                column_names = cols
                break
                
        if column_names is None:
            raise ValueError("Could not read column names from any file!")
            
    else:
        # Need to sample files to determine uniform time grid
        print("Determining uniform time grid from all files...")
        time_arrays = []
        column_names = None
        
        # First pass: quickly read only time columns for grid determination
        for filepath in tqdm(replica_files, desc="Reading time grids", unit="file"):
            time = read_time_column_only(filepath)
            if time is not None:
                time_arrays.append(time)
        
        # Read one file completely to get column names
        if time_arrays:
            print("Reading column names from first file...")
            for filepath in replica_files:
                _, _, cols = read_energy_file(filepath)
                if cols is not None:
                    column_names = cols
                    break
        
        if not time_arrays:
            raise ValueError("No valid files found!")
        
        uniform_grid = create_uniform_time_grid(time_arrays, grid_points, time_min, time_max, 
                                              time_step, fast_mode)
    
    # Initialize arrays for online averaging
    n_grid_points = len(uniform_grid)
    n_columns = len(column_names)
    running_mean = np.zeros((n_grid_points, n_columns))
    running_variance = np.zeros((n_grid_points, n_columns))
    count_per_point = np.zeros(n_grid_points, dtype=int)  # Track samples per time point
    
    # Process replicas one by one for memory efficiency
    print(f"\nProcessing all energy replicas...")
    successful_replicas = 0
    
    for filepath in tqdm(replica_files, desc="Processing energy replicas", unit="file"):
        # Read current replica
        time, data, cols = read_energy_file(filepath)
        
        if time is not None and data is not None:
            # Interpolate current replica
            interp_data, valid_mask = interpolate_to_time_grid(time, data, uniform_grid)
            
            if interp_data is not None and valid_mask is not None:
                # Check if we have any valid data points
                if np.any(valid_mask):
                    # Update running statistics
                    running_mean, running_variance, count_per_point = update_online_stats_with_mask(
                        running_mean, running_variance, count_per_point, interp_data, valid_mask
                    )
                    successful_replicas += 1
                else:
                    print(f"No valid data points for file: {filepath}")
            else:
                print(f"Failed to interpolate file: {filepath}")
        else:
            print(f"Failed to read file: {filepath}")
    
    if successful_replicas == 0:
        print("No successful interpolations!")
        return None, None, None, None, None
    
    # Calculate final standard deviation with variable sample counts
    std_data = np.zeros_like(running_mean)
    for i in range(n_grid_points):
        for j in range(n_columns):
            if count_per_point[i] > 1:
                std_data[i, j] = np.sqrt(running_variance[i, j] / (count_per_point[i] - 1))
            else:
                std_data[i, j] = 0.0
    
    print(f"Successfully averaged {successful_replicas} replicas")
    print(f"Sample count per time point: min={count_per_point.min()}, max={count_per_point.max()}, median={np.median(count_per_point):.0f}")
    
    return uniform_grid, running_mean, std_data, column_names, count_per_point

def calculate_energy_components_from_averaged(data, std_data, columns):
    """Calculate the energy components we want to plot from averaged data."""
    
    # Total molecular energy = molecular kinetic + molecular potential
    # Molecular potential = harmonic + lj + ewald_short + ewald_long
    molecular_potential = (data[:, columns['harmonic_energy']] + 
                          data[:, columns['lj_energy']] + 
                          data[:, columns['ewald_short_energy']] + 
                          data[:, columns['ewald_long_energy']])
    molecular_total = data[:, columns['molecular_kinetic_energy']] + molecular_potential
    
    # Calculate standard deviations for molecular potential
    molecular_potential_std = np.sqrt(
        std_data[:, columns['harmonic_energy']]**2 + 
        std_data[:, columns['lj_energy']]**2 + 
        std_data[:, columns['ewald_short_energy']]**2 + 
        std_data[:, columns['ewald_long_energy']]**2
    )
    molecular_total_std = np.sqrt(
        std_data[:, columns['molecular_kinetic_energy']]**2 + 
        molecular_potential_std**2
    )
    
    # Total cavity energy = cavity kinetic + cavity potential
    cavity_total = (data[:, columns['cavity_kinetic_energy']] + 
                   data[:, columns['cavity_total_potential_energy']])
    cavity_total_std = np.sqrt(
        std_data[:, columns['cavity_kinetic_energy']]**2 + 
        std_data[:, columns['cavity_total_potential_energy']]**2
    )
    
    # Total reservoir energy (already calculated in file)
    reservoir_total = data[:, columns['total_reservoir_energy']]
    reservoir_total_std = std_data[:, columns['total_reservoir_energy']]
    
    # Universe energy (already calculated in file - should be conserved)
    universe_total = data[:, columns['universe_total_energy']]
    universe_total_std = std_data[:, columns['universe_total_energy']]
    
    return (molecular_total, cavity_total, reservoir_total, universe_total,
            molecular_total_std, cavity_total_std, reservoir_total_std, universe_total_std)

def plot_averaged_energy_and_temperature(uniform_grid, averaged_data, std_data, columns, count_per_point, 
                                        num_replicas, fictive_temps=None, fictive_stds=None, output_file=None):
    """Create a multi-panel plot of averaged energy components, temperature, and fictive temperatures."""
    
    # Calculate energy components
    (molecular_total, cavity_total, reservoir_total, universe_total,
     molecular_total_std, cavity_total_std, reservoir_total_std, universe_total_std) = calculate_energy_components_from_averaged(
        averaged_data, std_data, columns)
    
    # Determine number of panels
    has_fictive = fictive_temps is not None and len(fictive_temps) > 0
    n_panels = 3 if has_fictive else 2
    
    # Create figure with subplots
    fig, axes = plt.subplots(n_panels, 1, figsize=(12, 5*n_panels), sharex=True)
    if n_panels == 1:
        axes = [axes]
    
    ax1, ax2 = axes[0], axes[1]
    if has_fictive:
        ax3 = axes[2]
    
    # Panel 1: Energy components without error bars
    time_ps = uniform_grid
    
    ax1.plot(time_ps, molecular_total, label='Total Molecular Energy', linewidth=2, color='blue')
    ax1.plot(time_ps, cavity_total, label='Total Cavity Energy', linewidth=2, color='red')
    ax1.plot(time_ps, reservoir_total, label='Total Reservoir Energy', linewidth=2, color='green')
    ax1.plot(time_ps, universe_total, label='Universe Energy (Conserved)', linewidth=2, color='black', linestyle='--')
    
    ax1.set_ylabel('Energy (Hartree)', fontsize=12)
    ax1.set_title(f'Averaged Energy Components vs Time ({num_replicas} replicas)', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Add text box with energy conservation info
    energy_drift = universe_total[-1] - universe_total[0]
    energy_std_avg = np.mean(universe_total_std)
    textstr = f'Universe Energy Drift: {energy_drift:.2e} Hartree\nAvg Std Dev: {energy_std_avg:.2e} Hartree\nReplicas: {num_replicas}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=9,
             verticalalignment='top', bbox=props)
    
    # Panel 2: Temperature without error bars
    temperature = averaged_data[:, columns['temperature']]
    temperature_std = std_data[:, columns['temperature']]
    
    ax2.plot(time_ps, temperature, label='System Temperature', linewidth=2, color='orange')
    
    ax2.set_ylabel('Temperature (K)', fontsize=12)
    ax2.set_title(f'Averaged Temperature vs Time ({num_replicas} replicas)', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # Add horizontal line at target temperature
    temp_mean = np.mean(temperature)
    temp_std_avg = np.mean(temperature_std)
    ax2.axhline(y=temp_mean, color='red', linestyle=':', alpha=0.7, 
                label=f'Mean: {temp_mean:.1f} ± {temp_std_avg:.1f} K')
    ax2.legend(fontsize=10)
    
    # Panel 3: Fictive temperatures (if available)
    if has_fictive:
        colors = {'total_PE': 'purple', 'harmonic': 'brown', 'lj': 'pink', 'coulombic': 'cyan'}
        component_names = {'total_PE': 'Total PE', 'harmonic': 'Harmonic', 'lj': 'Lennard-Jones', 'coulombic': 'Coulombic'}
        
        for comp_name, fictive_T in fictive_temps.items():
            if comp_name in fictive_stds:
                fictive_T_std = fictive_stds[comp_name]
                
                # Only plot non-NaN values
                valid_mask = ~np.isnan(fictive_T)
                if np.any(valid_mask):
                    color = colors.get(comp_name, 'gray')
                    display_name = component_names.get(comp_name, comp_name)
                    
                    ax3.plot(time_ps[valid_mask], fictive_T[valid_mask], 
                            label=f'{display_name} Fictive T', 
                            linewidth=2, color=color, alpha=0.8)
        
        # Add actual temperature for comparison
        ax3.plot(time_ps, temperature, label='Actual Temperature', 
                linewidth=2, color='black', linestyle='--', alpha=0.7)
        
        ax3.set_ylabel('Fictive Temperature (K)', fontsize=12)
        ax3.set_title(f'Fictive Temperatures vs Time ({num_replicas} replicas)\n(Interpolated values only - extrapolated excluded)', 
                     fontsize=14, fontweight='bold')
        ax3.legend(fontsize=10)
        ax3.grid(True, alpha=0.3)
        
        # Add info about fictive temperatures
        if fictive_temps:
            n_components = len([comp for comp, temps in fictive_temps.items() if np.any(~np.isnan(temps))])
            textstr = f'Fictive T Components: {n_components}\nExtrapolated values excluded'
            props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
            ax3.text(0.02, 0.98, textstr, transform=ax3.transAxes, fontsize=9,
                     verticalalignment='top', bbox=props)
    
    # Set x-label on the bottom panel
    if has_fictive:
        ax3.set_xlabel('Time (ps)', fontsize=12)
    else:
        ax2.set_xlabel('Time (ps)', fontsize=12)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save or show plot
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()
    
    return fig

def print_averaged_energy_statistics(uniform_grid, averaged_data, std_data, columns, count_per_point, num_replicas):
    """Print useful statistics about the averaged energy data."""
    
    (molecular_total, cavity_total, reservoir_total, universe_total,
     molecular_total_std, cavity_total_std, reservoir_total_std, universe_total_std) = calculate_energy_components_from_averaged(
        averaged_data, std_data, columns)
    
    temperature = averaged_data[:, columns['temperature']]
    temperature_std = std_data[:, columns['temperature']]
    
    print("\n" + "="*60)
    print(f"AVERAGED ENERGY STATISTICS ({num_replicas} replicas)")
    print("="*60)
    
    print(f"Simulation time: {uniform_grid[0]:.6f} - {uniform_grid[-1]:.6f} ps")
    print(f"Number of time points: {len(averaged_data)}")
    print(f"Sample count per point: min={count_per_point.min()}, max={count_per_point.max()}, median={np.median(count_per_point):.0f}")
    
    print(f"\nEnergy Components (Hartree) - Mean ± Std:")
    print(f"  Molecular Energy:")
    print(f"    Initial: {molecular_total[0]:.6f} ± {molecular_total_std[0]:.6f}")
    print(f"    Final:   {molecular_total[-1]:.6f} ± {molecular_total_std[-1]:.6f}")
    print(f"    Change:  {molecular_total[-1] - molecular_total[0]:.6f}")
    
    print(f"  Cavity Energy:")
    print(f"    Initial: {cavity_total[0]:.6f} ± {cavity_total_std[0]:.6f}")
    print(f"    Final:   {cavity_total[-1]:.6f} ± {cavity_total_std[-1]:.6f}")
    print(f"    Change:  {cavity_total[-1] - cavity_total[0]:.6f}")
    
    print(f"  Reservoir Energy:")
    print(f"    Initial: {reservoir_total[0]:.6f} ± {reservoir_total_std[0]:.6f}")
    print(f"    Final:   {reservoir_total[-1]:.6f} ± {reservoir_total_std[-1]:.6f}")
    print(f"    Change:  {reservoir_total[-1] - reservoir_total[0]:.6f}")
    
    print(f"  Universe Energy (should be conserved):")
    print(f"    Initial: {universe_total[0]:.6f} ± {universe_total_std[0]:.6f}")
    print(f"    Final:   {universe_total[-1]:.6f} ± {universe_total_std[-1]:.6f}")
    print(f"    Drift:   {universe_total[-1] - universe_total[0]:.2e}")
    print(f"    Avg Std: {np.mean(universe_total_std):.2e}")
    
    print(f"\nTemperature Statistics:")
    print(f"  Mean: {np.mean(temperature):.2f} ± {np.mean(temperature_std):.2f} K")
    print(f"  Range: {np.min(temperature):.2f} - {np.max(temperature):.2f} K")
    
    print("="*60)

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
        print("Fictive temperature calculation will be skipped.")
        return None, None

def create_energy_temperature_interpolators(temp_data, temp_columns):
    """Create interpolating functions U(T) for energy components."""
    if temp_data is None or temp_columns is None:
        return None
        
    interpolators = {}
    temperatures = temp_data[:, temp_columns['temperature']]
    
    # Create interpolators for different energy components
    # Map display names to actual column names in temperature components file
    energy_components = {
        'total_PE': 'total_PE_hartree',
        'harmonic': 'harmonic_hartree', 
        'lj': 'lj_hartree',
        'coulombic': 'coulombic_hartree'
    }
    
    print(f"Available temperature component columns: {list(temp_columns.keys())}")
    
    for display_name, col_name in energy_components.items():
        if col_name in temp_columns:
            print(f"Creating interpolator for {display_name} using column {col_name}")
            energies = temp_data[:, temp_columns[col_name]]
            
            # Sort by temperature for proper interpolation
            sort_idx = np.argsort(temperatures)
            T_sorted = temperatures[sort_idx]
            U_sorted = energies[sort_idx]
            
            # Create U(T) and T(U) interpolators
            U_of_T = interp1d(T_sorted, U_sorted, kind='linear', bounds_error=False, fill_value='extrapolate')
            T_of_U = interp1d(U_sorted, T_sorted, kind='linear', bounds_error=False, fill_value='extrapolate')
            
            # Create linear extrapolation functions for extreme cases
            # Use first and last few points for linear fits
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
                'low_linear': low_linear,
                'high_linear': high_linear,
                'energy_column': col_name
            }
        else:
            print(f"Skipping {display_name}: column {col_name} not found in temperature components")
    
    if not interpolators:
        print("Warning: No interpolators could be created!")
        return None
    
    print(f"Successfully created {len(interpolators)} interpolators: {list(interpolators.keys())}")
    return interpolators

def calculate_fictive_temperatures_averaged(averaged_data, std_data, columns, interpolators, 
                                          exclude_extrapolated=True):
    """Calculate fictive temperatures from averaged energy components."""
    if interpolators is None:
        return None, None
        
    fictive_temps = {}
    fictive_stds = {}
    
    for comp_name, interp_data in interpolators.items():
        # Map the temperature components column name to energy tracker column name
        temp_col = interp_data['energy_column']  # e.g., 'total_PE_hartree'
        
        # Handle different energy components
        if temp_col == 'total_PE_hartree':
            # Try to find total potential energy column
            if 'total_potential_energy' in columns:
                energies = averaged_data[:, columns['total_potential_energy']]
                energy_stds = std_data[:, columns['total_potential_energy']]
            else:
                # Calculate from molecular components
                mol_pe_cols = ['harmonic_energy', 'lj_energy', 'ewald_short_energy', 'ewald_long_energy']
                available_cols = [col for col in mol_pe_cols if col in columns]
                if available_cols:
                    energies = np.sum([averaged_data[:, columns[col]] for col in available_cols], axis=0)
                    energy_stds = np.sqrt(np.sum([std_data[:, columns[col]]**2 for col in available_cols], axis=0))
                else:
                    print(f"    Skipping {comp_name}: No suitable potential energy columns found")
                    continue
                    
        elif temp_col == 'harmonic_hartree' and 'harmonic_energy' in columns:
            energies = averaged_data[:, columns['harmonic_energy']]
            energy_stds = std_data[:, columns['harmonic_energy']]
            
        elif temp_col == 'lj_hartree' and 'lj_energy' in columns:
            energies = averaged_data[:, columns['lj_energy']]
            energy_stds = std_data[:, columns['lj_energy']]
            
        elif temp_col == 'coulombic_hartree':
            # Combine ewald components for coulombic energy
            coulomb_cols = ['ewald_short_energy', 'ewald_long_energy']
            available_coulomb = [col for col in coulomb_cols if col in columns]
            if available_coulomb:
                energies = np.sum([averaged_data[:, columns[col]] for col in available_coulomb], axis=0)
                energy_stds = np.sqrt(np.sum([std_data[:, columns[col]]**2 for col in available_coulomb], axis=0))
            else:
                print(f"    Skipping {comp_name}: No coulombic energy columns found")
                continue
        else:
            print(f"    Skipping {comp_name}: Unknown energy component {temp_col}")
            continue
            
        # Now we have energies and energy_stds for this component
        T_of_U = interp_data['T_of_U']
        U_range = interp_data['U_range']
        
        # Calculate fictive temperatures
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fictive_T = T_of_U(energies)
        
        # Only mask physically unreasonable temperatures
        reasonable_mask = (fictive_T > 0) & (fictive_T < 10000)  # 0-10000 K seems reasonable
        
        if exclude_extrapolated:
            # Only keep interpolated values (within original energy range)
            interpolated_mask = (energies >= U_range[0]) & (energies <= U_range[1])
            valid_mask = reasonable_mask & interpolated_mask
        else:
            valid_mask = reasonable_mask
        
        # Set invalid values to NaN
        fictive_T_clean = fictive_T.copy()
        fictive_T_clean[~valid_mask] = np.nan
        
        # Propagate uncertainty (approximate using numerical derivative)
        # dT/dU ≈ ΔT/ΔU at each point
        delta_U = 0.001 * np.abs(energies)  # Small energy perturbation
        delta_U[delta_U == 0] = 1e-6  # Avoid zero division
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            T_plus = T_of_U(energies + delta_U)
            T_minus = T_of_U(energies - delta_U)
            dT_dU = (T_plus - T_minus) / (2 * delta_U)
        
        # Propagate uncertainty: σ_T = |dT/dU| * σ_U
        fictive_T_std = np.abs(dT_dU) * energy_stds
        fictive_T_std[~valid_mask] = np.nan
        
        fictive_temps[comp_name] = fictive_T_clean
        fictive_stds[comp_name] = fictive_T_std
        
        # Print statistics
        n_valid = np.sum(valid_mask)
        n_total = len(fictive_T)
        if exclude_extrapolated:
            n_interpolated = np.sum((energies >= U_range[0]) & (energies <= U_range[1]) & reasonable_mask)
            print(f"  {comp_name}: {n_valid}/{n_total} valid fictive temperatures (interpolated only)")
        else:
            print(f"  {comp_name}: {n_valid}/{n_total} valid fictive temperatures (including extrapolated)")
    
    return fictive_temps, fictive_stds

def save_averaged_data(uniform_grid, averaged_data, std_data, column_names, count_per_point, 
                      output_path, format='csv'):
    """
    Save averaged energy data to file in various formats.
    
    Parameters:
    -----------
    uniform_grid : np.array
        Time grid
    averaged_data : np.array
        Averaged data (2D: time_points x columns)
    std_data : np.array
        Standard deviation data (2D: time_points x columns)
    column_names : list
        List of column names
    count_per_point : np.array
        Sample count per time point
    output_path : str
        Output file path (extension will be added if needed)
    format : str
        Output format: 'csv', 'npz', or 'both'
    
    Returns:
    --------
    saved_files : list
        List of saved file paths
    """
    output_path = Path(output_path)
    saved_files = []
    
    if format in ['csv', 'both']:
        # Save as CSV with comprehensive column headers
        csv_path = output_path.with_suffix('.csv')
        
        # Create column headers for averaged data
        avg_columns = [f"{col}_avg" for col in column_names]
        std_columns = [f"{col}_std" for col in column_names]
        
        # Combine all data
        all_data = np.column_stack([
            uniform_grid,  # Time column
            averaged_data,  # Averaged values
            std_data,      # Standard deviations
            count_per_point  # Sample counts
        ])
        
        # Create comprehensive header
        all_column_names = ['time'] + avg_columns + std_columns + ['sample_count']
        
        # Create DataFrame and save
        df = pd.DataFrame(all_data, columns=all_column_names)
        df.to_csv(csv_path, index=False, float_format='%.8e')
        saved_files.append(str(csv_path))
        print(f"Averaged data saved to CSV: {csv_path}")
    
    if format in ['npz', 'both']:
        # Save as NumPy compressed archive
        npz_path = output_path.with_suffix('.npz')
        
        np.savez_compressed(
            npz_path,
            time=uniform_grid,
            averaged_data=averaged_data,
            std_data=std_data,
            column_names=np.array(column_names),
            count_per_point=count_per_point,
            # Metadata
            num_replicas=np.max(count_per_point),
            time_range=[uniform_grid[0], uniform_grid[-1]],
            num_time_points=len(uniform_grid),
            num_columns=len(column_names)
        )
        saved_files.append(str(npz_path))
        print(f"Averaged data saved to NPZ: {npz_path}")
    
    return saved_files

def load_averaged_data(file_path):
    """
    Load previously saved averaged data.
    
    Parameters:
    -----------
    file_path : str
        Path to saved data file (.csv or .npz)
        
    Returns:
    --------
    uniform_grid : np.array
        Time grid
    averaged_data : np.array
        Averaged data
    std_data : np.array
        Standard deviation data
    column_names : list
        Column names
    count_per_point : np.array
        Sample count per time point
    """
    file_path = Path(file_path)
    
    if file_path.suffix == '.npz':
        # Load from NPZ
        data = np.load(file_path)
        return (data['time'], data['averaged_data'], data['std_data'], 
                data['column_names'].tolist(), data['count_per_point'])
    
    elif file_path.suffix == '.csv':
        # Load from CSV
        df = pd.read_csv(file_path)
        
        # Extract time
        uniform_grid = df['time'].values
        
        # Extract column names (remove _avg suffix)
        avg_columns = [col for col in df.columns if col.endswith('_avg')]
        column_names = [col[:-4] for col in avg_columns]  # Remove '_avg'
        
        # Extract averaged data
        averaged_data = df[avg_columns].values
        
        # Extract standard deviations
        std_columns = [f"{col}_std" for col in column_names]
        std_data = df[std_columns].values
        
        # Extract sample counts
        count_per_point = df['sample_count'].values.astype(int)
        
        return uniform_grid, averaged_data, std_data, column_names, count_per_point
    
    else:
        raise ValueError(f"Unsupported file format: {file_path.suffix}")

def main():
    parser = argparse.ArgumentParser(
        description='Plot averaged energy tracker output from multiple replicas with error bars',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('directory', type=str, 
                       help='Directory containing energy tracker files (prod-*_energy_tracker.txt)')
    parser.add_argument('--output', '-o', type=str, 
                       help='Output plot file (default: show plot)')
    parser.add_argument('--stats', action='store_true', 
                       help='Print detailed energy statistics')
    parser.add_argument('--grid-points', type=int, default=None,
                       help='Number of grid points for interpolation (default: median of replica lengths)')
    parser.add_argument('--time-min', type=float, default=None,
                       help='Minimum time for interpolation grid (ps). Default: auto-detect from replicas')
    parser.add_argument('--time-max', type=float, default=None,
                       help='Maximum time for interpolation grid (ps). Default: auto-detect from replicas')
    parser.add_argument('--time-step', type=float, default=None,
                       help='Time step for uniform grid (ps). Alternative to --grid-points')
    parser.add_argument('--save-data', type=str, default=None,
                       help='Save averaged data to file (specify base name, extensions added automatically)')
    parser.add_argument('--save-format', choices=['csv', 'npz', 'both'], default='csv',
                       help='Format for saving averaged data (default: csv)')
    parser.add_argument('--load-data', type=str, default=None,
                       help='Load previously saved averaged data instead of processing replicas')
    parser.add_argument('--fast', action='store_true',
                       help='Fast mode: use fewer grid points for quick preview (max 1000 points)')
    parser.add_argument('--fictive-temp', action='store_true',
                       help='Include fictive temperature panel (requires temperature components file)')
    parser.add_argument('--temp-components', type=str, default='potential_energy_components_vs_temperature.txt',
                       help='Temperature components file for fictive temperature calculation')
    parser.add_argument('--include-extrapolated', action='store_true',
                       help='Include extrapolated fictive temperature values (default: exclude)')
    
    args = parser.parse_args()
    
    # Check if we should load existing data or process replicas
    if args.load_data:
        print(f"Loading previously saved averaged data from: {args.load_data}")
        try:
            uniform_grid, averaged_data, std_data, column_names, count_per_point = load_averaged_data(args.load_data)
            print(f"Successfully loaded averaged data with {len(uniform_grid)} time points")
        except Exception as e:
            print(f"Error loading data: {e}")
            sys.exit(1)
    else:
        # Check if directory exists
        input_path = Path(args.directory)
        if not input_path.exists():
            print(f"Error: Directory '{args.directory}' not found")
            sys.exit(1)
        
        if not input_path.is_dir():
            print(f"Error: '{args.directory}' is not a directory")
            sys.exit(1)
        
        print(f"Processing energy tracker files in: {args.directory}")
        
        # Average data from multiple replicas with timing
        start_time = time.time()
        uniform_grid, averaged_data, std_data, column_names, count_per_point = average_energy_replicas(
            args.directory, args.grid_points, args.time_min, args.time_max, 
            args.time_step, args.fast)
        processing_time = time.time() - start_time
        
        print(f"Processing completed in {processing_time:.2f} seconds")
        
        if uniform_grid is None:
            print("Failed to process replica files!")
            sys.exit(1)
    
    # Create column index mapping
    columns = {name: i for i, name in enumerate(column_names)}
    
    num_replicas = np.max(count_per_point)
    print(f"Loaded averaged data from {num_replicas} replicas")
    print(f"Time range: {uniform_grid[0]:.6f} - {uniform_grid[-1]:.6f} ps")
    
    # Calculate fictive temperatures if requested
    fictive_temps, fictive_stds = None, None
    if args.fictive_temp:
        print("Loading temperature components for fictive temperature calculation...")
        temp_data, temp_columns = load_temperature_energy_components(args.temp_components)
        
        if temp_data is not None:
            print("Creating energy-temperature interpolators...")
            interpolators = create_energy_temperature_interpolators(temp_data, temp_columns)
            
            if interpolators:
                print("Calculating fictive temperatures...")
                fictive_temps, fictive_stds = calculate_fictive_temperatures_averaged(
                    averaged_data, std_data, columns, interpolators, 
                    exclude_extrapolated=not args.include_extrapolated)
                
                if fictive_temps:
                    print(f"Successfully calculated fictive temperatures for {len(fictive_temps)} components")
                else:
                    print("No fictive temperatures could be calculated")
            else:
                print("Could not create interpolators - fictive temperature calculation skipped")
        else:
            print("Could not load temperature components - fictive temperature calculation skipped")
    
    # Print statistics if requested
    if args.stats:
        print_averaged_energy_statistics(uniform_grid, averaged_data, std_data, columns, count_per_point, num_replicas)
    
    # Save averaged data if requested
    if args.save_data:
        print("Saving averaged data...")
        saved_files = save_averaged_data(uniform_grid, averaged_data, std_data, column_names, 
                                       count_per_point, args.save_data, args.save_format)
        print(f"Data saved to: {', '.join(saved_files)}")
    
    # Create plot
    plot_description = "Creating averaged plot"
    if fictive_temps:
        plot_description += " with fictive temperatures"
    print(f"{plot_description}...")
    fig = plot_averaged_energy_and_temperature(uniform_grid, averaged_data, std_data, columns, 
                                              count_per_point, num_replicas, fictive_temps, fictive_stds, args.output)
    
    if not args.output:
        print("Close the plot window to exit.")

if __name__ == '__main__':
    main() 