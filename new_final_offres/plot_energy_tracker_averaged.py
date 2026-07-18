#!/usr/bin/env python3
"""
Plot Averaged Energy Tracker Output from Multiple Replicas

This script finds all energy tracker files (prod-*_energy_tracker.txt) in a directory,
averages them using interpolation onto a uniform time grid, and creates a two-panel plot:
1. Energy components vs time (molecular, cavity, reservoir, universe) with error bars
2. Temperature vs time with error bars

This script is based on the existing averaging procedure from average_energy_replicas.py
and the plotting functionality from plot_energy_tracker_simple.py.

Usage:
    python plot_energy_tracker_averaged.py <directory>
    
Example:
    python plot_energy_tracker_averaged.py cavity_coupling_1eneg03_switch_1.0ps/
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
import glob
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

def create_uniform_time_grid(time_arrays, num_points=None):
    """
    Create a uniform time grid based on the common range of time values.
    
    Parameters:
    -----------
    time_arrays : list of np.arrays
        List of time arrays from all replicas
    num_points : int, optional
        Number of grid points. If None, uses the median number of points from replicas
        
    Returns:
    --------
    uniform_grid : np.array
        Uniform grid of time values
    """
    # Find the common range (intersection of all ranges)
    min_time = max([time_array.min() for time_array in time_arrays])
    max_time = min([time_array.max() for time_array in time_arrays])
    
    if num_points is None:
        # Use median number of points from all replicas
        num_points = int(np.median([len(time_array) for time_array in time_arrays]))
    
    return np.linspace(min_time, max_time, num_points)

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

def average_energy_replicas(folder_path, grid_points=None):
    """
    Average energy replica files using memory-efficient online averaging.
    
    Parameters:
    -----------
    folder_path : str
        Path to the folder containing replica files
    grid_points : int, optional
        Number of grid points for interpolation
        
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
    
    # Sample files to determine uniform time grid
    print("Determining uniform time grid from all files...")
    time_arrays = []
    column_names = None
    
    for filepath in tqdm(replica_files, desc="Reading time grids", unit="file"):
        time, data, cols = read_energy_file(filepath)
        if time is not None:
            time_arrays.append(time)
            if column_names is None:
                column_names = cols
    
    if not time_arrays:
        raise ValueError("No valid files found!")
    
    uniform_grid = create_uniform_time_grid(time_arrays, grid_points)
    print(f"Created uniform time grid: {len(uniform_grid)} points from {uniform_grid.min():.6f} to {uniform_grid.max():.6f} ps")
    
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
                                        num_replicas, output_file=None):
    """Create a two-panel plot of averaged energy components and temperature."""
    
    # Calculate energy components
    (molecular_total, cavity_total, reservoir_total, universe_total,
     molecular_total_std, cavity_total_std, reservoir_total_std, universe_total_std) = calculate_energy_components_from_averaged(
        averaged_data, std_data, columns)
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
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
    
    ax2.set_xlabel('Time (ps)', fontsize=12)
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
    
    args = parser.parse_args()
    
    # Check if directory exists
    input_path = Path(args.directory)
    if not input_path.exists():
        print(f"Error: Directory '{args.directory}' not found")
        sys.exit(1)
    
    if not input_path.is_dir():
        print(f"Error: '{args.directory}' is not a directory")
        sys.exit(1)
    
    print(f"Processing energy tracker files in: {args.directory}")
    
    # Average data from multiple replicas
    uniform_grid, averaged_data, std_data, column_names, count_per_point = average_energy_replicas(
        args.directory, args.grid_points)
    
    if uniform_grid is None:
        print("Failed to process replica files!")
        sys.exit(1)
    
    # Create column index mapping
    columns = {name: i for i, name in enumerate(column_names)}
    
    num_replicas = np.max(count_per_point)
    print(f"Loaded averaged data from {num_replicas} replicas")
    print(f"Time range: {uniform_grid[0]:.6f} - {uniform_grid[-1]:.6f} ps")
    
    # Print statistics if requested
    if args.stats:
        print_averaged_energy_statistics(uniform_grid, averaged_data, std_data, columns, count_per_point, num_replicas)
    
    # Create plot
    print("Creating averaged plot with error bars...")
    fig = plot_averaged_energy_and_temperature(uniform_grid, averaged_data, std_data, columns, 
                                              count_per_point, num_replicas, args.output)
    
    if not args.output:
        print("Close the plot window to exit.")

if __name__ == '__main__':
    main() 