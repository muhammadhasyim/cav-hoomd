#!/usr/bin/env python3
"""
Script to average multiple energy replica files (prod-*_energy_tracker.txt) by interpolating onto a uniform time grid.

Usage: python average_energy_replicas.py <folder_path>

The script will:
1. Find all files matching pattern prod-*_energy_tracker.txt in the given folder
2. Extract time and all energy/property columns from each file
3. Create a uniform time grid based on the common time range
4. Interpolate each replica's data onto this grid
5. Average the interpolated data across all replicas using online averaging
6. Save the results as prod_energy_tracker_averaged.txt in the same folder
"""

import argparse
import os
import glob
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import sys
from tqdm import tqdm
import re


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


def determine_time_grid_from_sample(replica_files, sample_size=50, grid_points=None):
    """
    Determine uniform time grid by sampling a subset of files to save memory.
    
    Parameters:
    -----------
    replica_files : list
        List of all replica file paths
    sample_size : int
        Number of files to sample for grid determination
    grid_points : int, optional
        Number of grid points for interpolation
        
    Returns:
    --------
    uniform_grid : np.array
        Uniform grid of time values
    column_names : list
        List of column names
    """
    # Sample files for grid determination
    sample_files = np.random.choice(replica_files, 
                                   size=min(sample_size, len(replica_files)), 
                                   replace=False)
    
    print(f"Sampling {len(sample_files)} files to determine uniform time grid...")
    
    time_arrays = []
    column_names = None
    
    for filepath in tqdm(sample_files, desc="Sampling files", unit="file"):
        time, data, cols = read_energy_file(filepath)
        if time is not None:
            time_arrays.append(time)
            if column_names is None:
                column_names = cols
    
    if not time_arrays:
        raise ValueError("No valid files found in sample!")
    
    uniform_grid = create_uniform_time_grid(time_arrays, grid_points)
    
    return uniform_grid, column_names


def create_uniform_time_grid_from_range(min_time, max_time, num_points=1000):
    """
    Create a uniform time grid from specified time range.
    
    Parameters:
    -----------
    min_time : float
        Minimum time value
    max_time : float
        Maximum time value
    num_points : int
        Number of grid points
        
    Returns:
    --------
    uniform_grid : np.array
        Uniform grid of time values
    """
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
    # Process each column
    for col_idx in range(new_data.shape[1]):
        # Only update points with valid data for this column
        valid_points = valid_mask & ~np.isnan(new_data[:, col_idx])
        
        if np.any(valid_points):
            # Update count for valid points
            count_per_point[valid_points] += 1
            
            # Update mean and variance for valid points
            delta = new_data[valid_points, col_idx] - mean[valid_points, col_idx]
            mean[valid_points, col_idx] += delta / count_per_point[valid_points]
            delta2 = new_data[valid_points, col_idx] - mean[valid_points, col_idx]
            variance[valid_points, col_idx] += delta * delta2
    
    return mean, variance, count_per_point


def find_energy_replica_files(folder_path):
    """
    Find all energy replica files in the specified folder.
    
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
    return glob.glob(pattern)


def average_energy_replicas(folder_path, grid_points=None, sample_size=50, max_time=None, min_time=None):
    """
    Average energy replica files using memory-efficient online averaging.
    
    Parameters:
    -----------
    folder_path : str
        Path to the folder containing replica files
    grid_points : int, optional
        Number of grid points for interpolation
    sample_size : int
        Number of files to sample for grid determination
    max_time : float, optional
        Maximum time for the uniform grid
    min_time : float, optional
        Minimum time for the uniform grid (default: 0.0)
    """
    # Find all energy replica files
    replica_files = find_energy_replica_files(folder_path)
    
    if not replica_files:
        print(f"No files matching pattern 'prod-*_energy_tracker.txt' found in {folder_path}")
        return False
    
    print(f"\nProcessing energy replica files in {folder_path}:")
    print(f"Found {len(replica_files)} replica files:")
    for f in sorted(replica_files):
        print(f"  - {os.path.basename(f)}")
    
    # Determine uniform time grid and get column names
    if max_time is not None:
        # Use user-specified time range
        if min_time is None:
            min_time = 0.0
        if grid_points is None:
            grid_points = 1000
        
        print(f"\nCreating uniform time grid from user-specified range...")
        print(f"Time range: {min_time:.6f} to {max_time:.6f} ps")
        print(f"Grid points: {grid_points}")
        
        uniform_grid = create_uniform_time_grid_from_range(min_time, max_time, grid_points)
        
        # Get column names from first file
        _, _, column_names = read_energy_file(replica_files[0])
        if column_names is None:
            print("Error: Could not read column names from first file")
            return False
    else:
        # Use sampling approach
        print(f"\nDetermining uniform time grid using {sample_size} sample files...")
        uniform_grid, column_names = determine_time_grid_from_sample(replica_files, sample_size, grid_points)
    
    if column_names is None:
        print("Error: Could not determine column names")
        return False
    
    print(f"Created uniform time grid with {len(uniform_grid)} points")
    print(f"Grid range: {uniform_grid.min():.6f} to {uniform_grid.max():.6f} ps")
    print(f"Columns to average: {len(column_names)}")
    
    # Initialize online statistics
    n_grid_points = len(uniform_grid)
    n_columns = len(column_names)
    running_mean = np.zeros((n_grid_points, n_columns))
    running_variance = np.zeros((n_grid_points, n_columns))
    count_per_point = np.zeros(n_grid_points, dtype=int)  # Track samples per time point
    
    # Process replicas one by one for memory efficiency
    print(f"\nProcessing all energy replicas...")
    successful_replicas = 0
    partial_replicas = 0
    
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
                    
                    # Count as successful if we have data for most of the grid
                    coverage = np.sum(valid_mask) / len(valid_mask)
                    if coverage > 0.8:  # More than 80% coverage
                        successful_replicas += 1
                    else:
                        partial_replicas += 1
                else:
                    print(f"No valid data points for file: {filepath}")
            else:
                print(f"Failed to interpolate file: {filepath}")
        else:
            print(f"Failed to read file: {filepath}")
    
    if successful_replicas == 0 and partial_replicas == 0:
        print("No successful interpolations!")
        return False
    
    # Calculate final standard deviation with variable sample counts
    std_data = np.zeros_like(running_mean)
    for i in range(n_grid_points):
        for j in range(n_columns):
            if count_per_point[i] > 1:
                std_data[i, j] = np.sqrt(running_variance[i, j] / (count_per_point[i] - 1))
            else:
                std_data[i, j] = 0.0
    
    print(f"Successfully averaged {successful_replicas} full replicas and {partial_replicas} partial replicas")
    print(f"Sample count per time point: min={count_per_point.min()}, max={count_per_point.max()}, median={np.median(count_per_point):.0f}")
    
    # Save the averaged result
    output_file = os.path.join(folder_path, "prod_energy_tracker_averaged.txt")
    
    # Create header similar to original files
    header = f"""# Averaged Energy tracking from {successful_replicas + partial_replicas} replicas
# Output period: uniform time grid
# Variable sample count per time point (min={count_per_point.min()}, max={count_per_point.max()})
# All energies in Hartree (atomic units)
# Column definitions:
#   time(ps): simulation time in picoseconds
#   timestep: simulation timestep number (averaged)
#   harmonic_energy: harmonic bond potential energy
#   lj_energy: Lennard-Jones potential energy
#   ewald_short_energy: short-range Coulomb energy
#   ewald_long_energy: long-range Coulomb energy
#   cavity_harmonic_energy: cavity harmonic potential energy
#   cavity_coupling_energy: cavity-molecule coupling energy
#   cavity_dipole_self_energy: dipole self-energy
#   cavity_total_potential_energy: total cavity potential energy
#   molecular_kinetic_energy: molecular kinetic energy
#   cavity_kinetic_energy: cavity kinetic energy
#   total_kinetic_energy: total kinetic energy
#   total_potential_energy: total potential energy
#   system_total_energy: total system energy (KE + PE)
#   molecular_reservoir_energy: molecular reservoir energy
#   cavity_reservoir_energy: cavity reservoir energy
#   total_reservoir_energy: total reservoir energy
#   universe_total_energy: universe total energy (system + reservoir) [CONSERVED]
#   temperature: kinetic temperature (K)
# NOTE: Standard deviations and sample counts are appended for each column
"""
    
    # Create extended column names (original + std + count)
    extended_column_names = []
    for col_name in column_names:
        extended_column_names.append(col_name)
    for col_name in column_names:
        extended_column_names.append(f"{col_name}_std")
    extended_column_names.append("sample_count")
    
    header += " ".join(extended_column_names)
    
    # Create output data: time, all columns, all std columns, sample count
    output_data = np.column_stack([
        running_mean,  # All averaged columns
        std_data,      # All standard deviation columns
        count_per_point.reshape(-1, 1)  # Sample count per time point
    ])
    
    # Save to file with progress bar
    print(f"\nSaving averaged energy data...")
    with tqdm(desc="Writing output", unit="lines") as pbar:
        np.savetxt(output_file, output_data, 
                   header=header, 
                   fmt='%.6f',
                   delimiter=' ')
        pbar.update(len(output_data))
    
    print(f"Averaged data saved to: {output_file}")
    print(f"Statistics:")
    print(f"  - Time range: {uniform_grid.min():.6f} to {uniform_grid.max():.6f} ps")
    print(f"  - Sample count range: {count_per_point.min()} to {count_per_point.max()}")
    print(f"  - Columns averaged: {len(column_names)}")
    print(f"  - Memory usage: Processed {successful_replicas + partial_replicas} replicas without storing all data in memory")
    
    return True


def main():
    """Main function to parse arguments and run the averaging."""
    parser = argparse.ArgumentParser(
        description="Average multiple energy replica files (prod-*_energy_tracker.txt) by interpolating onto a uniform time grid."
    )
    parser.add_argument(
        "folder_path", 
        type=str, 
        help="Path to the folder containing energy replica files"
    )
    parser.add_argument(
        "--grid-points", 
        type=int, 
        default=None,
        help="Number of grid points for interpolation (default: median of replica lengths or 1000 if max-time specified)"
    )
    parser.add_argument(
        "--sample-size", 
        type=int, 
        default=50,
        help="Number of files to sample for grid determination (default: 50, ignored if max-time specified)"
    )
    parser.add_argument(
        "--max-time", 
        type=float, 
        default=None,
        help="Maximum time for the uniform grid (ps). If specified, no sampling needed."
    )
    parser.add_argument(
        "--min-time", 
        type=float, 
        default=None,
        help="Minimum time for the uniform grid (ps). Default: 0.0 if max-time specified."
    )
    
    args = parser.parse_args()
    
    # Check if folder exists
    if not os.path.isdir(args.folder_path):
        print(f"Error: Folder '{args.folder_path}' does not exist!")
        sys.exit(1)
    
    # Run averaging
    success = average_energy_replicas(args.folder_path, args.grid_points, args.sample_size, args.max_time, args.min_time)
    
    if success:
        print("\nEnergy replica averaging completed successfully!")
    else:
        print("\nEnergy replica averaging failed!")
        sys.exit(1)


if __name__ == "__main__":
    main() 