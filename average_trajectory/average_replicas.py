#!/usr/bin/env python3
"""
Script to average multiple replica files (prod-*_ref*.txt) by interpolating onto a uniform grid.

Usage: python average_replicas.py <folder_path>

The script will:
1. Automatically detect all ref numbers (ref0, ref1, ref2, etc.) in the given folder
2. For each ref number, find all files matching pattern prod-*_refN.txt
3. Extract lag_time and field_autocorr data from each file
4. Create a uniform grid for lag_time based on the common range
5. Interpolate each replica's data onto this grid
6. Average the interpolated data across all replicas using online averaging
7. Save the results as prod_refN.txt (one file for each ref number) in the same folder
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


def read_replica_file(filepath):
    """
    Read a replica file and extract lag_time and field_autocorr data.
    
    Parameters:
    -----------
    filepath : str
        Path to the replica file
        
    Returns:
    --------
    lag_time : np.array
        Array of lag times
    field_autocorr : np.array
        Array of field autocorrelation values
    """
    try:
        # Read the file, skipping header lines that start with '#'
        data = pd.read_csv(filepath, sep='\s+', comment='#', header=None, 
                          names=['timestep', 'lag_time', 'field_autocorr'])
        
        return data['lag_time'].values, data['field_autocorr'].values
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return None, None


def create_uniform_grid(lag_times_list, num_points=None):
    """
    Create a uniform grid based on the common range of lag times.
    
    Parameters:
    -----------
    lag_times_list : list of np.arrays
        List of lag time arrays from all replicas
    num_points : int, optional
        Number of grid points. If None, uses the median number of points from replicas
        
    Returns:
    --------
    uniform_grid : np.array
        Uniform grid of lag times
    """
    # Find the common range (intersection of all ranges)
    min_lag = max([lag_times.min() for lag_times in lag_times_list])
    max_lag = min([lag_times.max() for lag_times in lag_times_list])
    
    if num_points is None:
        # Use median number of points from all replicas
        num_points = int(np.median([len(lag_times) for lag_times in lag_times_list]))
    
    return np.linspace(min_lag, max_lag, num_points)


def determine_grid_from_sample(replica_files, sample_size=50, grid_points=None):
    """
    Determine uniform grid by sampling a subset of files to save memory.
    
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
        Uniform grid of lag times
    """
    # Sample files for grid determination
    sample_files = np.random.choice(replica_files, 
                                   size=min(sample_size, len(replica_files)), 
                                   replace=False)
    
    print(f"Sampling {len(sample_files)} files to determine uniform grid...")
    
    lag_times_list = []
    for filepath in tqdm(sample_files, desc="Sampling files", unit="file"):
        lag_time, _ = read_replica_file(filepath)  # Only read lag_time, ignore field_autocorr
        if lag_time is not None:
            lag_times_list.append(lag_time)
    
    if not lag_times_list:
        raise ValueError("No valid files found in sample!")
    
    return create_uniform_grid(lag_times_list, grid_points)


def create_uniform_grid_from_range(min_time, max_time, num_points=1000):
    """
    Create a uniform grid from specified time range.
    
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
        Uniform grid of lag times
    """
    return np.linspace(min_time, max_time, num_points)


def interpolate_to_grid(lag_time, field_autocorr, uniform_grid):
    """
    Interpolate field_autocorr data onto a uniform grid.
    
    Parameters:
    -----------
    lag_time : np.array
        Original lag time array
    field_autocorr : np.array
        Original field autocorrelation array
    uniform_grid : np.array
        Target uniform grid
        
    Returns:
    --------
    interpolated_data : np.array
        Field autocorrelation values interpolated onto the uniform grid
    valid_mask : np.array
        Boolean mask indicating which grid points have valid interpolated data
    """
    try:
        # Find the range of the current file
        min_time = lag_time.min()
        max_time = lag_time.max()
        
        # Create mask for grid points within the file's range
        valid_mask = (uniform_grid >= min_time) & (uniform_grid <= max_time)
        
        # Initialize output array with NaN
        interpolated_data = np.full_like(uniform_grid, np.nan)
        
        if np.any(valid_mask):
            # Create interpolation function
            interp_func = interp1d(lag_time, field_autocorr, kind='linear', 
                                  bounds_error=False, fill_value=np.nan, assume_sorted=True)
            
            # Interpolate only for valid grid points
            interpolated_data[valid_mask] = interp_func(uniform_grid[valid_mask])
        
        return interpolated_data, valid_mask
    except Exception as e:
        print(f"Error during interpolation: {e}")
        return None, None


def update_online_stats_with_mask(mean, variance, count_per_point, new_value, valid_mask):
    """
    Update running mean and variance using Welford's online algorithm with variable sample counts.
    
    Parameters:
    -----------
    mean : np.array
        Current running mean
    variance : np.array  
        Current running variance
    count_per_point : np.array
        Current count of samples per time point
    new_value : np.array
        New data point to incorporate (may contain NaN)
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
    # Only update points with valid data
    valid_points = valid_mask & ~np.isnan(new_value)
    
    if np.any(valid_points):
        # Update count for valid points
        count_per_point[valid_points] += 1
        
        # Update mean and variance for valid points
        delta = new_value[valid_points] - mean[valid_points]
        mean[valid_points] += delta / count_per_point[valid_points]
        delta2 = new_value[valid_points] - mean[valid_points]
        variance[valid_points] += delta * delta2
    
    return mean, variance, count_per_point


def detect_ref_numbers(folder_path):
    """
    Detect all ref numbers present in the folder.
    
    Parameters:
    -----------
    folder_path : str
        Path to the folder containing replica files
        
    Returns:
    --------
    ref_numbers : list
        List of unique ref numbers found
    """
    # Find all files matching pattern prod-*_ref*.txt
    pattern = os.path.join(folder_path, "prod-*_ref*.txt")
    all_files = glob.glob(pattern)
    
    ref_numbers = set()
    ref_pattern = re.compile(r'_ref(\d+)\.txt$')
    
    for filepath in all_files:
        match = ref_pattern.search(filepath)
        if match:
            ref_numbers.add(int(match.group(1)))
    
    return sorted(ref_numbers)


def find_replica_files_for_ref(folder_path, ref_number):
    """
    Find all replica files for a specific ref number.
    
    Parameters:
    -----------
    folder_path : str
        Path to the folder containing replica files
    ref_number : int
        Reference number to find files for
        
    Returns:
    --------
    replica_files : list
        List of file paths for the specified ref number
    """
    pattern = os.path.join(folder_path, f"prod-*_ref{ref_number}.txt")
    return glob.glob(pattern)


def average_replicas_for_ref(folder_path, ref_number, grid_points=None, sample_size=50, max_time=None, min_time=None):
    """
    Average replica files for a specific ref number using memory-efficient online averaging.
    
    Parameters:
    -----------
    folder_path : str
        Path to the folder containing replica files
    ref_number : int
        Reference number to process
    grid_points : int, optional
        Number of grid points for interpolation
    sample_size : int
        Number of files to sample for grid determination
    max_time : float, optional
        Maximum time for the uniform grid
    min_time : float, optional
        Minimum time for the uniform grid (default: 0.0)
    """
    # Find all replica files for this ref number
    replica_files = find_replica_files_for_ref(folder_path, ref_number)
    
    if not replica_files:
        print(f"No files matching pattern 'prod-*_ref{ref_number}.txt' found in {folder_path}")
        return False
    
    print(f"\nProcessing ref{ref_number} files:")
    print(f"Found {len(replica_files)} replica files:")
    for f in sorted(replica_files):
        print(f"  - {os.path.basename(f)}")
    
    # Determine uniform grid
    if max_time is not None:
        # Use user-specified time range
        if min_time is None:
            min_time = 0.0
        if grid_points is None:
            grid_points = 1000
        
        print(f"\nCreating uniform grid from user-specified range...")
        print(f"Time range: {min_time:.6f} to {max_time:.6f} ps")
        print(f"Grid points: {grid_points}")
        
        uniform_grid = create_uniform_grid_from_range(min_time, max_time, grid_points)
    else:
        # Use sampling approach (fallback)
        print(f"\nDetermining uniform grid using {sample_size} sample files...")
        uniform_grid = determine_grid_from_sample(replica_files, sample_size, grid_points)
    
    print(f"Created uniform grid with {len(uniform_grid)} points")
    print(f"Grid range: {uniform_grid.min():.6f} to {uniform_grid.max():.6f} ps")
    
    # Initialize online statistics
    n_grid_points = len(uniform_grid)
    running_mean = np.zeros(n_grid_points)
    running_variance = np.zeros(n_grid_points)
    count_per_point = np.zeros(n_grid_points, dtype=int)  # Track samples per time point
    
    # Process replicas one by one for memory efficiency
    print(f"\nProcessing all ref{ref_number} replicas...")
    successful_replicas = 0
    partial_replicas = 0
    
    for filepath in tqdm(replica_files, desc=f"Processing ref{ref_number} replicas", unit="file"):
        # Read current replica
        lag_time, field_autocorr = read_replica_file(filepath)
        
        if lag_time is not None and field_autocorr is not None:
            # Interpolate current replica
            interp_data, valid_mask = interpolate_to_grid(lag_time, field_autocorr, uniform_grid)
            
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
        print(f"No successful interpolations for ref{ref_number}!")
        return False
    
    # Calculate final standard deviation with variable sample counts
    std_field_autocorr = np.zeros_like(running_mean)
    for i in range(n_grid_points):
        if count_per_point[i] > 1:
            std_field_autocorr[i] = np.sqrt(running_variance[i] / (count_per_point[i] - 1))
        else:
            std_field_autocorr[i] = 0.0
    
    print(f"Successfully averaged {successful_replicas} full replicas and {partial_replicas} partial replicas for ref{ref_number}")
    print(f"Sample count per time point: min={count_per_point.min()}, max={count_per_point.max()}, median={np.median(count_per_point):.0f}")
    
    # Save the averaged result
    output_file = os.path.join(folder_path, f"prod_ref{ref_number}.txt")
    
    # Create header similar to original files
    header = f"""# Density_correlation field autocorrelation (averaged from {successful_replicas + partial_replicas} replicas)
# Reference {ref_number} at t=0.000000 ps
# Output period: uniform grid
# Variable sample count per time point (min={count_per_point.min()}, max={count_per_point.max()})
# timestep lag_time(ps) field_autocorr field_autocorr_std sample_count"""
    
    # Create output data
    output_data = np.column_stack([
        np.arange(1, len(uniform_grid) + 1),  # timestep
        uniform_grid,                          # lag_time
        running_mean,                          # averaged field_autocorr
        std_field_autocorr,                    # standard deviation
        count_per_point                        # sample count per time point
    ])
    
    # Save to file with progress bar
    print(f"\nSaving averaged ref{ref_number} data...")
    with tqdm(desc="Writing output", unit="lines") as pbar:
        np.savetxt(output_file, output_data, 
                   header=header, 
                   fmt=['%d', '%.6f', '%.6f', '%.6f', '%d'],
                   delimiter=' ')
        pbar.update(len(output_data))
    
    print(f"Averaged data saved to: {output_file}")
    print(f"Statistics for ref{ref_number}:")
    print(f"  - Mean autocorr range: {running_mean.min():.2f} to {running_mean.max():.2f}")
    print(f"  - Average std dev: {std_field_autocorr.mean():.2f}")
    print(f"  - Sample count range: {count_per_point.min()} to {count_per_point.max()}")
    print(f"  - Memory usage: Processed {successful_replicas + partial_replicas} replicas without storing all data in memory")
    
    return True


def average_replicas(folder_path, grid_points=None, sample_size=50, max_time=None, min_time=None):
    """
    Main function to average replica files for all detected ref numbers using memory-efficient online averaging.
    
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
    # Detect all ref numbers in the folder
    ref_numbers = detect_ref_numbers(folder_path)
    
    if not ref_numbers:
        print(f"No files matching pattern 'prod-*_ref*.txt' found in {folder_path}")
        return
    
    print(f"Detected ref numbers: {ref_numbers}")
    print(f"Will process {len(ref_numbers)} reference sets")
    
    successful_refs = 0
    
    # Process each ref number separately
    for ref_number in ref_numbers:
        print(f"\n{'='*60}")
        print(f"Processing Reference {ref_number}")
        print(f"{'='*60}")
        
        success = average_replicas_for_ref(folder_path, ref_number, grid_points, sample_size, max_time, min_time)
        if success:
            successful_refs += 1
    
    # Final summary
    print(f"\n{'='*60}")
    print("PROCESSING COMPLETE")
    print(f"{'='*60}")
    print(f"Successfully processed {successful_refs}/{len(ref_numbers)} reference sets")
    
    if successful_refs > 0:
        print("Generated files:")
        for ref_number in ref_numbers:
            output_file = os.path.join(folder_path, f"prod_ref{ref_number}.txt")
            if os.path.exists(output_file):
                print(f"  - {output_file}")


def main():
    """Main function to parse arguments and run the averaging."""
    parser = argparse.ArgumentParser(
        description="Average multiple replica files (prod-*_ref*.txt) by interpolating onto a uniform grid. Automatically detects and processes all ref numbers found in the directory."
    )
    parser.add_argument(
        "folder_path", 
        type=str, 
        help="Path to the folder containing replica files"
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
    average_replicas(args.folder_path, args.grid_points, args.sample_size, args.max_time, args.min_time)


if __name__ == "__main__":
    main()
