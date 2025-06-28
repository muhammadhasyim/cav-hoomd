#!/usr/bin/env python3
"""
F(k,t) Processing and Analysis Script for Cavity MD Simulations

This script processes F(k,t) files from cavity MD simulations, creating master files
by averaging across all replicas and generating comparison plots.

Key Features:
- Configurable timestep resolution for F(k,t) data
- Processes all experiment directories (with optional parallelization)
- Generates master F(k,t) files by averaging across all replicas
- Creates comprehensive F(k,t) comparison plots for different coupling strengths

Usage:
    python process_fskt_only.py --base_dir ./
    python process_fskt_only.py --fskt_dt 1.0 --max_time 50.0
    python process_fskt_only.py --processes 4  # Use 4 parallel processes
"""
 
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import glob
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.interpolate import interp1d
import time
import natsort
from pathlib import Path
import argparse
import logging
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
import sys
import re

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def setup_worker_logging():
    """Set up logging for worker processes."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s'
    )

def extract_k_value_from_header(file_path):
    """
    Extract k-value from F(k,t) file header.
    
    Args:
        file_path: Path to the F(k,t) file
        
    Returns:
        String representation of k-value (e.g., "k=1.0") or None if not found
    """
    try:
        with open(file_path, 'r') as f:
            # Read first few lines to find k-value in header
            for i, line in enumerate(f):
                if i > 10:  # Limit search to first 10 lines
                    break
                if line.startswith('#') and 'kmag' in line.lower():
                    # Look for patterns like "kmag = 1.0" or "kmag=1.0"
                    match = re.search(r'kmag\s*=\s*([\d.]+)', line.lower())
                    if match:
                        k_val = float(match.group(1))
                        return f"k={k_val:.2f}"
        return None
    except Exception as e:
        logger.warning(f"Could not extract k-value from {file_path}: {e}")
        return None

def remove_duplicate_times(times, values, tolerance=1e-12):
    """
    Remove duplicate time values, keeping the first occurrence.
    
    Args:
        times: Array of time values
        values: Array of corresponding values
        tolerance: Tolerance for considering times as duplicates
        
    Returns:
        Tuple of (unique_times, unique_values)
    """
    if len(times) <= 1:
        return times, values
    
    # Convert to pandas for easier duplicate handling
    df = pd.DataFrame({'time': times, 'value': values})
    
    # Round times to handle floating point precision issues
    df['time_rounded'] = np.round(df['time'] / tolerance) * tolerance
    
    # Remove duplicates, keeping first occurrence
    df_unique = df.drop_duplicates(subset=['time_rounded'], keep='first')
    
    # Sort by time to ensure monotonic order
    df_unique = df_unique.sort_values('time')
    
    return df_unique['time'].values, df_unique['value'].values

@dataclass
class ExperimentConfig:
    """Configuration for a single experiment type."""
    name: str
    molecular_thermostat: str
    cavity_thermostat: str
    coupling_strength: str
    job_name: str  # 'prod'
    display_name: str

def process_fskt_files(exp_dir, job_name, max_time=None, dt=1.0):
    """
    Process F(k,t) files using configurable timestep and proper cumulative averaging.
    
    Args:
        exp_dir: Path to experiment directory
        job_name: Job name ('prod' or 'finq')
        max_time: Maximum time in ps for F(k,t) data (None = auto-detect from data)
        dt: Timestep in ps for interpolation
    """
    # Setup logging for this process
    process_logger = logging.getLogger(f"Worker-{mp.current_process().name}")
    
    exp_path = Path(exp_dir)
    
    # Find all F(k,t) files for this experiment - only process ref0 files
    fskt_pattern = f"{job_name}-*_ref0.txt"
    fskt_files = list(exp_path.glob(fskt_pattern))
    
    if not fskt_files:
        process_logger.warning(f"No F(k,t) ref0 files found in {exp_dir} with pattern {fskt_pattern}")
        return False
    
    # Group files by k-value (extracted from file headers)
    fskt_groups = {}
    for fskt_file in fskt_files:
        # Extract k-value from file header
        k_value = extract_k_value_from_header(fskt_file)
        if k_value is None:
            # Fallback: try to extract from filename or use default
            filename = fskt_file.name
            # Look for patterns like prod-1_ref0.txt, prod-2_ref0.txt, etc.
            match = re.search(r'prod-(\d+)_ref0\.txt', filename)
            if match:
                replica_num = match.group(1)
                k_value = f"k=1.00"  # Default k-value
                process_logger.info(f"Using default k-value for {filename}: {k_value}")
            else:
                k_value = "k=1.00"  # Final fallback
                process_logger.warning(f"Could not determine k-value for {filename}, using default: {k_value}")
        
        if k_value not in fskt_groups:
            fskt_groups[k_value] = []
        fskt_groups[k_value].append(str(fskt_file))
    
    process_logger.info(f"Processing {len(fskt_groups)} F(k,t) ref0 groups in {exp_dir}")
    success_count = 0
    
    # Process fskt groups (removing nested progress bar for cleaner parallel output)
    for k_value, file_list in fskt_groups.items():
        try:
            # Sort files naturally
            file_list = natsort.natsorted(file_list)
            
            # Auto-detect maximum time range from all files if not specified
            detected_max_time = max_time
            if max_time is None:
                all_max_times = []
                # Process time detection
                for fskt_file in file_list:
                    try:
                        temp_data = pd.read_csv(fskt_file, sep=r'\s+', comment='#', header=None)
                        if not temp_data.empty and len(temp_data.columns) >= 2:
                            # Skip header rows if present
                            first_row = temp_data.iloc[0]
                            if any(isinstance(val, str) and not val.replace('.', '').replace('-', '').replace('e', '').replace('+', '').isdigit() for val in first_row):
                                temp_data = temp_data.iloc[1:]
                            
                            temp_data = temp_data.apply(pd.to_numeric, errors='coerce').dropna()
                            if not temp_data.empty:
                                # For 5-column format: t0_step t0(ps) t_step t(ps) F(k,t)
                                # Use column 3 (index 3) which is t(ps)
                                if len(temp_data.columns) >= 4:
                                    times = temp_data.iloc[:, 3].values  # t(ps) column
                                else:
                                    times = temp_data.iloc[:, 0].values  # Fallback to first column
                                valid_times = times[times >= 0]
                                if len(valid_times) > 0:
                                    all_max_times.append(valid_times.max())
                    except Exception as e:
                        process_logger.warning(f"Error reading {fskt_file} for time detection: {e}")
                        continue
                
                if all_max_times:
                    # Use the maximum time found across all replicas to include all available data
                    detected_max_time = max(all_max_times)
                    process_logger.info(f"Auto-detected max time for {k_value}: {detected_max_time:.3f} ps from {len(all_max_times)} files")
                    process_logger.info(f"  Time range across files: {min(all_max_times):.3f} to {max(all_max_times):.3f} ps")
                else:
                    detected_max_time = 50.0  # fallback
                    process_logger.warning(f"Could not detect time range for {k_value}, using fallback: {detected_max_time} ps")
            
            # Set fixed linspace with configurable timestep
            if detected_max_time is not None:
                num_points = int(detected_max_time / dt) + 1
                uniform_times = np.linspace(0, detected_max_time, num_points)
                process_logger.info(f"F(k,t) group {k_value}: Using time grid: 0 to {detected_max_time:.3f} ps with {len(uniform_times)} points (dt = {dt} ps)")
            else:
                process_logger.error(f"Could not determine time range for {k_value}")
                continue
            
            # Initialize storage for cumulative data
            cumulative_fskt = np.zeros(len(uniform_times))
            sample_counts = np.zeros(len(uniform_times), dtype=int)
            processed_count = 0
            
            # Process individual files
            for fskt_file in file_list:
                try:
                    # Read F(k,t) data
                    data = pd.read_csv(fskt_file, sep=r'\s+', comment='#', header=None)
                    
                    if data.empty:
                        continue
                    
                    # Skip header rows if present
                    first_row = data.iloc[0]
                    if any(isinstance(val, str) and not val.replace('.', '').replace('-', '').replace('e', '').replace('+', '').isdigit() for val in first_row):
                        data = data.iloc[1:]
                    
                    # Convert to numeric and drop invalid rows
                    data = data.apply(pd.to_numeric, errors='coerce').dropna()
                    
                    if data.empty or len(data.columns) < 2:
                        continue
                    
                    # Handle both 5-column and 2-column formats
                    if len(data.columns) >= 5:
                        # 5-column format: t0_step t0(ps) t_step t(ps) F(k,t)
                        # Use columns 3 and 4 (t(ps) and F(k,t))
                        time_col = data.iloc[:, 3]  # t(ps)
                        fskt_col = data.iloc[:, 4]  # F(k,t)
                        data = pd.DataFrame({'t(ps)': time_col, 'F(k,t)': fskt_col})
                    else:
                        # 2-column format: use last two columns
                        data = data.iloc[:, -2:]
                        data.columns = ['t(ps)', 'F(k,t)']
                    
                    times = data['t(ps)'].values
                    fskt_vals = data['F(k,t)'].values
                    
                    # Filter to valid time range (now using detected max time)
                    valid_mask = (times >= 0) & (times <= detected_max_time)
                    if np.sum(valid_mask) < 3:
                        process_logger.warning(f"Insufficient valid time points in {fskt_file}: {np.sum(valid_mask)} points")
                        continue
                    
                    times = times[valid_mask]
                    fskt_vals = fskt_vals[valid_mask]
                    
                    # Log the actual time range used for this file
                    file_max_time = times.max() if len(times) > 0 else 0
                    process_logger.debug(f"File {Path(fskt_file).name}: using {len(times)} points, max time: {file_max_time:.3f} ps")
                    
                    # Remove duplicate time values before interpolation
                    unique_times, unique_fskt = remove_duplicate_times(times, fskt_vals)
                    
                    if len(unique_times) < 2:
                        process_logger.warning(f"Insufficient unique time points in {fskt_file}: {len(unique_times)} points")
                        continue
                    
                    # Create interpolation function
                    if len(unique_times) > 3:
                        interp_func = interp1d(unique_times, unique_fskt, kind='cubic',
                                             bounds_error=False, fill_value=np.nan)
                    else:
                        interp_func = interp1d(unique_times, unique_fskt, kind='linear',
                                             bounds_error=False, fill_value=np.nan)
                    
                    # Evaluate at uniform time points
                    interpolated_fskt = interp_func(uniform_times)
                    valid_mask = ~np.isnan(interpolated_fskt)
                    
                    # Update sample counts only for valid points
                    sample_counts[valid_mask] += 1
                    
                    # Update cumulative average using proper online algorithm
                    for i in range(len(uniform_times)):
                        if valid_mask[i]:
                            n = sample_counts[i]
                            if n == 1:
                                cumulative_fskt[i] = interpolated_fskt[i]
                            else:
                                # Online average update: new_avg = old_avg + (new_val - old_avg) / n
                                old_avg = cumulative_fskt[i]
                                new_val = interpolated_fskt[i]
                                cumulative_fskt[i] = old_avg + (new_val - old_avg) / n
                    
                    processed_count += 1
                        
                except Exception as e:
                    process_logger.warning(f"Error processing {fskt_file}: {e}")
                    continue
            
            # Save results if we have data
            if processed_count > 0:
                master_file = exp_path / f"master_fskt_{k_value}.txt"
                sample_count_file = exp_path / f"master_fskt_{k_value}_sample_counts.txt"
                
                # Only keep time points where we have at least one sample
                valid_points = sample_counts > 0
                clean_times = uniform_times[valid_points]
                clean_fskt = cumulative_fskt[valid_points]
                clean_sample_counts = sample_counts[valid_points]
                
                # Create output data
                output_data = pd.DataFrame({
                    'lag_time': clean_times,
                    'fskt': clean_fskt
                })
                
                # Save master data with comprehensive header
                with open(master_file, 'w') as f:
                    f.write(f"# Master F(k,t) file for {exp_path.name} - {k_value}\n")
                    f.write(f"# Cumulative average across {processed_count} replicas\n")
                    f.write(f"# Time grid: 0 to {detected_max_time:.3f} ps with {len(uniform_times)} points (dt = {dt} ps)\n")
                    f.write(f"# Valid data points: {len(clean_times)} (with at least 1 sample)\n")
                    f.write(f"# Sample count range: {clean_sample_counts.min()} to {clean_sample_counts.max()}\n")
                    f.write(f"# Generated by process_fskt_only.py\n")
                    f.write(f"# Columns: lag_time fskt\n")
                
                output_data.to_csv(master_file, sep='\t', index=False, mode='a', header=True, na_rep='nan')
                np.savetxt(sample_count_file, clean_sample_counts, fmt='%d')
                
                process_logger.info(f"Created master F(k,t) file: {master_file} ({processed_count} replicas)")
                process_logger.info(f"  Valid time points: {len(clean_times)}/{len(uniform_times)}")
                process_logger.info(f"  Time range: {clean_times.min():.3f} - {clean_times.max():.3f} ps")
                process_logger.info(f"  Sample count: min={clean_sample_counts.min()}, max={clean_sample_counts.max()}, mean={clean_sample_counts.mean():.1f}")
                success_count += 1
                
        except Exception as e:
            process_logger.error(f"Error processing F(k,t) group {k_value}: {e}")
            continue
    
    return success_count > 0

def process_single_experiment(args_tuple):
    """
    Wrapper function for processing a single experiment in parallel.
    
    Args:
        args_tuple: Tuple of (exp_dir, job_name, max_time, dt)
    
    Returns:
        Tuple of (exp_dir_name, success_bool)
    """
    exp_dir, job_name, max_time, dt = args_tuple
    
    # Set up logging for this worker process
    setup_worker_logging()
    process_logger = logging.getLogger(f"Worker-{mp.current_process().name}")
    
    try:
        process_logger.info(f"Starting processing of {exp_dir}")
        success = process_fskt_files(exp_dir, job_name, max_time, dt)
        
        if success:
            process_logger.info(f"✅ Successfully processed {exp_dir}")
        else:
            process_logger.warning(f"❌ Failed to process {exp_dir}")
            
        return (Path(exp_dir).name, success)
        
    except Exception as e:
        process_logger.error(f"Error processing {exp_dir}: {e}")
        return (Path(exp_dir).name, False)

def process_all_experiments(base_dir="./", max_time=None, dt=1.0, num_processes=None):
    """
    Process all experiment directories to generate master F(k,t) files.
    
    Args:
        base_dir: Base directory containing experiment folders
        max_time: Maximum time for F(k,t) data (ps) - None for auto-detection
        dt: Timestep for F(k,t) data interpolation (ps)
        num_processes: Number of parallel processes to use (None = auto-detect)
    """
    base_path = Path(base_dir)
    
    # Define expected experiment patterns based on current directory structure
    experiment_patterns = [
        "bussi_bussi_*finiteq*",
        "bussi_langevin_*finiteq*",
        "*_finiteq*",  # Catch any other finiteq experiments
        "bussi_bussi_no_finiteq",
        "bussi_langevin_no_finiteq_coupling_*",
        "*_no_finiteq*"  # Catch any other no_finiteq experiments
    ]
    
    # Find all experiment directories
    experiment_dirs = []
    for pattern in experiment_patterns:
        experiment_dirs.extend(base_path.glob(pattern))
    
    if not experiment_dirs:
        logger.error("No experiment directories found!")
        return False
    
    experiment_dirs = sorted(experiment_dirs)
    logger.info(f"Found {len(experiment_dirs)} experiment directories")
    
    if max_time is None:
        logger.info("Using auto-detection of time ranges from data for each experiment")
    else:
        logger.info(f"Using fixed time range: 0 to {max_time} ps for all experiments")
    
    # Determine number of processes to use
    if num_processes is None:
        num_processes = min(len(experiment_dirs), mp.cpu_count())
    
    logger.info(f"Using {num_processes} parallel processes")
    
    # Prepare arguments for parallel processing
    job_name = "prod"  # All current experiments use 'prod' job name
    process_args = [(str(exp_dir), job_name, max_time, dt) for exp_dir in experiment_dirs]
    
    if num_processes == 1:
        # Sequential processing (no multiprocessing overhead)
        logger.info("Running in sequential mode")
        results = []
        for args in tqdm(process_args, desc="Processing experiments"):
            results.append(process_single_experiment(args))
    else:
        # Parallel processing
        logger.info(f"Running in parallel mode with {num_processes} processes")
        
        # Use multiprocessing Pool
        with mp.Pool(processes=num_processes, initializer=setup_worker_logging) as pool:
            # Use tqdm to show progress
            results = list(tqdm(
                pool.imap(process_single_experiment, process_args),
                total=len(process_args),
                desc="Processing experiments"
            ))
    
    # Analyze results
    successful_experiments = [name for name, success in results if success]
    failed_experiments = [name for name, success in results if not success]
    
    total_success = len(successful_experiments)
    
    # Final summary
    logger.info(f"\n{'='*60}")
    logger.info("PROCESSING COMPLETE")
    logger.info(f"{'='*60}")
    logger.info(f"Total experiments processed: {len(experiment_dirs)}")
    logger.info(f"F(k,t) master files created: {total_success}")
    logger.info(f"Overall success rate: {100.0 * total_success / len(experiment_dirs):.1f}%")
    
    if successful_experiments:
        logger.info(f"Successful experiments: {', '.join(successful_experiments)}")
    
    if failed_experiments:
        logger.warning(f"Failed experiments: {', '.join(failed_experiments)}")
    
    return total_success > 0

def read_fskt_data(file_path):
    """Read F(k,t) data from file."""
    try:
        data = pd.read_csv(file_path, sep=r'\s+', comment='#')
        if data.empty:
            return None
        
        # Convert all data to numeric, coercing errors to NaN
        data = data.apply(pd.to_numeric, errors='coerce')
        
        # Drop any rows with NaN values
        data = data.dropna()
        
        if data.empty:
            logger.warning(f"No valid numeric data found in {file_path}")
            return None
        
        if 'lag_time' in data.columns and 'fskt' in data.columns:
            return {'lag_time': data['lag_time'].values, 'fskt': data['fskt'].values}
        elif data.shape[1] >= 2:
            return {'lag_time': data.iloc[:, 0].values, 'fskt': data.iloc[:, 1].values}
        else:
            return None
    except Exception as e:
        logger.error(f"Error reading F(k,t) data from {file_path}: {e}")
        return None

class FSKTAnalyzer:
    """Analyzes F(k,t) data from all experiments and creates comparison plots."""
    
    def __init__(self, base_dir="./", max_time=None):
        self.base_dir = Path(base_dir).resolve()
        self.max_time = max_time  # For F(k,t) analysis, None = auto-detect
        
        # Define experiment configurations based on current directory structure
        self.experiments = [
            # finiteq experiments
            ExperimentConfig("bussi_langevin_finiteq_coupling_0epos00", "bussi", "langevin", "0.0", "prod", 
                           "Bussi + Langevin (No Coupling, FiniteQ)"),
            ExperimentConfig("bussi_langevin_finiteq_coupling_1eneg04", "bussi", "langevin", "1e-04", "prod", 
                           "Bussi + Langevin (Coupling 1e-04, FiniteQ)"),
            ExperimentConfig("bussi_langevin_finiteq_coupling_2eneg04", "bussi", "langevin", "2e-04", "prod", 
                           "Bussi + Langevin (Coupling 2e-04, FiniteQ)"),
            ExperimentConfig("bussi_langevin_finiteq_coupling_3eneg04", "bussi", "langevin", "3e-04", "prod", 
                           "Bussi + Langevin (Coupling 3e-04, FiniteQ)"),
            ExperimentConfig("bussi_langevin_finiteq_coupling_4eneg04", "bussi", "langevin", "4e-04", "prod", 
                           "Bussi + Langevin (Coupling 4e-04, FiniteQ)"),
            ExperimentConfig("bussi_langevin_finiteq_coupling_5eneg04", "bussi", "langevin", "5e-04", "prod", 
                           "Bussi + Langevin (Coupling 5e-04, FiniteQ)"),
            ExperimentConfig("bussi_langevin_finiteq_coupling_1eneg03", "bussi", "langevin", "1e-03", "prod", 
                           "Bussi + Langevin (Coupling 1e-03, FiniteQ)"),
            # no_finiteq experiments
            ExperimentConfig("bussi_langevin_no_finiteq_coupling_0epos00", "bussi", "langevin", "0.0", "prod", 
                           "Bussi + Langevin (No Coupling)"),
            ExperimentConfig("bussi_langevin_no_finiteq_coupling_1eneg04", "bussi", "langevin", "1e-04", "prod", 
                           "Bussi + Langevin (Coupling 1e-04)"),
            ExperimentConfig("bussi_langevin_no_finiteq_coupling_2eneg04", "bussi", "langevin", "2e-04", "prod", 
                           "Bussi + Langevin (Coupling 2e-04)"),
            ExperimentConfig("bussi_langevin_no_finiteq_coupling_3eneg04", "bussi", "langevin", "3e-04", "prod", 
                           "Bussi + Langevin (Coupling 3e-04)"),
            ExperimentConfig("bussi_langevin_no_finiteq_coupling_4eneg04", "bussi", "langevin", "4e-04", "prod", 
                           "Bussi + Langevin (Coupling 4e-04)"),
            ExperimentConfig("bussi_langevin_no_finiteq_coupling_5eneg04", "bussi", "langevin", "5e-04", "prod", 
                           "Bussi + Langevin (Coupling 5e-04)"),
            ExperimentConfig("bussi_langevin_no_finiteq_coupling_1eneg03", "bussi", "langevin", "1e-03", "prod", 
                           "Bussi + Langevin (Coupling 1e-03)"),
        ]

    def auto_detect_time_ranges(self):
        """Auto-detect appropriate time ranges for each experiment individually."""
        experiment_time_ranges = {}
        
        for exp in self.experiments:
            fskt_files = self.find_fskt_files(exp)
            if fskt_files:
                try:
                    fskt_data = read_fskt_data(str(fskt_files[0]))
                    if fskt_data:
                        max_time = fskt_data['lag_time'].max()
                        experiment_time_ranges[exp.name] = max_time
                        logger.info(f"Auto-detected time range for {exp.name}: 0 to {max_time:.3f} ps")
                except:
                    pass
        
        # Use a reasonable default for experiments without data
        default_time = 50.0
        for exp in self.experiments:
            if exp.name not in experiment_time_ranges:
                experiment_time_ranges[exp.name] = default_time
                logger.warning(f"No data found for {exp.name}, using default time range: {default_time} ps")
        
        return experiment_time_ranges

    def find_fskt_files(self, exp):
        """Find master F(k,t) data files for an experiment."""
        exp_dir = self.base_dir / exp.name
        
        if not exp_dir.exists():
            return []
        
        # Look for master files with any k-value pattern
        fskt_files = list(exp_dir.glob("master_fskt_k*.txt"))
        valid_files = [f for f in fskt_files if f.stat().st_size > 0]
        return sorted(valid_files)

    def create_fskt_comparison(self, all_data):
        """Create F(k,t) comparison plot focusing on different coupling strengths."""
        logger.info("Generating F(k,t) coupling strength comparison plot...")
        
        # Group experiments by thermostat combination
        thermostat_groups = {
            'Bussi + Bussi': [],
            'Bussi + Langevin': []
        }
        
        for exp in self.experiments:
            data = all_data.get(exp.name)
            if not data or not data.get('valid', False) or not data['fskt_data']:
                continue
                
            # Group by molecular + cavity thermostat combination
            group_key = f"{exp.molecular_thermostat.title()} + {exp.cavity_thermostat.title()}"
            if group_key in thermostat_groups:
                thermostat_groups[group_key].append((exp, data))
        
        # Create comparison plot
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown']
        linestyles = ['-', '--', '-.', ':', '-', '--']
        
        for panel_idx, (group_name, group_experiments) in enumerate(thermostat_groups.items()):
            if panel_idx >= len(axes):
                break
                
            ax = axes[panel_idx]
            
            has_data = False
            
            # Plot all experiments in this group
            for i, (exp, data) in enumerate(group_experiments):
                fskt_data = data['fskt_data']
                
                if fskt_data:
                    time = fskt_data['lag_time']
                    fskt = fskt_data['fskt']
                    
                    color = colors[i % len(colors)]
                    linestyle = linestyles[i % len(linestyles)]
                    
                    # Create label from coupling strength
                    if exp.coupling_strength == "no_coupling":
                        label = "No Coupling (Bussi only)"
                    elif exp.coupling_strength == "0.0":
                        label = "No Coupling"
                    else:
                        label = f"Coupling {exp.coupling_strength}"
                    
                    ax.plot(time, fskt, color=color, linestyle=linestyle, 
                           alpha=0.8, linewidth=2.5, label=label)
                    has_data = True
            
            if not has_data:
                ax.text(0.5, 0.5, f'No data available\nfor {group_name}', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=12)
            
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('F(k,t)')
            ax.set_ylim(bottom=0.0)
            ax.set_title(f'{group_name}', fontweight='bold', fontsize=12)
            ax.grid(True, alpha=0.3)
            
            if has_data:
                ax.legend(fontsize=10)
        
        plt.suptitle('F(k,t) Comparison: Different Coupling Strengths\nBy Thermostat Combination', 
                     fontsize=16, fontweight='bold')
        plt.tight_layout()
        output_file = self.base_dir / 'fskt_coupling_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved F(k,t) coupling comparison: {output_file}")
        plt.close()

    def create_overall_coupling_comparison(self, all_data):
        """Create overall comparison plot of all coupling strengths."""
        logger.info("Generating overall coupling strength comparison plot...")
        
        # Get all valid experiments
        valid_experiments = []
        
        for exp in self.experiments:
            data = all_data.get(exp.name)
            if not data or not data.get('valid', False) or not data['fskt_data']:
                continue
            valid_experiments.append((exp, data))
        
        if not valid_experiments:
            logger.warning("No valid experiments found for overall comparison")
            return
        
        # Create single plot with all experiments
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(valid_experiments)))
        
        for i, (exp, data) in enumerate(valid_experiments):
            fskt_data = data['fskt_data']
            if fskt_data:
                time = fskt_data['lag_time']
                fskt = fskt_data['fskt']
                
                ax.plot(time, fskt, color=colors[i], linewidth=2.5, 
                        alpha=0.8, label=exp.display_name)
        
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('F(k,t)')
        ax.set_ylim(bottom=0.0)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=11)
        ax.set_title('F(k,t) Overall Comparison: All Coupling Strengths', 
                     fontweight='bold', fontsize=14)
        
        plt.tight_layout()
        
        output_file = self.base_dir / 'fskt_overall_coupling_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved overall coupling comparison: {output_file}")
        plt.close()

    def create_detailed_fskt_plot(self, all_data):
        """Create detailed F(k,t) plot showing all experiments."""
        logger.info("Generating detailed F(k,t) plot...")
        
        # Determine grid size based on number of experiments
        n_experiments = len(self.experiments)
        if n_experiments <= 4:
            rows, cols = 2, 2
        else:
            rows = int(np.ceil(n_experiments / 3))
            cols = 3
        
        fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 4*rows))
        if n_experiments == 1:
            axes = [axes]
        elif rows == 1 or cols == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        for i, exp in enumerate(self.experiments):
            if i >= len(axes):
                break
                
            ax = axes[i]
            data = all_data.get(exp.name)
            
            if not data or not data.get('valid', False):
                error_msg = data.get('error', 'Unknown error') if data else 'No data'
                ax.text(0.5, 0.5, f'No data available\nfor {exp.display_name}\n({error_msg})', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=10)
                ax.set_title(exp.display_name, fontweight='bold')
                continue
            
            fskt_data = data['fskt_data']
            if not fskt_data:
                ax.text(0.5, 0.5, f'No F(k,t) data\nfor {exp.display_name}', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=10)
                ax.set_title(exp.display_name, fontweight='bold')
                continue
            
            time = fskt_data['lag_time']
            fskt = fskt_data['fskt']
            
            # Use different colors for different coupling strengths
            if exp.coupling_strength == "no_coupling":
                color = 'blue'
            elif exp.coupling_strength == "0.0":
                color = 'blue'  # Same as no_coupling
            elif exp.coupling_strength == "1e-03":
                color = 'red'
            elif exp.coupling_strength == "1e-04":
                color = 'green'
            elif exp.coupling_strength == "5e-04":
                color = 'orange'
            else:
                color = 'purple'
                
            ax.plot(time, fskt, color=color, linewidth=2, alpha=0.8)
            
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('F(k,t)')
            ax.set_ylim(bottom=0.0)
            ax.set_title(exp.display_name, fontweight='bold', fontsize=10)
            ax.grid(True, alpha=0.3)
            
            # Add statistics
            if len(fskt) > 10:
                final_value = fskt[-1]
                min_value = np.min(fskt)
                max_value = np.max(fskt)
                
                stats_text = f"Final: {final_value:.4f}\nMin: {min_value:.4f}\nMax: {max_value:.4f}"
                ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                       verticalalignment='top', fontsize=8,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Remove unused axes
        for i in range(len(self.experiments), len(axes)):
            fig.delaxes(axes[i])
        
        plt.suptitle('Detailed F(k,t) Analysis - All Experiments', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        output_file = self.base_dir / "detailed_fskt_analysis.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved detailed F(k,t) analysis: {output_file}")
        plt.close()

    def create_single_experiment_coupling_comparison(self, all_data, experiment_type):
        """Create coupling strength comparison for a specific experiment type."""
        logger.info(f"Generating single experiment coupling comparison for: {experiment_type}")
        
        # Find experiments that match the specified type
        matching_experiments = []
        
        for exp in self.experiments:
            if experiment_type.lower() in exp.name.lower():
                data = all_data.get(exp.name)
                if data and data.get('valid', False) and data['fskt_data']:
                    matching_experiments.append((exp, data))
        
        if not matching_experiments:
            logger.warning(f"No valid data found for experiment type: {experiment_type}")
            return False
        
        # Create the comparison plot
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        colors = ['blue', 'red', 'green', 'orange', 'purple']
        has_data = False
        
        for i, (exp, data) in enumerate(matching_experiments):
            fskt_data = data['fskt_data']
            if fskt_data:
                time = fskt_data['lag_time']
                fskt = fskt_data['fskt']
                
                color = colors[i % len(colors)]
                
                # Create label from coupling strength
                if exp.coupling_strength == "no_coupling":
                    label = "No Coupling"
                elif exp.coupling_strength == "0.0":
                    label = "No Coupling"
                else:
                    label = f"Coupling {exp.coupling_strength}"
                
                ax.plot(time, fskt, color=color, linewidth=3, alpha=0.8, 
                       label=label, marker='o', markersize=2)
                has_data = True
        
        if not has_data:
            ax.text(0.5, 0.5, f'No valid F(k,t) data found\nfor {experiment_type}', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=14)
            logger.warning(f"No valid F(k,t) data found for {experiment_type}")
            return False
        
        # Customize the plot
        ax.set_xlabel('Time (ps)', fontsize=14)
        ax.set_ylabel('F(k,t)', fontsize=14)
        ax.set_ylim(bottom=0.0)
        
        title = f'{experiment_type.replace("_", " ").title()}\nCoupling Strength Comparison'
        ax.set_title(title, fontweight='bold', fontsize=16)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=12, loc='best')
        
        plt.tight_layout()
        
        # Save the plot
        safe_name = experiment_type.replace('_', '-')
        output_file = self.base_dir / f'fskt_single_{safe_name}_coupling_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved single experiment coupling comparison: {output_file}")
        plt.close()
        
        return True

    def run_analysis(self, single_experiment=None):
        """Run the complete F(k,t) analysis."""
        logger.info("Starting F(k,t) analysis...")
        
        # Auto-detect time ranges for each experiment if not specified
        if self.max_time is None:
            experiment_time_ranges = self.auto_detect_time_ranges()
            logger.info("Using individual time ranges for each experiment")
        else:
            # Use the same time range for all experiments
            experiment_time_ranges = {exp.name: self.max_time for exp in self.experiments}
            logger.info(f"Using fixed time range for all experiments: 0 to {self.max_time:.3f} ps")
        
        # Process all experiments
        all_data = {}
        valid_count = 0
        
        # Add detailed logging for each experiment
        logger.info("Checking experiment status:")
        for exp in self.experiments:
            exp_dir = self.base_dir / exp.name
            fskt_files = self.find_fskt_files(exp)
            
            if not exp_dir.exists():
                logger.warning(f"  {exp.name}: Directory does not exist")
            elif not fskt_files:
                logger.warning(f"  {exp.name}: No master F(k,t) files found")
            else:
                time_range = experiment_time_ranges.get(exp.name, 50.0)
                logger.info(f"  {exp.name}: Found {len(fskt_files)} F(k,t) files, time range: 0 to {time_range:.3f} ps")
        
        for exp in self.experiments:
            try:
                logger.info(f"Processing experiment: {exp.name}")
                
                fskt_files = self.find_fskt_files(exp)
                
                if not fskt_files:
                    all_data[exp.name] = {'valid': False, 'error': 'No F(k,t) files found'}
                    continue
                
                # Read F(k,t) data from first file (assuming k1 for primary analysis)
                raw_fskt_data = read_fskt_data(str(fskt_files[0]))
                if not raw_fskt_data:
                    all_data[exp.name] = {'valid': False, 'error': 'Failed to read F(k,t) data'}
                    continue
                
                # Apply experiment-specific time filter
                fskt_data = raw_fskt_data
                exp_max_time = experiment_time_ranges.get(exp.name)
                if exp_max_time is not None:
                    fskt_time = raw_fskt_data['lag_time']
                    fskt_values = raw_fskt_data['fskt']
                    
                    # Ensure data types are correct
                    if not np.issubdtype(fskt_time.dtype, np.number):
                        logger.error(f"  {exp.name}: lag_time data is not numeric, dtype: {fskt_time.dtype}")
                        logger.error(f"  Sample values: {fskt_time[:5]}")
                        all_data[exp.name] = {'valid': False, 'error': 'Non-numeric time data'}
                        continue
                    
                    if not np.issubdtype(fskt_values.dtype, np.number):
                        logger.error(f"  {exp.name}: fskt data is not numeric, dtype: {fskt_values.dtype}")
                        logger.error(f"  Sample values: {fskt_values[:5]}")
                        all_data[exp.name] = {'valid': False, 'error': 'Non-numeric fskt data'}
                        continue
                    
                    fskt_mask = fskt_time <= exp_max_time
                    
                    # Log how much data we're using
                    total_points = len(fskt_time)
                    used_points = np.sum(fskt_mask)
                    max_available_time = fskt_time.max() if len(fskt_time) > 0 else 0
                    
                    logger.info(f"  {exp.name}: Using {used_points}/{total_points} data points, "
                              f"available time range: 0 to {max_available_time:.3f} ps, "
                              f"using up to: {exp_max_time:.3f} ps")
                    
                    fskt_data = {
                        'lag_time': fskt_time[fskt_mask],
                        'fskt': fskt_values[fskt_mask]
                    }
                
                all_data[exp.name] = {
                    'fskt_data': fskt_data,
                    'valid': True,
                    'error': None
                }
                valid_count += 1
                logger.info(f"Successfully processed {exp.name}")
                
            except Exception as e:
                logger.error(f"Error processing {exp.name}: {e}")
                all_data[exp.name] = {'valid': False, 'error': str(e)}
        
        if valid_count == 0:
            logger.error("No experiment data could be processed!")
            return False
        
        logger.info(f"Successfully processed {valid_count}/{len(self.experiments)} experiments")
        
        # Create comparison plots
        if single_experiment:
            # Focus on single experiment comparison
            logger.info(f"Creating single experiment comparison for: {single_experiment}")
            single_success = self.create_single_experiment_coupling_comparison(all_data, single_experiment)
            if not single_success:
                logger.warning(f"Failed to create single experiment comparison for {single_experiment}")
        else:
            # Create all comparison plots
            self.create_fskt_comparison(all_data)
            self.create_overall_coupling_comparison(all_data)
            self.create_detailed_fskt_plot(all_data)
        
        logger.info("F(k,t) analysis complete!")
        return True

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Process F(k,t) files and generate analysis plots')
    parser.add_argument('--base_dir', default='./', help='Base directory containing experiment folders')
    parser.add_argument('--max_time', type=float, default=None,
                       help='Maximum time for F(k,t) data (ps) - if not specified, auto-detect from data')
    parser.add_argument('--fskt_dt', type=float, default=1.0,
                       help='Timestep for F(k,t) data interpolation (ps)')
    parser.add_argument('--skip_processing', action='store_true',
                       help='Skip master file processing (use existing files)')
    parser.add_argument('--single_experiment', type=str, default=None,
                       help='Focus on single experiment comparison (e.g., "bussi_langevin", "bussi_bussi")')
    parser.add_argument('--processes', type=int, default=None,
                       help='Number of parallel processes to use (None = auto-detect, 1 = sequential)')
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    logger.info("Starting F(k,t) processing and analysis...")
    logger.info(f"Base directory: {Path(args.base_dir).resolve()}")
    logger.info(f"F(k,t) dt: {args.fskt_dt} ps")
    
    if args.max_time is not None:
        logger.info(f"Max F(k,t) time: {args.max_time} ps (fixed)")
    else:
        logger.info("Max F(k,t) time: auto-detect from data")
    
    if args.single_experiment:
        logger.info(f"Single experiment mode: {args.single_experiment}")
    
    # Log parallelization settings
    if args.processes is not None:
        if args.processes == 1:
            logger.info("Parallelization: Sequential mode (1 process)")
        else:
            logger.info(f"Parallelization: {args.processes} processes")
    else:
        available_cpus = mp.cpu_count()
        logger.info(f"Parallelization: Auto-detect (available CPUs: {available_cpus})")
    
    success = True
    
    # Step 1: Process experiments to create master F(k,t) files (unless skipped)
    if not args.skip_processing:
        logger.info("\n" + "="*60)
        logger.info("STEP 1: PROCESSING EXPERIMENTS TO CREATE MASTER F(k,t) FILES")
        logger.info("="*60)
        
        success = process_all_experiments(
            base_dir=args.base_dir,
            max_time=args.max_time,
            dt=args.fskt_dt,
            num_processes=args.processes
        )
        
        if not success:
            logger.error("Failed to process experiments!")
            return 1
    else:
        logger.info("Skipping master file processing (using existing files)")
    
    # Step 2: Run F(k,t) analysis and generate plots
    logger.info("\n" + "="*60)
    logger.info("STEP 2: RUNNING F(k,t) ANALYSIS AND GENERATING PLOTS")
    logger.info("="*60)
    
    analyzer = FSKTAnalyzer(
        base_dir=args.base_dir,
        max_time=args.max_time
    )
    
    analysis_success = analyzer.run_analysis(single_experiment=args.single_experiment)
    
    total_time = time.time() - start_time
    
    if analysis_success:
        logger.info(f"\n" + "="*60)
        logger.info("F(k,t) PROCESSING AND ANALYSIS COMPLETE!")
        logger.info("="*60)
        logger.info(f"Total execution time: {total_time:.2f} seconds")
        logger.info("Generated files:")
        logger.info("  - master_fskt_k*.txt (in each experiment folder)")
        
        if args.single_experiment:
            logger.info(f"  - fskt_single_{args.single_experiment.replace('_', '-')}_coupling_comparison.png (single experiment coupling comparison)")
        else:
            logger.info("  - fskt_coupling_comparison.png (side-by-side coupling comparison by thermostat)")
            logger.info("  - fskt_overall_coupling_comparison.png (overall coupling strength comparison)")
            logger.info("  - detailed_fskt_analysis.png (individual experiment plots)")
        return 0
    else:
        logger.error("F(k,t) analysis failed!")
        return 1

if __name__ == "__main__":
    # Ensure multiprocessing works properly on all platforms
    mp.set_start_method('spawn', force=True)
    exit(main()) 
