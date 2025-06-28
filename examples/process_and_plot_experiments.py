#!/usr/bin/env python3
"""
Combined Experiment Processing and Analysis Script for Cavity MD Simulations

This script combines the functionality of:
1. process_experiment_master_files.py - generates master files from individual replicas
2. plot_multi_experiment_analysis.py - creates comprehensive analysis plots

It processes all experiment directories, creates master files using configurable timestep
interpolations, and then generates comparison plots and analysis.

Key Features:
- Configurable timestep resolution for energy/temperature vs F(k,t) data
- Processes all 12 thermostat experiments
- Generates master files by averaging across all replicas
- Creates comprehensive comparison plots and statistics
- Uses standardized reservoir energy column names (molecular_reservoir_energy, cavity_reservoir_energy)
  that combine values from all thermostat types (Langevin + Bussi)

Usage:
    python process_and_plot_experiments.py --base_dir ./
    python process_and_plot_experiments.py --energy_dt 0.01 --fskt_dt 1.0 --target_temp 100
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
import subprocess
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

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
    finite_q: bool
    job_name: str  # 'prod' or 'finq'
    display_name: str

def process_energy_files(exp_dir, job_name, max_time=100.0, dt=0.01):
    """
    Process energy contribution files using configurable timestep and proper cumulative averaging.
    
    Args:
        exp_dir: Path to experiment directory
        job_name: Job name ('prod' or 'finq')
        max_time: Maximum time in ps for energy data
        dt: Timestep in ps for interpolation
    """
    exp_path = Path(exp_dir)
    
    # Find all energy contribution files for this experiment
    energy_pattern = f"{job_name}-*_energy_contributions.txt"
    energy_files = list(exp_path.glob(energy_pattern))
    energy_files = natsort.natsorted([str(f) for f in energy_files])
    
    if not energy_files:
        logger.warning(f"No energy files found in {exp_dir} with pattern {energy_pattern}")
        return False
    
    logger.info(f"Processing {len(energy_files)} energy files in {exp_dir}")
    
    # Set fixed linspace with configurable timestep
    num_points = int(max_time / dt) + 1
    uniform_times = np.linspace(0, max_time, num_points)
    logger.info(f"Using fixed time grid: 0 to {max_time} ps with {len(uniform_times)} points (dt = {dt} ps)")
    
    # Initialize storage for cumulative data
    cumulative_data = None
    sample_counts = np.zeros(len(uniform_times), dtype=int)
    processed_count = 0
    
    for energy_file in energy_files:
        try:
            # Read energy data, skipping comment lines
            data = pd.read_csv(energy_file, sep=r'\s+', comment='#', header=None)
            
            # Skip header rows if present
            if not data.empty:
                first_row = data.iloc[0]
                if any(isinstance(val, str) and not val.replace('.', '').replace('-', '').replace('e', '').replace('+', '').isdigit() for val in first_row):
                    data = data.iloc[1:]
            
            # Convert to numeric and drop invalid rows
            data = data.apply(pd.to_numeric, errors='coerce').dropna()
            
            if data.empty:
                logger.warning(f"No valid numeric data in {energy_file}")
                continue
            
            # Define expected columns based on data width
            num_cols = len(data.columns)
            if num_cols >= 15:  # Lowered threshold from 16 to 15
                expected_columns = [
                    'time(ps)', 'timestep', 'harmonic_energy', 'lj_energy',
                    'coulomb_short_energy', 'coulomb_long_energy', 'cavity_harmonic_energy',
                    'cavity_coupling_energy', 'cavity_dipole_self_energy', 'cavity_total_potential_energy',
                    'molecular_kinetic_energy', 'cavity_mode_kinetic_energy', 'molecular_reservoir_energy',
                    'cavity_reservoir_energy'
                ]
                
                # Add additional columns that may be present depending on thermostat configuration
                additional_columns = [
                    'molecular_mttk_thermostat_energy', 'cavity_mttk_thermostat_energy',
                    'molecular_bussi_translational_energy', 'molecular_bussi_rotational_energy',
                    'cavity_bussi_translational_energy', 'cavity_bussi_rotational_energy',
                    'total_potential_energy', 'total_kinetic_energy', 'universe_total_energy'
                ]
                
                # Extend expected columns up to the actual number of columns
                full_expected = expected_columns + additional_columns
                if num_cols >= len(expected_columns):
                    data.columns = full_expected[:num_cols]
                else:
                    data.columns = expected_columns[:num_cols]
            else:
                # Fallback for files with fewer columns - still try to assign meaningful names
                base_columns = ['time(ps)', 'timestep']
                if num_cols <= 2:
                    data.columns = base_columns[:num_cols]
                else:
                    # Use generic names but try to identify key columns
                    data.columns = [f'col_{i}' for i in range(num_cols)]
                    # Try to rename first column to time
                    if num_cols > 0:
                        data.columns = ['time(ps)'] + [f'col_{i}' for i in range(1, num_cols)]
            
            # Calculate derived energy columns if components are available
            if all(col in data.columns for col in ['harmonic_energy', 'lj_energy', 'coulomb_short_energy', 'coulomb_long_energy', 'cavity_total_potential_energy']):
                molecular_potential = data['harmonic_energy'] + data['lj_energy'] + data['coulomb_short_energy'] + data['coulomb_long_energy']
                data['total_potential_energy'] = molecular_potential + data['cavity_total_potential_energy']
            
            if all(col in data.columns for col in ['molecular_kinetic_energy', 'cavity_mode_kinetic_energy']):
                data['total_kinetic_energy'] = data['molecular_kinetic_energy'] + data['cavity_mode_kinetic_energy']
            
            # Calculate universe_total_energy if both potential and kinetic are available
            if 'total_potential_energy' in data.columns and 'total_kinetic_energy' in data.columns:
                data['universe_total_energy'] = data['total_potential_energy'] + data['total_kinetic_energy']
            
            # Calculate molecular and cavity total energies for subsystem analysis
            if all(col in data.columns for col in ['harmonic_energy', 'lj_energy', 'coulomb_short_energy', 'coulomb_long_energy', 'molecular_kinetic_energy']):
                data['molecular_total_energy'] = (data['harmonic_energy'] + data['lj_energy'] + 
                                                data['coulomb_short_energy'] + data['coulomb_long_energy'] + 
                                                data['molecular_kinetic_energy'])
            
            if all(col in data.columns for col in ['cavity_total_potential_energy', 'cavity_mode_kinetic_energy']):
                data['cavity_total_energy'] = data['cavity_total_potential_energy'] + data['cavity_mode_kinetic_energy']
            
            # Get time column and data
            time_col = 'time(ps)' if 'time(ps)' in data.columns else data.columns[0]
            times = data[time_col].values
            
            # Filter to valid time range
            valid_mask = (times >= 0) & (times <= max_time)
            if np.sum(valid_mask) < 3:
                logger.warning(f"Insufficient valid time points in {energy_file}")
                continue
            
            times = times[valid_mask]
            data = data.iloc[valid_mask]
            
            # Initialize cumulative data structure on first successful file
            if cumulative_data is None:
                uniform_data = pd.DataFrame({time_col: uniform_times})
                for col in data.columns:
                    if col != time_col:
                        uniform_data[col] = np.zeros(len(uniform_times))
                cumulative_data = uniform_data.copy()
                logger.info(f"Initialized cumulative data structure with {len(uniform_times)} time points")
            
            # Process each column with interpolation
            file_valid_mask = np.zeros(len(uniform_times), dtype=bool)
            
            for col in data.columns:
                if col != time_col:
                    values = data[col].values
                    
                    # Remove duplicate time values before interpolation
                    unique_times, unique_values = remove_duplicate_times(times, values)
                    
                    if len(unique_times) < 2:
                        logger.warning(f"Insufficient unique time points for column {col} in {energy_file}")
                        continue
                    
                    # Create interpolation function
                    if len(unique_times) > 3:
                        interp_func = interp1d(unique_times, unique_values, kind='cubic', 
                                             bounds_error=False, fill_value=np.nan)
                    else:
                        interp_func = interp1d(unique_times, unique_values, kind='linear',
                                             bounds_error=False, fill_value=np.nan)
                    
                    # Evaluate at uniform time points
                    interpolated_values = interp_func(uniform_times)
                    column_valid_mask = ~np.isnan(interpolated_values)
                    file_valid_mask |= column_valid_mask
                    
                    # Update cumulative average
                    for i in range(len(uniform_times)):
                        if column_valid_mask[i]:
                            current_count = sample_counts[i] + 1
                            
                            if current_count == 1:
                                cumulative_data.iloc[i, cumulative_data.columns.get_loc(col)] = interpolated_values[i]
                            else:
                                old_avg = cumulative_data.iloc[i, cumulative_data.columns.get_loc(col)]
                                new_val = interpolated_values[i]
                                new_avg = old_avg + (new_val - old_avg) / current_count
                                cumulative_data.iloc[i, cumulative_data.columns.get_loc(col)] = new_avg
            
            # Update sample counts
            sample_counts[file_valid_mask] += 1
            processed_count += 1
            
            if processed_count % 10 == 0:
                logger.info(f"Processed {processed_count}/{len(energy_files)} energy files")
                
        except Exception as e:
            logger.error(f"Error processing {energy_file}: {e}")
            continue
    
    # Save results if we have data
    if cumulative_data is not None and processed_count > 0:
        master_file = exp_path / "master_energy_contributions.txt"
        sample_count_file = exp_path / "master_energy_sample_counts.txt"
        
        # Only keep time points where we have at least one sample
        valid_points = sample_counts > 0
        clean_data = cumulative_data[valid_points].copy()
        clean_sample_counts = sample_counts[valid_points]
        
        # Save master data
        with open(master_file, 'w') as f:
            f.write(f"# Master energy contributions file for {exp_path.name}\n")
            f.write(f"# Cumulative average across {processed_count} replicas\n")
            f.write(f"# Fixed time grid: 0 to {max_time} ps with {len(uniform_times)} points (dt = {dt} ps)\n")
            f.write(f"# Generated by process_and_plot_experiments.py\n")
            f.write(f"# Uses standardized reservoir energy columns (molecular_reservoir_energy, cavity_reservoir_energy)\n")
            f.write(f"# that combine values from all thermostat types (Langevin + Bussi)\n")
            f.write(f"# Derived energy columns calculated: total_potential_energy, total_kinetic_energy, universe_total_energy\n")
            f.write(f"# molecular_total_energy, cavity_total_energy\n")
        
        clean_data.to_csv(master_file, sep='\t', index=False, mode='a', header=True, na_rep='nan')
        np.savetxt(sample_count_file, clean_sample_counts, fmt='%d')
        
        logger.info(f"Created master energy file: {master_file} ({processed_count} replicas)")
        logger.info(f"  Valid time points: {len(clean_data)}/{len(uniform_times)}")
        logger.info(f"  Time range: {clean_data[time_col].min():.3f} - {clean_data[time_col].max():.3f} ps")
        return True
    
    return False

def process_fskt_files(exp_dir, job_name, max_time=50.0, dt=1.0):
    """
    Process F(k,t) files using configurable timestep and proper cumulative averaging.
    
    Args:
        exp_dir: Path to experiment directory
        job_name: Job name ('prod' or 'finq')
        max_time: Maximum time in ps for F(k,t) data
        dt: Timestep in ps for interpolation
    """
    exp_path = Path(exp_dir)
    
    # Find all F(k,t) files for this experiment
    fskt_pattern = f"{job_name}-*_fskt_k*.txt"
    fskt_files = list(exp_path.glob(fskt_pattern))
    
    if not fskt_files:
        logger.warning(f"No F(k,t) files found in {exp_dir} with pattern {fskt_pattern}")
        return False
    
    # Group files by k-value and reference
    fskt_groups = {}
    for fskt_file in fskt_files:
        filename = fskt_file.name
        parts = filename.split('_fskt_')[1].replace('.txt', '')
        
        if parts not in fskt_groups:
            fskt_groups[parts] = []
        fskt_groups[parts].append(str(fskt_file))
    
    logger.info(f"Processing {len(fskt_groups)} F(k,t) groups in {exp_dir}")
    success_count = 0
    
    for parts, file_list in fskt_groups.items():
        try:
            # Sort files naturally
            file_list = natsort.natsorted(file_list)
            
            # Set fixed linspace with configurable timestep
            num_points = int(max_time / dt) + 1
            uniform_times = np.linspace(0, max_time, num_points)
            logger.info(f"F(k,t) group {parts}: Using fixed time grid: 0 to {max_time} ps with {len(uniform_times)} points (dt = {dt} ps)")
            
            # Initialize storage for cumulative data
            cumulative_fskt = np.zeros(len(uniform_times))
            sample_counts = np.zeros(len(uniform_times), dtype=int)
            processed_count = 0
            
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
                    
                    # Use only first two columns: time and F(k,t)
                    data = data.iloc[:, :2]
                    data.columns = ['lag_time', 'fskt']
                    
                    times = data['lag_time'].values
                    fskt_vals = data['fskt'].values
                    
                    # Filter to valid time range
                    valid_mask = (times >= 0) & (times <= max_time)
                    if np.sum(valid_mask) < 3:
                        continue
                    
                    times = times[valid_mask]
                    fskt_vals = fskt_vals[valid_mask]
                    
                    # Remove duplicate time values before interpolation
                    unique_times, unique_fskt = remove_duplicate_times(times, fskt_vals)
                    
                    if len(unique_times) < 2:
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
                    
                    # Update sample counts
                    sample_counts[valid_mask] += 1
                    
                    # Update cumulative average
                    for i in range(len(uniform_times)):
                        if valid_mask[i]:
                            if sample_counts[i] == 1:
                                cumulative_fskt[i] = interpolated_fskt[i]
                            else:
                                old_avg = cumulative_fskt[i]
                                new_val = interpolated_fskt[i]
                                n = sample_counts[i]
                                new_avg = old_avg + (new_val - old_avg) / n
                                cumulative_fskt[i] = new_avg
                    
                    processed_count += 1
                        
                except Exception as e:
                    logger.warning(f"Error processing {fskt_file}: {e}")
                    continue
            
            # Save results if we have data
            if processed_count > 0:
                master_file = exp_path / f"master_fskt_{parts}.txt"
                sample_count_file = exp_path / f"master_fskt_{parts}_sample_counts.txt"
                
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
                
                # Save master data
                with open(master_file, 'w') as f:
                    f.write(f"# Master F(k,t) file for {exp_path.name} - {parts}\n")
                    f.write(f"# Cumulative average across {processed_count} replicas\n")
                    f.write(f"# Fixed time grid: 0 to {max_time} ps with {len(uniform_times)} points (dt = {dt} ps)\n")
                    f.write(f"# Generated by process_and_plot_experiments.py\n")
                    f.write(f"# Columns: lag_time fskt\n")
                
                output_data.to_csv(master_file, sep='\t', index=False, mode='a', header=True, na_rep='nan')
                np.savetxt(sample_count_file, clean_sample_counts, fmt='%d')
                
                logger.info(f"Created master F(k,t) file: {master_file} ({processed_count} replicas)")
                logger.info(f"  Valid time points: {len(output_data)}/{len(uniform_times)}")
                logger.info(f"  Time range: {clean_times.min():.3f} - {clean_times.max():.3f} ps")
                success_count += 1
                
        except Exception as e:
            logger.error(f"Error processing F(k,t) group {parts}: {e}")
            continue
    
    return success_count > 0

def process_all_experiments(base_dir="./", max_energy_time=100.0, max_fskt_time=50.0, 
                           energy_dt=0.01, fskt_dt=1.0):
    """
    Process all experiment directories to generate master files.
    
    Args:
        base_dir: Base directory containing experiment folders
        max_energy_time: Maximum time for energy data (ps)
        max_fskt_time: Maximum time for F(k,t) data (ps)
        energy_dt: Timestep for energy data interpolation (ps)
        fskt_dt: Timestep for F(k,t) data interpolation (ps)
    """
    base_path = Path(base_dir)
    
    # Define expected experiment patterns
    experiment_patterns = [
        "*_langevin_finiteq",
        "*_langevin_no_finiteq", 
        "*_bussi_finiteq",
        "*_bussi_no_finiteq"
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
    
    total_energy_success = 0
    total_fskt_success = 0
    
    for exp_dir in experiment_dirs:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing experiment: {exp_dir.name}")
        logger.info(f"{'='*60}")
        
        # Determine job name based on experiment type
        job_name = "prod"
        
        logger.info(f"Using job name: {job_name}")
        
        # Process energy files
        logger.info("Processing energy contribution files...")
        energy_success = process_energy_files(exp_dir, job_name, max_energy_time, energy_dt)
        if energy_success:
            total_energy_success += 1
            logger.info("✅ Energy files processed successfully")
        else:
            logger.warning("❌ Energy file processing failed")
        
        # Process F(k,t) files
        logger.info("Processing F(k,t) files...")
        fskt_success = process_fskt_files(exp_dir, job_name, max_fskt_time, fskt_dt)
        if fskt_success:
            total_fskt_success += 1
            logger.info("✅ F(k,t) files processed successfully")
        else:
            logger.warning("❌ F(k,t) file processing failed")
    
    # Final summary
    logger.info(f"\n{'='*60}")
    logger.info("PROCESSING COMPLETE")
    logger.info(f"{'='*60}")
    logger.info(f"Total experiments processed: {len(experiment_dirs)}")
    logger.info(f"Energy master files created: {total_energy_success}")
    logger.info(f"F(k,t) master files created: {total_fskt_success}")
    logger.info(f"Overall success rate: {100.0 * min(total_energy_success, total_fskt_success) / len(experiment_dirs):.1f}%")
    
    return total_energy_success > 0 and total_fskt_success > 0

# Analysis and plotting functions (simplified versions from plot_multi_experiment_analysis.py)

def read_energy_data(file_path):
    """Read energy data from file."""
    try:
        df = pd.read_csv(file_path, sep='\s+', comment='#')
        return df
    except Exception as e:
        logger.error(f"Error reading energy data from {file_path}: {e}")
        return None

def calculate_derived_energies(df):
    """Calculate derived energy quantities."""
    derived = {}
    
    # Use actual column names from the dataframe
    for col in df.columns:
        if col != df.columns[0]:  # Skip time column
            derived[col] = df[col]
    
    # Calculate derived quantities if possible
    if 'total_potential_energy' in df.columns and 'total_kinetic_energy' in df.columns:
        derived['universe_total_energy'] = df['total_potential_energy'] + df['total_kinetic_energy']
        derived['molecular_total_energy'] = derived.get('molecular_kinetic_energy', 0) + \
                                          (derived.get('harmonic_energy', 0) + 
                                           derived.get('lj_energy', 0) + 
                                           derived.get('coulomb_short_energy', 0) + 
                                           derived.get('coulomb_long_energy', 0))
        derived['cavity_total_energy'] = derived.get('cavity_mode_kinetic_energy', 0) + \
                                        derived.get('cavity_total_potential_energy', 0)
    
    return derived

def read_fskt_data(file_path):
    """Read F(k,t) data from file."""
    try:
        data = pd.read_csv(file_path, sep='\s+', comment='#')
        if data.empty:
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

class MultiExperimentAnalyzer:
    """Analyzes data from all experiments and creates comparison plots."""
    
    def __init__(self, base_dir="./", target_temp=100.0, max_energy_time=None, max_fskt_time=None):
        self.base_dir = Path(base_dir).resolve()
        self.target_temp = target_temp
        self.max_energy_time = max_energy_time  # For energy and temperature analysis, None = auto-detect
        self.max_fskt_time = max_fskt_time      # For F(k,t) analysis, None = auto-detect
        
        # Define all 12 experiment configurations
        self.experiments = [
            ExperimentConfig("bussi_langevin_finiteq", "bussi", "langevin", True, "prod", 
                           "Bussi + Langevin + Finite-q"),
            ExperimentConfig("bussi_langevin_no_finiteq", "bussi", "langevin", False, "finq", 
                           "Bussi + Langevin + No Finite-q"),
            ExperimentConfig("langevin_langevin_finiteq", "langevin", "langevin", True, "prod", 
                           "Langevin + Langevin + Finite-q"),
            ExperimentConfig("langevin_langevin_no_finiteq", "langevin", "langevin", False, "finq", 
                           "Langevin + Langevin + No Finite-q"),
            ExperimentConfig("none_langevin_finiteq", "none", "langevin", True, "prod", 
                           "None + Langevin + Finite-q"),
            ExperimentConfig("none_langevin_no_finiteq", "none", "langevin", False, "finq", 
                           "None + Langevin + No Finite-q"),
            ExperimentConfig("bussi_bussi_finiteq", "bussi", "bussi", True, "prod", 
                           "Bussi + Bussi + Finite-q"),
            ExperimentConfig("bussi_bussi_no_finiteq", "bussi", "bussi", False, "finq", 
                           "Bussi + Bussi + No Finite-q"),
            ExperimentConfig("langevin_bussi_finiteq", "langevin", "bussi", True, "prod", 
                           "Langevin + Bussi + Finite-q"),
            ExperimentConfig("langevin_bussi_no_finiteq", "langevin", "bussi", False, "finq", 
                           "Langevin + Bussi + No Finite-q"),
            ExperimentConfig("none_bussi_finiteq", "none", "bussi", True, "prod", 
                           "None + Bussi + Finite-q"),
            ExperimentConfig("none_bussi_no_finiteq", "none", "bussi", False, "finq", 
                           "None + Bussi + No Finite-q"),
        ]

    def auto_detect_time_ranges(self):
        """Auto-detect appropriate time ranges from available data."""
        energy_max_times = []
        fskt_max_times = []
        
        for exp in self.experiments:
            # Check energy file time range
            files = self.find_master_files(exp)
            if files['energy']:
                try:
                    df = read_energy_data(str(files['energy']))
                    if df is not None and not df.empty:
                        max_time = df.iloc[:, 0].max()
                        energy_max_times.append(max_time)
                except:
                    pass
            
            # Check F(k,t) file time range
            fskt_files = self.find_fskt_files(exp)
            if fskt_files:
                try:
                    fskt_data = read_fskt_data(str(fskt_files[0]))
                    if fskt_data:
                        max_time = fskt_data['lag_time'].max()
                        fskt_max_times.append(max_time)
                except:
                    pass
        
        # Use the minimum of the maximum times to ensure all experiments have data
        if energy_max_times:
            detected_energy_time = min(energy_max_times)
            logger.info(f"Auto-detected energy time range: 0 to {detected_energy_time:.3f} ps")
        else:
            detected_energy_time = 1.0  # fallback
            logger.warning("Could not auto-detect energy time range, using 1.0 ps")
        
        if fskt_max_times:
            detected_fskt_time = min(fskt_max_times)
            logger.info(f"Auto-detected F(k,t) time range: 0 to {detected_fskt_time:.3f} ps")
        else:
            detected_fskt_time = 50.0  # fallback
            logger.warning("Could not auto-detect F(k,t) time range, using 50.0 ps")
        
        return detected_energy_time, detected_fskt_time

    def find_master_files(self, exp):
        """Find master data files for an experiment."""
        exp_dir = self.base_dir / exp.name
        
        files = {
            'energy': exp_dir / "master_energy_contributions.txt",
        }
        
        if not files['energy'].exists():
            files['energy'] = None
        
        return files

    def find_fskt_files(self, exp):
        """Find master F(k,t) data files for an experiment."""
        exp_dir = self.base_dir / exp.name
        
        if not exp_dir.exists():
            return []
        
        fskt_files = list(exp_dir.glob("master_fskt_k*.txt"))
        valid_files = [f for f in fskt_files if f.stat().st_size > 0]
        return sorted(valid_files)

    def process_single_experiment(self, exp):
        """Process a single experiment and return analysis data."""
        logger.info(f"Processing experiment: {exp.name}")
        
        files = self.find_master_files(exp)
        
        if not files['energy']:
            return {'valid': False, 'error': 'Energy file not found'}
        
        try:
            # Read energy data
            df = read_energy_data(str(files['energy']))
            if df is None:
                return {'valid': False, 'error': 'Failed to read energy data'}
            
            # Calculate derived energies
            derived_energies = calculate_derived_energies(df)
            
            # Calculate temperature data from kinetic energies
            temperature_data = {}
            kB_au = 3.166811563e-6  # Boltzmann constant in atomic units
            
            if 'molecular_kinetic_energy' in df.columns:
                N_molecules = 500  # Approximate
                molecular_temp = (2.0 * df['molecular_kinetic_energy']) / (3.0 * N_molecules * kB_au)
                temperature_data['molecular_temperature'] = molecular_temp.values
            
            if 'cavity_mode_kinetic_energy' in df.columns:
                cavity_temp = (2.0 * df['cavity_mode_kinetic_energy']) / (3.0 * kB_au)
                temperature_data['cavity_temperature'] = cavity_temp.values
            
            # Process F(k,t) data
            fskt_files = self.find_fskt_files(exp)
            fskt_data = None
            if fskt_files:
                raw_fskt_data = read_fskt_data(str(fskt_files[0]))
                if raw_fskt_data and self.max_fskt_time is not None:
                    # Apply time filter to F(k,t) data using max_fskt_time
                    fskt_time = raw_fskt_data['lag_time']
                    fskt_values = raw_fskt_data['fskt']
                    fskt_mask = fskt_time <= self.max_fskt_time
                    fskt_data = {
                        'lag_time': fskt_time[fskt_mask],
                        'fskt': fskt_values[fskt_mask]
                    }
                else:
                    fskt_data = raw_fskt_data
            
            # Apply time filter if specified - use max_energy_time for energy and temperature data
            time = df.iloc[:, 0].values
            if self.max_energy_time is not None:
                mask = time <= self.max_energy_time
                time = time[mask]
                for key in derived_energies:
                    if hasattr(derived_energies[key], '__len__') and len(derived_energies[key]) == len(mask):
                        if hasattr(derived_energies[key], 'iloc'):
                            derived_energies[key] = derived_energies[key].iloc[mask].reset_index(drop=True)
                        else:
                            derived_energies[key] = derived_energies[key][mask]
                for key in temperature_data:
                    if len(temperature_data[key]) == len(mask):
                        temperature_data[key] = temperature_data[key][mask]
            
            return {
                'time': time,
                'energy_components': derived_energies,
                'temperature_data': temperature_data,
                'fskt_data': fskt_data,
                'valid': True,
                'error': None
            }
            
        except Exception as e:
            logger.error(f"Error processing {exp.name}: {e}")
            return {'valid': False, 'error': str(e)}

    def create_master_energy_comparison(self, all_data):
        """Create master energy comparison plot."""
        logger.info("Generating master energy comparison plot...")
        
        fig, axes = plt.subplots(4, 3, figsize=(18, 16))
        axes = axes.flatten()
        
        for i, exp in enumerate(self.experiments):
            ax = axes[i]
            data = all_data.get(exp.name)
            
            if not data or not data.get('valid', False):
                error_msg = data.get('error', 'Unknown error') if data else 'No data'
                ax.text(0.5, 0.5, f'No data available\nfor {exp.display_name}\n({error_msg})', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=10)
                ax.set_title(exp.display_name, fontweight='bold')
                continue
            
            time = data['time']
            derived_energies = data['energy_components']
            
            # Plot total energies for each subsystem
            legend_added = False
            
            # Universe total energy (molecular + cavity total)
            if 'universe_total_energy' in derived_energies:
                ax.plot(time, derived_energies['universe_total_energy'], 'k-', 
                       linewidth=2, label='Universe Total', alpha=0.9)
                legend_added = True
            elif 'total_energy' in derived_energies:
                ax.plot(time, derived_energies['total_energy'], 'k-', 
                       linewidth=2, label='Universe Total', alpha=0.9)
                legend_added = True
            
            # Molecular total energy (molecular potential + kinetic)
            if 'molecular_total_energy' in derived_energies:
                ax.plot(time, derived_energies['molecular_total_energy'], 'r-', 
                       linewidth=1.5, label='Molecular Total', alpha=0.8)
                legend_added = True
            
            # Cavity total energy (cavity potential + kinetic)
            if 'cavity_total_energy' in derived_energies:
                ax.plot(time, derived_energies['cavity_total_energy'], 'b-', 
                       linewidth=1.5, label='Cavity Total', alpha=0.8)
                legend_added = True
            
            # If we don't have the derived total energies, fall back to raw components
            if not legend_added:
                if 'total_potential_energy' in derived_energies:
                    ax.plot(time, derived_energies['total_potential_energy'], 'b-', 
                           linewidth=1.5, label='Total Potential', alpha=0.8)
                    legend_added = True
                
                if 'total_kinetic_energy' in derived_energies:
                    ax.plot(time, derived_energies['total_kinetic_energy'], 'r-', 
                           linewidth=1.5, label='Total Kinetic', alpha=0.8)
                    legend_added = True
            
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('Energy (a.u.)')
            ax.set_title(exp.display_name, fontweight='bold', fontsize=12)
            ax.grid(True, alpha=0.3)
            
            if legend_added:
                ax.legend(fontsize=10)
            
            # Calculate energy conservation statistics for universe total energy
            if 'universe_total_energy' in derived_energies:
                total_energy = derived_energies['universe_total_energy']
            elif 'total_energy' in derived_energies:
                total_energy = derived_energies['total_energy']
            else:
                total_energy = None
            
            if total_energy is not None and len(total_energy) > 1:
                energy_drift = total_energy.iloc[-1] - total_energy.iloc[0]
                drift_percent = abs(energy_drift / total_energy.iloc[0]) * 100 if total_energy.iloc[0] != 0 else 0
                
                # Add energy drift information to the plot
                drift_text = f'Energy Drift:\n{energy_drift:.6f} a.u.\n({drift_percent:.3f}%)'
                ax.text(0.02, 0.02, drift_text, transform=ax.transAxes, 
                       verticalalignment='bottom', fontsize=8,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        
        plt.suptitle('Multi-Experiment Energy Comparison\n(Total Energies: Universe, Molecular, and Cavity)', 
                    fontsize=20, fontweight='bold')
        plt.tight_layout()
        output_file = self.base_dir / "master_energy_comparison.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved master energy comparison: {output_file}")
        plt.close()

    def create_master_temperature_comparison(self, all_data):
        """Create master temperature comparison plot."""
        logger.info("Generating master temperature comparison plot...")
        
        fig, axes = plt.subplots(4, 3, figsize=(18, 16))
        axes = axes.flatten()
        
        for i, exp in enumerate(self.experiments):
            ax = axes[i]
            data = all_data.get(exp.name)
            
            if not data or not data.get('valid', False):
                error_msg = data.get('error', 'Unknown error') if data else 'No data'
                ax.text(0.5, 0.5, f'No data available\nfor {exp.display_name}\n({error_msg})', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=10)
                ax.set_title(exp.display_name, fontweight='bold')
                continue
            
            time = data['time']
            temperature_data = data['temperature_data']
            
            # Calculate mean temperatures from latter half of data
            half_point = len(time) // 2
            mean_temps = {}
            
            # Plot temperatures
            if 'molecular_temperature' in temperature_data:
                mol_temp = temperature_data['molecular_temperature']
                min_len = min(len(time), len(mol_temp))
                ax.plot(time[:min_len], mol_temp[:min_len], 'r-', 
                       linewidth=1.5, label='Molecular', alpha=0.8)
                
                # Calculate mean from latter half
                if min_len > half_point:
                    mean_mol_temp = np.mean(mol_temp[half_point:min_len])
                    mean_temps['molecular'] = mean_mol_temp
            
            if 'cavity_temperature' in temperature_data:
                cav_temp = temperature_data['cavity_temperature']
                min_len = min(len(time), len(cav_temp))
                ax.plot(time[:min_len], cav_temp[:min_len], 'b-', 
                       linewidth=1.5, label='Cavity', alpha=0.8)
                
                # Calculate mean from latter half
                if min_len > half_point:
                    mean_cav_temp = np.mean(cav_temp[half_point:min_len])
                    mean_temps['cavity'] = mean_cav_temp
            
            # Add target temperature line
            ax.axhline(y=self.target_temp, color='k', linestyle='--', alpha=0.7, 
                      label=f'Target ({self.target_temp:.1f} K)')
            
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('Temperature (K)')
            ax.set_title(exp.display_name, fontweight='bold', fontsize=12)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=10)
            
            # Add mean temperature text box
            if mean_temps:
                temp_text = "Mean T (latter half):\n"
                if 'molecular' in mean_temps:
                    temp_text += f"Molecular: {mean_temps['molecular']:.1f} K\n"
                if 'cavity' in mean_temps:
                    temp_text += f"Cavity: {mean_temps['cavity']:.1f} K"
                
                ax.text(0.02, 0.98, temp_text.strip(), transform=ax.transAxes, 
                       verticalalignment='top', fontsize=9,
                       bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.suptitle('Multi-Experiment Temperature Comparison', fontsize=20, fontweight='bold')
        plt.tight_layout()
        output_file = self.base_dir / "master_temperature_comparison.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved master temperature comparison: {output_file}")
        plt.close()

    def create_master_fskt_comparison(self, all_data):
        """Create master F(k,t) comparison plot."""
        logger.info("Generating master F(k,t) comparison plot...")
        
        # Group experiments by thermostat combination
        thermostat_groups = {
            'Bussi + Langevin': [],
            'Langevin + Langevin': [],
            'None + Langevin': [],
            'Bussi + Bussi': [],
            'Langevin + Bussi': [],
            'None + Bussi': []
        }
        
        for exp in self.experiments:
            data = all_data.get(exp.name)
            if not data or not data.get('valid', False) or not data['fskt_data']:
                continue
                
            # Group by molecular + cavity thermostat combination
            group_key = f"{exp.molecular_thermostat.title()} + {exp.cavity_thermostat.title()}"
            if group_key in thermostat_groups:
                thermostat_groups[group_key].append((exp, data))
        
        # Create 6-panel plot
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        colors = {'finite_q': 'blue', 'no_finite_q': 'red'}
        
        for panel_idx, (group_name, group_experiments) in enumerate(thermostat_groups.items()):
            if panel_idx >= len(axes):
                break
                
            ax = axes[panel_idx]
            
            if not group_experiments:
                ax.text(0.5, 0.5, f'No data available\nfor {group_name}', 
                       ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{group_name}')
                continue
            
            # Plot F(k,t) for each experiment in this group
            for exp, data in group_experiments:
                fskt_data = data['fskt_data']
                if not fskt_data:
                    continue
                
                time = fskt_data['lag_time']
                fskt = fskt_data['fskt']
                
                color = colors['finite_q'] if exp.finite_q else colors['no_finite_q']
                finite_label = 'Finite-q' if exp.finite_q else 'No Finite-q'
                
                ax.plot(time, fskt, color=color, alpha=0.8, linewidth=2, label=f'{finite_label}')
            
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('F(k,t)')
            ax.set_ylim(bottom=0.0)
            ax.set_title(f'{group_name}')
            ax.grid(True, alpha=0.3)
            ax.legend()
        
        # Remove unused axes
        for panel_idx in range(len(thermostat_groups), len(axes)):
            fig.delaxes(axes[panel_idx])
        
        plt.suptitle(f'F(k,t) Comparison by Thermostat Combination (T = {self.target_temp:.1f} K)', 
                     fontsize=16, fontweight='bold')
        plt.tight_layout()
        output_file = self.base_dir / 'master_fskt_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved master F(k,t) comparison: {output_file}")
        plt.close()

    def run_analysis(self):
        """Run the complete analysis."""
        logger.info("Starting multi-experiment analysis...")
        
        # Auto-detect time ranges if not specified
        if self.max_energy_time is None or self.max_fskt_time is None:
            detected_energy_time, detected_fskt_time = self.auto_detect_time_ranges()
            if self.max_energy_time is None:
                self.max_energy_time = detected_energy_time
            if self.max_fskt_time is None:
                self.max_fskt_time = detected_fskt_time
        
        logger.info(f"Using time ranges - Energy: 0 to {self.max_energy_time:.3f} ps, F(k,t): 0 to {self.max_fskt_time:.3f} ps")
        
        # Process all experiments
        all_data = {}
        valid_count = 0
        
        # Add detailed logging for each experiment
        logger.info("Checking experiment status:")
        for exp in self.experiments:
            exp_dir = self.base_dir / exp.name
            master_file = exp_dir / "master_energy_contributions.txt"
            
            if not exp_dir.exists():
                logger.warning(f"  {exp.name}: Directory does not exist")
            elif not master_file.exists():
                logger.warning(f"  {exp.name}: Master energy file missing - {master_file}")
            else:
                logger.info(f"  {exp.name}: Master energy file found - {master_file.stat().st_size} bytes")
        
        for exp in self.experiments:
            data = self.process_single_experiment(exp)
            all_data[exp.name] = data
            if data['valid']:
                valid_count += 1
                logger.info(f"Successfully processed {exp.name}")
            else:
                logger.warning(f"Failed to process {exp.name}: {data.get('error', 'Unknown error')}")
        
        if valid_count == 0:
            logger.error("No experiment data could be processed!")
            return False
        
        logger.info(f"Successfully processed {valid_count}/{len(self.experiments)} experiments")
        
        # Create master comparison plots
        self.create_master_energy_comparison(all_data)
        self.create_master_temperature_comparison(all_data)
        self.create_master_fskt_comparison(all_data)
        
        logger.info("Multi-experiment analysis complete!")
        return True

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Process experiment folders and generate analysis plots')
    parser.add_argument('--base_dir', default='./', help='Base directory containing experiment folders')
    parser.add_argument('--max_energy_time', type=float, default=100.0, 
                       help='Maximum time for energy data (ps)')
    parser.add_argument('--max_fskt_time', type=float, default=50.0,
                       help='Maximum time for F(k,t) data (ps)')
    parser.add_argument('--energy_dt', type=float, default=0.01,
                       help='Timestep for energy data interpolation (ps)')
    parser.add_argument('--fskt_dt', type=float, default=1.0,
                       help='Timestep for F(k,t) data interpolation (ps)')
    parser.add_argument('--target_temp', type=float, default=100.0,
                       help='Target temperature for analysis (K)')
    parser.add_argument('--skip_processing', action='store_true',
                       help='Skip master file processing (use existing files)')
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    logger.info("Starting combined experiment processing and analysis...")
    logger.info(f"Base directory: {Path(args.base_dir).resolve()}")
    logger.info(f"Energy dt: {args.energy_dt} ps, F(k,t) dt: {args.fskt_dt} ps")
    logger.info(f"Max energy time: {args.max_energy_time} ps")
    logger.info(f"Max F(k,t) time: {args.max_fskt_time} ps")
    logger.info(f"Target temperature: {args.target_temp} K")
    
    success = True
    
    # Step 1: Process experiments to create master files (unless skipped)
    if not args.skip_processing:
        logger.info("\n" + "="*60)
        logger.info("STEP 1: PROCESSING EXPERIMENTS TO CREATE MASTER FILES")
        logger.info("="*60)
        
        success = process_all_experiments(
            base_dir=args.base_dir,
            max_energy_time=args.max_energy_time,
            max_fskt_time=args.max_fskt_time,
            energy_dt=args.energy_dt,
            fskt_dt=args.fskt_dt
        )
        
        if not success:
            logger.error("Failed to process experiments!")
            return 1
    else:
        logger.info("Skipping master file processing (using existing files)")
    
    # Step 2: Run analysis and generate plots
    logger.info("\n" + "="*60)
    logger.info("STEP 2: RUNNING ANALYSIS AND GENERATING PLOTS")
    logger.info("="*60)
    
    analyzer = MultiExperimentAnalyzer(
        base_dir=args.base_dir,
        target_temp=args.target_temp,
        max_energy_time=args.max_energy_time,
        max_fskt_time=args.max_fskt_time
    )
    
    analysis_success = analyzer.run_analysis()
    
    total_time = time.time() - start_time
    
    if analysis_success:
        logger.info(f"\n" + "="*60)
        logger.info("PROCESSING AND ANALYSIS COMPLETE!")
        logger.info("="*60)
        logger.info(f"Total execution time: {total_time:.2f} seconds")
        logger.info("Generated files:")
        logger.info("  - master_energy_contributions.txt (in each experiment folder)")
        logger.info("  - master_fskt_k*.txt (in each experiment folder)")
        logger.info("  - master_energy_comparison.png")
        logger.info("  - master_temperature_comparison.png")
        logger.info("  - master_fskt_comparison.png")
        return 0
    else:
        logger.error("Analysis failed!")
        return 1

if __name__ == "__main__":
    exit(main()) 
