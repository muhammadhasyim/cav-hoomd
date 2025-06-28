#!/usr/bin/env python3
"""
Energy Processing and Analysis Script for Cavity MD Simulations

This script processes energy contribution files from cavity MD simulations, creating master files
by averaging across all replicas and generating comprehensive energy comparison plots.

Key Features:
- Configurable timestep resolution for energy data
- Processes all experiment directories or focus on a single experiment
- Generates master energy files by averaging across all replicas
- Creates comprehensive energy comparison plots
- Auto-detects time ranges from data to utilize all available information
- Single experiment analysis with detailed energy component breakdown

Usage:
    python process_energy_only.py --base_dir ./
    python process_energy_only.py --energy_dt 0.01 --max_time 100.0
    python process_energy_only.py --single_experiment bussi_langevin_finiteq
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

def process_energy_files(exp_dir, job_name, max_time=None, dt=0.01):
    """
    Process energy contribution files using configurable timestep and proper cumulative averaging.
    
    Args:
        exp_dir: Path to experiment directory
        job_name: Job name ('prod' or 'finq')
        max_time: Maximum time in ps for energy data (None = auto-detect from data)
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
    
    # Auto-detect maximum time range from all files if not specified
    detected_max_time = max_time
    if max_time is None:
        all_max_times = []
        for energy_file in energy_files:
            try:
                temp_data = pd.read_csv(energy_file, sep=r'\s+', comment='#', header=None)
                if not temp_data.empty:
                    # Skip header rows if present
                    first_row = temp_data.iloc[0]
                    if any(isinstance(val, str) and not val.replace('.', '').replace('-', '').replace('e', '').replace('+', '').isdigit() for val in first_row):
                        temp_data = temp_data.iloc[1:]
                    
                    temp_data = temp_data.apply(pd.to_numeric, errors='coerce').dropna()
                    if not temp_data.empty and len(temp_data.columns) > 0:
                        times = temp_data.iloc[:, 0].values
                        valid_times = times[times >= 0]
                        if len(valid_times) > 0:
                            all_max_times.append(valid_times.max())
            except Exception as e:
                logger.warning(f"Error reading {energy_file} for time detection: {e}")
                continue
        
        if all_max_times:
            # Use the maximum time found across all replicas to include all available data
            detected_max_time = max(all_max_times)
            logger.info(f"Auto-detected max time for energy data: {detected_max_time:.3f} ps from {len(all_max_times)} files")
        else:
            detected_max_time = 100.0  # fallback
            logger.warning(f"Could not detect time range for energy data, using fallback: {detected_max_time} ps")
    
    # Set fixed linspace with configurable timestep
    num_points = int(detected_max_time / dt) + 1
    uniform_times = np.linspace(0, detected_max_time, num_points)
    logger.info(f"Using time grid: 0 to {detected_max_time:.3f} ps with {len(uniform_times)} points (dt = {dt} ps)")
    
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
            if num_cols >= 15:
                expected_columns = [
                    'time(ps)', 'timestep', 'harmonic_energy', 'lj_energy',
                    'coulomb_short_energy', 'coulomb_long_energy', 'cavity_harmonic_energy',
                    'cavity_coupling_energy', 'cavity_dipole_self_energy', 'cavity_total_potential_energy',
                    'molecular_kinetic_energy', 'cavity_mode_kinetic_energy', 'molecular_reservoir_energy',
                    'cavity_reservoir_energy'
                ]
                
                additional_columns = [
                    'molecular_mttk_thermostat_energy', 'cavity_mttk_thermostat_energy',
                    'molecular_bussi_translational_energy', 'molecular_bussi_rotational_energy',
                    'cavity_bussi_translational_energy', 'cavity_bussi_rotational_energy',
                    'total_potential_energy', 'total_kinetic_energy', 'universe_total_energy'
                ]
                
                full_expected = expected_columns + additional_columns
                if num_cols >= len(expected_columns):
                    data.columns = full_expected[:num_cols]
                else:
                    data.columns = expected_columns[:num_cols]
            else:
                # Fallback for files with fewer columns
                base_columns = ['time(ps)', 'timestep']
                if num_cols <= 2:
                    data.columns = base_columns[:num_cols]
                else:
                    data.columns = [f'col_{i}' for i in range(num_cols)]
                    if num_cols > 0:
                        data.columns = ['time(ps)'] + [f'col_{i}' for i in range(1, num_cols)]
            
            # Calculate derived energy columns if components are available
            if all(col in data.columns for col in ['harmonic_energy', 'lj_energy', 'coulomb_short_energy', 'coulomb_long_energy', 'cavity_total_potential_energy']):
                molecular_potential = data['harmonic_energy'] + data['lj_energy'] + data['coulomb_short_energy'] + data['coulomb_long_energy']
                data['total_potential_energy'] = molecular_potential + data['cavity_total_potential_energy']
            
            if all(col in data.columns for col in ['molecular_kinetic_energy', 'cavity_mode_kinetic_energy']):
                data['total_kinetic_energy'] = data['molecular_kinetic_energy'] + data['cavity_mode_kinetic_energy']
            
            if 'total_potential_energy' in data.columns and 'total_kinetic_energy' in data.columns:
                data['universe_total_energy'] = data['total_potential_energy'] + data['total_kinetic_energy']
            
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
            valid_mask = (times >= 0) & (times <= detected_max_time)
            if np.sum(valid_mask) < 3:
                logger.warning(f"Insufficient valid time points in {energy_file}: {np.sum(valid_mask)} points")
                continue
            
            times = times[valid_mask]
            data = data.iloc[valid_mask]
            
            # Log the actual time range used for this file
            file_max_time = times.max() if len(times) > 0 else 0
            logger.debug(f"File {Path(energy_file).name}: using {len(times)} points, max time: {file_max_time:.3f} ps")
            
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
                        logger.warning(f"Insufficient unique time points for column {col} in {energy_file}: {len(unique_times)} points")
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
                    
                    # Update cumulative average using proper online algorithm
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
            
            # Update sample counts only for points where we have valid data
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
        
        # Save master data with comprehensive header
        with open(master_file, 'w') as f:
            f.write(f"# Master energy contributions file for {exp_path.name}\n")
            f.write(f"# Cumulative average across {processed_count} replicas\n")
            f.write(f"# Time grid: 0 to {detected_max_time:.3f} ps with {len(uniform_times)} points (dt = {dt} ps)\n")
            f.write(f"# Valid data points: {len(clean_data)} (with at least 1 sample)\n")
            f.write(f"# Sample count range: {clean_sample_counts.min()} to {clean_sample_counts.max()}\n")
            f.write(f"# Generated by process_energy_only.py\n")
            f.write(f"# Uses standardized reservoir energy columns (molecular_reservoir_energy, cavity_reservoir_energy)\n")
            f.write(f"# that combine values from all thermostat types (Langevin + Bussi)\n")
            f.write(f"# Derived energy columns calculated: total_potential_energy, total_kinetic_energy, universe_total_energy\n")
            f.write(f"# molecular_total_energy, cavity_total_energy\n")
        
        clean_data.to_csv(master_file, sep='\t', index=False, mode='a', header=True, na_rep='nan')
        np.savetxt(sample_count_file, clean_sample_counts, fmt='%d')
        
        logger.info(f"Created master energy file: {master_file} ({processed_count} replicas)")
        logger.info(f"  Valid time points: {len(clean_data)}/{len(uniform_times)}")
        logger.info(f"  Time range: {clean_data[time_col].min():.3f} - {clean_data[time_col].max():.3f} ps")
        logger.info(f"  Sample count: min={clean_sample_counts.min()}, max={clean_sample_counts.max()}, mean={clean_sample_counts.mean():.1f}")
        return True
    
    return False 

def process_all_experiments(base_dir="./", max_time=None, dt=0.01):
    """
    Process all experiment directories to generate master energy files.
    
    Args:
        base_dir: Base directory containing experiment folders
        max_time: Maximum time for energy data (ps) - None for auto-detection
        dt: Timestep for energy data interpolation (ps)
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
    
    if max_time is None:
        logger.info("Using auto-detection of time ranges from data for each experiment")
    else:
        logger.info(f"Using fixed time range: 0 to {max_time} ps for all experiments")
    
    total_success = 0
    
    for exp_dir in experiment_dirs:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing experiment: {exp_dir.name}")
        logger.info(f"{'='*60}")
        
        # Determine job name based on experiment type
        job_name = "prod"
        
        logger.info(f"Using job name: {job_name}")
        
        # Process energy files with auto-detection or fixed time range
        logger.info("Processing energy contribution files...")
        success = process_energy_files(exp_dir, job_name, max_time, dt)
        if success:
            total_success += 1
            logger.info("✅ Energy files processed successfully")
        else:
            logger.warning("❌ Energy file processing failed")
    
    # Final summary
    logger.info(f"\n{'='*60}")
    logger.info("PROCESSING COMPLETE")
    logger.info(f"{'='*60}")
    logger.info(f"Total experiments processed: {len(experiment_dirs)}")
    logger.info(f"Energy master files created: {total_success}")
    logger.info(f"Overall success rate: {100.0 * total_success / len(experiment_dirs):.1f}%")
    
    return total_success > 0

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

class EnergyAnalyzer:
    """Analyzes energy data from all experiments and creates comparison plots."""
    
    def __init__(self, base_dir="./", target_temp=100.0, max_time=None, single_experiment=None):
        self.base_dir = Path(base_dir).resolve()
        self.target_temp = target_temp
        self.max_time = max_time  # For energy analysis, None = auto-detect
        self.single_experiment = single_experiment  # Focus on single experiment if specified
        
        # Define experiment configurations
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

        # Filter experiments if single_experiment is specified
        if self.single_experiment:
            self.experiments = [exp for exp in self.experiments if exp.name == self.single_experiment]
            if not self.experiments:
                logger.error(f"Experiment '{self.single_experiment}' not found in experiment list")
                available_experiments = [exp.name for exp in self.__class__(base_dir, target_temp, max_time).experiments]
                logger.error(f"Available experiments: {', '.join(available_experiments)}")
                raise ValueError(f"Invalid experiment name: {self.single_experiment}")

    def auto_detect_time_ranges(self):
        """Auto-detect appropriate time ranges for each experiment individually."""
        experiment_time_ranges = {}
        
        for exp in self.experiments:
            files = self.find_master_files(exp)
            if files['energy']:
                try:
                    df = read_energy_data(str(files['energy']))
                    if df is not None and not df.empty:
                        max_time = df.iloc[:, 0].max()
                        experiment_time_ranges[exp.name] = max_time
                        logger.info(f"Auto-detected time range for {exp.name}: 0 to {max_time:.3f} ps")
                except:
                    pass
        
        # Use a reasonable default for experiments without data
        default_time = 100.0
        for exp in self.experiments:
            if exp.name not in experiment_time_ranges:
                experiment_time_ranges[exp.name] = default_time
                logger.warning(f"No data found for {exp.name}, using default time range: {default_time} ps")
        
        return experiment_time_ranges

    def find_master_files(self, exp):
        """Find master data files for an experiment."""
        exp_dir = self.base_dir / exp.name
        
        files = {
            'energy': exp_dir / "master_energy_contributions.txt",
        }
        
        if not files['energy'].exists():
            files['energy'] = None
        
        return files 

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
            
            # Apply time filter if specified - use experiment-specific time range
            time = df.iloc[:, 0].values
            
            return {
                'time': time,
                'energy_components': derived_energies,
                'temperature_data': temperature_data,
                'valid': True,
                'error': None
            }
            
        except Exception as e:
            logger.error(f"Error processing {exp.name}: {e}")
            return {'valid': False, 'error': str(e)}

    def create_energy_comparison(self, all_data):
        """Create energy comparison plot."""
        logger.info("Generating energy comparison plot...")
        
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
        output_file = self.base_dir / "energy_comparison.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved energy comparison: {output_file}")
        plt.close()

    def create_temperature_comparison(self, all_data):
        """Create temperature comparison plot."""
        logger.info("Generating temperature comparison plot...")
        
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
        output_file = self.base_dir / "temperature_comparison.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved temperature comparison: {output_file}")
        plt.close()

    def create_potential_kinetic_comparison(self, all_data):
        """Create potential and kinetic energy comparison plot."""
        logger.info("Generating potential and kinetic energy comparison plot...")
        
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
            
            # Plot potential and kinetic energies
            legend_added = False
            
            # Total potential energy
            if 'total_potential_energy' in derived_energies:
                ax.plot(time, derived_energies['total_potential_energy'], 'b-', 
                       linewidth=2, label='Total Potential', alpha=0.9)
                legend_added = True
            
            # Total kinetic energy
            if 'total_kinetic_energy' in derived_energies:
                ax.plot(time, derived_energies['total_kinetic_energy'], 'r-', 
                       linewidth=2, label='Total Kinetic', alpha=0.9)
                legend_added = True
            
            # Molecular kinetic energy (subset of total kinetic)
            if 'molecular_kinetic_energy' in derived_energies:
                ax.plot(time, derived_energies['molecular_kinetic_energy'], 'r--', 
                       linewidth=1.5, label='Molecular Kinetic', alpha=0.7)
                legend_added = True
            
            # Cavity kinetic energy (subset of total kinetic)
            if 'cavity_mode_kinetic_energy' in derived_energies:
                ax.plot(time, derived_energies['cavity_mode_kinetic_energy'], 'orange', 
                       linewidth=1.5, label='Cavity Kinetic', alpha=0.7, linestyle='--')
                legend_added = True
            
            # Major potential energy components
            if 'cavity_total_potential_energy' in derived_energies:
                ax.plot(time, derived_energies['cavity_total_potential_energy'], 'cyan', 
                       linewidth=1.5, label='Cavity Potential', alpha=0.7, linestyle=':')
                legend_added = True
            
            # Calculate molecular potential if components are available
            if all(col in derived_energies for col in ['harmonic_energy', 'lj_energy', 'coulomb_short_energy', 'coulomb_long_energy']):
                molecular_potential = (derived_energies['harmonic_energy'] + 
                                     derived_energies['lj_energy'] + 
                                     derived_energies['coulomb_short_energy'] + 
                                     derived_energies['coulomb_long_energy'])
                ax.plot(time, molecular_potential, 'navy', 
                       linewidth=1.5, label='Molecular Potential', alpha=0.7, linestyle=':')
                legend_added = True
            
            ax.set_xlabel('Time (ps)')
            ax.set_ylabel('Energy (a.u.)')
            ax.set_title(exp.display_name, fontweight='bold', fontsize=12)
            ax.grid(True, alpha=0.3)
            
            if legend_added:
                ax.legend(fontsize=9, loc='best')
            
            # Calculate energy statistics
            if 'total_potential_energy' in derived_energies and 'total_kinetic_energy' in derived_energies:
                pot_energy = derived_energies['total_potential_energy']
                kin_energy = derived_energies['total_kinetic_energy']
                
                if len(pot_energy) > 1 and len(kin_energy) > 1:
                    # Calculate mean values from latter half
                    half_point = len(time) // 2
                    mean_pot = np.mean(pot_energy.iloc[half_point:]) if hasattr(pot_energy, 'iloc') else np.mean(pot_energy[half_point:])
                    mean_kin = np.mean(kin_energy.iloc[half_point:]) if hasattr(kin_energy, 'iloc') else np.mean(kin_energy[half_point:])
                    
                    # Calculate standard deviations
                    std_pot = np.std(pot_energy.iloc[half_point:]) if hasattr(pot_energy, 'iloc') else np.std(pot_energy[half_point:])
                    std_kin = np.std(kin_energy.iloc[half_point:]) if hasattr(kin_energy, 'iloc') else np.std(kin_energy[half_point:])
                    
                    # Add statistics text box
                    stats_text = f"Mean ± Std (latter half):\nPotential: {mean_pot:.4f} ± {std_pot:.4f}\nKinetic: {mean_kin:.4f} ± {std_kin:.4f}"
                    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                           verticalalignment='top', fontsize=8,
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.suptitle('Multi-Experiment Potential and Kinetic Energy Comparison', 
                    fontsize=20, fontweight='bold')
        plt.tight_layout()
        output_file = self.base_dir / "potential_kinetic_comparison.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved potential and kinetic energy comparison: {output_file}")
        plt.close()

    def create_single_experiment_energy_analysis(self, all_data):
        """Create detailed energy analysis for a single experiment."""
        if len(self.experiments) != 1:
            logger.error("Single experiment analysis requires exactly one experiment")
            return
        
        exp = self.experiments[0]
        data = all_data.get(exp.name)
        
        if not data or not data.get('valid', False):
            logger.error(f"No valid data for experiment {exp.name}")
            return
        
        logger.info(f"Generating detailed energy analysis for {exp.display_name}...")
        
        time = data['time']
        derived_energies = data['energy_components']
        temperature_data = data['temperature_data']
        
        # Create a comprehensive 3-panel plot
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 16))
        
        # Panel 1: Total Energy Components
        ax1.set_title(f'Energy Components - {exp.display_name}', fontweight='bold', fontsize=14)
        
        # Universe total energy
        if 'universe_total_energy' in derived_energies:
            ax1.plot(time, derived_energies['universe_total_energy'], 'k-', 
                    linewidth=2.5, label='Universe Total', alpha=0.9)
        
        # Molecular and cavity total energies
        if 'molecular_total_energy' in derived_energies:
            ax1.plot(time, derived_energies['molecular_total_energy'], 'r-', 
                    linewidth=2, label='Molecular Total', alpha=0.8)
        
        if 'cavity_total_energy' in derived_energies:
            ax1.plot(time, derived_energies['cavity_total_energy'], 'b-', 
                    linewidth=2, label='Cavity Total', alpha=0.8)
        
        ax1.set_xlabel('Time (ps)')
        ax1.set_ylabel('Energy (a.u.)')
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=11)
        
        # Calculate energy drift statistics
        if 'universe_total_energy' in derived_energies:
            total_energy = derived_energies['universe_total_energy']
            if len(total_energy) > 1:
                energy_drift = total_energy.iloc[-1] - total_energy.iloc[0]
                drift_percent = abs(energy_drift / total_energy.iloc[0]) * 100 if total_energy.iloc[0] != 0 else 0
                
                drift_text = f'Energy Conservation:\nDrift: {energy_drift:.6f} a.u.\nRelative: {drift_percent:.4f}%'
                ax1.text(0.02, 0.98, drift_text, transform=ax1.transAxes, 
                        verticalalignment='top', fontsize=10,
                        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
        
        # Panel 2: Potential vs Kinetic Energy Breakdown
        ax2.set_title('Potential vs Kinetic Energy Breakdown', fontweight='bold', fontsize=14)
        
        # Total potential and kinetic
        if 'total_potential_energy' in derived_energies:
            ax2.plot(time, derived_energies['total_potential_energy'], 'b-', 
                    linewidth=2.5, label='Total Potential', alpha=0.9)
        
        if 'total_kinetic_energy' in derived_energies:
            ax2.plot(time, derived_energies['total_kinetic_energy'], 'r-', 
                    linewidth=2.5, label='Total Kinetic', alpha=0.9)
        
        # Kinetic energy components
        if 'molecular_kinetic_energy' in derived_energies:
            ax2.plot(time, derived_energies['molecular_kinetic_energy'], 'r--', 
                    linewidth=1.5, label='Molecular Kinetic', alpha=0.7)
        
        if 'cavity_mode_kinetic_energy' in derived_energies:
            ax2.plot(time, derived_energies['cavity_mode_kinetic_energy'], 'orange', 
                    linewidth=1.5, label='Cavity Kinetic', alpha=0.7, linestyle='--')
        
        # Potential energy components
        if 'cavity_total_potential_energy' in derived_energies:
            ax2.plot(time, derived_energies['cavity_total_potential_energy'], 'cyan', 
                    linewidth=1.5, label='Cavity Potential', alpha=0.7, linestyle=':')
        
        # Calculate molecular potential if components are available
        if all(col in derived_energies for col in ['harmonic_energy', 'lj_energy', 'coulomb_short_energy', 'coulomb_long_energy']):
            molecular_potential = (derived_energies['harmonic_energy'] + 
                                 derived_energies['lj_energy'] + 
                                 derived_energies['coulomb_short_energy'] + 
                                 derived_energies['coulomb_long_energy'])
            ax2.plot(time, molecular_potential, 'navy', 
                    linewidth=1.5, label='Molecular Potential', alpha=0.7, linestyle=':')
        
        ax2.set_xlabel('Time (ps)')
        ax2.set_ylabel('Energy (a.u.)')
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=10, loc='best')
        
        # Calculate and display energy statistics
        if 'total_potential_energy' in derived_energies and 'total_kinetic_energy' in derived_energies:
            pot_energy = derived_energies['total_potential_energy']
            kin_energy = derived_energies['total_kinetic_energy']
            
            if len(pot_energy) > 1 and len(kin_energy) > 1:
                # Calculate mean values from latter half
                half_point = len(time) // 2
                mean_pot = np.mean(pot_energy.iloc[half_point:]) if hasattr(pot_energy, 'iloc') else np.mean(pot_energy[half_point:])
                mean_kin = np.mean(kin_energy.iloc[half_point:]) if hasattr(kin_energy, 'iloc') else np.mean(kin_energy[half_point:])
                
                # Calculate standard deviations
                std_pot = np.std(pot_energy.iloc[half_point:]) if hasattr(pot_energy, 'iloc') else np.std(pot_energy[half_point:])
                std_kin = np.std(kin_energy.iloc[half_point:]) if hasattr(kin_energy, 'iloc') else np.std(kin_energy[half_point:])
                
                # Add statistics text box
                stats_text = f"Mean ± Std (latter half):\nPotential: {mean_pot:.4f} ± {std_pot:.4f}\nKinetic: {mean_kin:.4f} ± {std_kin:.4f}"
                ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes, 
                        verticalalignment='top', fontsize=10,
                        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        # Panel 3: Temperature Analysis
        ax3.set_title('Temperature Analysis', fontweight='bold', fontsize=14)
        
        # Calculate mean temperatures from latter half of data
        half_point = len(time) // 2
        mean_temps = {}
        
        # Plot temperatures
        if 'molecular_temperature' in temperature_data:
            mol_temp = temperature_data['molecular_temperature']
            min_len = min(len(time), len(mol_temp))
            ax3.plot(time[:min_len], mol_temp[:min_len], 'r-', 
                    linewidth=2, label='Molecular Temperature', alpha=0.8)
            
            # Calculate mean from latter half
            if min_len > half_point:
                mean_mol_temp = np.mean(mol_temp[half_point:min_len])
                std_mol_temp = np.std(mol_temp[half_point:min_len])
                mean_temps['molecular'] = (mean_mol_temp, std_mol_temp)
        
        if 'cavity_temperature' in temperature_data:
            cav_temp = temperature_data['cavity_temperature']
            min_len = min(len(time), len(cav_temp))
            ax3.plot(time[:min_len], cav_temp[:min_len], 'b-', 
                    linewidth=2, label='Cavity Temperature', alpha=0.8)
            
            # Calculate mean from latter half
            if min_len > half_point:
                mean_cav_temp = np.mean(cav_temp[half_point:min_len])
                std_cav_temp = np.std(cav_temp[half_point:min_len])
                mean_temps['cavity'] = (mean_cav_temp, std_cav_temp)
        
        # Add target temperature line
        ax3.axhline(y=self.target_temp, color='k', linestyle='--', alpha=0.7, linewidth=2,
                   label=f'Target Temperature ({self.target_temp:.1f} K)')
        
        ax3.set_xlabel('Time (ps)')
        ax3.set_ylabel('Temperature (K)')
        ax3.grid(True, alpha=0.3)
        ax3.legend(fontsize=11)
        
        # Add detailed temperature statistics
        if mean_temps:
            temp_text = "Temperature Statistics (latter half):\n"
            if 'molecular' in mean_temps:
                mean_val, std_val = mean_temps['molecular']
                temp_text += f"Molecular: {mean_val:.1f} ± {std_val:.1f} K\n"
            if 'cavity' in mean_temps:
                mean_val, std_val = mean_temps['cavity']
                temp_text += f"Cavity: {mean_val:.1f} ± {std_val:.1f} K"
            
            ax3.text(0.02, 0.98, temp_text.strip(), transform=ax3.transAxes, 
                    verticalalignment='top', fontsize=10,
                    bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
        
        plt.suptitle(f'Detailed Energy Analysis: {exp.display_name}', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        # Save the plot
        output_file = self.base_dir / f"single_energy_analysis_{exp.name}.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved single experiment energy analysis: {output_file}")
        plt.close()

    def run_analysis(self):
        """Run the complete energy analysis."""
        logger.info("Starting energy analysis...")
        
        # Auto-detect time ranges for each experiment if not specified
        if self.max_time is None:
            experiment_time_ranges = self.auto_detect_time_ranges()
            if self.single_experiment:
                logger.info(f"Using individual time range for experiment {self.single_experiment}")
            else:
                logger.info("Using individual time ranges for each experiment")
        else:
            # Use the same time range for all experiments
            experiment_time_ranges = {exp.name: self.max_time for exp in self.experiments}
            if self.single_experiment:
                logger.info(f"Using fixed time range for experiment {self.single_experiment}: 0 to {self.max_time:.3f} ps")
            else:
                logger.info(f"Using fixed time range for all experiments: 0 to {self.max_time:.3f} ps")
        
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
                time_range = experiment_time_ranges.get(exp.name, 100.0)
                logger.info(f"  {exp.name}: Master energy file found - {master_file.stat().st_size} bytes, time range: 0 to {time_range:.3f} ps")
        
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
        
        # Create plots based on analysis mode
        if self.single_experiment:
            # Single experiment mode - create detailed analysis
            self.create_single_experiment_energy_analysis(all_data)
        else:
            # Multi-experiment mode - create comparison plots
            self.create_energy_comparison(all_data)
            self.create_temperature_comparison(all_data)
            self.create_potential_kinetic_comparison(all_data)
        
        logger.info("Energy analysis complete!")
        return True

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Process energy files and generate analysis plots')
    parser.add_argument('--base_dir', default='./', help='Base directory containing experiment folders')
    parser.add_argument('--max_time', type=float, default=None,
                       help='Maximum time for energy data (ps) - if not specified, auto-detect from data')
    parser.add_argument('--energy_dt', type=float, default=0.01,
                       help='Timestep for energy data interpolation (ps)')
    parser.add_argument('--target_temp', type=float, default=100.0,
                       help='Target temperature for analysis (K)')
    parser.add_argument('--single_experiment', type=str, default=None,
                       help='Process only a single experiment (e.g., bussi_langevin_finiteq)')
    parser.add_argument('--skip_processing', action='store_true',
                       help='Skip master file processing (use existing files)')
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    logger.info("Starting energy processing and analysis...")
    logger.info(f"Base directory: {Path(args.base_dir).resolve()}")
    logger.info(f"Energy dt: {args.energy_dt} ps")
    logger.info(f"Target temperature: {args.target_temp} K")
    
    if args.single_experiment:
        logger.info(f"Single experiment mode: {args.single_experiment}")
    
    if args.max_time is not None:
        logger.info(f"Max energy time: {args.max_time} ps (fixed)")
    else:
        logger.info("Max energy time: auto-detect from data")
    
    success = True
    
    # Step 1: Process experiments to create master energy files (unless skipped)
    if not args.skip_processing:
        logger.info("\n" + "="*60)
        logger.info("STEP 1: PROCESSING EXPERIMENTS TO CREATE MASTER ENERGY FILES")
        logger.info("="*60)
        
        success = process_all_experiments(
            base_dir=args.base_dir,
            max_time=args.max_time,
            dt=args.energy_dt
        )
        
        if not success:
            logger.error("Failed to process experiments!")
            return 1
    else:
        logger.info("Skipping master file processing (using existing files)")
    
    # Step 2: Run energy analysis and generate plots
    logger.info("\n" + "="*60)
    logger.info("STEP 2: RUNNING ENERGY ANALYSIS AND GENERATING PLOTS")
    logger.info("="*60)
    
    analyzer = EnergyAnalyzer(
        base_dir=args.base_dir,
        target_temp=args.target_temp,
        max_time=args.max_time,
        single_experiment=args.single_experiment
    )
    
    analysis_success = analyzer.run_analysis()
    
    total_time = time.time() - start_time
    
    if analysis_success:
        logger.info(f"\n" + "="*60)
        logger.info("ENERGY PROCESSING AND ANALYSIS COMPLETE!")
        logger.info("="*60)
        logger.info(f"Total execution time: {total_time:.2f} seconds")
        logger.info("Generated files:")
        if args.single_experiment:
            logger.info(f"  - single_energy_analysis_{args.single_experiment}.png (detailed energy analysis)")
        else:
            logger.info("  - master_energy_contributions.txt (in each experiment folder)")
            logger.info("  - energy_comparison.png (total energy comparison)")
            logger.info("  - temperature_comparison.png (temperature analysis)")
            logger.info("  - potential_kinetic_comparison.png (potential vs kinetic energies)")
        return 0
    else:
        logger.error("Energy analysis failed!")
        return 1

if __name__ == "__main__":
    exit(main()) 