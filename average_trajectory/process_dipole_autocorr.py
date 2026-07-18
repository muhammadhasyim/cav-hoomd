#!/usr/bin/env python3
"""
Dipole Autocorrelation Processing and IR Spectrum Analysis Script for Cavity MD Simulations

This script processes dipole autocorrelation files from a single cavity MD simulation directory,
creating master files by averaging across all replicas and generating IR spectra using DCT.

Key Features:
- Configurable timestep resolution for autocorrelation data
- Processes a single specified experiment directory
- Generates master autocorrelation files by averaging across all replicas
- Creates IR spectra using Discrete Cosine Transform (DCT)
- Comprehensive analysis plots for autocorrelations and spectra

Usage:
    python process_dipole_autocorr.py --exp_dir ./cavity_coupling_1eneg03_switch_1.0ps
    python process_dipole_autocorr.py --exp_dir ./my_experiment --autocorr_dt 0.5 --max_time 100.0
    python process_dipole_autocorr.py --exp_dir ./test_run --skip_processing  # Only create plots
"""
 
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import glob
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.interpolate import interp1d
from scipy import fftpack
import time
import natsort
from pathlib import Path
import argparse
import logging
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from tqdm import tqdm
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants for IR spectrum calculation
BOLTZ = 1.38064852E-23  # m^2 kg s^-2 K^-1
LIGHTSPEED = 299792458.  # m s^-1
REDUCED_PLANCK = 1.05457180013E-34  # kg m^2 s^-1

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

def process_dipole_autocorr_files(exp_dir, job_name="prod", max_time=None, dt=1.0):
    """
    Process dipole autocorrelation files using configurable timestep and proper cumulative averaging.
    
    Args:
        exp_dir: Path to experiment directory
        job_name: Job name ('prod' or 'finq')
        max_time: Maximum time in fs for autocorrelation data (None = auto-detect from data)
        dt: Timestep in fs for interpolation
    
    Returns:
        Dictionary with processing results
    """
    exp_path = Path(exp_dir)
    
    if not exp_path.exists():
        logger.error(f"Experiment directory does not exist: {exp_dir}")
        return {'success': False, 'error': 'Directory does not exist'}
    
    # Find all dipole autocorrelation files - pattern: prod-i_dipole_autocorr_j.txt
    autocorr_pattern = f"{job_name}-*_dipole_autocorr_*.txt"
    autocorr_files = list(exp_path.glob(autocorr_pattern))
    
    if not autocorr_files:
        logger.warning(f"No dipole autocorrelation files found in {exp_dir} with pattern {autocorr_pattern}")
        return {'success': False, 'error': f'No files found with pattern {autocorr_pattern}'}
    
    # Group files by reference number (extracting j from prod-i_dipole_autocorr_j.txt)
    autocorr_groups = {}
    for autocorr_file in autocorr_files:
        filename = autocorr_file.name
        # Extract reference number from filename like prod-0_dipole_autocorr_1.txt -> ref1
        if '_dipole_autocorr_' in filename:
            ref_part = filename.split('_dipole_autocorr_')[1].replace('.txt', '')
            group_key = f"ref{ref_part}"
        else:
            # Fallback: use the whole filename without extension as group key
            group_key = filename.replace('.txt', '')
        
        if group_key not in autocorr_groups:
            autocorr_groups[group_key] = []
        autocorr_groups[group_key].append(str(autocorr_file))
    
    logger.info(f"Processing {len(autocorr_groups)} dipole autocorrelation reference groups in {exp_dir}")
    logger.info(f"Found reference groups: {list(autocorr_groups.keys())}")
    
    results = {'success': True, 'processed_groups': {}, 'total_groups': len(autocorr_groups)}
    success_count = 0
    
    # Process autocorr groups
    for group_key, file_list in tqdm(autocorr_groups.items(), desc="Processing autocorrelation groups"):
        # Sort files naturally
        file_list = natsort.natsorted(file_list)
        
        # Auto-detect maximum time range from all files if not specified
        detected_max_time = max_time
        if max_time is None:
            all_max_times = []
            # Process time detection with improved error handling
            for autocorr_file in file_list:
                try:
                    temp_data = pd.read_csv(autocorr_file, sep=r'\s+', comment='#', header=None)
                    if not temp_data.empty and len(temp_data.columns) >= 2:
                        # Skip header rows if present
                        first_row = temp_data.iloc[0]
                        if any(isinstance(val, str) and not val.replace('.', '').replace('-', '').replace('e', '').replace('+', '').isdigit() for val in first_row):
                            temp_data = temp_data.iloc[1:]
                        
                        temp_data = temp_data.apply(pd.to_numeric, errors='coerce').dropna()
                        if not temp_data.empty and len(temp_data) >= 10:  # Require at least 10 data points
                            times = temp_data.iloc[:, 0].values
                            valid_times = times[times >= 0]
                            if len(valid_times) >= 10:  # Require at least 10 valid time points
                                all_max_times.append(valid_times.max())
                                logger.debug(f"File {Path(autocorr_file).name}: {len(valid_times)} time points, max time: {valid_times.max():.3f} fs")
                except Exception as e:
                    logger.debug(f"Could not read {autocorr_file} for time detection: {e}")
                    continue
            
            if all_max_times:
                # Use a conservative approach: use the minimum time found to ensure all files have data
                detected_max_time = min(all_max_times)
                logger.info(f"Auto-detected max time for {group_key}: {detected_max_time:.3f} fs from {len(all_max_times)} valid files")
                logger.info(f"  Using conservative estimate (minimum across files) to handle incomplete data")
                logger.info(f"  Time range across files: {min(all_max_times):.3f} to {max(all_max_times):.3f} fs")
            else:
                detected_max_time = 1000.0  # Conservative fallback for running simulations (1 ps)
                logger.warning(f"Could not detect time range for {group_key} (files may be empty/incomplete), using conservative fallback: {detected_max_time} fs")
        
        # Set fixed linspace with configurable timestep
        if detected_max_time is not None:
            num_points = int(detected_max_time / dt) + 1
            uniform_times = np.linspace(0, detected_max_time, num_points)
            logger.info(f"Autocorrelation group {group_key}: Using time grid: 0 to {detected_max_time:.3f} fs with {len(uniform_times)} points (dt = {dt} fs)")
        else:
            logger.error(f"Could not determine time range for {group_key}")
            continue
        
        # Initialize storage for cumulative data
        cumulative_autocorr = np.zeros(len(uniform_times))
        sample_counts = np.zeros(len(uniform_times), dtype=int)
        processed_count = 0
        
        # Process each autocorrelation file in this group
        logger.info(f"Processing {group_key}: {len(file_list)} files...")
        for autocorr_file in file_list:
            try:
                # Read the autocorrelation data
                data = pd.read_csv(autocorr_file, sep=r'\s+', comment='#', header=None)
                
                if data.empty:
                    continue
                
                # Convert all data to numeric, coercing errors to NaN
                data = data.apply(pd.to_numeric, errors='coerce')
                
                # Drop any rows with NaN values
                data = data.dropna()
                
                if data.empty:
                    continue
                
                # Assume 2-column format: time(fs), autocorrelation
                if data.shape[1] >= 2:
                    times = data.iloc[:, 0].values  # First column is time (fs)
                    autocorr_vals = data.iloc[:, 1].values  # Second column is autocorrelation
                else:
                    continue
                
                # Filter to valid time range - be more lenient for incomplete files
                valid_mask = (times >= 0)
                if np.sum(valid_mask) < 3:
                    continue
                
                times = times[valid_mask]
                autocorr_vals = autocorr_vals[valid_mask]
                
                # Remove duplicate time values before interpolation
                unique_times, unique_autocorr = remove_duplicate_times(times, autocorr_vals)
                
                if len(unique_times) < 2:
                    continue
                
                # Check how much of our target grid this file can cover
                grid_coverage = np.sum((uniform_times >= unique_times.min()) & 
                                     (uniform_times <= unique_times.max()))
                
                if grid_coverage == 0:
                    continue
                
                # Interpolate to uniform grid within the file's time range
                valid_grid_mask = (uniform_times >= unique_times.min()) & (uniform_times <= unique_times.max())
                valid_grid_times = uniform_times[valid_grid_mask]
                
                if len(valid_grid_times) == 0:
                    continue
                
                # Perform interpolation
                interpolated_autocorr = np.interp(valid_grid_times, unique_times, unique_autocorr)
                
                # Update cumulative average using proper online algorithm
                for i, grid_idx in enumerate(np.where(valid_grid_mask)[0]):
                    n = sample_counts[grid_idx] + 1
                    if sample_counts[grid_idx] == 0:
                        cumulative_autocorr[grid_idx] = interpolated_autocorr[i]
                    else:
                        # Online average update: new_avg = old_avg + (new_val - old_avg) / n
                        old_avg = cumulative_autocorr[grid_idx]
                        new_val = interpolated_autocorr[i]
                        cumulative_autocorr[grid_idx] = old_avg + (new_val - old_avg) / n
                    sample_counts[grid_idx] = n
                
                processed_count += 1
                
            except Exception as e:
                logger.warning(f"Error processing {autocorr_file}: {e}")
                continue
        
        # Save results if we have data
        if processed_count > 0:
            master_file = exp_path / f"master_dipole_autocorr_{group_key}.txt"
            sample_count_file = exp_path / f"master_dipole_autocorr_{group_key}_sample_counts.txt"
            
            # Only keep time points where we have at least one sample
            valid_points = sample_counts > 0
            clean_times = uniform_times[valid_points]
            clean_autocorr = cumulative_autocorr[valid_points]
            clean_sample_counts = sample_counts[valid_points]
            
            if len(clean_times) == 0:
                logger.warning(f"No valid data points found for {group_key} after processing {processed_count} files")
                continue
            
            # Create output data
            output_data = pd.DataFrame({
                'time_fs': clean_times,
                'autocorr': clean_autocorr
            })
            
            # Calculate data coverage statistics
            total_possible_points = len(uniform_times)
            actual_data_points = len(clean_times)
            coverage_percent = 100.0 * actual_data_points / total_possible_points
            
            # Save master data with comprehensive header
            with open(master_file, 'w') as f:
                f.write(f"# Master dipole autocorrelation file for {exp_path.name} - {group_key}\n")
                f.write(f"# Cumulative average across {processed_count} replicas\n")
                f.write(f"# Time grid: 0 to {detected_max_time:.3f} fs with {len(uniform_times)} points (dt = {dt} fs)\n")
                f.write(f"# Valid data points: {len(clean_times)} ({coverage_percent:.1f}% coverage)\n")
                f.write(f"# Sample count range: {clean_sample_counts.min()} to {clean_sample_counts.max()}\n")
                f.write(f"# Note: Simulation may be incomplete - sample counts indicate data availability\n")
                f.write(f"# Generated by process_dipole_autocorr.py\n")
                f.write(f"# Columns: time_fs autocorr\n")
            
            output_data.to_csv(master_file, sep='\t', index=False, mode='a', header=True, na_rep='nan')
            np.savetxt(sample_count_file, clean_sample_counts, fmt='%d')
            
            logger.info(f"✓ {group_key}: {processed_count} files → {len(clean_times)} points ({coverage_percent:.1f}% coverage)")
            
            if coverage_percent < 50:
                logger.warning(f"  Low coverage for {group_key} ({coverage_percent:.1f}%) - simulation may still be running")
            
            results['processed_groups'][group_key] = {
                'master_file': str(master_file),
                'sample_count_file': str(sample_count_file),
                'processed_count': processed_count,
                'time_range': (clean_times.min(), clean_times.max()),
                'valid_points': len(clean_times),
                'coverage_percent': coverage_percent
            }
            success_count += 1
        else:
            logger.warning(f"✗ {group_key}: No files could be processed")
    
    results['success'] = success_count > 0
    results['success_count'] = success_count
    
    return results

def read_autocorr_data(file_path):
    """Read dipole autocorrelation data from file."""
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
        
        if 'time_fs' in data.columns and 'autocorr' in data.columns:
            return {'time_fs': data['time_fs'].values, 'autocorr': data['autocorr'].values}
        elif data.shape[1] >= 2:
            return {'time_fs': data.iloc[:, 0].values, 'autocorr': data.iloc[:, 1].values}
        else:
            return None
    except Exception as e:
        logger.error(f"Error reading autocorrelation data from {file_path}: {e}")
        return None

def calculate_ir_spectrum(time_fs, autocorr, temperature=300.0, max_freq_cm=4000.0):
    """
    Calculate IR spectrum from dipole autocorrelation function using DCT.
    
    Args:
        time_fs: Time array in femtoseconds
        autocorr: Autocorrelation function values
        temperature: Temperature in Kelvin
        max_freq_cm: Maximum frequency for plotting in cm^-1
    
    Returns:
        Dictionary with frequency and spectrum data
    """
    # Convert time from fs to seconds
    timestep = (time_fs[1] - time_fs[0]) * 1.0E-15
    
    # Calculate DCT (lineshape)
    lineshape = fftpack.dct(autocorr, type=1)[1:]
    lineshape_frequencies = np.linspace(0, 0.5/timestep, len(autocorr))[1:]
    lineshape_frequencies_wn = lineshape_frequencies / (100.0 * LIGHTSPEED)  # Convert to wavenumbers (cm^-1)
    
    # Calculate quantum correction and field description
    field_description = lineshape_frequencies * (1.0 - np.exp(-REDUCED_PLANCK * lineshape_frequencies / (BOLTZ * temperature)))
    quantum_correction = lineshape_frequencies / (1.0 - np.exp(-REDUCED_PLANCK * lineshape_frequencies / (BOLTZ * temperature)))
    
    # Calculate final spectrum with quantum correction
    spectra = lineshape * field_description
    spectra_qm = spectra * quantum_correction
    
    # Create mask for plotting range
    mask = (lineshape_frequencies_wn >= 0) & (lineshape_frequencies_wn <= max_freq_cm)
    
    return {
        'frequencies_cm': lineshape_frequencies_wn,
        'lineshape': lineshape,
        'field_description': field_description,
        'quantum_correction': quantum_correction,
        'spectra': spectra,
        'spectra_qm': spectra_qm,
        'plot_mask': mask
    }

class SingleDirectoryDipoleAutocorrAnalyzer:
    """Analyzes dipole autocorrelation data from a single experiment directory and creates IR spectra."""
    
    def __init__(self, exp_dir, max_time=None, temperature=300.0):
        self.exp_dir = Path(exp_dir).resolve()
        self.max_time = max_time
        self.temperature = temperature
        self.exp_name = self.exp_dir.name
        
        logger.info(f"Initialized analyzer for directory: {self.exp_dir}")

    def find_autocorr_files(self):
        """Find master dipole autocorrelation data files in the experiment directory."""
        if not self.exp_dir.exists():
            logger.error(f"Experiment directory does not exist: {self.exp_dir}")
            return []
        
        autocorr_files = list(self.exp_dir.glob("master_dipole_autocorr_ref*.txt"))
        valid_files = [f for f in autocorr_files if f.stat().st_size > 0]
        
        if valid_files:
            logger.info(f"Found {len(valid_files)} master dipole autocorrelation files")
            for f in valid_files:
                logger.info(f"  {f.name}")
        else:
            logger.warning("No master dipole autocorrelation files found")
        
        return sorted(valid_files)

    def create_autocorr_overview_plot(self, autocorr_files):
        """Create overview plot showing all autocorrelation functions from the directory."""
        logger.info("Generating dipole autocorrelation overview plot...")
        
        if not autocorr_files:
            logger.warning("No autocorrelation files to plot")
            return False
        
        # Create figure with appropriate size
        n_files = len(autocorr_files)
        if n_files <= 4:
            rows, cols = 2, 2
        elif n_files <= 6:
            rows, cols = 2, 3
        elif n_files <= 9:
            rows, cols = 3, 3
        else:
            rows = int(np.ceil(n_files / 4))
            cols = 4
        
        fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 3*rows))
        
        # Handle single subplot case
        if n_files == 1:
            axes = [axes]
        elif rows == 1 or cols == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        colors = plt.cm.tab10(np.linspace(0, 1, n_files))
        
        for i, autocorr_file in enumerate(autocorr_files):
            if i >= len(axes):
                break
                
            ax = axes[i]
            
            # Extract reference group from filename
            ref_group = autocorr_file.name.replace('master_dipole_autocorr_', '').replace('.txt', '')
            
            # Read and plot data
            autocorr_data = read_autocorr_data(str(autocorr_file))
            if autocorr_data:
                time_fs = autocorr_data['time_fs']
                autocorr = autocorr_data['autocorr']
                
                ax.plot(time_fs, autocorr, color=colors[i], linewidth=2, alpha=0.8)
                
                # Add statistics
                if len(autocorr) > 10:
                    initial_value = autocorr[0]
                    final_value = autocorr[-1]
                    min_value = np.min(autocorr)
                    max_value = np.max(autocorr)
                    
                    stats_text = f"Initial: {initial_value:.4f}\nFinal: {final_value:.4f}\nMin: {min_value:.4f}\nMax: {max_value:.4f}"
                    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                           verticalalignment='top', fontsize=8,
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            else:
                ax.text(0.5, 0.5, f'No valid data\nfor {ref_group}', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=10)
            
            ax.set_xlabel('Time (fs)', fontsize=10)
            ax.set_ylabel('Autocorrelation', fontsize=10)
            ax.set_title(f'{ref_group}', fontweight='bold', fontsize=11)
            ax.grid(True, alpha=0.3)
        
        # Remove unused axes
        for i in range(len(autocorr_files), len(axes)):
            fig.delaxes(axes[i])
        
        plt.suptitle(f'Dipole Autocorrelation Analysis - {self.exp_name}\nAll Reference Groups', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.exp_dir / "dipole_autocorr_overview.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved autocorrelation overview plot: {output_file}")
        plt.close()
        
        return True

    def create_ir_spectra_plot(self, autocorr_files):
        """Create IR spectra plot from autocorrelation functions."""
        logger.info("Generating IR spectra plot...")
        
        if not autocorr_files:
            logger.warning("No autocorrelation files to process for IR spectra")
            return False
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(autocorr_files)))
        has_data = False
        
        for i, autocorr_file in enumerate(autocorr_files):
            # Extract reference group from filename
            ref_group = autocorr_file.name.replace('master_dipole_autocorr_', '').replace('.txt', '')
            
            # Read autocorrelation data
            autocorr_data = read_autocorr_data(str(autocorr_file))
            if autocorr_data:
                time_fs = autocorr_data['time_fs']
                autocorr = autocorr_data['autocorr']
                
                # Plot autocorrelation function
                ax1.plot(time_fs, autocorr, color=colors[i], linewidth=2.5, 
                        alpha=0.8, label=ref_group)
                
                # Calculate and plot IR spectrum
                try:
                    ir_data = calculate_ir_spectrum(time_fs, autocorr, self.temperature)
                    frequencies_cm = ir_data['frequencies_cm']
                    spectra_qm = ir_data['spectra_qm']
                    mask = ir_data['plot_mask']
                    
                    ax2.plot(frequencies_cm[mask], spectra_qm[mask], color=colors[i], 
                            linewidth=2.5, alpha=0.8, label=ref_group)
                    
                    # Save IR spectrum data
                    spectrum_file = self.exp_dir / f"ir_spectrum_{ref_group}.txt"
                    spectrum_data = pd.DataFrame({
                        'frequency_cm': frequencies_cm,
                        'spectra_qm': spectra_qm
                    })
                    
                    with open(spectrum_file, 'w') as f:
                        f.write(f"# IR spectrum for {self.exp_name} - {ref_group}\n")
                        f.write(f"# Temperature: {self.temperature} K\n")
                        f.write(f"# Calculated using DCT from dipole autocorrelation\n")
                        f.write(f"# Columns: frequency_cm spectra_qm\n")
                    
                    spectrum_data.to_csv(spectrum_file, sep='\t', index=False, mode='a', header=True)
                    
                    has_data = True
                    
                except Exception as e:
                    logger.warning(f"Could not calculate IR spectrum for {ref_group}: {e}")
                    continue
        
        if not has_data:
            ax1.text(0.5, 0.5, f'No valid autocorrelation data found\nin {self.exp_name}', 
                    ha='center', va='center', transform=ax1.transAxes, fontsize=14)
            ax2.text(0.5, 0.5, f'No valid IR spectra calculated\nfor {self.exp_name}', 
                    ha='center', va='center', transform=ax2.transAxes, fontsize=14)
            logger.warning(f"No valid data found for {self.exp_name}")
            plt.close()
            return False
        
        # Customize autocorrelation plot
        ax1.set_xlabel('Time (fs)', fontsize=14)
        ax1.set_ylabel('Dipole Autocorrelation', fontsize=14)
        ax1.set_title(f'Master Dipole Autocorrelation Functions - {self.exp_name}', 
                     fontweight='bold', fontsize=16)
        ax1.grid(True, alpha=0.3)
        ax1.legend(fontsize=10, loc='best')
        
        # Customize IR spectrum plot
        ax2.set_xlabel(r'Frequency (cm$^{-1}$)', fontsize=14)
        ax2.set_ylabel('IR Intensity (arb. units)', fontsize=14)
        ax2.set_title(f'IR Spectra - {self.exp_name} (T = {self.temperature} K)', 
                     fontweight='bold', fontsize=16)
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=10, loc='best')
        ax2.invert_xaxis()  # Conventional IR plot orientation
        
        plt.tight_layout()
        
        # Save the plot
        output_file = self.exp_dir / 'dipole_autocorr_and_ir_spectra.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved combined autocorrelation and IR spectra plot: {output_file}")
        plt.close()
        
        return True

    def create_summary_report(self, autocorr_files, processing_results=None):
        """Create a summary report of the dipole autocorrelation and IR analysis."""
        logger.info("Generating dipole autocorrelation summary report...")
        
        report_file = self.exp_dir / "dipole_autocorr_analysis_summary.txt"
        
        with open(report_file, 'w') as f:
            f.write(f"Dipole Autocorrelation and IR Spectrum Analysis Summary Report\n")
            f.write(f"{'='*70}\n\n")
            f.write(f"Experiment Directory: {self.exp_dir}\n")
            f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Temperature: {self.temperature} K\n\n")
            
            if processing_results:
                f.write(f"Processing Results:\n")
                f.write(f"  Total groups processed: {processing_results.get('success_count', 0)}/{processing_results.get('total_groups', 0)}\n")
                f.write(f"  Processing success: {processing_results.get('success', False)}\n\n")
            
            f.write(f"Master Autocorrelation Files Found: {len(autocorr_files)}\n\n")
            
            if autocorr_files:
                f.write("File Details:\n")
                f.write("-" * 40 + "\n")
                
                for autocorr_file in autocorr_files:
                    ref_group = autocorr_file.name.replace('master_dipole_autocorr_', '').replace('.txt', '')
                    
                    # Read data to get statistics
                    autocorr_data = read_autocorr_data(str(autocorr_file))
                    if autocorr_data:
                        time_fs = autocorr_data['time_fs']
                        autocorr = autocorr_data['autocorr']
                        
                        f.write(f"\n{ref_group}:\n")
                        f.write(f"  File: {autocorr_file.name}\n")
                        f.write(f"  Data points: {len(time_fs)}\n")
                        f.write(f"  Time range: {time_fs.min():.3f} - {time_fs.max():.3f} fs\n")
                        f.write(f"  Autocorr range: {autocorr.min():.6f} - {autocorr.max():.6f}\n")
                        f.write(f"  Initial autocorr: {autocorr[0]:.6f}\n")
                        f.write(f"  Final autocorr: {autocorr[-1]:.6f}\n")
                        
                        # Check for sample count file
                        sample_file = autocorr_file.parent / f"{autocorr_file.stem}_sample_counts.txt"
                        if sample_file.exists():
                            try:
                                sample_counts = np.loadtxt(sample_file, dtype=int)
                                f.write(f"  Sample counts: {sample_counts.min()} - {sample_counts.max()} (mean: {sample_counts.mean():.1f})\n")
                            except:
                                f.write(f"  Sample count file found but could not read\n")
                        
                        # Check for IR spectrum file
                        spectrum_file = self.exp_dir / f"ir_spectrum_{ref_group}.txt"
                        if spectrum_file.exists():
                            f.write(f"  IR spectrum file: ir_spectrum_{ref_group}.txt\n")
                    else:
                        f.write(f"\n{ref_group}: Could not read data\n")
            else:
                f.write("No master autocorrelation files found.\n")
                f.write("Make sure to run processing first with --skip_processing=False\n")
        
        logger.info(f"Saved analysis summary: {report_file}")
        return True

    def run_analysis(self, processing_results=None):
        """Run the complete dipole autocorrelation and IR analysis for the single directory."""
        logger.info(f"Starting dipole autocorrelation and IR analysis for: {self.exp_name}")
        
        # Find autocorrelation files
        autocorr_files = self.find_autocorr_files()
        
        if not autocorr_files:
            logger.error("No master autocorrelation files found! Make sure to run processing first.")
            return False
        
        # Create analysis plots
        overview_success = self.create_autocorr_overview_plot(autocorr_files)
        ir_success = self.create_ir_spectra_plot(autocorr_files)
        
        # Create summary report
        summary_success = self.create_summary_report(autocorr_files, processing_results)
        
        if overview_success or ir_success:
            logger.info("Dipole autocorrelation and IR analysis complete!")
            return True
        else:
            logger.error("Failed to create analysis plots")
            return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Process dipole autocorrelation files and generate IR spectra from a single experiment directory')
    parser.add_argument('--exp_dir', required=True, help='Experiment directory to process')
    parser.add_argument('--job_name', default='prod', help='Job name for autocorrelation files (default: prod)')
    parser.add_argument('--max_time', type=float, default=None,
                       help='Maximum time for autocorrelation data (fs) - if not specified, auto-detect from data. Use smaller values for running simulations.')
    parser.add_argument('--autocorr_dt', type=float, default=1.0,
                       help='Timestep for autocorrelation data interpolation (fs)')
    parser.add_argument('--temperature', type=float, default=300.0,
                       help='Temperature for IR spectrum calculation (K)')
    parser.add_argument('--skip_processing', action='store_true',
                       help='Skip master file processing (use existing files)')
    parser.add_argument('--running_sim', action='store_true',
                       help='Optimized settings for running simulations (conservative time detection, better error handling)')
    
    args = parser.parse_args()
    
    start_time = time.time()
    
    # Validate experiment directory
    exp_dir = Path(args.exp_dir).resolve()
    if not exp_dir.exists():
        logger.error(f"Experiment directory does not exist: {exp_dir}")
        return 1
    
    if not exp_dir.is_dir():
        logger.error(f"Path is not a directory: {exp_dir}")
        return 1
    
    logger.info("Starting dipole autocorrelation processing and IR analysis...")
    logger.info(f"Experiment directory: {exp_dir}")
    logger.info(f"Job name: {args.job_name}")
    logger.info(f"Autocorrelation dt: {args.autocorr_dt} fs")
    logger.info(f"Temperature: {args.temperature} K")
    
    if args.running_sim:
        logger.info("Running simulation mode: Using conservative settings for incomplete data")
        if args.max_time is None:
            logger.info("Auto-detecting time range with conservative approach")
        else:
            logger.info(f"Max autocorrelation time: {args.max_time} fs (fixed for running simulation)")
    else:
        if args.max_time is not None:
            logger.info(f"Max autocorrelation time: {args.max_time} fs (fixed)")
        else:
            logger.info("Max autocorrelation time: auto-detect from data")
    
    processing_results = None
    
    # Step 1: Process autocorrelation files to create master files (unless skipped)
    if not args.skip_processing:
        logger.info("\n" + "="*60)
        logger.info("STEP 1: PROCESSING AUTOCORRELATION FILES TO CREATE MASTER FILES")
        logger.info("="*60)
        
        processing_results = process_dipole_autocorr_files(
            exp_dir=exp_dir,
            job_name=args.job_name,
            max_time=args.max_time,
            dt=args.autocorr_dt
        )
        
        if not processing_results['success']:
            logger.error(f"Failed to process autocorrelation files: {processing_results.get('error', 'Unknown error')}")
            return 1
        
        logger.info(f"Successfully processed {processing_results['success_count']} autocorrelation groups")
    else:
        logger.info("Skipping autocorrelation processing (using existing master files)")
    
    # Step 2: Run autocorrelation analysis and generate IR spectra
    logger.info("\n" + "="*60)
    logger.info("STEP 2: RUNNING AUTOCORRELATION ANALYSIS AND GENERATING IR SPECTRA")
    logger.info("="*60)
    
    analyzer = SingleDirectoryDipoleAutocorrAnalyzer(
        exp_dir=exp_dir,
        max_time=args.max_time,
        temperature=args.temperature
    )
    
    analysis_success = analyzer.run_analysis(processing_results)
    
    total_time = time.time() - start_time
    
    if analysis_success:
        logger.info(f"\n" + "="*60)
        logger.info("DIPOLE AUTOCORRELATION AND IR ANALYSIS COMPLETE!")
        logger.info("="*60)
        logger.info(f"Total execution time: {total_time:.2f} seconds")
        logger.info("Generated files:")
        logger.info(f"  - master_dipole_autocorr_ref*.txt (autocorrelation master files)")
        logger.info(f"  - ir_spectrum_ref*.txt (IR spectrum data files)")
        logger.info(f"  - dipole_autocorr_overview.png (individual autocorrelation plots)")
        logger.info(f"  - dipole_autocorr_and_ir_spectra.png (combined autocorr and IR plot)")
        logger.info(f"  - dipole_autocorr_analysis_summary.txt (analysis summary)")
        
        if processing_results and any(g.get('coverage_percent', 100) < 80 for g in processing_results.get('processed_groups', {}).values()):
            logger.info("\nNote: Some reference groups have low data coverage.")
            logger.info("This is normal for running simulations. Rerun when simulation completes for full coverage.")
            logger.info("Use --running_sim flag for optimized processing of incomplete data.")
        
        return 0
    else:
        logger.error("Dipole autocorrelation and IR analysis failed!")
        return 1

if __name__ == "__main__":
    exit(main()) 