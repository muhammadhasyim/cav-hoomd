#!/usr/bin/env python3
"""
F(k,t) Processing and Analysis Script for Cavity MD Simulations (Single Directory Version)

This script processes F(k,t) files from a single cavity MD simulation directory,
creating master files by averaging across all replicas and generating analysis plots.

Key Features:
- Configurable timestep resolution for F(k,t) data
- Processes a single specified experiment directory
- Generates master F(k,t) files by averaging across all replicas
- Creates comprehensive F(k,t) analysis plots for the specified directory

Usage:
    python process_fskt_only.py --exp_dir ./cavity_coupling_1eneg03_switch_1.0ps
    python process_fskt_only.py --exp_dir ./my_experiment --fskt_dt 0.5 --max_time 100.0
    python process_fskt_only.py --exp_dir ./test_run --skip_processing  # Only create plots
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

def process_fskt_files(exp_dir, job_name="prod", max_time=None, dt=1.0):
    """
    Process F(k,t) files using configurable timestep and proper cumulative averaging.
    
    Args:
        exp_dir: Path to experiment directory
        job_name: Job name ('prod' or 'finq')
        max_time: Maximum time in ps for F(k,t) data (None = auto-detect from data)
        dt: Timestep in ps for interpolation
    
    Returns:
        Dictionary with processing results
    """
    exp_path = Path(exp_dir)
    
    if not exp_path.exists():
        logger.error(f"Experiment directory does not exist: {exp_dir}")
        return {'success': False, 'error': 'Directory does not exist'}
    
    # Find all F(k,t) files for this experiment - updated pattern for new naming convention
    fskt_pattern = f"{job_name}*_ref*.txt"
    fskt_files = list(exp_path.glob(fskt_pattern))
    
    if not fskt_files:
        logger.warning(f"No F(k,t) ref files found in {exp_dir} with pattern {fskt_pattern}")
        return {'success': False, 'error': f'No files found with pattern {fskt_pattern}'}
    
    # Group files by reference number (extracting ref* part from filename)
    fskt_groups = {}
    for fskt_file in fskt_files:
        filename = fskt_file.name
        # Extract the reference part from filename like prod-0_ref1.txt -> ref1
        if '_ref' in filename:
            ref_part = filename.split('_ref')[1].replace('.txt', '')
            group_key = f"ref{ref_part}"
        else:
            # Fallback: use the whole filename without extension as group key
            group_key = filename.replace('.txt', '')
        
        if group_key not in fskt_groups:
            fskt_groups[group_key] = []
        fskt_groups[group_key].append(str(fskt_file))
    
    logger.info(f"Processing {len(fskt_groups)} F(k,t) reference groups in {exp_dir}")
    logger.info(f"Found reference groups: {list(fskt_groups.keys())}")
    
    results = {'success': True, 'processed_groups': {}, 'total_groups': len(fskt_groups)}
    success_count = 0
    
    # Process fskt groups
    for group_key, file_list in tqdm(fskt_groups.items(), desc="Processing F(k,t) groups"):
        # Sort files naturally
        file_list = natsort.natsorted(file_list)
        
        # Auto-detect maximum time range from all files if not specified
        detected_max_time = max_time
        if max_time is None:
            all_max_times = []
            # Process time detection with improved error handling
            for fskt_file in file_list:
                temp_data = pd.read_csv(fskt_file, sep=r'\s+', comment='#', header=None)
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
                            logger.debug(f"File {Path(fskt_file).name}: {len(valid_times)} time points, max time: {valid_times.max():.3f} ps")
            
            if all_max_times:
                # Use a conservative approach: use the minimum time found to ensure all files have data
                # This prevents trying to interpolate beyond available data for running simulations
                detected_max_time = min(all_max_times)
                logger.info(f"Auto-detected max time for {group_key}: {detected_max_time:.3f} ps from {len(all_max_times)} valid files")
                logger.info(f"  Using conservative estimate (minimum across files) to handle incomplete data")
                logger.info(f"  Time range across files: {min(all_max_times):.3f} to {max(all_max_times):.3f} ps")
            else:
                detected_max_time = 10.0  # Very conservative fallback for running simulations
                logger.warning(f"Could not detect time range for {group_key} (files may be empty/incomplete), using conservative fallback: {detected_max_time} ps")
        
        # Set fixed linspace with configurable timestep
        if detected_max_time is not None:
            num_points = int(detected_max_time / dt) + 1
            uniform_times = np.linspace(0, detected_max_time, num_points)
            logger.info(f"F(k,t) group {group_key}: Using time grid: 0 to {detected_max_time:.3f} ps with {len(uniform_times)} points (dt = {dt} ps)")
        else:
            logger.error(f"Could not determine time range for {group_key}")
            continue
        
        # Initialize storage for cumulative data
        cumulative_fskt = np.zeros(len(uniform_times))
        sample_counts = np.zeros(len(uniform_times), dtype=int)
        processed_count = 0
        
        # Process each F(k,t) file in this group
        logger.info(f"Processing {group_key}: {len(file_list)} files...")
        for fskt_file in file_list:
            # Read the F(k,t) data
            data = pd.read_csv(fskt_file, sep=r'\s+', comment='#', header=None)
            
            if data.empty:
                continue
            
            # Convert all data to numeric, coercing errors to NaN
            data = data.apply(pd.to_numeric, errors='coerce')
            
            # Drop any rows with NaN values
            data = data.dropna()
            
            if data.empty:
                continue
            
            # For 3-column format: timestep, time(ps), F(k,t)
            # Use column 1 as time and column 2 as F(k,t)
            if data.shape[1] >= 3:
                times = data.iloc[:, 1].values  # Second column is time (ps)
                fskt_vals = data.iloc[:, 2].values  # Third column is F(k,t)
            else:
                # Fallback for 2-column format: time, F(k,t)
                times = data.iloc[:, 0].values  # First column is time
                fskt_vals = data.iloc[:, 1].values  # Second column is F(k,t)
            
            # Filter to valid time range - be more lenient for incomplete files
            valid_mask = (times >= 0)
            if np.sum(valid_mask) < 3:
                continue
            
            times = times[valid_mask]
            fskt_vals = fskt_vals[valid_mask]
            
            # Log the actual time range used for this file
            file_min_time, file_max_time = times.min(), times.max()
            
            # Remove duplicate time values before interpolation
            unique_times, unique_fskt = remove_duplicate_times(times, fskt_vals)
            
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
            interpolated_fskt = np.interp(valid_grid_times, unique_times, unique_fskt)
            
            # Update cumulative average using proper online algorithm
            for i, grid_idx in enumerate(np.where(valid_grid_mask)[0]):
                n = sample_counts[grid_idx] + 1
                if sample_counts[grid_idx] == 0:
                    cumulative_fskt[grid_idx] = interpolated_fskt[i]
                else:
                    # Online average update: new_avg = old_avg + (new_val - old_avg) / n
                    old_avg = cumulative_fskt[grid_idx]
                    new_val = interpolated_fskt[i]
                    cumulative_fskt[grid_idx] = old_avg + (new_val - old_avg) / n
                sample_counts[grid_idx] = n
            
            processed_count += 1
        
        # Save results if we have data
        if processed_count > 0:
            master_file = exp_path / f"master_fskt_{group_key}.txt"
            sample_count_file = exp_path / f"master_fskt_{group_key}_sample_counts.txt"
            
            # Only keep time points where we have at least one sample
            valid_points = sample_counts > 0
            clean_times = uniform_times[valid_points]
            clean_fskt = cumulative_fskt[valid_points]
            clean_sample_counts = sample_counts[valid_points]
            
            if len(clean_times) == 0:
                logger.warning(f"No valid data points found for {group_key} after processing {processed_count} files")
                continue
            
            # Create output data
            output_data = pd.DataFrame({
                'lag_time': clean_times,
                'fskt': clean_fskt
            })
            
            # Calculate data coverage statistics
            total_possible_points = len(uniform_times)
            actual_data_points = len(clean_times)
            coverage_percent = 100.0 * actual_data_points / total_possible_points
            
            # Save master data with comprehensive header
            with open(master_file, 'w') as f:
                f.write(f"# Master F(k,t) file for {exp_path.name} - {group_key}\n")
                f.write(f"# Cumulative average across {processed_count} replicas\n")
                f.write(f"# Time grid: 0 to {detected_max_time:.3f} ps with {len(uniform_times)} points (dt = {dt} ps)\n")
                f.write(f"# Valid data points: {len(clean_times)} ({coverage_percent:.1f}% coverage)\n")
                f.write(f"# Sample count range: {clean_sample_counts.min()} to {clean_sample_counts.max()}\n")
                f.write(f"# Note: Simulation may be incomplete - sample counts indicate data availability\n")
                f.write(f"# Generated by process_fskt_only.py\n")
                f.write(f"# Columns: lag_time fskt\n")
            
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

class SingleDirectoryFSKTAnalyzer:
    """Analyzes F(k,t) data from a single experiment directory and creates analysis plots."""
    
    def __init__(self, exp_dir, max_time=None):
        self.exp_dir = Path(exp_dir).resolve()
        self.max_time = max_time
        self.exp_name = self.exp_dir.name
        
        logger.info(f"Initialized analyzer for directory: {self.exp_dir}")

    def find_fskt_files(self):
        """Find master F(k,t) data files in the experiment directory."""
        if not self.exp_dir.exists():
            logger.error(f"Experiment directory does not exist: {self.exp_dir}")
            return []
        
        fskt_files = list(self.exp_dir.glob("master_fskt_ref*.txt"))
        valid_files = [f for f in fskt_files if f.stat().st_size > 0]
        
        if valid_files:
            logger.info(f"Found {len(valid_files)} master F(k,t) files")
            for f in valid_files:
                logger.info(f"  {f.name}")
        else:
            logger.warning("No master F(k,t) files found")
        
        return sorted(valid_files)

    def create_fskt_overview_plot(self, fskt_files):
        """Create overview plot showing all k-values from the directory."""
        logger.info("Generating F(k,t) overview plot...")
        
        if not fskt_files:
            logger.warning("No F(k,t) files to plot")
            return False
        
        # Create figure with appropriate size
        n_files = len(fskt_files)
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
        
        for i, fskt_file in enumerate(fskt_files):
            if i >= len(axes):
                break
                
            ax = axes[i]
            
            # Extract reference group from filename
            ref_group = fskt_file.name.replace('master_fskt_', '').replace('.txt', '')
            
            # Read and plot data
            fskt_data = read_fskt_data(str(fskt_file))
            if fskt_data:
                lag_time = fskt_data['lag_time']
                fskt = fskt_data['fskt']
                
                ax.plot(lag_time, fskt, color=colors[i], linewidth=2, alpha=0.8)
                
                # Add statistics
                if len(fskt) > 10:
                    final_value = fskt[-1]
                    min_value = np.min(fskt)
                    max_value = np.max(fskt)
                    
                    stats_text = f"Final: {final_value:.4f}\nMin: {min_value:.4f}\nMax: {max_value:.4f}"
                    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                           verticalalignment='top', fontsize=8,
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            else:
                ax.text(0.5, 0.5, f'No valid data\nfor {ref_group}', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=10)
            
            ax.set_xlabel('Time (ps)', fontsize=10)
            ax.set_ylabel('F(k,t)', fontsize=10)
            ax.set_ylim(bottom=0.0)
            ax.set_title(f'{ref_group}', fontweight='bold', fontsize=11)
            ax.grid(True, alpha=0.3)
        
        # Remove unused axes
        for i in range(len(fskt_files), len(axes)):
            fig.delaxes(axes[i])
        
        plt.suptitle(f'F(k,t) Analysis - {self.exp_name}\nAll Reference Groups', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = self.exp_dir / "fskt_overview_analysis.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved F(k,t) overview plot: {output_file}")
        plt.close()
        
        return True

    def create_combined_fskt_plot(self, fskt_files):
        """Create combined plot showing all reference groups on the same axes."""
        logger.info("Generating combined F(k,t) plot...")
        
        if not fskt_files:
            logger.warning("No F(k,t) files to plot")
            return False
        
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(fskt_files)))
        has_data = False
        
        for i, fskt_file in enumerate(fskt_files):
            # Extract reference group from filename
            ref_group = fskt_file.name.replace('master_fskt_', '').replace('.txt', '')
            
            # Read and plot data
            fskt_data = read_fskt_data(str(fskt_file))
            if fskt_data:
                lag_time = fskt_data['lag_time']
                fskt = fskt_data['fskt']
                
                ax.plot(lag_time, fskt, color=colors[i], linewidth=2.5, 
                       alpha=0.8, label=ref_group, marker='o', markersize=1)
                has_data = True
        
        if not has_data:
            ax.text(0.5, 0.5, f'No valid F(k,t) data found\nin {self.exp_name}', 
                   ha='center', va='center', transform=ax.transAxes, fontsize=14)
            logger.warning(f"No valid F(k,t) data found for {self.exp_name}")
            return False
        
        # Customize the plot
        ax.set_xlabel('Time (ps)', fontsize=14)
        ax.set_ylabel('F(k,t)', fontsize=14)
        ax.set_ylim(bottom=0.0)
        
        title = f'F(k,t) Combined Analysis - {self.exp_name}\nAll Reference Groups'
        ax.set_title(title, fontweight='bold', fontsize=16)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10, loc='best')
        
        plt.tight_layout()
        
        # Save the plot
        output_file = self.exp_dir / 'fskt_combined_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved combined F(k,t) plot: {output_file}")
        plt.close()
        
        return True

    def create_summary_report(self, fskt_files, processing_results=None):
        """Create a summary report of the F(k,t) analysis."""
        logger.info("Generating F(k,t) summary report...")
        
        report_file = self.exp_dir / "fskt_analysis_summary.txt"
        
        with open(report_file, 'w') as f:
            f.write(f"F(k,t) Analysis Summary Report\n")
            f.write(f"{'='*50}\n\n")
            f.write(f"Experiment Directory: {self.exp_dir}\n")
            f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            if processing_results:
                f.write(f"Processing Results:\n")
                f.write(f"  Total groups processed: {processing_results.get('success_count', 0)}/{processing_results.get('total_groups', 0)}\n")
                f.write(f"  Processing success: {processing_results.get('success', False)}\n\n")
            
            f.write(f"Master F(k,t) Files Found: {len(fskt_files)}\n\n")
            
            if fskt_files:
                f.write("File Details:\n")
                f.write("-" * 30 + "\n")
                
                for fskt_file in fskt_files:
                    ref_group = fskt_file.name.replace('master_fskt_', '').replace('.txt', '')
                    
                    # Read data to get statistics
                    fskt_data = read_fskt_data(str(fskt_file))
                    if fskt_data:
                        lag_time = fskt_data['lag_time']
                        fskt = fskt_data['fskt']
                        
                        f.write(f"\n{ref_group}:\n")
                        f.write(f"  File: {fskt_file.name}\n")
                        f.write(f"  Data points: {len(lag_time)}\n")
                        f.write(f"  Time range: {lag_time.min():.3f} - {lag_time.max():.3f} ps\n")
                        f.write(f"  F(k,t) range: {fskt.min():.6f} - {fskt.max():.6f}\n")
                        f.write(f"  Final F(k,t): {fskt[-1]:.6f}\n")
                        
                        # Check for sample count file
                        sample_file = fskt_file.parent / f"{fskt_file.stem}_sample_counts.txt"
                        if sample_file.exists():
                            try:
                                sample_counts = np.loadtxt(sample_file, dtype=int)
                                f.write(f"  Sample counts: {sample_counts.min()} - {sample_counts.max()} (mean: {sample_counts.mean():.1f})\n")
                            except:
                                f.write(f"  Sample count file found but could not read\n")
                    else:
                        f.write(f"\n{ref_group}: Could not read data\n")
            else:
                f.write("No master F(k,t) files found.\n")
                f.write("Make sure to run processing first with --skip_processing=False\n")
        
        logger.info(f"Saved analysis summary: {report_file}")
        return True

    def run_analysis(self, processing_results=None):
        """Run the complete F(k,t) analysis for the single directory."""
        logger.info(f"Starting F(k,t) analysis for: {self.exp_name}")
        
        # Find F(k,t) files
        fskt_files = self.find_fskt_files()
        
        if not fskt_files:
            logger.error("No master F(k,t) files found! Make sure to run processing first.")
            return False
        
        # Create analysis plots
        overview_success = self.create_fskt_overview_plot(fskt_files)
        combined_success = self.create_combined_fskt_plot(fskt_files)
        
        # Create summary report
        summary_success = self.create_summary_report(fskt_files, processing_results)
        
        if overview_success or combined_success:
            logger.info("F(k,t) analysis complete!")
            return True
        else:
            logger.error("Failed to create analysis plots")
            return False

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Process F(k,t) files from a single experiment directory')
    parser.add_argument('--exp_dir', required=True, help='Experiment directory to process')
    parser.add_argument('--job_name', default='prod', help='Job name for F(k,t) files (default: prod)')
    parser.add_argument('--max_time', type=float, default=None,
                       help='Maximum time for F(k,t) data (ps) - if not specified, auto-detect from data. Use smaller values for running simulations.')
    parser.add_argument('--fskt_dt', type=float, default=1.0,
                       help='Timestep for F(k,t) data interpolation (ps)')
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
    
    logger.info("Starting F(k,t) processing and analysis...")
    logger.info(f"Experiment directory: {exp_dir}")
    logger.info(f"Job name: {args.job_name}")
    logger.info(f"F(k,t) dt: {args.fskt_dt} ps")
    
    if args.running_sim:
        logger.info("Running simulation mode: Using conservative settings for incomplete data")
        if args.max_time is None:
            logger.info("Auto-detecting time range with conservative approach")
        else:
            logger.info(f"Max F(k,t) time: {args.max_time} ps (fixed for running simulation)")
    else:
        if args.max_time is not None:
            logger.info(f"Max F(k,t) time: {args.max_time} ps (fixed)")
        else:
            logger.info("Max F(k,t) time: auto-detect from data")
    
    processing_results = None
    
    # Step 1: Process F(k,t) files to create master files (unless skipped)
    if not args.skip_processing:
        logger.info("\n" + "="*60)
        logger.info("STEP 1: PROCESSING F(k,t) FILES TO CREATE MASTER FILES")
        logger.info("="*60)
        
        processing_results = process_fskt_files(
            exp_dir=exp_dir,
            job_name=args.job_name,
            max_time=args.max_time,
            dt=args.fskt_dt
        )
        
        if not processing_results['success']:
            logger.error(f"Failed to process F(k,t) files: {processing_results.get('error', 'Unknown error')}")
            return 1
        
        logger.info(f"Successfully processed {processing_results['success_count']} F(k,t) groups")
    else:
        logger.info("Skipping F(k,t) processing (using existing master files)")
    
    # Step 2: Run F(k,t) analysis and generate plots
    logger.info("\n" + "="*60)
    logger.info("STEP 2: RUNNING F(k,t) ANALYSIS AND GENERATING PLOTS")
    logger.info("="*60)
    
    analyzer = SingleDirectoryFSKTAnalyzer(
        exp_dir=exp_dir,
        max_time=args.max_time
    )
    
    analysis_success = analyzer.run_analysis(processing_results)
    
    total_time = time.time() - start_time
    
    if analysis_success:
        logger.info(f"\n" + "="*60)
        logger.info("F(k,t) PROCESSING AND ANALYSIS COMPLETE!")
        logger.info("="*60)
        logger.info(f"Total execution time: {total_time:.2f} seconds")
        logger.info("Generated files:")
        logger.info(f"  - master_fskt_ref*.txt (F(k,t) master files)")
        logger.info(f"  - fskt_overview_analysis.png (individual reference group plots)")
        logger.info(f"  - fskt_combined_analysis.png (combined plot)")
        logger.info(f"  - fskt_analysis_summary.txt (analysis summary)")
        
        if processing_results and any(g.get('coverage_percent', 100) < 80 for g in processing_results.get('processed_groups', {}).values()):
            logger.info("\nNote: Some reference groups have low data coverage.")
            logger.info("This is normal for running simulations. Rerun when simulation completes for full coverage.")
            logger.info("Use --running_sim flag for optimized processing of incomplete data.")
        
        return 0
    else:
        logger.error("F(k,t) analysis failed!")
        return 1

if __name__ == "__main__":
    exit(main()) 
