#!/usr/bin/env python3
"""
Simple analysis script for coupling experiment data.
This script analyzes data where all replica files are in a single directory.
"""

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

def analyze_coupling_experiment(coupling_dir, max_time_ps=100.0):
    """Analyze a single coupling experiment directory."""
    coupling_dir = Path(coupling_dir)
    
    print(f"Analyzing coupling experiment: {coupling_dir}")
    
    # Find all energy files
    energy_files = list(coupling_dir.glob("prod-*_energy_contributions.txt"))
    print(f"Found {len(energy_files)} energy files")
    
    # Find all F(k,t) files
    fskt_files = list(coupling_dir.glob("prod-*_fskt_*.txt"))
    print(f"Found {len(fskt_files)} F(k,t) files")
    
    if len(energy_files) == 0:
        print("No energy files found!")
        return False
    
    # Process energy files
    print("Processing energy files...")
    all_energy_data = []
    
    for i, energy_file in enumerate(sorted(energy_files)[:20]):  # Limit to first 20 for speed
        try:
            print(f"  Processing {energy_file.name}...")
            
            # Read the file more carefully
            with open(energy_file, 'r') as f:
                lines = f.readlines()
            
            # Skip comment lines and find data
            data_lines = [line.strip() for line in lines if not line.startswith('#') and line.strip()]
            
            if not data_lines:
                print(f"    Warning: No data lines found")
                continue
            
            # Parse header from the comment lines
            header_line = None
            for line in lines:
                if line.startswith('# time(ps)'):
                    header_line = line[2:].strip()  # Remove '# '
                    break
            
            if not header_line:
                print(f"    Warning: No header found")
                continue
            
            # Split header - handle long headers that may wrap
            headers = header_line.split()
            
            # Parse data lines
            data_rows = []
            for line in data_lines:
                try:
                    values = [float(x) for x in line.split()]
                    if len(values) >= len(headers):
                        data_rows.append(values[:len(headers)])  # Take only as many as headers
                    elif len(values) >= 3:  # At least time and some data
                        # Pad with NaN if needed
                        padded_values = values + [np.nan] * (len(headers) - len(values))
                        data_rows.append(padded_values)
                except ValueError:
                    continue
            
            if not data_rows:
                print(f"    Warning: No valid data rows")
                continue
            
            # Create DataFrame
            df = pd.DataFrame(data_rows, columns=headers[:len(data_rows[0])])
            
            if df.empty or 'time(ps)' not in df.columns:
                print(f"    Warning: Invalid DataFrame structure")
                continue
            
            # Rename time column for consistency
            df = df.rename(columns={'time(ps)': 'time'})
            
            # Filter by max time
            df_filtered = df[df['time'] <= max_time_ps].copy()
            df_filtered['replica'] = i
            all_energy_data.append(df_filtered)
            print(f"    Loaded {len(df_filtered)} time points, columns: {list(df_filtered.columns[:5])}...")
            
        except Exception as e:
            print(f"    Error: {e}")
            import traceback
            traceback.print_exc()
    
    if not all_energy_data:
        print("No valid energy data found!")
        return False
    
    # Combine all energy data
    combined_energy = pd.concat(all_energy_data, ignore_index=True)
    print(f"Combined energy data: {len(combined_energy)} total points")
    print(f"Available columns: {list(combined_energy.columns)}")
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Coupling Experiment Analysis: {coupling_dir.name}', fontsize=14)
    
    # Find appropriate column names (handle variations)
    total_energy_col = None
    temperature_col = None
    cavity_energy_col = None
    
    for col in combined_energy.columns:
        if 'total' in col.lower() and 'energy' in col.lower():
            total_energy_col = col
        elif 'temperature' in col.lower():
            temperature_col = col
        elif 'cavity' in col.lower() and 'energy' in col.lower():
            cavity_energy_col = col
    
    if not total_energy_col:
        # Try alternative column names
        if 'universe_total_energy' in combined_energy.columns:
            total_energy_col = 'universe_total_energy'
        elif 'total_potential_energy' in combined_energy.columns:
            total_energy_col = 'total_potential_energy'
    
    # Plot 1: Total energy vs time
    if total_energy_col and total_energy_col in combined_energy.columns:
        for replica_id in combined_energy['replica'].unique():
            replica_data = combined_energy[combined_energy['replica'] == replica_id]
            axes[0,0].plot(replica_data['time'], replica_data[total_energy_col], 
                          alpha=0.7, linewidth=0.8, label=f'Replica {replica_id}')
        axes[0,0].set_xlabel('Time (ps)')
        axes[0,0].set_ylabel('Total Energy')
        axes[0,0].set_title(f'Total Energy vs Time ({total_energy_col})')
        axes[0,0].grid(True, alpha=0.3)
    else:
        axes[0,0].text(0.5, 0.5, 'Total energy column not found', 
                       ha='center', va='center', transform=axes[0,0].transAxes)
    
    # Plot 2: Potential energy vs time
    if 'total_potential_energy' in combined_energy.columns:
        for replica_id in combined_energy['replica'].unique():
            replica_data = combined_energy[combined_energy['replica'] == replica_id]
            axes[0,1].plot(replica_data['time'], replica_data['total_potential_energy'], 
                          alpha=0.7, linewidth=0.8, label=f'Replica {replica_id}')
        axes[0,1].set_xlabel('Time (ps)')
        axes[0,1].set_ylabel('Potential Energy')
        axes[0,1].set_title('Potential Energy vs Time')
        axes[0,1].grid(True, alpha=0.3)
    else:
        axes[0,1].text(0.5, 0.5, 'Potential energy column not found', 
                       ha='center', va='center', transform=axes[0,1].transAxes)
    
    # Plot 3: Cavity energy vs time
    cavity_cols = [col for col in combined_energy.columns if 'cavity' in col.lower()]
    if cavity_cols:
        cavity_col = cavity_cols[0]  # Use first cavity column found
        for replica_id in combined_energy['replica'].unique():
            replica_data = combined_energy[combined_energy['replica'] == replica_id]
            axes[1,0].plot(replica_data['time'], replica_data[cavity_col], 
                          alpha=0.7, linewidth=0.8, label=f'Replica {replica_id}')
        axes[1,0].set_xlabel('Time (ps)')
        axes[1,0].set_ylabel('Cavity Energy')
        axes[1,0].set_title(f'Cavity Energy vs Time ({cavity_col})')
        axes[1,0].grid(True, alpha=0.3)
    else:
        axes[1,0].text(0.5, 0.5, 'No cavity energy columns found', 
                       ha='center', va='center', transform=axes[1,0].transAxes)
    
    # Plot 4: Energy statistics
    if total_energy_col and total_energy_col in combined_energy.columns:
        # Calculate mean and std for each time point
        energy_stats = combined_energy.groupby('time')[total_energy_col].agg(['mean', 'std']).reset_index()
        axes[1,1].plot(energy_stats['time'], energy_stats['mean'], 'b-', linewidth=2, label='Mean')
        axes[1,1].fill_between(energy_stats['time'], 
                              energy_stats['mean'] - energy_stats['std'],
                              energy_stats['mean'] + energy_stats['std'],
                              alpha=0.3, label='±1 σ')
        axes[1,1].set_xlabel('Time (ps)')
        axes[1,1].set_ylabel('Total Energy')
        axes[1,1].set_title('Energy Statistics')
        axes[1,1].grid(True, alpha=0.3)
        axes[1,1].legend()
    else:
        axes[1,1].text(0.5, 0.5, 'Cannot calculate energy statistics', 
                       ha='center', va='center', transform=axes[1,1].transAxes)
    
    plt.tight_layout()
    
    # Save plot
    output_file = coupling_dir / f"{coupling_dir.name}_analysis.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved analysis plot: {output_file}")
    
    # Save summary statistics
    summary_file = coupling_dir / f"{coupling_dir.name}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Coupling Experiment Analysis Summary\n")
        f.write(f"====================================\n")
        f.write(f"Experiment: {coupling_dir.name}\n")
        f.write(f"Number of replicas analyzed: {len(all_energy_data)}\n")
        f.write(f"Time range: 0 - {max_time_ps} ps\n")
        f.write(f"Total data points: {len(combined_energy)}\n")
        f.write(f"Available columns: {list(combined_energy.columns)}\n\n")
        
        if total_energy_col and total_energy_col in combined_energy.columns:
            f.write(f"Total Energy Statistics ({total_energy_col}):\n")
            f.write(f"  Mean: {combined_energy[total_energy_col].mean():.6f}\n")
            f.write(f"  Std:  {combined_energy[total_energy_col].std():.6f}\n")
            f.write(f"  Min:  {combined_energy[total_energy_col].min():.6f}\n")
            f.write(f"  Max:  {combined_energy[total_energy_col].max():.6f}\n\n")
        
        if 'total_potential_energy' in combined_energy.columns:
            f.write(f"Potential Energy Statistics:\n")
            f.write(f"  Mean: {combined_energy['total_potential_energy'].mean():.6f}\n")
            f.write(f"  Std:  {combined_energy['total_potential_energy'].std():.6f}\n")
            f.write(f"  Min:  {combined_energy['total_potential_energy'].min():.6f}\n")
            f.write(f"  Max:  {combined_energy['total_potential_energy'].max():.6f}\n\n")
    
    print(f"Saved summary: {summary_file}")
    return True

def main():
    parser = argparse.ArgumentParser(description="Analyze coupling experiment data")
    parser.add_argument('coupling_dirs', nargs='+', help='Coupling experiment directories to analyze')
    parser.add_argument('--max-time', type=float, default=100.0, help='Maximum time to analyze (ps)')
    
    args = parser.parse_args()
    
    success_count = 0
    for coupling_dir in args.coupling_dirs:
        print(f"\n{'='*60}")
        if analyze_coupling_experiment(coupling_dir, args.max_time):
            success_count += 1
        print(f"{'='*60}")
    
    print(f"\nAnalyzed {success_count}/{len(args.coupling_dirs)} experiments successfully!")

if __name__ == '__main__':
    main() 