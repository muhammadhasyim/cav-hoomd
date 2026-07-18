#!/usr/bin/env python3
"""
Energy Analysis and Plotting Script for Cavity MD Simulations

This script reads the energy contributions file from cavity MD simulations
and creates comprehensive plots of different energy components including:
- Cavity total energy (kinetic + potential)
- Molecular total energy (kinetic + potential) 
- Individual potential and kinetic energies for both molecular and cavity systems
- Thermostat reservoir energies (Langevin, Bussi, or MTTK)
- Universe total energy (system + reservoir)

Temperature Analysis Features:
- Molecular system temperature time series (from kinetic energy file)
- Cavity particle temperature time series (from cavity mode file)
- Temperature statistics and comparison with target temperature
- Velocity distribution analysis with Maxwell-Boltzmann comparisons
- Energy values expressed in both atomic units and thermal energy units (kT)

The script generates multiple plots:
1. Main components plot (cavity total, molecular total, reservoirs, universe total)
2. Potential/kinetic components plot (individual P.E. and K.E. for molecular and cavity)
3. Overview plot (all main components together)
4. Temperature time series plots (if temperature data files are available)
5. Velocity histogram plots (if velocity analysis is requested)

Usage Examples:
    # Basic energy analysis at default temperature (100 K)
    python plot_energy_analysis.py --file prod-0_energy_contributions.txt
    
    # Analysis with custom temperature and time limit
    python plot_energy_analysis.py --file prod-0_energy_contributions.txt --target_temp 200 --max_time 100
    
    # Temperature-only analysis with statistics
    python plot_energy_analysis.py --file prod-0_energy_contributions.txt --temp_only --target_temp 77 --stats
    
    # Full analysis including velocity distributions
    python plot_energy_analysis.py --file prod-0_energy_contributions.txt --target_temp 298.15 --velocity_analysis --stats
    
    # Cavity-only velocity analysis
    python plot_energy_analysis.py --file prod-0_energy_contributions.txt --velocity_analysis --cavity_velocity_only
"""

import numpy as np
# Set matplotlib backend to non-interactive before importing pyplot
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless systems
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os
from pathlib import Path
import gsd.hoomd

def read_energy_data(filename):
    """
    Read energy contributions data from the simulation output file.
    Handles master files which may have mixed data types and inconsistent structure.
    
    Args:
        filename (str): Path to the energy contributions file
        
    Returns:
        pandas.DataFrame: DataFrame containing all energy data
    """
    try:
        # First, read the header line to get column names
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Find the header line (starts with # and contains column names)
        header_line = None
        for line in lines:
            if line.startswith('# time(ps)'):
                header_line = line.strip()
                break
        
        if header_line is None:
            print("Warning: Could not find header line, using standard column names")
            # Read without header and assign standard names based on replica file structure
            data = pd.read_csv(filename, sep=r'\s+', comment='#', header=None, dtype=str)
            
            # Standard column names from replica files (17 columns)
            standard_columns = [
                'time(ps)', 'timestep', 'harmonic_energy', 'lj_energy', 'coulomb_short_energy', 
                'coulomb_long_energy', 'cavity_harmonic_energy', 'cavity_coupling_energy', 
                'cavity_dipole_self_energy', 'cavity_total_potential_energy', 'molecular_kinetic_energy', 
                'cavity_mode_kinetic_energy', 'molecular_reservoir_energy', 'cavity_reservoir_energy',
                'total_potential_energy', 'total_kinetic_energy', 'total_energy'
            ]
            
            # Dynamically assign column names based on actual number of columns
            n_cols = data.shape[1]
            
            # Use as many standard names as we have columns, then generic names for extra columns
            if n_cols <= len(standard_columns):
                column_names = standard_columns[:n_cols]
            else:
                column_names = standard_columns + [f'extra_col_{i}' for i in range(len(standard_columns), n_cols)]
            
            data.columns = column_names
            print(f"  Assigned {n_cols} column names from standard structure")
        else:
            # Extract column names from header line
            column_names = header_line[2:].split()  # Remove '# ' and split
            
            # Read the data without the header, force string type first to handle mixed types
            data = pd.read_csv(filename, sep=r'\s+', comment='#', header=None, dtype=str)
            
            # Handle case where header has more/fewer names than actual columns
            n_cols = data.shape[1]
            if len(column_names) != n_cols:
                print(f"Warning: Header has {len(column_names)} names but data has {n_cols} columns")
                if n_cols < len(column_names):
                    column_names = column_names[:n_cols]
                else:
                    column_names.extend([f'col_{i}' for i in range(len(column_names), n_cols)])
            
            data.columns = column_names
        
        # Convert columns to numeric, handling any mixed types or non-numeric values
        for col in data.columns:
            data[col] = pd.to_numeric(data[col], errors='coerce')
        
        # Drop rows with NaN values that resulted from non-numeric conversion
        initial_rows = len(data)
        data = data.dropna()
        final_rows = len(data)
        
        if initial_rows != final_rows:
            print(f"Warning: Dropped {initial_rows - final_rows} rows with non-numeric data")
        
        # Check for minimum required columns (more flexible)
        required_columns = ['time(ps)', 'molecular_kinetic_energy', 'cavity_mode_kinetic_energy']
        missing_required = [col for col in required_columns if col not in data.columns]
        
        if missing_required:
            print(f"Warning: Missing critical columns: {missing_required}")
            # Still try to proceed if we have some energy data
            if 'molecular_kinetic_energy' not in data.columns and 'cavity_mode_kinetic_energy' not in data.columns:
                print("Error: No kinetic energy columns found")
                return None
        
        # Add total_potential_energy if missing but we have components
        if 'total_potential_energy' not in data.columns:
            potential_components = ['harmonic_energy', 'lj_energy', 'coulomb_short_energy', 
                                  'coulomb_long_energy', 'cavity_total_potential_energy']
            available_components = [col for col in potential_components if col in data.columns]
            
            if len(available_components) >= 3:  # Need at least some components
                print("Warning: total_potential_energy missing, calculating from components")
                data['total_potential_energy'] = sum(data[col] for col in available_components)
        
        # Add total_kinetic_energy if missing but we have components
        if 'total_kinetic_energy' not in data.columns:
            kinetic_components = ['molecular_kinetic_energy', 'cavity_mode_kinetic_energy']
            available_ke_components = [col for col in kinetic_components if col in data.columns]
            
            if available_ke_components:
                print("Warning: total_kinetic_energy missing, calculating from components")
                data['total_kinetic_energy'] = sum(data[col] for col in available_ke_components)
        
        print(f"Successfully loaded energy data from {filename}")
        print(f"Data shape: {data.shape}")
        print(f"Columns: {list(data.columns)}")
        if len(data) > 0:
            print(f"Time range: {data.iloc[0, 0]:.3f} - {data.iloc[-1, 0]:.3f} ps")
        return data
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        return None

def calculate_derived_energies(df):
    """
    Calculate the main energy components from the raw data, including individual
    potential and kinetic energies for molecular and cavity systems.
    Handles cases where total_energy column might be missing.
    
    Args:
        df (pd.DataFrame): Raw energy data
        
    Returns:
        dict: Dictionary containing all energy components including individual P.E. and K.E.
    """
    # Molecular potential energy components
    molecular_potential = (df['harmonic_energy'] + df['lj_energy'] + 
                          df['coulomb_short_energy'] + df['coulomb_long_energy'])
    
    # Molecular kinetic energy
    molecular_kinetic = df['molecular_kinetic_energy']
    
    # Molecular total energy (kinetic + potential)
    molecular_total = molecular_potential + molecular_kinetic
    
    # Cavity potential energy 
    cavity_potential = df['cavity_total_potential_energy']
    
    # Cavity kinetic energy
    cavity_kinetic = df['cavity_mode_kinetic_energy']
    
    # Cavity total energy (kinetic + potential)  
    cavity_total = cavity_kinetic + cavity_potential
    
    # Separate thermostat energy components (both Langevin reservoirs)
    molecular_reservoir = df['molecular_reservoir_energy'] if 'molecular_reservoir_energy' in df.columns else pd.Series([0.0] * len(df))
    cavity_reservoir = df['cavity_reservoir_energy'] if 'cavity_reservoir_energy' in df.columns else pd.Series([0.0] * len(df))
    
    # Total thermostat energy (for universe total calculation)
    total_thermostat_energy = molecular_reservoir + cavity_reservoir
    
    # Calculate or get total energy (system total)
    if 'total_energy' in df.columns:
        system_total_energy = df['total_energy']
    else:
        # Calculate total energy from potential + kinetic if column is missing
        print("Warning: 'total_energy' column not found, calculating from total_potential_energy + total_kinetic_energy")
        system_total_energy = df['total_potential_energy'] + df['total_kinetic_energy']
    
    # Universe total energy (system + all thermostat components)
    universe_total = system_total_energy + total_thermostat_energy
    
    return {
        # Individual molecular components
        'molecular_potential_energy': molecular_potential,
        'molecular_kinetic_energy': molecular_kinetic,
        'molecular_total_energy': molecular_total,
        
        # Individual cavity components
        'cavity_potential_energy': cavity_potential,
        'cavity_kinetic_energy': cavity_kinetic,
        'cavity_total_energy': cavity_total,
        
        # Thermostat components
        'molecular_reservoir_energy': molecular_reservoir,
        'cavity_reservoir_energy': cavity_reservoir,
        
        # Universe total
        'universe_total_energy': universe_total
    }

def create_energy_plots(df, derived_energies, output_prefix='energy_analysis', max_time=None, target_temp=100.0):
    """
    Create comprehensive energy plots including individual potential and kinetic energies.
    
    Args:
        df (pd.DataFrame): Raw energy data
        derived_energies (dict): Dictionary containing all energy components
        output_prefix (str): Prefix for output files
        max_time (float): Maximum time to plot (None for all data)
        target_temp (float): Target temperature in Kelvin for plot context
    """
    # Get time array from the dataframe
    time = df.iloc[:, 0].values  # First column is time
    
    # Apply time filter if specified
    if max_time is not None:
        mask = time <= max_time
        time = time[mask]
        # Apply mask to all energy arrays and reset index
        for key in derived_energies:
            derived_energies[key] = derived_energies[key][mask].reset_index(drop=True)
    
    # Create the main plot with the original five energy components
    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(20, 12))
    fig.suptitle(f'Energy Analysis: Main Components (T = {target_temp:.1f} K)', fontsize=16, fontweight='bold')
    
    # 1. Cavity Total Energy
    ax1.plot(time, derived_energies['cavity_total_energy'], 'b-', linewidth=2, label='Cavity Total')
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Energy (a.u.)')
    ax1.set_title('Cavity Total Energy\n(Kinetic + Potential)', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # 2. Molecular Total Energy
    ax2.plot(time, derived_energies['molecular_total_energy'], 'r-', linewidth=2, label='Molecular Total')
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Energy (a.u.)')
    ax2.set_title('Molecular Total Energy\n(Kinetic + Potential)', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # 3. Molecular Reservoir Energy
    ax3.plot(time, derived_energies['molecular_reservoir_energy'], 'g-', linewidth=2, label='Molecular Reservoir')
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('Energy (a.u.)')
    ax3.set_title('Molecular Thermostat Reservoir\n(Molecular Heat Exchange)', fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # 4. Cavity Reservoir Energy
    ax4.plot(time, derived_energies['cavity_reservoir_energy'], 'orange', linewidth=2, label='Cavity Reservoir')
    ax4.set_xlabel('Time (ps)')
    ax4.set_ylabel('Energy (a.u.)')
    ax4.set_title('Cavity Thermostat Reservoir\n(Cavity Heat Exchange)', fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    # 5. Universe Total Energy
    ax5.plot(time, derived_energies['universe_total_energy'], 'purple', linewidth=2, label='Universe Total')
    ax5.set_xlabel('Time (ps)')
    ax5.set_ylabel('Energy (a.u.)')
    ax5.set_title('Universe Total Energy\n(System + All Thermostats)', fontweight='bold')
    ax5.grid(True, alpha=0.3)
    ax5.legend()
    
    # 6. Leave empty or add comparison plot
    ax6.axis('off')  # Turn off the 6th subplot
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_main_components.png', dpi=300, bbox_inches='tight')
    print(f"Saved main component plot: {output_prefix}_main_components.png")
    
    # Create a new detailed plot showing individual potential and kinetic energies
    fig2, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig2.suptitle(f'Energy Analysis: Potential and Kinetic Components (T = {target_temp:.1f} K)', 
                  fontsize=16, fontweight='bold')
    
    # 1. Molecular Energy Components
    ax1.plot(time, derived_energies['molecular_potential_energy'], 'r-', linewidth=2, 
             label='Molecular Potential', alpha=0.8)
    ax1.plot(time, derived_energies['molecular_kinetic_energy'], 'b-', linewidth=2, 
             label='Molecular Kinetic', alpha=0.8)
    ax1.plot(time, derived_energies['molecular_total_energy'], 'k-', linewidth=2, 
             label='Molecular Total', alpha=0.9)
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Energy (a.u.)')
    ax1.set_title('Molecular System Energy Components\n(Harmonic + LJ + Coulomb + Kinetic)', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # 2. Cavity Energy Components
    ax2.plot(time, derived_energies['cavity_potential_energy'], 'r-', linewidth=2, 
             label='Cavity Potential', alpha=0.8)
    ax2.plot(time, derived_energies['cavity_kinetic_energy'], 'b-', linewidth=2, 
             label='Cavity Kinetic', alpha=0.8)
    ax2.plot(time, derived_energies['cavity_total_energy'], 'k-', linewidth=2, 
             label='Cavity Total', alpha=0.9)
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Energy (a.u.)')
    ax2.set_title('Cavity System Energy Components\n(Harmonic + Coupling + Self-energy + Kinetic)', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # 3. Potential Energy Comparison
    ax3.plot(time, derived_energies['molecular_potential_energy'], 'r-', linewidth=2, 
             label='Molecular Potential', alpha=0.8)
    ax3.plot(time, derived_energies['cavity_potential_energy'], 'b-', linewidth=2, 
             label='Cavity Potential', alpha=0.8)
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('Energy (a.u.)')
    ax3.set_title('Potential Energy Comparison\n(Molecular vs Cavity)', fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # 4. Kinetic Energy Comparison
    ax4.plot(time, derived_energies['molecular_kinetic_energy'], 'r-', linewidth=2, 
             label='Molecular Kinetic', alpha=0.8)
    ax4.plot(time, derived_energies['cavity_kinetic_energy'], 'b-', linewidth=2, 
             label='Cavity Kinetic', alpha=0.8)
    ax4.set_xlabel('Time (ps)')
    ax4.set_ylabel('Energy (a.u.)')
    ax4.set_title('Kinetic Energy Comparison\n(Molecular vs Cavity)', fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_potential_kinetic_components.png', dpi=300, bbox_inches='tight')
    print(f"Saved potential/kinetic component plot: {output_prefix}_potential_kinetic_components.png")
    
    # Create a summary plot with all components on the same axes (updated)
    fig3, ax = plt.subplots(1, 1, figsize=(14, 8))
    ax.plot(time, derived_energies['cavity_total_energy'], 'b-', linewidth=2, label='Cavity Total Energy')
    ax.plot(time, derived_energies['molecular_total_energy'], 'r-', linewidth=2, label='Molecular Total Energy')
    ax.plot(time, derived_energies['molecular_reservoir_energy'], 'g-', linewidth=2, label='Molecular Reservoir Energy')
    ax.plot(time, derived_energies['cavity_reservoir_energy'], 'orange', linewidth=2, label='Cavity Reservoir Energy')
    ax.plot(time, derived_energies['universe_total_energy'], 'purple', linewidth=2, label='Universe Total Energy')
    
    ax.set_xlabel('Time (ps)', fontsize=12)
    ax.set_ylabel('Energy (a.u.)', fontsize=12)
    ax.set_title(f'Energy Components Overview (T = {target_temp:.1f} K)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_overview.png', dpi=300, bbox_inches='tight')
    print(f"Saved overview plot: {output_prefix}_overview.png")
    
    return fig, fig2, fig3

def print_energy_statistics(derived_energies, max_time=None, target_temp=100.0):
    """
    Print statistical analysis of energy components including individual P.E. and K.E.
    
    Args:
        derived_energies (dict): Dictionary containing all energy components
        max_time (float): Maximum time to analyze (None for all data)
        target_temp (float): Target temperature in Kelvin for thermal energy calculations
    """
    print("\n" + "="*60)
    print("ENERGY STATISTICS ANALYSIS")
    print("="*60)
    
    # Apply time filter if specified
    if max_time is not None:
        print(f"Analysis limited to first {max_time} ps")
        # Note: time filtering already applied in main function
    
    # Calculate thermal energy reference
    kB = 3.167e-6  # Hartree/K
    kT = kB * target_temp
    print(f"Target temperature: {target_temp:.1f} K")
    print(f"Thermal energy (kT): {kT:.6f} a.u.")
    
    # Energy conservation analysis
    universe_initial = derived_energies['universe_total_energy'].iloc[0]
    universe_final = derived_energies['universe_total_energy'].iloc[-1]
    universe_drift = universe_final - universe_initial
    universe_std = np.std(derived_energies['universe_total_energy'])
    
    print(f"\nENERGY CONSERVATION:")
    print(f"  Initial universe energy: {universe_initial:.6f} a.u.")
    print(f"  Final universe energy:   {universe_final:.6f} a.u.")
    print(f"  Total energy drift:      {universe_drift:.6f} a.u. ({universe_drift/kT:.2f} kT)")
    print(f"  Relative drift:          {abs(universe_drift/universe_initial)*100:.4f}%")
    print(f"  Universe energy std:     {universe_std:.6f} a.u. ({universe_std/kT:.2f} kT)")
    
    # Component analysis - Individual molecular and cavity components
    print(f"\nINDIVIDUAL ENERGY COMPONENT AVERAGES:")
    component_names = {
        'molecular_potential_energy': 'Molecular Potential',
        'molecular_kinetic_energy': 'Molecular Kinetic',
        'molecular_total_energy': 'Molecular Total',
        'cavity_potential_energy': 'Cavity Potential',
        'cavity_kinetic_energy': 'Cavity Kinetic',
        'cavity_total_energy': 'Cavity Total',
        'molecular_reservoir_energy': 'Molecular Reservoir',
        'cavity_reservoir_energy': 'Cavity Reservoir',
        'universe_total_energy': 'Universe Total'
    }
    
    for component, name in component_names.items():
        if component in derived_energies:
            avg_val = np.mean(derived_energies[component])
            std_val = np.std(derived_energies[component])
            min_val = np.min(derived_energies[component])
            max_val = np.max(derived_energies[component])
            print(f"  {name:20s}: {avg_val:10.6f} ± {std_val:8.6f} a.u. ({avg_val/kT:6.1f} ± {std_val/kT:5.1f} kT)")
            print(f"  {'':<20s}  Range: {min_val:.6f} to {max_val:.6f} a.u. ({min_val/kT:.1f} to {max_val/kT:.1f} kT)")
    
    # Energy distribution analysis
    print(f"\nENERGY DISTRIBUTION ANALYSIS:")
    
    # Molecular system energy distribution
    mol_potential_avg = np.mean(derived_energies['molecular_potential_energy'])
    mol_kinetic_avg = np.mean(derived_energies['molecular_kinetic_energy'])
    mol_total_avg = np.mean(derived_energies['molecular_total_energy'])
    
    print(f"  Molecular System:")
    print(f"    Potential energy:     {mol_potential_avg:.6f} a.u. ({mol_potential_avg/kT:.1f} kT)")
    print(f"    Kinetic energy:       {mol_kinetic_avg:.6f} a.u. ({mol_kinetic_avg/kT:.1f} kT)")
    print(f"    Total energy:         {mol_total_avg:.6f} a.u. ({mol_total_avg/kT:.1f} kT)")
    print(f"    P.E./K.E. ratio:      {abs(mol_potential_avg/mol_kinetic_avg):.2f}")
    
    # Cavity system energy distribution
    cav_potential_avg = np.mean(derived_energies['cavity_potential_energy'])
    cav_kinetic_avg = np.mean(derived_energies['cavity_kinetic_energy'])
    cav_total_avg = np.mean(derived_energies['cavity_total_energy'])
    
    print(f"  Cavity System:")
    print(f"    Potential energy:     {cav_potential_avg:.6f} a.u. ({cav_potential_avg/kT:.1f} kT)")
    print(f"    Kinetic energy:       {cav_kinetic_avg:.6f} a.u. ({cav_kinetic_avg/kT:.1f} kT)")
    print(f"    Total energy:         {cav_total_avg:.6f} a.u. ({cav_total_avg/kT:.1f} kT)")
    if cav_kinetic_avg != 0:
        print(f"    P.E./K.E. ratio:      {abs(cav_potential_avg/cav_kinetic_avg):.2f}")
    else:
        print(f"    P.E./K.E. ratio:      N/A (zero kinetic energy)")
    
    # System comparison
    print(f"  System Comparison:")
    print(f"    Mol/Cav P.E. ratio:   {abs(mol_potential_avg/cav_potential_avg):.2f}")
    if cav_kinetic_avg != 0:
        print(f"    Mol/Cav K.E. ratio:   {abs(mol_kinetic_avg/cav_kinetic_avg):.2f}")
    print(f"    Mol/Cav total ratio:  {abs(mol_total_avg/cav_total_avg):.2f}")
    
    # Thermostat energy analysis (both reservoir components)
    mol_res_initial = derived_energies['molecular_reservoir_energy'].iloc[0]
    mol_res_final = derived_energies['molecular_reservoir_energy'].iloc[-1]
    mol_res_change = mol_res_final - mol_res_initial
    
    cav_res_initial = derived_energies['cavity_reservoir_energy'].iloc[0]
    cav_res_final = derived_energies['cavity_reservoir_energy'].iloc[-1]
    cav_res_change = cav_res_final - cav_res_initial
    
    print(f"\nTHERMOSTAT ENERGY ANALYSIS:")
    print(f"  Molecular Reservoir:")
    print(f"    Initial energy: {mol_res_initial:.6f} a.u. ({mol_res_initial/kT:.1f} kT)")
    print(f"    Final energy:   {mol_res_final:.6f} a.u. ({mol_res_final/kT:.1f} kT)")
    print(f"    Net change:     {mol_res_change:.6f} a.u. ({mol_res_change/kT:.1f} kT)")
    
    print(f"  Cavity Reservoir:")
    print(f"    Initial energy: {cav_res_initial:.6f} a.u. ({cav_res_initial/kT:.1f} kT)")
    print(f"    Final energy:   {cav_res_final:.6f} a.u. ({cav_res_final/kT:.1f} kT)")
    print(f"    Net change:     {cav_res_change:.6f} a.u. ({cav_res_change/kT:.1f} kT)")
    
    total_thermostat_change = mol_res_change + cav_res_change
    print(f"  Total thermostat change: {total_thermostat_change:.6f} a.u. ({total_thermostat_change/kT:.1f} kT)")
    if total_thermostat_change > 0:
        print(f"  → Net energy ADDED to thermostats (system cooled)")
    else:
        print(f"  → Net energy REMOVED from thermostats (system heated)")
    
    print("="*60)

def read_temperature_data(energy_file):
    """
    Read temperature data from kinetic energy and cavity mode files.
    
    Args:
        energy_file (str): Path to the energy contributions file (used to determine prefix)
        
    Returns:
        dict: Dictionary containing molecular and cavity temperature data
    """
    # Determine file prefix from energy file name
    # e.g., "prod-0_energy_contributions.txt" -> "prod-0"
    base_path = Path(energy_file)
    if '_energy_contributions.txt' in energy_file:
        prefix = energy_file.replace('_energy_contributions.txt', '')
    elif '_energy.txt' in energy_file:
        prefix = energy_file.replace('_energy.txt', '')
    else:
        # Try to extract prefix by removing common suffixes
        prefix = str(base_path.stem)
        for suffix in ['_energy_contributions', '_energy']:
            if prefix.endswith(suffix):
                prefix = prefix[:-len(suffix)]
                break
    
    temperature_data = {}
    
    # Read molecular temperature from kinetic energy file
    molecular_temp_file = f"{prefix}_kinetic_energy.txt"
    if os.path.exists(molecular_temp_file):
        try:
            # Read molecular kinetic energy file
            # Columns: time(ps) timestep kinetic_energy(a.u.) temperature(K) avg_vel_magnitude(a.u.)
            mol_data = pd.read_csv(molecular_temp_file, sep=r'\s+', comment='#', header=None)
            mol_data.columns = ['time_ps', 'timestep', 'kinetic_energy_au', 'temperature_K', 'avg_vel_magnitude_au']
            
            temperature_data['molecular'] = {
                'time': mol_data['time_ps'].values,
                'temperature': mol_data['temperature_K'].values,
                'kinetic_energy': mol_data['kinetic_energy_au'].values,
                'file': molecular_temp_file
            }
            print(f"Successfully loaded molecular temperature data from {molecular_temp_file}")
            print(f"  Data points: {len(mol_data)}")
            print(f"  Time range: {mol_data['time_ps'].iloc[0]:.3f} - {mol_data['time_ps'].iloc[-1]:.3f} ps")
            print(f"  Temperature range: {mol_data['temperature_K'].min():.1f} - {mol_data['temperature_K'].max():.1f} K")
        except Exception as e:
            print(f"Warning: Could not read molecular temperature file {molecular_temp_file}: {e}")
    else:
        print(f"Warning: Molecular temperature file {molecular_temp_file} not found")
    
    # Read cavity temperature from cavity mode file
    cavity_temp_file = f"{prefix}_cavity_mode.txt"
    if os.path.exists(cavity_temp_file):
        try:
            # Read cavity mode file - need to handle the header formatting issue
            # The header seems to have a line break in the middle, so we'll skip comment lines and read data
            cav_data = pd.read_csv(cavity_temp_file, sep=r'\s+', comment='#', header=None)
            
            # Check if we have the expected number of columns (should be 12)
            if cav_data.shape[1] >= 6:
                # Assign column names based on the expected structure
                cav_data.columns = ['time_ps', 'timestep', 'cavity_ke_au', 'cavity_pe_harmonic_au', 
                                   'cavity_total_energy_au', 'cavity_temperature_K'] + [f'col_{i}' for i in range(6, cav_data.shape[1])]
                
                temperature_data['cavity'] = {
                    'time': cav_data['time_ps'].values,
                    'temperature': cav_data['cavity_temperature_K'].values,
                    'kinetic_energy': cav_data['cavity_ke_au'].values,
                    'file': cavity_temp_file
                }
                print(f"Successfully loaded cavity temperature data from {cavity_temp_file}")
                print(f"  Data points: {len(cav_data)}")
                print(f"  Time range: {cav_data['time_ps'].iloc[0]:.3f} - {cav_data['time_ps'].iloc[-1]:.3f} ps")
                print(f"  Temperature range: {cav_data['cavity_temperature_K'].min():.1f} - {cav_data['cavity_temperature_K'].max():.1f} K")
            else:
                print(f"Warning: Cavity file {cavity_temp_file} has unexpected format (only {cav_data.shape[1]} columns)")
        except Exception as e:
            print(f"Warning: Could not read cavity temperature file {cavity_temp_file}: {e}")
    else:
        print(f"Note: Cavity temperature file {cavity_temp_file} not found (normal for non-cavity simulations)")
    
    return temperature_data

def create_temperature_plots(temperature_data, output_prefix='energy_analysis', max_time=None, target_temp=100.0):
    """
    Create temperature time series plots.
    
    Args:
        temperature_data (dict): Dictionary containing temperature data
        output_prefix (str): Prefix for output files
        max_time (float): Maximum time to plot (None for all data)
        target_temp (float): Target temperature in Kelvin for reference line
    """
    if not temperature_data:
        print("No temperature data available for plotting")
        return None
    
    # Determine number of subplots based on available data
    has_molecular = 'molecular' in temperature_data
    has_cavity = 'cavity' in temperature_data
    
    if has_molecular and has_cavity:
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 12))
        fig.suptitle('Temperature Time Series Analysis', fontsize=16, fontweight='bold')
    elif has_molecular or has_cavity:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
        fig.suptitle('Temperature Time Series Analysis', fontsize=16, fontweight='bold')
    else:
        print("No valid temperature data found")
        return None
    
    plot_idx = 0
    
    # Plot molecular temperature
    if has_molecular:
        mol_time = temperature_data['molecular']['time']
        mol_temp = temperature_data['molecular']['temperature']
        
        # Apply time filter if specified
        if max_time is not None:
            mask = mol_time <= max_time
            mol_time = mol_time[mask]
            mol_temp = mol_temp[mask]
        
        ax = [ax1, ax2, ax3][plot_idx] if has_molecular and has_cavity else [ax1, ax2][plot_idx]
        ax.plot(mol_time, mol_temp, 'r-', linewidth=1.5, label='Molecular Temperature')
        ax.axhline(y=target_temp, color='k', linestyle='--', alpha=0.7, label=f'Target ({target_temp} K)')
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('Temperature (K)')
        ax.set_title('Molecular System Temperature\n(O and N particles)', fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add statistics text using last half of data for consistency with analysis
        n_points = len(mol_temp)
        last_half_start = n_points // 2
        mol_temp_last_half = mol_temp[last_half_start:]
        mol_time_last_half = mol_time[last_half_start:]
        
        mol_mean = np.mean(mol_temp_last_half)
        mol_std = np.std(mol_temp_last_half)
        
        # Add vertical line to show where analysis period starts
        analysis_start_time = mol_time_last_half[0]
        ax.axvline(x=analysis_start_time, color='gray', linestyle=':', alpha=0.7, 
                  label=f'Analysis start ({analysis_start_time:.1f} ps)')
        
        ax.text(0.02, 0.98, f'Mean (last half): {mol_mean:.1f} ± {mol_std:.1f} K\nAnalysis: {analysis_start_time:.1f}-{mol_time[-1]:.1f} ps', 
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plot_idx += 1
    
    # Plot cavity temperature
    if has_cavity:
        cav_time = temperature_data['cavity']['time']
        cav_temp = temperature_data['cavity']['temperature']
        
        # Apply time filter if specified
        if max_time is not None:
            mask = cav_time <= max_time
            cav_time = cav_time[mask]
            cav_temp = cav_temp[mask]
        
        ax = [ax1, ax2, ax3][plot_idx] if has_molecular and has_cavity else [ax1, ax2][plot_idx]
        ax.plot(cav_time, cav_temp, 'b-', linewidth=1.5, label='Cavity Temperature')
        ax.axhline(y=target_temp, color='k', linestyle='--', alpha=0.7, label=f'Target ({target_temp} K)')
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('Temperature (K)')
        ax.set_title('Cavity Particle Temperature\n(Photon mode)', fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add statistics text using last half of data for consistency with analysis
        n_points = len(cav_temp)
        last_half_start = n_points // 2
        cav_temp_last_half = cav_temp[last_half_start:]
        cav_time_last_half = cav_time[last_half_start:]
        
        cav_mean = np.mean(cav_temp_last_half)
        cav_std = np.std(cav_temp_last_half)
        
        # Add vertical line to show where analysis period starts
        analysis_start_time = cav_time_last_half[0]
        ax.axvline(x=analysis_start_time, color='gray', linestyle=':', alpha=0.7, 
                  label=f'Analysis start ({analysis_start_time:.1f} ps)')
        
        ax.text(0.02, 0.98, f'Mean (last half): {cav_mean:.1f} ± {cav_std:.1f} K\nAnalysis: {analysis_start_time:.1f}-{cav_time[-1]:.1f} ps', 
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plot_idx += 1
    
    # Combined plot if both are available
    if has_molecular and has_cavity:
        mol_time = temperature_data['molecular']['time']
        mol_temp = temperature_data['molecular']['temperature']
        cav_time = temperature_data['cavity']['time']
        cav_temp = temperature_data['cavity']['temperature']
        
        # Apply time filter if specified
        if max_time is not None:
            mol_mask = mol_time <= max_time
            mol_time = mol_time[mol_mask]
            mol_temp = mol_temp[mol_mask]
            cav_mask = cav_time <= max_time
            cav_time = cav_time[cav_mask]
            cav_temp = cav_temp[cav_mask]
        
        ax3.plot(mol_time, mol_temp, 'r-', linewidth=1.5, label='Molecular Temperature', alpha=0.8)
        ax3.plot(cav_time, cav_temp, 'b-', linewidth=1.5, label='Cavity Temperature', alpha=0.8)
        ax3.axhline(y=target_temp, color='k', linestyle='--', alpha=0.7, label=f'Target ({target_temp} K)')
        ax3.set_xlabel('Time (ps)')
        ax3.set_ylabel('Temperature (K)')
        ax3.set_title('Temperature Comparison\n(Molecular vs Cavity)', fontweight='bold')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        
        # Add comparison statistics using last half of data
        mol_n_points = len(mol_temp)
        mol_last_half_start = mol_n_points // 2
        mol_temp_last_half = mol_temp[mol_last_half_start:]
        mol_time_last_half = mol_time[mol_last_half_start:]
        
        cav_n_points = len(cav_temp)
        cav_last_half_start = cav_n_points // 2
        cav_temp_last_half = cav_temp[cav_last_half_start:]
        cav_time_last_half = cav_time[cav_last_half_start:]
        
        mol_mean = np.mean(mol_temp_last_half)
        cav_mean = np.mean(cav_temp_last_half)
        
        # Add vertical lines to show where analysis periods start
        mol_analysis_start = mol_time_last_half[0]
        cav_analysis_start = cav_time_last_half[0]
        analysis_start = min(mol_analysis_start, cav_analysis_start)
        
        ax3.axvline(x=analysis_start, color='gray', linestyle=':', alpha=0.7, 
                   label=f'Analysis start ({analysis_start:.1f} ps)')
        
        temp_diff = mol_mean - cav_mean
        ax3.text(0.02, 0.98, f'Mean (last half):\nMolecular: {mol_mean:.1f} K\nCavity: {cav_mean:.1f} K\nDifference: {temp_diff:+.1f} K', 
                transform=ax3.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_temperature_timeseries.png', dpi=300, bbox_inches='tight')
    print(f"Saved temperature time series plot: {output_prefix}_temperature_timeseries.png")
    
    return fig

def print_temperature_statistics(temperature_data, max_time=None, target_temp=100.0):
    """
    Print statistical analysis of temperature data.
    Statistics are computed using only the last half of the data to exclude equilibration effects.
    
    Args:
        temperature_data (dict): Dictionary containing temperature data
        max_time (float): Maximum time to analyze (None for all data)
        target_temp (float): Target temperature in Kelvin
    """
    if not temperature_data:
        print("No temperature data available for analysis")
        return
    
    print("\n" + "="*60)
    print("TEMPERATURE STATISTICS ANALYSIS")
    print("="*60)
    print("(Statistics computed using last half of data to exclude equilibration)")
    
    if max_time is not None:
        print(f"Analysis limited to first {max_time} ps")
    
    print(f"Target temperature: {target_temp} K")
    
    # Analyze molecular temperature
    if 'molecular' in temperature_data:
        mol_temp = temperature_data['molecular']['temperature']
        mol_time = temperature_data['molecular']['time']
        
        # Apply time filter if specified
        if max_time is not None:
            mask = mol_time <= max_time
            mol_temp = mol_temp[mask]
            mol_time_filtered = mol_time[mask]
        else:
            mol_time_filtered = mol_time
        
        # Use only last half of data for statistics
        n_points = len(mol_temp)
        last_half_start = n_points // 2
        mol_temp_last_half = mol_temp[last_half_start:]
        
        mol_mean = np.mean(mol_temp_last_half)
        mol_std = np.std(mol_temp_last_half)
        mol_min = np.min(mol_temp_last_half)
        mol_max = np.max(mol_temp_last_half)
        mol_deviation = mol_mean - target_temp
        
        # Time range for last half analysis
        time_start_analysis = mol_time_filtered[last_half_start]
        time_end_analysis = mol_time_filtered[-1]
        
        print(f"\nMOLECULAR SYSTEM TEMPERATURE:")
        print(f"  Analysis period:      {time_start_analysis:.2f} - {time_end_analysis:.2f} ps (last half)")
        print(f"  Data points used:     {len(mol_temp_last_half)} / {n_points}")
        print(f"  Mean temperature:     {mol_mean:.2f} ± {mol_std:.2f} K")
        print(f"  Temperature range:    {mol_min:.2f} - {mol_max:.2f} K")
        print(f"  Deviation from target: {mol_deviation:+.2f} K ({mol_deviation/target_temp*100:+.1f}%)")
        print(f"  Coefficient of variation: {mol_std/mol_mean*100:.2f}%")
    
    # Analyze cavity temperature
    if 'cavity' in temperature_data:
        cav_temp = temperature_data['cavity']['temperature']
        cav_time = temperature_data['cavity']['time']
        
        # Apply time filter if specified
        if max_time is not None:
            mask = cav_time <= max_time
            cav_temp = cav_temp[mask]
            cav_time_filtered = cav_time[mask]
        else:
            cav_time_filtered = cav_time
        
        # Use only last half of data for statistics
        n_points = len(cav_temp)
        last_half_start = n_points // 2
        cav_temp_last_half = cav_temp[last_half_start:]
        
        cav_mean = np.mean(cav_temp_last_half)
        cav_std = np.std(cav_temp_last_half)
        cav_min = np.min(cav_temp_last_half)
        cav_max = np.max(cav_temp_last_half)
        cav_deviation = cav_mean - target_temp
        
        # Time range for last half analysis
        time_start_analysis = cav_time_filtered[last_half_start]
        time_end_analysis = cav_time_filtered[-1]
        
        print(f"\nCAVITY PARTICLE TEMPERATURE:")
        print(f"  Analysis period:      {time_start_analysis:.2f} - {time_end_analysis:.2f} ps (last half)")
        print(f"  Data points used:     {len(cav_temp_last_half)} / {n_points}")
        print(f"  Mean temperature:     {cav_mean:.2f} ± {cav_std:.2f} K")
        print(f"  Temperature range:    {cav_min:.2f} - {cav_max:.2f} K")
        print(f"  Deviation from target: {cav_deviation:+.2f} K ({cav_deviation/target_temp*100:+.1f}%)")
        print(f"  Coefficient of variation: {cav_std/cav_mean*100:.2f}%")
    
    # Compare temperatures if both are available
    if 'molecular' in temperature_data and 'cavity' in temperature_data:
        mol_temp = temperature_data['molecular']['temperature']
        mol_time = temperature_data['molecular']['time']
        cav_temp = temperature_data['cavity']['temperature']
        cav_time = temperature_data['cavity']['time']
        
        # Apply time filter if specified
        if max_time is not None:
            mol_mask = mol_time <= max_time
            mol_temp = mol_temp[mol_mask]
            cav_mask = cav_time <= max_time
            cav_temp = cav_temp[cav_mask]
        
        # Use only last half of data for comparison statistics
        mol_n_points = len(mol_temp)
        mol_last_half_start = mol_n_points // 2
        mol_temp_last_half = mol_temp[mol_last_half_start:]
        
        cav_n_points = len(cav_temp)
        cav_last_half_start = cav_n_points // 2
        cav_temp_last_half = cav_temp[cav_last_half_start:]
        
        mol_mean = np.mean(mol_temp_last_half)
        cav_mean = np.mean(cav_temp_last_half)
        temp_difference = mol_mean - cav_mean
        
        print(f"\nTEMPERATURE COMPARISON (last half of data):")
        print(f"  Molecular mean:       {mol_mean:.2f} K")
        print(f"  Cavity mean:          {cav_mean:.2f} K")
        print(f"  Difference (Mol-Cav): {temp_difference:+.2f} K")
        if abs(temp_difference) > 5.0:
            print(f"  → Significant temperature difference detected!")
        else:
            print(f"  → Temperatures are well equilibrated")
    
    print("="*60)

def maxwell_boltzmann_3d(v, mass, kT):
    """
    3D Maxwell-Boltzmann distribution for velocity magnitude.
    
    Args:
        v: velocity magnitude array
        mass: particle mass in atomic units
        kT: thermal energy (kB * T) in atomic units
        
    Returns:
        probability density
    """
    # Normalization constant
    norm = 4 * np.pi * (mass / (2 * np.pi * kT))**(3/2)
    # Exponential term
    exp_term = np.exp(-mass * v**2 / (2 * kT))
    return norm * v**2 * exp_term

def maxwell_boltzmann_1d(v, mass, kT):
    """
    1D Maxwell-Boltzmann distribution for velocity component.
    
    Args:
        v: velocity component array
        mass: particle mass in atomic units
        kT: thermal energy (kB * T) in atomic units
        
    Returns:
        probability density
    """
    # Normalization constant
    norm = np.sqrt(mass / (2 * np.pi * kT))
    # Exponential term
    exp_term = np.exp(-mass * v**2 / (2 * kT))
    return norm * exp_term

def read_gsd_velocities(gsd_filename, frame_indices=None, particle_types=None):
    """
    Read velocity data from GSD trajectory file.
    
    Args:
        gsd_filename (str): Path to GSD file
        frame_indices (list): List of frame indices to read (None for all frames)
        particle_types (list): List of particle types to include (None for all)
        
    Returns:
        dict: Dictionary containing velocity data organized by particle type
    """
    print(f"Reading GSD file: {gsd_filename}")
    
    velocity_data = {}
    
    try:
        with gsd.hoomd.open(gsd_filename, 'r') as f:
            n_frames = len(f)
            print(f"Total frames in GSD file: {n_frames}")
            
            if frame_indices is None:
                frame_indices = list(range(n_frames))
            else:
                # Ensure frame indices are valid
                frame_indices = [i for i in frame_indices if 0 <= i < n_frames]
            
            print(f"Reading {len(frame_indices)} frames")
            
            # Get particle type information from first frame
            first_frame = f[0]
            particle_typenames = first_frame.particles.types
            particle_typeids = first_frame.particles.typeid
            masses = first_frame.particles.mass
            
            print(f"Particle types: {particle_typenames}")
            print(f"Particle type IDs: {np.unique(particle_typeids)}")
            
            # Initialize velocity storage
            for i, typename in enumerate(particle_typenames):
                if particle_types is None or typename in particle_types:
                    velocity_data[typename] = {
                        'velocities': [],
                        'masses': [],
                        'times': [],
                        'typeid': i
                    }
            
            # Read velocities from selected frames
            for frame_idx in frame_indices:
                frame = f[frame_idx]
                
                # Extract time if available
                if hasattr(frame.configuration, 'step'):
                    # Estimate time from step (assuming dt ~ 0.001 ps)
                    time_ps = frame.configuration.step * 0.001
                else:
                    time_ps = frame_idx * 0.1  # Fallback estimate
                
                # Get velocities and convert from momentum
                # HOOMD stores momentum = mass * velocity
                momenta = frame.particles.velocity  # This is actually momentum in HOOMD
                particle_masses = frame.particles.mass
                
                # Convert momentum to velocity
                velocities = momenta / particle_masses[:, np.newaxis]
                
                # Group by particle type
                for typename in velocity_data.keys():
                    typeid = velocity_data[typename]['typeid']
                    mask = particle_typeids == typeid
                    
                    if np.any(mask):
                        type_velocities = velocities[mask]
                        type_masses = particle_masses[mask]
                        
                        velocity_data[typename]['velocities'].append(type_velocities)
                        velocity_data[typename]['masses'].append(type_masses)
                        velocity_data[typename]['times'].append(time_ps)
            
            # Convert lists to arrays
            for typename in velocity_data.keys():
                if velocity_data[typename]['velocities']:
                    velocity_data[typename]['velocities'] = np.array(velocity_data[typename]['velocities'])
                    velocity_data[typename]['masses'] = np.array(velocity_data[typename]['masses'])
                    velocity_data[typename]['times'] = np.array(velocity_data[typename]['times'])
                    
                    print(f"Type {typename}: {velocity_data[typename]['velocities'].shape} velocity samples")
        
        return velocity_data
    
    except Exception as e:
        print(f"Error reading GSD file {gsd_filename}: {e}")
        return {}

def calculate_velocity_statistics(velocities, masses, temperature_K=100.0):
    """
    Calculate velocity statistics and compare with Maxwell-Boltzmann theory.
    
    Args:
        velocities (np.array): Velocity array [frames, particles, 3]
        masses (np.array): Mass array [frames, particles]
        temperature_K (float): Target temperature in Kelvin
        
    Returns:
        dict: Dictionary containing velocity statistics
    """
    # Convert temperature to atomic units
    kB = 3.167e-6  # Hartree/K
    kT = kB * temperature_K
    
    # Flatten arrays across frames and particles
    vel_flat = velocities.reshape(-1, 3)  # [N_samples, 3]
    mass_flat = masses.reshape(-1)  # [N_samples]
    
    # Calculate velocity magnitudes
    vel_magnitudes = np.sqrt(np.sum(vel_flat**2, axis=1))
    
    # Calculate individual components
    vel_x = vel_flat[:, 0]
    vel_y = vel_flat[:, 1]
    vel_z = vel_flat[:, 2]
    
    # Calculate kinetic energy and temperature
    kinetic_energies = 0.5 * mass_flat * vel_magnitudes**2
    avg_kinetic_energy = np.mean(kinetic_energies)
    
    # Temperature from equipartition theorem: <KE> = (3/2) * kT
    measured_temp = (2.0/3.0) * avg_kinetic_energy / kB
    
    # Theoretical expectations for Maxwell-Boltzmann
    avg_mass = np.mean(mass_flat)
    theoretical_v_mean = np.sqrt(8 * kT / (np.pi * avg_mass))  # <v> for 3D MB
    theoretical_v_rms = np.sqrt(3 * kT / avg_mass)  # RMS velocity
    theoretical_v_most_probable = np.sqrt(2 * kT / avg_mass)  # Most probable velocity
    
    stats_dict = {
        'n_samples': len(vel_magnitudes),
        'avg_mass': avg_mass,
        'measured_temperature': measured_temp,
        'target_temperature': temperature_K,
        'temp_error_percent': 100 * (measured_temp - temperature_K) / temperature_K,
        'velocity_magnitudes': vel_magnitudes,
        'velocity_components': {'x': vel_x, 'y': vel_y, 'z': vel_z},
        'kinetic_energies': kinetic_energies,
        'avg_kinetic_energy': avg_kinetic_energy,
        'theoretical_v_mean': theoretical_v_mean,
        'theoretical_v_rms': theoretical_v_rms,
        'theoretical_v_most_probable': theoretical_v_most_probable,
        'measured_v_mean': np.mean(vel_magnitudes),
        'measured_v_rms': np.sqrt(np.mean(vel_magnitudes**2)),
        'measured_v_std': np.std(vel_magnitudes),
        'kT': kT,
        'masses': mass_flat
    }
    
    return stats_dict

def create_velocity_histograms(gsd_file, temperature_K=100.0, output_prefix='energy_analysis', 
                             cavity_only=False, molecular_only=False, max_frames=50):
    """
    Create velocity histogram plots from GSD trajectory file.
    
    Args:
        gsd_file (str): Path to GSD trajectory file
        temperature_K (float): Target temperature in Kelvin
        output_prefix (str): Prefix for output files
        cavity_only (bool): Plot only cavity particle velocities
        molecular_only (bool): Plot only molecular particle velocities
        max_frames (int): Maximum number of frames to analyze
    """
    # Check if GSD file exists
    if not os.path.exists(gsd_file):
        print(f"Warning: GSD file {gsd_file} not found. Skipping velocity analysis.")
        return
    
    print(f"\n{'='*60}")
    print("VELOCITY HISTOGRAM ANALYSIS")
    print(f"{'='*60}")
    
    # Read velocity data from last frames for better statistics
    with gsd.hoomd.open(gsd_file, 'r') as f:
        total_frames = len(f)
    
    # Use last frames for analysis (more equilibrated)
    start_frame = max(0, total_frames - max_frames)
    frame_indices = list(range(start_frame, total_frames))
    
    print(f"Analyzing {len(frame_indices)} frames from {start_frame} to {total_frames-1}")
    
    velocity_data = read_gsd_velocities(gsd_file, frame_indices=frame_indices)
    
    if not velocity_data:
        print("No velocity data available for analysis")
        return
    
    # Filter particle types based on options
    plot_types = []
    for typename in velocity_data.keys():
        if cavity_only and typename not in ['L']:
            continue
        if molecular_only and typename not in ['O', 'N']:
            continue
        if velocity_data[typename]['velocities'].size > 0:
            plot_types.append(typename)
    
    if not plot_types:
        print("No velocity data available for the specified particle types")
        return
    
    # Create figure with subplots
    n_types = len(plot_types)
    fig, axes = plt.subplots(n_types, 3, figsize=(18, 6*n_types))
    if n_types == 1:
        axes = axes.reshape(1, -1)
    
    fig.suptitle(f'Velocity Distribution Analysis (T = {temperature_K:.1f} K)', 
                 fontsize=16, fontweight='bold')
    
    colors = ['red', 'blue', 'green', 'orange', 'purple']
    
    for i, typename in enumerate(plot_types):
        velocities = velocity_data[typename]['velocities']
        masses = velocity_data[typename]['masses']
        
        if velocities.size == 0:
            continue
            
        # Calculate statistics
        stats = calculate_velocity_statistics(velocities, masses, temperature_K)
        
        # Plot 1: Velocity magnitude histogram with Maxwell-Boltzmann fit
        ax1 = axes[i, 0]
        vel_magnitudes = stats['velocity_magnitudes']
        
        # Create histogram
        counts, bins, _ = ax1.hist(vel_magnitudes, bins=50, density=True, alpha=0.7, 
                                  color=colors[i % len(colors)], 
                                  label=f'{typename} particles')
        
        # Plot theoretical Maxwell-Boltzmann distribution
        v_theory = np.linspace(0, np.max(vel_magnitudes), 200)
        mb_3d = maxwell_boltzmann_3d(v_theory, stats['avg_mass'], stats['kT'])
        ax1.plot(v_theory, mb_3d, 'k-', linewidth=2, 
                label=f'Maxwell-Boltzmann (T={temperature_K:.1f} K)')
        
        # Add vertical lines for characteristic velocities
        ax1.axvline(stats['theoretical_v_most_probable'], color='red', linestyle='--', 
                   alpha=0.8, label=f'Most probable: {stats["theoretical_v_most_probable"]:.3f}')
        ax1.axvline(stats['measured_v_mean'], color='green', linestyle='--', 
                   alpha=0.8, label=f'Measured mean: {stats["measured_v_mean"]:.3f}')
        
        ax1.set_xlabel('Velocity Magnitude (a.u.)')
        ax1.set_ylabel('Probability Density')
        ax1.set_title(f'{typename} Particles: Velocity Magnitude Distribution')
        ax1.legend(fontsize=8)
        ax1.grid(True, alpha=0.3)
        
        # Add statistics text
        stats_text = (f'Measured T: {stats["measured_temperature"]:.1f} K\n'
                     f'Error: {stats["temp_error_percent"]:.1f}%\n'
                     f'Samples: {stats["n_samples"]:,}')
        ax1.text(0.98, 0.98, stats_text, transform=ax1.transAxes, 
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Plot 2: Velocity component histograms
        ax2 = axes[i, 1]
        vel_components = stats['velocity_components']
        
        for j, (comp, vel_comp) in enumerate(vel_components.items()):
            ax2.hist(vel_comp, bins=40, density=True, alpha=0.5, 
                    label=f'{comp}-component', color=['red', 'green', 'blue'][j])
        
        # Plot theoretical 1D Maxwell-Boltzmann for comparison
        v_comp_range = np.linspace(-4*np.sqrt(stats['kT']/stats['avg_mass']), 
                                  4*np.sqrt(stats['kT']/stats['avg_mass']), 200)
        mb_1d = maxwell_boltzmann_1d(v_comp_range, stats['avg_mass'], stats['kT'])
        ax2.plot(v_comp_range, mb_1d, 'k-', linewidth=2, 
                label=f'1D Maxwell-Boltzmann')
        
        ax2.set_xlabel('Velocity Component (a.u.)')
        ax2.set_ylabel('Probability Density')
        ax2.set_title(f'{typename} Particles: Velocity Components')
        ax2.legend(fontsize=8)
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Kinetic energy distribution
        ax3 = axes[i, 2]
        kinetic_energies = stats['kinetic_energies']
        
        ax3.hist(kinetic_energies, bins=50, density=True, alpha=0.7, 
                color=colors[i % len(colors)], label=f'{typename} particles')
        
        # Theoretical kinetic energy distribution (3D)
        # KE distribution: P(E) ∝ sqrt(E) * exp(-E/kT)
        ke_max = np.max(kinetic_energies) if len(kinetic_energies) > 0 else stats['kT'] * 10
        ke_theory = np.linspace(0, ke_max, 200)
        ke_dist = np.sqrt(ke_theory) * np.exp(-ke_theory / stats['kT'])
        ke_dist /= np.trapz(ke_dist, ke_theory)  # Normalize
        
        ax3.plot(ke_theory, ke_dist, 'k-', linewidth=2, 
                label=f'Theoretical (3D)')
        
        # Add mean kinetic energy line
        ax3.axvline(stats['avg_kinetic_energy'], color='red', linestyle='--', 
                   alpha=0.8, label=f'Mean KE: {stats["avg_kinetic_energy"]:.4f}')
        ax3.axvline(1.5 * stats['kT'], color='green', linestyle='--', 
                   alpha=0.8, label=f'3/2 kT: {1.5 * stats["kT"]:.4f}')
        
        ax3.set_xlabel('Kinetic Energy (a.u.)')
        ax3.set_ylabel('Probability Density')
        ax3.set_title(f'{typename} Particles: Kinetic Energy Distribution')
        ax3.legend(fontsize=8)
        ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_velocity_histograms.png', dpi=300, bbox_inches='tight')
    print(f"Saved velocity histogram plot: {output_prefix}_velocity_histograms.png")
    
    # Print summary statistics
    print("\n" + "="*60)
    print("VELOCITY DISTRIBUTION ANALYSIS SUMMARY")
    print("="*60)
    
    for typename in plot_types:
        velocities = velocity_data[typename]['velocities']
        masses = velocity_data[typename]['masses']
        
        if velocities.size == 0:
            continue
            
        stats = calculate_velocity_statistics(velocities, masses, temperature_K)
        
        print(f"\n{typename.upper()} PARTICLES:")
        print(f"  Number of samples: {stats['n_samples']:,}")
        print(f"  Average mass: {stats['avg_mass']:.6f} a.u.")
        print(f"  Target temperature: {stats['target_temperature']:.1f} K")
        print(f"  Measured temperature: {stats['measured_temperature']:.1f} K")
        print(f"  Temperature error: {stats['temp_error_percent']:.1f}%")
        print(f"  Average kinetic energy: {stats['avg_kinetic_energy']:.6f} a.u.")
        print(f"  Expected kinetic energy (3/2 kT): {1.5 * stats['kT']:.6f} a.u.")
        print(f"  Velocity statistics:")
        print(f"    Measured mean: {stats['measured_v_mean']:.4f} a.u.")
        print(f"    Theoretical mean: {stats['theoretical_v_mean']:.4f} a.u.")
        print(f"    Measured RMS: {stats['measured_v_rms']:.4f} a.u.")
        print(f"    Theoretical RMS: {stats['theoretical_v_rms']:.4f} a.u.")
        print(f"    Most probable (theory): {stats['theoretical_v_most_probable']:.4f} a.u.")

def read_cavity_velocity_data(energy_file, max_frames=None):
    """
    Read cavity velocity data from cavity mode files.
    
    Args:
        energy_file (str): Path to the energy contributions file (used to determine prefix)
        max_frames (int): Maximum number of frames to read (None for all data)
        
    Returns:
        dict: Dictionary containing cavity velocity data and metadata
    """
    # Determine file prefix from energy file name
    base_path = Path(energy_file)
    if '_energy_contributions.txt' in energy_file:
        prefix = energy_file.replace('_energy_contributions.txt', '')
    elif '_energy.txt' in energy_file:
        prefix = energy_file.replace('_energy.txt', '')
    else:
        # Try to extract prefix by removing common suffixes
        prefix = str(base_path.stem)
        for suffix in ['_energy_contributions', '_energy']:
            if prefix.endswith(suffix):
                prefix = prefix[:-len(suffix)]
                break
    
    cavity_mode_file = f"{prefix}_cavity_mode.txt"
    
    if not os.path.exists(cavity_mode_file):
        print(f"Warning: Cavity mode file {cavity_mode_file} not found")
        return None
    
    try:
        # Read cavity mode file
        # Columns: time(ps) timestep cavity_ke(a.u.) cavity_pe_harmonic(a.u.) cavity_total_energy(a.u.) 
        #          cavity_temperature(K) cavity_pos_x cavity_pos_y cavity_pos_z cavity_vel_x cavity_vel_y cavity_vel_z
        print(f"Reading cavity velocity data from {cavity_mode_file}...")
        
        # Read the header to understand the format
        with open(cavity_mode_file, 'r') as f:
            lines = f.readlines()
        
        # Find the actual data start (skip multi-line header)
        data_start_line = 0
        for i, line in enumerate(lines):
            if not line.startswith('#') and line.strip():
                # Check if this line contains data (starts with a number)
                try:
                    float(line.split()[0])
                    data_start_line = i
                    break
                except (ValueError, IndexError):
                    continue
        
        print(f"  Data starts at line {data_start_line + 1}")
        
        # Read the data, skipping header lines
        data = pd.read_csv(cavity_mode_file, sep=r'\s+', skiprows=data_start_line, header=None)
        
        # Check if we have the expected number of columns (should be at least 12)
        if data.shape[1] < 12:
            print(f"Warning: Cavity mode file has unexpected format (only {data.shape[1]} columns, expected at least 12)")
            return None
        
        # Assign column names
        data.columns = ['time_ps', 'timestep', 'cavity_ke_au', 'cavity_pe_harmonic_au', 
                       'cavity_total_energy_au', 'cavity_temperature_K', 
                       'cavity_pos_x', 'cavity_pos_y', 'cavity_pos_z',
                       'cavity_vel_x', 'cavity_vel_y', 'cavity_vel_z'] + [f'col_{i}' for i in range(12, data.shape[1])]
        
        # Limit frames if requested
        if max_frames is not None and len(data) > max_frames:
            # Take evenly spaced frames
            indices = np.linspace(0, len(data)-1, max_frames, dtype=int)
            data = data.iloc[indices].reset_index(drop=True)
            print(f"  Subsampled to {max_frames} frames from {len(data)} total")
        
        # Extract velocity components
        velocities = np.column_stack([
            data['cavity_vel_x'].values,
            data['cavity_vel_y'].values, 
            data['cavity_vel_z'].values
        ])
        
        # Calculate velocity magnitudes
        vel_magnitudes = np.sqrt(np.sum(velocities**2, axis=1))
        
        # Get cavity particle mass (assuming it's 1.0 a.u. for photon mode)
        cavity_mass = 1.0  # a.u.
        
        # Calculate kinetic energies
        kinetic_energies = 0.5 * cavity_mass * vel_magnitudes**2
        
        velocity_data = {
            'time': data['time_ps'].values,
            'velocities': velocities,
            'vel_magnitudes': vel_magnitudes,
            'vel_x': data['cavity_vel_x'].values,
            'vel_y': data['cavity_vel_y'].values,
            'vel_z': data['cavity_vel_z'].values,
            'kinetic_energies': kinetic_energies,
            'temperatures': data['cavity_temperature_K'].values,
            'mass': cavity_mass,
            'n_frames': len(data),
            'file': cavity_mode_file
        }
        
        print(f"Successfully loaded cavity velocity data:")
        print(f"  Data points: {len(data)}")
        print(f"  Time range: {data['time_ps'].iloc[0]:.3f} - {data['time_ps'].iloc[-1]:.3f} ps")
        print(f"  Velocity magnitude range: {vel_magnitudes.min():.4f} - {vel_magnitudes.max():.4f} a.u.")
        print(f"  Temperature range: {data['cavity_temperature_K'].min():.1f} - {data['cavity_temperature_K'].max():.1f} K")
        
        return velocity_data
        
    except Exception as e:
        print(f"Error reading cavity velocity data from {cavity_mode_file}: {e}")
        return None

def create_cavity_velocity_histogram(velocity_data, target_temp=100.0, output_prefix='energy_analysis'):
    """
    Create cavity velocity histogram with Maxwell-Boltzmann comparison.
    
    Args:
        velocity_data (dict): Cavity velocity data from read_cavity_velocity_data
        target_temp (float): Target temperature in Kelvin
        output_prefix (str): Output prefix for plot files
    """
    if velocity_data is None:
        print("No cavity velocity data available for histogram analysis")
        return
    
    # Boltzmann constant in atomic units
    kB = 3.167e-6  # Hartree/K
    kT = kB * target_temp
    
    # Extract data and use only the last half to exclude equilibration
    n_frames_total = velocity_data['n_frames']
    last_half_start = n_frames_total // 2
    
    velocities = velocity_data['velocities'][last_half_start:]
    vel_magnitudes = velocity_data['vel_magnitudes'][last_half_start:]
    vel_x = velocity_data['vel_x'][last_half_start:]
    vel_y = velocity_data['vel_y'][last_half_start:]
    vel_z = velocity_data['vel_z'][last_half_start:]
    kinetic_energies = velocity_data['kinetic_energies'][last_half_start:]
    temperatures_last_half = velocity_data['temperatures'][last_half_start:]
    time_last_half = velocity_data['time'][last_half_start:]
    mass = velocity_data['mass']
    n_frames = len(velocities)  # Number of frames in last half
    
    # Calculate statistics
    stats = calculate_velocity_statistics(velocities, np.full(len(velocities), mass), target_temp)
    
    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f'Cavity Particle Velocity Distribution Analysis - Last Half Data (T = {target_temp:.1f} K)', 
                 fontsize=16, fontweight='bold')
    
    # Plot 1: Velocity magnitude histogram
    ax1.hist(vel_magnitudes, bins=50, density=True, alpha=0.7, color='blue', 
             label=f'Cavity particle (N={n_frames:,})')
    
    # Plot theoretical Maxwell-Boltzmann distribution
    v_max = np.max(vel_magnitudes) if len(vel_magnitudes) > 0 else 4*np.sqrt(kT/mass)
    v_theory = np.linspace(0, v_max, 200)
    mb_3d = maxwell_boltzmann_3d(v_theory, mass, kT)
    ax1.plot(v_theory, mb_3d, 'k-', linewidth=2, 
             label=f'Maxwell-Boltzmann (T={target_temp:.1f} K)')
    
    # Add characteristic velocities
    v_most_probable = np.sqrt(2 * kT / mass)
    v_mean_theory = np.sqrt(8 * kT / (np.pi * mass))
    v_rms_theory = np.sqrt(3 * kT / mass)
    
    ax1.axvline(v_most_probable, color='red', linestyle='--', alpha=0.8, 
               label=f'Most probable: {v_most_probable:.4f}')
    ax1.axvline(stats['measured_v_mean'], color='green', linestyle='--', alpha=0.8, 
               label=f'Measured mean: {stats["measured_v_mean"]:.4f}')
    
    ax1.set_xlabel('Velocity Magnitude (a.u.)')
    ax1.set_ylabel('Probability Density')
    ax1.set_title('Cavity Particle: Velocity Magnitude Distribution')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Add statistics text
    stats_text = (f'Measured T: {stats["measured_temperature"]:.1f} K\n'
                 f'Target T: {target_temp:.1f} K\n'
                 f'Error: {stats["temp_error_percent"]:.1f}%\n'
                 f'Samples: {n_frames:,} (last half)\n'
                 f'Total frames: {n_frames_total:,}')
    ax1.text(0.98, 0.98, stats_text, transform=ax1.transAxes, 
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Plot 2: Velocity component histograms
    components = {'x': vel_x, 'y': vel_y, 'z': vel_z}
    colors = ['red', 'green', 'blue']
    
    for i, (comp, vel_comp) in enumerate(components.items()):
        ax2.hist(vel_comp, bins=40, density=True, alpha=0.5, 
                label=f'{comp}-component', color=colors[i])
    
    # Plot theoretical 1D Maxwell-Boltzmann for comparison
    v_comp_max = max(np.max(np.abs(vel_x)), np.max(np.abs(vel_y)), np.max(np.abs(vel_z)))
    v_comp_range = np.linspace(-v_comp_max, v_comp_max, 200)
    mb_1d = maxwell_boltzmann_1d(v_comp_range, mass, kT)
    ax2.plot(v_comp_range, mb_1d, 'k-', linewidth=2, 
            label=f'1D Maxwell-Boltzmann')
    
    ax2.set_xlabel('Velocity Component (a.u.)')
    ax2.set_ylabel('Probability Density')
    ax2.set_title('Cavity Particle: Velocity Components')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Kinetic energy distribution
    ax3.hist(kinetic_energies, bins=50, density=True, alpha=0.7, color='purple', 
            label='Cavity particle')
    
    # Theoretical kinetic energy distribution (1D oscillator)
    # For a 1D harmonic oscillator: P(E) ∝ exp(-E/kT)
    ke_max = np.max(kinetic_energies) if len(kinetic_energies) > 0 else kT * 10
    ke_theory = np.linspace(0, ke_max, 200)
    ke_dist = np.exp(-ke_theory / kT)
    ke_dist /= np.trapz(ke_dist, ke_theory)  # Normalize
    
    ax3.plot(ke_theory, ke_dist, 'k-', linewidth=2, 
            label=f'Exponential (1D, T={target_temp:.1f} K)')
    
    # Add mean kinetic energy lines
    mean_ke = np.mean(kinetic_energies)
    ax3.axvline(mean_ke, color='red', linestyle='--', alpha=0.8, 
               label=f'Mean KE: {mean_ke:.4f}')
    ax3.axvline(0.5 * kT, color='green', linestyle='--', alpha=0.8, 
               label=f'kT/2: {0.5 * kT:.4f}')
    
    ax3.set_xlabel('Kinetic Energy (a.u.)')
    ax3.set_ylabel('Probability Density')
    ax3.set_title('Cavity Particle: Kinetic Energy Distribution')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Temperature time series (showing last half data)
    ax4.plot(time_last_half, temperatures_last_half, 'b-', linewidth=1, alpha=0.7)
    ax4.axhline(y=target_temp, color='k', linestyle='--', alpha=0.7, 
               label=f'Target ({target_temp:.1f} K)')
    
    # Calculate statistics for the last half data
    temp_mean = np.mean(temperatures_last_half)
    temp_std = np.std(temperatures_last_half)
    
    ax4.axhline(y=temp_mean, color='red', linestyle=':', alpha=0.7, 
               label=f'Mean (last half): {temp_mean:.1f} K')
    
    ax4.set_xlabel('Time (ps)')
    ax4.set_ylabel('Temperature (K)')
    ax4.set_title('Cavity Particle: Temperature Time Series (Last Half)')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    # Add statistics text
    temp_stats_text = (f'Mean: {temp_mean:.1f} ± {temp_std:.1f} K\n'
                      f'Deviation: {temp_mean - target_temp:+.1f} K\n'
                      f'Error: {(temp_mean - target_temp)/target_temp*100:+.1f}%')
    ax4.text(0.02, 0.98, temp_stats_text, transform=ax4.transAxes, 
            verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_cavity_velocity_histogram.png', dpi=300, bbox_inches='tight')
    print(f"Saved cavity velocity histogram: {output_prefix}_cavity_velocity_histogram.png")
    
    # Print summary statistics
    print("\n" + "="*60)
    print("CAVITY VELOCITY DISTRIBUTION ANALYSIS (LAST HALF DATA)")
    print("="*60)
    print(f"Total data points: {n_frames_total:,}")
    print(f"Analysis period: last {n_frames:,} samples (excluding equilibration)")
    print(f"Time range analyzed: {time_last_half[0]:.3f} - {time_last_half[-1]:.3f} ps")
    print(f"Cavity particle mass: {mass:.6f} a.u.")
    print(f"Target temperature: {target_temp:.1f} K")
    print(f"Measured temperature: {stats['measured_temperature']:.1f} K")
    print(f"Temperature error: {stats['temp_error_percent']:.1f}%")
    print(f"Average kinetic energy: {stats['avg_kinetic_energy']:.6f} a.u.")
    print(f"Expected kinetic energy (1/2 kT): {0.5 * kT:.6f} a.u.")
    print(f"Velocity statistics:")
    print(f"  Measured mean: {stats['measured_v_mean']:.4f} a.u.")
    print(f"  Theoretical mean: {stats['theoretical_v_mean']:.4f} a.u.")
    print(f"  Measured RMS: {stats['measured_v_rms']:.4f} a.u.")
    print(f"  Theoretical RMS: {stats['theoretical_v_rms']:.4f} a.u.")
    print(f"  Most probable (theory): {stats['theoretical_v_most_probable']:.4f} a.u.")
    print("="*60)

def main():
    """Main function to run the energy analysis."""
    parser = argparse.ArgumentParser(description='Analyze and plot energy data from cavity MD simulations including individual potential and kinetic energy components')
    parser.add_argument('--file', '-f', type=str, required=True,
                       help='Path to the energy contributions file')
    parser.add_argument('--output', '-o', type=str, default='energy_analysis',
                       help='Output prefix for plot files (default: energy_analysis)')
    parser.add_argument('--max_time', '-t', type=float, default=None,
                       help='Maximum time to plot in ps (default: all data)')
    parser.add_argument('--stats', '-s', action='store_true',
                       help='Print detailed energy statistics including P.E./K.E. breakdown')
    parser.add_argument('--target_temp', type=float, default=100.0,
                       help='Target temperature in Kelvin for reference (default: 100.0)')
    parser.add_argument('--temp_only', action='store_true',
                       help='Only create temperature plots (skip energy plots)')
    parser.add_argument('--velocity_analysis', action='store_true',
                       help='Include velocity histogram analysis from GSD file')
    parser.add_argument('--gsd_file', type=str, default=None,
                       help='Path to GSD trajectory file for velocity analysis (auto-detected if not provided)')
    parser.add_argument('--cavity_velocity_only', action='store_true',
                       help='Analyze only cavity particle velocities')
    parser.add_argument('--molecular_velocity_only', action='store_true',
                       help='Analyze only molecular particle velocities')
    parser.add_argument('--cavity_velocity_analysis', action='store_true',
                       help='Include cavity velocity histogram analysis from cavity mode files')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.file):
        print(f"Error: File {args.file} not found!")
        return 1
    
    # Read energy data
    df = read_energy_data(args.file)
    if df is None:
        return 1
    
    # Calculate derived energies
    derived_energies = calculate_derived_energies(df)
    
    # Print statistics if requested
    if args.stats:
        print_energy_statistics(derived_energies, args.max_time, args.target_temp)
    
    # Create plots
    if not args.temp_only:
        create_energy_plots(df, derived_energies, args.output, args.max_time, args.target_temp)
    
    # Read temperature data
    temperature_data = read_temperature_data(args.file)
    
    # Create temperature plots
    create_temperature_plots(temperature_data, args.output, args.max_time, args.target_temp)
    
    # Print temperature statistics
    print_temperature_statistics(temperature_data, args.max_time, args.target_temp)
    
    # Velocity analysis if requested
    if args.velocity_analysis:
        # Auto-detect GSD file if not provided
        gsd_file = args.gsd_file
        if gsd_file is None:
            # Try to find GSD file based on energy file name
            base_path = Path(args.file)
            if '_energy_contributions.txt' in args.file:
                gsd_prefix = args.file.replace('_energy_contributions.txt', '')
            elif '_energy.txt' in args.file:
                gsd_prefix = args.file.replace('_energy.txt', '')
            else:
                gsd_prefix = str(base_path.stem)
            
            gsd_file = f"{gsd_prefix}.gsd"
        
        create_velocity_histograms(
            gsd_file, 
            temperature_K=args.target_temp, 
            output_prefix=args.output,
            cavity_only=args.cavity_velocity_only,
            molecular_only=args.molecular_velocity_only
        )
    
    # Cavity velocity analysis if requested
    if args.cavity_velocity_analysis:
        # Read cavity velocity data
        cavity_velocity_data = read_cavity_velocity_data(args.file)
        
        # Create cavity velocity histogram
        create_cavity_velocity_histogram(cavity_velocity_data, args.target_temp, args.output)
    
    print(f"\nEnergy analysis complete!")
    print(f"Generated plots:")
    print(f"  - Main components: {args.output}_main_components.png")
    print(f"  - Potential/Kinetic breakdown: {args.output}_potential_kinetic_components.png") 
    print(f"  - Overview: {args.output}_overview.png")
    print(f"  - Temperature time series: {args.output}_temperature_timeseries.png")
    if args.velocity_analysis:
        print(f"  - Velocity histograms: {args.output}_velocity_histograms.png")
    if args.cavity_velocity_analysis:
        print(f"  - Cavity velocity histogram: {args.output}_cavity_velocity_histogram.png")
    return 0

if __name__ == '__main__':
    exit(main()) 
