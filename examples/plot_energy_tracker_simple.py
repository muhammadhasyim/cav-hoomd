#!/usr/bin/env python3
"""
Plot Energy Tracker Output (Simple Version)

This script reads the output from EnergyTracker and creates a two-panel plot:
1. Energy components vs time (molecular, cavity, reservoir, universe)
2. Temperature vs time

This version uses only numpy and matplotlib to avoid pandas compatibility issues.

Usage:
    python plot_energy_tracker_simple.py <energy_tracker_file>
    
Example:
    python plot_energy_tracker_simple.py cavity_coupling_1eneg03_switch_1.0ps/prod-0_energy_tracker.txt
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from pathlib import Path
from scipy.optimize import curve_fit

def load_energy_data(filename):
    """Load energy tracker data from file using numpy."""
    try:
        # First, read all lines to find the header
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Find the header line with column names (contains 'time(ps)' but is not a comment)
        header_line = None
        data_start_line = 0
        
        for i, line in enumerate(lines):
            line = line.strip()
            if not line.startswith('#') and 'time(ps)' in line:
                header_line = line
                data_start_line = i + 1
                break
        
        if header_line is None:
            raise ValueError("Could not find column headers in file")
        
        column_names = header_line.split()
        
        # Create a dictionary mapping column names to indices
        columns = {name: i for i, name in enumerate(column_names)}
        
        # Read the data starting from the line after the header
        data = np.loadtxt(filename, skiprows=data_start_line)
        
        return data, columns
        
    except Exception as e:
        print(f"Error loading data from {filename}: {e}")
        sys.exit(1)

def calculate_energy_components(data, columns):
    """Calculate the energy components we want to plot."""
    
    # Total molecular energy = molecular kinetic + molecular potential
    # Molecular potential = harmonic + lj + ewald_short + ewald_long
    molecular_potential = (data[:, columns['harmonic_energy']] + 
                          data[:, columns['lj_energy']] + 
                          data[:, columns['ewald_short_energy']] + 
                          data[:, columns['ewald_long_energy']])
    molecular_total = data[:, columns['molecular_kinetic_energy']] + molecular_potential
    
    # Total cavity energy = cavity kinetic + cavity potential
    cavity_total = (data[:, columns['cavity_kinetic_energy']] + 
                   data[:, columns['cavity_total_potential_energy']])
    
    # Total reservoir energy (already calculated in file)
    reservoir_total = data[:, columns['total_reservoir_energy']]
    
    # Universe energy (already calculated in file - should be conserved)
    universe_total = data[:, columns['universe_total_energy']]
    
    return molecular_total, cavity_total, reservoir_total, universe_total

def plot_energy_and_temperature(data, columns, output_file=None):
    """Create a two-panel plot of energy components and temperature."""
    
    # Calculate energy components
    molecular_total, cavity_total, reservoir_total, universe_total = calculate_energy_components(data, columns)
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
    
    # Panel 1: Energy components
    time_ps = data[:, columns['time(ps)']]
    
    ax1.plot(time_ps, molecular_total, label='Total Molecular Energy', linewidth=1, color='blue')
    ax1.plot(time_ps, cavity_total, label='Total Cavity Energy', linewidth=1, color='red')
    ax1.plot(time_ps, reservoir_total, label='Total Reservoir Energy', linewidth=1, color='green')
    ax1.plot(time_ps, universe_total, label='Universe Energy (Conserved)', linewidth=1, color='black', linestyle='-')
    # Plot total system energy
    ax1.plot(time_ps, molecular_total + cavity_total, label='Total System Energy', linewidth=1, color='purple', linestyle='-')
    
    # Fit exponential decay for time > 100 ps
    mask = time_ps > 100
    time_fit = time_ps[mask]
    energy_fit = (universe_total)[mask]  # Total system energy
    
    # Only do exponential fit if we have enough data points
    if len(time_fit) > 100:  # Require at least 100 points for meaningful fit
        # Shift time to start from 0
        t_shifted = time_fit - time_fit[0]
        
        # Define exponential decay function
        def exp_decay(t, A, tau, b):
            """Exponential decay to constant: A*exp(-t/tau) + b"""
            return A * np.exp(-t/tau) + b
        
        # Initial parameter guesses
        A_init = np.abs(energy_fit[0] - energy_fit[-1])  # Total energy change
        tau_init = (time_fit[-1] - time_fit[0]) / 3.0  # 1/3 of time range
        b_init = np.mean(energy_fit[-50:])  # Final energy level
        
        p0 = [A_init, tau_init, b_init]
        
        try:
            # Fit exponential decay
            popt, pcov = curve_fit(exp_decay, t_shifted, energy_fit, p0=p0)
            A, tau, b = popt
            
            print(f"\nExponential fit results:")
            print(f"Amplitude (A): {A:.2e} Hartree")
            print(f"Time constant (tau): {tau:.2f} ps")
            print(f"Final energy (b): {b:.6f} Hartree")
            
            # Linear fit for first 100 points
            t_linear = t_shifted#[#:200]
            e_linear = energy_fit#[:200]
            
            # Fit linear function E(t) = E0 + (dE/dt)*t
            slope, intercept = np.polyfit(t_linear, e_linear, 1)
            
            # Calculate linear time constant (time to decay by 1/e)
            linear_tau = (1-np.exp(-1))*intercept/slope if slope != 0 else np.inf
            
            print(f"\nLinear short-time fit results:")
            print(f"Initial slope (dE/dt): {slope:.2e} Hartree/ps")
            print(f"Linear time constant: {linear_tau:.2f} ps")
            
            # Generate fitted curves
            fit_line = exp_decay(t_shifted, A, tau, b)
            linear_fit = slope*t_linear + intercept
            
            # Plot exponential fit
            ax1.plot(time_fit, fit_line, '--', color='red', alpha=0.8, linewidth=2,
                     label=f'Exp Fit (τ={tau:.1f} ps)')
                     
            # Plot linear fit
            ax1.plot(time_fit, linear_fit, '--', color='green', alpha=0.8, linewidth=2,
                     label=f'Linear Fit (τ={linear_tau:.1f} ps)')
            
            # Create inset plot showing first 3 tau of decay
            # Position inset at upper right of main plot
            ax_inset = ax1.inset_axes([0.55, 0.55, 0.4, 0.4])
            
            # Plot data and fit for first 3 tau
            t_max_inset = time_fit[0] + 10*tau
            mask_inset = time_fit <= t_max_inset
            
            # Plot actual data
            ax_inset.plot(time_fit[mask_inset], energy_fit[mask_inset], 
                         color='black', linewidth=1, label='Data')
            
            # Plot exponential fit
            ax_inset.plot(time_fit[mask_inset], fit_line[mask_inset],
                         '--', color='red', alpha=0.8, linewidth=2, 
                         label='Exp Fit')
                         
            # Plot linear fit (only first 100 points)
            ax_inset.plot(time_fit[:100], linear_fit,
                         '--', color='green', alpha=0.8, linewidth=2,
                         label='Linear Fit')
            
            # Add grid and legend
            ax_inset.grid(True, alpha=0.3)
            ax_inset.legend(fontsize=8)
            
            # Set title showing zoom range
            ax_inset.set_title(f'First 10τ ({10*tau:.1f} ps)', fontsize=10)
            
            # Mark inset region on main plot
            ax1.indicate_inset_zoom(ax_inset, edgecolor="black")
        except Exception as e:
            print(f"\nFitting failed: {e}")
            print("Skipping fit analysis...")
            
        ax1.set_ylabel('Energy (Hartree)', fontsize=12)
        ax1.set_title('Energy Components vs Time', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='upper left', fontsize=9)
    else:
        print("\nWarning: Not enough data points after t=100 ps for exponential fit")
        tau = None  # Signal to later code that fit wasn't performed
   
    
    # Add text box with energy conservation info
    energy_drift = universe_total[-1] - universe_total[0]
    energy_std = np.std(universe_total)
    textstr = f'Universe Energy Drift: {energy_drift:.2e} Hartree\nStd Dev: {energy_std:.2e} Hartree'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=9,
             verticalalignment='top', bbox=props)
    
    # Panel 2: Temperature
    temperature = data[:, columns['temperature']]
    ax2.plot(time_ps, temperature, label='System Temperature', linewidth=1, color='orange', linestyle='-')
    
    ax2.set_xlabel('Time (ps)', fontsize=12)
    ax2.set_ylabel('Temperature (K)', fontsize=12)
    ax2.set_title('Temperature vs Time', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # Add horizontal line at target temperature if it's roughly constant
    temp_mean = np.mean(temperature)
    temp_std = np.std(temperature)
    ax2.axhline(y=temp_mean, color='red', linestyle=':', alpha=0.7, 
                label=f'Mean: {temp_mean:.1f} ± {temp_std:.1f} K')
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

def print_energy_statistics(data, columns):
    """Print useful statistics about the energy data."""
    
    molecular_total, cavity_total, reservoir_total, universe_total = calculate_energy_components(data, columns)
    time_ps = data[:, columns['time(ps)']]
    temperature = data[:, columns['temperature']]
    
    print("\n" + "="*60)
    print("ENERGY STATISTICS")
    print("="*60)
    
    print(f"Simulation time: {time_ps[0]:.6f} - {time_ps[-1]:.6f} ps")
    print(f"Number of data points: {len(data)}")
    
    print(f"\nEnergy Components (Hartree):")
    print(f"  Molecular Energy:")
    print(f"    Initial: {molecular_total[0]:.6f}")
    print(f"    Final:   {molecular_total[-1]:.6f}")
    print(f"    Change:  {molecular_total[-1] - molecular_total[0]:.6f}")
    
    print(f"  Cavity Energy:")
    print(f"    Initial: {cavity_total[0]:.6f}")
    print(f"    Final:   {cavity_total[-1]:.6f}")
    print(f"    Change:  {cavity_total[-1] - cavity_total[0]:.6f}")
    
    print(f"  Reservoir Energy:")
    print(f"    Initial: {reservoir_total[0]:.6f}")
    print(f"    Final:   {reservoir_total[-1]:.6f}")
    print(f"    Change:  {reservoir_total[-1] - reservoir_total[0]:.6f}")
    
    print(f"  Universe Energy (should be conserved):")
    print(f"    Initial: {universe_total[0]:.6f}")
    print(f"    Final:   {universe_total[-1]:.6f}")
    print(f"    Drift:   {universe_total[-1] - universe_total[0]:.2e}")
    print(f"    Std Dev: {np.std(universe_total):.2e}")
    
    print(f"\nTemperature Statistics:")
    print(f"  Mean: {np.mean(temperature):.2f} K")
    print(f"  Std:  {np.std(temperature):.2f} K")
    print(f"  Min:  {np.min(temperature):.2f} K")
    print(f"  Max:  {np.max(temperature):.2f} K")
    
    print("="*60)

def main():
    parser = argparse.ArgumentParser(
        description='Plot energy tracker output with energy components and temperature',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('input_file', type=str, 
                       help='Energy tracker output file (e.g., prod-0_energy_tracker.txt)')
    parser.add_argument('--output', '-o', type=str, 
                       help='Output plot file (default: show plot)')
    parser.add_argument('--stats', action='store_true', 
                       help='Print detailed energy statistics')
    
    args = parser.parse_args()
    
    # Check if input file exists
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file '{args.input_file}' not found")
        sys.exit(1)
    
    print(f"Loading energy data from: {args.input_file}")
    
    # Load data
    data, columns = load_energy_data(args.input_file)
    
    print(f"Loaded {len(data)} data points")
    print(f"Time range: {data[0, columns['time(ps)']]:.6f} - {data[-1, columns['time(ps)']]:.6f} ps")
    
    # Print statistics if requested
    if args.stats:
        print_energy_statistics(data, columns)
    
    # Create plot
    print("Creating plot...")
    fig = plot_energy_and_temperature(data, columns, args.output)
    
    if not args.output:
        print("Close the plot window to exit.")

if __name__ == '__main__':
    main() 