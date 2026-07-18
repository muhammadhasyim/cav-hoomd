#!/usr/bin/env python3
"""
Script to plot energy components for all detected coupling strength folders.
Creates a multi-panel figure with each panel showing:
1. Molecular total energy
2. Cavity total energy  
3. Cavity reservoir energy
4. Molecular reservoir energy
5. Universe energy

Usage: python plot_energy_components.py
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path


def extract_coupling_strength(folder_name):
    """
    Extract coupling strength from folder name for sorting.
    
    Parameters:
    -----------
    folder_name : str
        Folder name like 'cavity_coupling_1eneg04'
        
    Returns:
    --------
    float
        Coupling strength value for sorting
    """
    # Extract the coupling strength pattern
    match = re.search(r'cavity_coupling_(.+)', folder_name)
    if match:
        coupling_str = match.group(1)
        
        # Convert scientific notation patterns
        if 'eneg' in coupling_str:
            # Handle patterns like '1eneg04' -> 1e-04
            coupling_str = coupling_str.replace('eneg', 'e-')
        elif 'epos' in coupling_str:
            # Handle patterns like '1epos04' -> 1e+04
            coupling_str = coupling_str.replace('epos', 'e+')
        
        try:
            return float(coupling_str)
        except ValueError:
            return 0.0
    return 0.0


def read_averaged_energy_file(filepath):
    """
    Read averaged energy file and return time and energy data.
    
    Parameters:
    -----------
    filepath : str
        Path to the averaged energy file
        
    Returns:
    --------
    time : np.array
        Time array
    data : np.array
        2D array of energy data
    success : bool
        Whether reading was successful
    """
    try:
        # Read the file, skipping header lines
        data = np.loadtxt(filepath, comments='#')
        
        # Extract time (first column)
        time = data[:, 0]
        
        return time, data, True
        
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None, None, False


def plot_energy_components(base_dir='.'):
    """
    Plot energy components for all detected coupling strength folders.
    
    Parameters:
    -----------
    base_dir : str
        Base directory to search for coupling folders
    """
    # Find all coupling strength folders
    pattern = os.path.join(base_dir, "cavity_coupling_*")
    folders = glob.glob(pattern)
    
    if not folders:
        print("No coupling strength folders found!")
        return
    
    # Sort folders by coupling strength
    folders.sort(key=lambda x: extract_coupling_strength(os.path.basename(x)))
    
    # Check which folders have averaged energy files
    valid_folders = []
    for folder in folders:
        energy_file = os.path.join(folder, "prod_energy_tracker_averaged.txt")
        if os.path.exists(energy_file):
            valid_folders.append(folder)
        else:
            print(f"Warning: {folder} does not have averaged energy file")
    
    if not valid_folders:
        print("No folders with averaged energy files found!")
        return
    
    print(f"Found {len(valid_folders)} folders with averaged energy data:")
    for folder in valid_folders:
        coupling = extract_coupling_strength(os.path.basename(folder))
        print(f"  - {os.path.basename(folder)} (coupling = {coupling:.2e})")
    
    # Calculate subplot layout
    n_folders = len(valid_folders)
    n_cols = min(3, n_folders)  # Maximum 3 columns
    n_rows = (n_folders + n_cols - 1) // n_cols
    
    # Create figure and subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    
    # Ensure axes is always 2D for consistent indexing
    if n_folders == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # Colors for different energy components
    colors = {
        'molecular_total': '#1f77b4',      # Blue
        'cavity_total': '#ff7f0e',        # Orange
        'cavity_reservoir': '#2ca02c',    # Green
        'molecular_reservoir': '#d62728', # Red
        'universe': '#9467bd'             # Purple
    }
    
    # Process each folder
    for i, folder in enumerate(valid_folders):
        row = i // n_cols
        col = i % n_cols
        
        # Get the correct axis (axes is always 2D now)
        ax = axes[row][col]
        
        # Read energy data
        energy_file = os.path.join(folder, "prod_energy_tracker_averaged.txt")
        time, data, success = read_averaged_energy_file(energy_file)
        
        if not success or data is None or time is None:
            ax.text(0.5, 0.5, f'Error reading\n{os.path.basename(folder)}', 
                   ha='center', va='center', transform=ax.transAxes)
            continue
        
        # Extract energy components (adjust indices for 0-based indexing)
        # Molecular total energy = molecular KE + molecular PE (harmonic + LJ + ewald)
        molecular_total = data[:, 10] + data[:, 2] + data[:, 3] + data[:, 4] + data[:, 5]
        # That's: molecular_kinetic + harmonic + lj + ewald_short + ewald_long
        
        # Cavity total energy = cavity KE + cavity PE
        cavity_total = data[:, 11] + data[:, 9]  # cavity_kinetic + cavity_total_potential
        
        cavity_reservoir = data[:, 16]  # cavity_reservoir_energy (column 17)
        molecular_reservoir = data[:, 15]  # molecular_reservoir_energy (column 16)
        universe = data[:, 18]  # universe_total_energy (column 19)
        
        # Plot energy components
        ax.plot(time, molecular_total, label='Molecular Total', 
                color=colors['molecular_total'], linewidth=1.5)
        ax.plot(time, cavity_total, label='Cavity Total', 
                color=colors['cavity_total'], linewidth=1.5)
        ax.plot(time, cavity_reservoir, label='Cavity Reservoir', 
                color=colors['cavity_reservoir'], linewidth=1.5)
        ax.plot(time, molecular_reservoir, label='Molecular Reservoir', 
                color=colors['molecular_reservoir'], linewidth=1.5)
        ax.plot(time, universe, label='Universe', 
                color=colors['universe'], linewidth=1.5, linestyle='--')
        
        # Set labels and title
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('Energy (Hartree)')
        
        # Extract coupling strength for title
        coupling = extract_coupling_strength(os.path.basename(folder))
        ax.set_title(f'Coupling = {coupling:.2e}', fontsize=12, fontweight='bold')
        
        # Add legend (only for first subplot to avoid clutter)
        if i == 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        
        # Add grid for better readability
        ax.grid(True, alpha=0.3)
        
        # Set axis limits for better visualization
        ax.set_xlim(time.min(), time.max())
    
    # Hide empty subplots
    for i in range(n_folders, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row][col].set_visible(False)
    
    # Adjust layout
    plt.tight_layout()
    
    # Add overall title
    fig.suptitle('Energy Components vs Time for Different Coupling Strengths', 
                fontsize=16, fontweight='bold', y=0.98)
    
    # Adjust layout to accommodate suptitle
    plt.subplots_adjust(top=0.93)
    
    # Save figure
    output_file = "energy_components_all_couplings.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved as: {output_file}")
    
    # Also save as PDF
    output_pdf = "energy_components_all_couplings.pdf"
    plt.savefig(output_pdf, bbox_inches='tight')
    print(f"Figure also saved as: {output_pdf}")
    
    plt.show()


def main():
    """Main function to run the plotting script."""
    plot_energy_components()


if __name__ == "__main__":
    main() 