#!/usr/bin/env python3
"""
Script to plot Fourier transforms of energy components with ω² weighting.
Creates a multi-panel figure showing FT(E(t))*ω² vs frequency in cm⁻¹.

Usage: python plot_energy_ft.py
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
        Time array in ps
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


def calculate_ft_spectrum(time, energy, freq_max_cm=3500):
    """
    Calculate Fourier transform spectrum with ω² weighting.
    
    Parameters:
    -----------
    time : np.array
        Time array in ps
    energy : np.array
        Energy time series
    freq_max_cm : float
        Maximum frequency in cm⁻¹ for plotting
        
    Returns:
    --------
    freq_cm : np.array
        Frequency array in cm⁻¹
    spectrum : np.array
        FT(E(t))*ω² spectrum, normalized by area
    """
    # Convert time from ps to seconds
    time_s = time * 1e-12
    
    # Calculate time step
    dt = time_s[1] - time_s[0]
    
    # Remove DC component (mean)
    energy_fluct = energy - np.mean(energy)
    
    # Apply window function to reduce spectral leakage
    window = np.hanning(len(energy_fluct))
    energy_windowed = energy_fluct * window
    
    # Calculate FFT
    fft_energy = np.fft.fft(energy_windowed)
    
    # Calculate frequency array in Hz
    freqs_hz = np.fft.fftfreq(len(time), dt)
    
    # Take only positive frequencies
    positive_freqs = freqs_hz[:len(freqs_hz)//2]
    fft_positive = fft_energy[:len(fft_energy)//2]
    
    # Convert frequency to cm⁻¹
    # 1 cm⁻¹ = c × 100 Hz, where c is speed of light in m/s
    c_light = 2.998e8  # m/s
    freq_cm = positive_freqs / (c_light * 100)
    
    # Calculate power spectrum
    power_spectrum = np.abs(fft_positive)**2
    
    # Apply ω² weighting (ω = 2π × frequency in Hz)
    omega = 2 * np.pi * positive_freqs
    omega_squared = omega**2
    
    # Avoid division by zero at DC
    omega_squared[0] = 0
    
    # Apply ω² weighting
    weighted_spectrum = power_spectrum * omega_squared
    
    # Filter to desired frequency range
    freq_mask = (freq_cm >= 0) & (freq_cm <= freq_max_cm)
    freq_cm_filtered = freq_cm[freq_mask]
    spectrum_filtered = weighted_spectrum[freq_mask]
    
    # Normalize by area under the curve
    if len(spectrum_filtered) > 1:
        area = np.trapezoid(spectrum_filtered, freq_cm_filtered)
        if float(area) > 0:
            spectrum_normalized = spectrum_filtered / area
        else:
            spectrum_normalized = spectrum_filtered
    else:
        spectrum_normalized = spectrum_filtered
    
    return freq_cm_filtered, spectrum_normalized


def plot_energy_ft_spectra(base_dir='.'):
    """
    Plot FT spectra of energy components for all detected coupling strength folders.
    
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
        
        # Extract energy components (same as before)
        # Molecular total energy = molecular KE + molecular PE (harmonic + LJ + ewald)
        molecular_total = data[:, 10] + data[:, 2] + data[:, 3] + data[:, 4] + data[:, 5]
        
        # Cavity total energy = cavity KE + cavity PE
        cavity_total = data[:, 11] + data[:, 9]
        
        cavity_reservoir = data[:, 16]  # cavity_reservoir_energy
        molecular_reservoir = data[:, 15]  # molecular_reservoir_energy
        universe = data[:, 18]  # universe_total_energy
        
        # Calculate FT spectra for each energy component
        energies = {
            'molecular_total': molecular_total,
            'cavity_total': cavity_total,
            'cavity_reservoir': cavity_reservoir,
            'molecular_reservoir': molecular_reservoir,
            'universe': universe
        }
        
        labels = {
            'molecular_total': 'Molecular Total',
            'cavity_total': 'Cavity Total',
            'cavity_reservoir': 'Cavity Reservoir',
            'molecular_reservoir': 'Molecular Reservoir',
            'universe': 'Universe'
        }
        
        # Plot FT spectra
        for energy_type, energy_data in energies.items():
            try:
                freq_cm, spectrum = calculate_ft_spectrum(time, energy_data)
                if len(freq_cm) > 0 and len(spectrum) > 0:
                    ax.plot(freq_cm, spectrum, label=labels[energy_type], 
                           color=colors[energy_type], linewidth=1.5, alpha=0.8)
            except Exception as e:
                print(f"Error calculating FT for {energy_type} in {folder}: {e}")
                continue
        
        # Set labels and title
        ax.set_xlabel('Frequency (cm⁻¹)')
        ax.set_ylabel('FT(E(t))×ω² (normalized)')
        
        # Extract coupling strength for title
        coupling = extract_coupling_strength(os.path.basename(folder))
        ax.set_title(f'Coupling = {coupling:.2e}', fontsize=12, fontweight='bold')
        
        # Add legend (only for first subplot to avoid clutter)
        if i == 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
        
        # Add grid for better readability
        ax.grid(True, alpha=0.3)
        
        # Set axis limits
        ax.set_xlim(0, 3500)
        ax.set_ylim(bottom=0)
        
        # Use log scale for y-axis if needed
        ax.set_yscale('log')
    
    # Hide empty subplots
    for i in range(n_folders, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row][col].set_visible(False)
    
    # Adjust layout
    plt.tight_layout()
    
    # Add overall title
    fig.suptitle('Energy FT Spectra: FT(E(t))×ω² vs Frequency', 
                fontsize=16, fontweight='bold', y=0.98)
    
    # Adjust layout to accommodate suptitle
    plt.subplots_adjust(top=0.93)
    
    # Save figure
    output_file = "energy_ft_spectra_all_couplings.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved as: {output_file}")
    
    # Also save as PDF
    output_pdf = "energy_ft_spectra_all_couplings.pdf"
    plt.savefig(output_pdf, bbox_inches='tight')
    print(f"Figure also saved as: {output_pdf}")
    
    plt.show()


def main():
    """Main function to run the FT plotting script."""
    plot_energy_ft_spectra()


if __name__ == "__main__":
    main() 