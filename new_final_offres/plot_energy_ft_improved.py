#!/usr/bin/env python3
"""
Improved script to plot Fourier transforms of energy components.
Offers multiple analysis options to avoid ω² weighting issues.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path


def extract_coupling_strength(folder_name):
    """Extract coupling strength from folder name for sorting."""
    match = re.search(r'cavity_coupling_(.+)', folder_name)
    if match:
        coupling_str = match.group(1)
        if 'eneg' in coupling_str:
            coupling_str = coupling_str.replace('eneg', 'e-')
        elif 'epos' in coupling_str:
            coupling_str = coupling_str.replace('epos', 'e+')
        try:
            return float(coupling_str)
        except ValueError:
            return 0.0
    return 0.0


def read_averaged_energy_file(filepath):
    """Read averaged energy file and return time and energy data."""
    try:
        data = np.loadtxt(filepath, comments='#')
        time = data[:, 0]
        return time, data, True
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None, None, False


def calculate_improved_spectrum(time, energy, analysis_type="power_spectrum", freq_max_cm=3500):
    """
    Calculate improved spectral analysis avoiding ω² weighting issues.
    
    Parameters:
    -----------
    time : np.array
        Time array in ps
    energy : np.array
        Energy time series
    analysis_type : str
        Type of analysis: "power_spectrum", "velocity_autocorr"
    freq_max_cm : float
        Maximum frequency in cm⁻¹
        
    Returns:
    --------
    freq_cm : np.array
        Frequency array in cm⁻¹
    spectrum : np.array
        Spectral density
    """
    # Convert time from ps to seconds
    time_s = time * 1e-12
    dt = time_s[1] - time_s[0]
    
    # Remove DC component
    energy_fluct = energy - np.mean(energy)
    
    if analysis_type == "power_spectrum":
        # Standard power spectrum (no ω² weighting)
        window = np.hanning(len(energy_fluct))
        energy_windowed = energy_fluct * window
        
        fft_energy = np.fft.fft(energy_windowed)
        freqs_hz = np.fft.fftfreq(len(time), dt)
        
        positive_freqs = freqs_hz[:len(freqs_hz)//2]
        fft_positive = fft_energy[:len(fft_energy)//2]
        
        power_spectrum = np.abs(fft_positive)**2
        spectrum = power_spectrum
        
    elif analysis_type == "velocity_autocorr":
        # Calculate energy "velocity" (time derivative) and its power spectrum
        energy_velocity = np.gradient(energy_fluct, time_s)
        
        window = np.hanning(len(energy_velocity))
        velocity_windowed = energy_velocity * window
        
        fft_velocity = np.fft.fft(velocity_windowed)
        freqs_hz = np.fft.fftfreq(len(time), dt)
        
        positive_freqs = freqs_hz[:len(freqs_hz)//2]
        fft_positive = fft_velocity[:len(fft_velocity)//2]
        
        # Power spectrum of velocity ~ ω² × power spectrum of position
        velocity_power = np.abs(fft_positive)**2
        spectrum = velocity_power
        
    else:
        raise ValueError(f"Unknown analysis type: {analysis_type}")
    
    # Convert frequency to cm⁻¹
    c_light = 2.998e8  # m/s
    freq_cm = positive_freqs / (c_light * 100)
    
    # Filter to desired frequency range
    freq_mask = (freq_cm >= 0) & (freq_cm <= freq_max_cm)
    freq_cm_filtered = freq_cm[freq_mask]
    spectrum_filtered = spectrum[freq_mask]
    
    # Normalize by area under the curve
    if len(spectrum_filtered) > 1 and np.trapezoid(spectrum_filtered, freq_cm_filtered) > 0:
        area = np.trapezoid(spectrum_filtered, freq_cm_filtered)
        spectrum_normalized = spectrum_filtered / area
    else:
        spectrum_normalized = spectrum_filtered
    
    return freq_cm_filtered, spectrum_normalized


def plot_improved_energy_spectra(base_dir='.', analysis_type="power_spectrum"):
    """
    Plot improved spectral analysis of energy components.
    
    Parameters:
    -----------
    base_dir : str
        Base directory to search for coupling folders
    analysis_type : str
        Type of spectral analysis
    """
    # Find valid folders
    pattern = os.path.join(base_dir, "cavity_coupling_*")
    folders = glob.glob(pattern)
    
    if not folders:
        print("No coupling strength folders found!")
        return
    
    folders.sort(key=lambda x: extract_coupling_strength(os.path.basename(x)))
    
    valid_folders = []
    for folder in folders:
        energy_file = os.path.join(folder, "prod_energy_tracker_averaged.txt")
        if os.path.exists(energy_file):
            valid_folders.append(folder)
    
    if not valid_folders:
        print("No folders with averaged energy files found!")
        return
    
    print(f"Analysis type: {analysis_type}")
    print(f"Found {len(valid_folders)} folders with averaged energy data:")
    for folder in valid_folders:
        coupling = extract_coupling_strength(os.path.basename(folder))
        print(f"  - {os.path.basename(folder)} (coupling = {coupling:.2e})")
    
    # Calculate layout
    n_folders = len(valid_folders)
    n_cols = min(3, n_folders)
    n_rows = (n_folders + n_cols - 1) // n_cols
    
    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    
    if n_folders == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # Colors for energy components
    colors = {
        'molecular_total': '#1f77b4',
        'cavity_total': '#ff7f0e',
        'cavity_reservoir': '#2ca02c',
        'molecular_reservoir': '#d62728',
        'universe': '#9467bd'
    }
    
    # Process each folder
    for i, folder in enumerate(valid_folders):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row][col]
        
        # Read energy data
        energy_file = os.path.join(folder, "prod_energy_tracker_averaged.txt")
        time, data, success = read_averaged_energy_file(energy_file)
        
        if not success or data is None or time is None:
            ax.text(0.5, 0.5, f'Error reading\n{os.path.basename(folder)}', 
                   ha='center', va='center', transform=ax.transAxes)
            continue
        
        # Extract energy components
        molecular_total = data[:, 10] + data[:, 2] + data[:, 3] + data[:, 4] + data[:, 5]
        cavity_total = data[:, 11] + data[:, 9]
        cavity_reservoir = data[:, 16]
        molecular_reservoir = data[:, 15]
        universe = data[:, 18]
        
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
        
        # Plot spectra for each energy component
        for energy_type, energy_data in energies.items():
            try:
                freq_cm, spectrum = calculate_improved_spectrum(time, energy_data, analysis_type)
                if len(freq_cm) > 0 and len(spectrum) > 0:
                    ax.semilogy(freq_cm, spectrum, label=labels[energy_type], 
                               color=colors[energy_type], linewidth=1.5, alpha=0.8)
            except Exception as e:
                print(f"Error calculating spectrum for {energy_type} in {folder}: {e}")
                continue
        
        # Set labels and title
        ax.set_xlabel('Frequency (cm⁻¹)')
        
        if analysis_type == "power_spectrum":
            ax.set_ylabel('Power Spectral Density (normalized)')
        elif analysis_type == "velocity_autocorr":
            ax.set_ylabel('dE/dt Power Spectrum (normalized)')
        
        coupling = extract_coupling_strength(os.path.basename(folder))
        ax.set_title(f'Coupling = {coupling:.2e}', fontsize=12, fontweight='bold')
        
        if i == 0:
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
        
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 3500)
        ax.set_ylim(bottom=1e-6)  # Set reasonable lower limit for log scale
    
    # Hide empty subplots
    for i in range(n_folders, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row][col].set_visible(False)
    
    plt.tight_layout()
    
    # Title based on analysis type
    if analysis_type == "power_spectrum":
        title = 'Energy Power Spectral Density vs Frequency'
    elif analysis_type == "velocity_autocorr":
        title = 'Energy Velocity Power Spectrum vs Frequency'
    
    fig.suptitle(title, fontsize=16, fontweight='bold', y=0.98)
    plt.subplots_adjust(top=0.93)
    
    # Save figure
    output_file = f"energy_spectra_{analysis_type}_all_couplings.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved as: {output_file}")
    
    output_pdf = f"energy_spectra_{analysis_type}_all_couplings.pdf"
    plt.savefig(output_pdf, bbox_inches='tight')
    print(f"Figure also saved as: {output_pdf}")
    
    plt.show()


def main():
    """Main function with different analysis options."""
    import sys
    
    if len(sys.argv) > 1:
        analysis_type = sys.argv[1]
    else:
        analysis_type = "power_spectrum"  # Default
    
    valid_types = ["power_spectrum", "velocity_autocorr"]
    if analysis_type not in valid_types:
        print(f"Invalid analysis type. Choose from: {valid_types}")
        return
    
    plot_improved_energy_spectra(analysis_type=analysis_type)


if __name__ == "__main__":
    main() 