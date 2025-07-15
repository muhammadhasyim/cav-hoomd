#!/usr/bin/env python3
"""
Convert Dipole Autocorrelation to IR Spectrum
Based on the DCT calculation sample code

Takes a pre-calculated autocorrelation function and converts it to IR absorption spectrum.

Usage:
    python autocorr_to_ir.py

Input file format:
    # Header comments (optional)
    # timestep lag_time(ps) C(t)
    10888 0.000000 982.757471
    10889 0.000298 982.740510
    ...

Author: Based on provided DCT calculation sample
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
import argparse
import sys
import os

# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================

# Input/output files
input_autocorr_file = 'cavity_coupling_1eneg03_switch_10.0ps/prod-499_dipole_autocorr_1.txt'
output_data_file = 'IR-data.txt'

# Calculation parameters  
T = 300.0  # Temperature in Kelvin
zoom_wavenum = 4000.0  # Maximum wavenumber for zoomed plots (cm^-1)
interp_dt_fs = 0.1  # Interpolation timestep in femtoseconds for uniform spacing

# Physical constants (SI units)
BOLTZ = 1.38064852e-23  # Boltzmann constant (J/K)
LIGHTSPEED = 299792458.0  # Speed of light (m/s)
REDUCED_PLANCK = 1.05457180013e-34  # Reduced Planck constant (J·s)

def print_header():
    """Print script information."""
    print("=" * 60)
    print("Convert Dipole Autocorrelation to IR Spectrum")
    print("=" * 60)
    print(f"Temperature: {T:.1f} K")
    print(f"Input file: {input_autocorr_file}")
    print(f"Zoom range: 0-{zoom_wavenum:.0f} cm⁻¹")
    print("-" * 60)

def load_autocorrelation(filename):
    """Load autocorrelation function from file and interpolate to uniform spacing."""
    print(f"Loading autocorrelation from {filename}...")
    
    if not os.path.exists(filename):
        print(f"ERROR: File '{filename}' not found!")
        print("Expected format:")
        print("# Optional header")
        print("# timestep lag_time(ps) C(t)")
        print("10888 0.000000 982.757471")
        print("10889 0.000298 982.740510")
        print("...")
        sys.exit(1)
    
    try:
        # Read the file and automatically skip comment lines starting with #
        data = np.loadtxt(filename, comments='#')
        
        # Extract columns: timestep, lag_time(ps), C(t)
        if data.shape[1] != 3:
            raise ValueError(f"Expected 3 columns (timestep, lag_time(ps), C(t)), found {data.shape[1]}")
        
        timestep = data[:, 0]
        time_ps = data[:, 1]  # Time in picoseconds
        autocorr = data[:, 2]  # Autocorrelation C(t)
        
        # Convert time from picoseconds to femtoseconds
        time_fs = time_ps * 1000.0
        
        print(f"Loaded {len(autocorr)} autocorrelation points")
        print(f"Time range: {time_fs[0]:.2f} - {time_fs[-1]:.2f} fs")
        print(f"Time range (original): {time_ps[0]:.4f} - {time_ps[-1]:.4f} ps")
        if len(time_fs) > 1:
            print(f"Original timestep range: {np.min(np.diff(time_fs)):.4f} - {np.max(np.diff(time_fs)):.4f} fs")
        
        # Interpolate to uniform spacing
        print(f"Interpolating to uniform spacing with dt = {interp_dt_fs:.4f} fs...")
        time_uniform = np.arange(time_fs[0], time_fs[-1] + interp_dt_fs, interp_dt_fs)
        autocorr_interp = np.interp(time_uniform, time_fs, autocorr)
        
        print(f"Interpolated to {len(autocorr_interp)} points")
        print(f"Uniform timestep: {interp_dt_fs:.4f} fs")
        
        return time_uniform, autocorr_interp
    except Exception as e:
        print(f"ERROR loading autocorrelation: {e}")
        print("File format should be:")
        print("# Header lines (optional)")
        print("# timestep lag_time(ps) C(t)")
        print("10888 0.000000 982.757471")
        print("10889 0.000298 982.740510")
        print("...")
        sys.exit(1)

def calculate_ir_spectrum(time, autocorr):
    """
    Calculate IR absorption spectrum from autocorrelation function.
    
    This calculates the final IR spectrum with ω² prefactor.
    """
    print("Calculating IR spectrum from autocorrelation...")
    
    # Convert time from fs to seconds
    timestep = (time[1] - time[0]) * 1e-15
    
    # Calculate FFT of autocorrelation function
    # Use DCT type 1 as in the original code
    lineshape = fftpack.dct(autocorr, type=1)[1:]  # Skip DC component
    
    # Calculate frequency arrays
    lineshape_frequencies = np.linspace(0, 0.5/timestep, len(autocorr))[1:]
    lineshape_frequencies_wn = lineshape_frequencies / (100.0 * LIGHTSPEED)  # Convert to cm^-1
    
    # Calculate quantum correction factor
    quantum_correction = lineshape_frequencies / (1.0 - np.exp(-REDUCED_PLANCK * lineshape_frequencies / (BOLTZ * T)))
    
    # Calculate IR spectrum with ω² prefactor
    # IR intensity ∝ ω² × |F[μ(t)·μ(0)]|² × quantum correction
    omega_squared = (2.0 * np.pi * lineshape_frequencies)**2
    ir_spectrum = lineshape * omega_squared * quantum_correction
    
    print(f"Spectrum calculated for {len(lineshape_frequencies_wn)} frequency points")
    print(f"Frequency range: {lineshape_frequencies_wn[0]:.1f} - {lineshape_frequencies_wn[-1]:.1f} cm⁻¹")
    
    return lineshape_frequencies_wn, ir_spectrum

def save_spectrum_data(frequencies, ir_spectrum, filename):
    """Save IR spectrum data to file."""
    print(f"Saving IR spectrum data to {filename}...")
    
    header_line = 'Frequency(cm^-1) IR_Spectrum'
    
    data = np.column_stack((frequencies, ir_spectrum))
    np.savetxt(filename, data, header=header_line, fmt='%.6f %.6e')

def create_plots(time, autocorr, frequencies, ir_spectrum):
    """Create and save IR spectrum plot."""
    print("Creating IR spectrum plot...")
    
    # Create mask for zoom limit (upper frequency limit)
    mask = frequencies <= zoom_wavenum
    
    # Plot IR spectrum with zoom limit
    plt.figure(figsize=(10, 6))
    plt.plot(frequencies[mask], ir_spectrum[mask], 'b-', linewidth=1.5)
    plt.axvline(x=1560, color='k', linestyle='--', linewidth=1.5)
    plt.axvline(x=2331, color='k', linestyle='--', linewidth=1.5)
    plt.xlabel(r'Frequency (cm$^{-1}$)', fontsize=12)
    plt.ylabel(r'IR Intensity (arb. units)', fontsize=12)
    plt.title(f'IR Absorption Spectrum (T = {T:.0f} K)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, zoom_wavenum)
    plt.tight_layout()
    plt.savefig('IR-spectrum.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("IR spectrum plot saved as IR-spectrum.png")

def main():
    """Main calculation workflow."""
    print_header()
    
    # Load autocorrelation function
    time, autocorr = load_autocorrelation(input_autocorr_file)
    
    # Calculate IR spectrum
    spectrum_data = calculate_ir_spectrum(time, autocorr)
    frequencies, ir_spectrum = spectrum_data
    
    # Save spectrum data
    save_spectrum_data(frequencies, ir_spectrum, output_data_file)
    
    # Create plots
    create_plots(time, autocorr, frequencies, ir_spectrum)
    
    # Print summary
    print("-" * 60)
    print("CALCULATION COMPLETE!")
    print(f"Autocorrelation points: {len(autocorr)}")
    print(f"Spectrum points: {len(frequencies)}")
    print(f"Frequency range: {frequencies[0]:.1f} - {frequencies[-1]:.1f} cm⁻¹")
    print(f"Zoom range: 0 - {zoom_wavenum:.1f} cm⁻¹")
    print(f"Temperature: {T:.1f} K")
    print(f"Interpolation timestep: {interp_dt_fs:.4f} fs")
    print()
    print("Output files created:")
    print(f"  {output_data_file} - IR spectrum data")
    print("  IR-spectrum.png - IR absorption spectrum plot")
    print("=" * 60)

if __name__ == "__main__":
    # Command line argument support
    parser = argparse.ArgumentParser(description='Convert autocorrelation to IR spectrum')
    parser.add_argument('--input', '-i', 
                       default='cavity_coupling_1eneg03_switch_10.0ps/prod-499_dipole_autocorr_1.txt', 
                       help='Input autocorrelation file')
    parser.add_argument('--temp', '-T', type=float, default=300.0, help='Temperature in K')
    parser.add_argument('--zoom', '-z', type=float, default=4000.0, 
                       help='Max frequency for plot (cm^-1)')
    parser.add_argument('--dt', type=float, default=0.1, 
                       help='Interpolation timestep in fs')
    
    args = parser.parse_args()
    
    # Update parameters from command line
    input_autocorr_file = args.input
    T = args.temp
    zoom_wavenum = args.zoom
    interp_dt_fs = args.dt
    
    try:
        main()
    except KeyboardInterrupt:
        print("\nCalculation interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1) 
