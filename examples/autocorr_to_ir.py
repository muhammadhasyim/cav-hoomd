#!/usr/bin/env python3
"""
Convert Dipole Autocorrelation to IR Spectrum with Multi-Replica Averaging
Based on the DCT calculation sample code

Takes multiple pre-calculated autocorrelation functions with pattern:
prod-i_dipole_autocorr_j.txt (i=replica, j=reference) and averages them
to produce master averaged autocorrelation and IR spectrum.

Usage:
    python autocorr_to_ir.py --pattern "prod-*_dipole_autocorr_*.txt"

Input file format:
    # Header comments (optional)
    # timestep lag_time(ps) C(t)
    10888 0.000000 982.757471
    10889 0.000298 982.740510
    ...

Author: Based on provided DCT calculation sample, extended for multi-replica averaging
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack
import argparse
import sys
import os
import glob

# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================

# Global parameters
input_pattern = 'cavity_coupling_*/*_dipole_autocorr_*.txt'
T = 300.0  # Temperature in Kelvin
zoom_wavenum = 4000.0  # Max wavenumber for plots (cm^-1)
interp_dt_ps = 0.1  # Interpolation timestep in picoseconds for uniform spacing
output_dir = 'averaged_analysis'

# Physical constants (SI units)
BOLTZ = 1.38064852e-23  # Boltzmann constant (J/K)
LIGHTSPEED = 299792458.0  # Speed of light (m/s)
REDUCED_PLANCK = 1.05457180013e-34  # Reduced Planck constant (J·s)

def print_header():
    """Print script information."""
    print("=" * 60)
    print("Convert Multiple Dipole Autocorrelations to IR Spectrum")
    print("Multi-Replica and Multi-Reference Averaging")
    print("=" * 60)
    print(f"Temperature: {T:.1f} K")
    print(f"Input pattern: {input_pattern}")
    print(f"Zoom range: 0-{zoom_wavenum:.0f} cm⁻¹")
    print("-" * 60)

def find_autocorr_files(pattern):
    """Find all autocorrelation files matching the pattern."""
    print(f"Searching for files matching pattern: {pattern}")
    
    files = glob.glob(pattern)
    files.sort()  # Sort for consistent ordering
    
    if not files:
        print(f"ERROR: No files found matching pattern '{pattern}'!")
        print("Expected file naming pattern: prod-i_dipole_autocorr_j.txt")
        print("  where i = replica number, j = reference number")
        sys.exit(1)
    
    print(f"Found {len(files)} autocorrelation files:")
    for i, file in enumerate(files):
        print(f"  {i+1:3d}: {file}")
    
    return files

def load_single_autocorr(filename):
    """Load a single autocorrelation function from file."""
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
        
        return time_fs, autocorr
    except Exception as e:
        print(f"ERROR loading {filename}: {e}")
        return None, None

def load_and_average_autocorrelations(filenames):
    """Load multiple autocorrelation files and compute average."""
    print(f"\nLoading and averaging {len(filenames)} autocorrelation files...")
    
    all_autocorrs = []
    all_times = []
    successful_loads = 0
    
    # Load all files
    for i, filename in enumerate(filenames):
        print(f"  Loading {i+1}/{len(filenames)}: {os.path.basename(filename)}")
        
        time_fs, autocorr = load_single_autocorr(filename)
        if time_fs is not None and autocorr is not None:
            all_times.append(time_fs)
            all_autocorrs.append(autocorr)
            successful_loads += 1
        else:
            print(f"    WARNING: Failed to load {filename}")
    
    if successful_loads == 0:
        print("ERROR: No files loaded successfully!")
        sys.exit(1)
    
    print(f"Successfully loaded {successful_loads}/{len(filenames)} files")
    
    # Find common time range (intersection of all time ranges)
    min_time = max(times[0] for times in all_times)
    max_time = min(times[-1] for times in all_times)
    
    print(f"Common time range: {min_time:.2f} - {max_time:.2f} fs")
    
    # Create uniform time grid
    interp_dt_fs = interp_dt_ps * 1000.0  # Convert ps to fs
    time_uniform = np.arange(min_time, max_time + interp_dt_fs, interp_dt_fs)
    print(f"Interpolating to uniform grid with {len(time_uniform)} points (dt = {interp_dt_ps:.4f} ps = {interp_dt_fs:.1f} fs)")
    
    # Interpolate all autocorrelations to uniform grid
    interpolated_autocorrs = []
    for i, (time_fs, autocorr) in enumerate(zip(all_times, all_autocorrs)):
        # Only interpolate within the common time range
        mask = (time_fs >= min_time) & (time_fs <= max_time)
        if np.sum(mask) < 2:
            print(f"    WARNING: File {i+1} has insufficient data in common time range")
            continue
            
        autocorr_interp = np.interp(time_uniform, time_fs[mask], autocorr[mask])
        interpolated_autocorrs.append(autocorr_interp)
    
    if len(interpolated_autocorrs) == 0:
        print("ERROR: No autocorrelations could be interpolated to common grid!")
        sys.exit(1)
    
    # Average all interpolated autocorrelations
    autocorr_matrix = np.array(interpolated_autocorrs)
    averaged_autocorr = np.mean(autocorr_matrix, axis=0)
    autocorr_std = np.std(autocorr_matrix, axis=0)
    
    print(f"Averaged {len(interpolated_autocorrs)} autocorrelations")
    print(f"Average C(0): {averaged_autocorr[0]:.3f} ± {autocorr_std[0]:.3f}")
    print(f"Average C(final): {averaged_autocorr[-1]:.3f} ± {autocorr_std[-1]:.3f}")
    
    # Calculate IR spectrum
    print("Calculating IR spectrum for averaged autocorrelation...")
    frequencies, ir_spectrum = calculate_ir_spectrum(time_uniform, averaged_autocorr)
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    return time_uniform, averaged_autocorr, autocorr_std, frequencies, ir_spectrum

def calculate_ir_spectrum(time, autocorr):
    """
    Calculate IR absorption spectrum from autocorrelation function.
    
    This calculates the final IR spectrum with ω² prefactor.
    """
    print("\nCalculating IR spectrum from averaged autocorrelation...")
    
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

def save_averaged_autocorr(time, autocorr, autocorr_std, filename):
    """Save averaged autocorrelation data to file."""
    print(f"\nSaving averaged autocorrelation to {filename}...")
    
    # Convert time back to ps for output
    time_ps = time / 1000.0
    
    header_lines = [
        'Averaged Dipole Autocorrelation Function',
        f'Temperature: {T:.1f} K',
        f'Interpolation timestep: {interp_dt_ps:.4f} ps',
        f'Data points: {len(time)}',
        'time(fs) C(t) C_std(t)'
    ]
    header = '\n'.join([f'# {line}' for line in header_lines])
    
    data = np.column_stack((time_ps, autocorr, autocorr_std))
    np.savetxt(filename, data, header=header, fmt='%.6f %.6e %.6e')

def save_ir_spectrum(frequencies, ir_spectrum, filename):
    """Save IR spectrum data to file."""
    print(f"Saving IR spectrum to {filename}...")
    
    header_lines = [
        'IR Absorption Spectrum from Averaged Autocorrelation',
        f'Temperature: {T:.1f} K',
        'Columns: Frequency(cm^-1) IR_Intensity(arb.units)'
    ]
    header = '\n'.join([f'# {line}' for line in header_lines])
    
    data = np.column_stack((frequencies, ir_spectrum))
    np.savetxt(filename, data, header=header, fmt='%.6f %.6e')

def create_plots(time, autocorr, autocorr_std, frequencies, ir_spectrum):
    """Create and save plots for autocorrelation and IR spectrum."""
    print("\nCreating plots...")
    
    # Convert time to ps for plotting
    time_ps = time / 1000.0
    
    # Plot averaged autocorrelation with error bars
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(time_ps, autocorr, 'b-', linewidth=1.5, label='Averaged C(t)')
    plt.fill_between(time_ps, autocorr - autocorr_std, autocorr + autocorr_std, 
                     alpha=0.3, color='blue', label='±1σ')
    plt.xlabel('Lag Time (ps)', fontsize=12)
    plt.ylabel('Autocorrelation C(t)', fontsize=12)
    plt.title('Averaged Dipole Autocorrelation', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot IR spectrum with zoom limit
    plt.subplot(1, 2, 2)
    mask = frequencies <= zoom_wavenum
    plt.plot(frequencies[mask], ir_spectrum[mask], 'r-', linewidth=1.5)
    plt.axvline(x=1560, color='k', linestyle='--', linewidth=1.0, alpha=0.7, label='1560 cm⁻¹')
    plt.axvline(x=2331, color='k', linestyle='--', linewidth=1.0, alpha=0.7, label='2331 cm⁻¹')
    plt.xlabel(r'Frequency (cm$^{-1}$)', fontsize=12)
    plt.ylabel(r'IR Intensity (arb. units)', fontsize=12)
    plt.title(f'IR Spectrum (T = {T:.0f} K)', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim(0, zoom_wavenum)
    
    plt.tight_layout()
    plt.savefig('averaged_autocorr_and_IR.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Combined plot saved as averaged_autocorr_and_IR.png")

def main():
    """Main calculation workflow."""
    print_header()
    
    # Find all autocorrelation files
    autocorr_files = find_autocorr_files(input_pattern)
    
    # Load and average autocorrelations
    time, averaged_autocorr, autocorr_std, frequencies, ir_spectrum = load_and_average_autocorrelations(autocorr_files)
    
    # Save master data files
    save_averaged_autocorr(time, averaged_autocorr, autocorr_std, os.path.join(output_dir, 'averaged_dipole_autocorr.txt'))
    save_ir_spectrum(frequencies, ir_spectrum, os.path.join(output_dir, 'averaged_IR_spectrum.txt'))
    
    # Create plots
    create_plots(time, averaged_autocorr, autocorr_std, frequencies, ir_spectrum)
    
    # Print summary
    print("\n" + "=" * 60)
    print("MULTI-REPLICA AVERAGING COMPLETE!")
    print("=" * 60)
    print(f"Input files processed: {len(autocorr_files)}")
    print(f"Averaged autocorr points: {len(averaged_autocorr)}")
    print(f"IR spectrum points: {len(frequencies)}")
    print(f"Time range: {time[0]/1000:.4f} - {time[-1]/1000:.4f} ps")
    print(f"Frequency range: {frequencies[0]:.1f} - {frequencies[-1]:.1f} cm⁻¹")
    print(f"Temperature: {T:.1f} K")
    print(f"Interpolation timestep: {interp_dt_ps:.4f} ps")
    print()
    print("Master output files created:")
    print(f"  {os.path.join(output_dir, 'averaged_dipole_autocorr.txt')} - Averaged dipole autocorrelation")
    print(f"  {os.path.join(output_dir, 'averaged_IR_spectrum.txt')} - IR spectrum from averaged data")
    print("  averaged_autocorr_and_IR.png - Combined plots")
    print("=" * 60)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Average dipole autocorrelations and convert to IR spectrum')
    parser.add_argument('--pattern', '-p', type=str, 
                       default='cavity_coupling_*/*_dipole_autocorr_*.txt', 
                       help='File pattern for autocorrelation files')
    parser.add_argument('--temp', '-T', type=float, default=300.0, help='Temperature in K')
    parser.add_argument('--zoom', '-z', type=float, default=4000.0, 
                       help='Max frequency for plot (cm^-1)')
    parser.add_argument('--dt', type=float, default=0.1, 
                       help='Interpolation timestep in ps (picoseconds)')
    parser.add_argument('--output-dir', default='averaged_analysis',
                       help='Output directory for averaged results')
    
    args = parser.parse_args()
    
    # Update parameters from command line
    input_pattern = args.pattern
    T = args.temp
    zoom_wavenum = args.zoom
    interp_dt_ps = args.dt  # Now in picoseconds
    output_dir = args.output_dir
    
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
