#!/usr/bin/env python3
"""
Script to plot energy power spectral density using FFT.
Creates individual figures for each coupling strength with zero coupling as background reference.
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Set non-interactive backend before importing pyplot


import re
from pathlib import Path
from scipy.signal import find_peaks
from scipy.interpolate import interp1d


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


def calculate_proper_energy_spectrum(time, energy, T=100.0, freq_max_cm=None):
    """
    Calculate power spectral density using FFT.
    
    Parameters:
    -----------
    time : np.array
        Time array in ps
    energy : np.array
        Energy time series
    T : float
        Temperature in K (kept for compatibility)
    freq_max_cm : float
        Maximum frequency in cm⁻¹
        
    Returns:
    --------
    freq_cm : np.array
        Frequency array in cm⁻¹
    spectrum : np.array
        Normalized power spectral density
    """
    # Constants
    boltz = 1.38064852E-23  # J/K
    lightspeed = 299792458.0  # m/s
    reduced_planck = 1.05457180013E-34  # J·s
    
    # Convert time from ps to seconds
    time_s = time * 1e-12
    dt = time_s[1] - time_s[0]
    total_time = time[-1] - time[0]  # Total simulation time in ps
    
    print(f"dt: {dt*1e12} ps")
    print(f"Total simulation time: {total_time} ps")
    
    # Calculate minimum resolvable frequency
    freq_resolution_hz = 1.0 / (total_time * 1e-12)  # in Hz
    freq_resolution_cm = freq_resolution_hz / (lightspeed * 100)  # in cm^-1
    print(f"Frequency resolution: {freq_resolution_cm:.1f} cm^-1")
    # Remove DC component
    energy_fluct = energy #- np.mean(energy)
    
    # Apply windowing
    window = np.hanning(len(energy_fluct))
    energy_windowed = energy_fluct * window
    
    # Calculate FFT for proper power spectral density
    fft_result = np.fft.fft(energy_windowed)
    
    # Calculate power spectral density
    # Factor of 2 accounts for negative frequencies (one-sided spectrum)
    # Normalize by sampling frequency and window power
    window_power = np.sum(window**2)
    power_spectrum = 2.0 * np.abs(fft_result)**2 / ((1/dt) * window_power)
    
    # Calculate frequencies for FFT
    freqs = np.fft.fftfreq(len(energy_windowed), dt)
    
    # Take only positive frequencies (one-sided spectrum)
    positive_mask = freqs > 0  # Exclude DC component (freq = 0)
    positive_freqs = freqs[positive_mask]
    power_spectrum = power_spectrum[positive_mask]
    
    # Convert frequency to cm⁻¹
    freq_cm = positive_freqs / (lightspeed * 100)
    
    # No frequency weighting - use raw power spectrum
    weighted_spectrum = power_spectrum
    
    # Filter to desired frequency range
    # Set minimum frequency based on simulation time resolution
    freq_min_cm = freq_resolution_cm
    
    if freq_max_cm is not None:
        freq_mask = (freq_cm >= freq_min_cm) & (freq_cm <= freq_max_cm)
        freq_cm_filtered = freq_cm[freq_mask]
        weighted_spectrum_filtered = weighted_spectrum[freq_mask]
    else:
        # Use full frequency range but exclude poorly resolved low frequencies
        freq_mask = freq_cm >= freq_min_cm
        freq_cm_filtered = freq_cm[freq_mask]
        weighted_spectrum_filtered = weighted_spectrum[freq_mask]
    
    print(f"Frequency range: {freq_cm_filtered.min():.1f} to {freq_cm_filtered.max():.1f} cm^-1")
    
    # Normalize by area under the curve
    if len(weighted_spectrum_filtered) > 1:
        area = np.trapezoid(weighted_spectrum_filtered, freq_cm_filtered)
        if area > 0:
            spectrum_normalized = weighted_spectrum_filtered / area
        else:
            spectrum_normalized = weighted_spectrum_filtered
    else:
        spectrum_normalized = weighted_spectrum_filtered
    
    return freq_cm_filtered, spectrum_normalized


def read_peaks_reference(filepath='peaks_reference.txt'):
    """
    Read the peaks reference file containing polariton and B-B frequencies.
    
    Returns:
    --------
    reference_data : dict
        Dictionary with coupling strengths as keys and frequency data as values
    """
    if not os.path.exists(filepath):
        print(f"Warning: {filepath} not found. Reference frequencies will not be plotted.")
        return None
    
    try:
        reference_data = {}
        
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Skip header lines starting with #
        data_lines = [line for line in lines if not line.strip().startswith('#') and line.strip()]
        
        for line in data_lines[1:]:  # Skip the column header line
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                try:
                    epsilon = float(parts[1])  # Coupling strength
                    peak2_freq = float(parts[2]) if parts[2] != '-1.0' else None  # Lower polariton
                    peak3_freq = float(parts[4]) if parts[4] != '-1.0' else None  # Upper polariton
                    peak4_freq = float(parts[6]) if parts[6] != '-1.0' else None  # B-B frequency
                    
                    reference_data[epsilon] = {
                        'lower_polariton': peak2_freq,
                        'upper_polariton': peak3_freq,
                        'bb_frequency': peak4_freq
                    }
                except (ValueError, IndexError):
                    continue
        
        print(f"Loaded reference frequencies for {len(reference_data)} coupling strengths")
        return reference_data
        
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None


def interpolate_reference_frequencies(reference_data, target_coupling):
    """
    Interpolate reference frequencies for a target coupling strength.
    
    Parameters:
    -----------
    reference_data : dict
        Reference frequency data
    target_coupling : float
        Target coupling strength to interpolate for
        
    Returns:
    --------
    frequencies : dict
        Interpolated frequencies, None if interpolation fails
    """
    if reference_data is None or len(reference_data) < 2:
        return None
    
    # Extract coupling strengths and frequencies
    couplings = np.array(sorted(reference_data.keys()))
    lower_pol_freqs = []
    upper_pol_freqs = []
    bb_freqs = []
    
    for coupling in couplings:
        data = reference_data[coupling]
        lower_pol_freqs.append(data['lower_polariton'])
        upper_pol_freqs.append(data['upper_polariton'])
        bb_freqs.append(data['bb_frequency'])
    
    result = {}
    
    # Interpolate each frequency type
    for freq_name, freq_list in [('lower_polariton', lower_pol_freqs), 
                                ('upper_polariton', upper_pol_freqs),
                                ('bb_frequency', bb_freqs)]:
        # Remove None values for interpolation
        valid_mask = np.array([f is not None for f in freq_list])
        if np.sum(valid_mask) >= 2:  # Need at least 2 points for interpolation
            valid_couplings = couplings[valid_mask]
            valid_freqs = np.array([f for f in freq_list if f is not None])
            
            # Check if target coupling is within range
            if valid_couplings.min() <= target_coupling <= valid_couplings.max():
                try:
                    if len(valid_freqs) >= 4:
                        # Use cubic interpolation if we have enough points
                        interp_func = interp1d(valid_couplings, valid_freqs, kind='cubic', 
                                             bounds_error=False, fill_value=np.nan)
                    else:
                        # Use linear interpolation for fewer points
                        interp_func = interp1d(valid_couplings, valid_freqs, kind='linear',
                                             bounds_error=False, fill_value=np.nan)
                    
                    interpolated_freq = interp_func(target_coupling)
                    if not np.isnan(interpolated_freq):
                        result[freq_name] = float(interpolated_freq)
                    else:
                        result[freq_name] = None
                except Exception as e:
                    print(f"Interpolation failed for {freq_name}: {e}")
                    result[freq_name] = None
            else:
                result[freq_name] = None
        else:
            result[freq_name] = None
    
    # Calculate derived frequencies BEFORE doubling (using original values)
    if result.get('lower_polariton') is not None and result.get('upper_polariton') is not None:
        result['rabi_frequency'] = result['upper_polariton'] - result['lower_polariton']  # ωUP - ωLP (not doubled)
        result['sum_frequency'] = result['lower_polariton'] + result['upper_polariton']    # ωLP + ωUP (not doubled)
    else:
        result['rabi_frequency'] = None
        result['sum_frequency'] = None
    
    # Double only the base reference frequencies for plotting
    for freq_type in ['lower_polariton', 'upper_polariton', 'bb_frequency']:
        if result.get(freq_type) is not None:
            result[freq_type] = 2.0 * result[freq_type]
    
    return result


def plot_reference_lines(ax, frequencies):
    """
    Plot reference frequency lines on the axis.
    
    Parameters:
    -----------
    ax : matplotlib axis
        The axis to plot on
    frequencies : dict
        Dictionary of interpolated frequencies
    """
    if frequencies is None:
        return
    
    # Define line styles and colors for different frequency types
    line_styles = {
        'lower_polariton': {'color': '#e74c3c', 'linestyle': ':', 'linewidth': 2, 'alpha': 0.8},
        'upper_polariton': {'color': '#e74c3c', 'linestyle': '--', 'linewidth': 2, 'alpha': 0.8},
        'bb_frequency': {'color': '#8e44ad', 'linestyle': '-', 'linewidth': 2, 'alpha': 0.8},
        'rabi_frequency': {'color': '#27ae60', 'linestyle': '-.', 'linewidth': 2, 'alpha': 0.8},
        'sum_frequency': {'color': '#f39c12', 'linestyle': '-', 'linewidth': 2, 'alpha': 0.8}
    }
    
    labels = {
        'lower_polariton': r'$2\omega_{LP}$ (2× Lower Polariton)',
        'upper_polariton': r'$2\omega_{UP}$ (2× Upper Polariton)',
        'bb_frequency': '2× B-B Frequency',
        'rabi_frequency': r'$\Omega_{Rabi} = \omega_{UP} - \omega_{LP}$',
        'sum_frequency': r'$\omega_{LP} + \omega_{UP}$'
    }
    
    reference_lines_added = []
    
    for freq_type, freq_value in frequencies.items():
        if freq_value is not None and 0 <= freq_value <= 5000:  # Within plot range
            style = line_styles.get(freq_type, {'color': 'gray', 'linestyle': '-', 'linewidth': 1})
            label = labels.get(freq_type, freq_type)
            
            ax.axvline(x=freq_value, label=label, **style)
            reference_lines_added.append((freq_type, freq_value))
    
    return reference_lines_added


def plot_energy_spectral_density_individual(base_dir='.', T=100.0):
    """
    Plot energy power spectral density using FFT, creating individual figures.
    Each figure shows one coupling strength with zero coupling as background reference.
    
    Parameters:
    -----------
    base_dir : str
        Base directory to search for coupling folders
    T : float
        Temperature in K (kept for compatibility)
    """
    # Read reference frequencies
    reference_data = read_peaks_reference()
    
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
    
    print(f"Temperature: {T} K")
    print(f"Found {len(valid_folders)} folders with averaged energy data:")
    for folder in valid_folders:
        coupling = extract_coupling_strength(os.path.basename(folder))
        print(f"  - {os.path.basename(folder)} (coupling = {coupling:.2e})")
    
    # Find zero coupling folder for background reference
    zero_coupling_folder = None
    zero_coupling_data = None
    
    for folder in valid_folders:
        coupling = extract_coupling_strength(os.path.basename(folder))
        if abs(coupling) < 1e-12:  # Essentially zero
            zero_coupling_folder = folder
            # Read zero coupling data once
            energy_file = os.path.join(folder, "prod_energy_tracker_averaged.txt")
            time_zero, data_zero, success_zero = read_averaged_energy_file(energy_file)
            if success_zero and data_zero is not None:
                molecular_total_zero = data_zero[:, 10] + data_zero[:, 2] + data_zero[:, 3] + data_zero[:, 4] + data_zero[:, 5]
                cavity_total_zero = data_zero[:, 11] + data_zero[:, 9]
                zero_coupling_data = {
                    'time': time_zero,
                    'molecular_total': molecular_total_zero,
                    'cavity_total': cavity_total_zero
                }
            break
    
    if zero_coupling_folder is None:
        print("Warning: No zero coupling folder found for background reference!")
    else:
        print(f"Using {os.path.basename(zero_coupling_folder)} as background reference")
    
    # Colors for energy components - in-cavity (coupled) spectra in black
    colors = {
        'molecular_total': '#000000',  # Black for in-cavity
        'cavity_total': '#000000'     # Black for in-cavity
    }
    
    # Colors for cavity-free (zero coupling) spectra in grey
    background_colors = {
        'molecular_total': '#808080',  # Grey for cavity-free
        'cavity_total': '#808080'     # Grey for cavity-free
    }
    
    labels = {
        'molecular_total': 'Molecular Total',
        'cavity_total': 'Cavity Total'
    }
    
    # Process each folder to create individual figures
    for folder in valid_folders:
        coupling = extract_coupling_strength(os.path.basename(folder))
        
        # Skip zero coupling since it will be used as background
        if abs(coupling) < 1e-12:
            continue
            
        print(f"\nCreating figure for coupling = {coupling:.2e}")
        
        # Interpolate reference frequencies for this coupling
        interpolated_freqs = interpolate_reference_frequencies(reference_data, coupling)
        
        # Create individual figure
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))  # Increased size to accommodate legend
        
        # Read current coupling data
        energy_file = os.path.join(folder, "prod_energy_tracker_averaged.txt")
        time, data, success = read_averaged_energy_file(energy_file)
        
        if not success or data is None or time is None:
            print(f"Error reading data for {os.path.basename(folder)}")
            plt.close(fig)
            continue
        
        # Extract energy components
        molecular_total = data[:, 10] + data[:, 2] + data[:, 3] + data[:, 4] + data[:, 5]
        cavity_total = data[:, 11] + data[:, 9]
        
        current_energies = {
            'molecular_total': molecular_total,
            'cavity_total': cavity_total
        }
        
        # Plot zero coupling in background first (if available)
        if zero_coupling_data is not None:
            for energy_type in ['molecular_total', 'cavity_total']:
                try:
                    freq_cm_zero, spectrum_zero = calculate_proper_energy_spectrum(
                        zero_coupling_data['time'], zero_coupling_data[energy_type], T)
                    if len(freq_cm_zero) > 0 and len(spectrum_zero) > 0:
                        # Molecular: bold solid line, Cavity: dashed line
                        linestyle = '-' if energy_type == 'molecular_total' else '--'
                        linewidth = 3.0 if energy_type == 'molecular_total' else 2.5  # Bold for molecular
                        ax.plot(freq_cm_zero, spectrum_zero, 
                               color=background_colors[energy_type], linewidth=linewidth, alpha=0.8,
                               linestyle=linestyle, label=f'{labels[energy_type]} (Zero Coupling)')
                except Exception as e:
                    print(f"Error calculating zero coupling spectrum for {energy_type}: {e}")
        
        # Store current spectra data for peak detection
        spectra_data = {}
        
        # Plot current coupling spectra in foreground
        for energy_type, energy_data in current_energies.items():
            try:
                freq_cm, spectrum = calculate_proper_energy_spectrum(time, energy_data, T)
                if len(freq_cm) > 0 and len(spectrum) > 0:
                    # Molecular: bold solid line, Cavity: dashed line
                    linestyle = '-' if energy_type == 'molecular_total' else '--'
                    linewidth = 3.0 if energy_type == 'molecular_total' else 2.0  # Bold for molecular
                    ax.plot(freq_cm, spectrum, label=f'{labels[energy_type]} (Coupling = {coupling:.2e})', 
                               color=colors[energy_type], linewidth=linewidth, alpha=0.9, linestyle=linestyle)
                    spectra_data[energy_type] = (freq_cm, spectrum)
            except Exception as e:
                print(f"Error calculating spectrum for {energy_type}: {e}")
                continue
        
        # Plot reference frequency lines
        reference_lines_added = plot_reference_lines(ax, interpolated_freqs)
        
        # Peak detection and annotation for current coupling only
        peak_text_lines = []
        
        for energy_type, (freq_cm, spectrum) in spectra_data.items():
            # Find peaks with minimum height threshold and minimum distance
            min_height = 0.00025  # Fixed minimum height threshold
            min_distance = max(1, len(spectrum) // 50)  # Minimum distance between peaks
            
            peaks, peak_properties = find_peaks(spectrum, 
                                              height=min_height,
                                              distance=min_distance,
                                              prominence=np.max(spectrum) * 0.05)
            
            if len(peaks) > 0:
                # Get peak heights and sort by height to find top 10
                peak_heights = spectrum[peaks]
                peak_freqs = freq_cm[peaks]
                
                # Sort peaks by height (descending)
                sorted_indices = np.argsort(peak_heights)[::-1]
                top_peaks = sorted_indices[:10]  # Top 10 peaks
                
                # Add peak annotations
                for j, peak_idx in enumerate(top_peaks):
                    peak_freq = peak_freqs[peak_idx]
                    peak_height = peak_heights[peak_idx]
                    
                    # Add vertical line at peak location
                    ax.axvline(x=peak_freq, color=colors[energy_type], 
                             linestyle='-.', alpha=0.5, linewidth=0.8)
                    
                    # Add text annotation
                    ax.annotate(f'{peak_freq:.0f}', 
                              xy=(peak_freq, peak_height),
                              xytext=(5, 5), textcoords='offset points',
                              fontsize=9, color=colors[energy_type],
                              bbox=dict(boxstyle='round,pad=0.2', 
                                      facecolor='white', alpha=0.8, edgecolor=colors[energy_type]),
                              ha='left')
                
                # Add peak information to text summary
                energy_label = labels[energy_type]
                peak_text_lines.append(f"{energy_label} Top 10 Peaks:")
                
                # Sort the top peaks by frequency (lowest to highest) for display
                top_peak_freqs = peak_freqs[top_peaks]
                freq_sorted_indices = np.argsort(top_peak_freqs)
                
                for j, sort_idx in enumerate(freq_sorted_indices):
                    peak_idx = top_peaks[sort_idx]
                    peak_freq = peak_freqs[peak_idx]
                    peak_text_lines.append(f"  {j+1}. {peak_freq:.0f} cm⁻¹")
        
        # Add reference frequency information to text summary
        if interpolated_freqs and any(v is not None for v in interpolated_freqs.values()):
            peak_text_lines.append("\nReference Frequencies:")
            for freq_type, freq_value in interpolated_freqs.items():
                if freq_value is not None:
                    freq_labels = {
                        'lower_polariton': '2× Lower Polariton',
                        'upper_polariton': '2× Upper Polariton', 
                        'bb_frequency': '2× B-B Frequency',
                        'rabi_frequency': 'Rabi Frequency',
                        'sum_frequency': 'ωLP + ωUP'
                    }
                    label = freq_labels.get(freq_type, freq_type)
                    peak_text_lines.append(f"  {label}: {freq_value:.0f} cm⁻¹")
        
        # Set labels and title
        ax.set_xlabel('Frequency (cm⁻¹)', fontsize=12)
        ax.set_ylabel('Power Spectral Density\n(normalized)', fontsize=12)
        ax.set_title(f'Energy Power Spectral Density - Coupling = {coupling:.2e}\n(T = {T} K, FFT-based)', 
                    fontsize=14, fontweight='bold')
        
        # Add legend with better positioning
        ax.legend(fontsize=9, loc='upper right', bbox_to_anchor=(1.0, 1.0))
        
        ax.grid(True, alpha=0.3)
        
        # Hardcoded x-limits
        ax.set_xlim(0, 5000)
        ax.set_ylim(bottom=0)  # Set y-axis to start at 0
        
        # Add peak summary text box
        if peak_text_lines:
            peak_text = '\n'.join(peak_text_lines)
            ax.text(0.02, 0.98, peak_text, transform=ax.transAxes,
                   fontsize=8, verticalalignment='top',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8))
        
        plt.tight_layout()
        
        # Save individual figure
        coupling_str = f"{coupling:.2e}".replace('e-0', 'e-').replace('e+0', 'e+')
        output_file = f"energy_spectral_density_coupling_{coupling_str}_T{T}K_FFT_with_refs.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Figure saved as: {output_file}")
        
        output_pdf = f"energy_spectral_density_coupling_{coupling_str}_T{T}K_FFT_with_refs.pdf"
        plt.savefig(output_pdf, bbox_inches='tight')
        print(f"Figure also saved as: {output_pdf}")
        
        plt.close(fig)  # Close figure to free memory
    
    print(f"\nCompleted generating individual figures for all coupling strengths!")


def main():
    """Main function."""
    import sys
    
    # Default temperature
    T = 100.0
    
    # Allow temperature to be specified as command line argument
    if len(sys.argv) > 1:
        try:
            T = float(sys.argv[1])
        except ValueError:
            print("Invalid temperature. Using default T = 300 K")
            T = 100.0
    
    print(f"Using temperature: {T} K")
    plot_energy_spectral_density_individual(T=T)


if __name__ == "__main__":
    main() 