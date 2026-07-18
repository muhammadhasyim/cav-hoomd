#!/usr/bin/env python3
"""
Debug script to analyze energy FT calculation step by step.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def debug_energy_ft():
    """Debug the energy FT calculation."""
    
    # Read one of the energy files for debugging
    energy_file = "cavity_coupling_0epos00/prod_energy_tracker_averaged.txt"
    
    if not os.path.exists(energy_file):
        print(f"File {energy_file} not found!")
        return
    
    # Read data
    data = np.loadtxt(energy_file, comments='#')
    time = data[:, 0]  # time in ps
    
    # Calculate molecular total energy
    molecular_total = data[:, 10] + data[:, 2] + data[:, 3] + data[:, 4] + data[:, 5]
    
    print(f"Time range: {time.min():.6f} to {time.max():.6f} ps")
    print(f"Number of points: {len(time)}")
    print(f"Time step: {time[1] - time[0]:.6f} ps")
    print(f"Energy range: {molecular_total.min():.6f} to {molecular_total.max():.6f}")
    print(f"Energy mean: {molecular_total.mean():.6f}")
    print(f"Energy std: {molecular_total.std():.6f}")
    
    # Create debugging plots
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Plot 1: Time series
    axes[0,0].plot(time, molecular_total)
    axes[0,0].set_xlabel('Time (ps)')
    axes[0,0].set_ylabel('Molecular Total Energy')
    axes[0,0].set_title('Energy vs Time')
    axes[0,0].grid(True)
    
    # Plot 2: Energy fluctuations (remove mean)
    energy_fluct = molecular_total - np.mean(molecular_total)
    axes[0,1].plot(time, energy_fluct)
    axes[0,1].set_xlabel('Time (ps)')
    axes[0,1].set_ylabel('Energy Fluctuations')
    axes[0,1].set_title('Energy Fluctuations (mean removed)')
    axes[0,1].grid(True)
    
    # Calculate FFT
    time_s = time * 1e-12  # Convert to seconds
    dt = time_s[1] - time_s[0]
    
    # Apply window
    window = np.hanning(len(energy_fluct))
    energy_windowed = energy_fluct * window
    
    # FFT
    fft_energy = np.fft.fft(energy_windowed)
    freqs_hz = np.fft.fftfreq(len(time), dt)
    
    # Positive frequencies only
    positive_freqs = freqs_hz[:len(freqs_hz)//2]
    fft_positive = fft_energy[:len(fft_energy)//2]
    
    # Convert to cm⁻¹
    c_light = 2.998e8  # m/s
    freq_cm = positive_freqs / (c_light * 100)
    
    print(f"Frequency range: {freq_cm.min():.2f} to {freq_cm.max():.2f} cm⁻¹")
    print(f"Frequency resolution: {freq_cm[1] - freq_cm[0]:.2f} cm⁻¹")
    
    # Plot 3: Power spectrum (no weighting)
    power_spectrum = np.abs(fft_positive)**2
    
    # Filter to 0-3500 cm⁻¹ range
    freq_mask = (freq_cm >= 0) & (freq_cm <= 3500)
    freq_plot = freq_cm[freq_mask]
    power_plot = power_spectrum[freq_mask]
    
    axes[0,2].semilogy(freq_plot, power_plot)
    axes[0,2].set_xlabel('Frequency (cm⁻¹)')
    axes[0,2].set_ylabel('Power Spectrum')
    axes[0,2].set_title('Unweighted Power Spectrum')
    axes[0,2].set_xlim(0, 3500)
    axes[0,2].grid(True)
    
    # Plot 4: ω² weighting
    omega = 2 * np.pi * positive_freqs
    omega_squared = omega**2
    omega_squared[0] = 0  # Avoid division by zero
    
    omega_plot = omega_squared[freq_mask]
    axes[1,0].plot(freq_plot, omega_plot)
    axes[1,0].set_xlabel('Frequency (cm⁻¹)')
    axes[1,0].set_ylabel('ω²')
    axes[1,0].set_title('ω² Weighting Factor')
    axes[1,0].set_xlim(0, 3500)
    axes[1,0].grid(True)
    
    # Plot 5: ω² weighted spectrum
    weighted_spectrum = power_spectrum * omega_squared
    weighted_plot = weighted_spectrum[freq_mask]
    
    axes[1,1].semilogy(freq_plot, weighted_plot)
    axes[1,1].set_xlabel('Frequency (cm⁻¹)')
    axes[1,1].set_ylabel('FT(E(t))×ω²')
    axes[1,1].set_title('ω² Weighted Spectrum')
    axes[1,1].set_xlim(0, 3500)
    axes[1,1].grid(True)
    
    # Plot 6: Normalized spectrum
    area = np.trapezoid(weighted_plot, freq_plot)
    if area > 0:
        normalized_spectrum = weighted_plot / area
    else:
        normalized_spectrum = weighted_plot
    
    axes[1,2].semilogy(freq_plot, normalized_spectrum)
    axes[1,2].set_xlabel('Frequency (cm⁻¹)')
    axes[1,2].set_ylabel('Normalized FT(E(t))×ω²')
    axes[1,2].set_title('Final Normalized Spectrum')
    axes[1,2].set_xlim(0, 3500)
    axes[1,2].grid(True)
    
    plt.tight_layout()
    plt.savefig('debug_energy_ft.png', dpi=300, bbox_inches='tight')
    print("\nDebug plots saved as: debug_energy_ft.png")
    
    # Print some statistics
    print(f"\nSpectrum statistics:")
    print(f"Power spectrum max: {power_plot.max():.2e}")
    print(f"Weighted spectrum max: {weighted_plot.max():.2e}")
    print(f"Area under weighted curve: {area:.2e}")
    
    # Check where most power is
    max_power_idx = np.argmax(power_plot)
    max_weighted_idx = np.argmax(weighted_plot)
    print(f"Peak unweighted power at: {freq_plot[max_power_idx]:.2f} cm⁻¹")
    print(f"Peak weighted power at: {freq_plot[max_weighted_idx]:.2f} cm⁻¹")
    
    plt.show()

if __name__ == "__main__":
    debug_energy_ft() 