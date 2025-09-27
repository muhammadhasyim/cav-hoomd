#!/usr/bin/env python3
"""
Test script to verify IR spectrum calculation matches reference code.
"""

import numpy as np
from scipy import fftpack

# Physical constants (SI units) - same as in our main code
BOLTZ = 1.38064852e-23  # Boltzmann constant (J/K)
LIGHTSPEED = 299792458.0  # Speed of light (m/s)
REDUCED_PLANCK = 1.05457180013e-34  # Reduced Planck constant (J·s)

# Test parameters
T = 300.0  # Temperature in Kelvin

def test_ir_calculation():
    """Test our IR calculation against the reference method."""
    print("Testing IR spectrum calculation method")
    print("=" * 50)
    
    # Create a simple test autocorrelation function
    time_fs = np.linspace(0, 1000, 101)  # 0 to 1000 fs, 101 points
    autocorr = np.exp(-time_fs / 100.0) * np.cos(2*np.pi*time_fs/50.0)  # Decaying oscillation
    
    print(f"Test autocorrelation: {len(time_fs)} points from {time_fs[0]:.1f} to {time_fs[-1]:.1f} fs")
    print(f"Temperature: {T:.1f} K")
    
    # Calculate using our method
    timestep = (time_fs[1] - time_fs[0]) * 1e-15  # Convert to seconds
    
    # Calculate FFT of autocorrelation function
    lineshape = fftpack.dct(autocorr, type=1)[1:]  # Skip DC component
    
    # Calculate frequency arrays
    lineshape_frequencies = np.linspace(0, 0.5/timestep, len(autocorr))[1:]
    lineshape_frequencies_wn = lineshape_frequencies / (100.0 * LIGHTSPEED)  # Convert to cm^-1
    
    # Our new method (following reference code)
    field_description = lineshape_frequencies * (1.0 - np.exp(-REDUCED_PLANCK * lineshape_frequencies / (BOLTZ * T)))
    quantum_correction = lineshape_frequencies / (1.0 - np.exp(-REDUCED_PLANCK * lineshape_frequencies / (BOLTZ * T)))
    spectra = lineshape * field_description
    spectra_qm = spectra * quantum_correction
    
    # Old method (for comparison)
    omega_squared_old = (2.0 * np.pi * lineshape_frequencies)**2
    quantum_correction_old = lineshape_frequencies / (1.0 - np.exp(-REDUCED_PLANCK * lineshape_frequencies / (BOLTZ * T)))
    ir_spectrum_old = lineshape * omega_squared_old * quantum_correction_old
    
    # Mathematical verification: field_description * quantum_correction should equal omega^2
    omega_squared_new = field_description * quantum_correction
    omega_squared_theoretical = lineshape_frequencies**2
    
    print(f"\nFrequency range: {lineshape_frequencies_wn[0]:.1f} to {lineshape_frequencies_wn[-1]:.1f} cm^-1")
    print(f"Number of frequency points: {len(lineshape_frequencies_wn)}")
    
    # Compare some values
    print(f"\nComparison at first few frequency points:")
    print(f"{'Freq (cm-1)':<12} {'Lineshape':<12} {'Field_desc':<12} {'Quantum_corr':<12} {'Spectra':<12} {'Spectra_qm':<12}")
    print("-" * 78)
    for i in range(min(5, len(lineshape_frequencies_wn))):
        print(f"{lineshape_frequencies_wn[i]:<12.1f} {lineshape[i]:<12.3e} {field_description[i]:<12.3e} "
              f"{quantum_correction[i]:<12.3e} {spectra[i]:<12.3e} {spectra_qm[i]:<12.3e}")
    
    # Verify mathematical relationship
    ratio_check = omega_squared_new / omega_squared_theoretical
    print(f"\nMathematical verification:")
    print(f"field_description * quantum_correction should equal omega^2")
    print(f"Ratio check (should be ~1.0): min={ratio_check.min():.6f}, max={ratio_check.max():.6f}")
    
    # Compare new vs old method
    ratio_new_old = spectra_qm / ir_spectrum_old
    print(f"\nComparison: new method / old method")
    print(f"Ratio range: {ratio_new_old.min():.3e} to {ratio_new_old.max():.3e}")
    print(f"The old method used (2π*ω)^2, new method uses ω^2")
    print(f"Expected ratio should be ~1/(2π)^2 = {1/(2*np.pi)**2:.6f}")
    
    expected_ratio = 1.0 / (2*np.pi)**2
    print(f"Actual average ratio: {np.mean(ratio_new_old):.6f}")
    print(f"Matches expected ratio: {np.allclose(ratio_new_old, expected_ratio, rtol=1e-10)}")

if __name__ == "__main__":
    test_ir_calculation() 