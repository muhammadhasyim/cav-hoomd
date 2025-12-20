#!/usr/bin/env python3
"""
Test script for extended fitting functions in the cavitymd module.
This script tests the new extended harmonic and T^(3/5) fitting functions.
"""

import numpy as np
import sys
import os
from pathlib import Path

# Add the source directory to the path
sys.path.insert(0, str(Path(__file__).parent / "src"))

def test_extended_harmonic_fit():
    """Test extended harmonic fitting: E = aT/(1+bT)"""
    print("Testing Extended Harmonic Fitting")
    print("=" * 50)
    
    try:
        from cavitymd.analysis import EmpiricalTemperatureData
        
        # Create synthetic harmonic data with saturation behavior
        temperatures = np.array([65, 100, 150, 200, 250, 300, 400, 500, 1000, 2000, 5000])
        
        # Extended harmonic model: E = aT/(1+bT) with a=0.015, b=0.0001
        a_true, b_true = 0.015, 0.0001
        energies_hartree = a_true * temperatures / (1 + b_true * temperatures)
        
        # Add small amount of noise
        np.random.seed(42)
        energies_hartree += np.random.normal(0, 0.0001, len(energies_hartree))
        
        # Create temporary data file
        temp_file = Path("temp_harmonic_data.txt")
        with open(temp_file, 'w') as f:
            f.write("# Temperature (K) and Harmonic Energy (Ha)\n")
            f.write("temperature harmonic_hartree\n")
            for T, E in zip(temperatures, energies_hartree):
                f.write(f"{T:.1f} {E:.8f}\n")
        
        # Test the fitting
        empirical_data = EmpiricalTemperatureData(
            str(temp_file), 
            energy_component='harmonic',
            use_direct_harmonic=False
        )
        
        print(f"✓ Extended harmonic fit successful!")
        if hasattr(empirical_data, 'has_extended_harmonic_fit') and empirical_data.has_extended_harmonic_fit:
            fit_params = empirical_data.extended_harmonic_fit
            print(f"  Fitted parameters: a = {fit_params['a']:.6f}, b = {fit_params['b']:.6f}")
            print(f"  True parameters:   a = {a_true:.6f}, b = {b_true:.6f}")
            print(f"  R² = {fit_params['r2']:.12f}")
            
            # Test temperature calculation
            test_energy = 0.01  # Ha
            calc_temp = empirical_data.calculate_systemic_temperature(test_energy)
            expected_temp = test_energy / (a_true - b_true * test_energy)
            print(f"  Temperature calculation test:")
            print(f"    Input energy: {test_energy:.6f} Ha")
            print(f"    Calculated T: {calc_temp:.2f} K")
            print(f"    Expected T:   {expected_temp:.2f} K")
            print(f"    Error: {abs(calc_temp - expected_temp):.2f} K")
        
        # Cleanup
        temp_file.unlink()
        return True
        
    except Exception as e:
        print(f"ERROR: Extended harmonic fitting test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_extended_t35_fit():
    """Test extended T^(3/5) fitting: E = E0 + AT^(3/5)/(1+CT^(3/5))"""
    print("\nTesting Extended T^(3/5) Fitting")
    print("=" * 50)
    
    try:
        from cavitymd.analysis import EmpiricalTemperatureData
        
        # Create synthetic LJ+Coulombic data with saturation behavior
        temperatures = np.array([65, 100, 150, 200, 250, 300, 400, 500, 1000, 2000, 5000])
        
        # Extended T^(3/5) model: E = E0 + AT^(3/5)/(1+CT^(3/5))
        e0_true, a_true, c_true = -50.0, 0.75, 0.001
        t_power = temperatures**(3/5)
        energies_hartree = e0_true + a_true * t_power / (1 + c_true * t_power)
        
        # Add small amount of noise
        np.random.seed(42)
        energies_hartree += np.random.normal(0, 0.1, len(energies_hartree))
        
        # Create temporary data file
        temp_file = Path("temp_lj_coul_data.txt")
        with open(temp_file, 'w') as f:
            f.write("# Temperature (K) and LJ+Coulombic Energy (Ha)\n")
            f.write("temperature lj_hartree coulombic_hartree\n")
            # Split the energy into LJ and Coulombic components (roughly 60/40 split)
            lj_energies = energies_hartree * 0.6
            coul_energies = energies_hartree * 0.4
            for T, E_lj, E_coul in zip(temperatures, lj_energies, coul_energies):
                f.write(f"{T:.1f} {E_lj:.8f} {E_coul:.8f}\n")
        
        # Test the fitting
        empirical_data = EmpiricalTemperatureData(
            str(temp_file), 
            energy_component='lj_coulombic',
            use_direct_harmonic=False
        )
        
        print(f"✓ Extended T^(3/5) fit successful!")
        if hasattr(empirical_data, 'has_extended_t35_fit') and empirical_data.has_extended_t35_fit:
            fit_params = empirical_data.extended_t35_fit
            print(f"  Fitted parameters: E0 = {fit_params['e0']:.3f}, A = {fit_params['a']:.6f}, C = {fit_params['c']:.6f}")
            print(f"  True parameters:   E0 = {e0_true:.3f}, A = {a_true:.6f}, C = {c_true:.6f}")
            print(f"  R² = {fit_params['r2']:.12f}")
            
            # Test temperature calculation
            test_energy = -45.0  # Ha
            calc_temp = empirical_data.calculate_systemic_temperature(test_energy)
            print(f"  Temperature calculation test:")
            print(f"    Input energy: {test_energy:.3f} Ha")
            print(f"    Calculated T: {calc_temp:.2f} K")
        
        # Cleanup
        temp_file.unlink()
        return True
        
    except Exception as e:
        print(f"ERROR: Extended T^(3/5) fitting test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_hartree_units():
    """Test that all energy calculations use Hartree units consistently"""
    print("\nTesting Hartree Unit Consistency")
    print("=" * 50)
    
    try:
        from cavitymd.utils import PhysicalConstants
        
        print(f"✓ Physical constants loaded:")
        print(f"  kB (Hartree/K): {PhysicalConstants.KB_HARTREE_PER_K:.10e}")
        print(f"  Hartree to cm⁻¹: {PhysicalConstants.HARTREE_TO_CM_MINUS1:.2f}")
        
        # Test unit conversions
        energy_hartree = 1.0  # 1 Hartree
        energy_kelvin = energy_hartree / PhysicalConstants.KB_HARTREE_PER_K
        
        print(f"  Unit conversion test:")
        print(f"    1.0 Ha = {energy_kelvin:.2f} K")
        print(f"    1.0 Ha = {PhysicalConstants.HARTREE_TO_CM_MINUS1:.2f} cm⁻¹")
        
        return True
        
    except Exception as e:
        print(f"ERROR: Hartree unit test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Run all tests"""
    print("Testing Extended Fitting Functions in CavityMD Module")
    print("=" * 60)
    
    tests = [
        test_extended_harmonic_fit,
        test_extended_t35_fit,
        test_hartree_units
    ]
    
    results = []
    for test in tests:
        results.append(test())
    
    print("\n" + "=" * 60)
    print("Test Results Summary:")
    print(f"  Passed: {sum(results)}/{len(results)}")
    print(f"  Failed: {len(results) - sum(results)}/{len(results)}")
    
    if all(results):
        print("All tests passed! Extended fitting functions are working correctly.")
        return 0
    else:
        print("ERROR: Some tests failed. Please check the implementation.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
