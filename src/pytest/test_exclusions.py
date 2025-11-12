#!/usr/bin/env python3
"""
Test script to demonstrate data exclusion and time shifting functionality
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def create_test_autocorr_data():
    """Create synthetic autocorrelation data for testing."""
    # Create test data
    time_ps = np.linspace(0, 10, 1000)  # 0 to 10 ps
    
    # Synthetic autocorrelation: exponential decay with oscillations
    autocorr = np.exp(-time_ps/2.0) * (1 + 0.2 * np.cos(2*np.pi*time_ps/0.5))
    
    # Add some noise to the beginning and end (to simulate artifacts)
    noise_start = np.random.normal(0, 0.1, 50)
    noise_end = np.random.normal(0, 0.05, 50)
    
    autocorr[:50] += noise_start
    autocorr[-50:] += noise_end
    
    return time_ps, autocorr

def save_test_file(time_ps, autocorr, filename):
    """Save test autocorrelation data in the expected format."""
    # Convert time to ps and create timestep column
    timestep = np.arange(len(time_ps))
    
    header = "# Test autocorrelation data\n# timestep lag_time(ps) C(t)"
    data = np.column_stack((timestep, time_ps, autocorr))
    
    np.savetxt(filename, data, header=header, fmt='%d %.6f %.6f')
    print(f"Created test file: {filename}")

def main():
    """Create test files and demonstrate usage."""
    print("Creating test autocorrelation files...")
    
    # Create test directory
    test_dir = "test_exclusions"
    os.makedirs(test_dir, exist_ok=True)
    
    # Create multiple test files (simulating replicas and references)
    for replica in range(2):
        for ref in range(2):
            time_ps, autocorr = create_test_autocorr_data()
            
            # Add slight variations between replicas/references
            autocorr *= (1 + 0.05 * np.random.normal())
            
            filename = f"{test_dir}/prod-{replica}_dipole_autocorr_{ref}.txt"
            save_test_file(time_ps, autocorr, filename)
    
    print(f"\nTest files created in {test_dir}/")
    print("\nNow you can test the exclusion functionality with:")
    print(f"python examples/autocorr_to_ir_individual.py --directory {test_dir}")
    print(f"python examples/autocorr_to_ir_individual.py --directory {test_dir} --exclude-initial 0.1 --exclude-final 0.05")
    
    print("\nThe second command will:")
    print("1. Remove the first 10% of data (noisy beginning)")
    print("2. Remove the last 5% of data (noisy end)")
    print("3. Shift the time array to start from t=0")
    print("4. Calculate IR spectra from the cleaned, shifted data")

if __name__ == "__main__":
    main() 