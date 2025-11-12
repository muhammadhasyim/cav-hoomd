#!/usr/bin/env python3
"""
Test script for the new DualIndependentTemperatureFeedback controller.

This script demonstrates the dual independent temperature feedback controller
where cavity and molecular baths are controlled independently with different
target signals and parameters.

Example usage:
    python test_dual_controller.py --enable-dual-feedback \\
        --dual-cavity-method harmonic_equipartition \\
        --dual-molecular-method lj_coulombic_kinetic \\
        --dual-cavity-target 80.0 \\
        --dual-molecular-target 120.0 \\
        --dual-cavity-time-constant 5.0 \\
        --dual-molecular-time-constant 10.0
"""

import subprocess
import sys
from pathlib import Path

def run_dual_controller_test():
    """Run a test simulation with the dual independent controller."""
    
    # Test parameters
    test_args = [
        'python', '18_unified_cavity_dynamics.py',
        '--coupling', '1e-3',
        '--frequency', '1560.0',
        '--temperature', '100.0',
        '--runtime', '50.0',  # Short run for testing
        '--replicas', '1',
        '--molecular-bath', 'bussi',
        '--cavity-bath', 'langevin',
        '--enable-energy-tracker',
        '--enable-temp-tracker',
        
        # Enable dual independent feedback controller
        '--enable-dual-feedback',
        '--dual-cavity-method', 'harmonic_equipartition',
        '--dual-molecular-method', 'lj_coulombic_kinetic',
        '--dual-cavity-target', '80.0',
        '--dual-molecular-target', '120.0',
        '--dual-cavity-time-constant', '5.0',
        '--dual-molecular-time-constant', '10.0',
        '--dual-turn-on-time', '5.0',  # Start control after 5 ps
        '--dual-update-interval', '0.1',
        
        # Output settings
        '--console-output-period', '1.0',
        '--output-dir', 'test_dual_controller',
    ]
    
    print("=" * 60)
    print("TESTING DUAL INDEPENDENT TEMPERATURE FEEDBACK CONTROLLER")
    print("=" * 60)
    print()
    print("Configuration:")
    print("  Cavity bath: harmonic_equipartition → 80.0K (τ=5.0ps)")
    print("  Molecular bath: lj_coulombic_kinetic → 120.0K (τ=10.0ps)")
    print("  Runtime: 50 ps")
    print("  Turn-on time: 5 ps")
    print()
    print("This test demonstrates independent control of cavity and molecular")
    print("bath temperatures using different signal types and targets.")
    print()
    
    try:
        # Run the simulation
        print("Running simulation...")
        result = subprocess.run(test_args, check=True, capture_output=True, text=True)
        
        print("Simulation completed successfully!")
        print()
        print("Output files created:")
        print("  - test_dual_controller/dual_feedback_replica_0.csv")
        print("  - test_dual_controller/temperature_tracker_replica_0.csv")
        print("  - test_dual_controller/energy_tracker_replica_0.csv")
        print()
        print("To analyze the results:")
        print("  1. Check dual_feedback_replica_0.csv for controller behavior")
        print("  2. Check temperature_tracker_replica_0.csv for all temperature signals")
        print("  3. Plot the data to see independent temperature evolution")
        
        return True
        
    except subprocess.CalledProcessError as e:
        print("ERROR: Simulation failed!")
        print(f"Error: {e}")
        if e.stdout:
            print("STDOUT:")
            print(e.stdout)
        if e.stderr:
            print("STDERR:")
            print(e.stderr)
        return False
    except Exception as e:
        print(f"ERROR: Unexpected error: {e}")
        return False

def create_simple_analysis_script():
    """Create a simple script to analyze the dual controller output."""
    
    analysis_script = '''#!/usr/bin/env python3
"""
Simple analysis script for dual controller output.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def plot_dual_controller_results():
    """Plot the results from dual controller test."""
    
    # Read dual controller data
    dual_file = Path("test_dual_controller/dual_feedback_replica_0.csv")
    if not dual_file.exists():
        print(f"Error: {dual_file} not found")
        return
    
    df_dual = pd.read_csv(dual_file)
    
    # Create plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Dual Independent Temperature Controller Results', fontsize=16)
    
    # Plot 1: Temperature signals
    ax1 = axes[0, 0]
    ax1.plot(df_dual['time_ps'], df_dual['cavity_temp'], 'b-', label='Cavity signal')
    ax1.plot(df_dual['time_ps'], df_dual['molecular_temp'], 'r-', label='Molecular signal')
    ax1.axhline(y=80.0, color='b', linestyle='--', alpha=0.7, label='Cavity target (80K)')
    ax1.axhline(y=120.0, color='r', linestyle='--', alpha=0.7, label='Molecular target (120K)')
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Temperature (K)')
    ax1.set_title('Temperature Signals')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Bath temperatures (controller outputs)
    ax2 = axes[0, 1]
    ax2.plot(df_dual['time_ps'], df_dual['cavity_output'], 'b-', label='Cavity bath')
    ax2.plot(df_dual['time_ps'], df_dual['molecular_output'], 'r-', label='Molecular bath')
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Bath Temperature (K)')
    ax2.set_title('Controller Outputs (Bath Temperatures)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Error signals
    ax3 = axes[1, 0]
    ax3.plot(df_dual['time_ps'], df_dual['cavity_error'], 'b-', label='Cavity error')
    ax3.plot(df_dual['time_ps'], df_dual['molecular_error'], 'r-', label='Molecular error')
    ax3.axhline(y=0, color='k', linestyle='-', alpha=0.5)
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('Error (K)')
    ax3.set_title('Control Errors')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Effective temperatures
    ax4 = axes[1, 1]
    ax4.plot(df_dual['time_ps'], df_dual['cavity_effective'], 'b-', label='Cavity effective')
    ax4.plot(df_dual['time_ps'], df_dual['molecular_effective'], 'r-', label='Molecular effective')
    ax4.axhline(y=80.0, color='b', linestyle='--', alpha=0.7, label='Cavity target')
    ax4.axhline(y=120.0, color='r', linestyle='--', alpha=0.7, label='Molecular target')
    ax4.set_xlabel('Time (ps)')
    ax4.set_ylabel('Effective Temperature (K)')
    ax4.set_title('Effective Temperatures')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('test_dual_controller_results.png', dpi=300, bbox_inches='tight')
    print("Plot saved as: test_dual_controller_results.png")
    plt.show()

if __name__ == "__main__":
    plot_dual_controller_results()
'''
    
    with open('/home/mh7373/GitRepos/cav-hoomd/analyze_dual_controller.py', 'w') as f:
        f.write(analysis_script)
    
    print("Created analysis script: analyze_dual_controller.py")

def main():
    """Main test function."""
    print("Dual Independent Temperature Controller Test")
    print("=" * 50)
    
    # Check if we have the required input file
    input_file = Path("molecular-0.gsd")
    if not input_file.exists():
        print("WARNING: molecular-0.gsd not found in current directory")
        print("   The simulation may fail without input coordinates")
        print("   Please ensure you have a valid input GSD file")
        print()
    
    # Create analysis script
    create_simple_analysis_script()
    print()
    
    # Run the test
    success = run_dual_controller_test()
    
    if success:
        print()
        print("Test completed successfully!")
        print()
        print("Next steps:")
        print("  1. Run: python analyze_dual_controller.py")
        print("  2. Examine the generated plots to verify independent control")
        print("  3. Check CSV files for detailed time series data")
        
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
