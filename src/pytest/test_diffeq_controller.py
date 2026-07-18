#!/usr/bin/env python3
"""
Test script for DiffEqController

This script demonstrates the differential equation temperature controller:
dT_bath(t)/dt = -(T_bath(t) - T_signal(t)) / τ

The controller makes the bath temperature exponentially approach the signal temperature
with a time constant τ.
"""

import subprocess
import sys
import os

def run_diffeq_controller_test():
    """Run a test of the DiffEqController with composite coupling."""
    
    print("Testing DiffEqController with Composite Coupling")
    print("=" * 60)
    print()
    print("📋 Test Configuration:")
    print("  • Controller: DiffEqController")
    print("  • Signal: lj_coulombic_bath (LJ+Coulombic fictive temperature)")
    print("  • Time constant: 2.0 ps")
    print("  • Differential equation: dT_bath/dt = -(T_bath - T_signal)/τ")
    print("  • Coupling: Composite (sinusoid + square wave)")
    print("  • Runtime: 100 ps")
    print()
    
    # Command to run the DiffEq controller test
    cmd = [
        "python3", "18_unified_cavity_dynamics.py",
        
        # Basic simulation parameters
        "--molecular-bath", "bussi",
        "--cavity-bath", "langevin", 
        "--coupling-type", "composite",
        "--coupling", "5e-4",
        
        # Composite coupling parameters
        "--composite-sinusoid-amplitude", "2e-4",
        "--composite-sinusoid-period", "2.0",
        "--composite-sinusoid-phase", "0.0",
        "--composite-square-amplitude", "3e-4", 
        "--composite-square-period", "20.0",
        "--composite-square-duty-cycle", "0.1",
        "--composite-square-start-time", "10.0",
        "--composite-square-stop-time", "80.0",
        "--composite-sinusoid-start-time", "80.0",
        
        # DiffEq controller parameters
        "--enable-diffeq-controller",
        "--diffeq-temperature-method", "lj_coulombic_bath",
        "--diffeq-time-constant", "2.0",
        "--diffeq-turn-on-time", "5.0",
        "--diffeq-update-interval", "0.01",
        "--diffeq-apply-to", "molecular",
        "--diffeq-T-min", "1.0",
        
        # System parameters
        "--temperature", "300.0",
        "--frequency", "1560.0",
        "--molecular-tau", "1.0",
        "--cavity-tau", "1.0",
        "--runtime", "100.0",
        
        # Tracking and analysis
        "--enable-energy-tracker",
        "--enable-temp-tracker",
        "--temp-tracker-empirical-data-file", "/home/mh7373/GitRepos/cav-hoomd/final_nodiss_cavitymd/potential_energy_components_vs_temperature.txt",
        
        # Hardware and output
        "--device", "GPU",
        "--seed", "42",
        "--replicas", "0",
        "--energy-output-period-ps", "0.01",
        "--console-output-period-ps", "1.0",
        "--error-tolerance", "5.0",
        "--initial-fraction", "1e-6",
        "--input-gsd", "final_nodiss_periodic/init-0.gsd",
        "--time-constant-ps", "10.0",
        "--enable-dynamic-coupling-detection",
        "--coupling-change-threshold", "1e-5"
    ]
    
    print("Running DiffEq controller test...")
    print("Command:", " ".join(cmd))
    print()
    
    try:
        # Run the command
        result = subprocess.run(cmd, capture_output=False, text=True, cwd=os.getcwd())
        
        if result.returncode == 0:
            print("DiffEq controller test completed successfully!")
        else:
            print(f"ERROR: DiffEq controller test failed with return code {result.returncode}")
            return False
            
    except KeyboardInterrupt:
        print("\nTest interrupted by user")
        return False
    except Exception as e:
        print(f"ERROR: Error running DiffEq controller test: {e}")
        return False
    
    return True

if __name__ == "__main__":
    success = run_diffeq_controller_test()
    sys.exit(0 if success else 1)
