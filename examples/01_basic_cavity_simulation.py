#!/usr/bin/env python3
"""
Example 1: Basic Cavity MD Simulation

This example demonstrates the most basic cavity MD simulation setup using
the existing infrastructure and init-0.gsd system file.

Key features:
- Uses the existing molecular system from init-0.gsd  
- Simple cavity-molecule coupling
- Short simulation to demonstrate basic functionality
- Uses the CavityMDSimulation framework for proper setup
"""

import sys
from pathlib import Path

# Import from the modular hoomd.cavitymd package
from hoomd.cavitymd import CavityMDSimulation

def run_basic_cavity_simulation():
    """Run a basic cavity MD simulation using the existing CavityMDSimulation class."""
    
    print("="*60)
    print("BASIC CAVITY MD SIMULATION")
    print("="*60)
    
    # Create output directory
    output_dir = Path("basic_cavity_example")
    output_dir.mkdir(exist_ok=True)
    
    print(f"Output directory: {output_dir}")
    print(f"Using system file: ../init-0.gsd")
    print(f"System: 500 particles (O and N types)")
    
    # Set up basic simulation parameters
    coupling_strength = 1e-3    # Weak coupling for demonstration
    runtime_ps = 10.0          # Short 10 ps simulation
    replica = 1                # Single replica
    frame = -1                 # Use last frame from gsd file
    
    print(f"Coupling strength: {coupling_strength}")
    print(f"Runtime: {runtime_ps} ps")
    print(f"Replica: {replica}")
    
    # Create and run the simulation
    print("\nSetting up cavity MD simulation...")
    
    sim = CavityMDSimulation(
        job_dir=str(output_dir),
        replica=replica,
        freq=1560.0,                    # Cavity frequency (cm^-1)
        couplstr=coupling_strength,     # Coupling strength  
        incavity=True,                  # Enable cavity coupling
        runtime_ps=runtime_ps,          # Simulation time
        input_gsd='../init-0.gsd',      # Use existing system file
        frame=frame,                    # Use last frame
        name='basic_prod',              # Output file prefix
        error_tolerance=0.01,           # Small error tolerance for adaptive timestep
        temperature=100.0,              # Temperature (K)
        molecular_thermostat='bussi',   # BussiReservoir thermostat for molecules
        cavity_thermostat='langevin',   # Langevin thermostat for cavity
        cavity_damping_factor=1.0,      # Moderate damping
        use_brownian_overdamped=True,   # Use overdamped dynamics
        add_cavity_particle=True,       # Add cavity particle
        finite_q=True,                  # Use finite Q model
        molecular_thermostat_tau=5.0,   # Molecular thermostat time constant (ps)
        cavity_thermostat_tau=5.0,      # Cavity thermostat time constant (ps)
        log_to_file=False,              # Don't log to file for this basic example
        log_to_console=True,            # Log to console
        log_level='INFO',               # Information level logging
        enable_fkt=False,               # Disable F(k,t) for simplicity
        max_energy_output_time_ps=None, # No energy output limit
        enable_energy_tracking=True,    # Track energy components
        dt_fs=None,                     # Use adaptive timestep
        device='CPU',                   # Use CPU
        gpu_id=0                        # GPU ID (not used with CPU)
    )
    
    print("Running cavity MD simulation...")
    print("-" * 40)
    
    # Run the simulation
    exit_code = sim.run()
    
    print("-" * 40)
    print(f"Simulation completed with exit code: {exit_code}")
    
    if exit_code == 0:
        print("✅ Basic cavity MD simulation completed successfully!")
        print(f"\nOutput files are in: {output_dir}/")
        print("Files created:")
        
        # List output files
        for file in output_dir.glob("*"):
            if file.is_file():
                size_kb = file.stat().st_size / 1024
                print(f"  - {file.name} ({size_kb:.1f} KB)")
        
        print(f"\nThis simulation demonstrated:")
        print(f"  - Loading a molecular system from GSD file")
        print(f"  - Setting up cavity-molecule coupling") 
        print(f"  - Running a short cavity MD simulation")
        print(f"  - Using the CavityMDSimulation framework")
        
        return True
    else:
        print("❌ Simulation failed!")
        return False

if __name__ == "__main__":
    try:
        success = run_basic_cavity_simulation()
        if success:
            print(f"\n🎉 Basic cavity simulation example completed successfully!")
        else:
            print(f"\n💥 Basic cavity simulation example failed!")
            sys.exit(1)
    except Exception as e:
        print(f"\n❌ Error running basic cavity simulation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1) 