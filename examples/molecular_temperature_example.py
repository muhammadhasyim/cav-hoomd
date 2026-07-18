#!/usr/bin/env python3
"""
Example script showing how to use DiatomicMolecularTemperatures tracker.

This demonstrates how to calculate translational, rotational, and vibrational
temperatures for a system of diatomic molecules.
"""

import sys
sys.path.insert(0, '../src')

from cavitymd.molecular_temperatures import DiatomicMolecularTemperatures
from cavitymd.analysis import ElapsedTimeTracker
import hoomd
import numpy as np


def main():
    """Run example simulation with molecular temperature tracking."""
    
    # Initialize HOOMD
    device = hoomd.device.CPU()
    sim = hoomd.Simulation(device=device, seed=42)
    
    # Load an existing GSD file with diatomic molecules
    # Replace this with your actual GSD file path
    gsd_file = "molecular-0.gsd"  # Your molecular system file
    
    try:
        sim.create_state_from_gsd(filename=gsd_file)
    except:
        print(f"Error: Could not load {gsd_file}")
        print("Please provide a valid GSD file with diatomic molecules")
        return
    
    # Setup integrator and thermostat
    dt = 0.0001  # timestep in reduced units
    integrator = hoomd.md.Integrator(dt=dt)
    
    # Add thermostat
    langevin = hoomd.md.methods.Langevin(
        filter=hoomd.filter.All(),
        kT=0.0009505  # 300 K in Hartree (300 * 3.167e-6)
    )
    integrator.methods.append(langevin)
    
    # Add bond forces (harmonic bonds for diatomic molecules)
    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params['O-O'] = dict(k=2*0.36602, r0=2.281655158)
    harmonic.params['N-N'] = dict(k=2*0.71625, r0=2.0743522177)
    integrator.forces.append(harmonic)
    
    # Add pair forces (Lennard-Jones)
    cell = hoomd.md.nlist.Cell(buffer=1.0, exclusions=('bond',))
    lj = hoomd.md.pair.LJ(nlist=cell, mode='shift')
    lj.params[('O', 'O')] = dict(epsilon=0.00016685201, sigma=6.230426584)
    lj.r_cut[('O', 'O')] = 15.0
    lj.params[('N', 'N')] = dict(epsilon=0.000083426, sigma=5.48277488)
    lj.r_cut[('N', 'N')] = 15.0
    lj.params[('N', 'O')] = dict(epsilon=0.00025027802, sigma=4.9832074319)
    lj.r_cut[('N', 'O')] = 15.0
    integrator.forces.append(lj)
    
    sim.operations.integrator = integrator
    
    # Setup time tracker
    time_tracker = ElapsedTimeTracker(
        timestep=dt,
        output_file="elapsed_time.txt",
        output_period_ps=1.0
    )
    sim.operations.writers.append(time_tracker)
    
    # Setup molecular temperature tracker
    mol_temp_tracker = DiatomicMolecularTemperatures(
        simulation=sim,
        time_tracker=time_tracker,
        output_period_ps=1.0,  # Output every 1 ps
        output_file="molecular_temperatures.csv",
        debug=True  # Enable debug output
    )
    sim.operations.writers.append(mol_temp_tracker)
    
    # Run simulation
    print("\nRunning simulation with molecular temperature tracking...")
    print("Output will be written to: molecular_temperatures.csv\n")
    
    # Run for a short equilibration
    sim.run(10000)  # 1 ps
    
    # Print current temperatures
    print("\nCurrent Molecular Temperatures:")
    print(f"  Translational:  {mol_temp_tracker.translational_temp:.2f} K")
    print(f"  Rotational:     {mol_temp_tracker.rotational_temp:.2f} K")
    print(f"  Vibrational:    {mol_temp_tracker.vibrational_temp:.2f} K")
    print(f"  Total Kinetic:  {mol_temp_tracker.total_kinetic_temp:.2f} K")
    
    print("\nData saved to: molecular_temperatures.csv")
    print("\nColumn definitions:")
    print("  T_trans: Translational temperature (center of mass motion)")
    print("  T_rot: Rotational temperature (rotation around COM)")
    print("  T_vib: Vibrational temperature (motion along bond)")
    print("  T_kinetic_total: Total kinetic temperature")
    print("  T_*_O2: Temperatures for O-O molecules only")
    print("  T_*_N2: Temperatures for N-N molecules only")


if __name__ == "__main__":
    main()

