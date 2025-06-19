#!/usr/bin/env python3
"""
Example usage of the BussiReservoir thermostat plugin.

This script demonstrates how to:
1. Set up a simple LJ liquid simulation
2. Use the BussiReservoir thermostat to track reservoir energy
3. Log and analyze the reservoir energy over time
"""

import hoomd
import hoomd.bussi_reservoir as bussi_res
import numpy as np
import matplotlib.pyplot as plt

def create_lj_liquid():
    """Create a simple Lennard-Jones liquid system."""
    # Simulation parameters
    device = hoomd.device.CPU()
    simulation = hoomd.Simulation(device=device, seed=42)
    
    # System parameters
    N = 500  # Number of particles
    L = 15.0  # Box size
    density = N / L**3
    
    # Create initial configuration
    snapshot = hoomd.Snapshot(device.communicator)
    if snapshot.communicator.rank == 0:
        box = [L, L, L, 0, 0, 0]
        snapshot.configuration.box = box
        snapshot.particles.N = N
        
        # Random positions
        snapshot.particles.position[:] = np.random.uniform(-L/2, L/2, (N, 3))
        
        # Maxwell-Boltzmann velocities at temperature kT=1.5
        kT = 1.5
        snapshot.particles.velocity[:] = np.random.normal(0, np.sqrt(kT), (N, 3))
        
        # Single particle type
        snapshot.particles.typeid[:] = 0
        snapshot.particles.types = ['A']
    
    simulation.create_state_from_snapshot(snapshot)
    return simulation

def main():
    """Main simulation function."""
    print("Setting up LJ liquid simulation with BussiReservoir thermostat...")
    
    # Create simulation
    simulation = create_lj_liquid()
    
    # Create BussiReservoir thermostat
    kT = 1.5
    dt = 0.002
    tau = dt * 50  # Thermostat time constant
    
    bussi = bussi_res.thermostats.BussiReservoir(kT=kT, tau=tau)
    
    # Set up integration method
    nve = hoomd.md.methods.ConstantVolume(
        filter=hoomd.filter.All(),
        thermostat=bussi
    )
    
    integrator = hoomd.md.Integrator(dt=dt, methods=[nve])
    simulation.operations.integrator = integrator
    
    # Add Lennard-Jones potential
    lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4))
    lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)
    lj.r_cut[('A', 'A')] = 2.5
    simulation.operations.integrator.forces.append(lj)
    
    # Set up logging
    logger = hoomd.logging.Logger()
    
    # Add thermodynamic quantities
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    simulation.operations.computes.append(thermo)
    logger.add(thermo, quantities=['kinetic_energy', 'potential_energy', 'kinetic_temperature'])
    
    # Add reservoir energy quantities
    logger.add(bussi, quantities=[
        'total_reservoir_energy',
        'reservoir_energy_translational',
        'reservoir_energy_rotational',
        'instantaneous_reservoir_total'
    ])
    
    # Set up file output
    table_writer = hoomd.write.Table(
        output=open('reservoir_energy.log', 'w'),
        trigger=hoomd.trigger.Periodic(100),
        logger=logger
    )
    simulation.operations.writers.append(table_writer)
    
    print(f"Initial parameters:")
    print(f"  Temperature: {kT}")
    print(f"  Time step: {dt}")
    print(f"  Thermostat tau: {tau}")
    print(f"  Number of particles: {simulation.state.N_particles}")
    
    # Initialize simulation
    print("\nInitializing simulation...")
    simulation.run(0)
    
    # Equilibration phase
    print("Running equilibration (5000 steps)...")
    simulation.run(5000)
    
    print(f"After equilibration:")
    print(f"  Total reservoir energy: {bussi.total_reservoir_energy:.4f}")
    print(f"  Translational component: {bussi.reservoir_energy_translational:.4f}")
    print(f"  Rotational component: {bussi.reservoir_energy_rotational:.4f}")
    
    # Reset counters for production run
    print("\nResetting reservoir energy counters...")
    bussi.reset_reservoir_energy()
    
    # Production phase
    print("Running production (10000 steps)...")
    n_steps = 10000
    
    # Collect data during production
    reservoir_energies = []
    timesteps = []
    
    for i in range(100):  # 100 measurements
        simulation.run(100)  # Run 100 steps between measurements
        reservoir_energies.append(bussi.total_reservoir_energy)
        timesteps.append((i + 1) * 100)
    
    print(f"\nProduction run completed!")
    print(f"Final reservoir energies:")
    print(f"  Total: {bussi.total_reservoir_energy:.4f}")
    print(f"  Translational: {bussi.reservoir_energy_translational:.4f}")
    print(f"  Rotational: {bussi.reservoir_energy_rotational:.4f}")
    
    # Calculate some statistics
    avg_temp = thermo.kinetic_temperature
    total_energy = thermo.kinetic_energy + thermo.potential_energy
    
    print(f"\nFinal thermodynamic state:")
    print(f"  Temperature: {avg_temp:.4f} (target: {kT})")
    print(f"  Total energy: {total_energy:.4f}")
    print(f"  Kinetic energy: {thermo.kinetic_energy:.4f}")
    print(f"  Potential energy: {thermo.potential_energy:.4f}")
    
    # Plot results if matplotlib is available
    try:
        plt.figure(figsize=(10, 6))
        
        plt.subplot(1, 2, 1)
        plt.plot(timesteps, reservoir_energies, 'b-', linewidth=2)
        plt.xlabel('Time steps')
        plt.ylabel('Cumulative Reservoir Energy')
        plt.title('Reservoir Energy vs Time')
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 2, 2)
        # Calculate instantaneous rate of energy dump
        energy_rates = np.diff(reservoir_energies)
        plt.plot(timesteps[1:], energy_rates, 'r-', alpha=0.7)
        plt.xlabel('Time steps')
        plt.ylabel('Energy Rate (per 100 steps)')
        plt.title('Rate of Energy Dump to Reservoir')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('reservoir_energy_analysis.png', dpi=150, bbox_inches='tight')
        print(f"\nPlot saved as 'reservoir_energy_analysis.png'")
        
    except ImportError:
        print("\nMatplotlib not available - skipping plots")
    
    print(f"\nData saved to 'reservoir_energy.log'")
    print("Analysis complete!")

if __name__ == '__main__':
    main() 