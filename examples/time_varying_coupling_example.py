#!/usr/bin/env python3
"""
Time-Varying Coupling Strength Example
======================================

This example demonstrates how to use time-varying coupling strength in cavity molecular dynamics
simulations using hoomd.variant. The coupling strength g(t) can vary over time according to 
different patterns: constant, ramp, cycle, or power functions.

Usage:
    python time_varying_coupling_example.py --variant constant --coupling 1e-3
    python time_varying_coupling_example.py --variant ramp --coupling 1e-3 --final_coupling 2e-3 --ramp_time 10000
    python time_varying_coupling_example.py --variant cycle --coupling 1e-3 --final_coupling 2e-3 --cycle_time 5000
"""

import hoomd
import hoomd.variant
import numpy as np
import argparse
import matplotlib.pyplot as plt
import sys
import os

# Add the cav-hoomd module to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cavitymd.forces import CavityForce
# Use HOOMD's standard Langevin integrator


def create_test_system():
    """Create a simple test system: water molecule + cavity mode."""
    
    # Create simple water molecule
    snapshot = hoomd.Snapshot()
    
    if hoomd.communicator.rank == 0:
        # Water molecule: O-H-H
        snapshot.particles.N = 4  # O, H, H, L (cavity mode)
        snapshot.particles.types = ['O', 'H', 'L']
        
        # Positions: water molecule + cavity particle at origin
        snapshot.particles.position = [
            [0.0, 0.0, 0.0],      # O at origin
            [0.757, 0.586, 0.0],  # H 
            [-0.757, 0.586, 0.0], # H
            [0.0, 0.0, 0.0]       # L (cavity mode) at origin
        ]
        
        # Masses: realistic masses for water + light photon mass
        snapshot.particles.mass = [15.999, 1.008, 1.008, 1.0]  # O, H, H, L
        
        # Charges: partial charges for water
        snapshot.particles.charge = [-0.834, 0.417, 0.417, 0.0]  # O, H, H, L
        
        # Velocities: small random velocities
        np.random.seed(42)
        snapshot.particles.velocity = np.random.normal(0, 0.1, (4, 3))
        snapshot.particles.velocity[3] = [0, 0, 0]  # Cavity mode starts at rest
        
        # Type IDs
        snapshot.particles.typeid = [0, 1, 1, 2]  # O, H, H, L
        
        # Box
        snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
    
    return snapshot


def setup_variant_coupling(variant_type, base_coupling, **kwargs):
    """Setup different types of time-varying coupling strength."""
    
    if variant_type == "constant":
        coupling_variant = hoomd.variant.Constant(base_coupling)
        description = f"Constant coupling: g = {base_coupling:.3e} a.u."
        
    elif variant_type == "ramp":
        final_coupling = kwargs.get('final_coupling', 2 * base_coupling)
        ramp_time = kwargs.get('ramp_time', 10000)
        start_time = kwargs.get('start_time', 0)
        
        coupling_variant = hoomd.variant.Ramp(
            A=base_coupling,
            B=final_coupling, 
            t_start=start_time,
            t_ramp=ramp_time
        )
        description = f"Ramp coupling: g = {base_coupling:.3e} → {final_coupling:.3e} a.u. over {ramp_time} steps"
        
    elif variant_type == "cycle":
        final_coupling = kwargs.get('final_coupling', 2 * base_coupling)
        cycle_time = kwargs.get('cycle_time', 5000)
        hold_time = kwargs.get('hold_time', 1000)
        start_time = kwargs.get('start_time', 0)
        
        coupling_variant = hoomd.variant.Cycle(
            A=base_coupling,
            B=final_coupling,
            t_start=start_time,
            t_A=hold_time,
            t_AB=cycle_time,
            t_B=hold_time,
            t_BA=cycle_time
        )
        description = f"Cycle coupling: g = {base_coupling:.3e} ↔ {final_coupling:.3e} a.u., cycle = {cycle_time} steps"
        
    elif variant_type == "power":
        final_coupling = kwargs.get('final_coupling', 2 * base_coupling)
        ramp_time = kwargs.get('ramp_time', 10000)
        power = kwargs.get('power', 0.5)
        start_time = kwargs.get('start_time', 0)
        
        coupling_variant = hoomd.variant.Power(
            A=base_coupling,
            B=final_coupling,
            power=power,
            t_start=start_time,
            t_ramp=ramp_time
        )
        description = f"Power coupling: g = {base_coupling:.3e} → {final_coupling:.3e} a.u., power = {power}"
        
    else:
        raise ValueError(f"Unknown variant type: {variant_type}")
    
    return coupling_variant, description


def main():
    parser = argparse.ArgumentParser(description="Time-varying coupling strength example")
    parser.add_argument("--variant", choices=["constant", "ramp", "cycle", "power"], 
                       default="ramp", help="Type of coupling variation")
    parser.add_argument("--coupling", type=float, default=1e-3,
                       help="Base coupling strength in a.u.")
    parser.add_argument("--final_coupling", type=float, default=None,
                       help="Final coupling strength (for ramp/cycle/power)")
    parser.add_argument("--ramp_time", type=int, default=10000,
                       help="Ramp duration in timesteps")
    parser.add_argument("--cycle_time", type=int, default=5000,
                       help="Cycle duration in timesteps")
    parser.add_argument("--hold_time", type=int, default=1000,
                       help="Hold time in cycle variant")
    parser.add_argument("--power", type=float, default=0.5,
                       help="Power for power variant")
    parser.add_argument("--start_time", type=int, default=0,
                       help="Start time for variant")
    parser.add_argument("--steps", type=int, default=20000,
                       help="Number of simulation steps")
    parser.add_argument("--dt", type=float, default=0.001,
                       help="Timestep size in a.u.")
    parser.add_argument("--temp", type=float, default=100.0,
                       help="Temperature in K")
    parser.add_argument("--freq", type=float, default=2000.0,
                       help="Cavity frequency in cm^-1")
    parser.add_argument("--plot", action="store_true",
                       help="Generate plots of coupling strength vs time")
    
    args = parser.parse_args()
    
    # Convert units
    kB_au = 3.1668e-6  # Boltzmann constant in a.u./K
    kT = args.temp * kB_au
    
    # Convert cavity frequency from cm^-1 to a.u.
    omegac_au = args.freq * 4.556e-6  # cm^-1 to a.u.
    
    # Set default final coupling if not specified
    if args.final_coupling is None:
        args.final_coupling = 2 * args.coupling
    
    print("Time-Varying Coupling Strength Example")
    print("=====================================")
    print(f"Base coupling: {args.coupling:.3e} a.u.")
    print(f"Temperature: {args.temp} K (kT = {kT:.6f} a.u.)")
    print(f"Cavity frequency: {args.freq} cm^-1 ({omegac_au:.6f} a.u.)")
    print()
    
    # Setup coupling variant
    coupling_variant, description = setup_variant_coupling(
        args.variant, 
        args.coupling,
        final_coupling=args.final_coupling,
        ramp_time=args.ramp_time,
        cycle_time=args.cycle_time,
        hold_time=args.hold_time,
        power=args.power,
        start_time=args.start_time
    )
    
    print(f"Coupling variant: {description}")
    print()
    
    # Create simulation
    cpu = hoomd.device.CPU()
    sim = hoomd.Simulation(device=cpu, seed=42)
    
    # Create system
    snapshot = create_test_system()
    sim.create_state_from_snapshot(snapshot)
    
    print(f"System created: {sim.state.N_particles} particles")
    print(f"Particle types: {sim.state.particle_types}")
    
    # Setup integrator
    integrator = hoomd.md.Integrator(dt=args.dt, integrate_rotational_dof=True)
    
    # Thermostat for molecular particles
    molecular_filter = hoomd.filter.Type(['O', 'H'])
    molecular_langevin = hoomd.md.methods.Langevin(
        filter=molecular_filter,
        kT=kT,
        default_gamma=1.0,
        tally_reservoir_energy=True
    )
    
    # Thermostat for cavity mode
    cavity_filter = hoomd.filter.Type(['L'])
    cavity_langevin = hoomd.md.methods.Langevin(
        filter=cavity_filter,
        kT=kT,
        default_gamma=0.1,  # Lighter damping for cavity mode
        tally_reservoir_energy=True
    )
    
    integrator.methods = [molecular_langevin, cavity_langevin]
    
    # Setup cavity force with time-varying coupling
    cavity_force = CavityForce(
        kvector=[0, 0, 1],
        couplstr=coupling_variant,  # Time-varying coupling strength
        omegac=omegac_au,
        phmass=1.0,
        damping_ratio=0.01
    )
    
    integrator.forces = [cavity_force]
    sim.operations.integrator = integrator
    
    print(f"Integrator setup complete")
    print(f"Initial coupling strength: {coupling_variant(0):.6f} a.u.")
    
    # Setup logging
    logger = hoomd.logging.Logger()
    logger.add(sim, quantities=['timestep'])
    logger.add(cavity_force, quantities=['harmonic_energy', 'coupling_energy', 'dipole_self_energy', 'total_cavity_energy'])
    logger.add(molecular_langevin, quantities=['reservoir_energy'])
    logger.add(cavity_langevin, quantities=['reservoir_energy'])
    
    # Custom logger for coupling strength
    class CouplingLogger:
        def __init__(self, variant):
            self.variant = variant
            
        @property
        def coupling_strength(self):
            return self.variant(sim.timestep)
    
    coupling_logger = CouplingLogger(coupling_variant)
    logger.add(coupling_logger, quantities=['coupling_strength'])
    
    # Setup file writer
    gsd_writer = hoomd.write.GSD(
        filename=f"time_varying_coupling_{args.variant}.gsd",
        trigger=hoomd.trigger.Periodic(100),
        mode='wb'
    )
    sim.operations.writers.append(gsd_writer)
    
    # Setup table writer for analysis
    table_writer = hoomd.write.Table(
        output=open(f"time_varying_coupling_{args.variant}.log", 'w'),
        trigger=hoomd.trigger.Periodic(50),
        logger=logger
    )
    sim.operations.writers.append(table_writer)
    
    print("Starting simulation...")
    print()
    
    # Run simulation
    sim.run(0)  # Initialize
    
    # Store coupling values for plotting
    if args.plot:
        coupling_values = []
        time_values = []
    
    # Run in chunks to monitor progress
    chunk_size = args.steps // 10
    for i in range(10):
        current_step = (i + 1) * chunk_size
        sim.run(chunk_size)
        
        current_coupling = coupling_variant(sim.timestep)
        total_energy = cavity_force.total_cavity_energy
        
        print(f"Step {sim.timestep:6d}: g = {current_coupling:.6f} a.u., "
              f"E_total = {total_energy:.6f} a.u.")
        
        if args.plot:
            coupling_values.append(current_coupling)
            time_values.append(sim.timestep)
    
    print()
    print(f"Simulation completed: {sim.timestep} steps")
    print(f"Final coupling strength: {coupling_variant(sim.timestep):.6f} a.u.")
    
    # Generate plots if requested
    if args.plot:
        plt.figure(figsize=(12, 4))
        
        # Plot 1: Coupling strength vs time
        plt.subplot(1, 2, 1)
        full_time = np.arange(0, args.steps + 1, 100)
        full_coupling = [coupling_variant(t) for t in full_time]
        
        plt.plot(full_time, full_coupling, 'b-', linewidth=2, label='Coupling strength')
        plt.scatter(time_values, coupling_values, color='red', s=30, alpha=0.7, label='Simulation points')
        plt.xlabel('Time step')
        plt.ylabel('Coupling strength g(t) [a.u.]')
        plt.title(f'Time-varying coupling: {args.variant}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Plot 2: Coupling strength range
        plt.subplot(1, 2, 2)
        plt.plot(full_time, full_coupling, 'b-', linewidth=2)
        plt.fill_between(full_time, full_coupling, alpha=0.3)
        plt.xlabel('Time step')
        plt.ylabel('Coupling strength g(t) [a.u.]')
        plt.title('Coupling strength variation')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plot_filename = f"coupling_variation_{args.variant}.png"
        plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
        print(f"Plot saved: {plot_filename}")
        plt.show()


if __name__ == "__main__":
    main() 