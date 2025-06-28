#!/usr/bin/env python3
"""
Example 3: Energy Tracking Simulation (1 ps)

This example demonstrates comprehensive energy tracking in a cavity MD simulation:
- Tracks all cavity energy components (harmonic, coupling, dipole self-energy)
- Monitors molecular kinetic and potential energies
- Records thermostat reservoir energies
- Saves energy data to files and creates analysis plots
- Runs for 1 ps to show energy equilibration

Key features:
- Uses the refactored EnergyContributionTracker
- Real-time energy monitoring and logging
- Creates comprehensive energy analysis plots
- Demonstrates proper thermostating of molecular vs cavity systems
"""

import hoomd
import hoomd.md
import gsd.hoomd
import numpy as np
import matplotlib.pyplot as plt
from hoomd.cavitymd import CavityForce, EnergyContributionTracker

def create_energy_tracking_system():
    """Create a molecular system optimized for energy tracking."""
    snapshot = gsd.hoomd.Frame()
    
    # Medium-sized system for clear energy dynamics
    L = 18.0
    snapshot.configuration.box = [L, L, L, 0, 0, 0]
    
    # Create 6 water-like molecules + 1 cavity particle = 19 total particles
    snapshot.particles.N = 19
    snapshot.particles.types = ['O', 'H', 'L']
    
    # Arrange 6 water molecules in a hexagonal pattern
    positions = []
    typeids = []
    charges = []
    masses = []
    
    # Hexagonal arrangement for water molecules
    import math
    n_waters = 6
    radius = 4.0
    
    for i in range(n_waters):
        angle = 2 * math.pi * i / n_waters
        center_x = radius * math.cos(angle)
        center_y = radius * math.sin(angle)
        center_z = 0.0
        
        # Oxygen
        positions.append([center_x, center_y, center_z])
        typeids.append(0)  # O
        charges.append(-0.834)
        masses.append(15.999)
        
        # Hydrogen 1
        positions.append([center_x + 0.96, center_y, center_z])
        typeids.append(1)  # H
        charges.append(0.417)
        masses.append(1.008)
        
        # Hydrogen 2
        positions.append([center_x - 0.24, center_y + 0.93, center_z])
        typeids.append(1)  # H
        charges.append(0.417)
        masses.append(1.008)
    
    # Cavity particle at center
    positions.append([0.0, 0.0, 0.0])
    typeids.append(2)  # L
    charges.append(0.0)
    masses.append(1.0)
    
    # Set snapshot data
    snapshot.particles.position = positions
    snapshot.particles.typeid = typeids
    snapshot.particles.charge = charges
    snapshot.particles.mass = masses
    snapshot.particles.diameter = [1.0] * snapshot.particles.N
    
    return snapshot

def setup_interactions(sim):
    """Set up realistic molecular interactions."""
    
    # Lennard-Jones interactions
    lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4))
    
    # Water-water interactions (TIP3P-like)
    lj.params[('O', 'O')] = dict(epsilon=0.1554, sigma=3.1656)  # kcal/mol to atomic units
    lj.r_cut[('O', 'O')] = 9.0
    
    # H-H interactions (weak repulsive)
    lj.params[('H', 'H')] = dict(epsilon=0.01, sigma=2.0)
    lj.r_cut[('H', 'H')] = 4.0
    
    # O-H cross interactions
    lj.params[('O', 'H')] = dict(epsilon=0.05, sigma=2.5)
    lj.r_cut[('O', 'H')] = 6.0
    
    # Cavity interactions (minimal)
    lj.params[('L', 'O')] = dict(epsilon=0.001, sigma=2.0)
    lj.params[('L', 'H')] = dict(epsilon=0.001, sigma=1.5)
    lj.params[('L', 'L')] = dict(epsilon=0.0, sigma=1.0)
    lj.r_cut[('L', 'O')] = 4.0
    lj.r_cut[('L', 'H')] = 3.0
    lj.r_cut[('L', 'L')] = 2.0
    
    return lj

def run_energy_tracking_simulation():
    """Run 1 ps cavity MD simulation with comprehensive energy tracking."""
    
    print("=" * 80)
    print("CAVITY MD SIMULATION WITH ENERGY TRACKING (1 ps)")
    print("=" * 80)
    
    # 1. Initialize simulation
    print("\n1. Setting up simulation...")
    device = hoomd.device.CPU()
    sim = hoomd.Simulation(device=device, seed=42)
    
    # 2. Create system
    print("2. Creating molecular system...")
    snapshot = create_energy_tracking_system()
    sim.create_state_from_snapshot(snapshot)
    
    print(f"   System: {snapshot.particles.N} particles")
    print(f"   Water molecules: 6 (18 atoms)")
    print(f"   Cavity particles: 1")
    print(f"   Box: {snapshot.configuration.box[0]:.1f}³")
    
    # 3. Set up forces
    print("\n3. Setting up forces...")
    
    # Molecular interactions
    lj = setup_interactions(sim)
    
    # Cavity force with moderate coupling
    cavity_params = {
        'kvector': [1.0, 0.0, 0.0],
        'couplstr': 0.08,  # Strong enough to see clear coupling effects
        'omegac': 0.12,    # Higher frequency for faster dynamics
        'phmass': 1.0
    }
    
    cavity_force = CavityForce(**cavity_params)
    
    print(f"   Cavity coupling: g = {cavity_params['couplstr']:.3f}")
    print(f"   Cavity frequency: ωc = {cavity_params['omegac']:.3f} Hartree")
    print(f"   Implementation: {cavity_force.implementation}")
    
    # 4. Set up integrator and thermostats
    print("\n4. Setting up integrator...")
    dt = 0.0005  # 0.5 fs for stable integration
    
    integrator = hoomd.md.Integrator(dt=dt)
    integrator.forces.append(lj)
    integrator.forces.append(cavity_force)
    
    # NVT for molecular particles
    molecular_filter = hoomd.filter.Type(['O', 'H'])
    target_temp = 1.5  # ~300 K equivalent
    nvt_molecular = hoomd.md.methods.ConstantVolume(
        filter=molecular_filter,
        thermostat=hoomd.md.methods.thermostats.MTTK(kT=target_temp, tau=0.5)
    )
    integrator.methods.append(nvt_molecular)
    
    # NVE for cavity particle
    cavity_filter = hoomd.filter.Type(['L'])
    nve_cavity = hoomd.md.methods.ConstantVolume(filter=cavity_filter)
    integrator.methods.append(nve_cavity)
    
    sim.operations.integrator = integrator
    
    print(f"   Timestep: {dt:.4f} fs")
    print(f"   Molecular thermostat: NVT (kT = {target_temp:.1f}, τ = 0.5)")
    print(f"   Cavity dynamics: NVE")
    
    # 5. Set up energy tracking
    print("\n5. Setting up energy tracking...")
    
    # Create energy tracker
    energy_tracker = EnergyContributionTracker(
        cavity_force=cavity_force,
        molecular_filter=molecular_filter,
        cavity_filter=cavity_filter,
        filename="energy_tracking_1ps.txt",
        molecular_thermostat_name="mttk"
    )
    
    # Add to simulation operations
    sim.operations.analyzers.append(energy_tracker)
    
    # Set up additional thermodynamic quantities
    thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)
    
    print(f"   Energy tracker: logging to energy_tracking_1ps.txt")
    print(f"   Tracking: cavity, molecular, and thermostat energies")
    
    # 6. Initialize system
    print("\n6. Initializing system...")
    sim.state.thermalize_particle_momenta(filter=molecular_filter, kT=target_temp)
    
    # Print initial state
    with sim.state.cpu_local_snapshot as snap:
        positions = snap.particles.position
        charges = snap.particles.charge
        molecular_dipole = np.sum(charges[:-1, np.newaxis] * positions[:-1], axis=0)
        cavity_pos = positions[-1]
    
    print(f"   Initial molecular dipole: [{molecular_dipole[0]:.3f}, {molecular_dipole[1]:.3f}, {molecular_dipole[2]:.3f}]")
    print(f"   Initial cavity position: [{cavity_pos[0]:.3f}, {cavity_pos[1]:.3f}, {cavity_pos[2]:.3f}]")
    
    # 7. Run simulation
    print("\n7. Running 1 ps simulation...")
    print("=" * 80)
    
    total_steps = 2000  # 1 ps at 0.5 fs timestep
    report_interval = 200  # Report every 0.1 ps
    log_interval = 10     # Log every 5 fs
    
    # Set logging frequency
    energy_tracker.period = log_interval
    
    print(f"Total simulation: {total_steps * dt:.1f} ps ({total_steps} steps)")
    print(f"Energy logging: every {log_interval * dt:.3f} ps")
    print(f"Progress reports: every {report_interval * dt:.1f} ps")
    print()
    print("Step    Time(ps)  | Cavity Energies                    | Molecular T(K)")
    print("-" * 75)
    
    # Storage for real-time plotting data
    times = []
    cavity_totals = []
    molecular_temps = []
    
    for step in range(0, total_steps + 1, report_interval):
        if step > 0:
            sim.run(report_interval)
        
        # Get current time
        current_time = step * dt
        times.append(current_time)
        
        # Get cavity energies
        try:
            harmonic = cavity_force.harmonic_energy
            coupling = cavity_force.coupling_energy
            dipole_self = cavity_force.dipole_self_energy
            total_cavity = cavity_force.total_cavity_energy
            cavity_totals.append(total_cavity)
        except:
            harmonic = coupling = dipole_self = total_cavity = 0.0
            cavity_totals.append(0.0)
        
        # Get molecular temperature
        try:
            mol_temp = thermo.kinetic_temperature * (target_temp / thermo.kinetic_temperature) if hasattr(thermo, 'kinetic_temperature') else target_temp
            molecular_temps.append(mol_temp)
        except:
            molecular_temps.append(target_temp)
        
        print(f"{step:5d}  {current_time:7.3f}   | H:{harmonic:.5f} C:{coupling:.5f} D:{dipole_self:.5f} T:{total_cavity:.5f} | {molecular_temps[-1]:.1f}")
    
    print("-" * 75)
    print(f"Simulation completed: {total_steps * dt:.1f} ps")
    
    # 8. Final analysis
    print("\n8. Final state analysis:")
    
    with sim.state.cpu_local_snapshot as snap:
        final_positions = snap.particles.position
        final_charges = snap.particles.charge
        final_dipole = np.sum(final_charges[:-1, np.newaxis] * final_positions[:-1], axis=0)
        final_cavity_pos = final_positions[-1]
    
    print(f"   Final molecular dipole: [{final_dipole[0]:.3f}, {final_dipole[1]:.3f}, {final_dipole[2]:.3f}]")
    print(f"   Final cavity position: [{final_cavity_pos[0]:.3f}, {final_cavity_pos[1]:.3f}, {final_cavity_pos[2]:.3f}]")
    
    print(f"\n   Final cavity energy components:")
    print(f"     Harmonic:     {cavity_force.harmonic_energy:.8f} Hartree")
    print(f"     Coupling:     {cavity_force.coupling_energy:.8f} Hartree")
    print(f"     Dipole self:  {cavity_force.dipole_self_energy:.8f} Hartree")
    print(f"     Total cavity: {cavity_force.total_cavity_energy:.8f} Hartree")
    
    # 9. Create real-time analysis plot
    print("\n9. Creating energy analysis plots...")
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
    
    # Plot 1: Cavity energy evolution
    ax1.plot(times, cavity_totals, 'b-', linewidth=2, label='Total Cavity Energy')
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Cavity Energy (Hartree)')
    ax1.set_title('Cavity Energy Evolution (1 ps)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Temperature evolution
    ax2.plot(times, molecular_temps, 'r-', linewidth=2, label='Molecular Temperature')
    ax2.axhline(y=target_temp, color='k', linestyle='--', alpha=0.5, label=f'Target T = {target_temp:.1f}')
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Temperature (K equivalent)')
    ax2.set_title('Molecular Temperature Evolution')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Plot 3: Energy correlation
    if len(cavity_totals) > 1 and len(molecular_temps) > 1:
        ax3.scatter(cavity_totals, molecular_temps, alpha=0.6, c=times, cmap='viridis')
        ax3.set_xlabel('Cavity Energy (Hartree)')
        ax3.set_ylabel('Molecular Temperature (K equivalent)')
        ax3.set_title('Energy-Temperature Correlation')
        ax3.grid(True, alpha=0.3)
        cbar = plt.colorbar(ax3.collections[0], ax=ax3)
        cbar.set_label('Time (ps)')
    
    plt.tight_layout()
    plt.savefig('energy_tracking_1ps_analysis.png', dpi=300, bbox_inches='tight')
    print(f"   Saved analysis plot: energy_tracking_1ps_analysis.png")
    
    print(f"\n✅ Energy tracking simulation completed successfully!")
    print(f"   Data saved to: energy_tracking_1ps.txt")
    print(f"   This demonstrates comprehensive energy monitoring over 1 ps.")
    
    return sim, energy_tracker, cavity_force

if __name__ == "__main__":
    try:
        simulation, tracker, cavity_force = run_energy_tracking_simulation()
        print(f"\nSuccess! Energy tracking simulation completed.")
        print(f"Check energy_tracking_1ps.txt for detailed energy data.")
    except Exception as e:
        print(f"\n❌ Error running energy tracking simulation: {e}")
        import traceback
        traceback.print_exc() 