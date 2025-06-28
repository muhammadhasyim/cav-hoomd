#!/usr/bin/env python3
"""
Example 4: Structure and Dipole Tracking Simulation (10 ps)

This example demonstrates advanced analysis capabilities in cavity MD:
- Tracks intermediate scattering function F(k,t) for structural analysis
- Monitors dipole moment autocorrelation function
- Records cavity mode dynamics over longer timescales
- Runs for 10 ps to capture relevant correlation timescales
- Creates comprehensive structure analysis plots

Key features:
- Uses DensityCorrelationTracker for F(k,t) calculation
- Uses DipoleAutocorrelation for dipole dynamics
- Long-time simulation to see decorrelation
- Advanced analysis of cavity-molecular coupling effects
"""

import hoomd
import hoomd.md
import gsd.hoomd
import numpy as np
import matplotlib.pyplot as plt
from hoomd.cavitymd import (
    CavityForce, DensityCorrelationTracker, 
    DipoleAutocorrelation, CavityModeTracker
)

def create_structure_tracking_system():
    """Create a larger molecular system for structure analysis."""
    snapshot = gsd.hoomd.Frame()
    
    # Larger system for meaningful structure factors
    L = 25.0
    snapshot.configuration.box = [L, L, L, 0, 0, 0]
    
    # Create 12 water-like molecules + 1 cavity particle = 37 total particles
    snapshot.particles.N = 37
    snapshot.particles.types = ['O', 'H', 'L']
    
    # Arrange molecules in a 3x4 grid with some randomization
    positions = []
    typeids = []
    charges = []
    masses = []
    
    import math
    np.random.seed(42)  # For reproducible arrangements
    
    # Grid arrangement with slight randomization
    n_waters = 12
    grid_x, grid_y = 4, 3
    spacing = 4.5
    
    for i in range(n_waters):
        grid_i = i % grid_x
        grid_j = i // grid_x
        
        # Base position on grid
        base_x = (grid_i - (grid_x-1)/2) * spacing
        base_y = (grid_j - (grid_y-1)/2) * spacing
        base_z = 0.0
        
        # Add small random displacement
        offset_x = np.random.uniform(-0.5, 0.5)
        offset_y = np.random.uniform(-0.5, 0.5)
        offset_z = np.random.uniform(-0.2, 0.2)
        
        center_x = base_x + offset_x
        center_y = base_y + offset_y
        center_z = base_z + offset_z
        
        # Add random rotation to water orientation
        angle = np.random.uniform(0, 2*math.pi)
        
        # Oxygen
        positions.append([center_x, center_y, center_z])
        typeids.append(0)  # O
        charges.append(-0.834)
        masses.append(15.999)
        
        # Hydrogen 1 (rotated)
        h1_x = center_x + 0.96 * math.cos(angle)
        h1_y = center_y + 0.96 * math.sin(angle)
        positions.append([h1_x, h1_y, center_z])
        typeids.append(1)  # H
        charges.append(0.417)
        masses.append(1.008)
        
        # Hydrogen 2 (rotated)
        h2_angle = angle + 2.09  # ~120 degrees
        h2_x = center_x + 0.96 * math.cos(h2_angle)
        h2_y = center_y + 0.96 * math.sin(h2_angle)
        positions.append([h2_x, h2_y, center_z])
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

def setup_enhanced_interactions(sim):
    """Set up realistic molecular interactions for structure analysis."""
    
    # Lennard-Jones interactions with more realistic parameters
    lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.5))
    
    # Enhanced water-water interactions
    lj.params[('O', 'O')] = dict(epsilon=0.1554, sigma=3.1656)
    lj.r_cut[('O', 'O')] = 10.0
    
    # H-H interactions (prevent overlap)
    lj.params[('H', 'H')] = dict(epsilon=0.02, sigma=2.0)
    lj.r_cut[('H', 'H')] = 5.0
    
    # O-H cross interactions
    lj.params[('O', 'H')] = dict(epsilon=0.08, sigma=2.5)
    lj.r_cut[('O', 'H')] = 7.0
    
    # Cavity interactions (very weak)
    lj.params[('L', 'O')] = dict(epsilon=0.002, sigma=2.5)
    lj.params[('L', 'H')] = dict(epsilon=0.001, sigma=2.0)
    lj.params[('L', 'L')] = dict(epsilon=0.0, sigma=1.0)
    lj.r_cut[('L', 'O')] = 5.0
    lj.r_cut[('L', 'H')] = 4.0
    lj.r_cut[('L', 'L')] = 2.0
    
    return lj

def run_structure_tracking_simulation():
    """Run 10 ps cavity MD simulation with structure and dipole tracking."""
    
    print("=" * 80)
    print("CAVITY MD SIMULATION WITH STRUCTURE TRACKING (10 ps)")
    print("=" * 80)
    
    # 1. Initialize simulation
    print("\n1. Setting up simulation...")
    device = hoomd.device.CPU()
    sim = hoomd.Simulation(device=device, seed=42)
    
    # 2. Create system
    print("2. Creating molecular system...")
    snapshot = create_structure_tracking_system()
    sim.create_state_from_snapshot(snapshot)
    
    print(f"   System: {snapshot.particles.N} particles")
    print(f"   Water molecules: 12 (36 atoms)")
    print(f"   Cavity particles: 1")
    print(f"   Box: {snapshot.configuration.box[0]:.1f}³")
    print(f"   Density: {snapshot.particles.N / snapshot.configuration.box[0]**3:.4f} particles/Å³")
    
    # 3. Set up forces
    print("\n3. Setting up forces...")
    
    # Enhanced molecular interactions
    lj = setup_enhanced_interactions(sim)
    
    # Cavity force optimized for longer simulations
    cavity_params = {
        'kvector': [1.0, 0.0, 0.0],
        'couplstr': 0.06,  # Moderate coupling for stable long simulation
        'omegac': 0.08,    # Lower frequency for slower cavity dynamics
        'phmass': 1.0
    }
    
    cavity_force = CavityForce(**cavity_params)
    
    print(f"   Cavity coupling: g = {cavity_params['couplstr']:.3f}")
    print(f"   Cavity frequency: ωc = {cavity_params['omegac']:.3f} Hartree")
    print(f"   Cavity period: ~{2*np.pi/cavity_params['omegac']:.1f} a.u. ≈ {2*np.pi/cavity_params['omegac']*0.024:.1f} fs")
    print(f"   Implementation: {cavity_force.implementation}")
    
    # 4. Set up integrator and thermostats
    print("\n4. Setting up integrator...")
    dt = 0.0005  # 0.5 fs timestep
    
    integrator = hoomd.md.Integrator(dt=dt)
    integrator.forces.append(lj)
    integrator.forces.append(cavity_force)
    
    # NVT for molecular particles (lower temperature for stable structure)
    molecular_filter = hoomd.filter.Type(['O', 'H'])
    target_temp = 1.0  # Lower temperature for clearer structure
    nvt_molecular = hoomd.md.methods.ConstantVolume(
        filter=molecular_filter,
        thermostat=hoomd.md.methods.thermostats.MTTK(kT=target_temp, tau=1.0)
    )
    integrator.methods.append(nvt_molecular)
    
    # NVE for cavity particle
    cavity_filter = hoomd.filter.Type(['L'])
    nve_cavity = hoomd.md.methods.ConstantVolume(filter=cavity_filter)
    integrator.methods.append(nve_cavity)
    
    sim.operations.integrator = integrator
    
    print(f"   Timestep: {dt:.4f} fs")
    print(f"   Molecular thermostat: NVT (kT = {target_temp:.1f}, τ = 1.0)")
    print(f"   Cavity dynamics: NVE")
    
    # 5. Set up structure and dipole tracking
    print("\n5. Setting up analysis trackers...")
    
    # Structure analysis - F(k,t)
    # Multiple k-values for comprehensive analysis
    k_values = [0.5, 1.0, 1.5, 2.0]  # Å⁻¹
    
    structure_trackers = []
    for k_val in k_values:
        tracker = DensityCorrelationTracker(
            molecular_filter=molecular_filter,
            k_vector=[k_val, 0.0, 0.0],
            filename=f"structure_k_{k_val:.1f}_10ps.txt",
            reference_period=50  # Reference every 25 fs
        )
        structure_trackers.append(tracker)
        sim.operations.analyzers.append(tracker)
    
    # Dipole autocorrelation
    dipole_tracker = DipoleAutocorrelation(
        molecular_filter=molecular_filter,
        filename="dipole_autocorr_10ps.txt",
        reference_period=50
    )
    sim.operations.analyzers.append(dipole_tracker)
    
    # Cavity mode tracking
    cavity_tracker = CavityModeTracker(
        cavity_force=cavity_force,
        cavity_filter=cavity_filter,
        filename="cavity_mode_10ps.txt"
    )
    sim.operations.analyzers.append(cavity_tracker)
    
    print(f"   Structure tracking: F(k,t) for k = {k_values} Å⁻¹")
    print(f"   Dipole tracking: autocorrelation function")
    print(f"   Cavity tracking: mode position and velocity")
    print(f"   Reference period: {50 * dt:.3f} ps")
    
    # Set tracker periods
    for tracker in structure_trackers:
        tracker.period = 20  # Every 10 fs
    dipole_tracker.period = 20
    cavity_tracker.period = 20
    
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
    print(f"   Initial dipole magnitude: {np.linalg.norm(molecular_dipole):.3f}")
    
    # 7. Run simulation
    print("\n7. Running 10 ps simulation...")
    print("=" * 80)
    
    total_steps = 20000  # 10 ps at 0.5 fs timestep
    report_interval = 2000  # Report every 1 ps
    
    print(f"Total simulation: {total_steps * dt:.1f} ps ({total_steps} steps)")
    print(f"Analysis logging: every {20 * dt:.3f} ps")
    print(f"Progress reports: every {report_interval * dt:.1f} ps")
    print()
    print("Step     Time(ps)  | Cavity Position        | Dipole Magnitude | Status")
    print("-" * 80)
    
    # Storage for real-time analysis
    times = []
    cavity_positions = []
    dipole_magnitudes = []
    cavity_energies = []
    
    for step in range(0, total_steps + 1, report_interval):
        if step > 0:
            sim.run(report_interval)
        
        # Get current time
        current_time = step * dt
        times.append(current_time)
        
        # Get current state
        with sim.state.cpu_local_snapshot as snap:
            positions = snap.particles.position
            charges = snap.particles.charge
            
            # Calculate current dipole
            current_dipole = np.sum(charges[:-1, np.newaxis] * positions[:-1], axis=0)
            dipole_mag = np.linalg.norm(current_dipole)
            dipole_magnitudes.append(dipole_mag)
            
            # Cavity position
            cavity_pos = positions[-1]
            cavity_positions.append(cavity_pos.copy())
        
        # Get cavity energy
        try:
            cavity_energy = cavity_force.total_cavity_energy
            cavity_energies.append(cavity_energy)
        except:
            cavity_energy = 0.0
            cavity_energies.append(0.0)
        
        # Progress status
        progress = step / total_steps * 100
        status = f"{progress:5.1f}%"
        
        print(f"{step:6d}  {current_time:8.3f}   | [{cavity_pos[0]:6.3f}, {cavity_pos[1]:6.3f}, {cavity_pos[2]:6.3f}] | {dipole_mag:14.3f} | {status}")
    
    print("-" * 80)
    print(f"Simulation completed: {total_steps * dt:.1f} ps")
    
    # 8. Final analysis
    print("\n8. Final structure analysis:")
    
    # Final dipole
    final_dipole = current_dipole
    initial_dipole = molecular_dipole
    
    print(f"   Initial dipole: [{initial_dipole[0]:.3f}, {initial_dipole[1]:.3f}, {initial_dipole[2]:.3f}] (mag: {np.linalg.norm(initial_dipole):.3f})")
    print(f"   Final dipole:   [{final_dipole[0]:.3f}, {final_dipole[1]:.3f}, {final_dipole[2]:.3f}] (mag: {np.linalg.norm(final_dipole):.3f})")
    
    # Cavity displacement
    initial_cavity = np.array([0.0, 0.0, 0.0])
    final_cavity = cavity_pos
    cavity_displacement = np.linalg.norm(final_cavity - initial_cavity)
    
    print(f"   Cavity displacement: {cavity_displacement:.3f} Å")
    print(f"   Average cavity energy: {np.mean(cavity_energies):.6f} ± {np.std(cavity_energies):.6f} Hartree")
    
    # 9. Create comprehensive analysis plots
    print("\n9. Creating analysis plots...")
    
    # Convert times and positions to numpy arrays
    times = np.array(times)
    cavity_positions = np.array(cavity_positions)
    dipole_magnitudes = np.array(dipole_magnitudes)
    cavity_energies = np.array(cavity_energies)
    
    fig = plt.figure(figsize=(15, 12))
    
    # Plot 1: Cavity trajectory (2D projection)
    ax1 = plt.subplot(3, 3, 1)
    plt.plot(cavity_positions[:, 0], cavity_positions[:, 1], 'b-', alpha=0.7, linewidth=1)
    plt.scatter(cavity_positions[0, 0], cavity_positions[0, 1], c='green', s=100, marker='o', label='Start', zorder=5)
    plt.scatter(cavity_positions[-1, 0], cavity_positions[-1, 1], c='red', s=100, marker='s', label='End', zorder=5)
    plt.xlabel('X Position (Å)')
    plt.ylabel('Y Position (Å)')
    plt.title('Cavity Particle Trajectory (XY)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.axis('equal')
    
    # Plot 2: Cavity position vs time
    ax2 = plt.subplot(3, 3, 2)
    plt.plot(times, cavity_positions[:, 0], 'r-', label='X', alpha=0.8)
    plt.plot(times, cavity_positions[:, 1], 'g-', label='Y', alpha=0.8)
    plt.plot(times, cavity_positions[:, 2], 'b-', label='Z', alpha=0.8)
    plt.xlabel('Time (ps)')
    plt.ylabel('Position (Å)')
    plt.title('Cavity Position vs Time')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 3: Dipole magnitude evolution
    ax3 = plt.subplot(3, 3, 3)
    plt.plot(times, dipole_magnitudes, 'purple', linewidth=2)
    plt.xlabel('Time (ps)')
    plt.ylabel('Dipole Magnitude')
    plt.title('Molecular Dipole Evolution')
    plt.grid(True, alpha=0.3)
    
    # Plot 4: Cavity energy evolution
    ax4 = plt.subplot(3, 3, 4)
    plt.plot(times, cavity_energies, 'orange', linewidth=2)
    plt.xlabel('Time (ps)')
    plt.ylabel('Cavity Energy (Hartree)')
    plt.title('Cavity Energy Evolution')
    plt.grid(True, alpha=0.3)
    
    # Plot 5: Position-energy correlation
    ax5 = plt.subplot(3, 3, 5)
    cavity_displacement_mag = np.linalg.norm(cavity_positions, axis=1)
    plt.scatter(cavity_displacement_mag, cavity_energies, alpha=0.6, c=times, cmap='viridis')
    plt.xlabel('Cavity Displacement (Å)')
    plt.ylabel('Cavity Energy (Hartree)')
    plt.title('Position-Energy Correlation')
    plt.grid(True, alpha=0.3)
    cbar = plt.colorbar(ax5.collections[0], ax=ax5)
    cbar.set_label('Time (ps)')
    
    # Plot 6: Dipole-cavity correlation
    ax6 = plt.subplot(3, 3, 6)
    plt.scatter(dipole_magnitudes, cavity_displacement_mag, alpha=0.6, c=times, cmap='plasma')
    plt.xlabel('Dipole Magnitude')
    plt.ylabel('Cavity Displacement (Å)')
    plt.title('Dipole-Cavity Coupling')
    plt.grid(True, alpha=0.3)
    cbar = plt.colorbar(ax6.collections[0], ax=ax6)
    cbar.set_label('Time (ps)')
    
    # Plot 7: Histograms of cavity position
    ax7 = plt.subplot(3, 3, 7)
    plt.hist(cavity_positions[:, 0], bins=20, alpha=0.7, label='X', density=True)
    plt.hist(cavity_positions[:, 1], bins=20, alpha=0.7, label='Y', density=True)
    plt.xlabel('Position (Å)')
    plt.ylabel('Probability Density')
    plt.title('Cavity Position Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Plot 8: Energy histogram
    ax8 = plt.subplot(3, 3, 8)
    plt.hist(cavity_energies, bins=20, alpha=0.7, color='orange', density=True)
    plt.xlabel('Cavity Energy (Hartree)')
    plt.ylabel('Probability Density')
    plt.title('Cavity Energy Distribution')
    plt.grid(True, alpha=0.3)
    
    # Plot 9: System overview
    ax9 = plt.subplot(3, 3, 9)
    # Plot both dipole and cavity displacement on same plot (normalized)
    dipole_norm = dipole_magnitudes / np.max(dipole_magnitudes)
    cavity_norm = cavity_displacement_mag / np.max(cavity_displacement_mag)
    plt.plot(times, dipole_norm, 'purple', label='Dipole (norm)', linewidth=2)
    plt.plot(times, cavity_norm, 'blue', label='Cavity (norm)', linewidth=2)
    plt.xlabel('Time (ps)')
    plt.ylabel('Normalized Magnitude')
    plt.title('System Dynamics Overview')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('structure_tracking_10ps_analysis.png', dpi=300, bbox_inches='tight')
    print(f"   Saved comprehensive analysis: structure_tracking_10ps_analysis.png")
    
    print(f"\n✅ Structure tracking simulation completed successfully!")
    print(f"   Data files created:")
    for k_val in k_values:
        print(f"     - structure_k_{k_val:.1f}_10ps.txt (F(k,t) analysis)")
    print(f"     - dipole_autocorr_10ps.txt (dipole correlation)")
    print(f"     - cavity_mode_10ps.txt (cavity dynamics)")
    print(f"   This demonstrates advanced structure and correlation analysis over 10 ps.")
    
    return sim, structure_trackers, dipole_tracker, cavity_tracker

if __name__ == "__main__":
    try:
        simulation, struct_trackers, dipole_tracker, cavity_tracker = run_structure_tracking_simulation()
        print(f"\nSuccess! Structure tracking simulation completed.")
        print(f"Check the generated data files for detailed correlation analysis.")
    except Exception as e:
        print(f"\n❌ Error running structure tracking simulation: {e}")
        import traceback
        traceback.print_exc() 