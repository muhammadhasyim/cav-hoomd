#!/usr/bin/env python3
"""
FDR Temperature Estimator Example

This example demonstrates how to use the fluctuation-dissipation ratio (FDR)
based effective temperature estimator for cavity-coupled molecular dynamics.

The example shows:
1. Setting up a basic molecular system
2. Configuring cavity coupling
3. Integrating FDR temperature monitoring
4. Running a non-equilibrium simulation with temperature tracking
5. Analyzing the results

Scientific Context:
- Demonstrates mode-specific effective temperature measurement
- Shows non-equilibrium thermalization dynamics
- Validates FDR theory in practice

Author: Computational Chemistry Team
Date: 2025
"""

import numpy as np
import hoomd
import matplotlib.pyplot as plt
from pathlib import Path
import logging

# Import CavityMD modules
import sys
sys.path.insert(0, '../src')

from cavitymd.fdr_integration import (
    FDRTemperatureMonitor, 
    CavityFDRAnalyzer,
    create_cavity_fdr_monitor
)
from cavitymd.simulation import CavityMDSimulation
from cavitymd.forces import CavityForce


def setup_logging():
    """Setup logging for the example."""
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('fdr_example.log'),
            logging.StreamHandler()
        ]
    )
    
    return logging.getLogger(__name__)


def create_molecular_system():
    """Create a simple molecular system for demonstration."""
    
    # Simple harmonic oscillator system (like diatomic molecule)
    # This provides a well-defined vibrational mode for testing
    
    N_particles = 2
    box_size = 10.0
    
    # Create simulation snapshot
    snapshot = hoomd.Snapshot()
    snapshot.particles.N = N_particles
    snapshot.particles.types = ['A']
    snapshot.configuration.box = [box_size, box_size, box_size, 0, 0, 0]
    
    # Position particles for harmonic oscillator
    positions = np.array([
        [-0.5, 0.0, 0.0],  # Particle 1
        [0.5, 0.0, 0.0]    # Particle 2 
    ])
    snapshot.particles.position[:] = positions
    
    # Set masses and charges
    masses = np.array([1.0, 1.0])  # Equal masses
    charges = np.array([1.0, -1.0])  # Oppositely charged
    snapshot.particles.mass[:] = masses
    snapshot.particles.charge[:] = charges
    
    # Velocities (thermal distribution)
    T_initial = 400.0  # K (start hot)
    k_B = 8.314e-3  # kJ/mol/K, approximate for simplicity
    v_thermal = np.sqrt(k_B * T_initial / masses[0])
    
    velocities = np.random.normal(0, v_thermal, (N_particles, 3))
    snapshot.particles.velocity[:] = velocities
    
    return snapshot, masses, charges


def setup_cavity_simulation(snapshot, cavity_freq_cm=2000.0):
    """Setup cavity-coupled molecular dynamics simulation."""
    
    # Initialize HOOMD simulation
    cpu = hoomd.device.CPU()
    sim = hoomd.Simulation(device=cpu, seed=42)
    sim.create_state_from_snapshot(snapshot)
    
    # Harmonic bond between particles (spring constant tuned to cavity frequency)
    bond_force = hoomd.md.bond.Harmonic()
    bond_force.params['A-A'] = dict(k=100.0, r0=1.0)  # Spring constant and equilibrium distance
    
    # Add bonds
    sim.state.bonds.types = ['A-A']
    bond_snapshot = sim.state.get_snapshot()
    bond_snapshot.bonds.N = 1
    bond_snapshot.bonds.group[0] = [0, 1]  # Bond between particles 0 and 1
    bond_snapshot.bonds.typeid[0] = 0
    sim.state.set_snapshot(bond_snapshot)
    
    # Lennard-Jones interactions (weak, to avoid complications)
    lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4))
    lj.params[('A', 'A')] = dict(epsilon=0.1, sigma=1.0)
    lj.r_cut[('A', 'A')] = 2.5
    
    # Cavity coupling force
    cavity_force = CavityForce(
        coupling_strength=1e-3,  # Weak coupling
        cavity_frequency_cm=cavity_freq_cm,
        cavity_mode_volume=1.0
    )
    
    # Thermostat (Langevin for bath coupling)
    langevin = hoomd.md.methods.Langevin(
        filter=hoomd.filter.All(), 
        kT=300.0 * 8.314e-3,  # Bath temperature in simulation units
        gamma=1.0  # Friction coefficient
    )
    
    # Integrator
    dt = 0.001  # ps
    integrator = hoomd.md.Integrator(dt=dt, methods=[langevin], forces=[bond_force, lj])
    sim.operations.integrator = integrator
    
    return sim, dt


def run_fdr_temperature_example():
    """Run the complete FDR temperature example."""
    
    logger = setup_logging()
    logger.info("Starting FDR Temperature Estimator Example")
    
    # 1. Create molecular system
    logger.info("Creating molecular system...")
    snapshot, masses, charges = create_molecular_system()
    
    # 2. Setup cavity simulation
    cavity_freq_cm = 2000.0  # cm^-1 - typical molecular vibration
    logger.info(f"Setting up cavity simulation (ν = {cavity_freq_cm} cm^-1)...")
    sim, dt = setup_cavity_simulation(snapshot, cavity_freq_cm)
    
    # 3. Create FDR temperature monitor
    logger.info("Initializing FDR temperature monitor...")
    fdr_monitor = create_cavity_fdr_monitor(
        cavity_frequency_cm=cavity_freq_cm,
        simulation_dt=dt,
        observable_type='dipole',
        axis='x',  # Oscillation along x-axis
        output_file='fdr_temperature_trajectory.dat',
        log_interval=500
    )
    
    # 4. Attach monitor to simulation
    fdr_monitor.attach_to_simulation(sim)
    
    # 5. Equilibration and calibration
    logger.info("Running equilibration for FDR calibration...")
    T_bath = 300.0  # K
    fdr_monitor.calibrate_equilibrium(T_cal=T_bath, n_steps=5000)
    
    # 6. Production run with FDR monitoring
    logger.info("Starting production run with FDR monitoring...")
    
    n_production_steps = 20000
    update_interval = 10  # Update FDR every 10 steps
    
    T_eff_trajectory = []
    timesteps = []
    
    for step in range(n_production_steps):
        # Run simulation
        sim.run(update_interval)
        
        # Update FDR temperature
        T_eff, diagnostics = fdr_monitor.update()
        
        # Store results
        T_eff_trajectory.append(T_eff)
        timesteps.append(sim.timestep)
        
        # Progress reporting
        if (step + 1) % 1000 == 0:
            recent_T = fdr_monitor.get_recent_temperature(100)
            progress = (step + 1) / n_production_steps * 100
            logger.info(f"Progress: {progress:.1f}% | Recent T_eff: {recent_T:.1f} K")
            
    # 7. Analysis
    logger.info("Analyzing results...")
    
    analyzer = CavityFDRAnalyzer(fdr_monitor)
    
    # Temperature statistics
    temp_stats = fdr_monitor.get_temperature_statistics()
    logger.info(f"Temperature Statistics:")
    logger.info(f"  Mean: {temp_stats.get('mean', 0):.1f} ± {temp_stats.get('std', 0):.1f} K")
    logger.info(f"  Range: {temp_stats.get('min', 0):.1f} - {temp_stats.get('max', 0):.1f} K")
    
    # Thermalization analysis
    therm_results = analyzer.analyze_thermalization_dynamics()
    if therm_results:
        logger.info(f"Thermalization Analysis:")
        logger.info(f"  Equilibrium T: {therm_results.get('equilibrium_temperature', 0):.1f} K")
        logger.info(f"  Relaxation time: {therm_results.get('relaxation_time', 0):.2f} ps")
        
    # Comparison with equipartition
    comparison = analyzer.compare_with_equipartition(T_bath)
    if comparison:
        logger.info(f"FDR vs Equipartition:")
        logger.info(f"  FDR T_eff: {comparison['fdr_temperature']:.1f} K")
        logger.info(f"  Bath T: {comparison['bath_temperature']:.1f} K") 
        logger.info(f"  Deviation: {comparison['relative_deviation']*100:.1f}%")
        logger.info(f"  Non-equilibrium: {comparison['is_non_equilibrium']}")
        
    # 8. Plotting
    logger.info("Creating plots...")
    create_analysis_plots(timesteps, T_eff_trajectory, dt, T_bath, fdr_monitor)
    
    logger.info("FDR Temperature Example completed successfully!")
    
    return fdr_monitor, analyzer


def create_analysis_plots(timesteps, T_eff_trajectory, dt, T_bath, fdr_monitor):
    """Create analysis plots for the FDR temperature example."""
    
    # Convert to time in ps
    times_ps = np.array(timesteps) * dt
    T_eff_array = np.array(T_eff_trajectory)
    
    # Remove NaN values for plotting
    valid_mask = ~np.isnan(T_eff_array)
    times_valid = times_ps[valid_mask]
    T_eff_valid = T_eff_array[valid_mask]
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('FDR Effective Temperature Analysis', fontsize=14)
    
    # 1. Temperature trajectory
    ax1 = axes[0, 0]
    ax1.plot(times_valid, T_eff_valid, 'b-', alpha=0.7, linewidth=1)
    ax1.axhline(T_bath, color='r', linestyle='--', label=f'Bath T = {T_bath} K')
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('T_eff (K)')
    ax1.set_title('FDR Temperature Trajectory')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Temperature distribution
    ax2 = axes[0, 1]
    ax2.hist(T_eff_valid, bins=30, alpha=0.7, density=True, color='skyblue', edgecolor='black')
    ax2.axvline(T_bath, color='r', linestyle='--', label=f'Bath T = {T_bath} K')
    ax2.axvline(np.mean(T_eff_valid), color='g', linestyle='-', label=f'Mean = {np.mean(T_eff_valid):.1f} K')
    ax2.set_xlabel('T_eff (K)')
    ax2.set_ylabel('Probability Density')
    ax2.set_title('Temperature Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Moving average
    ax3 = axes[1, 0]
    window = min(100, len(T_eff_valid) // 10)
    if window > 1:
        moving_avg = np.convolve(T_eff_valid, np.ones(window)/window, mode='valid')
        moving_times = times_valid[window//2:len(moving_avg)+window//2]
        ax3.plot(moving_times, moving_avg, 'g-', linewidth=2, label=f'Moving average (n={window})')
        ax3.plot(times_valid, T_eff_valid, 'b-', alpha=0.3, linewidth=0.5, label='Raw data')
    else:
        ax3.plot(times_valid, T_eff_valid, 'b-', linewidth=1)
        
    ax3.axhline(T_bath, color='r', linestyle='--', label=f'Bath T = {T_bath} K')
    ax3.set_xlabel('Time (ps)')
    ax3.set_ylabel('T_eff (K)')
    ax3.set_title('Smoothed Temperature Evolution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. Diagnostics from recent data
    ax4 = axes[1, 1]
    if len(fdr_monitor.diagnostics_history) > 0:
        recent_diag = fdr_monitor.diagnostics_history[-min(1000, len(fdr_monitor.diagnostics_history)):]
        snr_values = [d.snr for d in recent_diag]
        gamma_values = [d.gamma for d in recent_diag]
        
        ax4_twin = ax4.twinx()
        
        line1 = ax4.plot(snr_values, 'b-', label='SNR')
        ax4.set_ylabel('Signal-to-Noise Ratio', color='b')
        ax4.tick_params(axis='y', labelcolor='b')
        
        line2 = ax4_twin.plot(gamma_values, 'r-', label='γ (1/ps)')
        ax4_twin.set_ylabel('Damping γ (1/ps)', color='r')
        ax4_twin.tick_params(axis='y', labelcolor='r')
        
        ax4.set_xlabel('Recent Steps')
        ax4.set_title('FDR Diagnostics')
        
        # Combined legend
        lines = line1 + line2
        labels = [l.get_label() for l in lines]
        ax4.legend(lines, labels, loc='upper left')
    else:
        ax4.text(0.5, 0.5, 'No diagnostics available', 
                ha='center', va='center', transform=ax4.transAxes)
        ax4.set_title('FDR Diagnostics')
    
    plt.tight_layout()
    plt.savefig('fdr_temperature_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig('fdr_temperature_analysis.pdf', bbox_inches='tight')
    plt.show()
    
    print(f"Analysis plots saved as fdr_temperature_analysis.png/pdf")


if __name__ == "__main__":
    try:
        fdr_monitor, analyzer = run_fdr_temperature_example()
        print("\n✅ FDR Temperature Example completed successfully!")
        print("📁 Check output files:")
        print("   - fdr_temperature_trajectory.dat")
        print("   - fdr_temperature_analysis.png")
        print("   - fdr_example.log")
        
    except Exception as e:
        print(f"\n❌ Example failed with error: {e}")
        import traceback
        traceback.print_exc()
