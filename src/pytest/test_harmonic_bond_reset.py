#!/usr/bin/env python3
"""
Test script for HarmonicBondReset class.

This demonstrates how to use the one-time harmonic bond reset
for diatomic molecules.
"""

import hoomd
import numpy as np
from cavitymd.harmonic_bond_reset import HarmonicBondReset


def create_simple_diatomic_system():
    """
    Create a simple test system with a few diatomic molecules.
    
    Returns
    -------
    snapshot : hoomd.Snapshot
        System snapshot with diatomics
    """
    # Create snapshot
    snap = hoomd.Snapshot()
    snap.communicator = hoomd.communicator.Communicator()
    
    # Simple cubic box
    snap.configuration.box = [10, 10, 10, 0, 0, 0]
    
    # Create 4 diatomic molecules (8 particles)
    N = 8
    snap.particles.N = N
    snap.particles.types = ['A']
    
    # Initialize positions (4 dimers along x-axis)
    positions = []
    for i in range(4):
        x = i * 2.5
        positions.append([x, 0, 0])  # Atom 1
        positions.append([x + 1.2, 0, 0])  # Atom 2 (bond length ~1.2)
    
    snap.particles.position[:] = positions
    snap.particles.typeid[:] = [0] * N
    snap.particles.mass[:] = [1.0] * N
    
    # Random velocities
    np.random.seed(123)
    snap.particles.velocity[:] = np.random.randn(N, 3) * 0.1
    
    # Create bonds
    snap.bonds.N = 4
    snap.bonds.types = ['bond']
    snap.bonds.typeid[:] = [0, 0, 0, 0]
    snap.bonds.group[:] = [[0, 1], [2, 3], [4, 5], [6, 7]]
    
    return snap


def main():
    """Run test of harmonic bond reset."""
    
    print("="*60)
    print("HARMONIC BOND RESET TEST")
    print("="*60)
    
    # Create simulation
    device = hoomd.device.CPU()
    sim = hoomd.Simulation(device=device, seed=42)
    
    # Create initial state
    snap = create_simple_diatomic_system()
    sim.create_state_from_snapshot(snap)
    
    # Define bond potential
    harmonic = hoomd.md.bond.Harmonic()
    harmonic.params['bond'] = dict(k=10000.0, r0=1.0)
    
    # Integrator
    integrator = hoomd.md.Integrator(dt=0.001)
    nvt = hoomd.md.methods.NVT(
        filter=hoomd.filter.All(),
        kT=1.0,
        tau=1.0
    )
    integrator.forces.append(harmonic)
    integrator.methods.append(nvt)
    sim.operations.integrator = integrator
    
    # Create harmonic bond reset action
    print("\nCreating HarmonicBondReset...")
    bond_reset = HarmonicBondReset(
        bond_type='bond',
        spring_constant=10000.0,
        equilibrium_length=1.0,
        temperature=2.0,  # Higher temperature for visibility
        kB=1.0,
        seed=42
    )
    print(f"  {bond_reset}")
    
    # Add as custom updater (check every 100 steps)
    reset_updater = hoomd.update.CustomUpdater(
        action=bond_reset,
        trigger=hoomd.trigger.Periodic(100)
    )
    sim.operations.updaters.append(reset_updater)
    
    # Run initial equilibration
    print("\n" + "="*60)
    print("INITIAL EQUILIBRATION (500 steps)")
    print("="*60)
    sim.run(500)
    
    # Get initial bond statistics
    snap_before = sim.state.get_snapshot()
    if snap_before.communicator.rank == 0:
        bond_lengths_before = []
        for i in range(snap_before.bonds.N):
            idx1, idx2 = snap_before.bonds.group[i]
            r1 = snap_before.particles.position[idx1]
            r2 = snap_before.particles.position[idx2]
            bond_lengths_before.append(np.linalg.norm(r2 - r1))
        
        print(f"\nBond lengths BEFORE reset:")
        print(f"  Mean: {np.mean(bond_lengths_before):.6f}")
        print(f"  Std:  {np.std(bond_lengths_before):.6f}")
        print(f"  Min:  {np.min(bond_lengths_before):.6f}")
        print(f"  Max:  {np.max(bond_lengths_before):.6f}")
    
    # TRIGGER THE RESET
    print("\n" + "="*60)
    print("TRIGGERING RESET")
    print("="*60)
    bond_reset.enabled = True
    print(f"bond_reset.enabled = {bond_reset.enabled}")
    print(f"bond_reset.has_reset = {bond_reset.has_reset}")
    
    # Run a few steps to allow reset to happen
    sim.run(200)
    
    # Get statistics after reset
    snap_after = sim.state.get_snapshot()
    if snap_after.communicator.rank == 0:
        bond_lengths_after = []
        for i in range(snap_after.bonds.N):
            idx1, idx2 = snap_after.bonds.group[i]
            r1 = snap_after.particles.position[idx1]
            r2 = snap_after.particles.position[idx2]
            bond_lengths_after.append(np.linalg.norm(r2 - r1))
        
        print(f"\nBond lengths AFTER reset:")
        print(f"  Mean: {np.mean(bond_lengths_after):.6f}")
        print(f"  Std:  {np.std(bond_lengths_after):.6f}")
        print(f"  Min:  {np.min(bond_lengths_after):.6f}")
        print(f"  Max:  {np.max(bond_lengths_after):.6f}")
        
        # Expected thermal distribution width
        sigma_x_expected = np.sqrt(bond_reset.kB * bond_reset.T / bond_reset.K)
        print(f"\nExpected thermal width σ_x = √(k_B T / K) = {sigma_x_expected:.6f}")
    
    print(f"\nFinal state:")
    print(f"  bond_reset.enabled = {bond_reset.enabled}")
    print(f"  bond_reset.has_reset = {bond_reset.has_reset}")
    
    # Try to trigger again (should not reset)
    print("\n" + "="*60)
    print("ATTEMPTING SECOND RESET (should be ignored)")
    print("="*60)
    bond_reset.enabled = True
    sim.run(100)
    
    print("\n" + "="*60)
    print("TEST COMPLETED")
    print("="*60)


if __name__ == '__main__':
    main()

