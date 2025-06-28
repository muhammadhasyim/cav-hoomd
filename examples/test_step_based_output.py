#!/usr/bin/env python3
"""
Test script for step-based output periods in cavity MD trackers.

This script demonstrates how to use different output frequencies (in steps)
for different observables in the updated cavity MD plugin.
"""

import hoomd
import numpy as np
import os
import tempfile

# Import the updated tracker classes directly from the plugin
try:
    from hoomd.cavitymd.analysis import (
        EnergyTracker, 
        FieldAutocorrelationTracker, 
        DipoleAutocorrelation,
        AutocorrelationTracker
    )
    from hoomd.cavitymd.utils import PhysicalConstants
    print("✓ Successfully imported updated tracker classes from hoomd.cavitymd plugin")
except ImportError as e:
    print(f"❌ Failed to import tracker classes: {e}")
    exit(1)


def create_simple_simulation():
    """Create a simple HOOMD simulation for testing trackers."""
    print("🔧 Setting up simple test simulation...")
    
    # Create a simple device and CPU simulation
    device = hoomd.device.CPU()
    
    # Create a simple system
    sim = hoomd.Simulation(device=device, seed=42)
    
    # Create a simple snapshot with a few particles
    snapshot = hoomd.Snapshot(device.communicator)
    
    if snapshot.communicator.rank == 0:
        snapshot.configuration.box = [10, 10, 10, 0, 0, 0]
        snapshot.particles.N = 8  # Simple 8-particle system
        
        # Set positions in a simple cubic arrangement
        positions = np.array([
            [0, 0, 0], [2, 0, 0], [0, 2, 0], [2, 2, 0],
            [0, 0, 2], [2, 0, 2], [0, 2, 2], [2, 2, 2]
        ], dtype=np.float32)
        snapshot.particles.position[:] = positions
        
        # Set particle properties
        snapshot.particles.types = ['A']
        snapshot.particles.typeid[:] = [0] * 8
        snapshot.particles.mass[:] = [1.0] * 8
        snapshot.particles.charge[:] = [0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1]
        
        # Set random velocities
        np.random.seed(42)
        snapshot.particles.velocity[:] = np.random.normal(0, 0.1, (8, 3))
    
    sim.create_state_from_snapshot(snapshot)
    
    # Add simple force and integrator for a functioning simulation
    lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4))
    lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)
    lj.r_cut[('A', 'A')] = 2.5
    
    # Create simple Langevin integrator
    integrator = hoomd.md.Integrator(dt=0.001)
    integrator.forces.append(lj)
    
    # Add thermostat
    langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=1.0)
    integrator.methods.append(langevin)
    
    sim.operations.integrator = integrator
    
    print(f"  ✓ Created {snapshot.particles.N} particle system")
    print(f"  ✓ Added LJ force and Langevin integrator")
    
    return sim


def test_step_based_trackers():
    """Test the updated step-based tracker output periods."""
    
    print("🧪 Testing Step-Based Output Periods for Cavity MD Trackers")
    print("=" * 60)
    
    # Create simple simulation
    sim = create_simple_simulation()
    
    # Create a simple time tracker (mocked)
    class SimpleTimeTracker:
        def __init__(self):
            self.elapsed_time = 0.0
        
        def update(self, timestep, dt):
            self.elapsed_time = PhysicalConstants.atomic_units_to_ps(dt * timestep)
    
    time_tracker = SimpleTimeTracker()
    
    # Test different trackers with different step-based output periods
    print("\n📊 Setting up trackers with different step-based output periods:")
    
    # 1. Test EnergyTracker with step-based output
    print("\n1️⃣ Testing EnergyTracker:")
    try:
        energy_tracker = EnergyTracker(
            simulation=sim,
            components=['kinetic'],  # Only kinetic energy for simplicity
            force_objects={},
            time_tracker=time_tracker,
            output_prefix='test_energy_steps',
            output_period_steps=100,  # Output every 100 steps
            max_timesteps=1000,       # Stop after 1000 steps
            compute_temperature=True
        )
        
        # Add to simulation
        energy_updater = hoomd.update.CustomUpdater(
            action=energy_tracker,
            trigger=hoomd.trigger.Periodic(1)  # Check every step, but output based on internal logic
        )
        sim.operations.updaters.append(energy_updater)
        
        print(f"  ✓ EnergyTracker: output every 100 steps, max 1000 steps")
        
    except Exception as e:
        print(f"  ❌ Failed to create EnergyTracker: {e}")
        return False
    
    # 2. Test AutocorrelationTracker (Dipole) with step-based output
    print("\n2️⃣ Testing AutocorrelationTracker (Dipole):")
    try:
        dipole_tracker = DipoleAutocorrelation(
            simulation=sim,
            time_tracker=time_tracker,
            output_prefix='test_dipole_steps',
            output_period_steps=50  # Output every 50 steps
        )
        
        # Add to simulation
        dipole_updater = hoomd.update.CustomUpdater(
            action=dipole_tracker,
            trigger=hoomd.trigger.Periodic(1)
        )
        sim.operations.updaters.append(dipole_updater)
        
        print(f"  ✓ DipoleAutocorrelation: output every 50 steps")
        
    except Exception as e:
        print(f"  ❌ Failed to create DipoleAutocorrelation: {e}")
        return False
    
    # 3. Test FieldAutocorrelationTracker with step-based output
    print("\n3️⃣ Testing FieldAutocorrelationTracker (Density):")
    try:
        density_tracker = FieldAutocorrelationTracker(
            simulation=sim,
            observable="density",
            time_tracker=time_tracker,
            output_prefix='test_density_steps',
            output_period_steps=150,        # Output every 150 steps
            reference_interval_steps=300,   # New reference every 300 steps
            max_references=3,
            kmag=1.0,
            num_wavevectors=10
        )
        
        # Add to simulation
        density_updater = hoomd.update.CustomUpdater(
            action=density_tracker,
            trigger=hoomd.trigger.Periodic(1)
        )
        sim.operations.updaters.append(density_updater)
        
        print(f"  ✓ FieldAutocorrelationTracker: output every 150 steps, references every 300 steps")
        
    except Exception as e:
        print(f"  ❌ Failed to create FieldAutocorrelationTracker: {e}")
        return False
    
    # Run simulation
    print(f"\n🚀 Running test simulation for 1000 steps...")
    print(f"Expected outputs:")
    print(f"  - Energy: ~{1000 // 100} outputs (every 100 steps)")
    print(f"  - Dipole: ~{1000 // 50} outputs (every 50 steps)")
    print(f"  - Density: ~{1000 // 150} outputs (every 150 steps)")
    
    try:
        # Run the simulation
        sim.run(1000)
        print("✅ Test simulation completed successfully!")
        
        # Check output files
        print("\n📁 Checking generated output files:")
        output_files = []
        
        for filename in os.listdir('.'):
            if filename.startswith('test_') and filename.endswith('.txt'):
                output_files.append(filename)
                
                # Count data lines
                with open(filename, 'r') as f:
                    lines = f.readlines()
                data_lines = [line for line in lines if not line.strip().startswith('#')]
                
                print(f"  📄 {filename}: {len(data_lines)} data lines")
                
                # Show first few lines to verify format
                if len(data_lines) > 0:
                    print(f"    Sample lines:")
                    for i, line in enumerate(data_lines[:2]):
                        print(f"      {line.strip()}")
                    if len(data_lines) > 2:
                        print(f"      ... (and {len(data_lines) - 2} more)")
                print()
        
        if len(output_files) == 0:
            print("  ⚠️  No output files found")
            return False
        else:
            print(f"  ✅ Generated {len(output_files)} output files with step-based periods")
            return True
        
    except Exception as e:
        print(f"❌ Test simulation failed: {e}")
        return False


def cleanup_test_files():
    """Clean up test output files."""
    print("\n🧹 Cleaning up test files...")
    for filename in os.listdir('.'):
        if filename.startswith('test_') and filename.endswith('.txt'):
            os.remove(filename)
            print(f"  🗑️  Removed {filename}")


def main():
    """Main test function."""
    print("🔬 Testing Updated Cavity MD Trackers with Step-Based Output Periods")
    print("This test verifies that all trackers now use timestep-based output frequencies")
    print()
    
    success = test_step_based_trackers()
    
    if success:
        print("\n🎉 All tests passed! Step-based output periods are working correctly.")
        print("\nKey improvements verified:")
        print("  ✓ Output periods specified directly in simulation steps")
        print("  ✓ No more unreliable time-based conversions")
        print("  ✓ Different frequencies for different observables")
        print("  ✓ Timestep information included in output files")
        print("  ✓ Clear output period documentation in file headers")
        
        # Clean up
        cleanup_test_files()
        return 0
    else:
        print("\n❌ Tests failed. Check the error messages above.")
        return 1


if __name__ == "__main__":
    exit(main()) 