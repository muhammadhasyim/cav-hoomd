#!/usr/bin/env python3
"""
Example 2: Force Validation - Python vs C++ CPU vs CUDA GPU

This example validates that all three implementations of the cavity force
produce identical results:
1. Pure Python implementation (CavityForcePython)
2. C++ CPU implementation 
3. CUDA GPU implementation

Key features:
- Creates identical test systems for each implementation
- Compares forces on each particle
- Compares energy components
- Reports numerical differences and validates they're within tolerance
- Tests both simple and complex molecular configurations
"""

import hoomd
import hoomd.md
import gsd.hoomd
import numpy as np
import matplotlib.pyplot as plt
from hoomd.cavitymd import CavityForce

def create_validation_system(complexity="simple"):
    """Create test system for force validation.
    
    Args:
        complexity: "simple" or "complex" - determines system size
    """
    snapshot = gsd.hoomd.Frame()
    
    if complexity == "simple":
        # Simple system: 2 charged particles + 1 cavity
        L = 10.0
        snapshot.configuration.box = [L, L, L, 0, 0, 0]
        
        snapshot.particles.N = 3
        snapshot.particles.types = ['A', 'L']
        snapshot.particles.typeid = [0, 0, 1]  # A, A, L
        
        positions = [
            [2.0, 0.0, 0.0],   # Positive charge
            [-2.0, 0.0, 0.0],  # Negative charge
            [0.5, 0.5, 0.0]    # Cavity particle (off-center)
        ]
        
        charges = [1.0, -1.0, 0.0]
        masses = [1.0, 1.0, 1.0]
        
    else:  # complex
        # Complex system: 4 water molecules + 1 cavity
        L = 20.0
        snapshot.configuration.box = [L, L, L, 0, 0, 0]
        
        snapshot.particles.N = 13  # 4 waters (12 atoms) + 1 cavity
        snapshot.particles.types = ['O', 'H', 'L']
        
        # 4 water molecules arranged in a square
        positions = []
        typeids = []
        charges = []
        masses = []
        
        water_positions = [
            [3.0, 3.0, 0.0],   # Water 1 center
            [-3.0, 3.0, 0.0],  # Water 2 center
            [3.0, -3.0, 0.0],  # Water 3 center
            [-3.0, -3.0, 0.0]  # Water 4 center
        ]
        
        for i, water_center in enumerate(water_positions):
            # Oxygen
            positions.append(water_center)
            typeids.append(0)  # O
            charges.append(-0.834)
            masses.append(15.999)
            
            # Hydrogen 1
            positions.append([water_center[0] + 0.96, water_center[1], water_center[2]])
            typeids.append(1)  # H
            charges.append(0.417)
            masses.append(1.008)
            
            # Hydrogen 2
            positions.append([water_center[0] - 0.24, water_center[1] + 0.93, water_center[2]])
            typeids.append(1)  # H
            charges.append(0.417)
            masses.append(1.008)
        
        # Cavity particle at center
        positions.append([0.0, 0.0, 0.0])
        typeids.append(2)  # L
        charges.append(0.0)
        masses.append(1.0)
        
        snapshot.particles.typeid = typeids
        snapshot.particles.charge = charges
        snapshot.particles.mass = masses
    
    snapshot.particles.position = positions
    snapshot.particles.diameter = [1.0] * snapshot.particles.N
    
    return snapshot

def run_single_validation(snapshot, implementation, device_type="CPU"):
    """Run cavity force calculation with a specific implementation."""
    
    # Create simulation
    if device_type == "GPU":
        try:
            device = hoomd.device.GPU()
        except:
            print(f"   GPU not available, skipping GPU test")
            return None
    else:
        device = hoomd.device.CPU()
    
    sim = hoomd.Simulation(device=device, seed=42)
    sim.create_state_from_snapshot(snapshot)
    
    # Cavity parameters
    cavity_params = {
        'kvector': [1.0, 0.0, 0.0],
        'couplstr': 0.1,
        'omegac': 0.05,
        'phmass': 1.0,
        'force_python': (implementation == "python")
    }
    
    cavity_force = CavityForce(**cavity_params)
    
    # Simple integrator (just need to compute forces once)
    integrator = hoomd.md.Integrator(dt=0.001)
    integrator.forces.append(cavity_force)
    
    # NVE for all particles
    nve = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
    integrator.methods.append(nve)
    
    sim.operations.integrator = integrator
    
    # Run one step to compute forces
    sim.run(1)
    
    # Extract forces and energies
    with sim.state.cpu_local_snapshot as snap:
        forces = snap.particles.net_force.copy()
        positions = snap.particles.position.copy()
        charges = snap.particles.charge.copy()
    
    # Get energy components
    energy_components = {
        'harmonic': cavity_force.harmonic_energy,
        'coupling': cavity_force.coupling_energy,
        'dipole_self': cavity_force.dipole_self_energy,
        'total': cavity_force.total_cavity_energy
    }
    
    return {
        'implementation': cavity_force.implementation,
        'forces': forces,
        'positions': positions,
        'charges': charges,
        'energies': energy_components,
        'device': device_type
    }

def compare_results(results, tolerance=1e-10):
    """Compare results from different implementations."""
    
    implementations = list(results.keys())
    print(f"\n📊 COMPARISON RESULTS")
    print(f"Implementations: {implementations}")
    print(f"Tolerance: {tolerance:.0e}")
    print("=" * 80)
    
    # Compare forces
    print(f"\n🔍 FORCE COMPARISON:")
    reference = results[implementations[0]]
    
    for i, impl in enumerate(implementations[1:], 1):
        if results[impl] is None:
            continue
            
        ref_forces = reference['forces']
        test_forces = results[impl]['forces']
        
        # Calculate differences
        force_diff = test_forces - ref_forces
        max_diff = np.max(np.abs(force_diff))
        rms_diff = np.sqrt(np.mean(force_diff**2))
        
        print(f"\n{reference['implementation']} vs {results[impl]['implementation']}:")
        print(f"  Maximum force difference: {max_diff:.2e}")
        print(f"  RMS force difference:     {rms_diff:.2e}")
        
        if max_diff < tolerance:
            print(f"  ✅ Forces match within tolerance")
        else:
            print(f"  ❌ Forces differ by more than tolerance")
            
            # Show detailed per-particle differences
            print(f"  Per-particle force differences:")
            for j in range(len(force_diff)):
                fx, fy, fz = force_diff[j]
                magnitude = np.sqrt(fx**2 + fy**2 + fz**2)
                print(f"    Particle {j}: ΔF = [{fx:.2e}, {fy:.2e}, {fz:.2e}], |ΔF| = {magnitude:.2e}")
    
    # Compare energies
    print(f"\n🔍 ENERGY COMPARISON:")
    for i, impl in enumerate(implementations[1:], 1):
        if results[impl] is None:
            continue
            
        ref_energies = reference['energies']
        test_energies = results[impl]['energies']
        
        print(f"\n{reference['implementation']} vs {results[impl]['implementation']}:")
        for component in ['harmonic', 'coupling', 'dipole_self', 'total']:
            ref_val = ref_energies[component]
            test_val = test_energies[component]
            diff = abs(test_val - ref_val)
            
            print(f"  {component:12s}: {ref_val:.8e} vs {test_val:.8e}, diff = {diff:.2e}")
            
            if diff < tolerance:
                print(f"                   ✅ Match")
            else:
                print(f"                   ❌ Differ")
    
    # Overall assessment
    print(f"\n📋 OVERALL ASSESSMENT:")
    all_match = True
    
    for i, impl in enumerate(implementations[1:], 1):
        if results[impl] is None:
            continue
            
        # Check forces
        force_diff = results[impl]['forces'] - reference['forces']
        max_force_diff = np.max(np.abs(force_diff))
        
        # Check energies
        max_energy_diff = 0
        for component in ['harmonic', 'coupling', 'dipole_self', 'total']:
            diff = abs(results[impl]['energies'][component] - reference['energies'][component])
            max_energy_diff = max(max_energy_diff, diff)
        
        if max_force_diff < tolerance and max_energy_diff < tolerance:
            print(f"  ✅ {results[impl]['implementation']} matches {reference['implementation']}")
        else:
            print(f"  ❌ {results[impl]['implementation']} differs from {reference['implementation']}")
            all_match = False
    
    return all_match

def plot_force_comparison(results):
    """Create plots comparing forces from different implementations."""
    
    implementations = [k for k, v in results.items() if v is not None]
    if len(implementations) < 2:
        print("Not enough implementations for plotting")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Force magnitude comparison
    ax = axes[0, 0]
    for impl in implementations:
        forces = results[impl]['forces']
        force_mags = np.sqrt(np.sum(forces[:, :3]**2, axis=1))  # Only x,y,z components
        ax.plot(force_mags, 'o-', label=results[impl]['implementation'], alpha=0.7)
    
    ax.set_xlabel('Particle Index')
    ax.set_ylabel('Force Magnitude')
    ax.set_title('Force Magnitude Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Force component comparison (x-component)
    ax = axes[0, 1]
    for impl in implementations:
        forces = results[impl]['forces']
        ax.plot(forces[:, 0], 'o-', label=results[impl]['implementation'], alpha=0.7)
    
    ax.set_xlabel('Particle Index')
    ax.set_ylabel('Force X-Component')
    ax.set_title('X-Force Component Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Energy component comparison
    ax = axes[1, 0]
    components = ['harmonic', 'coupling', 'dipole_self', 'total']
    x_pos = np.arange(len(components))
    
    for i, impl in enumerate(implementations):
        energies = [results[impl]['energies'][comp] for comp in components]
        ax.bar(x_pos + i*0.25, energies, 0.25, label=results[impl]['implementation'], alpha=0.7)
    
    ax.set_xlabel('Energy Component')
    ax.set_ylabel('Energy')
    ax.set_title('Energy Component Comparison')
    ax.set_xticks(x_pos + 0.125)
    ax.set_xticklabels(components)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Force difference heatmap (if we have reference)
    ax = axes[1, 1]
    if len(implementations) >= 2:
        ref_forces = results[implementations[0]]['forces']
        test_forces = results[implementations[1]]['forces']
        force_diff = test_forces - ref_forces
        
        im = ax.imshow(force_diff[:, :3].T, aspect='auto', cmap='RdBu_r')
        ax.set_xlabel('Particle Index')
        ax.set_ylabel('Force Component (x, y, z)')
        ax.set_title(f'Force Differences\n({implementations[1]} - {implementations[0]})')
        plt.colorbar(im, ax=ax)
    
    plt.tight_layout()
    plt.savefig('force_validation_comparison.png', dpi=300, bbox_inches='tight')
    print(f"   Saved force comparison plot: force_validation_comparison.png")

def run_force_validation():
    """Run complete force validation."""
    
    print("=" * 80)
    print("CAVITY FORCE VALIDATION")
    print("Testing Python vs C++ CPU vs CUDA GPU implementations")
    print("=" * 80)
    
    # Test both simple and complex systems
    for complexity in ["simple", "complex"]:
        print(f"\n🧪 TESTING {complexity.upper()} SYSTEM")
        print("=" * 40)
        
        # Create test system
        snapshot = create_validation_system(complexity)
        print(f"Created {complexity} system:")
        print(f"  Particles: {snapshot.particles.N}")
        print(f"  Types: {snapshot.particles.types}")
        print(f"  Box: {snapshot.configuration.box[0]:.1f}³")
        
        # Test all implementations
        results = {}
        
        # Python implementation
        print(f"\n1. Testing Python implementation...")
        results['python'] = run_single_validation(snapshot, "python", "CPU")
        if results['python']:
            print(f"   ✅ Python implementation: {results['python']['implementation']}")
        
        # C++ CPU implementation
        print(f"\n2. Testing C++ CPU implementation...")
        results['cpp'] = run_single_validation(snapshot, "cpp", "CPU")
        if results['cpp']:
            print(f"   ✅ C++ implementation: {results['cpp']['implementation']}")
        
        # CUDA GPU implementation
        print(f"\n3. Testing CUDA GPU implementation...")
        results['cuda'] = run_single_validation(snapshot, "cpp", "GPU")
        if results['cuda']:
            print(f"   ✅ CUDA implementation: {results['cuda']['implementation']}")
        
        # Compare results
        valid_results = {k: v for k, v in results.items() if v is not None}
        if len(valid_results) >= 2:
            all_match = compare_results(valid_results)
            
            # Create plots
            plot_force_comparison(valid_results)
            
            if all_match:
                print(f"\n🎉 SUCCESS: All implementations match for {complexity} system!")
            else:
                print(f"\n⚠️  WARNING: Implementations differ for {complexity} system!")
        else:
            print(f"\n❌ Not enough implementations available for comparison")
        
        print("\n" + "=" * 40)

if __name__ == "__main__":
    try:
        run_force_validation()
        print(f"\n✅ Force validation completed!")
    except Exception as e:
        print(f"\n❌ Error during force validation: {e}")
        import traceback
        traceback.print_exc() 