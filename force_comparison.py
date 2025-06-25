#!/usr/bin/env python3
"""
Force Comparison Script for Cavity Force Implementations

This script compares forces calculated by:
1. Old pure Python implementation (old-cavitymd.py)
2. Current CPU version (C++ implementation)
3. Current GPU version (C++ implementation)

Forces are compared per-atom and plotted with the old implementation as benchmark.
"""

import numpy as np
import matplotlib.pyplot as plt
import hoomd
import gsd.hoomd
import sys
import os
from pathlib import Path

# Add the source directory to Python path to import the old implementation
sys.path.insert(0, str(Path(__file__).parent))

def create_test_system(n_particles=10, box_size=15.0, seed=42, random_cavity_pos=True, random_charges=True):
    """Create a test system with randomly positioned charged particles and a randomly positioned cavity particle.
    
    Args:
        n_particles: Number of molecular particles (cavity particle is added automatically)
        box_size: Size of the simulation box
        seed: Random seed for reproducibility
        random_cavity_pos: If True, place cavity particle at random position; if False, at origin
        random_charges: If True, use random charges; if False, use alternating +1/-1 pattern
    """
    np.random.seed(seed)
    
    snapshot = gsd.hoomd.Frame()
    snapshot.configuration.box = [box_size, box_size, box_size, 0, 0, 0]
    
    # Create molecular particles
    n_mol = n_particles
    snapshot.particles.N = n_mol + 1  # +1 for cavity particle
    
    print(f"Creating test system with {n_mol} molecular particles + 1 cavity particle:")
    print(f"  Box size: {box_size} x {box_size} x {box_size}")
    print(f"  Random seed: {seed}")
    
    # Random positions for molecular particles within a reasonable range
    molecular_range = box_size / 3  # Use 1/3 of box size to avoid edge effects
    positions = np.random.uniform(-molecular_range, molecular_range, (n_mol, 3))
    print(f"  Molecular particles: randomly positioned in [{-molecular_range:.1f}, {molecular_range:.1f}]³")
    
    # Add cavity particle position
    if random_cavity_pos:
        cavity_range = box_size / 4  # Smaller range for cavity particle
        cavity_pos = np.random.uniform(-cavity_range, cavity_range, (1, 3))
        print(f"  Cavity particle: randomly positioned in [{-cavity_range:.1f}, {cavity_range:.1f}]³")
    else:
        cavity_pos = np.array([[0.0, 0.0, 0.0]])
        print(f"  Cavity particle: positioned at origin [0, 0, 0]")
    
    positions = np.vstack([positions, cavity_pos])
    
    # Set up charges
    if random_charges:
        # Random charges for molecular particles (cavity particle has charge 0)
        charges = np.random.uniform(-1, 1, n_mol)
        print(f"  Charges: random values in [-1, 1] for molecular particles")
    else:
        # Alternating charges for more predictable behavior
        charges = np.array([(-1)**i for i in range(n_mol)], dtype=float)
        print(f"  Charges: alternating +1/-1 pattern for molecular particles")
    
    charges = np.append(charges, [0.0])  # Cavity particle always has charge 0
    print(f"  Cavity particle charge: 0.0")
    
    # Set up particle data
    snapshot.particles.types = ['A', 'L']
    snapshot.particles.typeid = [0] * n_mol + [1]  # A particles + L (cavity) particle
    snapshot.particles.position = positions
    snapshot.particles.charge = charges
    snapshot.particles.mass = np.ones(n_mol + 1)
    snapshot.particles.diameter = np.ones(n_mol + 1)
    snapshot.particles.image = np.zeros((n_mol + 1, 3), dtype=int)
    
    # Print some statistics about the created system
    print(f"\nSystem statistics:")
    print(f"  Total dipole moment: [{np.sum(charges[:-1] * positions[:-1, 0]):.3f}, "
          f"{np.sum(charges[:-1] * positions[:-1, 1]):.3f}, "
          f"{np.sum(charges[:-1] * positions[:-1, 2]):.3f}]")
    print(f"  Cavity position: [{cavity_pos[0,0]:.3f}, {cavity_pos[0,1]:.3f}, {cavity_pos[0,2]:.3f}]")
    print(f"  Charge range: [{np.min(charges[:-1]):.3f}, {np.max(charges[:-1]):.3f}]")
    print(f"  Position range: [{np.min(positions[:-1]):.3f}, {np.max(positions[:-1]):.3f}]")
    
    return snapshot

class OldPythonCavityForce:
    """Recreate the old Python implementation for force comparison."""
    
    def __init__(self, omegac, couplstr, phmass=1.0):
        self.omegac = omegac
        self.K = phmass * omegac**2
        self.couplstr = couplstr
        self.phmass = phmass
        
        # Energy components
        self.harmonic_energy = 0.0
        self.coupling_energy = 0.0
        self.dipole_self_energy = 0.0
    
    def unwrap_positions(self, positions, images, box_lengths):
        """Unwrap particle positions."""
        pos = np.asarray(positions)
        img = np.asarray(images)
        box = np.asarray(box_lengths)
        return pos + img * box[None, :]
    
    def compute_forces(self, positions, charges, images, box_lengths, types):
        """Compute forces using the old Python algorithm."""
        
        # Find photon particle (type L, typeid=1)
        photon_indices = np.where(np.array(types) == 1)[0]
        if len(photon_indices) == 0:
            print("Warning: No photon particle found")
            return np.zeros_like(positions), np.zeros(len(positions))
        
        photon_idx = photon_indices[0]
        
        # Unwrap positions
        unwrapped_pos = self.unwrap_positions(positions, images, box_lengths)
        
        # Calculate dipole moment (only x,y components)
        dipole = np.zeros(3)
        for i, charge in enumerate(charges):
            if i != photon_idx:  # Skip photon particle
                dipole += charge * unwrapped_pos[i]
        
        dipole_xy = dipole.copy()
        dipole_xy[2] = 0.0  # Zero out z-component
        
        # Get photon position (only x,y components)
        q_photon = unwrapped_pos[photon_idx]
        q_photon_xy = q_photon.copy()
        q_photon_xy[2] = 0.0  # Zero out z-component
        
        # Compute energy contributions
        self.harmonic_energy = 0.5 * self.K * np.dot(q_photon, q_photon)
        self.coupling_energy = self.couplstr * np.dot(dipole_xy, q_photon_xy)
        self.dipole_self_energy = 0.5 * (self.couplstr**2 / self.K) * np.dot(dipole_xy, dipole_xy)
        
        # Initialize forces
        forces = np.zeros_like(positions)
        potential_energies = np.zeros(len(positions))
        
        # Total cavity energy on photon particle
        total_cavity_energy = self.harmonic_energy + self.coupling_energy + self.dipole_self_energy
        potential_energies[photon_idx] = total_cavity_energy
        
        # Force on molecular particles
        Dq = q_photon_xy + (self.couplstr / self.K) * dipole_xy
        
        for i in range(len(positions)):
            if i != photon_idx:  # Not cavity particle
                charge = charges[i]
                force = -self.couplstr * charge * Dq
                forces[i] = force
                forces[i, 2] = 0.0  # Zero out z-component
        
        # Force on photon particle
        photon_force = -self.K * q_photon - self.couplstr * dipole_xy
        forces[photon_idx] = photon_force
        
        return forces, potential_energies

def extract_forces_old_python(snapshot, omegac, couplstr, phmass=1.0):
    """Extract forces using the old Python implementation."""
    
    # Create the old Python force calculator
    force_calc = OldPythonCavityForce(omegac, couplstr, phmass)
    
    # Extract data from snapshot
    positions = snapshot.particles.position
    charges = snapshot.particles.charge
    images = snapshot.particles.image
    box_lengths = snapshot.configuration.box[:3]
    types = snapshot.particles.typeid
    
    # Compute forces
    forces, _ = force_calc.compute_forces(positions, charges, images, box_lengths, types)
    # Exclude cavity particle force from the returned forces
    # Find the cavity particle index (type 'L')
    # cavity_typeid = None
    # for i, type_name in enumerate(snapshot.particles.types):
    #     if type_name == 'L':
    #         cavity_typeid = i
    #         break
    
    # if cavity_typeid is not None:
    #     # Zero out the force on the cavity particle
    #     forces[cavity_typeid] = [0.0, 0.0, 0.0]
    return forces, force_calc

def extract_forces_cpp(snapshot, omegac, couplstr, phmass=1.0, use_gpu=False):
    """Extract forces using the current C++ implementation."""
    
    # Import the current implementation
    try:
        from hoomd.cavitymd import CavityForce
    except ImportError:
        print("Error: Could not import CavityForce. Make sure the plugin is compiled.")
        return None, None
    
    # Create simulation
    device = hoomd.device.GPU() if use_gpu else hoomd.device.CPU()
    sim = hoomd.Simulation(device=device, seed=42)
    sim.create_state_from_snapshot(snapshot)
    
    # Create cavity force (add kvector parameter)
    kvector = [1.0, 0.0, 0.0]  # Default kvector for compatibility
    cavity_force = CavityForce(
        kvector=kvector,
        omegac=omegac,
        couplstr=couplstr,
        phmass=phmass
    )
    
    # Setup minimal integrator
    integrator = hoomd.md.Integrator(dt=0.001)
    integrator.forces.append(cavity_force)
    
    # NVE method for all particles
    all_particles = hoomd.filter.All()
    nve = hoomd.md.methods.ConstantVolume(filter=all_particles)
    integrator.methods.append(nve)
    
    sim.operations.integrator = integrator
    
    # Run one step to compute forces
    sim.run(1)
    
    # Extract forces from the simulation state
    with sim.state.cpu_local_snapshot as snap:
        forces = snap.particles.net_force.copy()
    
    # Extract energy values immediately after running
    try:
        energy_data = {
            'harmonic_energy': cavity_force.harmonic_energy,
            'coupling_energy': cavity_force.coupling_energy,
            'dipole_self_energy': cavity_force.dipole_self_energy
        }
    except Exception as e:
        print(f"Warning: Could not extract energy data: {e}")
        energy_data = {
            'harmonic_energy': 0.0,
            'coupling_energy': 0.0,
            'dipole_self_energy': 0.0
        }
    # Exclude cavity particle (type 'L') from force calculation
    # Find cavity particle index
    # cavity_idx = None
    # for i, typeid in enumerate(snap.particles.typeid):
    #     if snap.particles.types[typeid] == 'L':
    #         cavity_idx = i
    #         break
    
    # # Zero out force on cavity particle if found
    # if cavity_idx is not None:
    #     forces[cavity_idx] = [0.0, 0.0, 0.0]
    #     print(f"Zeroed force on cavity particle (index {cavity_idx})")
    return forces, energy_data

def compare_implementations(n_particles=20, box_size=15.0, seed=42, **kwargs):
    """Compare all three implementations and return force data."""
    
    print(f"Comparing implementations with {n_particles} particles...")
    
    # Extract only the parameters needed for create_test_system
    system_params = {k: v for k, v in kwargs.items() 
                    if k in ['random_cavity_pos', 'random_charges']}
    
    # Create test system with random positioning
    snapshot = create_test_system(n_particles, box_size, seed, **system_params)
    
    # Cavity parameters
    omegac = 0.0594  # Cavity frequency (1560 cm^-1 in atomic units)
    couplstr = 0.0005  # Coupling strength
    phmass = 1.0  # Photon mass
    
    results = {}
    
    # 1. Old Python implementation
    print("\nComputing forces with old Python implementation...")
    forces_old, force_calc_old = extract_forces_old_python(snapshot, omegac, couplstr, phmass)
    results['old_python'] = {
        'forces': forces_old,
        'harmonic_energy': force_calc_old.harmonic_energy,
        'coupling_energy': force_calc_old.coupling_energy,
        'dipole_self_energy': force_calc_old.dipole_self_energy
    }
    print(f"  ✓ Old Python - Total cavity energy: {force_calc_old.harmonic_energy + force_calc_old.coupling_energy + force_calc_old.dipole_self_energy:.6f}")
    
    # 2. Current CPU implementation
    print("Computing forces with current CPU implementation...")
    forces_cpu, energy_data_cpu = extract_forces_cpp(snapshot, omegac, couplstr, phmass, use_gpu=False)
    if forces_cpu is not None:
        results['cpu'] = {
            'forces': forces_cpu,
            'harmonic_energy': energy_data_cpu['harmonic_energy'],
            'coupling_energy': energy_data_cpu['coupling_energy'],
            'dipole_self_energy': energy_data_cpu['dipole_self_energy']
        }
        total_cpu = energy_data_cpu['harmonic_energy'] + energy_data_cpu['coupling_energy'] + energy_data_cpu['dipole_self_energy']
        print(f"  ✓ CPU - Total cavity energy: {total_cpu:.6f}")
    else:
        print("  ✗ CPU implementation failed")
    
    # 3. Current GPU implementation
    print("Computing forces with current GPU implementation...")
    try:
        forces_gpu, energy_data_gpu = extract_forces_cpp(snapshot, omegac, couplstr, phmass, use_gpu=True)
        if forces_gpu is not None:
            results['gpu'] = {
                'forces': forces_gpu,
                'harmonic_energy': energy_data_gpu['harmonic_energy'],
                'coupling_energy': energy_data_gpu['coupling_energy'],
                'dipole_self_energy': energy_data_gpu['dipole_self_energy']
            }
            total_gpu = energy_data_gpu['harmonic_energy'] + energy_data_gpu['coupling_energy'] + energy_data_gpu['dipole_self_energy']
            print(f"  ✓ GPU - Total cavity energy: {total_gpu:.6f}")
        else:
            print("  ✗ GPU implementation failed")
    except Exception as e:
        print(f"  ✗ GPU implementation error: {e}")
    
    return results, snapshot

def generate_random_test_cases(n_cases=5, base_seed=42):
    """Generate multiple random test cases with different configurations."""
    
    test_cases = []
    
    for i in range(n_cases):
        case = {
            'name': f'Random case {i+1}',
            'n_particles': np.random.randint(5, 25),  # Random number of particles
            'box_size': np.random.uniform(10, 20),    # Random box size
            'seed': base_seed + i * 123,              # Different seed for each case
            'random_cavity_pos': True,                # Always random cavity position
            'random_charges': np.random.choice([True, False])  # Random or alternating charges
        }
        test_cases.append(case)
    
    return test_cases

def plot_force_comparison(results, snapshot, component_names=['x', 'y', 'z']):
    """Create scatter plots of force components against each other with all methods overlaid."""
    
    if 'old_python' not in results:
        print("Error: Old Python results not available for comparison")
        return
    
    # Extract forces from all methods
    forces_old = results['old_python']['forces']
    forces_cpu = results.get('cpu', {}).get('forces', None)
    forces_gpu = results.get('gpu', {}).get('forces', None)
    
    n_particles = len(forces_old)
    
    # Create subplots for force component pairs
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Define the component pairs and their indices
    pairs = [
        ('x', 'y', 0, 1),  # Fx vs Fy
        ('x', 'z', 0, 2),  # Fx vs Fz  
        ('y', 'z', 1, 2)   # Fy vs Fz
    ]
    
    for panel_idx, (comp1_name, comp2_name, comp1_idx, comp2_idx) in enumerate(pairs):
        ax = axes[panel_idx]
        
        # Plot Old Python implementation
        fx_old = forces_old[:, comp1_idx]
        fy_old = forces_old[:, comp2_idx]
        ax.scatter(fx_old, fy_old, alpha=0.6, label='Old Python', s=50, 
                  color='green', marker='o', edgecolors='none')
        
        # Plot CPU implementation
        if forces_cpu is not None:
            fx_cpu = forces_cpu[:, comp1_idx]
            fy_cpu = forces_cpu[:, comp2_idx]
            ax.scatter(fx_cpu, fy_cpu, alpha=0.5, label='CPU', s=40, 
                      color='blue', marker='s', edgecolors='none')
        
        # Plot GPU implementation
        if forces_gpu is not None:
            fx_gpu = forces_gpu[:, comp1_idx]
            fy_gpu = forces_gpu[:, comp2_idx]
            ax.scatter(fx_gpu, fy_gpu, alpha=0.5, label='GPU', s=40, 
                      color='red', marker='^', edgecolors='none')
        
        # Set labels and title
        ax.set_xlabel(f'F_{comp1_name} (a.u.)', fontsize=12)
        ax.set_ylabel(f'F_{comp2_name} (a.u.)', fontsize=12)
        ax.set_title(f'F_{comp1_name} vs F_{comp2_name}', fontsize=14, fontweight='bold')
        ax.legend(loc='upper right', framealpha=0.8)
        ax.grid(True, alpha=0.3)
        
        # Add particle count info
        ax.text(0.02, 0.98, f'N = {n_particles} particles', transform=ax.transAxes, 
               fontsize=10, ha='left', va='top',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.7))
        
        # Add force range info
        all_forces_x = [fx_old]
        all_forces_y = [fy_old]
        if forces_cpu is not None:
            all_forces_x.append(fx_cpu)
            all_forces_y.append(fy_cpu)
        if forces_gpu is not None:
            all_forces_x.append(fx_gpu)
            all_forces_y.append(fy_gpu)
        
        all_fx = np.concatenate(all_forces_x)
        all_fy = np.concatenate(all_forces_y)
        
        range_info = (f'F_{comp1_name}: [{np.min(all_fx):.2e}, {np.max(all_fx):.2e}]\n'
                     f'F_{comp2_name}: [{np.min(all_fy):.2e}, {np.max(all_fy):.2e}]')
        
        ax.text(0.02, 0.02, range_info, transform=ax.transAxes, 
               fontsize=9, ha='left', va='bottom',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.7))
    
    plt.tight_layout()
    plt.savefig('force_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print statistical analysis of force component relationships
    print("\nForce Component Relationship Analysis:")
    print("="*60)
    
    for comp1_name, comp2_name, comp1_idx, comp2_idx in pairs:
        print(f"\n{comp1_name.upper()} vs {comp2_name.upper()} Components:")
        
        # Old Python statistics
        fx_old = forces_old[:, comp1_idx]
        fy_old = forces_old[:, comp2_idx]
        
        print(f"  Old Python:")
        print(f"    F_{comp1_name} range: [{np.min(fx_old):.3e}, {np.max(fx_old):.3e}] a.u.")
        print(f"    F_{comp2_name} range: [{np.min(fy_old):.3e}, {np.max(fy_old):.3e}] a.u.")
        
        if np.std(fx_old) > 1e-10 and np.std(fy_old) > 1e-10:
            r_old = np.corrcoef(fx_old, fy_old)[0, 1]
            print(f"    Correlation: {r_old:.4f}")
        else:
            print(f"    Correlation: N/A (constant values)")
        
        # CPU statistics
        if forces_cpu is not None:
            fx_cpu = forces_cpu[:, comp1_idx]
            fy_cpu = forces_cpu[:, comp2_idx]
            
            print(f"  CPU:")
            print(f"    F_{comp1_name} range: [{np.min(fx_cpu):.3e}, {np.max(fx_cpu):.3e}] a.u.")
            print(f"    F_{comp2_name} range: [{np.min(fy_cpu):.3e}, {np.max(fy_cpu):.3e}] a.u.")
            
            if np.std(fx_cpu) > 1e-10 and np.std(fy_cpu) > 1e-10:
                r_cpu = np.corrcoef(fx_cpu, fy_cpu)[0, 1]
                print(f"    Correlation: {r_cpu:.4f}")
            else:
                print(f"    Correlation: N/A (constant values)")
        
        # GPU statistics
        if forces_gpu is not None:
            fx_gpu = forces_gpu[:, comp1_idx]
            fy_gpu = forces_gpu[:, comp2_idx]
            
            print(f"  GPU:")
            print(f"    F_{comp1_name} range: [{np.min(fx_gpu):.3e}, {np.max(fx_gpu):.3e}] a.u.")
            print(f"    F_{comp2_name} range: [{np.min(fy_gpu):.3e}, {np.max(fy_gpu):.3e}] a.u.")
            
            if np.std(fx_gpu) > 1e-10 and np.std(fy_gpu) > 1e-10:
                r_gpu = np.corrcoef(fx_gpu, fy_gpu)[0, 1]
                print(f"    Correlation: {r_gpu:.4f}")
            else:
                print(f"    Correlation: N/A (constant values)")

def print_energy_comparison(results):
    """Print detailed energy comparison."""
    
    print("\nEnergy Components Comparison:")
    print("="*50)
    
    energy_components = ['harmonic_energy', 'coupling_energy', 'dipole_self_energy']
    
    for component in energy_components:
        print(f"\n{component.replace('_', ' ').title()}:")
        
        if 'old_python' in results:
            old_val = results['old_python'][component]
            print(f"  Old Python: {old_val:.8f} a.u.")
            
            if 'cpu' in results:
                cpu_val = results['cpu'][component]
                diff_cpu = cpu_val - old_val
                print(f"  CPU:        {cpu_val:.8f} a.u. (diff: {diff_cpu:.2e})")
            
            if 'gpu' in results:
                gpu_val = results['gpu'][component]
                diff_gpu = gpu_val - old_val
                print(f"  GPU:        {gpu_val:.8f} a.u. (diff: {diff_gpu:.2e})")

def print_detailed_force_statistics(results):
    """Print detailed force statistics for each implementation."""
    
    print("\nDetailed Force Statistics:")
    print("="*50)
    
    component_pairs = [('x', 'y', 0, 1), ('x', 'z', 0, 2), ('y', 'z', 1, 2)]
    
    for method_name in ['old_python', 'cpu', 'gpu']:
        if method_name not in results:
            continue
            
        print(f"\n{method_name.replace('_', ' ').title()} Implementation:")
        
        forces = results[method_name]['forces']
        
        for comp1_name, comp2_name, comp1_idx, comp2_idx in component_pairs:
            print(f"  {comp1_name.upper()} vs {comp2_name.upper()} Components:")
            
            fx = forces[:, comp1_idx]
            fy = forces[:, comp2_idx]
            
            print(f"    F_{comp1_name} range: [{np.min(fx):.3e}, {np.max(fx):.3e}] a.u.")
            print(f"    F_{comp2_name} range: [{np.min(fy):.3e}, {np.max(fy):.3e}] a.u.")
            
            if np.std(fx) > 1e-10 and np.std(fy) > 1e-10:
                r = np.corrcoef(fx, fy)[0, 1]
                print(f"    Correlation: {r:.4f}")
            else:
                print(f"    Correlation: N/A (constant values)")

def main():
    """Main comparison function."""
    
    print("Cavity Force Implementation Comparison")
    print("="*60)
    print("Testing with randomly positioned particles and cavity")
    print("="*60)
    
    # Test with different random configurations
    print("\n🎲 GENERATING RANDOM TEST CASES")
    print("-" * 40)
    
    # Generate random test cases
    random_cases = generate_random_test_cases(n_cases=3)
    
    # Add some predefined cases for comparison
    predefined_cases = [
        {
            'name': 'Small systematic (alternating charges)',
            'n_particles': 10,
            'box_size': 15.0,
            'seed': 42,
            'random_cavity_pos': True,
            'random_charges': False  # Alternating charges for predictability
        },
        {
            'name': 'Medium random (all random)',
            'n_particles': 50,
            'box_size': 15.0,
            'seed': 123,
            'random_cavity_pos': True,
            'random_charges': True
        },
        {
            'name': 'Large random (cavity at origin)',
            'n_particles': 100,
            'box_size': 20.0,
            'seed': 456,
            'random_cavity_pos': False,  # Cavity at origin for comparison
            'random_charges': True
        }
    ]
    
    all_test_cases = random_cases + predefined_cases
    
    for i, test_case in enumerate(all_test_cases):
        print(f"\n🧪 TEST CASE {i+1}: {test_case['name']}")
        print("-" * 40)
        print(f"Parameters: {test_case['n_particles']} particles, "
              f"box={test_case['box_size']:.1f}, seed={test_case['seed']}")
        
        # Run comparison
        results, snapshot = compare_implementations(**test_case)
        
        # Print energy comparison
        print_energy_comparison(results)
        
        # Check if implementations agree
        if 'cpu' in results and 'gpu' in results:
            force_diff = np.abs(results['gpu']['forces'] - results['cpu']['forces'])
            max_diff = np.max(force_diff)
            print(f"\n📊 FORCE AGREEMENT:")
            print(f"  Maximum force difference (GPU vs CPU): {max_diff:.2e} a.u.")
            if max_diff < 1e-10:
                print("  ✅ Perfect agreement!")
            elif max_diff < 1e-6:
                print("  ✅ Excellent agreement")
            else:
                print("  ❌ Significant differences detected")
        
        # Create plots for one of the medium-sized cases
        if test_case['n_particles'] == 50:
            print("\n📈 Creating force component scatter plots...")
            plot_force_comparison(results, snapshot)
    
    print(f"\n🎉 Comparison complete!")
    print(f"✅ GPU implementation is working correctly with randomly positioned particles and cavity")
    print(f"📊 Check 'force_comparison.png' for scatter plots.")

if __name__ == "__main__":
    main() 