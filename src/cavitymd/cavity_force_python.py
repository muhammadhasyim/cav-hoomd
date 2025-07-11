# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Pure Python implementation of cavity force computation."""

import hoomd
import numpy as np
from .utils import unwrap_positions


class CavityForcePython(hoomd.md.force.Custom):
    """
    Pure Python implementation of the cavity force.
    
    This serves as a fallback when the C++/CUDA implementations are not available.
    Implements the same physics as the optimized versions but with Python performance.
    
    Now includes dissipation support for cavity mode damping.
    """
    
    def __init__(self, kvector, couplstr, omegac, phmass=1.0, dissipation=0.0):
        super().__init__(aniso=False)
        
        # Store parameters
        self.kvector = np.array(kvector)
        self.couplstr = couplstr
        self.omegac = omegac
        self.phmass = phmass
        self.dissipation = dissipation
        self.K = phmass * omegac**2
        
        # Initialize energy components
        self.harmonic_energy = 0.0
        self.coupling_energy = 0.0
        self.dipole_self_energy = 0.0
        
        print(f"CavityForcePython initialized:")
        print(f"  Coupling strength: {self.couplstr:.6f} a.u.")
        print(f"  Cavity frequency: {self.omegac:.6f} a.u.")
        print(f"  Photon mass: {self.phmass:.6f} a.u.")
        print(f"  Spring constant K: {self.K:.6f} a.u.")
        print(f"  Dissipation rate: {self.dissipation:.6f} a.u.")
    
    @property
    def total_cavity_energy(self):
        """Total cavity energy: sum of all components"""
        return self.harmonic_energy + self.coupling_energy + self.dipole_self_energy
    
    def set_forces(self, timestep):
        """
        Compute cavity forces and potential energy.
        
        This is called by HOOMD at each timestep to compute forces.
        Implements the cavity-molecule interaction Hamiltonian with dissipation.
        """
        with self._state.cpu_local_snapshot as snap:
            try:
                # Find cavity particle (type ID = 1)
                cavity_mask = snap.particles.typeid == 1
                cavity_indices = np.where(cavity_mask)[0]
                
                if len(cavity_indices) == 0:
                    print("Warning: No cavity particle found (typeid=1)")
                    self._zero_forces_and_energy()
                    return
                elif len(cavity_indices) > 1:
                    print(f"Warning: Multiple cavity particles found ({len(cavity_indices)}), using first one")
                
                cavity_idx = cavity_indices[0]
                
                # Get unwrapped positions
                box_lengths = np.array([
                    snap.global_box.L[0],
                    snap.global_box.L[1], 
                    snap.global_box.L[2]
                ])
                
                unwrapped_positions = unwrap_positions(
                    snap.particles.position,
                    snap.particles.image,
                    box_lengths
                )
                
                # Calculate molecular dipole moment (only x,y components)
                dipole_moment = np.dot(snap.particles.charge, unwrapped_positions)
                dipole_xy = dipole_moment.copy()
                dipole_xy[2] = 0.0  # Zero out z-component
                
                # Get cavity particle position and velocity
                cavity_position = unwrapped_positions[cavity_idx]
                cavity_xy = cavity_position.copy()
                cavity_xy[2] = 0.0  # Zero out z-component
                
                # Get cavity particle velocity for dissipation
                cavity_velocity = snap.particles.velocity[cavity_idx]
                
                # Compute energy components
                # 1. Harmonic energy: (1/2) * K * |q|²
                self.harmonic_energy = 0.5 * self.K * np.dot(cavity_position, cavity_position)
                
                # 2. Coupling energy: g * (q_xy · d_xy)
                self.coupling_energy = self.couplstr * np.dot(cavity_xy, dipole_xy)
                
                # 3. Dipole self-energy: (g²/2K) * |d_xy|²
                self.dipole_self_energy = 0.5 * (self.couplstr**2 / self.K) * np.dot(dipole_xy, dipole_xy)
                
                # Total cavity energy
                total_energy = self.total_cavity_energy
                
                # Compute forces
                with self.cpu_local_force_arrays as arrays:
                    # Initialize arrays
                    arrays.force[:] = 0.0
                    arrays.potential_energy[:] = 0.0
                    
                    # DO NOT assign energy to particle potential energy - prevents double-counting
                    # Energy is accessed directly through force object attributes
                    arrays.potential_energy[cavity_idx] = 0.0
                    
                    # Force on molecular particles: F_i = -g * q_i * [q_xy + (g/K) * d_xy]
                    # Only x,y components contribute
                    force_factor = cavity_xy + (self.couplstr / self.K) * dipole_xy
                    num_particles = len(snap.particles.position)
                    for i in range(num_particles):
                        if i != cavity_idx:  # Skip cavity particle
                            charge = snap.particles.charge[i]
                            force_molecular = -self.couplstr * charge * force_factor
                            arrays.force[i] = force_molecular
                    
                    # Force on cavity particle: F_cavity = -K * q - g * d_xy - γ * v_cavity
                    # Add dissipation term: -γ * v_cavity (friction force)
                    force_cavity = (-self.K * cavity_position - 
                                   self.couplstr * dipole_xy - 
                                   self.dissipation * cavity_velocity)
                    arrays.force[cavity_idx] = force_cavity
                    
            except Exception as e:
                print(f"Error in cavity force calculation: {e}")
                self._zero_forces_and_energy()
    
    def _zero_forces_and_energy(self):
        """Zero out all forces and energies."""
        with self.cpu_local_force_arrays as arrays:
            arrays.force[:] = 0.0
            arrays.potential_energy[:] = 0.0
        
        self.harmonic_energy = 0.0
        self.coupling_energy = 0.0
        self.dipole_self_energy = 0.0 