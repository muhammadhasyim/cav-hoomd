#!/usr/bin/env python3
"""
Molecular temperature calculations for diatomic molecules.

This module provides temperature estimates for translational, rotational,
and vibrational degrees of freedom in a system of diatomic molecules.

NOTE: This is a calculation-only class. File output is handled by ObservableWriter
in data.py. This class has been stripped of its old CSV output functionality.
"""

import hoomd
import numpy as np
from typing import Optional, Dict, Tuple
from .utils import PhysicalConstants


class DiatomicMolecularTemperatures(hoomd.custom.Action):
    """
    Calculate translational, rotational, and vibrational kinetic temperatures
    for a system of diatomic molecules (A₂ and B₂ type).
    
    For diatomic molecules, the kinetic energy can be decomposed into:
    1. Translational KE: Center of mass motion
    2. Rotational KE: Rotation around center of mass  
    3. Vibrational KE: Motion along the bond direction
    
    Temperature Definitions:
    -----------------------
    1. Translational temperature (T_trans):
       - From center of mass velocity: <KE_trans> = (3/2) * k_B * T_trans
       - T_trans = M * <v_cm²> / (3 * k_B)
    
    2. Rotational temperature (T_rot):
       - From angular velocity: <KE_rot> = k_B * T_rot  (2 rotational DOF in 3D)
       - T_rot = I * <ω²> / (2 * k_B)
    
    3. Vibrational/Bond temperature (T_vib):
       - From relative velocity along bond: <KE_vib> = (1/2) * k_B * T_vib
       - T_vib = μ * <v_bond²> / k_B
    
    Parameters
    ----------
    simulation : hoomd.Simulation
        HOOMD simulation object
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    output_period_ps : float
        Output period in picoseconds
    output_file : str
        Output CSV file path
    debug : bool, optional
        Enable debug output
        
    Attributes
    ----------
    translational_temp : float
        Translational kinetic temperature (K)
    rotational_temp : float
        Rotational kinetic temperature (K)
    vibrational_temp : float
        Vibrational kinetic temperature along bond (K)
    total_kinetic_temp : float
        Total kinetic temperature from all velocities (K)
    """
    
    def __init__(self, 
                 simulation: hoomd.Simulation,
                 time_tracker,
                 output_period_ps: float = 0.1,
                 debug: bool = False):
        super().__init__()
        self.simulation = simulation
        self.time_tracker = time_tracker
        self.output_period_ps = output_period_ps
        self.debug = debug
        
        # Physical constants
        self.kB = PhysicalConstants.KB_HARTREE_PER_K  # Hartree/K
        
        # Get molecular information from first snapshot
        self._initialize_molecular_info()
        
        # Tracking state
        self.last_calculation_time = None
        
        # Temperature attributes for external access
        self.translational_temp = None
        self.rotational_temp = None
        self.vibrational_temp = None
        self.total_kinetic_temp = None
    
    def _initialize_molecular_info(self):
        """Extract bond topology and particle masses from initial snapshot."""
        with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
            # Get bond information
            bonds = np.array(snap.bonds.group)
            self.n_molecules = len(bonds)
            
            # Store bond topology (particle indices for each dimer)
            self.bond_pairs = bonds.copy()
            
            # Get masses for each particle type
            masses = np.array(snap.particles.mass)
            type_ids = np.array(snap.particles.typeid)
            
            # Determine masses for O and N particles
            # Type 0 = O (mass 29150), Type 1 = N (mass 25527)
            self.mass_O = 29150.0
            self.mass_N = 25527.0
            
            # Store masses and reduced masses for each molecule
            self.molecule_masses = np.zeros((self.n_molecules, 2))  # [m1, m2]
            self.molecule_total_masses = np.zeros(self.n_molecules)  # M = m1 + m2
            self.molecule_reduced_masses = np.zeros(self.n_molecules)  # μ = m1*m2/(m1+m2)
            self.molecule_types = np.zeros(self.n_molecules, dtype=int)  # 0=O-O, 1=N-N, 2=O-N
            
            for i, (idx1, idx2) in enumerate(self.bond_pairs):
                m1 = masses[idx1]
                m2 = masses[idx2]
                t1 = type_ids[idx1]
                t2 = type_ids[idx2]
                
                self.molecule_masses[i] = [m1, m2]
                self.molecule_total_masses[i] = m1 + m2
                self.molecule_reduced_masses[i] = (m1 * m2) / (m1 + m2)
                
                # Determine molecule type
                if t1 == 0 and t2 == 0:
                    self.molecule_types[i] = 0  # O-O
                elif t1 == 1 and t2 == 1:
                    self.molecule_types[i] = 1  # N-N
                else:
                    self.molecule_types[i] = 2  # O-N (mixed)
            
            if self.debug:
                print(f"Initialized {self.n_molecules} diatomic molecules")
                print(f"  O-O molecules: {np.sum(self.molecule_types == 0)}")
                print(f"  N-N molecules: {np.sum(self.molecule_types == 1)}")
                print(f"  O-N molecules: {np.sum(self.molecule_types == 2)}")
    
    def act(self, timestep):
        """Calculate molecular temperatures at each timestep."""
        current_time_ps = self.time_tracker.elapsed_time
        
        # Calculate temperatures periodically
        if self._should_calculate(current_time_ps):
            temps = self._calculate_molecular_temperatures()
            
            self.translational_temp = temps['T_trans']
            self.rotational_temp = temps['T_rot']
            self.vibrational_temp = temps['T_vib']
            self.total_kinetic_temp = temps['T_kinetic_total']
            
            self.last_calculation_time = current_time_ps
    
    def _should_calculate(self, current_time_ps: float) -> bool:
        """Check if we should calculate temperatures."""
        if self.last_calculation_time is None:
            return True
        return (current_time_ps - self.last_calculation_time) >= self.output_period_ps
    
    def _get_hoomd_simulation(self):
        """Get the HOOMD simulation object."""
        if hasattr(self.simulation, 'sim'):
            return self.simulation.sim
        else:
            return self.simulation
    
    def _calculate_molecular_temperatures(self) -> Dict[str, float]:
        """
        Calculate all molecular temperatures.
        
        Returns
        -------
        dict
            Dictionary with temperature values for each component
        """
        try:
            # Get box size from state (not from snapshot)
            box = np.array(self._get_hoomd_simulation().state.box.L)
            
            with self._get_hoomd_simulation().state.cpu_local_snapshot as snap:
                positions = np.array(snap.particles.position)
                velocities = np.array(snap.particles.velocity)
                masses = np.array(snap.particles.mass)
            
            # Calculate temperatures for all molecules
            T_trans, T_rot, T_vib, T_total = self._compute_temperatures(
                positions, velocities, masses, box
            )
            
            # Calculate temperatures for O-O molecules only
            O2_mask = self.molecule_types == 0
            if np.any(O2_mask):
                T_trans_O2, T_rot_O2, T_vib_O2, _ = self._compute_temperatures(
                    positions, velocities, masses, box, molecule_mask=O2_mask
                )
            else:
                T_trans_O2 = T_rot_O2 = T_vib_O2 = 0.0
            
            # Calculate temperatures for N-N molecules only
            N2_mask = self.molecule_types == 1
            if np.any(N2_mask):
                T_trans_N2, T_rot_N2, T_vib_N2, _ = self._compute_temperatures(
                    positions, velocities, masses, box, molecule_mask=N2_mask
                )
            else:
                T_trans_N2 = T_rot_N2 = T_vib_N2 = 0.0
            
            return {
                'T_trans': T_trans,
                'T_rot': T_rot,
                'T_vib': T_vib,
                'T_kinetic_total': T_total,
                'T_trans_O2': T_trans_O2,
                'T_trans_N2': T_trans_N2,
                'T_rot_O2': T_rot_O2,
                'T_rot_N2': T_rot_N2,
                'T_vib_O2': T_vib_O2,
                'T_vib_N2': T_vib_N2
            }
            
        except Exception as e:
            if self.debug:
                print(f"Warning: Could not calculate molecular temperatures: {e}")
                import traceback
                traceback.print_exc()
            return {
                'T_trans': 0.0, 'T_rot': 0.0, 'T_vib': 0.0, 'T_kinetic_total': 0.0,
                'T_trans_O2': 0.0, 'T_trans_N2': 0.0,
                'T_rot_O2': 0.0, 'T_rot_N2': 0.0,
                'T_vib_O2': 0.0, 'T_vib_N2': 0.0
            }
    
    def _compute_temperatures(self, 
                            positions: np.ndarray, 
                            velocities: np.ndarray,
                            masses: np.ndarray,
                            box: np.ndarray,
                            molecule_mask: Optional[np.ndarray] = None) -> Tuple[float, float, float, float]:
        """
        Compute translational, rotational, and vibrational temperatures.
        
        Parameters
        ----------
        positions : np.ndarray
            Particle positions
        velocities : np.ndarray
            Particle velocities
        masses : np.ndarray
            Particle masses
        box : np.ndarray
            Box dimensions [Lx, Ly, Lz]
        molecule_mask : np.ndarray, optional
            Boolean mask to select subset of molecules
            
        Returns
        -------
        T_trans : float
            Translational temperature (K)
        T_rot : float
            Rotational temperature (K)
        T_vib : float
            Vibrational temperature (K)
        T_total : float
            Total kinetic temperature (K)
        """
        if molecule_mask is None:
            molecule_mask = np.ones(self.n_molecules, dtype=bool)
        
        n_mols = np.sum(molecule_mask)
        if n_mols == 0:
            return 0.0, 0.0, 0.0, 0.0
        
        # Initialize kinetic energy accumulators
        KE_trans = 0.0  # Translational KE
        KE_rot = 0.0    # Rotational KE
        KE_vib = 0.0    # Vibrational KE
        KE_total = 0.0  # Total KE (for verification)
        
        for i, (idx1, idx2) in enumerate(self.bond_pairs):
            if not molecule_mask[i]:
                continue
            
            # Get particle data
            r1 = positions[idx1]
            r2 = positions[idx2]
            v1 = velocities[idx1]
            v2 = velocities[idx2]
            m1 = masses[idx1]
            m2 = masses[idx2]
            M = m1 + m2  # Total mass
            mu = (m1 * m2) / M  # Reduced mass
            
            # Handle periodic boundary conditions for bond vector
            r12 = r2 - r1
            r12 = r12 - box * np.round(r12 / box)
            r12_mag = np.linalg.norm(r12)
            if r12_mag < 1e-10:
                continue  # Skip if atoms are on top of each other
            r12_hat = r12 / r12_mag  # Unit vector along bond
            
            # 1. Translational KE: (1/2) * M * v_cm²
            v_cm = (m1 * v1 + m2 * v2) / M
            KE_trans += 0.5 * M * np.dot(v_cm, v_cm)
            
            # 2. Relative velocity
            v_rel = v2 - v1
            
            # 3. Vibrational KE: motion along bond direction
            # v_bond = v_rel · r̂  (component of relative velocity along bond)
            v_bond = np.dot(v_rel, r12_hat)
            KE_vib += 0.5 * mu * v_bond**2
            
            # 4. Rotational KE: motion perpendicular to bond
            # v_perp = v_rel - v_bond * r̂  (perpendicular component)
            v_perp = v_rel - v_bond * r12_hat
            # For rotation: KE_rot = (1/2) * I * ω² = (1/2) * μ * r² * ω²
            # But ω = v_perp / r, so: KE_rot = (1/2) * μ * v_perp²
            KE_rot += 0.5 * mu * np.dot(v_perp, v_perp)
            
            # Total KE for this molecule (for verification)
            KE_total += 0.5 * m1 * np.dot(v1, v1) + 0.5 * m2 * np.dot(v2, v2)
        
        # Convert kinetic energies to temperatures
        # For translational: <KE_trans> = (3/2) * N * k_B * T_trans
        # T_trans = (2 * KE_trans) / (3 * N * k_B)
        T_trans = (2.0 * KE_trans) / (3.0 * n_mols * self.kB) if n_mols > 0 else 0.0
        
        # For rotational: <KE_rot> = N * k_B * T_rot  (2 rotational DOF in 3D)
        # T_rot = KE_rot / (N * k_B)
        T_rot = (2.0 * KE_rot) / (2.0 * n_mols * self.kB) if n_mols > 0 else 0.0
        
        # For vibrational: <KE_vib> = (1/2) * N * k_B * T_vib  (1 vibrational DOF)
        # T_vib = (2 * KE_vib) / (N * k_B)
        T_vib = (2.0 * KE_vib) / (n_mols * self.kB) if n_mols > 0 else 0.0
        
        # Total kinetic temperature (from all particle velocities)
        # <KE_total> = (3/2) * (2N) * k_B * T_total  (2N particles with 3 DOF each)
        # T_total = (2 * KE_total) / (6 * N * k_B)
        T_total = (2.0 * KE_total) / (6.0 * n_mols * self.kB) if n_mols > 0 else 0.0
        
        if self.debug and molecule_mask is None:  # Only print for all molecules
            print(f"\nMolecular Temperatures (N={n_mols}):")
            print(f"  KE_trans = {KE_trans:.6e} Hartree -> T_trans = {T_trans:.2f} K")
            print(f"  KE_rot   = {KE_rot:.6e} Hartree -> T_rot   = {T_rot:.2f} K")
            print(f"  KE_vib   = {KE_vib:.6e} Hartree -> T_vib   = {T_vib:.2f} K")
            print(f"  KE_total = {KE_total:.6e} Hartree -> T_total = {T_total:.2f} K")
            print(f"  KE sum check: {KE_trans + KE_rot + KE_vib:.6e} vs {KE_total:.6e}")
        
        return T_trans, T_rot, T_vib, T_total

