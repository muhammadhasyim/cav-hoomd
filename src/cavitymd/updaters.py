# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Custom updaters for cavity MD simulations."""

import hoomd
from hoomd.operation import Updater
from hoomd.variant import Variant
from hoomd.trigger import Trigger, On
import numpy as np
from typing import Optional

# Try to import C++ compiled module
try:
    from . import _cavitymd
    _cpp_available = True
except ImportError:
    _cpp_available = False

class CavityParticleDisplacer(Updater):
    """Displaces the cavity particle to its new equilibrium position.

    Args:
        trigger (hoomd.trigger.Trigger): The trigger to activate this updater.
        couplstr (hoomd.variant.Variant): The coupling strength variant.
        omegac (float): The cavity frequency.
        phmass (float): The photon mass.
    """

    def __init__(self, trigger, couplstr, omegac, phmass=1.0):
        super().__init__(trigger)
        
        # Store parameters for later use when attached to simulation
        self._couplstr = couplstr
        self._omegac = omegac
        self._phmass = phmass
        self._cpp_obj = None

        if not _cpp_available:
            raise RuntimeError("C++ implementation of CavityParticleDisplacer not available.")

    def _attach_hook(self):
        """Called when the updater is attached to a simulation."""
        # Now we can create the C++ object since we have access to the simulation
        if _cpp_available and self._cpp_obj is None:
            self._cpp_obj = _cavitymd.CavityParticleDisplacer(
                self._simulation.state._cpp_sys_def,
                self.trigger,
                self._couplstr,
                self._omegac,
                self._phmass
            )


class HarmonicBondReset(hoomd.custom.Action):
    """
    Perform one-time thermal reset of harmonic bond stretch DOF.
    
    This class samples bond lengths and relative velocities from the canonical
    distribution at a target temperature, preserving COM motion and molecular
    orientation. The reset is performed exactly once when triggered.
    
    Parameters
    ----------
    bond_type : str
        The bond type name to reset (e.g., 'bond')
    spring_constant : float
        Harmonic spring constant K (in HOOMD units)
    equilibrium_length : float
        Equilibrium bond length r0 (in HOOMD units)
    temperature : float
        Target temperature for internal stretch DOF (in HOOMD units)
    kB : float, optional
        Boltzmann constant (default: 1.0 for reduced units)
    seed : int, optional
        Random seed for reproducibility (default: 42)
    
    Attributes
    ----------
    enabled : bool
        Set to True to trigger the reset on next call
    has_reset : bool
        True if reset has been performed
    
    Notes
    -----
    The reset preserves:
    - Center-of-mass position
    - Center-of-mass velocity
    - Bond orientation (axis direction)
    
    Only the internal stretch DOF (bond length and relative velocity along
    bond axis) are resampled from the thermal distribution.
    
    Mathematical Details
    --------------------
    For bond extension x = r - r0 from equilibrium:
        x ~ N(0, σ_x²)  where σ_x² = k_B T / K
    
    For relative velocity u along bond axis:
        u ~ N(0, σ_u²)  where σ_u² = k_B T / μ
    
    where μ = m1*m2/(m1+m2) is the reduced mass.
    
    COM-preserving updates:
        r1' = r1 - (m2/M) * Δr * ê
        r2' = r2 + (m1/M) * Δr * ê
        v1' = v1 - (m2/M) * δu * ê
        v2' = v2 + (m1/M) * δu * ê
    
    where ê is the bond axis unit vector, M = m1 + m2.
    
    Examples
    --------
    >>> # Create reset action
    >>> bond_reset = HarmonicBondReset(
    ...     bond_type='bond',
    ...     spring_constant=10000.0,
    ...     equilibrium_length=1.0,
    ...     temperature=300.0,
    ...     kB=1.0
    ... )
    >>> 
    >>> # Add to simulation as a custom updater
    >>> reset_updater = hoomd.update.CustomUpdater(
    ...     action=bond_reset,
    ...     trigger=hoomd.trigger.Periodic(1000)  # Check every 1000 steps
    ... )
    >>> sim.operations.updaters.append(reset_updater)
    >>> 
    >>> # Trigger the reset at any time
    >>> bond_reset.enabled = True
    >>> 
    >>> # Reset happens on next update, then bond_reset.has_reset becomes True
    """
    
    def __init__(
        self,
        bond_type: str = None,
        spring_constant: float = None,
        equilibrium_length: float = None,
        temperature: float = None,
        bond_params: dict = None,
        kB: float = 1.0,
        seed: int = 42
    ):
        """
        Initialize the harmonic bond reset action.
        
        Can handle either a single bond type or multiple bond types.
        
        Parameters
        ----------
        bond_type : str, optional
            Single bond type name (for single bond type mode)
        spring_constant : float, optional
            Spring constant K for single bond type
        equilibrium_length : float, optional
            Equilibrium length r0 for single bond type
        temperature : float
            Target temperature for internal stretch DOF
        bond_params : dict, optional
            Dictionary mapping bond type names to {'K': ..., 'r0': ...}
            Example: {'OO': {'K': 0.73204, 'r0': 2.281655158},
                      'NN': {'K': 1.4325, 'r0': 2.0743522177}}
        kB : float, optional
            Boltzmann constant (default: 1.0)
        seed : int, optional
            Random seed (default: 42)
        """
        # Handle single bond type vs. multiple bond types
        if bond_params is not None:
            # Multiple bond types mode
            self.bond_params = bond_params
            self.multi_bond_mode = True
            self.bond_type = None  # Not used in multi-bond mode
            self.K = None
            self.r0 = None
        else:
            # Single bond type mode (backward compatibility)
            if bond_type is None or spring_constant is None or equilibrium_length is None:
                raise ValueError("Must provide either bond_params dict OR (bond_type, spring_constant, equilibrium_length)")
            self.bond_type = bond_type
            self.K = spring_constant
            self.r0 = equilibrium_length
            self.bond_params = {bond_type: {'K': spring_constant, 'r0': equilibrium_length}}
            self.multi_bond_mode = False
        
        self.T = temperature
        self.kB = kB
        self.seed = seed
        
        # Control flags
        self.enabled = False
        self.has_reset = False
        
        # Initialize random number generator
        self.rng = np.random.RandomState(seed)
        
        # Empirical data for harmonic fictive temperature calculation
        self.empirical_data = None
        
        # Validation
        for btype, params in self.bond_params.items():
            if params['K'] <= 0:
                raise ValueError(f"Spring constant must be positive for {btype}, got {params['K']}")
            if params['r0'] <= 0:
                raise ValueError(f"Equilibrium length must be positive for {btype}, got {params['r0']}")
        if self.T is not None and self.T < 0:
            raise ValueError(f"Temperature must be non-negative, got {self.T}")
    
    def set_empirical_data(self, empirical_data):
        """
        Set empirical temperature data for harmonic fictive temperature calculation.
        
        Parameters
        ----------
        empirical_data : EmpiricalTemperatureData
            Empirical data object for harmonic energy component
        """
        self.empirical_data = empirical_data
    
    def act(self, timestep: int):
        """
        Perform the bond reset if enabled and not yet performed.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep
        """
        # Only reset once when enabled
        if not self.enabled or self.has_reset:
            return
        
        # Auto-detect bond parameters if we have a placeholder
        if hasattr(self, '_sim_obj') and 'placeholder' in self.bond_params:
            self._auto_detect_bond_params()
        
        # Get snapshot (use get_snapshot() method)
        snap = self._state.get_snapshot()
        
        if snap.communicator.rank == 0:
            print(f"\n{'='*60}")
            print(f"HARMONIC BOND RESET at timestep {timestep}")
            print(f"{'='*60}")
            if self.multi_bond_mode:
                print(f"Bond types: {list(self.bond_params.keys())}")
                for btype, params in self.bond_params.items():
                    print(f"  {btype}: K={params['K']:.6f}, r0={params['r0']:.6f}")
            else:
                print(f"Bond type: {self.bond_type}")
                print(f"Spring constant K: {self.K}")
                print(f"Equilibrium length r0: {self.r0}")
            print(f"Temperature T: {self.T}")
            print(f"Boltzmann constant kB: {self.kB}")
            
            # Calculate harmonic energy BEFORE reset
            E_harmonic_before = self._calculate_harmonic_energy(snap)
            T_fictive_before = self._harmonic_energy_to_temperature(E_harmonic_before, snap)
            
            print(f"\nBEFORE RESET:")
            print(f"  Total harmonic energy: {E_harmonic_before:.6f} Hartree")
            print(f"  Harmonic fictive temperature: {T_fictive_before:.2f} K")
            
            # Perform the reset
            n_bonds_reset = self._reset_bonds(snap)
            
            # Calculate harmonic energy AFTER reset
            E_harmonic_after = self._calculate_harmonic_energy(snap)
            T_fictive_after = self._harmonic_energy_to_temperature(E_harmonic_after, snap)
            
            print(f"\nAFTER RESET:")
            print(f"  Total harmonic energy: {E_harmonic_after:.6f} Hartree")
            print(f"  Harmonic fictive temperature: {T_fictive_after:.2f} K")
            print(f"  Target temperature: {self.T:.2f} K")
            print(f"  Energy change: {(E_harmonic_after - E_harmonic_before):.6f} Hartree")
            
            print(f"\n✓ Reset {n_bonds_reset} bonds")
            print(f"{'='*60}\n")
        
        # Set the modified snapshot back to the state
        self._state.set_snapshot(snap)
        
        # Mark as completed
        self.has_reset = True
        self.enabled = False  # Auto-disable after reset
    
    def _auto_detect_bond_params(self):
        """Auto-detect bond parameters from simulation force field."""
        try:
            # Find harmonic bond force
            bond_force = None
            for f in self._sim_obj.sim.operations.integrator.forces:
                if hasattr(f, 'params') and len(f.params) > 0:
                    # Check if this looks like a bond force (has 'k' and 'r0' parameters)
                    first_key = list(f.params.keys())[0]
                    if 'k' in f.params[first_key] and 'r0' in f.params[first_key]:
                        bond_force = f
                        break
            
            if bond_force is None:
                raise ValueError("Could not find harmonic bond force in simulation")
            
            # Extract parameters for all bond types
            bond_params_dict = {}
            for bond_type in bond_force.params.keys():
                params = bond_force.params[bond_type]
                bond_params_dict[bond_type] = {
                    'K': params['k'],
                    'r0': params['r0']
                }
            
            if not bond_params_dict:
                raise ValueError("No bond types found in harmonic bond force")
            
            # Update bond parameters
            self.bond_params = bond_params_dict
            self.multi_bond_mode = True
            
            print(f"\n✓ Auto-detected {len(bond_params_dict)} bond type(s) from force field")
            for btype, params in bond_params_dict.items():
                print(f"  {btype}: K={params['K']:.6f}, r0={params['r0']:.6f}")
            print()
                    
        except Exception as e:
            raise ValueError(f"Could not auto-detect bond parameters from simulation. Error: {e}")
    
    def _reset_bonds(self, snap) -> int:
        """
        Reset all bonds of the specified type.
        
        Parameters
        ----------
        snap : hoomd.Snapshot
            System snapshot
        
        Returns
        -------
        int
            Number of bonds reset
        """
        n_bonds_reset = 0
        bonds_reset_by_type = {btype: 0 for btype in self.bond_params.keys()}
        
        # Get bond data
        bonds = snap.bonds
        positions = snap.particles.position
        velocities = snap.particles.velocity
        masses = snap.particles.mass
        box = snap.configuration.box  # CRITICAL: Get box for PBC handling
        
        # Iterate over all bonds
        for bond_idx in range(bonds.N):
            # Get bond type name
            bond_type_name = bonds.types[bonds.typeid[bond_idx]]
            
            # Check if this bond type should be reset
            if bond_type_name not in self.bond_params:
                continue
            
            # Get particle indices
            i, j = bonds.group[bond_idx]
            
            # Get particle data
            r1 = positions[i]
            r2 = positions[j]
            v1 = velocities[i]
            v2 = velocities[j]
            m1 = masses[i]
            m2 = masses[j]
            
            # Get bond-specific parameters
            K = self.bond_params[bond_type_name]['K']
            r0 = self.bond_params[bond_type_name]['r0']
            
            # Reset this bond using Option B (thermal sampling)
            # CRITICAL: Pass box for PBC handling
            r1_new, r2_new, v1_new, v2_new = self._thermal_reset_bond(
                r1, r2, v1, v2, m1, m2, K, r0, box
            )
            
            # Update positions and velocities
            positions[i] = r1_new
            positions[j] = r2_new
            velocities[i] = v1_new
            velocities[j] = v2_new
            
            n_bonds_reset += 1
            bonds_reset_by_type[bond_type_name] += 1
        
        # Print breakdown by type
        if n_bonds_reset > 0:
            print(f"\nBonds reset by type:")
            for btype, count in bonds_reset_by_type.items():
                if count > 0:
                    print(f"  {btype}: {count} bonds")
        
        return n_bonds_reset
    
    def _thermal_reset_bond(
        self,
        r1: np.ndarray,
        r2: np.ndarray,
        v1: np.ndarray,
        v2: np.ndarray,
        m1: float,
        m2: float,
        K: float,
        r0: float,
        box: np.ndarray
    ) -> tuple:
        """
        Thermally reset a single bond using Option B (exact thermal sampling).
        
        CRITICAL: Properly handles periodic boundary conditions using minimum image convention.
        
        Parameters
        ----------
        r1, r2 : np.ndarray
            Positions of atoms 1 and 2 (shape: (3,))
        v1, v2 : np.ndarray
            Velocities of atoms 1 and 2 (shape: (3,))
        m1, m2 : float
            Masses of atoms 1 and 2
        K : float
            Spring constant
        r0 : float
            Equilibrium bond length
        box : np.ndarray
            Box dimensions [Lx, Ly, Lz, ...] (from snap.configuration.box)
        
        Returns
        -------
        r1_new, r2_new, v1_new, v2_new : np.ndarray
            Updated positions and velocities
        """
        # Calculate mass properties
        M = m1 + m2
        mu = (m1 * m2) / M
        
        # CRITICAL: Apply minimum image convention for bond vector (PBC!)
        # This ensures we get the shortest distance between atoms
        bond_vec = r2 - r1
        Lx, Ly, Lz = box[:3]
        bond_vec[0] -= Lx * np.round(bond_vec[0] / Lx)
        bond_vec[1] -= Ly * np.round(bond_vec[1] / Ly)
        bond_vec[2] -= Lz * np.round(bond_vec[2] / Lz)
        
        # Current bond length (using minimum image)
        r_current = np.linalg.norm(bond_vec)
        
        # Bond axis unit vector (from minimum image)
        e = bond_vec / r_current
        
        # ========== POSITION UPDATE (sample bond extension) ==========
        # Sample extension from thermal distribution
        # K is in atomic units (Hartree/Bohr^2), kB in Hartree/K, T in K
        sigma_x = np.sqrt(self.kB * self.T / K)
        x = self.rng.normal(0, sigma_x)
        
        # New bond length
        r_new = r0 + x
        
        # CRITICAL: With PBC, we CANNOT reliably preserve COM position from wrapped coordinates!
        # Strategy: Keep particle 1 fixed, place particle 2 at correct bond length
        # This avoids COM calculation from wrapped coords (which fails at boundaries)
        # COM velocity is still preserved by the velocity update below
        r1_new = r1.copy()
        r2_new = r1 + r_new * e  # Place r2 at distance r_new along bond axis from r1
        
        # Manually wrap r2_new back into the periodic box
        # HOOMD expects positions to be within [-L/2, L/2] for each dimension
        r2_new[0] = r2_new[0] - Lx * np.floor(r2_new[0] / Lx + 0.5)
        r2_new[1] = r2_new[1] - Ly * np.floor(r2_new[1] / Ly + 0.5)
        r2_new[2] = r2_new[2] - Lz * np.floor(r2_new[2] / Lz + 0.5)
        
        # ========== VELOCITY UPDATE (sample relative velocity) ==========
        # Current relative velocity along bond axis
        u_cur = np.dot(v2 - v1, e)
        
        # Sample new relative velocity from thermal distribution
        # CRITICAL: Masses in HOOMD are in AMU, but atomic units use electron masses
        # Must convert to ensure correct thermal velocity distribution
        from .utils import PhysicalConstants
        mu_atomic = mu * PhysicalConstants.AMU_TO_ELECTRON_MASS  # Convert to electron masses
        sigma_u = np.sqrt(self.kB * self.T / mu_atomic)
        u_star = self.rng.normal(0, sigma_u)
        
        # Change in relative velocity
        du = u_star - u_cur
        
        # COM-preserving velocity update along bond axis
        v1_new = v1 - (m2 / M) * du * e
        v2_new = v2 + (m1 / M) * du * e
        
        return r1_new, r2_new, v1_new, v2_new
    
    def _calculate_harmonic_energy(self, snap):
        """
        Calculate total harmonic bond potential energy.
        
        Parameters
        ----------
        snap : gsd.hoomd.Snapshot
            System snapshot
            
        Returns
        -------
        float
            Total harmonic energy in Hartree
        """
        bonds = snap.bonds.group
        positions = snap.particles.position
        box = snap.configuration.box
        
        E_total = 0.0
        
        for bond_idx in range(len(bonds)):
            i, j = bonds[bond_idx]
            bond_type_idx = snap.bonds.typeid[bond_idx]
            bond_type = snap.bonds.types[bond_type_idx]
            
            # Get bond parameters
            if self.multi_bond_mode:
                if bond_type not in self.bond_params:
                    continue
                K = self.bond_params[bond_type]['K']
                r0 = self.bond_params[bond_type]['r0']
            else:
                K = self.K
                r0 = self.r0
            
            # Calculate bond vector (with PBC)
            r1 = positions[i]
            r2 = positions[j]
            
            # Apply minimum image convention for periodic boundaries
            dr = r2 - r1
            Lx, Ly, Lz = box[:3]
            dr[0] -= Lx * np.round(dr[0] / Lx)
            dr[1] -= Ly * np.round(dr[1] / Ly)
            dr[2] -= Lz * np.round(dr[2] / Lz)
            
            r = np.linalg.norm(dr)
            
            # Harmonic potential: E = 0.5 * K * (r - r0)^2
            E_bond = 0.5 * K * (r - r0)**2
            E_total += E_bond
        
        return E_total
    
    def _harmonic_energy_to_temperature(self, E_harmonic, snap):
        """
        Convert harmonic bond energy to fictive temperature using empirical data.
        
        CRITICAL: Must exclude cavity particle from N_particles count!
        Empirical data was calibrated with 500 molecular particles.
        
        Parameters
        ----------
        E_harmonic : float
            Total harmonic bond POTENTIAL energy in Hartree
        snap : gsd.hoomd.Snapshot
            System snapshot (to count particles, excluding cavity)
            
        Returns
        -------
        float
            Fictive temperature in Kelvin
        """
        if E_harmonic <= 0:
            return 0.0
        
        # Use empirical data if available (preferred method)
        if self.empirical_data is not None:
            try:
                # CRITICAL: Count only MOLECULAR particles (exclude cavity particle type 'X')
                # The empirical data was calibrated with 500 molecular particles
                N_total = len(snap.particles.position)
                N_molecular = N_total - 1  # Assume 1 cavity particle (type 'X')
                
                temperature_K = self.empirical_data.calculate_systemic_temperature(
                    E_harmonic,
                    num_particles=N_molecular
                )
                return temperature_K
            except Exception as e:
                print(f"Warning: Empirical temperature calculation failed: {e}, using direct method")
        
        # Fallback to direct calculation: T = 2*E/(N_bonds*kB)
        # For harmonic potential energy: E_PE = N_bonds * (1/2) * kB * T
        N_bonds = len(snap.bonds.group)
        if N_bonds == 0:
            return 0.0
        
        temperature_K = 2.0 * E_harmonic / (N_bonds * self.kB)
        return temperature_K
    
    def reset_flag(self):
        """Reset the has_reset flag to allow another reset."""
        self.has_reset = False
    
    def __repr__(self):
        if self.multi_bond_mode:
            bond_info = f"bond_types={list(self.bond_params.keys())}"
        else:
            bond_info = f"bond_type='{self.bond_type}', K={self.K}, r0={self.r0}"
        return (f"HarmonicBondReset({bond_info}, T={self.T}, "
                f"enabled={self.enabled}, has_reset={self.has_reset})") 