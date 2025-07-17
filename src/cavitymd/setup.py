"""
Specialized setup classes for cavity molecular dynamics simulations.

This module provides organized setup functionality for forces, thermostats,
and particle management, following HOOMD-blue best practices.
"""

import logging
import numpy as np
from typing import List, Tuple, Dict, Optional, Any
import hoomd
import hoomd.md
import gsd.hoomd
from hoomd.bussi_reservoir.thermostats import BussiReservoir as Bussi

from .validation import CavitySimulationParams
from .utils import PhysicalConstants, unwrap_positions
from .forces import CavityForce
from .variants import StepVariant


def convert_to_numpy_array(array_data):
    """
    Convert array data to NumPy array, handling both CPU (NumPy) and GPU (CuPy) arrays.
    
    Parameters
    ----------
    array_data : array-like
        Array data from local snapshot (could be NumPy or CuPy)
        
    Returns
    -------
    np.ndarray
        NumPy array
    """
    if hasattr(array_data, 'get'):  # CuPy array (GPU)
        return array_data.get()
    else:  # NumPy array (CPU) or other array-like
        return np.asarray(array_data)


class DeviceSetup:
    """Handles HOOMD device configuration and setup."""
    
    @staticmethod
    def create_device(device_type: str, gpu_id: int = 0) -> hoomd.device.Device:
        """
        Create and configure HOOMD device.
        
        Parameters
        ----------
        device_type : str
            Device type ('CPU' or 'GPU')
        gpu_id : int
            GPU ID for GPU devices
            
        Returns
        -------
        hoomd.device.Device
            Configured HOOMD device
            
        Raises
        ------
        ValueError
            If device configuration is invalid
        """
        device_type = device_type.upper()
        
        if device_type == 'CPU':
            device = hoomd.device.CPU()
            logging.debug("Using CPU device")
        elif device_type == 'GPU':
            try:
                # Try different GPU initialization methods for different HOOMD versions
                try:
                    # First try with gpu_ids parameter (newer HOOMD versions)
                    device = hoomd.device.GPU(gpu_ids=[gpu_id])
                except TypeError:
                    # Fall back to older parameter name
                    device = hoomd.device.GPU(gpu_id=gpu_id)
                logging.info(f"Using GPU device {gpu_id}")
            except Exception as e:
                logging.error(f"Failed to initialize GPU {gpu_id}: {e}")
                logging.info("Falling back to CPU device")
                device = hoomd.device.CPU()
        else:
            raise ValueError(f"Invalid device type: {device_type}")
        
        return device
    
    @staticmethod
    def generate_simulation_seed(seed: Optional[int], replica: int) -> int:
        """
        Generate deterministic simulation seed.
        
        Parameters
        ----------
        seed : Optional[int]
            User-specified seed (if None, use replica-based)
        replica : int
            Replica number for deterministic seed generation
            
        Returns
        -------
        int
            Simulation seed
        """
        if seed is not None:
            logging.info(f"Using user-specified seed: {seed}")
            return seed
        else:
            replica_seed = hash(str(replica)) % (2**31)
            logging.info(f"Using replica-based seed: {replica_seed} (replica {replica})")
            return replica_seed


class ParticleSetup:
    """Handles cavity particle creation and validation."""
    
    @staticmethod
    def check_cavity_particle_exists(snapshot: Any) -> Tuple[bool, int]:
        """
        Check if cavity particle already exists in snapshot.
        
        Parameters
        ----------
        snapshot : Any
            GSD snapshot to check
            
        Returns
        -------
        Tuple[bool, int]
            (exists, count) where exists is True if cavity particles found,
            and count is the number of cavity particles
        """
        if 'L' not in snapshot.particles.types:
            return False, 0
        
        # Count cavity particles (type 'L', typeid = 2)
        cavity_count = np.sum(snapshot.particles.typeid == 2)
        return cavity_count > 0, cavity_count
    
    @staticmethod
    def create_cavity_particle(snapshot: Any, params: CavitySimulationParams) -> Any:
        """
        Add a new cavity particle to the simulation snapshot.
        
        Parameters
        ----------
        snapshot : Any
            Original snapshot
        params : CavitySimulationParams
            Simulation parameters
            
        Returns
        -------
        Any
            Modified snapshot with cavity particle
        """
        logging.debug("Creating new cavity particle and adding to system...")
        
        # Set numpy seed for reproducible cavity particle positioning
        if params.seed is not None:
            np.random.seed(params.seed + 2)
            logging.debug(f"Using seed {params.seed + 2} for cavity particle positioning")
        
        # Calculate dipole moment and photon position
        positions = unwrap_positions(snapshot.particles.position, 
                                   snapshot.particles.image, 
                                   snapshot.configuration.box[:3])
        dipmom = np.einsum('i,ij->j', snapshot.particles.charge, positions)
        
        # Get physical constants
        constants = params.get_physical_constants()
        omegac = constants['omegac_au']
        
        # Determine coupling strength for initial placement
        initial_couplstr = 0.0 if params.switch_time_ps is not None else params.coupling_strength
        
        # Calculate photon mass and spring constant
        phmass = 1.0  # Photon mass in atomic units
        K = phmass * omegac * omegac
        
        # Calculate new equilibrium position
        q_eq = -(initial_couplstr / K) * dipmom
        
        # Add some randomness to avoid exact zero position
        if np.allclose(q_eq, 0.0):
            q_eq = np.random.normal(0.0, 0.1, 3)
        
        newpos = q_eq
        
        # Wrap position and get image flags
        def wrap_position(x, L):
            image_flags = np.floor((x + L/2) / L)
            wrapped_position = x - image_flags * L
            return wrapped_position, image_flags.astype(int)
        
        newpos, image_flags = wrap_position(newpos, np.array(snapshot.configuration.box[:3]))
        
        # Add photon particle
        if 'L' not in snapshot.particles.types:
            snapshot.particles.types.append('L')
        
        snapshot.particles.N += 1
        snapshot.particles.typeid = np.append(snapshot.particles.typeid, [2])
        snapshot.particles.position = np.append(snapshot.particles.position, [newpos], axis=0)
        snapshot.particles.charge = np.append(snapshot.particles.charge, [0.0])
        snapshot.particles.mass = np.append(snapshot.particles.mass, [phmass])
        snapshot.particles.diameter = np.append(snapshot.particles.diameter, [1.0])
        snapshot.particles.image = np.vstack([snapshot.particles.image, image_flags])
        
        # Set additional particle properties
        if hasattr(snapshot.particles, "body"):
            snapshot.particles.body = np.append(snapshot.particles.body, [-1], axis=0)
        
        if hasattr(snapshot.particles, "orientation"):
            snapshot.particles.orientation = np.append(
                snapshot.particles.orientation, [[1.0, 0.0, 0.0, 0.0]], axis=0
            )
        
        if hasattr(snapshot.particles, "moment_inertia"):
            snapshot.particles.moment_inertia = np.vstack([
                snapshot.particles.moment_inertia, np.zeros((1, 3))
            ])
        
        if hasattr(snapshot.particles, "velocity"):
            snapshot.particles.velocity = np.vstack([
                snapshot.particles.velocity, np.zeros((1, 3))
            ])
        
        if hasattr(snapshot.particles, "angmom"):
            snapshot.particles.angmom = np.vstack([
                snapshot.particles.angmom, np.zeros((1, 4))
            ])
        
        logging.debug(f"Cavity particle added at position {newpos}")
        return snapshot
    
    @staticmethod
    def validate_cavity_particle(sim: hoomd.Simulation) -> None:
        """
        Validate that cavity particle exists and is properly configured.
        
        Parameters
        ----------
        sim : hoomd.Simulation
            HOOMD simulation object
            
        Raises
        ------
        ValueError
            If cavity particle validation fails
        """
        # Check particle types
        particle_types = sim.state.particle_types
        if 'L' not in particle_types:
            raise ValueError("Cavity particle type 'L' not found in simulation")
        
        # Check particle count - use device-aware local snapshot access
        # Detect device type from simulation
        device = sim.device
        if hasattr(device, 'gpu_ids') or 'GPU' in str(type(device)):
            local_snapshot = sim.state.gpu_local_snapshot
        else:
            local_snapshot = sim.state.cpu_local_snapshot
            
        with local_snapshot as snap:
            # Convert to numpy array for robust handling across CPU/GPU
            typeid_array = convert_to_numpy_array(snap.particles.typeid)
            
            if 2 not in typeid_array:
                raise ValueError("No cavity particles found in simulation")
            
            cavity_count = np.sum(typeid_array == 2)
            if cavity_count != 1:
                raise ValueError(f"Expected exactly 1 cavity particle, found {cavity_count}")
            
            cavity_indices = np.where(typeid_array == 2)[0]
            cavity_index = cavity_indices[0]
            cavity_position = snap.particles.position[cavity_index]
            logging.info(f"Cavity particle validated at index {cavity_index}, position {cavity_position}")


class ForceSetup:
    """Handles force field setup and configuration."""
    
    @staticmethod
    def create_cavity_force(params: CavitySimulationParams, time_tracker=None) -> Optional[CavityForce]:
        """
        Create cavity force with proper configuration.
        
        Parameters
        ----------
        params : CavitySimulationParams
            Simulation parameters
        time_tracker : Optional
            Time tracker for time-varying parameters
            
        Returns
        -------
        Optional[CavityForce]
            Configured cavity force (None if not cavity simulation)
        """
        if not params.incavity:
            return None
        
        constants = params.get_physical_constants()
        omegac = constants['omegac_au']
        
        # Create time-varying coupling and dissipation if switch_time is specified
        if params.switch_time_ps is not None:
            if time_tracker is None:
                raise ValueError("Time tracker required for time-varying parameters")
            
            # Create step variants for instantaneous switching
            coupling_variant = StepVariant(
                target_value=params.coupling_strength,
                switch_time_ps=params.switch_time_ps,
                time_tracker=time_tracker
            )
            
            dissipation_variant = StepVariant(
                target_value=params.dissipation,
                switch_time_ps=params.switch_time_ps,
                time_tracker=time_tracker
            )
            
            logging.info("Using time-varying cavity force:")
            logging.info(f"  Switch time: {params.switch_time_ps} ps")
            logging.info(f"  Coupling: 0 → {params.coupling_strength} a.u.")
            logging.info(f"  Dissipation: 0 → {params.dissipation} a.u.")
            
            # Create cavity force with variants
            cavityforce = CavityForce(
                kvector=np.array([0, 0, 1]),
                couplstr=coupling_variant,
                omegac=omegac,
                dissipation=dissipation_variant
            )
        else:
            # Use constant values
            logging.info("Using constant cavity force:")
            logging.info(f"  Coupling: {params.coupling_strength} a.u.")
            logging.info(f"  Dissipation: {params.dissipation} a.u.")
            
            cavityforce = CavityForce(
                kvector=np.array([0, 0, 1]),
                couplstr=params.coupling_strength,
                omegac=omegac,
                dissipation=params.dissipation
            )
        
        return cavityforce
    
    @staticmethod
    def create_harmonic_bonds() -> hoomd.md.bond.Harmonic:
        """
        Create harmonic bond forces for the molecular system.
        
        Returns
        -------
        hoomd.md.bond.Harmonic
            Configured harmonic bond force
        """
        harmonic = hoomd.md.bond.Harmonic()
        harmonic.params['O-O'] = dict(k=2*0.36602, r0=2.281655158)
        harmonic.params['N-N'] = dict(k=2*0.71625, r0=2.0743522177)
        logging.debug("Harmonic bond forces configured")
        return harmonic
    
    @staticmethod
    def create_lennard_jones_force(params: CavitySimulationParams, rcut: float = 15.0) -> hoomd.md.pair.LJ:
        """
        Create Lennard-Jones pair force.
        
        Parameters
        ----------
        params : CavitySimulationParams
            Simulation parameters
        rcut : float
            Cutoff distance
            
        Returns
        -------
        hoomd.md.pair.LJ
            Configured LJ force
        """
        # Create neighbor list
        cell = hoomd.md.nlist.Cell(buffer=1.0, exclusions=('bond',))
        
        # Create LJ force
        lj = hoomd.md.pair.LJ(nlist=cell, mode='shift')
        lj.params[('O', 'O')] = dict(epsilon=0.00016685201, sigma=6.230426584)
        lj.r_cut[('O', 'O')] = rcut
        lj.params[('N', 'N')] = dict(epsilon=0.000083426, sigma=5.48277488)
        lj.r_cut[('N', 'N')] = rcut
        lj.params[('N', 'O')] = dict(epsilon=0.00025027802, sigma=4.9832074319)
        lj.r_cut[('N', 'O')] = rcut
        
        # Disable pair interactions with cavity particle
        if params.incavity:
            for particle_type in ['N', 'O']:
                lj.params[('L', particle_type)] = dict(epsilon=0.0, sigma=1.0)
                lj.r_cut[('L', particle_type)] = 0.0
                lj.params[(particle_type, 'L')] = dict(epsilon=0.0, sigma=1.0)
                lj.r_cut[(particle_type, 'L')] = 0.0
            lj.params[('L', 'L')] = dict(epsilon=0.0, sigma=1.0)
            lj.r_cut[('L', 'L')] = 0.0
        
        logging.debug(f"Lennard-Jones forces configured with rcut = {rcut} a.u.")
        return lj
    
    @staticmethod
    def create_coulomb_forces(rcut: float = 15.0) -> Tuple[hoomd.md.long_range.pppm.Coulomb, hoomd.md.pair.Ewald]:
        """
        Create long-range Coulomb forces using PPPM method.
        
        Parameters
        ----------
        rcut : float
            Cutoff distance
            
        Returns
        -------
        Tuple[hoomd.md.long_range.pppm.Coulomb, hoomd.md.pair.Ewald]
            Short-range and long-range Coulomb forces
        """
        # Create neighbor list
        cell = hoomd.md.nlist.Cell(buffer=1.0, exclusions=('bond',))
        
        # Create PPPM forces
        numpoints = 32
        order = 6
        short, long = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
            nlist=cell,
            resolution=[numpoints, numpoints, numpoints],
            order=order,
            r_cut=rcut,
            alpha=0.0
        )
        
        logging.debug(f"PPPM Coulomb forces configured with rcut = {rcut} a.u.")
        return short, long
    
    @staticmethod
    def create_all_forces(params: CavitySimulationParams, time_tracker=None, rcut: float = 15.0) -> List[hoomd.md.force.Force]:
        """
        Create all forces for the simulation.
        
        Parameters
        ----------
        params : CavitySimulationParams
            Simulation parameters
        time_tracker : Optional
            Time tracker for time-varying parameters
        rcut : float
            Cutoff distance
            
        Returns
        -------
        List[hoomd.md.force.Force]
            List of all configured forces
        """
        forces = []
        
        # Add cavity force if requested
        if params.incavity:
            cavity_force = ForceSetup.create_cavity_force(params, time_tracker)
            if cavity_force is not None:
                forces.append(cavity_force)
        
        # Add harmonic bonds
        forces.append(ForceSetup.create_harmonic_bonds())
        
        # Add Lennard-Jones forces
        forces.append(ForceSetup.create_lennard_jones_force(params, rcut))
        
        # Add Coulomb forces
        short, long = ForceSetup.create_coulomb_forces(rcut)
        forces.extend([short, long])
        
        logging.debug(f"Total forces configured: {len(forces)}")
        return forces


class ThermostatSetup:
    """Handles thermostat configuration and setup."""
    
    @staticmethod
    def create_molecular_thermostat(params: CavitySimulationParams) -> hoomd.md.methods.Method:
        """
        Create molecular thermostat method.
        
        Parameters
        ----------
        params : CavitySimulationParams
            Simulation parameters
            
        Returns
        -------
        hoomd.md.methods.Method
            Configured molecular thermostat method
        """
        constants = params.get_physical_constants()
        kT = constants['kT']
        molecular_tau_au = constants['molecular_tau_au']
        
        molecular_filter = hoomd.filter.Type(['O', 'N'])
        
        if params.molecular_thermostat.lower() == 'bussi':
            logging.debug("Configuring molecular Bussi thermostat")
            molecular_bussi = Bussi(kT=kT, tau=molecular_tau_au)
            molecular_method = hoomd.md.methods.ConstantVolume(
                filter=molecular_filter, thermostat=molecular_bussi
            )
            logging.debug(f"Molecular Bussi: T={params.temperature:.1f} K, tau={params.molecular_thermostat_tau:.3f} ps")
            
        elif params.molecular_thermostat.lower() == 'langevin':
            logging.debug("Configuring molecular Langevin thermostat")
            molecular_gamma = PhysicalConstants.gamma_from_tau_ps(params.molecular_thermostat_tau)
            molecular_method = hoomd.md.methods.Langevin(
                filter=molecular_filter, kT=kT, 
                default_gamma=molecular_gamma, tally_reservoir_energy=True
            )
            logging.debug(f"Molecular Langevin: T={params.temperature:.1f} K, gamma={molecular_gamma:.6f} a.u.^-1")
            
        elif params.molecular_thermostat.lower() == 'none':
            logging.debug("Configuring molecular NVE (no thermostat)")
            molecular_method = hoomd.md.methods.ConstantVolume(filter=molecular_filter)
            
        else:
            raise ValueError(f"Invalid molecular thermostat: {params.molecular_thermostat}")
        
        return molecular_method
    
    @staticmethod
    def create_cavity_thermostat(params: CavitySimulationParams) -> Optional[hoomd.md.methods.Method]:
        """
        Create cavity thermostat method.
        
        Parameters
        ----------
        params : CavitySimulationParams
            Simulation parameters
            
        Returns
        -------
        Optional[hoomd.md.methods.Method]
            Configured cavity thermostat method (None if no cavity)
        """
        if not params.incavity:
            return None
        
        constants = params.get_physical_constants()
        kT = constants['kT']
        cavity_tau_au = constants['cavity_tau_au']
        
        cavity_filter = hoomd.filter.Type(['L'])
        
        if params.cavity_thermostat.lower() == 'bussi':
            logging.debug("Configuring cavity Bussi thermostat")
            cavity_bussi = Bussi(kT=kT, tau=cavity_tau_au)
            cavity_method = hoomd.md.methods.ConstantVolume(
                filter=cavity_filter, thermostat=cavity_bussi
            )
            logging.debug(f"Cavity Bussi: T={params.temperature:.1f} K, tau={params.cavity_thermostat_tau:.3f} ps")
            
        elif params.cavity_thermostat.lower() == 'langevin':
            logging.debug("Configuring cavity Langevin thermostat")
            base_gamma = PhysicalConstants.gamma_from_tau_ps(params.cavity_thermostat_tau)
            cavity_gamma = params.cavity_damping_factor * base_gamma
            cavity_method = hoomd.md.methods.Langevin(
                filter=cavity_filter, kT=kT,
                default_gamma=cavity_gamma, tally_reservoir_energy=True
            )
            logging.debug(f"Cavity Langevin: T={params.temperature:.1f} K, gamma={cavity_gamma:.6f} a.u.^-1")
            logging.debug(f"  Damping factor: {params.cavity_damping_factor:.1f}x")
            
        elif params.cavity_thermostat.lower() == 'none':
            logging.debug("Configuring cavity NVE (no thermostat)")
            cavity_method = hoomd.md.methods.ConstantVolume(filter=cavity_filter)
            
        else:
            raise ValueError(f"Invalid cavity thermostat: {params.cavity_thermostat}")
        
        return cavity_method
    
    @staticmethod
    def create_all_methods(params: CavitySimulationParams) -> Tuple[List[hoomd.md.methods.Method], Dict[str, Any]]:
        """
        Create all integration methods for the simulation.
        
        Parameters
        ----------
        params : CavitySimulationParams
            Simulation parameters
            
        Returns
        -------
        Tuple[List[hoomd.md.methods.Method], Dict[str, Any]]
            List of integration methods and thermostat references
        """
        methods = []
        thermostat_refs = {}
        
        # Create molecular method
        molecular_method = ThermostatSetup.create_molecular_thermostat(params)
        methods.append(molecular_method)
        
        # Store thermostat references for energy tracking
        if params.molecular_thermostat.lower() == 'bussi':
            thermostat_refs['molecular_bussi'] = molecular_method.thermostat
        elif params.molecular_thermostat.lower() == 'langevin':
            thermostat_refs['molecular_langevin'] = molecular_method
        
        # Create cavity method if needed
        if params.incavity:
            cavity_method = ThermostatSetup.create_cavity_thermostat(params)
            if cavity_method is not None:
                methods.append(cavity_method)
                
                # Store thermostat references
                if params.cavity_thermostat.lower() == 'bussi':
                    thermostat_refs['cavity_bussi'] = cavity_method.thermostat
                elif params.cavity_thermostat.lower() == 'langevin':
                    thermostat_refs['cavity_langevin'] = cavity_method
        
        logging.debug(f"Total integration methods: {len(methods)}")
        return methods, thermostat_refs


class ThermalizationSetup:
    """Handles particle velocity thermalization."""
    
    @staticmethod
    def thermalize_system(sim: hoomd.Simulation, params: CavitySimulationParams) -> None:
        """
        Initialize particle velocities using Maxwell-Boltzmann distribution.
        
        Parameters
        ----------
        sim : hoomd.Simulation
            HOOMD simulation object
        params : CavitySimulationParams
            Simulation parameters
        """
        if not params.restart_velocities:
            logging.debug("Velocity restart disabled - keeping existing velocities")
            return
        
        constants = params.get_physical_constants()
        kT = constants['kT']
        
        # Set numpy seed for reproducible thermalization
        if params.seed is not None:
            np.random.seed(params.seed + 1)
            logging.debug(f"Using seed {params.seed + 1} for thermalization")
        
        molecular_filter = hoomd.filter.Type(['O', 'N'])
        
        if params.incavity:
            # Thermalize molecular particles using HOOMD State
            sim.state.thermalize_particle_momenta(kT=kT, filter=molecular_filter)
            logging.debug("Thermalized molecular particles using HOOMD State")
            
            # Manual thermalization for cavity particle
            ThermalizationSetup._thermalize_cavity_particle_manually(sim, kT)
        else:
            # No cavity particle, thermalize all particles
            sim.state.thermalize_particle_momenta(kT=kT, filter=hoomd.filter.All())
            logging.debug("Thermalized all molecular particles")
        
        logging.debug("Thermalization completed")
    
    @staticmethod
    def _thermalize_cavity_particle_manually(sim: hoomd.Simulation, kT: float) -> None:
        """
        Manually thermalize cavity particle using Maxwell-Boltzmann distribution.
        
        Parameters
        ----------
        sim : hoomd.Simulation
            HOOMD simulation object
        kT : float
            Thermal energy in atomic units
        """
        # Use device-aware local snapshot access
        device = sim.device
        if hasattr(device, 'gpu_ids') or 'GPU' in str(type(device)):
            local_snapshot = sim.state.gpu_local_snapshot
        else:
            local_snapshot = sim.state.cpu_local_snapshot
            
        with local_snapshot as snap:
            # Convert to numpy array for robust handling across CPU/GPU
            typeid_array = convert_to_numpy_array(snap.particles.typeid)
            # Find cavity particle
            cavity_mask = typeid_array == 2
            if not np.any(cavity_mask):
                logging.warning("No cavity particle found for thermalization")
                return
            
            cavity_indices = np.where(cavity_mask)[0]
            
            cavity_idx = cavity_indices[0]
            cavity_mass = snap.particles.mass[cavity_idx]
            
            # Generate Maxwell-Boltzmann velocities
            sigma = np.sqrt(kT / cavity_mass)
            cavity_velocities = np.random.normal(0.0, sigma, 3)
            
            # Follow CuPy guidelines for CPU/GPU agnostic code
            # Get the appropriate array module based on snapshot arrays
            try:
                import cupy as cp
                xp = cp.get_array_module(snap.particles.velocity)
                HAS_CUPY_LOCAL = True
            except ImportError:
                xp = np
                HAS_CUPY_LOCAL = False
            
            # Set cavity particle velocities (HOOMD snapshots ALWAYS expect NumPy arrays)
            snap.particles.velocity[cavity_idx] = cavity_velocities
            
            logging.debug(f"Manually thermalized cavity particle: mass={cavity_mass:.3f}, sigma={sigma:.6f}")
            logging.debug(f"Cavity velocities: {cavity_velocities}")


class TimestepSetup:
    """Handles timestep calculation and configuration."""
    
    @staticmethod
    def calculate_initial_timestep(params: CavitySimulationParams) -> float:
        """
        Calculate initial timestep based on parameters.
        
        Parameters
        ----------
        params : CavitySimulationParams
            Simulation parameters
            
        Returns
        -------
        float
            Initial timestep in atomic units
        """
        if params.dt_fs is not None:
            # Use user-specified timestep
            dt_au = PhysicalConstants.fs_to_atomic_units(params.dt_fs)
            logging.debug(f"Using user-specified timestep: {params.dt_fs} fs = {dt_au:.6f} a.u.")
        else:
            # Use default timestep
            dt_ps = 0.0001  # 0.1 fs
            dt_au = PhysicalConstants.ps_to_atomic_units(dt_ps)
            logging.debug(f"Using default timestep: {dt_ps*1000:.1f} fs = {dt_au:.6f} a.u.")
        
        return dt_au
    
    @staticmethod
    def compute_optimal_timestep(sim: hoomd.Simulation, error_tolerance: float) -> float:
        """
        Compute optimal timestep based on forces and error tolerance.
        
        Parameters
        ----------
        sim : hoomd.Simulation
            HOOMD simulation object
        error_tolerance : float
            Error tolerance for timestep computation
            
        Returns
        -------
        float
            Optimal timestep in atomic units
        """
        # Initialize forces by running one step
        sim.run(1)
        
        # Collect forces and masses
        particle_data = sim.state.get_snapshot().particles
        masses = np.array(particle_data.mass)
        n_particles = len(masses)
        
        # Initialize total forces array
        total_forces = np.zeros((n_particles, 3))
        
        # Sum forces from all force objects
        for force in sim.operations.integrator.forces:
            if hasattr(force, 'forces'):
                particle_forces = force.forces
                if particle_forces is not None and len(particle_forces) == n_particles:
                    total_forces += np.asarray(particle_forces)
                else:
                    logging.debug(f"Force object {type(force).__name__} has no forces or wrong size")
            else:
                logging.debug(f"Force object {type(force).__name__} does not have forces attribute")
        
        # Calculate optimal timestep
        force_norm = np.array([np.linalg.norm(f) for f in total_forces])
        force_mass_sum = np.sum(force_norm / masses)
        
        if force_mass_sum > 0:
            optimal_dt = np.sqrt(error_tolerance / force_mass_sum)
            # Apply maximum timestep limit to prevent numerical instability
            max_dt = 5.0  # a.u. (~121 fs, appropriate for molecular systems)
            if optimal_dt > max_dt:
                logging.warning(f"Computed timestep {optimal_dt:.6f} a.u. exceeds maximum limit")
                logging.warning(f"Clamping to maximum timestep: {max_dt:.6f} a.u.")
                optimal_dt = max_dt
            logging.debug(f"Computed optimal timestep: {optimal_dt:.6f} a.u.")
            return optimal_dt
        else:
            logging.warning("Zero force detected, using default timestep")
            return PhysicalConstants.fs_to_atomic_units(0.1)  # 0.1 fs default 
    
    @staticmethod
    def compute_optimal_timestep_without_run(sim: hoomd.Simulation, error_tolerance: float) -> float:
        """
        Compute optimal timestep based on forces and error tolerance.
        
        This method assumes that forces have already been initialized by previous sim.run() calls.
        Use this when reservoir energy initialization and force initialization have already been done.
        
        Parameters
        ----------
        sim : hoomd.Simulation
            HOOMD simulation object
        error_tolerance : float
            Error tolerance for timestep computation
            
        Returns
        -------
        float
            Optimal timestep in atomic units
        """
        # Collect forces and masses (forces assumed already initialized)
        particle_data = sim.state.get_snapshot().particles
        masses = np.array(particle_data.mass)
        n_particles = len(masses)
        
        # Initialize total forces array
        total_forces = np.zeros((n_particles, 3))
        
        # Sum forces from all force objects
        for force in sim.operations.integrator.forces:
            if hasattr(force, 'forces'):
                particle_forces = force.forces
                if particle_forces is not None and len(particle_forces) == n_particles:
                    total_forces += np.asarray(particle_forces)
                else:
                    logging.debug(f"Force object {type(force).__name__} has no forces or wrong size")
            else:
                logging.debug(f"Force object {type(force).__name__} does not have forces attribute")
        
        # Calculate optimal timestep
        force_norm = np.array([np.linalg.norm(f) for f in total_forces])
        force_mass_sum = np.sum(force_norm / masses)
        
        if force_mass_sum > 0:
            optimal_dt = np.sqrt(error_tolerance / force_mass_sum)
            # Apply maximum timestep limit to prevent numerical instability
            max_dt = 5.0  # a.u. (~121 fs, appropriate for molecular systems)
            if optimal_dt > max_dt:
                logging.warning(f"Computed timestep {optimal_dt:.6f} a.u. exceeds maximum limit")
                logging.warning(f"Clamping to maximum timestep: {max_dt:.6f} a.u.")
                optimal_dt = max_dt
            logging.debug(f"Computed optimal timestep (without run): {optimal_dt:.6f} a.u.")
            return optimal_dt
        else:
            logging.warning("Zero force detected, using default timestep")
            return PhysicalConstants.fs_to_atomic_units(0.1)  # 0.1 fs default 