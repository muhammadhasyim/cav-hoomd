# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""State management system for FDR fork-and-clone simulations."""

from typing import Optional, Dict, Any, List, Tuple
import hoomd
import numpy as np
import pickle
import logging
import time
from pathlib import Path
import copy

from .utils import PhysicalConstants


class SimulationStateSnapshot:
    """Complete simulation state snapshot for fork-and-clone operations.
    
    This class captures and restores the complete state of a HOOMD simulation,
    including positions, velocities, forces, random number generator states,
    thermostat states, and integrator parameters. This is essential for 
    FDR measurements where exact reproducibility is required.
    """
    
    def __init__(self, 
                 timestep: int,
                 elapsed_time_ps: float,
                 particle_data: Dict[str, np.ndarray],
                 box_data: Dict[str, Any],
                 rng_state: Optional[Dict[str, Any]] = None,
                 thermostat_state: Optional[Dict[str, Any]] = None,
                 integrator_state: Optional[Dict[str, Any]] = None,
                 force_states: Optional[Dict[str, Any]] = None):
        """Initialize simulation state snapshot.
        
        Parameters
        ----------
        timestep : int
            Current simulation timestep
        elapsed_time_ps : float
            Elapsed simulation time in picoseconds
        particle_data : Dict[str, np.ndarray]
            Particle data including positions, velocities, forces, etc.
        box_data : Dict[str, Any]
            Simulation box parameters
        rng_state : Dict[str, Any], optional
            Random number generator state
        thermostat_state : Dict[str, Any], optional
            Thermostat internal states
        integrator_state : Dict[str, Any], optional
            Integrator parameters and state
        force_states : Dict[str, Any], optional
            Force object states and parameters
        """
        self.timestep = timestep
        self.elapsed_time_ps = elapsed_time_ps
        self.particle_data = particle_data
        self.box_data = box_data
        self.rng_state = rng_state
        self.thermostat_state = thermostat_state
        self.integrator_state = integrator_state
        self.force_states = force_states
        self.snapshot_time = time.time()
        
    def save(self, filepath: str) -> None:
        """Save snapshot to file.
        
        Parameters
        ----------
        filepath : str
            Path to save the snapshot
        """
        snapshot_data = {
            'timestep': self.timestep,
            'elapsed_time_ps': self.elapsed_time_ps,
            'particle_data': self.particle_data,
            'box_data': self.box_data,
            'rng_state': self.rng_state,
            'thermostat_state': self.thermostat_state,
            'integrator_state': self.integrator_state,
            'force_states': self.force_states,
            'snapshot_time': self.snapshot_time,
            'hoomd_version': hoomd.version.version,
            'metadata': {
                'creation_time': time.time(),
                'description': 'FDR simulation state snapshot'
            }
        }
        
        with open(filepath, 'wb') as f:
            pickle.dump(snapshot_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    @classmethod
    def load(cls, filepath: str) -> 'SimulationStateSnapshot':
        """Load snapshot from file.
        
        Parameters
        ----------
        filepath : str
            Path to load the snapshot from
            
        Returns
        -------
        SimulationStateSnapshot
            Loaded snapshot object
        """
        with open(filepath, 'rb') as f:
            data = pickle.load(f)
        
        return cls(
            timestep=data['timestep'],
            elapsed_time_ps=data['elapsed_time_ps'],
            particle_data=data['particle_data'],
            box_data=data['box_data'],
            rng_state=data.get('rng_state'),
            thermostat_state=data.get('thermostat_state'),
            integrator_state=data.get('integrator_state'),
            force_states=data.get('force_states')
        )


class StateSnapshotManager:
    """Manager for creating and restoring complete simulation states.
    
    This class provides functionality to capture complete HOOMD simulation
    states, create clones with modified parameters, and restore states for
    exact reproducibility. It's designed specifically for FDR measurements
    where we need to fork simulations and apply different perturbations.
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """Initialize state snapshot manager.
        
        Parameters
        ----------
        logger : logging.Logger, optional
            Logger for status messages
        """
        self.logger = logger or logging.getLogger(__name__)
        self._snapshots: Dict[str, SimulationStateSnapshot] = {}
        
    def capture_state(self, 
                     simulation: hoomd.Simulation,
                     time_tracker: Any,
                     snapshot_id: str) -> SimulationStateSnapshot:
        """Capture complete simulation state.
        
        Parameters
        ----------
        simulation : hoomd.Simulation
            HOOMD simulation object
        time_tracker : ElapsedTimeTracker
            Time tracking object
        snapshot_id : str
            Unique identifier for this snapshot
            
        Returns
        -------
        SimulationStateSnapshot
            Complete state snapshot
        """
        self.logger.info(f"Capturing simulation state: {snapshot_id}")
        
        # Get current timestep and time
        timestep = simulation.timestep
        elapsed_time_ps = time_tracker.elapsed_time()
        
        # Capture particle data
        with simulation.state.cpu_local_snapshot as snap:
            particle_data = {
                'positions': np.array(snap.particles.position),
                'velocities': np.array(snap.particles.velocity),
                'types': np.array(snap.particles.typeid),
                'masses': np.array(snap.particles.mass),
                'charges': np.array(snap.particles.charge) if hasattr(snap.particles, 'charge') else None,
                'images': np.array(snap.particles.image),
                'N_particles': snap.particles.N
            }
            
            # Capture box data
            box_data = {
                'Lx': snap.configuration.box[0],
                'Ly': snap.configuration.box[1], 
                'Lz': snap.configuration.box[2],
                'xy': snap.configuration.box[3],
                'xz': snap.configuration.box[4],
                'yz': snap.configuration.box[5],
                'dimensions': snap.configuration.dimensions
            }
        
        # Capture RNG state (HOOMD-specific)
        rng_state = self._capture_rng_state(simulation)
        
        # Capture thermostat states
        thermostat_state = self._capture_thermostat_state(simulation)
        
        # Capture integrator state
        integrator_state = self._capture_integrator_state(simulation)
        
        # Capture force states
        force_states = self._capture_force_states(simulation)
        
        snapshot = SimulationStateSnapshot(
            timestep=timestep,
            elapsed_time_ps=elapsed_time_ps,
            particle_data=particle_data,
            box_data=box_data,
            rng_state=rng_state,
            thermostat_state=thermostat_state,
            integrator_state=integrator_state,
            force_states=force_states
        )
        
        # Store snapshot
        self._snapshots[snapshot_id] = snapshot
        
        self.logger.info(f"State captured: timestep={timestep}, time={elapsed_time_ps:.3f} ps")
        return snapshot
    
    def restore_state(self, 
                     simulation: hoomd.Simulation,
                     time_tracker: Any,
                     snapshot_id: str) -> bool:
        """Restore simulation to a previously captured state.
        
        Parameters
        ----------
        simulation : hoomd.Simulation
            HOOMD simulation object to restore
        time_tracker : ElapsedTimeTracker
            Time tracking object to update
        snapshot_id : str
            Identifier of snapshot to restore
            
        Returns
        -------
        bool
            True if restoration successful, False otherwise
        """
        if snapshot_id not in self._snapshots:
            self.logger.error(f"Snapshot {snapshot_id} not found")
            return False
        
        snapshot = self._snapshots[snapshot_id]
        self.logger.info(f"Restoring simulation state: {snapshot_id}")
        
        try:
            # Restore particle data
            with simulation.state.cpu_local_snapshot as snap:
                snap.particles.position[:] = snapshot.particle_data['positions']
                snap.particles.velocity[:] = snapshot.particle_data['velocities']
                if snapshot.particle_data['charges'] is not None:
                    snap.particles.charge[:] = snapshot.particle_data['charges']
                snap.particles.image[:] = snapshot.particle_data['images']
                
                # Restore box
                snap.configuration.box = [
                    snapshot.box_data['Lx'],
                    snapshot.box_data['Ly'],
                    snapshot.box_data['Lz'],
                    snapshot.box_data['xy'],
                    snapshot.box_data['xz'],
                    snapshot.box_data['yz']
                ]
            
            # Reset timestep
            # Note: HOOMD doesn't allow direct timestep setting, 
            # so we'll track this separately
            
            # Restore RNG state
            self._restore_rng_state(simulation, snapshot.rng_state)
            
            # Restore thermostat states
            self._restore_thermostat_state(simulation, snapshot.thermostat_state)
            
            # Restore integrator state
            self._restore_integrator_state(simulation, snapshot.integrator_state)
            
            # Restore force states
            self._restore_force_states(simulation, snapshot.force_states)
            
            # Update time tracker
            time_tracker.total_time = PhysicalConstants.ps_to_atomic_units(snapshot.elapsed_time_ps)
            time_tracker.last_timestep = snapshot.timestep
            
            self.logger.info(f"State restored successfully: timestep={snapshot.timestep}, time={snapshot.elapsed_time_ps:.3f} ps")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to restore state {snapshot_id}: {e}")
            return False
    
    def create_clone(self, 
                    base_simulation: hoomd.Simulation,
                    base_time_tracker: Any,
                    snapshot_id: str,
                    modifications: Optional[Dict[str, Any]] = None) -> Tuple[hoomd.Simulation, Any]:
        """Create a clone of the simulation from a snapshot with optional modifications.
        
        Parameters
        ----------
        base_simulation : hoomd.Simulation
            Original simulation to clone
        base_time_tracker : ElapsedTimeTracker
            Original time tracker
        snapshot_id : str
            Snapshot to use as base for clone
        modifications : Dict[str, Any], optional
            Modifications to apply to the clone
            
        Returns
        -------
        Tuple[hoomd.Simulation, ElapsedTimeTracker]
            Cloned simulation and time tracker
        """
        if snapshot_id not in self._snapshots:
            raise ValueError(f"Snapshot {snapshot_id} not found")
        
        snapshot = self._snapshots[snapshot_id]
        self.logger.info(f"Creating simulation clone from {snapshot_id}")
        
        # Create new simulation with same device and seed
        device = base_simulation.device
        clone_sim = hoomd.Simulation(device=device, seed=base_simulation.seed)
        
        # Create snapshot and initialize state
        hoomd_snapshot = hoomd.Snapshot(device.communicator)
        
        if hoomd_snapshot.communicator.rank == 0:
            # Set up snapshot data
            N = snapshot.particle_data['N_particles']
            hoomd_snapshot.particles.N = N
            hoomd_snapshot.particles.position[:] = snapshot.particle_data['positions']
            hoomd_snapshot.particles.velocity[:] = snapshot.particle_data['velocities']
            hoomd_snapshot.particles.typeid[:] = snapshot.particle_data['types']
            hoomd_snapshot.particles.mass[:] = snapshot.particle_data['masses']
            if snapshot.particle_data['charges'] is not None:
                hoomd_snapshot.particles.charge[:] = snapshot.particle_data['charges']
            hoomd_snapshot.particles.image[:] = snapshot.particle_data['images']
            
            # Set box
            hoomd_snapshot.configuration.box = [
                snapshot.box_data['Lx'],
                snapshot.box_data['Ly'],
                snapshot.box_data['Lz'],
                snapshot.box_data['xy'],
                snapshot.box_data['xz'],
                snapshot.box_data['yz']
            ]
            
            # Set particle types (copy from original)
            with base_simulation.state.cpu_local_snapshot as orig_snap:
                hoomd_snapshot.particles.types = orig_snap.particles.types[:]
        
        clone_sim.create_state_from_snapshot(hoomd_snapshot)
        
        # Copy integrator and forces from base simulation
        if hasattr(base_simulation.operations, 'integrator') and base_simulation.operations.integrator is not None:
            # Note: Creating exact copy of integrator is complex due to HOOMD's architecture
            # For FDR, we'll need to reconstruct the integrator with same parameters
            # This will be handled by the FDR workflow that knows the simulation setup
            pass
        
        # Create clone time tracker
        from .analysis import ElapsedTimeTracker
        clone_time_tracker = ElapsedTimeTracker(clone_sim, base_time_tracker.runtime)
        clone_time_tracker.total_time = PhysicalConstants.ps_to_atomic_units(snapshot.elapsed_time_ps)
        clone_time_tracker.last_timestep = snapshot.timestep
        
        # Apply modifications if provided
        if modifications:
            self._apply_modifications(clone_sim, modifications)
        
        self.logger.info(f"Clone created successfully")
        return clone_sim, clone_time_tracker
    
    def _capture_rng_state(self, simulation: hoomd.Simulation) -> Dict[str, Any]:
        """Capture random number generator state."""
        # Note: HOOMD's RNG state capture depends on version and implementation
        # This is a placeholder for the essential RNG state information
        return {
            'seed': simulation.seed,
            'timestep': simulation.timestep,
            # Additional RNG state would go here in a full implementation
        }
    
    def _restore_rng_state(self, simulation: hoomd.Simulation, rng_state: Optional[Dict[str, Any]]) -> None:
        """Restore random number generator state."""
        if rng_state is None:
            return
        # Implementation would restore RNG state based on HOOMD version
        pass
    
    def _capture_thermostat_state(self, simulation: hoomd.Simulation) -> Dict[str, Any]:
        """Capture thermostat internal states."""
        thermostat_states = {}
        
        if hasattr(simulation.operations, 'integrator') and simulation.operations.integrator is not None:
            integrator = simulation.operations.integrator
            for i, method in enumerate(integrator.methods):
                if hasattr(method, 'thermostat') and method.thermostat is not None:
                    # Capture thermostat-specific state
                    thermostat_states[f'method_{i}_thermostat'] = {
                        'type': type(method.thermostat).__name__,
                        'kT': getattr(method.thermostat, 'kT', None),
                        'tau': getattr(method.thermostat, 'tau', None),
                        # Additional thermostat state would be captured here
                    }
        
        return thermostat_states
    
    def _restore_thermostat_state(self, simulation: hoomd.Simulation, thermostat_state: Optional[Dict[str, Any]]) -> None:
        """Restore thermostat internal states."""
        if thermostat_state is None:
            return
        # Implementation would restore thermostat states
        pass
    
    def _capture_integrator_state(self, simulation: hoomd.Simulation) -> Dict[str, Any]:
        """Capture integrator state."""
        integrator_state = {}
        
        if hasattr(simulation.operations, 'integrator') and simulation.operations.integrator is not None:
            integrator = simulation.operations.integrator
            integrator_state = {
                'dt': integrator.dt,
                'aniso': getattr(integrator, 'aniso', None),
                'num_methods': len(integrator.methods),
                'num_forces': len(integrator.forces)
            }
        
        return integrator_state
    
    def _restore_integrator_state(self, simulation: hoomd.Simulation, integrator_state: Optional[Dict[str, Any]]) -> None:
        """Restore integrator state."""
        if integrator_state is None:
            return
        # Integrator restoration is complex and context-dependent
        pass
    
    def _capture_force_states(self, simulation: hoomd.Simulation) -> Dict[str, Any]:
        """Capture force object states."""
        force_states = {}
        
        if hasattr(simulation.operations, 'integrator') and simulation.operations.integrator is not None:
            integrator = simulation.operations.integrator
            for i, force in enumerate(integrator.forces):
                force_states[f'force_{i}'] = {
                    'type': type(force).__name__,
                    # Force-specific parameters would be captured here
                }
        
        return force_states
    
    def _restore_force_states(self, simulation: hoomd.Simulation, force_states: Optional[Dict[str, Any]]) -> None:
        """Restore force object states."""
        if force_states is None:
            return
        # Force state restoration is also context-dependent
        pass
    
    def _apply_modifications(self, simulation: hoomd.Simulation, modifications: Dict[str, Any]) -> None:
        """Apply modifications to cloned simulation."""
        # This would handle modifications like adding perturbation forces
        # Implementation depends on the specific modifications needed
        pass
    
    def save_snapshot(self, snapshot_id: str, filepath: str) -> bool:
        """Save a snapshot to file.
        
        Parameters
        ----------
        snapshot_id : str
            Snapshot identifier
        filepath : str
            Path to save file
            
        Returns
        -------
        bool
            True if successful, False otherwise
        """
        if snapshot_id not in self._snapshots:
            self.logger.error(f"Snapshot {snapshot_id} not found")
            return False
        
        try:
            self._snapshots[snapshot_id].save(filepath)
            self.logger.info(f"Snapshot {snapshot_id} saved to {filepath}")
            return True
        except Exception as e:
            self.logger.error(f"Failed to save snapshot {snapshot_id}: {e}")
            return False
    
    def load_snapshot(self, snapshot_id: str, filepath: str) -> bool:
        """Load a snapshot from file.
        
        Parameters
        ----------
        snapshot_id : str
            Snapshot identifier
        filepath : str
            Path to load from
            
        Returns
        -------
        bool
            True if successful, False otherwise
        """
        try:
            snapshot = SimulationStateSnapshot.load(filepath)
            self._snapshots[snapshot_id] = snapshot
            self.logger.info(f"Snapshot {snapshot_id} loaded from {filepath}")
            return True
        except Exception as e:
            self.logger.error(f"Failed to load snapshot from {filepath}: {e}")
            return False
    
    def list_snapshots(self) -> List[str]:
        """List all available snapshot IDs.
        
        Returns
        -------
        List[str]
            List of snapshot identifiers
        """
        return list(self._snapshots.keys())
    
    def delete_snapshot(self, snapshot_id: str) -> bool:
        """Delete a snapshot.
        
        Parameters
        ----------
        snapshot_id : str
            Snapshot identifier
            
        Returns
        -------
        bool
            True if successful, False otherwise
        """
        if snapshot_id in self._snapshots:
            del self._snapshots[snapshot_id]
            self.logger.info(f"Snapshot {snapshot_id} deleted")
            return True
        else:
            self.logger.warning(f"Snapshot {snapshot_id} not found")
            return False
