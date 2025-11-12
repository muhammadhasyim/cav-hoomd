#!/usr/bin/env python3
"""
Unified Observable Data Output System (data.py)

Provides efficient HDF5-based output for all simulation observables with
support for concurrent read access (SWMR mode) enabling live monitoring.

This module replaces the old text-file-based trackers (EnergyTracker, 
TemperatureTracker, DiatomicMolecularTemperatures) with a unified HDF5 system.

Key Features:
- Single HDF5 file for all observables (energies, temperatures, order parameters)
- Hierarchical organization with groups
- SWMR mode for live read access during writing
- Efficient chunked storage with compression
- Thread-safe operations
- ~10-100x smaller file sizes vs text files
- Fast random access to any observable
"""

import h5py
import numpy as np
import hoomd
from typing import Dict, List, Optional, Any
from pathlib import Path
import threading


class ObservableWriter(hoomd.custom.Action):
    """
    Unified HDF5 writer for all simulation observables.
    
    Supports concurrent read access via SWMR mode, allowing real-time
    monitoring of simulation progress without interrupting the run.
    
    Organizes data hierarchically:
    - /time: simulation time array
    - /energies/: all energy components
    - /temperatures/: all temperature measurements
    - /order_parameters/dipole: dipole moment components
    - /order_parameters/density: density at wavevectors
    
    Parameters
    ----------
    output_file : str or Path
        HDF5 output file path
    time_tracker : ElapsedTimeTracker
        Time tracker for accurate timing
    output_period_ps : float
        Output period in picoseconds
    chunk_size : int, optional
        Number of data points per chunk (affects I/O performance). Default: 1000
    compression : str, optional
        Compression algorithm ('gzip', 'lzf', or None). Default: 'gzip'
    compression_level : int, optional
        Compression level (0-9 for gzip). Default: 4
    enable_swmr : bool, optional
        Enable SWMR mode for concurrent read access. Default: True
    flush_interval : int, optional
        Number of writes before flushing to disk. Default: 10
        
    Attributes
    ----------
    file : h5py.File
        Open HDF5 file handle
    datasets : dict
        Dictionary of HDF5 datasets for efficient access
    write_count : int
        Number of writes since last flush
        
    Methods
    -------
    add_energy_tracker(energy_tracker)
        Register energy tracker for output
    add_temperature_tracker(temperature_tracker)
        Register temperature tracker for output
    add_molecular_temperature_tracker(molecular_temp_tracker)
        Register molecular temperature tracker for output
    add_dipole_tracker(dipole_tracker)
        Register dipole moment tracker for output
    add_density_tracker(density_tracker)
        Register density correlation tracker for output
    act(timestep)
        Write data at each timestep (called by HOOMD)
    close()
        Close HDF5 file and finalize
        
    Examples
    --------
    >>> from cavitymd.data import ObservableWriter
    >>> 
    >>> # Create writer
    >>> writer = ObservableWriter(
    ...     output_file='observables.h5',
    ...     time_tracker=time_tracker,
    ...     output_period_ps=0.1,
    ...     enable_swmr=True
    ... )
    >>> 
    >>> # Register trackers
    >>> writer.add_energy_tracker(energy_tracker)
    >>> writer.add_temperature_tracker(temp_tracker)
    >>> writer.add_dipole_tracker(dipole_tracker)
    >>> 
    >>> # Add to simulation
    >>> updater = hoomd.update.CustomUpdater(
    ...     action=writer,
    ...     trigger=hoomd.trigger.Periodic(100)
    ... )
    >>> sim.operations.updaters.append(updater)
    >>> 
    >>> # Run simulation
    >>> sim.run(100000)
    >>> 
    >>> # Close writer
    >>> writer.close()
    
    Notes
    -----
    - **SWMR Mode**: Enables live monitoring with external scripts while simulation runs
    - **Read Example**: `with h5py.File('observables.h5', 'r', swmr=True) as f: ...`
    - **Performance**: Adjust chunk_size based on output frequency and memory
    - **Compression**: Use 'gzip' for best compression, 'lzf' for speed, None for fastest I/O
    
    See Also
    --------
    EnergyTracker : Energy component tracking
    TemperatureTracker : Temperature measurement tracking
    DipoleMomentFDRTracker : Dipole moment tracking
    DensityCorrelationTracker : Density correlation tracking
    """
    
    def __init__(self,
                 output_file: str,
                 time_tracker,
                 output_period_ps: float,
                 chunk_size: int = 1000,
                 compression: Optional[str] = 'gzip',
                 compression_level: int = 4,
                 enable_swmr: bool = True,
                 flush_interval: int = 10):
        
        super().__init__()
        
        self.output_file = Path(output_file)
        self.time_tracker = time_tracker
        self.output_period_ps = output_period_ps
        self.chunk_size = chunk_size
        self.compression = compression
        self.compression_level = compression_level if compression == 'gzip' else None
        self.enable_swmr = enable_swmr
        self.flush_interval = flush_interval
        
        # Trackers
        self.energy_tracker = None
        self.temperature_tracker = None
        self.molecular_temperature_tracker = None
        self.dipole_tracker = None
        self.density_tracker = None
        self.controllers = {}  # name -> controller object
        
        # State
        self.last_output_time = None  # Use None to mark first output
        self.write_count = 0
        self.current_index = 0
        self.lock = threading.Lock()
        
        # Datasets dictionary
        self.datasets = {}
        
        # Initialize HDF5 file
        self._initialize_hdf5_file()
        
        print(f"ObservableWriter initialized:")
        print(f"  Output file: {self.output_file}")
        print(f"  Output period: {self.output_period_ps} ps")
        print(f"  Chunk size: {self.chunk_size}")
        print(f"  Compression: {self.compression}")
        print(f"  SWMR mode: {self.enable_swmr}")
        print(f"  Flush interval: {self.flush_interval}")
    
    def _initialize_hdf5_file(self):
        """Initialize HDF5 file with hierarchical structure."""
        # Create file with latest HDF5 format for better SWMR support
        self.file = h5py.File(
            self.output_file,
            'w',
            libver='latest'
        )
        
        # Store metadata as attributes
        self.file.attrs['format_version'] = '1.0'
        self.file.attrs['output_period_ps'] = self.output_period_ps
        self.file.attrs['creation_time'] = np.bytes_(np.datetime64('now').astype(str), 'utf-8')
        
        # Create hierarchical group structure
        self.file.create_group('energies')
        self.file.create_group('temperatures')
        self.file['temperatures'].create_group('molecular')  # For molecular temperature decomposition
        self.file.create_group('controllers')
        self.file.create_group('order_parameters')
        self.file['order_parameters'].create_group('dipole')
        self.file['order_parameters'].create_group('density')
        
        # Create time dataset (shared by all observables)
        self._create_resizable_dataset(
            '/', 'time', 
            description='Simulation time in picoseconds',
            units='ps'
        )
        
        print(f"✓ HDF5 file structure created: {self.output_file}")
    
    def _create_resizable_dataset(self, group_path: str, name: str, 
                                  shape: tuple = (1,),
                                  description: str = '', 
                                  units: str = ''):
        """
        Create a resizable dataset with chunking and optional compression.
        
        Parameters
        ----------
        group_path : str
            HDF5 group path (e.g., '/energies', '/temperatures')
        name : str
            Dataset name
        shape : tuple
            Initial shape (use (1,) for scalar time series, (1, N) for vector)
        description : str
            Human-readable description
        units : str
            Physical units
        """
        group = self.file[group_path]
        
        # Determine maxshape (unlimited first dimension for time)
        if len(shape) == 1:
            maxshape = (None,)
            chunks = (self.chunk_size,)
        else:
            maxshape = (None,) + shape[1:]
            chunks = (self.chunk_size,) + shape[1:]
        
        # Create dataset
        ds = group.create_dataset(
            name,
            shape=shape,
            maxshape=maxshape,
            chunks=chunks,
            dtype=np.float64,
            compression=self.compression,
            compression_opts=self.compression_level
        )
        
        # Add metadata
        ds.attrs['description'] = description
        ds.attrs['units'] = units
        
        # Store reference
        full_path = f"{group_path}/{name}" if group_path != '/' else name
        self.datasets[full_path] = ds
        
        return ds
    
    def add_energy_tracker(self, energy_tracker):
        """
        Register energy tracker and create corresponding datasets.
        
        Parameters
        ----------
        energy_tracker : EnergyTracker
            Energy tracker instance
        """
        self.energy_tracker = energy_tracker
        
        # Create energy datasets
        energy_components = [
            ('harmonic', 'Harmonic bond potential energy', 'Hartree'),
            ('lj', 'Lennard-Jones potential energy', 'Hartree'),
            ('ewald_short', 'Short-range Coulomb energy', 'Hartree'),
            ('ewald_long', 'Long-range Coulomb energy', 'Hartree'),
            ('cavity_harmonic', 'Cavity harmonic potential energy', 'Hartree'),
            ('cavity_coupling', 'Cavity-molecule coupling energy', 'Hartree'),
            ('cavity_dipole_self', 'Dipole self-energy', 'Hartree'),
            ('cavity_total_potential', 'Total cavity potential energy', 'Hartree'),
            ('molecular_kinetic', 'Molecular kinetic energy', 'Hartree'),
            ('cavity_kinetic', 'Cavity kinetic energy', 'Hartree'),
            ('total_kinetic', 'Total kinetic energy', 'Hartree'),
            ('total_potential', 'Total potential energy', 'Hartree'),
            ('system_total', 'Total system energy (KE + PE)', 'Hartree'),
            ('molecular_reservoir', 'Molecular reservoir energy', 'Hartree'),
            ('cavity_reservoir', 'Cavity reservoir energy', 'Hartree'),
            ('total_reservoir', 'Total reservoir energy', 'Hartree'),
            ('universe_total', 'Universe total energy (conserved)', 'Hartree'),
            ('temperature', 'Kinetic temperature', 'K'),
        ]
        
        for name, desc, units in energy_components:
            self._create_resizable_dataset('/energies', name, description=desc, units=units)
        
        print(f"✓ Energy tracker registered ({len(energy_components)} components)")
    
    def add_temperature_tracker(self, temperature_tracker):
        """
        Register temperature tracker and create corresponding datasets.
        
        Parameters
        ----------
        temperature_tracker : TemperatureTracker
            Temperature tracker instance
        """
        self.temperature_tracker = temperature_tracker
        
        # Create temperature datasets
        temp_components = [
            ('kinetic', 'Kinetic temperature', 'K'),
            ('harmonic_fictive', 'Harmonic fictive temperature (empirical)', 'K'),
            ('lj_coul_fictive', 'LJ+Coulombic fictive temperature (empirical)', 'K'),
            ('cavity_bath', 'Cavity bath temperature', 'K'),
            ('molecular_bath', 'Molecular bath temperature', 'K'),
            ('harmonic_equipartition', 'Harmonic equipartition temperature', 'K'),
        ]
        
        for name, desc, units in temp_components:
            self._create_resizable_dataset('/temperatures', name, description=desc, units=units)
        
        print(f"✓ Temperature tracker registered ({len(temp_components)} components)")
        
        # If molecular tracking is enabled, create molecular temperature datasets
        if hasattr(temperature_tracker, 'track_molecular') and temperature_tracker.track_molecular:
            molecular_components = [
                ('translational', 'Translational temperature', 'K'),
                ('rotational', 'Rotational temperature', 'K'),
                ('vibrational', 'Vibrational temperature', 'K'),
                ('total_kinetic', 'Total kinetic temperature', 'K'),
                ('translational_O2', 'O2 translational temperature', 'K'),
                ('translational_N2', 'N2 translational temperature', 'K'),
                ('rotational_O2', 'O2 rotational temperature', 'K'),
                ('rotational_N2', 'N2 rotational temperature', 'K'),
                ('vibrational_O2', 'O2 vibrational temperature', 'K'),
                ('vibrational_N2', 'N2 vibrational temperature', 'K'),
            ]
            
            for name, desc, units in molecular_components:
                self._create_resizable_dataset('/temperatures/molecular', name, description=desc, units=units)
            
            print(f"✓ Molecular temperature tracking enabled ({len(molecular_components)} components)")
    
    def add_molecular_temperature_tracker(self, molecular_temp_tracker):
        """
        Register molecular temperature tracker and create corresponding datasets.
        
        Parameters
        ----------
        molecular_temp_tracker : DiatomicMolecularTemperatures
            Molecular temperature tracker instance
        """
        self.molecular_temperature_tracker = molecular_temp_tracker
        
        # Create molecular temperature datasets
        molecular_components = [
            ('translational', 'Translational temperature', 'K'),
            ('rotational', 'Rotational temperature', 'K'),
            ('vibrational', 'Vibrational temperature', 'K'),
            ('total_kinetic', 'Total kinetic temperature', 'K'),
            ('translational_O2', 'O2 translational temperature', 'K'),
            ('translational_N2', 'N2 translational temperature', 'K'),
            ('rotational_O2', 'O2 rotational temperature', 'K'),
            ('rotational_N2', 'N2 rotational temperature', 'K'),
            ('vibrational_O2', 'O2 vibrational temperature', 'K'),
            ('vibrational_N2', 'N2 vibrational temperature', 'K'),
        ]
        
        for name, desc, units in molecular_components:
            self._create_resizable_dataset('/temperatures/molecular', name, description=desc, units=units)
        
        print(f"✓ Molecular temperature tracker registered ({len(molecular_components)} components)")
    
    def add_dipole_tracker(self, dipole_tracker):
        """
        Register dipole moment tracker and create corresponding datasets.
        
        Parameters
        ----------
        dipole_tracker : DipoleMomentFDRTracker
            Dipole moment tracker instance
        """
        self.dipole_tracker = dipole_tracker
        
        # Create dipole datasets (3D vector + magnitude)
        self._create_resizable_dataset(
            '/order_parameters/dipole', 'components',
            shape=(1, 3),
            description='Dipole moment components (x, y, z)',
            units='e·Å'
        )
        self._create_resizable_dataset(
            '/order_parameters/dipole', 'magnitude',
            description='Dipole moment magnitude',
            units='e·Å'
        )
        
        print(f"✓ Dipole moment tracker registered")
    
    def add_density_tracker(self, density_tracker, num_wavevectors: int):
        """
        Register density correlation tracker and create corresponding datasets.
        
        Parameters
        ----------
        density_tracker : DensityCorrelationTracker
            Density correlation tracker instance
        num_wavevectors : int
            Number of wavevectors being tracked
        """
        self.density_tracker = density_tracker
        self.num_wavevectors = num_wavevectors
        
        # Create density datasets (complex array for each wavevector)
        self._create_resizable_dataset(
            '/order_parameters/density', 'wavevectors',
            shape=(num_wavevectors, 3),
            description='Wavevector list (kx, ky, kz)',
            units='Å⁻¹'
        )
        self._create_resizable_dataset(
            '/order_parameters/density', 'rho_k_real',
            shape=(1, num_wavevectors),
            description='Real part of density at wavevectors',
            units='dimensionless'
        )
        self._create_resizable_dataset(
            '/order_parameters/density', 'rho_k_imag',
            shape=(1, num_wavevectors),
            description='Imaginary part of density at wavevectors',
            units='dimensionless'
        )
        
        # Write wavevectors (static, only once)
        if hasattr(density_tracker, 'wavevectors'):
            self.datasets['/order_parameters/density/wavevectors'][:] = density_tracker.wavevectors
        
        print(f"✓ Density tracker registered ({num_wavevectors} wavevectors)")
    
    def add_controller(self, name: str, controller):
        """
        Register a controller for HDF5 output.
        
        Parameters
        ----------
        name : str
            Controller name (e.g., 'diffeq', 'simple_setpoint')
        controller : object
            Controller instance with control_history attribute
        """
        self.controllers[name] = controller
        print(f"✓ Controller '{name}' registered for HDF5 output")
    
    def act(self, timestep):
        """
        Write all observable data to HDF5 file.
        
        Called by HOOMD at each timestep according to the trigger.
        """
        current_time_ps = self.time_tracker.elapsed_time
        
        # DEBUG: Log first few calls
        if self.write_count < 3 or timestep % 100000 == 0:
            print(f"[HDF5 Writer DEBUG] timestep={timestep}, current_time_ps={current_time_ps:.3f}, last_output_time={self.last_output_time}, write_count={self.write_count}")
        
        # Check if it's time to output
        if not self._should_output(current_time_ps):
            if self.write_count < 3:
                print(f"[HDF5 Writer DEBUG] Skipping output: not time yet (current={current_time_ps:.3f}, last={self.last_output_time}, period={self.output_period_ps})")
            return
        
        with self.lock:
            # DEBUG: Confirm we're writing
            if self.write_count < 3:
                print(f"[HDF5 Writer DEBUG] Writing data point {self.current_index} at time {current_time_ps:.3f} ps")
            
            # Resize all datasets
            new_size = self.current_index + 1
            for ds_path, ds in self.datasets.items():
                if ds_path == '/order_parameters/density/wavevectors':
                    continue  # Skip static wavevector dataset
                if len(ds.shape) == 1:
                    ds.resize((new_size,))
                else:
                    ds.resize((new_size,) + ds.shape[1:])
            
            # Write time
            self.datasets['time'][self.current_index] = current_time_ps
            
            # Write energy data
            if self.energy_tracker is not None:
                self._write_energy_data(self.current_index)
            
            # Write temperature data
            if self.temperature_tracker is not None:
                self._write_temperature_data(self.current_index)
            
            # Write molecular temperature data
            if self.molecular_temperature_tracker is not None:
                self._write_molecular_temperature_data(self.current_index)
            
            # Write dipole data
            if self.dipole_tracker is not None:
                self._write_dipole_data(self.current_index)
            
            # Write density data
            if self.density_tracker is not None:
                self._write_density_data(self.current_index)
            
            # Write controller data
            for controller_name, controller in self.controllers.items():
                self._write_controller_data(self.current_index, controller_name, controller)
            
            # Increment index and write count
            self.current_index += 1
            self.write_count += 1
            self.last_output_time = current_time_ps
            
            # Periodic flush for SWMR and data safety
            if self.write_count % self.flush_interval == 0:
                self.file.flush()
                # Enable SWMR mode after first flush (if not already enabled)
                if self.enable_swmr and not self.file.swmr_mode:
                    try:
                        self.file.swmr_mode = True
                        print("✓ SWMR mode enabled - file can now be read concurrently")
                    except Exception as e:
                        print(f"Warning: Could not enable SWMR mode: {e}")
    
    def _write_energy_data(self, index: int):
        """Write energy components at given index."""
        et = self.energy_tracker
        
        energy_map = {
            '/energies/harmonic': et.current_harmonic_energy,
            '/energies/lj': et.current_lj_energy,
            '/energies/ewald_short': et.current_ewald_short_energy,
            '/energies/ewald_long': et.current_ewald_long_energy,
            '/energies/cavity_harmonic': et.current_cavity_harmonic_energy,
            '/energies/cavity_coupling': et.current_cavity_coupling_energy,
            '/energies/cavity_dipole_self': et.current_cavity_dipole_self_energy,
            '/energies/cavity_total_potential': et.current_cavity_total_potential_energy,
            '/energies/molecular_kinetic': et.current_molecular_kinetic_energy,
            '/energies/cavity_kinetic': et.current_cavity_kinetic_energy,
            '/energies/total_kinetic': et.current_total_kinetic_energy,
            '/energies/total_potential': et.current_total_potential_energy,
            '/energies/system_total': et.current_system_total_energy,
            '/energies/molecular_reservoir': et.current_molecular_reservoir_energy,
            '/energies/cavity_reservoir': et.current_cavity_reservoir_energy,
            '/energies/total_reservoir': et.current_total_reservoir_energy,
            '/energies/universe_total': et.current_universe_total_energy,
            '/energies/temperature': et.current_temperature,
        }
        
        for ds_path, value in energy_map.items():
            if ds_path in self.datasets:
                self.datasets[ds_path][index] = value
    
    def _write_temperature_data(self, index: int):
        """Write temperature components at given index."""
        tt = self.temperature_tracker
        
        temp_map = {
            '/temperatures/kinetic': tt.kinetic_temperature,
            '/temperatures/harmonic_fictive': tt.harmonic_fictive_temperature,
            '/temperatures/lj_coul_fictive': tt.lj_coulombic_fictive_temperature,
            '/temperatures/cavity_bath': tt.cavity_bath_temperature,
            '/temperatures/molecular_bath': tt.molecular_bath_temperature,
            '/temperatures/harmonic_equipartition': tt.harmonic_equipartition_temperature,
        }
        
        for ds_path, value in temp_map.items():
            if ds_path in self.datasets:
                self.datasets[ds_path][index] = value
        
        # Write molecular temperatures if tracking is enabled
        if hasattr(tt, 'track_molecular') and tt.track_molecular:
            mol_temps = tt._calculate_molecular_temperatures()
            
            molecular_map = {
                '/temperatures/molecular/translational': mol_temps.get('T_trans', 0.0),
                '/temperatures/molecular/rotational': mol_temps.get('T_rot', 0.0),
                '/temperatures/molecular/vibrational': mol_temps.get('T_vib', 0.0),
                '/temperatures/molecular/total_kinetic': mol_temps.get('T_kinetic_total', 0.0),
                '/temperatures/molecular/translational_O2': mol_temps.get('T_trans_O2', 0.0),
                '/temperatures/molecular/translational_N2': mol_temps.get('T_trans_N2', 0.0),
                '/temperatures/molecular/rotational_O2': mol_temps.get('T_rot_O2', 0.0),
                '/temperatures/molecular/rotational_N2': mol_temps.get('T_rot_N2', 0.0),
                '/temperatures/molecular/vibrational_O2': mol_temps.get('T_vib_O2', 0.0),
                '/temperatures/molecular/vibrational_N2': mol_temps.get('T_vib_N2', 0.0),
            }
            
            for ds_path, value in molecular_map.items():
                if ds_path in self.datasets:
                    self.datasets[ds_path][index] = value
    
    def _write_molecular_temperature_data(self, index: int):
        """Write molecular temperature components at given index."""
        mt = self.molecular_temperature_tracker
        
        # Get current temperatures
        temps = mt._calculate_molecular_temperatures() if hasattr(mt, '_calculate_molecular_temperatures') else {}
        
        molecular_map = {
            '/temperatures/molecular/translational': temps.get('T_trans', 0.0),
            '/temperatures/molecular/rotational': temps.get('T_rot', 0.0),
            '/temperatures/molecular/vibrational': temps.get('T_vib', 0.0),
            '/temperatures/molecular/total_kinetic': temps.get('T_kinetic_total', 0.0),
            '/temperatures/molecular/translational_O2': temps.get('T_trans_O2', 0.0),
            '/temperatures/molecular/translational_N2': temps.get('T_trans_N2', 0.0),
            '/temperatures/molecular/rotational_O2': temps.get('T_rot_O2', 0.0),
            '/temperatures/molecular/rotational_N2': temps.get('T_rot_N2', 0.0),
            '/temperatures/molecular/vibrational_O2': temps.get('T_vib_O2', 0.0),
            '/temperatures/molecular/vibrational_N2': temps.get('T_vib_N2', 0.0),
        }
        
        for ds_path, value in molecular_map.items():
            if ds_path in self.datasets:
                self.datasets[ds_path][index] = value
    
    def _write_dipole_data(self, index: int):
        """Write dipole moment at given index."""
        # Get dipole moment from tracker
        if hasattr(self.dipole_tracker, '_calculate_total_dipole_moment'):
            dipole = self.dipole_tracker._calculate_total_dipole_moment()
            magnitude = np.linalg.norm(dipole)
            
            self.datasets['/order_parameters/dipole/components'][index, :] = dipole
            self.datasets['/order_parameters/dipole/magnitude'][index] = magnitude
    
    def _write_density_data(self, index: int):
        """Write density at wavevectors at given index."""
        # Get density data from tracker
        if hasattr(self.density_tracker, '_compute_density_at_wavevectors'):
            rho_k = self.density_tracker._compute_density_at_wavevectors()
            
            self.datasets['/order_parameters/density/rho_k_real'][index, :] = np.real(rho_k)
            self.datasets['/order_parameters/density/rho_k_imag'][index, :] = np.imag(rho_k)
    
    def _write_controller_data(self, index: int, controller_name: str, controller):
        """Write controller data at given index."""
        if not hasattr(controller, 'control_history') or len(controller.control_history) == 0:
            return  # No data to write
        
        # Get the latest control history entry
        latest_entry = controller.control_history[-1]
        
        # Create controller subgroup if it doesn't exist
        controller_group_path = f'/controllers/{controller_name}'
        if controller_group_path not in self.file:
            self.file.create_group(controller_group_path)
        
        # Write all fields from the history entry
        for key, value in latest_entry.items():
            if value is None:
                continue  # Skip None values
            
            ds_path = f'{controller_group_path}/{key}'
            if ds_path not in self.datasets:
                # Create new dataset
                self._create_resizable_dataset(
                    controller_group_path, key,
                    description=f'{controller_name} {key}',
                    units=''
                )
            
            # Write value
            self.datasets[ds_path][index] = value
    
    def _should_output(self, current_time_ps: float) -> bool:
        """Check if we should output data."""
        # Always output the first time
        if self.last_output_time is None:
            return True
        # Then check time difference
        return (current_time_ps - self.last_output_time) >= self.output_period_ps
    
    def close(self):
        """Close HDF5 file and finalize."""
        if hasattr(self, 'file') and self.file is not None:
            self.file.flush()
            self.file.close()
            print(f"✓ HDF5 file closed: {self.output_file}")
            print(f"  Total data points written: {self.current_index}")
    
    def __del__(self):
        """Cleanup on deletion."""
        self.close()


# Convenience function for reading HDF5 files with SWMR
def read_observables_live(hdf5_file: str) -> h5py.File:
    """
    Open HDF5 observables file for live reading (SWMR mode).
    
    Parameters
    ----------
    hdf5_file : str
        Path to HDF5 observables file
        
    Returns
    -------
    h5py.File
        Open HDF5 file handle in SWMR read mode
        
    Examples
    --------
    >>> import h5py
    >>> from cavitymd.data import read_observables_live
    >>> 
    >>> # Open file for live reading
    >>> with read_observables_live('observables.h5') as f:
    ...     time = f['time'][:]
    ...     temp = f['temperatures/kinetic'][:]
    ...     energy = f['energies/system_total'][:]
    ...     print(f"Current time: {time[-1]:.2f} ps")
    ...     print(f"Current temperature: {temp[-1]:.2f} K")
    """
    return h5py.File(hdf5_file, 'r', swmr=True)

