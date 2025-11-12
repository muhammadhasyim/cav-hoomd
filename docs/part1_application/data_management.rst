==================
Data Management
==================

Modern HDF5-based data output system for efficient, scalable observable tracking.

.. contents:: In this Section
   :local:
   :depth: 2

Overview
========

Cavity HOOMD provides a comprehensive **HDF5-based data management system** that replaces the old text-file-based trackers with a unified, efficient solution.

**Key Benefits:**

- **10-100x smaller files** compared to text output
- **Live monitoring** via SWMR (Single Writer Multiple Reader) mode
- **Hierarchical organization** with logical grouping
- **Fast random access** to any observable
- **Thread-safe operations** for concurrent access
- **Automatic compression** with configurable levels
- **Single file** for all observables

Migration from Text Files
==========================

**Old Approach (Text Files):**

.. code-block:: python

   # Multiple separate files
   energy_tracker = EnergyTracker(output_file='energy.txt')
   temp_tracker = TemperatureTracker(output_file='temperature.txt')
   dipole_tracker = DipoleTracker(output_file='dipole.txt')
   
   # Large files, slow I/O, no live access

**New Approach (HDF5):**

.. code-block:: python

   from cavitymd.data import ObservableWriter
   
   # Single unified file
   writer = ObservableWriter(
       output_file='observables.h5',
       time_tracker=time_tracker,
       output_period_ps=0.1,
       enable_swmr=True  # Enable live monitoring
   )
   
   # Register all trackers
   writer.add_energy_tracker(energy_tracker)
   writer.add_temperature_tracker(temp_tracker)
   writer.add_dipole_tracker(dipole_tracker)
   
   # Automatically writes to single HDF5 file

Using ObservableWriter
=======================

Basic Setup
-----------

.. code-block:: python

   from cavitymd.data import ObservableWriter
   from cavitymd.analysis import ElapsedTimeTracker, EnergyTracker, TemperatureTracker
   
   # Setup time tracking
   time_tracker = ElapsedTimeTracker(sim)
   
   # Create writer with SWMR mode
   writer = ObservableWriter(
       output_file='sim_data.h5',
       time_tracker=time_tracker,
       output_period_ps=0.1,        # Write every 0.1 ps
       chunk_size=1000,              # 1000 points per chunk
       compression='gzip',           # Enable compression
       compression_level=4,          # Moderate compression
       enable_swmr=True,             # Allow concurrent reads
       flush_interval=10             # Flush every 10 writes
   )
   
   # Register trackers
   energy_tracker = EnergyTracker(sim, forces={'cavity': cavity_force})
   temp_tracker = TemperatureTracker(sim, time_tracker)
   
   writer.add_energy_tracker(energy_tracker)
   writer.add_temperature_tracker(temp_tracker)
   
   # Add to simulation
   sim.operations.updaters.append(hoomd.update.CustomUpdater(
       action=writer,
       trigger=hoomd.trigger.Periodic(1)
   ))
   
   # Run simulation
   sim.run(1000000)
   
   # Close file
   writer.close()

File Organization
-----------------

HDF5 files are organized hierarchically:

.. code-block:: text

   observables.h5
   ├── time [dataset]                      # Simulation time (ps)
   ├── timestep [dataset]                  # HOOMD timestep
   │
   ├── /energies                           # Energy group
   │   ├── kinetic [dataset]               # Molecular kinetic
   │   ├── potential [dataset]             # Molecular potential
   │   ├── cavity_photon [dataset]         # Photon energy
   │   ├── interaction [dataset]           # Molecule-cavity
   │   ├── reservoir [dataset]             # Thermostat reservoir
   │   └── total [dataset]                 # Total energy
   │
   ├── /temperatures                       # Temperature group
   │   ├── kinetic [dataset]               # Kinetic temperature
   │   ├── translational [dataset]         # Translational
   │   ├── rotational [dataset]            # Rotational
   │   ├── vibrational [dataset]           # Vibrational
   │   ├── harmonic_fictive [dataset]      # Harmonic fictive
   │   └── lj_coulombic [dataset]          # LJ+Coulombic fictive
   │
   └── /order_parameters                   # Order parameters group
       ├── /dipole                         # Dipole group
       │   ├── x [dataset]                 # X component
       │   ├── y [dataset]                 # Y component
       │   ├── z [dataset]                 # Z component
       │   └── magnitude [dataset]         # |μ|
       │
       └── /density                        # Density correlation
           ├── Fkt_real [dataset]          # Real part
           └── Fkt_imag [dataset]          # Imaginary part

Reading Data
============

Real-Time Monitoring
--------------------

**Read data while simulation runs:**

.. code-block:: python

   from cavitymd.data import ObservableReader
   import matplotlib.pyplot as plt
   
   # Open file in read mode (simulation still running)
   reader = ObservableReader('observables.h5', swmr_mode=True)
   
   # Get latest data
   time = reader.get_time()
   energy = reader.get_observable('energies/total')
   
   # Plot live data
   plt.plot(time, energy)
   plt.xlabel('Time (ps)')
   plt.ylabel('Total Energy (a.u.)')
   plt.show()
   
   reader.close()

**Automatic refresh for live updates:**

.. code-block:: python

   import time
   
   reader = ObservableReader('observables.h5', swmr_mode=True)
   
   while simulation_running:
       # Get latest data
       time_data = reader.get_time()
       temp = reader.get_observable('temperatures/kinetic')
       
       # Update plot
       update_live_plot(time_data, temp)
       
       # Wait before next update
       time.sleep(1.0)
   
   reader.close()

Post-Processing
---------------

**Load complete dataset after simulation:**

.. code-block:: python

   from cavitymd.data import ObservableReader
   import numpy as np
   
   reader = ObservableReader('observables.h5')
   
   # Get all time points
   time = reader.get_time()
   
   # Get all observables
   kinetic_energy = reader.get_observable('energies/kinetic')
   potential_energy = reader.get_observable('energies/potential')
   temperature = reader.get_observable('temperatures/kinetic')
   
   # Compute derived quantities
   total_energy = kinetic_energy + potential_energy
   
   # Statistical analysis
   mean_temp = np.mean(temperature)
   std_temp = np.std(temperature)
   
   print(f"Mean temperature: {mean_temp:.2f} ± {std_temp:.2f} K")
   
   reader.close()

**Time-windowed analysis:**

.. code-block:: python

   reader = ObservableReader('observables.h5')
   
   # Get time array
   time = reader.get_time()
   
   # Find indices for time window
   mask = (time >= 10.0) & (time <= 50.0)  # 10-50 ps
   
   # Get data in window
   energy = reader.get_observable('energies/total')
   energy_window = energy[mask]
   
   # Analyze equilibrated region
   mean_energy = np.mean(energy_window)
   
   reader.close()

Direct HDF5 Access
------------------

**For advanced users:**

.. code-block:: python

   import h5py
   
   with h5py.File('observables.h5', 'r') as f:
       # Browse structure
       print("Available groups:", list(f.keys()))
       print("Energy datasets:", list(f['energies'].keys()))
       
       # Direct access
       time = f['time'][:]
       temperature = f['temperatures/kinetic'][:]
       
       # Access attributes
       units = f['energies/kinetic'].attrs['units']
       description = f['energies/kinetic'].attrs['description']

Performance Optimization
========================

Chunk Size Selection
--------------------

**Chunk size affects I/O performance:**

.. code-block:: python

   # Small systems, frequent output
   writer = ObservableWriter(
       ...,
       chunk_size=500,      # Smaller chunks
       output_period_ps=0.01  # Frequent writes
   )
   
   # Large systems, infrequent output
   writer = ObservableWriter(
       ...,
       chunk_size=5000,     # Larger chunks
       output_period_ps=1.0   # Less frequent writes
   )

**General guidelines:**

- Chunk size ≈ 1000-10000 data points
- Match to expected dataset size
- Balance between I/O efficiency and memory

Compression Settings
--------------------

**Trade-off between file size and speed:**

.. code-block:: python

   # Maximum compression (slower, smallest files)
   writer = ObservableWriter(
       ...,
       compression='gzip',
       compression_level=9
   )
   
   # Balanced compression (recommended)
   writer = ObservableWriter(
       ...,
       compression='gzip',
       compression_level=4
   )
   
   # Fast compression
   writer = ObservableWriter(
       ...,
       compression='lzf'  # Fast but less compression
   )
   
   # No compression (fastest, largest files)
   writer = ObservableWriter(
       ...,
       compression=None
   )

**Typical compression ratios:**

- Energy data: 5-10x reduction
- Temperature data: 5-10x reduction
- Trajectory data: 50-100x reduction

Flush Control
-------------

**Control write frequency:**

.. code-block:: python

   # Frequent flushing (more disk I/O, safer)
   writer = ObservableWriter(
       ...,
       flush_interval=5  # Flush every 5 writes
   )
   
   # Infrequent flushing (less disk I/O, faster)
   writer = ObservableWriter(
       ...,
       flush_interval=50  # Flush every 50 writes
   )

**Guidelines:**

- Short simulations: flush_interval = 5-10
- Long simulations: flush_interval = 20-50
- Higher values = better performance, more risk if crash

Best Practices
==============

File Management
---------------

**1. Use descriptive filenames:**

.. code-block:: python

   filename = f'sim_T{temperature}K_lambda{coupling}_replica{replica}.h5'
   writer = ObservableWriter(output_file=filename, ...)

**2. Organize by experiment:**

.. code-block:: python

   from pathlib import Path
   
   output_dir = Path('experiment_01/observables')
   output_dir.mkdir(parents=True, exist_ok=True)
   
   output_file = output_dir / f'replica_{replica}.h5'
   writer = ObservableWriter(output_file=output_file, ...)

**3. Always close files:**

.. code-block:: python

   try:
       writer = ObservableWriter(...)
       sim.run(1000000)
   finally:
       writer.close()  # Ensures data is flushed

Error Handling
--------------

**Robust file handling:**

.. code-block:: python

   from pathlib import Path
   
   output_file = Path('observables.h5')
   
   # Check if file exists
   if output_file.exists():
       # Backup or remove
       backup = output_file.with_suffix('.h5.bak')
       output_file.rename(backup)
   
   # Create new file
   writer = ObservableWriter(output_file=output_file, ...)

Metadata
--------

**Store simulation metadata:**

.. code-block:: python

   import h5py
   
   # After closing writer, add metadata
   with h5py.File('observables.h5', 'a') as f:
       f.attrs['temperature'] = 100.0
       f.attrs['coupling'] = 0.001
       f.attrs['frequency'] = 2000.0
       f.attrs['n_particles'] = 1000
       f.attrs['runtime_ps'] = 1000.0
       f.attrs['timestep'] = 0.001
       f.attrs['date'] = str(datetime.now())

Migration Example
=================

**Complete migration from old to new system:**

.. code-block:: python

   # OLD CODE (text-based)
   # energy_tracker = EnergyTracker(sim, output_file='energy.txt')
   # temp_tracker = TemperatureTracker(sim, output_file='temperature.txt')
   # sim.operations.updaters.append(energy_tracker)
   # sim.operations.updaters.append(temp_tracker)
   
   # NEW CODE (HDF5-based)
   from cavitymd.data import ObservableWriter
   
   # Create unified writer
   writer = ObservableWriter(
       output_file='observables.h5',
       time_tracker=time_tracker,
       output_period_ps=0.1,
       enable_swmr=True
   )
   
   # Register trackers (no need to specify output files)
   energy_tracker = EnergyTracker(sim, forces={'cavity': cavity_force})
   temp_tracker = TemperatureTracker(sim, time_tracker)
   
   writer.add_energy_tracker(energy_tracker)
   writer.add_temperature_tracker(temp_tracker)
   
   # Single updater instead of multiple
   sim.operations.updaters.append(hoomd.update.CustomUpdater(
       action=writer,
       trigger=hoomd.trigger.Periodic(1)
   ))

Analysis Tools Integration
==========================

**Pandas integration:**

.. code-block:: python

   import pandas as pd
   from cavitymd.data import ObservableReader
   
   reader = ObservableReader('observables.h5')
   
   # Create DataFrame
   df = pd.DataFrame({
       'time': reader.get_time(),
       'kinetic_energy': reader.get_observable('energies/kinetic'),
       'potential_energy': reader.get_observable('energies/potential'),
       'temperature': reader.get_observable('temperatures/kinetic')
   })
   
   # Analysis with pandas
   print(df.describe())
   df.plot(x='time', y='temperature')
   
   reader.close()

**NumPy integration:**

.. code-block:: python

   import numpy as np
   from cavitymd.data import ObservableReader
   
   reader = ObservableReader('observables.h5')
   
   # Load as NumPy arrays
   data = {
       'time': reader.get_time(),
       'energy': reader.get_observable('energies/total'),
       'temperature': reader.get_observable('temperatures/kinetic')
   }
   
   # NumPy operations
   mean_energy = np.mean(data['energy'])
   energy_fft = np.fft.fft(data['energy'])
   
   reader.close()

Troubleshooting
===============

Common Issues
-------------

**1. File locked / already open:**

.. code-block:: python

   # Ensure previous file is closed
   if 'writer' in locals():
       writer.close()
   
   # Create new writer
   writer = ObservableWriter(...)

**2. SWMR mode not working:**

.. code-block:: python

   # Ensure HDF5 library supports SWMR
   import h5py
   print(f"HDF5 version: {h5py.version.hdf5_version}")
   # Need HDF5 >= 1.10.0
   
   # Enable explicitly
   writer = ObservableWriter(..., enable_swmr=True)

**3. Large file sizes:**

.. code-block:: python

   # Increase compression
   writer = ObservableWriter(
       ...,
       compression='gzip',
       compression_level=9  # Maximum compression
   )
   
   # Reduce output frequency
   writer = ObservableWriter(
       ...,
       output_period_ps=1.0  # Less frequent output
   )

**4. Slow I/O performance:**

.. code-block:: python

   # Optimize chunk size
   writer = ObservableWriter(
       ...,
       chunk_size=10000,  # Larger chunks
       flush_interval=50   # Less frequent flushing
   )

Next Steps
==========

- :doc:`analysis_tools` for analysis tracker documentation
- :doc:`../part3_advanced/performance` for optimization strategies
- :doc:`running_simulations` for basic simulation setup
- :doc:`../api/index` for API reference

