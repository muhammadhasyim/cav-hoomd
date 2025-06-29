========
Analysis
========

This module provides comprehensive analysis and tracking capabilities for cavity MD simulations.

.. currentmodule:: hoomd.cavitymd

Energy Tracking
===============

EnergyTracker
-------------

.. autoclass:: EnergyTracker
   :members:
   :show-inheritance:

   Comprehensive energy tracking with reservoir monitoring and conservation analysis.

   .. rubric:: Key Features

   - Individual force component energies
   - Thermostat reservoir energies  
   - Kinetic energy tracking
   - Total energy conservation monitoring
   - Configurable output periods

Cavity Mode Analysis
====================

CavityModeTracker
-----------------

.. autoclass:: CavityModeTracker
   :members:
   :show-inheritance:

   Tracks cavity mode properties including position, velocity, and energy components.

   .. rubric:: Tracked Quantities

   - Cavity mode position and velocity
   - Kinetic and potential energies
   - Mode amplitude and phase
   - Coupling interaction strength

Correlation Functions
=====================

FieldAutocorrelationTracker
---------------------------

.. autoclass:: FieldAutocorrelationTracker
   :members:
   :show-inheritance:

   Computes field autocorrelation functions including F(k,t) density correlations.

   .. rubric:: Supported Observables

   - Density correlation functions F(k,t)
   - Custom field observables
   - Multiple wavevector sampling
   - Reference frame management

AutocorrelationTracker
----------------------

.. autoclass:: AutocorrelationTracker
   :members:
   :show-inheritance:

   Base class for autocorrelation function calculations.

DipoleAutocorrelation
---------------------

.. autoclass:: DipoleAutocorrelation
   :members:
   :show-inheritance:

   Specialized tracker for molecular dipole moment autocorrelations.

Utility Classes
===============

Status
------

.. autoclass:: Status
   :members:
   :show-inheritance:

   Performance monitoring and status reporting.

ElapsedTimeTracker
------------------

.. autoclass:: ElapsedTimeTracker
   :members:
   :show-inheritance:

   Tracks simulation time progress and handles runtime termination.

TimestepFormatter
-----------------

.. autoclass:: TimestepFormatter
   :members:
   :show-inheritance:

   Formats timestep information for display (atomic units to femtoseconds).

Usage Examples
==============

Energy Conservation Monitoring
-------------------------------

.. code-block:: python

   from hoomd.cavitymd import EnergyTracker

   # Set up comprehensive energy tracking
   energy_tracker = EnergyTracker(
       simulation=sim,
       components=['harmonic', 'lj', 'ewald_short', 'ewald_long', 'cavity'],
       force_objects={
           'cavity': cavity_force,
           'lj': lj_force,
           'harmonic': harmonic_force
       },
       thermostat_objects={
           'langevin_molecular': molecular_method,
           'langevin_cavity': cavity_method
       },
       time_tracker=time_tracker,
       output_prefix='energy_data',
       output_period_steps=100,
       track_reservoirs=True
   )

   # Add to simulation
   energy_updater = hoomd.update.CustomUpdater(
       action=energy_tracker,
       trigger=hoomd.trigger.Periodic(100)
   )
   sim.operations.updaters.append(energy_updater)

F(k,t) Density Correlations
----------------------------

.. code-block:: python

   from hoomd.cavitymd import FieldAutocorrelationTracker

   # Set up F(k,t) calculations  
   fkt_tracker = FieldAutocorrelationTracker(
       simulation=sim,
       observable="density_correlation",
       time_tracker=time_tracker,
       output_period_steps=1000,
       output_prefix='fkt_data',
       reference_interval_steps=5000,
       max_references=20,
       kmag=1.0,                    # Wavevector magnitude
       num_wavevectors=100          # Number of k-vectors to sample
   )

   # Add to simulation
   fkt_updater = hoomd.update.CustomUpdater(
       action=fkt_tracker,
       trigger=hoomd.trigger.Periodic(1)
   )
   sim.operations.updaters.append(fkt_updater)

Cavity Mode Monitoring
-----------------------

.. code-block:: python

   from hoomd.cavitymd import CavityModeTracker

   # Track cavity mode properties
   cavity_tracker = CavityModeTracker(
       simulation=sim,
       cavityforce=cavity_force,
       time_tracker=time_tracker,
       output_prefix='cavity_mode',
       output_period_steps=10
   )

   # Add to simulation
   cavity_updater = hoomd.update.CustomUpdater(
       action=cavity_tracker,
       trigger=hoomd.trigger.Periodic(10)
   )
   sim.operations.updaters.append(cavity_updater)

Custom Analysis Pipeline
------------------------

.. code-block:: python

   # Create comprehensive analysis pipeline
   def setup_analysis(sim, cavity_force, time_tracker):
       """Set up complete analysis pipeline."""
       
       # Energy tracking
       energy_tracker = EnergyTracker(
           simulation=sim,
           components=['cavity', 'molecular'],
           track_reservoirs=True,
           output_period_steps=100
       )
       
       # Cavity mode tracking
       cavity_tracker = CavityModeTracker(
           simulation=sim,
           cavityforce=cavity_force,
           time_tracker=time_tracker
       )
       
       # F(k,t) correlations
       fkt_tracker = FieldAutocorrelationTracker(
           simulation=sim,
           observable="density_correlation",
           kmag=1.0,
           num_wavevectors=50
       )
       
       # Add all trackers
       trackers = [energy_tracker, cavity_tracker, fkt_tracker]
       for i, tracker in enumerate(trackers):
           updater = hoomd.update.CustomUpdater(
               action=tracker,
               trigger=hoomd.trigger.Periodic(10 * (i + 1))
           )
           sim.operations.updaters.append(updater)
       
       return trackers

   # Use in simulation
   trackers = setup_analysis(sim, cavity_force, time_tracker)

Data Analysis
=============

Reading Output Files
--------------------

The analysis classes generate several types of output files:

.. code-block:: python

   import numpy as np
   import pandas as pd

   # Read energy tracking data
   energy_data = pd.read_csv('energy_data.txt', delimiter='\t')
   
   # Plot energy conservation
   import matplotlib.pyplot as plt
   
   plt.figure(figsize=(10, 6))
   plt.plot(energy_data['time_ps'], energy_data['total_energy'])
   plt.xlabel('Time (ps)')
   plt.ylabel('Total Energy (Hartree)')
   plt.title('Energy Conservation')
   plt.show()

   # Check energy conservation
   total_energy = energy_data['total_energy']
   energy_drift = (total_energy.iloc[-1] - total_energy.iloc[0]) / total_energy.iloc[0]
   print(f"Relative energy drift: {energy_drift:.2e}")

F(k,t) Analysis
---------------

.. code-block:: python

   # Read F(k,t) correlation data
   fkt_data = np.loadtxt('fkt_data.txt')
   
   # Extract time and correlation values
   times = fkt_data[:, 0]  # Time values
   correlations = fkt_data[:, 1:]  # F(k,t) for different k-vectors
   
   # Plot autocorrelation decay
   plt.figure(figsize=(10, 6))
   for i in range(min(5, correlations.shape[1])):
       plt.plot(times, correlations[:, i], label=f'k-vector {i+1}')
   
   plt.xlabel('Time (ps)')
   plt.ylabel('F(k,t)')
   plt.title('Density Correlation Functions')
   plt.legend()
   plt.yscale('log')
   plt.show()

   # Fit exponential decay
   def exponential_decay(t, A, tau):
       return A * np.exp(-t / tau)
   
   from scipy.optimize import curve_fit
   
   # Fit first correlation function
   popt, pcov = curve_fit(exponential_decay, times, correlations[:, 0])
   decay_time = popt[1]
   print(f"Decay time: {decay_time:.2f} ps")

Cavity Mode Analysis
--------------------

.. code-block:: python

   # Read cavity mode data
   cavity_data = pd.read_csv('cavity_mode.txt', delimiter='\t')
   
   # Plot cavity mode trajectory
   fig, axes = plt.subplots(2, 2, figsize=(12, 8))
   
   # Position components
   axes[0, 0].plot(cavity_data['time_ps'], cavity_data['position_x'])
   axes[0, 0].set_title('Cavity X Position')
   axes[0, 1].plot(cavity_data['time_ps'], cavity_data['position_y']) 
   axes[0, 1].set_title('Cavity Y Position')
   
   # Energy components
   axes[1, 0].plot(cavity_data['time_ps'], cavity_data['kinetic_energy'])
   axes[1, 0].set_title('Cavity Kinetic Energy')
   axes[1, 1].plot(cavity_data['time_ps'], cavity_data['potential_energy'])
   axes[1, 1].set_title('Cavity Potential Energy')
   
   plt.tight_layout()
   plt.show()

   # Calculate mode properties
   amplitude = np.sqrt(cavity_data['position_x']**2 + cavity_data['position_y']**2)
   phase = np.arctan2(cavity_data['position_y'], cavity_data['position_x'])
   
   print(f"Average amplitude: {amplitude.mean():.4f}")
   print(f"Amplitude std: {amplitude.std():.4f}")

Performance Monitoring
======================

The analysis framework includes built-in performance monitoring:

.. code-block:: python

   from hoomd.cavitymd import Status

   # Create status monitor
   status = Status(sim, runtime_ps=1000.0, time_tracker=time_tracker)

   # Access performance metrics
   print(f"Simulation progress: {status.progress:.1%}")
   print(f"Performance: {status.ns_per_day:.2f} ns/day")
   print(f"ETA: {status.eta}")

The Status class provides:

- Real-time progress monitoring
- Performance metrics (ns/day)
- Estimated time to completion
- Memory usage tracking
- GPU utilization (when available) 