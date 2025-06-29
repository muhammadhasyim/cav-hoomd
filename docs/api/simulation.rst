==========
Simulation
==========

This module provides simulation management utilities for cavity molecular dynamics.

.. currentmodule:: hoomd.cavitymd

Adaptive Timestep Control
==========================

AdaptiveTimestepUpdater
-----------------------

.. autoclass:: AdaptiveTimestepUpdater
   :members:
   :show-inheritance:

   The `AdaptiveTimestepUpdater` automatically adjusts the simulation timestep based on 
   force magnitudes and error tolerances to maintain numerical stability while maximizing 
   computational efficiency.

   .. rubric:: Key Features

   - Automatic timestep adaptation based on force magnitudes
   - Configurable error tolerance and time constants
   - Integration with cavity thermostats and molecular dynamics
   - Real-time monitoring of simulation stability

Usage Examples
==============

Basic Adaptive Timestep Setup
------------------------------

.. code-block:: python

   from hoomd.cavitymd import AdaptiveTimestepUpdater

   # Create adaptive timestep updater
   adaptive_action = AdaptiveTimestepUpdater(
       state=simulation.state,
       integrator=simulation.operations.integrator,
       error_tolerance=0.01,           # Target error tolerance
       time_constant_ps=50.0,          # Adaptation time constant
       adaptiveerror=True,             # Enable adaptive error
       cavity_damping_factor=1.0,      # Cavity damping
       molecular_thermostat_tau=5.0,   # Molecular tau (ps)
       cavity_thermostat_tau=5.0       # Cavity tau (ps)
   )

   # Add to simulation as updater
   import hoomd
   adaptive_updater = hoomd.update.CustomUpdater(
       action=adaptive_action,
       trigger=hoomd.trigger.Periodic(100)  # Update every 100 steps
   )
   simulation.operations.updaters.append(adaptive_updater)

Advanced Configuration
----------------------

.. code-block:: python

   # Advanced adaptive timestep with time tracking
   from hoomd.cavitymd import ElapsedTimeTracker
   
   # Create time tracker for runtime control
   time_tracker = ElapsedTimeTracker(simulation, runtime_ps=1000.0)
   
   # Create adaptive timestep updater with time tracking
   adaptive_action = AdaptiveTimestepUpdater(
       state=simulation.state,
       integrator=simulation.operations.integrator,
       error_tolerance=0.005,          # Tighter tolerance
       time_constant_ps=25.0,          # Faster adaptation
       initial_fraction=1e-4,          # Conservative start
       adaptiveerror=True,
       cavity_damping_factor=2.0,      # Enhanced damping
       molecular_thermostat_tau=10.0,
       cavity_thermostat_tau=2.0,
       time_tracker=time_tracker       # Optional time tracking
   )

Integration with Energy Tracking
---------------------------------

.. code-block:: python

   # Complete simulation setup with adaptive timestep
   import hoomd
   from hoomd.cavitymd import (
       AdaptiveTimestepUpdater, ElapsedTimeTracker, 
       EnergyTracker, CavityModeTracker
   )

   # Setup simulation (assume simulation object exists)
   # ... simulation setup code ...

   # Create time tracker
   time_tracker = ElapsedTimeTracker(simulation, runtime_ps=500.0)

   # Create adaptive timestep updater
   adaptive_action = AdaptiveTimestepUpdater(
       state=simulation.state,
       integrator=simulation.operations.integrator,
       error_tolerance=0.01,
       time_constant_ps=50.0,
       time_tracker=time_tracker
   )

   # Create energy tracker for monitoring
   energy_tracker = EnergyTracker(
       simulation=simulation,
       components=['harmonic', 'lj', 'cavity'],
       output_prefix="energy_log",
       output_period_steps=100
   )

   # Add all updaters to simulation
   simulation.operations.updaters.extend([
       hoomd.update.CustomUpdater(
           action=time_tracker, 
           trigger=hoomd.trigger.Periodic(1)
       ),
       hoomd.update.CustomUpdater(
           action=adaptive_action, 
           trigger=hoomd.trigger.Periodic(100)
       ),
       hoomd.update.CustomUpdater(
           action=energy_tracker, 
           trigger=hoomd.trigger.Periodic(100)
       )
   ])

   # Run simulation - will automatically terminate when runtime_ps is reached
   # and maintain optimal timestep throughout
   simulation.run(999999999)  # Large number, limited by time_tracker

See Also
========

* :class:`~hoomd.cavitymd.ElapsedTimeTracker` - For runtime control
* :class:`~hoomd.cavitymd.EnergyTracker` - For energy monitoring
* :class:`~hoomd.cavitymd.CavityModeTracker` - For cavity analysis

.. note::
   
   **High-Level Simulation Class**: For a complete simulation framework that 
   includes `AdaptiveTimestepUpdater` and other features automatically, see 
   the `CavityMDSimulation` class in ``examples/05_advanced_run.py``. This 
   example class demonstrates how to integrate all cavity MD components into 
   a comprehensive simulation workflow. 