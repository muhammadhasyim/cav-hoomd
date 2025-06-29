==========
Simulation
==========

This module provides high-level simulation management classes and utilities.

.. currentmodule:: hoomd.cavitymd

High-Level Simulation Classes
=============================

CavityMDSimulation
------------------

.. autoclass:: CavityMDSimulation
   :members:
   :show-inheritance:

   .. rubric:: Main Methods

   .. autosummary::
      :nosignatures:

      ~CavityMDSimulation.__init__
      ~CavityMDSimulation.run
      ~CavityMDSimulation.setup_logging
      ~CavityMDSimulation.cleanup

   .. rubric:: Setup Methods

   .. autosummary::
      :nosignatures:

      ~CavityMDSimulation.setup_simulation
      ~CavityMDSimulation.setup_integrator
      ~CavityMDSimulation.setup_force_parameters
      ~CavityMDSimulation.setup_thermostat_parameters
      ~CavityMDSimulation.setup_trackers_and_loggers
      ~CavityMDSimulation.setup_output_writers

   .. rubric:: Utility Methods

   .. autosummary::
      :nosignatures:

      ~CavityMDSimulation.create_cavity_particle
      ~CavityMDSimulation.thermalize_system
      ~CavityMDSimulation.compute_and_set_optimal_timestep

Adaptive Timestep Control
==========================

AdaptiveTimestepUpdater
-----------------------

.. autoclass:: AdaptiveTimestepUpdater
   :members:
   :show-inheritance:

Usage Examples
==============

Basic Simulation Setup
-----------------------

.. code-block:: python

   from hoomd.cavitymd import CavityMDSimulation

   # Create a basic cavity MD simulation
   sim = CavityMDSimulation(
       job_dir="./simulation",
       replica=1,
       freq=2000.0,                    # Cavity frequency (cm⁻¹)
       couplstr=1e-3,                  # Coupling strength
       incavity=True,                  # Enable cavity coupling
       runtime_ps=1000.0,              # Runtime in picoseconds
       temperature=300.0,              # Temperature in Kelvin
       molecular_thermostat='bussi',   # Molecular thermostat
       cavity_thermostat='langevin'    # Cavity thermostat
   )

   # Run the simulation
   exit_code = sim.run()

Advanced Configuration
----------------------

.. code-block:: python

   # Advanced simulation with custom parameters
   sim = CavityMDSimulation(
       job_dir="./advanced_sim",
       replica=1,
       freq=1800.0,
       couplstr=5e-4,
       incavity=True,
       runtime_ps=2000.0,
       temperature=250.0,
       
       # Thermostat configuration
       molecular_thermostat='bussi',
       cavity_thermostat='langevin',
       molecular_thermostat_tau=10.0,  # ps
       cavity_thermostat_tau=2.0,      # ps
       cavity_damping_factor=2.0,
       
       # Advanced features
       finite_q=True,                  # Allow finite-q modes
       error_tolerance=0.005,          # Adaptive timestep
       enable_energy_tracking=True,    # Energy conservation
       enable_fkt=True,                # F(k,t) calculations
       
       # Output control
       energy_output_period_ps=0.05,
       gsd_output_period_ps=25.0,
       
       # Device configuration
       device='GPU',
       gpu_id=0
   )

   # Run with comprehensive logging
   exit_code = sim.run()

Parameter Sweeps
----------------

.. code-block:: python

   import numpy as np

   # Parameter sweep over coupling strengths
   coupling_values = [1e-4, 5e-4, 1e-3, 5e-3]
   temperatures = [200, 250, 300]

   for i, (coupling, temp) in enumerate(
       itertools.product(coupling_values, temperatures)
   ):
       sim = CavityMDSimulation(
           job_dir=f"./sweep_{i}",
           replica=1,
           freq=2000.0,
           couplstr=coupling,
           incavity=True,
           runtime_ps=500.0,
           temperature=temp,
           molecular_thermostat='bussi',
           cavity_thermostat='langevin'
       )
       
       print(f"Running coupling={coupling}, T={temp}K")
       sim.run()

No-Cavity Control Simulations
------------------------------

.. code-block:: python

   # Run control simulation without cavity
   control_sim = CavityMDSimulation(
       job_dir="./control",
       replica=1,
       freq=2000.0,      # Ignored when incavity=False
       couplstr=0.0,     # Ignored when incavity=False
       incavity=False,   # Disable cavity coupling
       runtime_ps=1000.0,
       temperature=300.0,
       molecular_thermostat='bussi'
       # cavity_thermostat not needed for no-cavity runs
   )

   control_sim.run()

Logging and Output Control
==========================

Logging Configuration
---------------------

.. code-block:: python

   sim = CavityMDSimulation(
       # ... other parameters ...
       
       # Logging options
       log_to_file=True,
       log_to_console=True,
       log_level='INFO',
       custom_log_file='my_simulation.log'
   )

Output Period Control
---------------------

.. code-block:: python

   sim = CavityMDSimulation(
       # ... other parameters ...
       
       # Fine-grained output control
       energy_output_period_ps=0.1,     # High-frequency energy data
       fkt_output_period_ps=1.0,        # F(k,t) calculations
       gsd_output_period_ps=50.0,       # Trajectory snapshots
       console_output_period_ps=5.0     # Console updates
   )

Error Handling
==============

The simulation framework provides robust error handling:

.. code-block:: python

   try:
       sim = CavityMDSimulation(...)
       exit_code = sim.run()
       
       if exit_code == 0:
           print("Simulation completed successfully")
       else:
           print("Simulation failed")
           
   except Exception as e:
       print(f"Setup error: {e}")

The simulation automatically handles:

- GPU initialization failures (falls back to CPU)
- Invalid parameter combinations
- File I/O errors
- Memory allocation issues
- HOOMD internal errors

Performance Considerations
==========================

Device Selection
----------------

.. code-block:: python

   # GPU acceleration (recommended for large systems)
   sim = CavityMDSimulation(
       device='GPU',
       gpu_id=0,  # Use first GPU
       # ... other parameters ...
   )

   # CPU for smaller systems or debugging
   sim = CavityMDSimulation(
       device='CPU',
       # ... other parameters ...
   )

Timestep Control
----------------

.. code-block:: python

   # Adaptive timestep (recommended)
   sim = CavityMDSimulation(
       error_tolerance=0.01,  # Adaptive control
       # ... other parameters ...
   )

   # Fixed timestep
   sim = CavityMDSimulation(
       error_tolerance=0.0,   # Disable adaptive
       dt_fs=1.0,            # Fixed 1 fs timestep
       # ... other parameters ...
   )

Memory Optimization
-------------------

.. code-block:: python

   # For large systems, limit output frequency
   sim = CavityMDSimulation(
       # ... other parameters ...
       
       # Reduce output frequency to save memory/disk space
       gsd_output_period_ps=100.0,      # Less frequent snapshots
       energy_output_period_ps=1.0,     # Less frequent energy data
       max_energy_output_time_ps=500.0  # Limit energy output duration
   ) 