==================
Advanced Features
==================

.. note::
   This guide covers advanced features of Cavity HOOMD. For basic usage, see :doc:`basic_usage`.

Adaptive Timestep Control
==========================

The :class:`~hoomd.cavitymd.AdaptiveTimestepUpdater` automatically adjusts timesteps based on force magnitudes and error tolerances:

.. code-block:: python

   from hoomd.cavitymd import AdaptiveTimestepUpdater

   adaptive_action = AdaptiveTimestepUpdater(
       state=simulation.state,
       integrator=simulation.operations.integrator,
       error_tolerance=0.01,
       time_constant_ps=50.0
   )

Energy Tracking
===============

Comprehensive energy monitoring with reservoir energy tracking:

.. code-block:: python

   from hoomd.cavitymd import EnergyTracker

   energy_tracker = EnergyTracker(
       simulation=simulation,
       components=['harmonic', 'lj', 'cavity'],
       output_prefix="energy_log",
       track_reservoirs=True
   )

Correlation Analysis
====================

F(k,t) density correlation functions:

.. code-block:: python

   from hoomd.cavitymd import FieldAutocorrelationTracker

   fkt_tracker = FieldAutocorrelationTracker(
       simulation=simulation,
       observable="density_correlation",
       kmag=1.0,
       num_wavevectors=50
   )

GPU Acceleration
================

Optimized GPU kernels for cavity forces and large system sizes:

.. code-block:: python

   # Enable GPU acceleration
   device = hoomd.device.GPU()
   simulation = hoomd.Simulation(device=device)

Custom Thermostats
==================

Integration with Bussi reservoir thermostats for enhanced temperature control.

.. note::
   **Coming Soon**: More detailed documentation for each advanced feature is being developed. 
   For now, see the API reference and examples for implementation details. 