========
Analysis
========

Post-Simulation Analysis
========================

Cavity HOOMD provides built-in analysis tools for understanding cavity-coupled dynamics.

Energy Analysis
===============

Monitor energy conservation and component breakdown:

.. code-block:: python

   from hoomd.cavitymd import EnergyTracker
   
   # During simulation
   energy_tracker = EnergyTracker(
       simulation=simulation,
       components=['harmonic', 'lj', 'cavity'],
       track_reservoirs=True
   )

   # Post-analysis
   data = np.loadtxt("energy_log.txt")
   total_energy = data[:, -1]  # Last column is total energy
   print(f"Energy drift: {np.std(total_energy):.6f}")

Correlation Functions
=====================

F(k,t) density correlation analysis:

.. code-block:: python

   from hoomd.cavitymd import FieldAutocorrelationTracker
   
   fkt_tracker = FieldAutocorrelationTracker(
       simulation=simulation,
       observable="density_correlation",
       kmag=1.0
   )

Trajectory Analysis
===================

Analyze saved GSD trajectories:

.. code-block:: python

   import gsd.hoomd
   
   with gsd.hoomd.open("trajectory.gsd", "r") as traj:
       for frame in traj:
           positions = frame.particles.position
           # Analyze positions, compute observables
           
.. note::
   
   Detailed analysis examples and post-processing scripts are available 
   in the examples directory. 