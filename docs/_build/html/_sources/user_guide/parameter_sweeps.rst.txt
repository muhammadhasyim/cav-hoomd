================
Parameter Sweeps
================

Overview
========

Parameter sweeps allow systematic exploration of cavity-coupling effects across different system parameters. This guide shows how to efficiently run multiple simulations with varying parameters.

Basic Parameter Sweep
======================

Here's a simple example sweeping over coupling strengths:

.. code-block:: python

   import numpy as np
   from hoomd.cavitymd import CavityForce

   # Define parameter range
   coupling_values = np.logspace(-4, -2, 10)  # 10 values from 1e-4 to 1e-2

   for i, coupling in enumerate(coupling_values):
       # Setup simulation with current parameters
       cavity_force = CavityForce(
           kvector=[0, 0, 1],
           couplstr=coupling,
           omegac=0.1
       )
       
       # Run simulation with unique output
       run_simulation(coupling, output_dir=f"coupling_{coupling:.1e}")

Multi-Parameter Sweeps
======================

Sweep over multiple parameters simultaneously:

.. code-block:: python

   import itertools

   coupling_values = [1e-4, 1e-3, 1e-2]
   temperature_values = [200, 300, 400]  # Kelvin
   frequency_values = [1500, 2000, 2500]  # cm⁻¹

   # Generate all combinations
   for coupling, temp, freq in itertools.product(
       coupling_values, temperature_values, frequency_values
   ):
       run_simulation(coupling=coupling, temperature=temp, frequency=freq)

HPC Integration
===============

For high-performance computing systems with job arrays:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --array=1-100
   #SBATCH --time=24:00:00

   # Use SLURM_ARRAY_TASK_ID to select parameters
   python run_parameter_sweep.py --task-id $SLURM_ARRAY_TASK_ID

.. note::
   
   See ``examples/05_advanced_run.py`` for a complete parameter sweep implementation 
   with SLURM integration and automatic parameter selection. 