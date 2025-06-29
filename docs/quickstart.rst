===========
Quick Start
===========

Get up and running with Cavity HOOMD in just a few minutes! This guide will walk you through installation and your first cavity-coupled simulation.

Installation
============

**Option 1: Install from PyPI (Recommended)**

.. code-block:: bash

   pip install cavity-hoomd

**Option 2: Install from Source**

.. code-block:: bash

   git clone https://github.com/yourusername/cavity-hoomd.git
   cd cavity-hoomd
   pip install -e .

**Verify Installation**

.. code-block:: python

   import hoomd
   from hoomd.cavitymd import CavityForce, CavityMDSimulation
   print("Cavity HOOMD installed successfully!")

Your First Simulation
=====================

Let's run a simple cavity-coupled molecular dynamics simulation:

.. code-block:: python

   from hoomd.cavitymd import CavityMDSimulation
   
   # Create a cavity MD simulation
   sim = CavityMDSimulation(
       job_dir="./my_first_sim",
       replica=1,
       freq=2000.0,        # Cavity frequency in cm⁻¹
       couplstr=1e-3,      # Coupling strength
       incavity=True,      # Enable cavity coupling
       runtime_ps=100.0,   # Short runtime for quick test
       temperature=300.0,  # Temperature in K
       molecular_thermostat='bussi',
       cavity_thermostat='langevin',
       input_gsd='init.gsd'  # You'll need an initial structure
   )
   
   # Run the simulation
   exit_code = sim.run()
   
   if exit_code == 0:
       print("Simulation completed successfully!")
   else:
       print("Simulation failed - check the logs")

.. note::

   You'll need an initial GSD file containing your molecular system. See :doc:`user_guide/basic_usage` for details on preparing initial structures.

Understanding the Output
========================

Your simulation will generate several output files:

.. code-block:: text

   my_first_sim/
   ├── prod-1.gsd              # Trajectory file
   ├── prod-1-energy.txt       # Energy tracking data
   ├── prod-1-cavity_mode.txt  # Cavity mode properties
   ├── prod-1-fkt.txt          # F(k,t) correlation data
   └── prod-1.log              # Simulation log

**Quick Analysis**

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt
   
   # Read energy data
   energy_data = pd.read_csv('my_first_sim/prod-1-energy.txt', delimiter='\t')
   
   # Plot total energy conservation
   plt.figure(figsize=(10, 6))
   plt.plot(energy_data['time_ps'], energy_data['total_energy'])
   plt.xlabel('Time (ps)')
   plt.ylabel('Total Energy (Hartree)')
   plt.title('Energy Conservation')
   plt.show()
   
   # Check energy drift
   energy_drift = (energy_data['total_energy'].iloc[-1] - 
                   energy_data['total_energy'].iloc[0]) / energy_data['total_energy'].iloc[0]
   print(f"Relative energy drift: {energy_drift:.2e}")

Key Parameters Explained
========================

**Cavity Parameters**

- **freq**: Cavity frequency in cm⁻¹ (typical range: 1000-4000)
- **couplstr**: Coupling strength (typical range: 1e-5 to 1e-2)
- **incavity**: Whether to enable cavity coupling (True/False)

**Simulation Parameters**

- **runtime_ps**: Simulation time in picoseconds
- **temperature**: Temperature in Kelvin
- **molecular_thermostat**: Thermostat for molecules ('bussi', 'langevin', 'none')
- **cavity_thermostat**: Thermostat for cavity ('bussi', 'langevin', 'none')

**Common Combinations**

.. tabs::

   .. tab:: Weak Coupling Study

      .. code-block:: python

         sim = CavityMDSimulation(
             # ... other parameters ...
             freq=2000.0,
             couplstr=1e-4,          # Weak coupling
             molecular_thermostat='bussi',
             cavity_thermostat='langevin'
         )

   .. tab:: Strong Coupling Study

      .. code-block:: python

         sim = CavityMDSimulation(
             # ... other parameters ...
             freq=1800.0,
             couplstr=1e-2,          # Strong coupling
             molecular_thermostat='bussi',
             cavity_thermostat='bussi'
         )

   .. tab:: Control (No Cavity)

      .. code-block:: python

         sim = CavityMDSimulation(
             # ... other parameters ...
             incavity=False,         # Disable cavity
             molecular_thermostat='bussi'
         )

Common First-Time Issues
========================

**Issue: "No GSD file found"**

.. code-block:: python

   # Solution: Create or specify correct path to initial structure
   sim = CavityMDSimulation(
       # ... other parameters ...
       input_gsd='path/to/your/initial_structure.gsd'
   )

**Issue: "GPU not available"**

.. code-block:: python

   # Solution: Use CPU explicitly
   sim = CavityMDSimulation(
       # ... other parameters ...
       device='CPU'
   )

**Issue: "Energy not conserved"**

.. code-block:: python

   # Solution: Enable adaptive timestep
   sim = CavityMDSimulation(
       # ... other parameters ...
       error_tolerance=0.01,  # Enable adaptive timestep
       enable_energy_tracking=True
   )

Next Steps
==========

🎉 **Congratulations!** You've run your first cavity-coupled simulation.

**What to do next:**

1. **Explore the Results**: Use the analysis examples above to understand your simulation output

2. **Read the User Guide**: Check out :doc:`user_guide/basic_usage` for detailed explanations

3. **Try Parameter Sweeps**: Learn how to explore parameter spaces systematically in :doc:`user_guide/parameter_sweeps`

4. **Understand the Theory**: Read :doc:`theory` for the scientific background

5. **Optimize Performance**: See :doc:`user_guide/performance` for getting the best performance

**Useful Resources:**

- :doc:`examples/index` - Complete working examples
- :doc:`api/index` - Detailed API reference  
- :doc:`user_guide/analysis` - Data analysis techniques
- `GitHub Repository <https://github.com/yourusername/cavity-hoomd>`_ - Source code and issues

**Having Problems?**

- Check the :doc:`user_guide/index` for detailed troubleshooting
- Look at the :doc:`examples/index` for working code
- Post an issue on `GitHub <https://github.com/yourusername/cavity-hoomd/issues>`_

Quick Reference
===============

**Most Common Commands**

.. code-block:: python

   # Basic simulation
   from hoomd.cavitymd import CavityMDSimulation
   sim = CavityMDSimulation(job_dir="./sim", replica=1, freq=2000, couplstr=1e-3, incavity=True, runtime_ps=1000, temperature=300)
   sim.run()

   # Parameter sweep
   for coupling in [1e-4, 1e-3, 1e-2]:
       sim = CavityMDSimulation(job_dir=f"./sweep_{coupling}", replica=1, couplstr=coupling, ...)
       sim.run()

   # Analysis
   import pandas as pd
   data = pd.read_csv('sim/prod-1-energy.txt', delimiter='\t')
   print(f"Energy drift: {(data['total_energy'].iloc[-1] - data['total_energy'].iloc[0]) / data['total_energy'].iloc[0]:.2e}")

**Key File Locations**

- **Trajectory**: `{job_dir}/prod-{replica}.gsd`
- **Energy data**: `{job_dir}/prod-{replica}-energy.txt`
- **Cavity mode**: `{job_dir}/prod-{replica}-cavity_mode.txt`
- **Logs**: `{job_dir}/prod-{replica}.log` 