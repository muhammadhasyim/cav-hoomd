===========
Quick Start
===========

This guide will get you running cavity-coupled molecular dynamics simulations using the **05_advanced_run.py** script in just a few minutes.

Installation
============

.. code-block:: bash

   # Install from source
   git clone https://github.com/yourusername/cavity-hoomd.git
   cd cavity-hoomd
   pip install -e .

**Verify Installation**

.. code-block:: bash

   python 05_advanced_run.py --help

If the help message appears, you're ready to go!

Your First Simulation
=====================

The easiest way to run a cavity simulation is using the command line:

.. code-block:: bash

   # Run a basic cavity simulation
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000

This command will:
- Use Bussi thermostat for molecules, Langevin for cavity
- Set cavity coupling strength to 1e-3
- Run for 1000 ps
- Use default temperature (100 K) and frequency (2000 cm⁻¹)

**Expected Output**

The simulation creates an output directory and generates several files:

.. code-block:: text

   bussi_langevin_finiteq_coupling_1e-03/
   ├── prod-1.gsd              # Trajectory file
   ├── prod-1-energy.txt       # Energy tracking data
   ├── prod-1-cavity_mode.txt  # Cavity mode properties  
   ├── prod-1.log              # Simulation log

Available Experiments
=====================

The script provides four pre-configured experiment types:

.. code-block:: bash

   # Bussi (molecules) + Langevin (cavity) with finite-q mode
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000

   # Same but with q=0 mode
   python 05_advanced_run.py --experiment bussi_langevin_no_finiteq --coupling 1e-3 --runtime 1000

   # Langevin for both molecules and cavity
   python 05_advanced_run.py --experiment langevin_langevin --coupling 1e-3 --runtime 1000

   # Bussi for both molecules and cavity
   python 05_advanced_run.py --experiment bussi_bussi --coupling 1e-3 --runtime 1000

Control Simulations (No Cavity)
================================

Run without cavity coupling for comparison:

.. code-block:: bash

   # Molecular-only simulation (no cavity)
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --no-cavity --runtime 1000

Parameter Sweeps
================

Explore multiple parameter values in one command:

.. code-block:: bash

   # Sweep over coupling strengths
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3,1e-4,1e-5 --runtime 500

   # Sweep over temperatures
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --temperature 100,200,300 --runtime 500

   # Sweep over frequencies
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --frequency 1800,2000,2200 --runtime 500

Multiple Replicas
=================

Run multiple independent replicas:

.. code-block:: bash

   # Run replicas 1-5
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --replicas "1-5" --runtime 1000

   # Run specific replicas
   python 05_advanced_run.py --experiment bussi_langevin_finiteq \
       --coupling 1e-3 --replicas "1,3,5,7" --runtime 1000

Common Options
==============

**Basic Parameters**

.. code-block:: bash

   --experiment          # Experiment type (required)
   --coupling           # Coupling strength(s)
   --temperature        # Temperature in K (default: 100)
   --frequency          # Cavity frequency in cm⁻¹ (default: 2000)
   --runtime            # Simulation time in ps (default: 500)
   --no-cavity          # Disable cavity (control simulation)

**Advanced Features**

.. code-block:: bash

   --enable-energy-tracker    # Detailed energy tracking
   --enable-fkt              # F(k,t) correlation analysis
   --fixed-timestep          # Use fixed timestep instead of adaptive
   --device GPU              # Use GPU acceleration

**Output Control**

.. code-block:: bash

   --energy-output-period-ps 0.1     # Energy output frequency
   --gsd-output-period-ps 50.0       # Trajectory output frequency
   --console-output-period-ps 1.0    # Console update frequency

Quick Analysis
==============

**Check Energy Conservation**

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt

   # Read energy data
   data = pd.read_csv('bussi_langevin_finiteq_coupling_1e-03/prod-1-energy.txt', delimiter='\t')

   # Plot total energy
   plt.figure(figsize=(10, 6))
   plt.plot(data['time_ps'], data['total_energy'])
   plt.xlabel('Time (ps)')
   plt.ylabel('Total Energy (Hartree)')
   plt.title('Energy Conservation')
   plt.show()

   # Calculate energy drift
   drift = (data['total_energy'].iloc[-1] - data['total_energy'].iloc[0]) / data['total_energy'].iloc[0]
   print(f"Relative energy drift: {drift:.2e}")

**Analyze Cavity Mode**

.. code-block:: python

   # Read cavity mode data
   cavity_data = pd.read_csv('bussi_langevin_finiteq_coupling_1e-03/prod-1-cavity_mode.txt', delimiter='\t')

   # Plot cavity amplitude
   plt.figure(figsize=(10, 6))
   plt.plot(cavity_data['time_ps'], cavity_data['amplitude'])
   plt.xlabel('Time (ps)')
   plt.ylabel('Cavity Amplitude')
   plt.title('Cavity Mode Dynamics')
   plt.show()

HPC Usage (SLURM)
=================

The script automatically detects SLURM environments:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --array=1-10
   #SBATCH --time=24:00:00

   # This will automatically use SLURM_ARRAY_TASK_ID as replica number
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 5000

Common Issues & Solutions
=========================

**"No GSD file found"**

The script looks for ``init-0.gsd`` in the parent directory. Make sure you have an initial structure file:

.. code-block:: bash

   # Check if file exists
   ls ../init-0.gsd

**"GPU not available"**

Force CPU usage:

.. code-block:: bash

   python 05_advanced_run.py --experiment bussi_langevin_finiteq --device CPU --coupling 1e-3 --runtime 1000

**Simulation runs too slowly**

Use larger timesteps or reduce output frequency:

.. code-block:: bash

   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --fixed-timestep --timestep 2.0 --gsd-output-period-ps 100.0

Next Steps
==========

🎉 **Congratulations!** You've run your first cavity-coupled simulation.

**Learn More:**

1. **Explore Parameters**: Try different coupling strengths, temperatures, and frequencies
2. **Compare Experiments**: Run different thermostat combinations  
3. **Analyze Results**: Use the generated data files for deeper analysis
4. **Scale Up**: Run parameter sweeps and multiple replicas

**Resources:**

- :doc:`user_guide/index` - Complete usage guide
- :doc:`theory` - Scientific background
- `GitHub <https://github.com/yourusername/cavity-hoomd>`_ - Source code and issues

**Need Help?**

- Check the :doc:`user_guide/troubleshooting` section
- Post an issue on `GitHub <https://github.com/yourusername/cavity-hoomd/issues>`_
- Read the complete help: ``python 05_advanced_run.py --help`` 