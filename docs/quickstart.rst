===========
Quick Start
===========

Get up and running with cavity-coupled molecular dynamics in just a few minutes using the **05_advanced_run.py** script.

Installation
============

.. code-block:: bash

   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd
   ./build_install.sh

**Verify Installation**

.. code-block:: bash

   python examples/05_advanced_run.py --help

If the help message appears, you're ready to go!

Your First Simulation
=====================

Run a basic cavity simulation:

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000

This will:
- Set cavity coupling strength to 1e-3
- Run for 1000 ps
- Use default temperature (100 K) and frequency (2000 cm⁻¹)
- Create output files with trajectory and energy data

**Expected Output**

The simulation creates files like:

.. code-block:: text

   cavity_coupling_1e-03/
   ├── prod-1.gsd              # Trajectory file
   ├── prod-1-energy.txt       # Energy tracking data
   ├── prod-1-cavity_mode.txt  # Cavity mode properties  
   └── prod-1.log              # Simulation log

Control Simulation
==================

Run without cavity coupling for comparison:

.. code-block:: bash

   python examples/05_advanced_run.py --no-cavity --runtime 1000

Common Usage Examples
=====================

**Different Parameters**

.. code-block:: bash

   # High temperature
   python examples/05_advanced_run.py --coupling 1e-3 --temperature 300 --runtime 1000

   # Different frequency
   python examples/05_advanced_run.py --coupling 1e-3 --frequency 1800 --runtime 1000

   # Strong coupling
   python examples/05_advanced_run.py --coupling 1e-2 --runtime 1000

**Multiple Replicas**

.. code-block:: bash

   # Run replicas 1-5
   python examples/05_advanced_run.py --coupling 1e-3 --replicas "1-5" --runtime 1000

**Advanced Features**

.. code-block:: bash

   # Enable detailed tracking
   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000 \
       --enable-energy-tracker --enable-fkt

   # Use GPU acceleration
   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000 --device GPU

   # Different thermostat combinations
   python examples/05_advanced_run.py --molecular-bath bussi --cavity-bath langevin \
       --coupling 1e-3 --runtime 1000

Key Options
===========

**Basic Parameters:**
- ``--coupling`` - Coupling strength (e.g., 1e-3)
- ``--temperature`` - Temperature in K (default: 100)
- ``--frequency`` - Cavity frequency in cm⁻¹ (default: 2000)
- ``--runtime`` - Simulation time in ps (default: 500)
- ``--no-cavity`` - Run without cavity (control simulation)

**Thermostat Options:**
- ``--molecular-bath`` - Molecular thermostat: bussi, langevin, none (default: bussi)
- ``--cavity-bath`` - Cavity thermostat: bussi, langevin, none (default: langevin)
- ``--finite-q`` - Enable finite-q cavity mode

**Advanced:**
- ``--replicas`` - Run multiple replicas (e.g., "1-5")
- ``--device GPU`` - Use GPU acceleration
- ``--enable-energy-tracker`` - Detailed energy tracking
- ``--enable-fkt`` - F(k,t) correlation analysis

Jupyter Notebook
================

For interactive usage, see the Jupyter notebook:

.. code-block:: bash

   jupyter notebook examples/05_advanced_run.ipynb

This notebook shows the same simulation setup with interactive analysis.

Quick Analysis
==============

**Check Results with Python**

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt

   # Read energy data
   data = pd.read_csv('cavity_coupling_1e-03/prod-1-energy.txt', delimiter='\t')

   # Plot energy over time
   plt.plot(data['time_ps'], data['total_energy'])
   plt.xlabel('Time (ps)')
   plt.ylabel('Total Energy (Hartree)')
   plt.show()

**Check Energy Conservation**

.. code-block:: python

   # Calculate energy drift
   drift = (data['total_energy'].iloc[-1] - data['total_energy'].iloc[0])
   drift_percent = drift / data['total_energy'].iloc[0] * 100
   print(f"Energy drift: {drift_percent:.3f}%")

Next Steps
==========

- Try different thermostat combinations
- Run multiple replicas to explore system behavior
- Compare cavity vs no-cavity simulations
- Use the Jupyter notebook for interactive analysis
- Check ``python examples/05_advanced_run.py --help`` for all options 