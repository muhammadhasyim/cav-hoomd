===========
Quick Start
===========

Get up and running with cavity-coupled molecular dynamics in just a few minutes using the **05_advanced_run.py** script with support for time-varying coupling and advanced analysis features.

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

**Basic Cavity Simulation**

Run a basic cavity simulation:

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000

This will:
- Set cavity coupling strength to 1e-3
- Run for 1000 ps
- Use default temperature (100 K) and frequency (2000 cm⁻¹)
- Create output files with trajectory and energy data

**Time-Varying Coupling Simulation**

Run a simulation with coupling that switches on at a specific time:

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 1.0 --runtime 1000

This will:
- Start with zero coupling (molecules only)
- Switch coupling to 1e-3 at t = 1.0 ps
- For finite-q mode, automatically displaces cavity particle to new equilibrium
- Monitor energy conservation during the switching process

**Expected Output**

The simulation creates organized output directories:

.. code-block:: text

   cavity_coupling_1e-03/               # Constant coupling
   ├── prod-1.gsd                      # Trajectory file
   ├── prod-1-energy.txt               # Energy tracking data
   ├── prod-1-cavity_mode.txt          # Cavity mode properties  
   └── prod-1.log                      # Simulation log
   
   cavity_coupling_1e-03_switch_1.0ps/ # Time-varying coupling
   ├── prod-1.gsd                      # Trajectory file
   ├── prod-1-energy.txt               # Energy tracking data
   ├── prod-1-cavity_mode.txt          # Cavity mode properties
   └── prod-1.log                      # Simulation log

Control Simulation
==================

Run without cavity coupling for comparison:

.. code-block:: bash

   python examples/05_advanced_run.py --no-cavity --runtime 1000

Advanced Features
=================

**Time-Varying Coupling**

The plugin now supports sophisticated time-varying coupling experiments:

.. code-block:: bash

   # Coupling switches from 0 to 1e-3 at t = 2.0 ps
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 2.0 --runtime 1000

   # Finite-q mode with coupling switching
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 1.0 --finite-q --runtime 1000

The time-varying coupling feature:
- Uses HOOMD variants for precise time control
- Automatically handles cavity particle displacement for finite-q modes
- Provides energy conservation monitoring during transitions
- Works with adaptive timestep algorithms

**Energy Conservation Testing**

Monitor energy conservation with detailed tracking:

.. code-block:: bash

   # Enable comprehensive energy tracking
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 1.0 \
       --enable-energy-tracker --runtime 1000

   # Test energy conservation in time-varying systems
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 1.0 \
       --enable-energy-tracker --molecular-bath none --cavity-bath none --runtime 1000

**GPU Acceleration**

Take advantage of optimized GPU kernels:

.. code-block:: bash

   # Use GPU for high-performance simulations
   python examples/05_advanced_run.py --coupling 1e-3 --device GPU --runtime 1000

   # GPU with time-varying coupling
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 1.0 \
       --device GPU --runtime 1000

Common Usage Examples
=====================

**Different Parameters**

.. code-block:: bash

   # High temperature
   python examples/05_advanced_run.py --coupling 1e-3 --temperature 300 --runtime 1000

   # Different frequency
   python examples/05_advanced_run.py --coupling 1e-3 --frequency 1800 --runtime 1000

   # Strong coupling with switching
   python examples/05_advanced_run.py --coupling 1e-2 --switch-time 0.5 --runtime 1000

**Multiple Replicas**

.. code-block:: bash

   # Run replicas 1-5 with constant coupling
   python examples/05_advanced_run.py --coupling 1e-3 --replicas "1-5" --runtime 1000

   # Run replicas with time-varying coupling
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 1.0 \
       --replicas "1-5" --runtime 1000

**Advanced Analysis**

.. code-block:: bash

   # Enable all analysis features
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 1.0 \
       --enable-energy-tracker --enable-fkt --runtime 1000

   # Adaptive timestep with detailed monitoring
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 1.0 \
       --enable-energy-tracker --runtime 1000

   # Different thermostat combinations
   python examples/05_advanced_run.py --molecular-bath bussi --cavity-bath langevin \
       --coupling 1e-3 --switch-time 1.0 --runtime 1000

Key Options
===========

**Basic Parameters:**
- ``--coupling`` - Coupling strength (e.g., 1e-3)
- ``--temperature`` - Temperature in K (default: 100)
- ``--frequency`` - Cavity frequency in cm⁻¹ (default: 2000)
- ``--runtime`` - Simulation time in ps (default: 500)
- ``--no-cavity`` - Run without cavity (control simulation)

**Time-Varying Parameters:**
- ``--switch-time`` - Time to switch coupling on in ps (enables time-varying mode)
- ``--damping-ratio`` - Cavity damping ratio (default: 0.0)

**Thermostat Options:**
- ``--molecular-bath`` - Molecular thermostat: bussi, langevin, none (default: bussi)
- ``--cavity-bath`` - Cavity thermostat: bussi, langevin, none (default: langevin)
- ``--finite-q`` - Enable finite-q cavity mode

**Analysis Options:**
- ``--enable-energy-tracker`` - Detailed energy component tracking
- ``--enable-fkt`` - F(k,t) correlation analysis
- ``--fixed-timestep`` - Use fixed timestep instead of adaptive

**Performance Options:**
- ``--device GPU`` - Use GPU acceleration
- ``--replicas`` - Run multiple replicas (e.g., "1-5")

**Output Control:**
- ``--energy-output-period`` - Energy output frequency in ps (default: 0.1)
- ``--gsd-output-period`` - Trajectory output frequency in ps (default: 50.0)
- ``--console-output-period`` - Console output frequency in ps (default: 1.0)

Jupyter Notebook
================

For interactive usage and analysis, see the Jupyter notebook:

.. code-block:: bash

   jupyter notebook examples/05_advanced_run.ipynb

This notebook shows:
- Interactive parameter exploration
- Real-time analysis of time-varying coupling
- Energy conservation monitoring
- Cavity particle displacement visualization

Quick Analysis
==============

**Check Results with Python**

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt

   # Read energy data
   data = pd.read_csv('cavity_coupling_1e-03_switch_1.0ps/prod-1-energy.txt', delimiter='\t')

   # Plot energy over time
   plt.figure(figsize=(10, 6))
   plt.plot(data['time_ps'], data['total_energy'], label='Total Energy')
   plt.axvline(x=1.0, color='red', linestyle='--', label='Coupling Switch')
   plt.xlabel('Time (ps)')
   plt.ylabel('Total Energy (Hartree)')
   plt.legend()
   plt.show()

**Check Energy Conservation**

.. code-block:: python

   # Calculate energy drift
   drift = (data['total_energy'].iloc[-1] - data['total_energy'].iloc[0])
   drift_percent = drift / data['total_energy'].iloc[0] * 100
   print(f"Energy drift: {drift_percent:.6f}%")

   # Check energy conservation around switch time
   switch_data = data[(data['time_ps'] >= 0.5) & (data['time_ps'] <= 1.5)]
   pre_switch = switch_data[switch_data['time_ps'] < 1.0]['total_energy'].mean()
   post_switch = switch_data[switch_data['time_ps'] > 1.0]['total_energy'].mean()
   print(f"Energy change at switch: {post_switch - pre_switch:.6f} Hartree")

**Analyze Cavity Mode**

.. code-block:: python

   # Read cavity mode data
   cavity_data = pd.read_csv('cavity_coupling_1e-03_switch_1.0ps/prod-1-cavity_mode.txt', delimiter='\t')

   # Plot cavity position over time
   plt.figure(figsize=(10, 6))
   plt.plot(cavity_data['time_ps'], cavity_data['cavity_position_x'], label='X Position')
   plt.plot(cavity_data['time_ps'], cavity_data['cavity_position_y'], label='Y Position')
   plt.axvline(x=1.0, color='red', linestyle='--', label='Coupling Switch')
   plt.xlabel('Time (ps)')
   plt.ylabel('Cavity Position (a.u.)')
   plt.legend()
   plt.show()

Understanding Time-Varying Coupling
====================================

**Physical Interpretation**

Time-varying coupling simulates experimental scenarios where:
- Cavity coupling is switched on suddenly (pump-probe experiments)
- System equilibrates to new cavity-coupled state
- Energy redistribution between molecular and cavity modes occurs

**Key Observations**

1. **Energy Conservation**: Total energy should be conserved during switching
2. **Cavity Displacement**: For finite-q modes, cavity particle jumps to new equilibrium
3. **Relaxation Dynamics**: System equilibrates to new steady state
4. **Coupling Strength**: Determines magnitude of energy redistribution

**Best Practices**

- Use energy conservation tests to validate simulations
- Monitor cavity particle displacement for finite-q modes
- Compare with control simulations (no cavity coupling)
- Use appropriate thermostat combinations for your system

Next Steps
==========

- Explore different thermostat combinations
- Test time-varying coupling with different switch times
- Run multiple replicas for statistical analysis
- Compare cavity vs no-cavity simulations
- Use GPU acceleration for larger systems
- Check the comprehensive help: ``python examples/05_advanced_run.py --help`` 