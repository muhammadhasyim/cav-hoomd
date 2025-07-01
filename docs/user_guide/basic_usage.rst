===========
Basic Usage
===========

This guide covers the essential usage patterns for **05_advanced_run.py**, the main script for running cavity MD simulations.

Prerequisites
=============

Before starting, ensure you have:

* Cavity HOOMD installed
* An initial GSD file (usually ``init-0.gsd``)
* Basic familiarity with command-line interfaces

Command-Line Interface
======================

The script provides a comprehensive command-line interface:

.. code-block:: bash

   python 05_advanced_run.py [OPTIONS]

**View all options:**

.. code-block:: bash

   python 05_advanced_run.py --help

Essential Options
================

**Required Parameters**

.. code-block:: bash

   --experiment EXPERIMENT_TYPE    # Choose experiment type (required)

**Basic Parameters**

.. code-block:: bash

   --coupling COUPLING            # Coupling strength (e.g., 1e-3)
   --temperature TEMPERATURE      # Temperature in K (default: 100)
   --frequency FREQUENCY          # Cavity frequency in cm⁻¹ (default: 2000)
   --runtime RUNTIME             # Simulation time in ps (default: 500)

**Control Options**

.. code-block:: bash

   --no-cavity                   # Disable cavity (molecular-only simulation)
   --replicas "1-5"              # Run multiple replicas
   --device GPU                  # Use GPU acceleration

Experiment Types
================

Choose from four pre-configured experiment types:

**bussi_langevin_finiteq** (Recommended)

.. code-block:: bash

   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000

- Bussi thermostat for molecules (deterministic)
- Langevin thermostat for cavity (stochastic)
- Finite-q cavity mode (allows momentum)

**bussi_langevin_no_finiteq**

.. code-block:: bash

   python 05_advanced_run.py --experiment bussi_langevin_no_finiteq --coupling 1e-3 --runtime 1000

- Same as above but q=0 cavity mode
- Good for comparison studies

**langevin_langevin**

.. code-block:: bash

   python 05_advanced_run.py --experiment langevin_langevin --coupling 1e-3 --runtime 1000

- Langevin thermostat for both molecules and cavity
- Fully stochastic dynamics

**bussi_bussi**

.. code-block:: bash

   python 05_advanced_run.py --experiment bussi_bussi --coupling 1e-3 --runtime 1000

- Bussi thermostat for both molecules and cavity
- Deterministic temperature control

Basic Simulation Examples
=========================

**Your First Simulation**

.. code-block:: bash

   # Simple cavity simulation
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000

**Control Simulation**

.. code-block:: bash

   # Same parameters but without cavity
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --no-cavity --runtime 1000

**Different Coupling Strengths**

.. code-block:: bash

   # Weak coupling
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-4 --runtime 1000

   # Strong coupling  
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-2 --runtime 1000

**Temperature Studies**

.. code-block:: bash

   # Low temperature
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --temperature 50 --runtime 1000

   # High temperature
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --temperature 300 --runtime 1000

**Frequency Studies**

.. code-block:: bash

   # Low frequency
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --frequency 1500 --runtime 1000

   # High frequency
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --frequency 3000 --runtime 1000

Advanced Features
=================

**Energy Tracking**

.. code-block:: bash

   # Enable detailed energy tracking
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --enable-energy-tracker

**F(k,t) Correlation Analysis**

.. code-block:: bash

   # Enable F(k,t) correlation tracking
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --enable-fkt --fkt-kmag 1.0

**Fixed Timestep**

.. code-block:: bash

   # Use fixed timestep instead of adaptive
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --fixed-timestep --timestep 1.0

**GPU Acceleration**

.. code-block:: bash

   # Use GPU acceleration
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --device GPU --gpu-id 0

Output Control
==============

**Output Periods**

.. code-block:: bash

   # Control output frequencies (all in ps)
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --energy-output-period-ps 0.1 \
       --gsd-output-period-ps 50.0 \
       --console-output-period-ps 1.0

**Logging Options**

.. code-block:: bash

   # Enable file and console logging
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --log-to-file --log-to-console

Understanding Output Files
==========================

Each simulation creates a directory with output files:

.. code-block:: text

   bussi_langevin_finiteq_coupling_1e-03/
   ├── prod-1.gsd              # Trajectory file
   ├── prod-1-energy.txt       # Energy tracking data
   ├── prod-1-cavity_mode.txt  # Cavity mode properties
   ├── prod-1-fkt.txt          # F(k,t) correlation data (if enabled)
   └── prod-1.log              # Simulation log

**Key Output Files:**

- **prod-1.gsd**: Trajectory data for visualization and analysis
- **prod-1-energy.txt**: Energy conservation monitoring
- **prod-1-cavity_mode.txt**: Cavity amplitude and phase evolution
- **prod-1.log**: Detailed simulation log with performance metrics

Thermostat Parameters
====================

**Thermostat Time Constants**

.. code-block:: bash

   # Adjust thermostat coupling times (in ps)
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000 \
       --molecular-tau 5.0 --cavity-tau 1.0

**Typical Values:**

- **Molecular tau**: 5-10 ps (moderate coupling to heat bath)
- **Cavity tau**: 1-5 ps (faster coupling for cavity mode)

Best Practices
==============

**Coupling Strength Guidelines**

- **Weak coupling**: 1e-4 to 1e-3 (linear response regime)
- **Strong coupling**: 1e-3 to 1e-2 (nonlinear effects)
- **Very strong**: > 1e-2 (may require shorter timesteps)

**Runtime Guidelines**

- **Test runs**: 100-500 ps
- **Production runs**: 1000-10000 ps  
- **Long dynamics**: 10000+ ps

**Temperature Guidelines**

- **Low temperature**: 50-100 K (reduced thermal motion)
- **Room temperature**: 200-300 K (realistic conditions)
- **High temperature**: 300+ K (enhanced dynamics)

**Frequency Guidelines**

- **Infrared**: 1000-2000 cm⁻¹
- **Mid-infrared**: 2000-3000 cm⁻¹
- **Near-infrared**: 3000+ cm⁻¹

Common Patterns
===============

**Systematic Study**

.. code-block:: bash

   # Study different coupling strengths
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-4 --runtime 1000
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-2 --runtime 1000

**Comparing Thermostats**

.. code-block:: bash

   # Compare different thermostat combinations
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000
   python 05_advanced_run.py --experiment langevin_langevin --coupling 1e-3 --runtime 1000
   python 05_advanced_run.py --experiment bussi_bussi --coupling 1e-3 --runtime 1000

**Cavity vs Control**

.. code-block:: bash

   # Compare cavity vs no-cavity
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --no-cavity --runtime 1000

Next Steps
==========

* Learn about :doc:`parameter_sweeps` for automated parameter exploration
* Explore :doc:`analysis` for understanding output files
* Check :doc:`troubleshooting` for common issues and solutions 