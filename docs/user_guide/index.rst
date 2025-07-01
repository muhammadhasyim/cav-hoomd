==========
User Guide
==========

This guide covers everything you need to know about using **05_advanced_run.py** for cavity-coupled molecular dynamics simulations.

.. toctree::
   :maxdepth: 2

   basic_usage
   parameter_sweeps
   analysis
   troubleshooting

Overview
========

**05_advanced_run.py** provides a complete command-line interface for cavity MD simulations. It handles:

- All simulation setup and configuration
- Multiple thermostat combinations
- Parameter sweeps and replica management  
- Comprehensive analysis and output
- HPC and SLURM integration

Quick Navigation
================

.. grid:: 2
    :gutter: 3

    .. grid-item-card:: 🎯 **Basic Usage**
        :link: basic_usage
        :link-type: doc

        Learn the fundamental command-line options and run
        your first simulations.

    .. grid-item-card:: 🔄 **Parameter Sweeps**
        :link: parameter_sweeps
        :link-type: doc

        Explore parameter spaces efficiently with automated
        sweeps and multiple replicas.

    .. grid-item-card:: 📊 **Analysis**
        :link: analysis
        :link-type: doc

        Understand the output files and analyze your
        simulation results.

    .. grid-item-card:: 🔧 **Troubleshooting**
        :link: troubleshooting
        :link-type: doc

        Common issues and solutions for cavity MD
        simulations.

Learning Path
=============

**New to Cavity HOOMD?** Follow this path:

1. **Basic Usage** (:doc:`basic_usage`) - Learn the essential commands
2. **Parameter Sweeps** (:doc:`parameter_sweeps`) - Run systematic studies  
3. **Analysis** (:doc:`analysis`) - Understand your results
4. **Troubleshooting** (:doc:`troubleshooting`) - Solve common problems

Common Workflows
================

**Single Simulation**

.. code-block:: bash

   # Basic cavity simulation
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000

**Parameter Study**

.. code-block:: bash

   # Coupling strength sweep
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3,1e-4,1e-5 --runtime 500

**Control vs Cavity**

.. code-block:: bash

   # With cavity
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000
   
   # Without cavity (control)
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --no-cavity --runtime 1000

**Multiple Replicas**

.. code-block:: bash

   # Run 5 independent replicas
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --replicas "1-5" --runtime 1000

Experiment Types
================

The script provides four pre-configured experiments:

**bussi_langevin_finiteq**
   - Bussi thermostat for molecules
   - Langevin thermostat for cavity  
   - Finite-q cavity mode
   - **Recommended for most studies**

**bussi_langevin_no_finiteq**
   - Same as above but q=0 cavity mode
   - Use for comparison with finite-q

**langevin_langevin**
   - Langevin thermostat for both molecules and cavity
   - Good for stochastic dynamics studies

**bussi_bussi**
   - Bussi thermostat for both molecules and cavity
   - Deterministic temperature control

Getting Help
============

**Command-line help:**

.. code-block:: bash

   python 05_advanced_run.py --help

**Documentation:**

- :doc:`basic_usage` - Command-line reference
- :doc:`parameter_sweeps` - Advanced usage patterns
- :doc:`analysis` - Output file formats
- :doc:`troubleshooting` - Problem solving

**Support:**

- `GitHub Issues <https://github.com/yourusername/cavity-hoomd/issues>`_ - Bug reports
- `GitHub Discussions <https://github.com/yourusername/cavity-hoomd/discussions>`_ - Questions 