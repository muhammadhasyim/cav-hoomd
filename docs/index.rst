===============================
Cavity HOOMD User Guide
===============================

**Cavity HOOMD** enables molecular dynamics simulations with optical cavity coupling using HOOMD-blue.

This package extends HOOMD-blue to simulate molecules interacting with optical cavity modes, enabling study of light-matter interactions, polariton chemistry, and cavity-modified molecular dynamics.

.. image:: _static/logo.png
   :width: 200px
   :align: center
   :alt: Cavity HOOMD Logo

Quick Links
===========

- **New Users:** Start with :doc:`part1_application/getting_started`
- **Running Simulations:** See :doc:`part1_application/running_simulations`
- **Theory Background:** Read :doc:`part2_theory/introduction`
- **API Reference:** Browse :doc:`api/index`
- **Installation Help:** :doc:`installation`

What Can You Do?
================

**Core Capabilities:**

- Simulate cavity-molecule coupling with time-dependent coupling strength
- Study strong coupling effects and polariton formation
- Analyze energy flow between molecular and photonic degrees of freedom
- Use advanced control systems for temperature and parameter management
- GPU-accelerated simulations for large systems

**Example Applications:**

- Vibrational strong coupling in molecular ensembles
- Cavity-modified chemical dynamics
- Non-equilibrium switching experiments
- Polariton-mediated energy transfer
- Mode-selective chemistry under strong coupling

Quick Example
=============

Get started in minutes:

.. code-block:: bash

   # Install
   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd
   ./build_install.sh

   # Run your first cavity simulation
   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000

That's it! See :doc:`part1_application/getting_started` for detailed instructions.

Documentation Structure
=======================

This documentation is organized into three main parts, following the structure of OpenMM:

**Part I: The Cavity HOOMD Application Layer**
  User-focused tutorials and guides for running simulations, from basic examples to advanced time-varying coupling experiments.

**Part II: Theory and Implementation**
  Mathematical foundations, equations of motion, energy conservation, and physical interpretation of cavity-molecule coupling.

**Part III: Advanced Features and Analysis**
  Advanced controllers, FDR temperature measurement, molecular temperature decomposition, correlation analysis, and performance optimization.

----

.. toctree::
   :maxdepth: 2
   :caption: Part I: Application Layer
   :hidden:
   
   part1_application/getting_started
   part1_application/running_simulations
   part1_application/advanced_simulations
   part1_application/analysis_tools
   part1_application/time_varying_coupling

.. toctree::
   :maxdepth: 2
   :caption: Part II: Theory and Implementation
   :hidden:
   
   part2_theory/introduction
   part2_theory/cavity_forces
   part2_theory/time_varying_coupling
   part2_theory/thermostats
   part2_theory/energy_conservation
   part2_theory/strong_coupling

.. toctree::
   :maxdepth: 2
   :caption: Part III: Advanced Features
   :hidden:
   
   part3_advanced/controllers
   part3_advanced/fdr_temperature
   part3_advanced/molecular_temperatures
   part3_advanced/correlation_analysis
   part3_advanced/performance

.. toctree::
   :maxdepth: 1
   :caption: Reference
   :hidden:
   
   installation
   reference
   api/index
   license

Help and Support
================

- **GitHub Issues:** https://github.com/muhammadhasyim/cav-hoomd/issues
- **Discussions:** https://github.com/muhammadhasyim/cav-hoomd/discussions
- **Documentation:** You're reading it!

Citation
========

If you use Cavity HOOMD in your research, please cite:

.. code-block:: bibtex

   @software{cavity_hoomd_2025,
     title={Cavity HOOMD: Molecular Dynamics with Optical Cavity Coupling},
     author={Development Team},
     year={2025},
     url={https://github.com/muhammadhasyim/cav-hoomd}
   }

License
=======

Cavity HOOMD is released under the BSD 3-Clause License. See :doc:`license` for details.

Indices and Tables
==================

* :ref:`genindex`
* :ref:`search`
