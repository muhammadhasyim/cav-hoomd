===============================
Cavity HOOMD Documentation
===============================

**Cavity HOOMD** enables molecular dynamics simulations with optical cavity coupling using HOOMD-blue.

This package provides advanced tools for studying light-matter interactions in molecular systems through cavity quantum electrodynamics (QED) simulations with comprehensive analysis and time-varying coupling capabilities.

Quick Start
===========

Get started quickly with a basic cavity simulation:

.. code-block:: bash

   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd
   ./build_install.sh
   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000

For time-varying coupling experiments:

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 1.0 --runtime 1000

That's it! See the :doc:`installation` and :doc:`quickstart` guides for more details.

What It Does
============

Cavity HOOMD implements cavity-molecule coupling through the single-mode Hamiltonian:

.. math::

   H = \frac{1}{2} K_\lambda \tilde{q}_{0,\lambda}^2 + \tilde{\varepsilon}_{0,\lambda} \tilde{q}_{0,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda} + \frac{\tilde{\varepsilon}_{0,\lambda}^2}{2K_\lambda} \left(\sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}\right)^2

Where:

- :math:`\tilde{q}_{0,\lambda}` is the normalized cavity mode coordinate for polarization λ (x,y)
- :math:`d_{ng,\lambda}` is the dipole moment component of molecule n in direction λ  
- :math:`\tilde{\varepsilon}_{0,\lambda}` is the effective coupling strength for the fundamental mode
- :math:`K_\lambda = m_{0,\lambda}\omega_{0,\lambda}^2` is the cavity spring constant

This implements the single-mode approximation (κ = 0) with transverse polarizations.

**NEW: Time-Varying Coupling**

The plugin now supports time-varying coupling strength and dissipation using HOOMD variants:

.. math::

   g(t) = \begin{cases}
   0 & \text{if } t < t_{\text{switch}} \\
   g_{\text{target}} & \text{if } t \geq t_{\text{switch}}
   \end{cases}

For complete theoretical background, see :doc:`theory`.

Features
========

**Core Capabilities**

- **Time-varying coupling**: Step-function coupling switching for dynamic simulations
- **Multiple thermostats**: Bussi and Langevin thermostats for molecular and cavity degrees of freedom
- **GPU acceleration**: Optimized CUDA kernels for high-performance simulations  
- **Flexible coupling**: Support for both q=0 and finite-q cavity modes
- **HOOMD integration**: Built as native HOOMD-blue plugins with full HOOMD ecosystem compatibility

**Advanced Analysis**

- **Energy tracking**: Comprehensive energy component monitoring (harmonic, coupling, dipole self-energy)
- **Correlation functions**: F(k,t) density correlation analysis
- **Cavity mode monitoring**: Real-time cavity particle position and velocity tracking
- **Performance tracking**: Detailed performance metrics and benchmarking
- **Adaptive timestep**: Automatic timestep optimization for stability and performance

**Simulation Framework**

- **Smart particle handling**: Automatic detection and validation of cavity particles
- **Replica management**: Built-in support for multiple independent simulations
- **SLURM integration**: Native support for high-performance computing environments
- **Comprehensive validation**: Parameter validation following scientific best practices
- **Modular architecture**: Clean separation of forces, analysis, and simulation components

**Time-Varying Experiments**

- **Coupling switching**: Instantaneous coupling activation at specified times
- **Dissipation ramping**: Time-dependent cavity damping
- **Finite-q displacement**: Automatic cavity particle repositioning for finite-q modes
- **Energy conservation**: Rigorous energy conservation testing for time-varying systems

**Output and Visualization**

- **Multiple output formats**: GSD trajectories, energy logs, cavity mode data
- **Real-time monitoring**: Console output with customizable verbosity
- **Analysis tools**: Built-in plotting and analysis utilities
- **Jupyter integration**: Interactive notebook examples for analysis

Documentation
=============

.. toctree::
   :maxdepth: 2

   installation
   quickstart
   theory
   api/index
   license

Help and Support
================

- **Installation help**: :doc:`installation`
- **Usage examples**: :doc:`quickstart`  
- **Issues**: https://github.com/muhammadhasyim/cav-hoomd/issues
- **Discussions**: https://github.com/muhammadhasyim/cav-hoomd/discussions

Citation
========

.. code-block:: bibtex

   @software{cavity_hoomd_2025,
     title={Cavity HOOMD: Molecular Dynamics with Optical Cavity Coupling},
     author={Development Team},
     year={2025},
     url={https://github.com/muhammadhasyim/cav-hoomd}
   } 