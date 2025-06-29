===============================
Cavity HOOMD Documentation
===============================

.. image:: _static/logo.png
   :width: 200px
   :align: center

**Cavity HOOMD** is a powerful molecular dynamics simulation framework that extends HOOMD-blue with advanced cavity-coupling capabilities. It enables researchers to simulate molecular systems coupled to optical cavity modes, providing insights into cavity-mediated chemical reactions, polariton dynamics, and light-matter interactions.

.. grid:: 2
    :gutter: 3

    .. grid-item-card:: 🚀 **Quick Start**
        :link: quickstart
        :link-type: doc

        Get up and running with Cavity HOOMD in minutes. Learn the basics
        of setting up and running your first cavity-coupled simulation.

    .. grid-item-card:: 📚 **User Guide**
        :link: user_guide/index
        :link-type: doc

        Comprehensive tutorials and guides for using Cavity HOOMD effectively.
        From simple simulations to advanced parameter sweeps.

    .. grid-item-card:: 🔧 **API Reference**
        :link: api/index
        :link-type: doc

        Complete API documentation for all classes, functions, and modules
        in the Cavity HOOMD framework.

    .. grid-item-card:: 🧪 **Examples**
        :link: examples/index
        :link-type: doc

        Real-world examples and Jupyter notebooks demonstrating various
        simulation scenarios and analysis techniques.

Key Features
============

🔬 **Advanced Cavity Coupling**
   - Optical cavity modes with tunable frequency and coupling strength
   - Support for both finite-q and q=0 cavity modes
   - Dipole-cavity interaction with self-energy corrections

⚡ **High Performance**
   - GPU acceleration with CUDA support
   - Optimized C++ kernels with Python fallbacks
   - Adaptive timestep control for efficiency

🎛️ **Flexible Thermostats**
   - Bussi thermostat for accurate temperature control
   - Langevin dynamics for stochastic sampling
   - Independent thermostats for molecular and cavity degrees of freedom

📊 **Comprehensive Analysis**
   - Real-time energy tracking and conservation monitoring
   - F(k,t) density correlation functions
   - Cavity mode amplitude and phase tracking
   - Extensive logging and performance metrics

🔄 **Robust Simulation Framework**
   - Parameter sweep capabilities
   - Automatic error handling and recovery
   - Checkpoint and restart functionality
   - Integration with SLURM for HPC workflows

Installation
============

.. code-block:: bash

   # Install from PyPI (recommended)
   pip install cavity-hoomd

   # Or install from source
   git clone https://github.com/yourusername/cavity-hoomd.git
   cd cavity-hoomd
   pip install -e .

Quick Example
=============

Here's a simple example of setting up a cavity-coupled molecular dynamics simulation:

.. code-block:: python

   import hoomd
   from hoomd.cavitymd import CavityForce, CavityMDSimulation

   # Create a cavity MD simulation
   sim = CavityMDSimulation(
       job_dir="./output",
       replica=1,
       freq=2000.0,        # Cavity frequency in cm⁻¹
       couplstr=1e-3,      # Coupling strength
       incavity=True,      # Enable cavity coupling
       runtime_ps=1000.0,  # Simulation time in ps
       temperature=300.0,  # Temperature in K
       molecular_thermostat='bussi',
       cavity_thermostat='langevin'
   )

   # Run the simulation
   sim.run()

What's New
==========

.. admonition:: Version 1.0.0 - Latest Release
   :class: tip

   - **New**: Adaptive timestep control with automatic error tolerance
   - **New**: Comprehensive energy tracking with reservoir monitoring
   - **New**: F(k,t) density correlation analysis
   - **Improved**: GPU acceleration with CUDA kernels
   - **Improved**: Enhanced parameter sweep capabilities
   - **Fixed**: Energy conservation in cavity-coupled systems

Scientific Background
=====================

Cavity HOOMD implements the theoretical framework for cavity quantum electrodynamics (cQED) in classical molecular dynamics simulations. The core Hamiltonian includes:

.. math::

   H = H_{\text{mol}} + \frac{1}{2}\hbar\omega_c (a^\dagger a + \frac{1}{2}) + g \vec{d} \cdot \vec{E}_{\text{cav}}

Where:
- :math:`H_{\text{mol}}` is the molecular Hamiltonian
- :math:`\omega_c` is the cavity frequency
- :math:`g` is the coupling strength
- :math:`\vec{d}` is the molecular dipole moment
- :math:`\vec{E}_{\text{cav}}` is the cavity electric field

This enables the study of:
- Cavity-mediated chemical reactions
- Polariton formation and dynamics
- Collective strong coupling effects
- Light-matter energy transfer

Contents
========

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart
   theory

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/index
   user_guide/basic_usage
   user_guide/advanced_features
   user_guide/parameter_sweeps
   user_guide/analysis
   user_guide/performance

.. toctree::
   :maxdepth: 2
   :caption: Examples

   examples/index
   examples/basic_simulation
   examples/energy_tracking
   examples/parameter_sweep
   examples/analysis_tutorial

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/index
   api/forces
   api/simulation
   api/analysis
   api/utils

.. toctree::
   :maxdepth: 2
   :caption: Development

   development/index
   development/contributing
   development/testing
   development/documentation

.. toctree::
   :maxdepth: 1
   :caption: Reference

   changelog
   license
   references

Support & Community
===================

.. grid:: 2
    :gutter: 3

    .. grid-item-card:: 🐛 **Bug Reports**
        :link: https://github.com/yourusername/cavity-hoomd/issues
        :link-type: url

        Found a bug? Report it on our GitHub issue tracker.

    .. grid-item-card:: 💬 **Discussions**
        :link: https://github.com/yourusername/cavity-hoomd/discussions
        :link-type: url

        Ask questions and discuss features with the community.

    .. grid-item-card:: 📖 **Documentation**
        :link: https://cavity-hoomd.readthedocs.io
        :link-type: url

        Complete documentation with examples and tutorials.

    .. grid-item-card:: 📧 **Contact**
        :link: mailto:cavity-hoomd@example.com
        :link-type: url

        Get in touch with the development team.

Citation
========

If you use Cavity HOOMD in your research, please cite:

.. code-block:: bibtex

   @article{cavity_hoomd_2025,
     title={Cavity HOOMD: A Framework for Cavity-Coupled Molecular Dynamics},
     author={Development Team},
     journal={Journal of Computational Chemistry},
     year={2025},
     volume={XX},
     pages={XXX-XXX},
     doi={10.1002/jcc.XXXXX}
   }

License
=======

Cavity HOOMD is released under the BSD 3-Clause License. See the :doc:`license` page for details.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search` 