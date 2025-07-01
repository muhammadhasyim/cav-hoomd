===============================
Cavity HOOMD Documentation
===============================

.. image:: _static/logo.png
   :width: 200px
   :align: center

**Cavity HOOMD** is a molecular dynamics simulation framework that extends HOOMD-blue with cavity-coupling capabilities. It enables researchers to simulate molecular systems coupled to optical cavity modes.

The easiest way to run simulations is using the **05_advanced_run.py** script, which provides a complete command-line interface for cavity MD simulations.

.. grid:: 2
    :gutter: 3

    .. grid-item-card:: 🚀 **Quick Start**
        :link: quickstart
        :link-type: doc

        Get up and running with Cavity HOOMD using 05_advanced_run.py.
        Run your first simulation in minutes.

    .. grid-item-card:: 📚 **User Guide**
        :link: user_guide/index
        :link-type: doc

        Learn how to use 05_advanced_run.py for different simulation scenarios,
        parameter sweeps, and analysis.

Key Features
============

🔬 **Complete Command-Line Interface**
   - Single script handles all simulation types
   - Built-in parameter sweeps and replica management
   - Automatic output organization

⚡ **Advanced Simulation Capabilities**
   - Multiple thermostat combinations (Bussi, Langevin)
   - Adaptive timestep control
   - GPU and CPU support

📊 **Comprehensive Analysis**
   - Real-time energy tracking
   - F(k,t) density correlation functions
   - Cavity mode monitoring

🔄 **Production-Ready Features**
   - SLURM integration for HPC
   - Robust error handling
   - Detailed logging and performance metrics

Quick Example
=============

Run a cavity-coupled simulation:

.. code-block:: bash

   # Basic cavity simulation
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3 --runtime 1000

   # Parameter sweep over coupling strengths
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --coupling 1e-3,1e-4,1e-5 --runtime 500

   # Run without cavity (control simulation)
   python 05_advanced_run.py --experiment bussi_langevin_finiteq --no-cavity --runtime 1000

Available Experiments
=====================

The script provides four pre-configured experiment types:

- **bussi_langevin_finiteq**: Bussi thermostat for molecules, Langevin for cavity (finite-q mode)
- **bussi_langevin_no_finiteq**: Same as above but q=0 mode  
- **langevin_langevin**: Langevin thermostat for both molecules and cavity
- **bussi_bussi**: Bussi thermostat for both molecules and cavity

Installation
============

.. code-block:: bash

   # Install from source
   git clone https://github.com/yourusername/cavity-hoomd.git
   cd cavity-hoomd
   pip install -e .

Scientific Background
=====================

Cavity HOOMD implements cavity quantum electrodynamics (cQED) in classical MD simulations:

.. math::

   H = H_{\text{mol}} + \frac{1}{2}\hbar\omega_c (a^\dagger a + \frac{1}{2}) + g \vec{d} \cdot \vec{E}_{\text{cav}}

This enables study of cavity-mediated chemical reactions, polariton dynamics, and light-matter interactions.

Contents
========

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/index

.. toctree::
   :maxdepth: 1
   :caption: Reference

   theory
   api/index
   license

Citation
========

If you use Cavity HOOMD in your research, please cite:

.. code-block:: bibtex

   @article{cavity_hoomd_2025,
     title={Cavity HOOMD: A Framework for Cavity-Coupled Molecular Dynamics},
     author={Development Team},
     journal={Journal of Computational Chemistry},
     year={2025}
   }

Support
=======

- **Documentation**: https://cavity-hoomd.readthedocs.io
- **Issues**: https://github.com/yourusername/cavity-hoomd/issues
- **Discussions**: https://github.com/yourusername/cavity-hoomd/discussions

License
=======

Cavity HOOMD is released under the BSD 3-Clause License.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search` 