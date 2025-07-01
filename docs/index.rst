===============================
Cavity HOOMD Documentation
===============================

**Cavity HOOMD** enables molecular dynamics simulations with optical cavity coupling using HOOMD-blue.

This package provides tools for studying light-matter interactions in molecular systems through cavity quantum electrodynamics (QED) simulations.

Quick Start
===========

Get started quickly with a basic cavity simulation:

.. code-block:: bash

   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd
   ./build_install.sh
   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000

That's it! See the :doc:`installation` and :doc:`quickstart` guides for more details.

What It Does
============

Cavity HOOMD implements cavity-molecule coupling through the Hamiltonian:

.. math::

   H = \frac{1}{2} K q^2 + g \vec{q} \cdot \vec{d} + \frac{g^2}{2K} d^2

Where :math:`q` is the cavity mode position, :math:`d` is the molecular dipole moment, 
:math:`g` is the coupling strength, and :math:`K` is the cavity spring constant.

For theoretical background, see :doc:`theory`.

Features
========

- **Multiple thermostats**: Bussi and Langevin thermostats for molecular and cavity degrees of freedom
- **GPU acceleration**: CUDA support for high-performance simulations  
- **Advanced analysis**: Energy tracking, correlation functions, and cavity mode monitoring
- **Flexible coupling**: Support for both q=0 and finite-q cavity modes
- **HOOMD integration**: Built as native HOOMD-blue plugins

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