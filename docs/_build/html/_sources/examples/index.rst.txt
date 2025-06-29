========
Examples
========

This section contains practical examples demonstrating various features of Cavity HOOMD.

Available Examples
==================

Basic Examples
--------------

* **Basic Simulation** - Minimal cavity-coupled MD setup
* **Energy Tracking** - Monitoring energy conservation  
* **Parameter Sweep** - Systematic parameter exploration

Analysis Examples  
-----------------

* **Correlation Analysis** - F(k,t) and dipole correlations
* **Trajectory Analysis** - Post-processing GSD files

Example Files
=============

Complete example scripts are available in the ``examples/`` directory:

* ``01_basic_cavity.py`` - Basic cavity coupling
* ``02_energy_tracking.py`` - Energy conservation
* ``03_parameter_sweep.py`` - Parameter exploration  
* ``04_correlation_analysis.py`` - Correlation functions
* ``05_advanced_run.py`` - Complete simulation framework

Running Examples
================

.. code-block:: bash

   cd examples/
   python 01_basic_cavity.py

Each example includes detailed comments explaining the setup and analysis.

.. note::
   
   Make sure you have a molecular system GSD file (``init-0.gsd``) in the 
   examples directory before running the scripts. 