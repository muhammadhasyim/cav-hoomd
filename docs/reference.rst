===============
Quick Reference
===============

This page provides quick reference tables, command-line flags, units, and a glossary for Cavity HOOMD.

.. contents:: In this Reference
   :local:
   :depth: 2

Command-Line Flags
==================

05_advanced_run.py Options
---------------------------

**Basic Parameters:**

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Flag
     - Default
     - Description
   * - ``--coupling VALUE``
     - 1e-3
     - Coupling strength (atomic units)
   * - ``--temperature VALUE``
     - 100
     - Temperature in Kelvin
   * - ``--frequency VALUE``
     - 2000
     - Cavity frequency in cm⁻¹
   * - ``--runtime VALUE``
     - 500
     - Simulation time in ps

**Thermostats:**

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Flag
     - Default
     - Description
   * - ``--molecular-bath``
     - bussi
     - Molecular thermostat: bussi, langevin, none
   * - ``--cavity-bath``
     - langevin
     - Cavity thermostat: bussi, langevin, none

**Modes:**

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Flag
     - Default
     - Description
   * - ``--finite-q``
     - False
     - Enable finite-q cavity mode
   * - ``--no-cavity``
     - False
     - Run without cavity (control)

**Time-Varying:**

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Flag
     - Default
     - Description
   * - ``--switch-time VALUE``
     - None
     - Time to switch coupling on (ps)

**Execution:**

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Flag
     - Default
     - Description
   * - ``--device``
     - auto
     - Device: GPU, CPU, auto
   * - ``--replicas``
     - 1
     - Replica IDs (e.g., "1-5" or "1,3,5")

Parameter Tables
================

Typical Values
--------------

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Parameter
     - Range
     - Notes
   * - Coupling strength
     - 10⁻⁵ - 10⁻³
     - Atomic units; strong coupling > 10⁻³
   * - Temperature
     - 50-500 K
     - Typical molecular simulations
   * - Cavity frequency
     - 1000-4000 cm⁻¹
     - Match molecular vibrations
   * - Timestep
     - 0.0005-0.002 ps
     - Smaller for high frequencies
   * - Runtime
     - 100-10000 ps
     - Depends on observable

Unit Conversions
================

Energy
------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Unit
     - Conversion
   * - 1 Hartree
     - 27.2 eV = 627 kcal/mol = 4.36×10⁻¹⁸ J
   * - 1 eV
     - 0.0367 Hartree = 23.06 kcal/mol
   * - 1 kcal/mol
     - 0.00159 Hartree = 0.043 eV

Temperature
-----------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Temperature
     - Reduced Units
   * - 50 K
     - 0.00016 Hartree/k_B
   * - 100 K
     - 0.00032 Hartree/k_B
   * - 300 K
     - 0.00095 Hartree/k_B

Length
------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Unit
     - Conversion
   * - 1 Bohr
     - 0.529 Å = 0.0529 nm
   * - 1 Å
     - 1.889 Bohr
   * - 1 nm
     - 18.89 Bohr

Time
----

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Unit
     - Conversion
   * - 1 atomic time unit
     - 0.024 fs = 2.42×10⁻¹⁷ s
   * - 1 fs
     - 41.34 atomic time units
   * - 1 ps
     - 41,341 atomic time units

Frequency
---------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Unit
     - Conversion
   * - 1 cm⁻¹
     - 4.556×10⁻⁶ Hartree/ℏ = 0.030 meV
   * - 1000 cm⁻¹
     - 0.00456 Hartree/ℏ = 124 meV
   * - 2000 cm⁻¹
     - 0.00916 Hartree/ℏ = 248 meV

Output Files
============

File Formats
------------

**trajectory.gsd:**

- Binary GSD format
- Contains: positions, velocities, forces, box
- Read with: ``gsd.hoomd.open()``

**energy_tracker.txt:**

- Tab-separated values
- Columns: time, kinetic, potential, cavity, coupling, self, total
- Units: ps for time, reduced units for energy

**cavity_mode.txt:**

- Tab-separated values
- Columns: time, q_x, q_y, v_x, v_y
- Units: ps for time, reduced units for positions/velocities

**molecular_temps.csv:**

- Comma-separated values
- Columns: time, T_trans, T_rot, T_vib, T_total
- Units: ps for time, Kelvin for temperatures

Glossary
========

.. glossary::

   Cavity Mode
      Quantized electromagnetic field mode confined in optical cavity

   Coupling Strength
      Parameter g controlling interaction between molecular dipoles and cavity field

   Collective Coupling
      Enhanced coupling :math:`g_{\text{eff}} = g\sqrt{N}` for N molecules

   Finite-q Mode
      Cavity mode with non-zero wave vector, spatial phase variation

   FDR
      Fluctuation-Dissipation Ratio, measures effective temperature

   Polariton
      Hybrid light-matter quasi-particle from strong coupling

   Rabi Splitting
      Energy splitting :math:`\Omega_R` between upper and lower polaritons

   q=0 Mode
      Uniform cavity mode with zero wave vector, all molecules couple in phase

   Self-Energy
      Dipole-dipole interaction energy term :math:`g^2 D^2/2K`

   Strong Coupling
      Regime where :math:`\Omega_R > \sqrt{\gamma\kappa}`, polaritons form

Symbol Reference
================

Mathematical Symbols
--------------------

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Symbol
     - Meaning
   * - :math:`g`
     - Coupling strength parameter
   * - :math:`\omega_c`
     - Cavity angular frequency
   * - :math:`q_\lambda`
     - Cavity mode coordinate (polarization λ)
   * - :math:`D_\lambda`
     - Total molecular dipole moment (direction λ)
   * - :math:`d_n`
     - Dipole moment of molecule n
   * - :math:`K`
     - Cavity spring constant :math:`m\omega^2`
   * - :math:`\tilde{\varepsilon}`
     - Effective coupling strength
   * - :math:`N`
     - Number of molecules
   * - :math:`\Omega_R`
     - Rabi splitting / Rabi frequency
   * - :math:`\gamma`
     - Molecular dephasing rate
   * - :math:`\kappa`
     - Cavity photon loss rate
   * - :math:`\tau`
     - Thermostat coupling time
   * - :math:`T`
     - Temperature
   * - :math:`k_B`
     - Boltzmann constant

Citations
=========

If you use Cavity HOOMD, please cite:

.. code-block:: bibtex

   @software{cavity_hoomd_2025,
     title={Cavity HOOMD: Molecular Dynamics with Optical Cavity Coupling},
     author={Development Team},
     year={2025},
     url={https://github.com/muhammadhasyim/cav-hoomd}
   }

Also cite HOOMD-blue:

.. code-block:: bibtex

   @article{anderson2020hoomd,
     title={HOOMD-blue: A Python package for high-performance molecular dynamics and hard particle Monte Carlo simulations},
     author={Anderson, Joshua A and Glaser, Jens and Glotzer, Sharon C},
     journal={Computational Materials Science},
     volume={173},
     pages={109363},
     year={2020},
     publisher={Elsevier}
   }

Related Documentation
=====================

- :doc:`part1_application/getting_started` - Getting started guide
- :doc:`part2_theory/introduction` - Theory introduction
- :doc:`api/index` - API reference
- :doc:`installation` - Installation instructions

