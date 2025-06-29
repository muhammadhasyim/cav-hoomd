==============
API Reference
==============

This section provides comprehensive API documentation for all modules, classes, and functions in Cavity HOOMD.

.. toctree::
   :maxdepth: 2

   forces
   simulation
   analysis
   utils

Overview
========

The Cavity HOOMD API is organized into several key modules:

**Forces** (:doc:`forces`)
   Classes for cavity-matter interactions and force calculations

**Simulation** (:doc:`simulation`)
   High-level simulation management and control

**Analysis** (:doc:`analysis`)
   Data analysis, tracking, and correlation functions

**Utilities** (:doc:`utils`)
   Helper functions, constants, and utility classes

Module Structure
================

.. autosummary::
   :toctree: _autosummary
   :recursive:

   hoomd.cavitymd.forces
   hoomd.cavitymd.simulation
   hoomd.cavitymd.analysis
   hoomd.cavitymd.utils

Quick Reference
===============

Most Common Classes
-------------------

.. autosummary::
   :nosignatures:

   hoomd.cavitymd.CavityForce
   hoomd.cavitymd.AdaptiveTimestepUpdater
   hoomd.cavitymd.EnergyTracker
   hoomd.cavitymd.CavityModeTracker
   hoomd.cavitymd.FieldAutocorrelationTracker

Key Functions
-------------

.. autosummary::
   :nosignatures:

   hoomd.cavitymd.unwrap_positions
   hoomd.cavitymd.PhysicalConstants

Constants and Units
===================

Cavity HOOMD uses atomic units internally with conversion utilities for common units:

- **Energy**: Hartree (atomic units)
- **Length**: Bohr (atomic units)
- **Time**: Atomic time units (≈ 24.2 attoseconds)
- **Temperature**: Kelvin
- **Frequency**: cm⁻¹ for cavity frequencies

See :class:`hoomd.cavitymd.PhysicalConstants` for conversion factors. 