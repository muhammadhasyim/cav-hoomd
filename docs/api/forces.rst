======
Forces
======

This module provides force implementations for cavity-matter interactions.

.. currentmodule:: hoomd.cavitymd

Core Force Classes
==================

CavityForce
-----------

.. autoclass:: CavityForce
   :members:
   :inherited-members:
   :show-inheritance:

   .. rubric:: Methods

   .. autosummary::
      :nosignatures:

      ~CavityForce.__init__
      ~CavityForce.harmonic_energy
      ~CavityForce.coupling_energy  
      ~CavityForce.dipole_self_energy
      ~CavityForce.total_cavity_energy

   .. rubric:: Properties

   .. autosummary::
      :nosignatures:

      ~CavityForce.implementation
      ~CavityForce.energy
      ~CavityForce.forces

Python Implementation
=====================

.. automodule:: hoomd.cavitymd.cavity_force_python
   :members:
   :show-inheritance:

Usage Examples
==============

Basic Cavity Force Setup
-------------------------

.. code-block:: python

   import numpy as np
   from hoomd.cavitymd import CavityForce

   # Create a cavity force with specified parameters
   cavity_force = CavityForce(
       kvector=[0, 0, 1],      # Wave vector direction
       couplstr=1e-3,          # Coupling strength
       omegac=0.092,           # Cavity frequency (atomic units)
       phmass=1.0              # Photon mass
   )

   # Add to integrator forces
   forces = [cavity_force, ...]
   integrator = hoomd.md.Integrator(dt=0.001, forces=forces)

Energy Components
-----------------

The cavity force provides access to individual energy components:

.. code-block:: python

   # After running simulation
   harmonic_energy = cavity_force.harmonic_energy
   coupling_energy = cavity_force.coupling_energy  
   dipole_self_energy = cavity_force.dipole_self_energy
   total_energy = cavity_force.total_cavity_energy

   print(f"Harmonic: {harmonic_energy:.6f}")
   print(f"Coupling: {coupling_energy:.6f}")
   print(f"Self-energy: {dipole_self_energy:.6f}")
   print(f"Total: {total_energy:.6f}")

Implementation Selection
------------------------

The force automatically selects the best available implementation:

.. code-block:: python

   # Check which implementation is being used
   print(f"Using {cavity_force.implementation} implementation")

   # Force Python implementation for debugging
   cavity_force_python = CavityForce(
       kvector=[0, 0, 1],
       couplstr=1e-3,
       omegac=0.092,
       force_python=True  # Force Python implementation
   )

Theory
======

Hamiltonian
-----------

The cavity force implements the interaction Hamiltonian:

.. math::

   H_{cav} = \frac{1}{2}K q^2 + g \vec{q} \cdot \vec{d} + \frac{g^2}{2K} d^2

where:

- :math:`K = m_{\text{ph}} \omega_c^2` is the cavity spring constant
- :math:`q` is the cavity mode coordinate 
- :math:`\vec{d}` is the molecular dipole moment
- :math:`g` is the coupling strength

Energy Components
-----------------

**Harmonic Energy**

.. math::

   E_{\text{harm}} = \frac{1}{2}K q^2

The harmonic oscillator energy of the cavity mode.

**Coupling Energy**

.. math::

   E_{\text{coup}} = g \vec{q} \cdot \vec{d}

The direct dipole-cavity interaction energy.

**Dipole Self-Energy**  

.. math::

   E_{\text{self}} = \frac{g^2}{2K} d^2

The dipole self-energy correction due to cavity coupling.

Forces
------

The forces on particles are derived from the energy gradients:

.. math::

   \vec{F}_i = -\frac{\partial H_{cav}}{\partial \vec{r}_i}

For molecular particles:

.. math::

   \vec{F}_i = -g q_i \vec{q} - \frac{g^2}{K} q_i \vec{d}

For the cavity particle:

.. math::

   \vec{F}_{cav} = -K \vec{q} - g \vec{d} 