======
Theory
======

Overview
========

Cavity HOOMD implements cavity quantum electrodynamics (cQED) within classical molecular dynamics simulations. This approach enables the study of light-matter interactions in the strong coupling regime, where molecular excitations hybridize with cavity photons to form polaritons.

Theoretical Framework
=====================

Cavity-Coupled Hamiltonian
---------------------------

The total Hamiltonian for the cavity-molecule system is:

.. math::

   H = H_{\text{mol}} + H_{\text{cav}} + H_{\text{int}}

where:

* :math:`H_{\text{mol}}` is the molecular Hamiltonian
* :math:`H_{\text{cav}}` is the cavity Hamiltonian  
* :math:`H_{\text{int}}` is the light-matter interaction

Molecular Hamiltonian
---------------------

The molecular system follows standard molecular dynamics:

.. math::

   H_{\text{mol}} = \sum_i \frac{|\mathbf{p}_i|^2}{2m_i} + V(\{\mathbf{r}_i\})

where :math:`\mathbf{p}_i` and :math:`\mathbf{r}_i` are momenta and positions of particle :math:`i`, and :math:`V` is the inter-particle potential.

Cavity Hamiltonian
-------------------

The cavity mode is treated as a quantum harmonic oscillator:

.. math::

   H_{\text{cav}} = \hbar\omega_c a^\dagger a + \frac{1}{2}\hbar\omega_c

where :math:`\omega_c` is the cavity frequency and :math:`a^\dagger, a` are creation/annihilation operators.

Light-Matter Interaction
-------------------------

In the electric dipole approximation:

.. math::

   H_{\text{int}} = -\boldsymbol{\mu} \cdot \mathbf{E}_{\text{cav}}

where :math:`\boldsymbol{\mu} = \sum_i q_i \mathbf{r}_i` is the total molecular dipole moment and :math:`\mathbf{E}_{\text{cav}}` is the cavity electric field.

For a single-mode cavity:

.. math::

   \mathbf{E}_{\text{cav}} = i\mathbf{e} E_0 (a^\dagger - a)

where :math:`\mathbf{e}` is the polarization vector and :math:`E_0 = \sqrt{\frac{\hbar\omega_c}{2\epsilon_0 V_{\text{cav}}}}`.

Classical Implementation
========================

Cavity Mode Dynamics
---------------------

The cavity mode is treated as a classical harmonic oscillator with coordinate :math:`q` and momentum :math:`p`:

.. math::

   H_{\text{cav}} = \frac{p^2}{2} + \frac{1}{2}\omega_c^2 q^2

The equations of motion are:

.. math::

   \dot{q} &= p \\
   \dot{p} &= -\omega_c^2 q - g \boldsymbol{\mu} \cdot \mathbf{e}

where :math:`g` is the coupling strength.

Molecular Dynamics
------------------

The molecular particles experience the usual forces plus the cavity coupling:

.. math::

   m_i \ddot{\mathbf{r}}_i = \mathbf{F}_i^{\text{mol}} + q_i g \mathbf{e} q_{\text{cav}}

where :math:`\mathbf{F}_i^{\text{mol}}` are the intermolecular forces.

Strong Coupling Regime
======================

Rabi Splitting
--------------

In the strong coupling regime, the cavity-molecule system exhibits Rabi splitting:

.. math::

   \Omega_R = 2g\sqrt{N}

where :math:`N` is the number of molecules and :math:`g` is the single-molecule coupling strength.

Polariton Formation
-------------------

The eigenstates become hybridized polaritons:

.. math::

   |P_\pm\rangle = \frac{1}{\sqrt{2}}(|0,e\rangle \pm |1,g\rangle)

with energies :math:`E_\pm = \frac{\omega_c + \omega_m}{2} \pm \frac{\Omega_R}{2}`.

Thermostats and Equilibration
==============================

Molecular Thermostats
---------------------

The molecular degrees of freedom are thermostated using:

* **Bussi thermostat**: Stochastic velocity rescaling
* **Langevin thermostat**: Friction and random forces

Cavity Thermostats
------------------

The cavity mode can be thermostated to maintain thermal equilibrium:

.. math::

   \dot{p} = -\omega_c^2 q - g \boldsymbol{\mu} \cdot \mathbf{e} - \gamma p + \sqrt{2\gamma k_B T} \xi(t)

where :math:`\gamma` is the damping coefficient and :math:`\xi(t)` is white noise.

Finite-q Effects
=================

Wave Vector Dependence
----------------------

For finite wave vector :math:`\mathbf{k} \neq 0`, the interaction becomes spatially dependent:

.. math::

   H_{\text{int}} = -g \sum_i q_i e^{i\mathbf{k} \cdot \mathbf{r}_i} (a^\dagger + a)

This enables study of spatially correlated phenomena and collective excitations.

Observables and Analysis
========================

Energy Components
-----------------

The total energy is conserved and consists of:

.. math::

   E_{\text{total}} = E_{\text{kinetic}} + E_{\text{potential}} + E_{\text{cavity}} + E_{\text{coupling}}

Correlation Functions
---------------------

Important correlation functions include:

* **Density correlation**: :math:`F(k,t) = \langle \rho_{\mathbf{k}}(t) \rho_{-\mathbf{k}}(0) \rangle`
* **Dipole correlation**: :math:`C_{\mu}(t) = \langle \boldsymbol{\mu}(t) \cdot \boldsymbol{\mu}(0) \rangle`
* **Cavity correlation**: :math:`C_q(t) = \langle q(t) q(0) \rangle`

Applications
============

Chemical Reactions
------------------

Cavity coupling can modify reaction rates and pathways through:

* Vibrational strong coupling
* Reaction coordinate modification
* Collective effects

Phase Transitions
-----------------

Strong coupling can induce new phases:

* Polariton condensation
* Modified phase diagram
* Collective ordering

Spectroscopy
------------

The framework enables simulation of:

* Absorption spectra
* Raman scattering
* Polariton dispersion

References
==========

Key theoretical works:

1. Jaynes, E. T. & Cummings, F. W. Comparison of quantum and semiclassical radiation theories. *Proc. IEEE* **51**, 89-109 (1963).

2. Hutchison, J. A., Schwartz, T., Genet, C., Devaux, E. & Ebbesen, T. W. Modifying chemical landscapes by coupling to vacuum fields. *Angew. Chem. Int. Ed.* **51**, 1592-1596 (2012).

3. Flick, J., Ruggenthaler, M., Appel, H. & Rubio, A. Atoms and molecules in cavities: from weak to strong coupling. *Proc. Natl. Acad. Sci. USA* **114**, 3026-3034 (2017).

4. Galego, J., Garcia-Vidal, F. J. & Feist, J. Cavity-induced modifications of molecular structure in the strong-coupling regime. *Phys. Rev. X* **5**, 041022 (2015). 