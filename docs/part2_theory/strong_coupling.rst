================
Strong Coupling
================

This section describes the strong coupling regime where cavity and molecular excitations hybridize to form polaritons.

.. contents:: In this Section
   :local:
   :depth: 2

Defining Strong Coupling
=========================

Collective Enhancement
----------------------

For N molecules coupling to a single cavity mode:

.. math::

   g_{\text{eff}} = g_{\text{single}} \sqrt{N}

Where:

- :math:`g_{\text{single}}`: Single-molecule coupling strength
- :math:`N`: Number of molecules
- :math:`g_{\text{eff}}`: Effective collective coupling

**Example:**

- N = 512 molecules
- :math:`g_{\text{single}} = 10^{-3}`
- :math:`g_{\text{eff}} = 10^{-3} \times \sqrt{512} = 0.0226`

This :math:`\sqrt{N}` enhancement is key to achieving strong coupling.

Strong Coupling Criterion
--------------------------

**Condition:**

.. math::

   \Omega_R > \sqrt{\gamma_{\text{mol}} \kappa_{\text{cav}}}

Where:

- :math:`\Omega_R = 2g_{\text{eff}}`: Rabi splitting
- :math:`\gamma_{\text{mol}}`: Molecular dephasing rate
- :math:`\kappa_{\text{cav}}`: Cavity photon loss rate

**Physical meaning:** Coherent energy exchange faster than decoherence.

Polariton Formation
===================

Hybrid Light-Matter States
---------------------------

In strong coupling, cavity photons and molecular excitations hybridize:

.. math::

   |\text{polariton}\rangle = \alpha |\text{photon}\rangle + \beta |\text{molecule}\rangle

Creating two new eigenstates:

- **Upper polariton (UP):** :math:`E_+ = E_0 + \frac{\Omega_R}{2}`
- **Lower polariton (LP):** :math:`E_- = E_0 - \frac{\Omega_R}{2}`

Where :math:`E_0` is the uncoupled resonance energy.

Rabi Splitting
--------------

**Energy splitting:**

.. math::

   \Omega_R = 2g_{\text{eff}} = 2g_{\text{single}}\sqrt{N}

**Observable in:**

- Transmission/reflection spectra
- Photoluminescence
- Time-resolved dynamics

**Typical experimental values:**

- Vibrational strong coupling: 50-200 cm⁻¹
- Electronic strong coupling: 100-1000 meV

Cavity-Modified Dynamics
=========================

Modified Energy Landscape
-------------------------

Strong coupling modifies the potential energy surface. Effective potentials include cavity contribution.

Rate Modifications
------------------

Chemical reaction rates can be modified under strong coupling through:

1. **Energy redistribution**
2. **Modified transition states**
3. **Collective effects**

This is an active area of research ("polariton chemistry").

See :doc:`../part1_application/running_simulations` for simulating strong coupling regimes.

