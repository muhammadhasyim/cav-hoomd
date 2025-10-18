===========
Thermostats
===========

This section describes the thermostats used in Cavity HOOMD for temperature control of both molecular and photonic degrees of freedom.

.. contents:: In this Section
   :local:
   :depth: 2

Overview
========

Thermostats maintain system temperature by coupling to an external heat bath. For cavity simulations, we need thermostats for:

1. **Molecular degrees of freedom:** Control molecular temperature
2. **Photonic degrees of freedom:** Represent cavity dissipation

Cavity HOOMD supports three main thermostats:

- **Bussi-Donadio-Parrinello (Bussi):** Stochastic velocity rescaling
- **Langevin:** Friction and random force
- **None:** Microcanonical (NVE) ensemble

Bussi-Donadio-Parrinello Thermostat
====================================

Theory
------

The Bussi thermostat [Bussi2007]_ generates correct canonical ensemble statistics through stochastic velocity rescaling.

**Algorithm:**

At each step, rescale velocities by factor :math:`\alpha`:

.. math::

   \vec{v}_i \to \alpha \vec{v}_i

where :math:`\alpha` is chosen to drive kinetic energy toward target value while preserving canonical distribution.

**Kinetic energy evolution:**

.. math::

   K(t+\Delta t) = K(t) + (K_{\text{target}} - K(t))\left(1 - e^{-\Delta t/\tau}\right) + 2\sqrt{\frac{K(t)K_{\text{target}}}{N_f}} \left(1 - e^{-\Delta t/2\tau}\right) R

Where:

- :math:`K(t)`: Current kinetic energy
- :math:`K_{\text{target}} = \frac{1}{2}N_f k_B T`: Target kinetic energy
- :math:`N_f`: Number of degrees of freedom
- :math:`\tau`: Coupling time constant
- :math:`R`: Gaussian random variable

**Properties:**

- ✓ Exact canonical ensemble
- ✓ Fast equilibration
- ✓ Minimal perturbation to dynamics
- ✓ Works for small systems

Parameters
----------

**Coupling time τ:**

.. math::

   \tau = \frac{\text{relaxation time}}{k_B T/E_{\text{typ}}}

**Typical values:**

- τ = 0.1 ps: Very strong coupling (fast equilibration)
- τ = 1.0 ps: Standard (recommended)
- τ = 10 ps: Weak coupling (near NVE)

**Choosing τ:**

- **Fast equilibration:** τ ~ 0.1-1.0 ps
- **Dynamics studies:** τ ~ 1-10 ps
- **Near-NVE:** τ > 10 ps

Usage in Cavity HOOMD
----------------------

**For molecules:**

.. code-block:: python

   from hoomd.bussi_reservoir import BussiReservoir
   
   bussi = BussiReservoir(
       kT=0.01,  # Temperature in reduced units (100 K)
       tau=1.0   # Coupling time in ps
   )
   
   molecular_filter = hoomd.filter.Type(['O', 'N'])
   integrator.methods.append(
       hoomd.md.methods.ConstantVolume(
           filter=molecular_filter,
           thermostat=bussi
       )
   )

**For cavity:**

.. code-block:: python

   cavity_filter = hoomd.filter.Type(['Cavity'])
   integrator.methods.append(
       hoomd.md.methods.ConstantVolume(
           filter=cavity_filter,
           thermostat=bussi
       )
   )

Langevin Thermostat
===================

Theory
------

The Langevin equation adds friction and random force:

.. math::

   m\ddot{\vec{r}} = \vec{F} - \gamma m\dot{\vec{r}} + \vec{F}_{\text{random}}

Where:

- :math:`\vec{F}`: Deterministic force
- :math:`\gamma`: Friction coefficient
- :math:`\vec{F}_{\text{random}}`: Random force satisfying fluctuation-dissipation theorem

**Fluctuation-dissipation:**

.. math::

   \langle F_{\text{random},i}(t) F_{\text{random},j}(t') \rangle = 2\gamma m k_B T \delta_{ij} \delta(t-t')

**Integration schemes:**

Various algorithms exist (BBK, BAOAB, etc.). HOOMD uses optimized schemes.

Parameters
----------

**Friction coefficient γ:**

.. math::

   \gamma = \frac{1}{\tau}

where τ is the damping timescale.

**Typical values:**

- γ = 0.1 ps⁻¹: Weak damping (near NVE)
- γ = 1.0 ps⁻¹: Moderate damping
- γ = 10 ps⁻¹: Strong damping (overdamped)

**Physical interpretation:**

- τ = 1/γ is the momentum relaxation time
- After time τ, velocity correlation decays to e⁻¹ ≈ 37%

Usage in Cavity HOOMD
----------------------

**For molecules:**

.. code-block:: python

   molecular_filter = hoomd.filter.Type(['O', 'N'])
   langevin = hoomd.md.methods.Langevin(
       filter=molecular_filter,
       kT=0.01,         # Temperature
       default_gamma=1.0  # Friction coefficient
   )
   integrator.methods.append(langevin)

**For cavity (represents dissipation):**

.. code-block:: python

   cavity_filter = hoomd.filter.Type(['Cavity'])
   cavity_langevin = hoomd.md.methods.Langevin(
       filter=cavity_filter,
       kT=0.01,
       default_gamma=0.1  # Lower γ for cavity
   )
   integrator.methods.append(cavity_langevin)

**Physical meaning for cavity:**

γ_cavity represents photon loss rate (cavity linewidth κ).

Microcanonical (NVE) Ensemble
==============================

Theory
------

Without thermostat, total energy is conserved:

.. math::

   E_{\text{total}} = \text{constant}

**Equations of motion:**

Standard Hamiltonian dynamics with energy-conserving integrator (Velocity Verlet).

**Properties:**

- ✓ Exact energy conservation (within numerical precision)
- ✓ Deterministic dynamics
- ✓ No artificial damping
- ✗ Temperature may drift if not initialized properly
- ✗ No temperature control

Usage in Cavity HOOMD
----------------------

**Setup:**

Don't add any thermostat. Just use ConstantVolume with no thermostat parameter:

.. code-block:: python

   all_particles = hoomd.filter.All()
   integrator.methods.append(
       hoomd.md.methods.ConstantVolume(filter=all_particles)
   )

**When to use:**

- Energy conservation tests
- Short-time dynamics (< 10 ps)
- Microcanonical ensemble studies
- Validating implementation

Hybrid Thermostat Strategies
=============================

Molecular + Cavity Combinations
--------------------------------

Different thermostat combinations serve different purposes:

**1. Bussi (molecules) + Langevin (cavity)** [Recommended]

.. code-block:: bash

   --molecular-bath bussi --cavity-bath langevin

**Rationale:**

- Molecules: Canonical ensemble, fast equilibration
- Cavity: Physical dissipation (represents photon loss)

**Use for:** Most production simulations

**2. Bussi (molecules) + Bussi (cavity)**

.. code-block:: bash

   --molecular-bath bussi --cavity-bath bussi

**Rationale:**

- Both subsystems in canonical ensemble
- Useful for equilibrium thermodynamics

**Use for:** Studying equilibrium properties

**3. None (molecules) + None (cavity)**

.. code-block:: bash

   --molecular-bath none --cavity-bath none

**Rationale:**

- Pure NVE dynamics
- No energy exchange with bath
- Total energy conserved

**Use for:** Energy conservation testing, short-time dynamics

**4. Langevin (molecules) + Langevin (cavity)**

.. code-block:: bash

   --molecular-bath langevin --cavity-bath langevin

**Rationale:**

- Both damped and stochastic
- Strongly dissipative

**Use for:** High-friction systems, overdamped dynamics

Energy Flow Considerations
---------------------------

**With thermostats:**

.. math::

   \frac{dE_{\text{total}}}{dt} = P_{\text{thermostat}}

where :math:`P_{\text{thermostat}}` is power injected/removed by thermostat.

**Energy balance:**

- Thermostats maintain temperature by absorbing/injecting energy
- Total energy NOT conserved (as expected for canonical ensemble)
- Temperature IS controlled

**Monitoring:**

Track both energy and temperature to ensure proper thermostat function.

Temperature Measurement
=======================

Kinetic Temperature Definition
-------------------------------

**Standard definition:**

.. math::

   T_{\text{kinetic}} = \frac{2\langle E_{\text{kinetic}} \rangle}{N_f k_B}

Where :math:`N_f` is number of degrees of freedom.

**For molecules:**

.. math::

   T_{\text{mol}} = \frac{2}{3Nk_B} \sum_{i=1}^{N} \frac{1}{2}m_i v_i^2

**For cavity:**

.. math::

   T_{\text{cavity}} = \frac{1}{N_{\lambda}k_B} \sum_{\lambda} \frac{1}{2}m_{0,\lambda}v_{0,\lambda}^2

where :math:`N_\lambda` is number of polarizations.

Equipartition Check
-------------------

**At equilibrium:**

Each quadratic degree of freedom should have energy :math:`\frac{1}{2}k_B T`.

**Check:**

.. math::

   \langle E_{\text{kinetic}} \rangle \stackrel{?}{=} \frac{1}{2}N_f k_B T

Deviation > 5% indicates:

- Non-equilibrium state
- Thermostat malfunction
- Insufficient equilibration time

Equilibration Protocol
=======================

Standard Procedure
------------------

**1. Initial configuration (t=0):**

- Generate or load initial positions
- Assign velocities from Maxwell-Boltzmann distribution at target T

**2. Equilibration phase (0 < t < t_eq):**

- Use strong thermostat (Bussi with τ=0.1-1.0 ps)
- Monitor temperature, energy, pressure
- Run until properties stabilize (typically 10-100 ps)

**3. Production phase (t > t_eq):**

- Switch to production thermostat settings
- Begin data collection
- Run for desired simulation time

Equilibration Checks
---------------------

**Temperature stability:**

.. math::

   |\langle T \rangle - T_{\text{target}}| < 0.01 T_{\text{target}}

**Energy fluctuations:**

.. math::

   \frac{\sigma_E}{\langle E \rangle} \approx \sqrt{\frac{2k_B T^2}{C_V}}

Should match expected fluctuations for canonical ensemble.

**Autocorrelation:**

Check that observables decorrelate:

.. math::

   C(t) = \frac{\langle O(t')O(t'+t) \rangle - \langle O \rangle^2}{\langle O^2 \rangle - \langle O \rangle^2} \to 0

as :math:`t \to \infty`.

Special Considerations for Cavity Systems
==========================================

Cavity Lifetime
----------------

**Physical cavity lifetime:**

.. math::

   \tau_{\text{cavity}} = \frac{1}{\kappa}

where κ is photon loss rate.

**Langevin thermostat maps:**

.. math::

   \gamma_{\text{cavity}} \leftrightarrow \kappa

**Typical experimental values:**

- High-Q cavity: κ ~ 0.01 ps⁻¹ (τ ~ 100 ps)
- Medium-Q: κ ~ 0.1 ps⁻¹ (τ ~ 10 ps)
- Low-Q: κ ~ 1.0 ps⁻¹ (τ ~ 1 ps)

Thermalization Timescales
--------------------------

**Molecular thermalization:**

.. math::

   \tau_{\text{mol}} \sim \frac{1}{\gamma_{\text{mol}}}

**Cavity thermalization:**

.. math::

   \tau_{\text{cavity}} \sim \frac{1}{\gamma_{\text{cavity}}}

**Energy exchange time:**

.. math::

   \tau_{\text{exchange}} \sim \frac{1}{g \sqrt{N}}

**Hierarchy:**

Typically :math:`\tau_{\text{exchange}} < \tau_{\text{mol}} < \tau_{\text{cavity}}`

Next Sections
=============

Continue to:

- :doc:`energy_conservation` for energy analysis with thermostats
- :doc:`strong_coupling` for polariton dynamics
- :doc:`../part1_application/running_simulations` for practical thermostat selection

.. [Bussi2007] G. Bussi, D. Donadio, and M. Parrinello, J. Chem. Phys. 126, 014101 (2007)

