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

- Exact canonical ensemble
- Fast equilibration
- Minimal perturbation to dynamics
- Works for small systems

Detailed Derivation
~~~~~~~~~~~~~~~~~~~

The Bussi-Parrinello thermostat is derived by analyzing standard Langevin dynamics where each particle experiences independent friction and noise:

.. math::

   \dot{\mathbf{R}}_i = \frac{\mathbf{P}_i}{M_i}, \quad \dot{\mathbf{P}}_i = -\nabla_{\mathbf{R}_i} V - \gamma_\mathrm{b} \mathbf{P}_i + \sqrt{2M_i\gamma_\mathrm{b}k_BT}\,\boldsymbol{\eta}_i(t)

where :math:`\boldsymbol{\eta}_i(t)` are independent white noise processes with :math:`\langle\boldsymbol{\eta}_i(t)\boldsymbol{\eta}_j(t')\rangle = \boldsymbol{\delta}_{ij}\delta(t-t')`.

**Energy Evolution:**

Calculating the time derivative of total energy :math:`H(t)` using the Ito chain rule and summing over all noise terms yields:

.. math::

   \dot{H}(t) = -\frac{K(t) - \bar{K}}{\tau_\mathrm{b}} + 2\sqrt{\frac{K(t)\bar{K}}{3N\tau_\mathrm{b}}}\,\eta(t)

where :math:`\bar{K} = \frac{3}{2}Nk_BT`, :math:`\tau_\mathrm{b} = (2\gamma_\mathrm{b})^{-1}`, and :math:`\eta(t)` is a single global white noise emerging from the sum of :math:`3N` independent noise terms.

**Key Observation:**

Only the total kinetic energy :math:`K` appears in the energy evolution. We can design an alternative stochastic force :math:`\mathbf{G}_i` that depends globally on :math:`K` and reproduces the same rate of energy exchange while minimizing perturbations to individual particle trajectories.

**Minimizing Disturbance:**

Quantify the disturbance on kinetic energy as:

.. math::

   D = \sum_{i=1}^N \frac{|\mathbf{G}_i|^2}{2M_i}

and minimize it subject to the constraint that :math:`\dot{H}(t)` matches the Langevin result. Geometrically, the constraint forces :math:`\mathbf{G}_i` to be parallel to each momentum vector:

.. math::

   \mathbf{G}_i = \mathbf{P}_i \left[-\gamma_K + \sqrt{2D_K}\,\eta(t)\right]

where the scalar functions :math:`\gamma_K` and :math:`D_K` depend on the total kinetic energy :math:`K` and are determined by matching the energy evolution:

.. math::

   \dot{H} = -2\gamma_K K + 2D_K K + 2\sqrt{2D_K}\,K\,\eta(t)

Enforcing equality with the Langevin energy equation gives:

.. math::

   \gamma_K = \frac{1}{2\tau_\mathrm{b}}\left(1 - \left(1 - \frac{1}{3N}\right)\frac{\bar{K}}{K}\right)

.. math::

   D_K = \frac{\bar{K}}{3N\tau_\mathrm{b}K}

**Equations of Motion:**

.. math::

   \dot{\mathbf{R}}_i = \frac{\mathbf{P}_i}{M_i}, \quad \dot{\mathbf{P}}_i = -\nabla_{\mathbf{R}_i} V - \gamma_K \mathbf{P}_i + \sqrt{2D_K}\,\mathbf{P}_i\eta(t)

**Physical Interpretation:**

This thermostat samples the canonical ensemble exactly while preserving short-time velocity autocorrelation functions, which are critical for maintaining correct microscopic dynamics in supercooled liquids. For a single degree of freedom (:math:`N=1`), it reduces to standard Langevin dynamics, whereas in many-particle systems, the global coupling reduces dynamical artifacts from bath coupling.

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

Caldeira-Leggett Cavity Bath Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Even in high-quality optical cavities, photons are lost due to radiative losses through imperfect mirrors. We model this dissipation through coupling with a thermal electromagnetic environment outside the cavity.

**Bath Hamiltonian:**

Following the Caldeira-Leggett approach, couple the cavity mode coordinate :math:`q_\mathrm{c}` (frequency :math:`\omega_\mathrm{c}`, unit mass) linearly to a continuum of bath oscillators:

.. math::

   H_\mathrm{b} = \sum_b \left[ \frac{p_b^2}{2} + \frac{\omega_b^2}{2} \left( x_b - \frac{c_b}{\omega_b^2} q_\alpha \right)^2 \right]

where :math:`(x_b, p_b)` are bath oscillator coordinates with frequencies :math:`\omega_b` and coupling constants :math:`c_b`.

**Generalized Langevin Equation:**

Integrating out the bath degrees of freedom yields:

.. math::

   \ddot{q}_\alpha(t) + \omega_\alpha^2 q_\alpha(t) + \int_0^t dt'\,G(t-t')\,\dot{q}_\alpha(t') = \eta_\alpha(t)

where the memory kernel and noise correlation are related through the spectral density:

.. math::

   J(\omega) = \frac{\pi}{2} \sum_b \frac{c_b^2}{\omega_b} \delta(\omega - \omega_b)

with the friction kernel and noise correlations:

.. math::

   G(t) = \frac{2}{\pi} \int_0^\infty d\omega\,\frac{J(\omega)}{\omega}\cos(\omega t), \quad \langle\eta_\mathrm{c}(t)\eta_\mathrm{c}(t')\rangle = k_BT\,G(t-t')

**Markovian Limit:**

In the Ohmic, Markovian limit where :math:`J(\omega) \approx \gamma_\mathrm{c}\omega` with a high-frequency cutoff, the memory kernel becomes instantaneous, :math:`G(t) \approx 2\gamma_\mathrm{c}\delta(t)`, yielding standard Langevin dynamics.

**Cavity Lifetime:**

The lifetime of the cavity :math:`\tau_\mathrm{c}` is related to the damping rate through :math:`\gamma_\mathrm{c} = 1/(2\tau_\mathrm{c})`. When coupling is turned on, :math:`q_\mathrm{c} \to q_\mathrm{c} + \lambda\mu_\alpha/\omega_\mathrm{c}`, the damping and noise statistics remain unchanged while the deterministic force acquires the coupling term.

**Physical Interpretation:**

This model captures the essential physics of cavity dissipation: photons leak out through mirrors, and thermal photons leak in from the environment. The balance between these processes maintains the cavity at temperature :math:`T`.

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

