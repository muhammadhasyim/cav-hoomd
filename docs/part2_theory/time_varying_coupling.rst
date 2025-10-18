=======================
Time-Varying Coupling
=======================

This section provides the theoretical foundation for time-dependent cavity-molecule coupling in Cavity HOOMD.

.. contents:: In this Section
   :local:
   :depth: 2

.. note::

   For practical implementation, see :doc:`../part1_application/time_varying_coupling`.

Motivation and Physical Context
================================

Why Time-Dependent Coupling?
-----------------------------

In experiments, cavity coupling can be dynamically controlled through:

1. **Cavity tuning:** Piezo-electric mirrors change cavity length
2. **Pump-probe:** External laser activates coupling
3. **Mode switching:** Transition between cavity modes
4. **Detuning:** Frequency sweeps across resonance

Time-varying coupling in simulations enables studying:

- Non-equilibrium dynamics
- Sudden quench experiments  
- Pump-probe spectroscopy
- Relaxation processes
- Energy redistribution

Experimental Analogues
----------------------

**Typical Experimental Protocols:**

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Experiment
     - Implementation
     - Timescale
   * - Pump-probe
     - Laser pulse activates coupling
     - fs to ps
   * - Cavity tuning
     - Piezo changes cavity length
     - µs to ms
   * - Mode switching
     - Electronic control
     - ns to µs
   * - Detuning sweep
     - Frequency scan
     - ms to s

**Cavity HOOMD simulates the molecular response** to these protocols.

Mathematical Formulation
========================

Time-Dependent Hamiltonian
---------------------------

The Hamiltonian with time-varying coupling:

.. math::

   H(t) = H_{\text{mol}} + \sum_{\lambda} \left[ \frac{1}{2}K_\lambda \tilde{q}_{0,\lambda}^2 + g(t)\tilde{\varepsilon}_{0,\lambda}^{(0)} \tilde{q}_{0,\lambda} D_\lambda + \frac{g(t)^2(\tilde{\varepsilon}_{0,\lambda}^{(0)})^2}{2K_\lambda} D_\lambda^2 \right]

Where:

- :math:`g(t)`: Time-dependent coupling function
- :math:`\tilde{\varepsilon}_{0,\lambda}^{(0)}`: Base coupling strength (constant)
- All other symbols as in :doc:`cavity_forces`

**Key Point:** The coupling strength :math:`g(t)` multiplies both the linear coupling and quadratic self-energy terms.

Modified Equations of Motion
-----------------------------

**Nuclear Motion:**

.. math::

   M_{nj}\ddot{R}_{nj} = F_{nj}^{(0)} - \sum_{\lambda} \left( g(t)\tilde{\varepsilon}_{0,\lambda}^{(0)}\tilde{q}_{0,\lambda} + \frac{g(t)^2(\tilde{\varepsilon}_{0,\lambda}^{(0)})^2}{K_\lambda} D_\lambda \right) \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}

**Photonic Mode Dynamics:**

.. math::

   m_{0,\lambda}\ddot{\tilde{q}}_{0,\lambda} = -K_\lambda \tilde{q}_{0,\lambda} - g(t)\tilde{\varepsilon}_{0,\lambda}^{(0)} D_\lambda

**Time Derivative Effects:**

When :math:`g(t)` changes, an additional "quench force" appears if we consider the full time-dependent Hamiltonian formulation. However, for practical implementations with smooth or instantaneous changes, this effect is negligible or absorbed into the redefinition of equilibrium positions.

Step Function Protocol
======================

Mathematical Definition
-----------------------

The step function is the simplest time-varying protocol:

.. math::

   g(t) = \begin{cases}
   g_{\text{initial}} & \text{if } t < t_{\text{switch}} \\
   g_{\text{final}} & \text{if } t \geq t_{\text{switch}}
   \end{cases}

**Common case:** :math:`g_{\text{initial}} = 0`, :math:`g_{\text{final}} = g_0` (coupling activation)

**Implementation:**

In discrete timesteps:

.. math::

   g(n\Delta t) = \begin{cases}
   g_{\text{initial}} & \text{if } n < n_{\text{switch}} \\
   g_{\text{final}} & \text{if } n \geq n_{\text{switch}}
   \end{cases}

where :math:`n_{\text{switch}} = \lfloor t_{\text{switch}}/\Delta t \rfloor`.

Energy Redistribution at Switch
--------------------------------

**Before switch** (t < t_switch, g=0):

.. math::

   E_{\text{total}}(t^-) = E_{\text{kinetic}}^{\text{mol}} + E_{\text{potential}}

No cavity or coupling energy.

**After switch** (t ≥ t_switch, g=g_0):

.. math::

   E_{\text{total}}(t^+) = E_{\text{kinetic}}^{\text{mol}} + E_{\text{potential}} + E_{\text{cavity}} + E_{\text{coupling}} + E_{\text{self}}

**Energy Conservation:**

.. math::

   E_{\text{total}}(t^+) = E_{\text{total}}(t^-)

Energy is **conserved** during the switch but **redistributes** among degrees of freedom.

**Physical Interpretation:**

1. Molecular kinetic/potential energy decreases
2. Cavity mode gains energy
3. Coupling and self-energy appear
4. Total energy unchanged (in NVE ensemble)

Cavity Particle Displacement
-----------------------------

**In q=0 mode:**

Cavity particles remain at origin; only forces change.

**In finite-q mode:**

Cavity particle position must jump to new equilibrium:

.. math::

   \vec{q}_{\text{eq}}(t^+) = -\frac{g_0}{K} \vec{D}_{\text{total}}(t_{\text{switch}})

**Reason:** The equilibrium position changes from :math:`\vec{q}=0` (g=0) to :math:`\vec{q}=-\frac{g}{K}\vec{D}` (g≠0).

**Implementation:**

At t = t_switch:

1. Calculate total dipole :math:`\vec{D}_{\text{total}}`
2. Compute new equilibrium :math:`\vec{q}_{\text{new}}`
3. Instantly move cavity particles to :math:`\vec{q}_{\text{new}}`
4. Cavity velocity remains unchanged (momentum conservation)

**Energy accounting:**

The position jump conserves total energy:

.. math::

   \Delta E_{\text{cavity}} + \Delta E_{\text{coupling}} + \Delta E_{\text{self}} = 0

Smooth Coupling Protocols
==========================

Exponential Approach
--------------------

**Formula:**

.. math::

   g(t) = g_{\text{final}} + (g_{\text{initial}} - g_{\text{final}}) e^{-(t-t_0)/\tau}

for :math:`t \geq t_0`.

Where:

- :math:`t_0`: Start time
- :math:`\tau`: Characteristic timescale

**Properties:**

- Smooth transition
- Asymptotically approaches :math:`g_{\text{final}}`
- Rate controlled by :math:`\tau`

**Use cases:**

- Gradual cavity activation
- Coupling decay (lossy cavity)
- Smooth detuning

**Characteristic times:**

- :math:`t_{1/2} = \tau \ln(2) \approx 0.69\tau` (half-life)
- :math:`t_{90\%} \approx 2.3\tau` (90% completion)

Sinusoidal Modulation
----------------------

**Formula:**

.. math::

   g(t) = g_0 + A\sin(2\pi f t + \phi)

Where:

- :math:`g_0`: DC offset
- :math:`A`: Modulation amplitude
- :math:`f`: Modulation frequency
- :math:`\phi`: Initial phase

**Physical realization:**

- Oscillating cavity length
- Frequency modulation
- Periodic pump pulses

**Floquet regime:**

When :math:`f \sim \omega_{\text{vib}}`, the system enters the Floquet (driven) regime with:

- Parametric resonances
- Sideband generation
- Modified selection rules

Square Wave (Periodic On-Off)
------------------------------

**Formula:**

.. math::

   g(t) = \begin{cases}
   g_{\text{on}} & \text{if } \mod(t, T) < T \times \text{duty} \\
   g_{\text{off}} & \text{otherwise}
   \end{cases}

Where:

- :math:`T`: Period
- duty: Duty cycle (0-1)

**Example:** duty=0.5 gives 50% on, 50% off.

**Applications:**

- Pulsed excitation
- Stroboscopic measurement
- Time-resolved spectroscopy

Energy Conservation Theory
===========================

Total Energy Functional
-----------------------

**For time-dependent g(t):**

.. math::

   E_{\text{total}}(t) = \sum_{n,j} \frac{1}{2}M_{nj}\dot{R}_{nj}^2 + V(\{R\}) + \sum_{\lambda} \left[ \frac{1}{2}m_\lambda \dot{\tilde{q}}_\lambda^2 + \frac{1}{2}K_\lambda \tilde{q}_\lambda^2 + g(t)\tilde{\varepsilon}\tilde{q}_\lambda D_\lambda + \frac{g(t)^2\tilde{\varepsilon}^2}{2K_\lambda} D_\lambda^2 \right]

**Time derivative:**

.. math::

   \frac{dE_{\text{total}}}{dt} = \frac{\partial E}{\partial t}\bigg|_{\text{explicit}} + \sum_i \frac{\partial E}{\partial q_i}\dot{q}_i + \sum_i \frac{\partial E}{\partial \dot{q}_i}\ddot{q}_i

For Hamiltonian dynamics with time-independent H: :math:`\frac{dE}{dt} = 0`

For time-dependent H: :math:`\frac{dE}{dt} = \frac{\partial H}{\partial t}`

Numerical Energy Conservation
------------------------------

**In practice:**

.. math::

   \frac{dE_{\text{total}}}{dt} = \frac{\partial g(t)}{\partial t} \left[ \tilde{\varepsilon}\tilde{q}_\lambda D_\lambda + \frac{g(t)\tilde{\varepsilon}^2}{K_\lambda} D_\lambda^2 \right]

**For step function:** :math:`\frac{\partial g}{\partial t} = 0` except at t=t_switch where it's a delta function.

**Energy conservation check:**

.. math::

   \Delta E_{\text{rel}} = \frac{|E(t) - E(t_0)|}{|E(t_0)|} < \epsilon_{\text{tol}}

Typically :math:`\epsilon_{\text{tol}} \approx 10^{-4}` for good simulations.

Adiabatic vs Sudden Limits
===========================

Adiabatic Limit
---------------

**Condition:**

.. math::

   \tau_{\text{change}} \gg \frac{2\pi}{\omega_{\text{vib}}}

The system adjusts instantaneously to changing coupling.

**Behavior:**

- System remains in instantaneous ground state
- No excess heating
- Reversible process

**Cavity mode response:**

.. math::

   \tilde{q}_\lambda(t) \approx -\frac{g(t)}{K_\lambda} D_\lambda(t)

Tracks equilibrium instantaneously.

Sudden Quench Limit
--------------------

**Condition:**

.. math::

   \tau_{\text{change}} \ll \frac{2\pi}{\omega_{\text{vib}}}

Parameters change faster than system can respond.

**Behavior:**

- Non-equilibrium excitation
- Energy redistribution
- Transient dynamics
- Heating possible

**Cavity mode response:**

Oscillates around new equilibrium with initial conditions from old equilibrium.

Intermediate Regime
-------------------

**When:**

.. math::

   \tau_{\text{change}} \sim \frac{2\pi}{\omega_{\text{vib}}}

**Behavior:**

- Partial excitation
- Depends sensitively on details
- Rich dynamics
- May exhibit resonances

Applications and Observables
=============================

Non-Equilibrium Spectroscopy
-----------------------------

**Pump-probe signal:**

.. math::

   S(t, \Delta t) = \langle O(t+\Delta t) O(t) \rangle

where O is an observable (dipole, energy, etc.).

**Time-resolved measurements:**

- Transient absorption
- Fluorescence upconversion
- Two-dimensional spectroscopy

Relaxation Dynamics
-------------------

**After sudden quench:**

.. math::

   O(t) - O_{\text{eq}} \propto e^{-t/\tau_{\text{relax}}}

Measure :math:`\tau_{\text{relax}}` to understand:

- Equilibration rates
- Energy dissipation pathways
- Coupling to thermostats

Energy Transfer Rates
---------------------

**From molecules to cavity:**

.. math::

   \frac{dE_{\text{cavity}}}{dt} = -g(t)\tilde{\varepsilon} \dot{\tilde{q}}_\lambda D_\lambda

**From cavity to molecules:**

.. math::

   \frac{dE_{\text{mol}}}{dt} = -\frac{dE_{\text{cavity}}}{dt} - \frac{dE_{\text{coupling}}}{dt} - \frac{dE_{\text{self}}}{dt}

**Net flow:** Depends on :math:`g(t)` protocol and initial conditions.

Practical Considerations
========================

Timestep Requirements
---------------------

**For smooth changes:**

.. math::

   \Delta t \ll \min\left( \frac{1}{\omega_{\text{max}}}, \tau_{\text{change}} \right)

**For step changes:**

.. math::

   \Delta t \ll \frac{1}{\omega_{\text{max}}}

No additional constraint from step itself (delta function in time).

Thermostat Effects
------------------

**With thermostats:**

- Energy is not conserved globally
- Thermostats may absorb/inject energy
- Distinguish physical vs numerical effects

**Recommendation:**

Test energy conservation in NVE first, then add thermostats.

Numerical Stability
-------------------

**Potential issues:**

1. **Large forces at switch:** Use gradual ramp-up
2. **Resonant driving:** Reduce timestep
3. **Cavity particle jump:** Ensure proper implementation
4. **Energy drift:** Monitor and reduce timestep if needed

**Best practices:**

- Start with small g(t) changes
- Verify energy conservation
- Compare different protocols
- Use adaptive timestep when available

Next Sections
=============

Continue to:

- :doc:`energy_conservation` for detailed energy analysis
- :doc:`thermostats` for temperature control with time-varying coupling
- :doc:`strong_coupling` for polariton dynamics

**Related practical guides:**

- :doc:`../part1_application/time_varying_coupling` for implementation
- :doc:`../part1_application/advanced_simulations` for custom protocols

