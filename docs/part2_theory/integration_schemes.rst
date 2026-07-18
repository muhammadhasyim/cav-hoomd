=====================
Integration Schemes
=====================

This section describes the numerical integration methods used in Cavity HOOMD for propagating both cavity and molecular degrees of freedom, including adaptive timestepping for handling disparate timescales and sudden protocol changes.

.. contents:: In this Section
   :local:
   :depth: 2

Overview
========

Cavity-coupled molecular dynamics requires careful numerical integration due to:

1. **Multiple timescales**: Femtosecond cavity oscillations to picosecond molecular motion
2. **Stochastic thermostats**: Both deterministic and random forces
3. **Time-varying coupling**: Sudden or gradual changes in coupling strength
4. **Energy conservation**: Maintaining accuracy while handling stiff equations

Cavity HOOMD employs **split-operator schemes** that separate deterministic propagation from stochastic thermostatting, combined with **adaptive timestepping** to handle protocol changes.

Velocity Verlet for Cavity Modes
=================================

Basic Algorithm
---------------

The cavity mode follows Langevin dynamics with canonical coordinates :math:`(q_\alpha, p_\alpha)` and unit mass. The equations of motion are:

.. math::

   \dot{q}_\alpha = p_\alpha, \quad \dot{p}_\alpha = f_\alpha(t) - \gamma_\mathrm{c} p_\alpha + \xi_\alpha(t)

where:

- :math:`f_\alpha(t) = -\omega_\alpha^2(q_\alpha + \lambda\mu_\alpha/\omega_\mathrm{c})`: Effective force (harmonic + coupling)
- :math:`\gamma_\mathrm{c} = 1/(2\tau_\mathrm{c})`: Friction coefficient
- :math:`\xi_\alpha(t)`: Gaussian white noise with :math:`\langle\xi_\alpha^2\rangle = 6\gamma_\mathrm{c}k_BT/\Delta t`

Split-Operator Scheme
----------------------

**First Half-Step (Deterministic):**

Advance momentum by half timestep:

.. math::

   p_\alpha(t + \Delta t/2) = p_\alpha(t) + \frac{1}{2}f_\alpha(t)\Delta t

**Position Update:**

Update position using half-step momentum:

.. math::

   q_\alpha(t + \Delta t) = q_\alpha(t) + p_\alpha(t + \Delta t/2)\Delta t

**Second Half-Step (Stochastic):**

Modify force with friction and noise:

.. math::

   f_\alpha(t + \Delta t) \to f_\alpha(t + \Delta t) - \gamma_\mathrm{c} p_\alpha(t + \Delta t/2) + \xi_\alpha(t)

Complete momentum update:

.. math::

   p_\alpha(t + \Delta t) = p_\alpha(t + \Delta t/2) + \frac{1}{2}f_\alpha(t + \Delta t)\Delta t

Noise Variance
~~~~~~~~~~~~~~

The random force variance satisfies the fluctuation-dissipation theorem:

.. math::

   \langle\xi_\alpha^2\rangle = \frac{6\gamma_\mathrm{c}k_BT}{\Delta t}

**Implementation:**

Sample :math:`\xi_\alpha \sim \mathcal{N}(0, \sqrt{6\gamma_\mathrm{c}k_BT/\Delta t})` from a normal distribution at each timestep.

Physical Interpretation
-----------------------

**Properties:**

- Maintains cavity at temperature :math:`T`
- Thermalizes kinetic energy at rate :math:`\gamma_\mathrm{c}`
- Represents photon loss (cavity lifetime :math:`\tau_\mathrm{c} = 1/(2\gamma_\mathrm{c})`)
- Typical: :math:`\gamma_\mathrm{c} = 0.5` ps⁻¹ for ultrafast energy exchange

**Time-Reversible Core:**

Without thermostat (:math:`\gamma_\mathrm{c} = 0`, :math:`\xi = 0`), the scheme is time-reversible and symplectic.

Velocity Verlet for Molecules
==============================

Two-Step Propagator
-------------------

Molecular coordinates :math:`(\mathbf{R}_i, \mathbf{P}_i)` with mass :math:`M_i` evolve with Bussi-Parrinello thermostat. Each timestep :math:`\Delta t` consists of:

**First Half-Step:**

Advance momenta by half timestep:

.. math::

   \mathbf{P}_i(t + \Delta t/2) = \mathbf{P}_i(t) + \frac{1}{2}\mathbf{F}_i(t)\Delta t

where :math:`\mathbf{F}_i(t)` is the total force (molecular + cavity coupling).

**First Velocity Rescaling:**

Rescale momenta by stochastic factor :math:`\alpha_\upsilon` computed from current kinetic energy:

.. math::

   K_{t+\Delta t/2} = \sum_i \frac{|\mathbf{P}_i(t + \Delta t/2)|^2}{2M_i}

The scaling factor is:

.. math::

   \alpha_\upsilon = \left[e^{-\Delta t/\tau_\mathrm{b}} + \frac{T}{2K_{t+\Delta t/2}}(1 - e^{-\Delta t/\tau_\mathrm{b}})(R_\Gamma + R_N^2) + 2R_N\sqrt{\frac{T}{2K_{t+\Delta t/2}}e^{-\Delta t/\tau_\mathrm{b}}(1 - e^{-\Delta t/\tau_\mathrm{b}})}\right]^{1/2}

where :math:`R_\Gamma` and :math:`R_N` are independent standard normal random variables.

Apply rescaling:

.. math::

   \mathbf{P}_i(t + \Delta t/2) \to \alpha_\upsilon \mathbf{P}_i(t + \Delta t/2)

**Position Update:**

Update positions:

.. math::

   \mathbf{R}_i(t + \Delta t) = \mathbf{R}_i(t) + \frac{\mathbf{P}_i(t + \Delta t/2)}{M_i}\Delta t

**Force Recalculation:**

Compute forces at new positions :math:`\mathbf{F}_i(t + \Delta t)`.

**Second Half-Step:**

Advance momenta:

.. math::

   \mathbf{P}_i(t + \Delta t) = \mathbf{P}_i(t + \Delta t/2) + \frac{1}{2}\mathbf{F}_i(t + \Delta t)\Delta t

**Second Velocity Rescaling:**

Apply a freshly computed rescaling factor :math:`\alpha'_\upsilon`:

.. math::

   \mathbf{P}_i(t + \Delta t) \to \alpha'_\upsilon \mathbf{P}_i(t + \Delta t)

Bussi-Parrinello Implementation
--------------------------------

**Key Parameters:**

- :math:`\tau_\mathrm{b} = 1.0` ps: Thermostat time constant
- Two rescalings per timestep ensure proper canonical sampling
- Rescaling factors depend on instantaneous kinetic energy

**Advantages:**

- Exactly samples canonical ensemble
- Preserves short-time velocity autocorrelation
- Minimal perturbation to microscopic dynamics
- Ideal for supercooled liquids where dynamics are critical

Force Composition
-----------------

The total force on particle :math:`i` is:

.. math::

   \mathbf{F}_i = \mathbf{F}_i^\mathrm{mol} + \mathbf{F}_i^\mathrm{cavity}

**Molecular Force:**

.. math::

   \mathbf{F}_i^\mathrm{mol} = -\nabla_{\mathbf{R}_i} V

Includes bonds, angles, dihedrals, Lennard-Jones, Coulomb, etc.

**Cavity Force:**

.. math::

   \mathbf{F}_i^\mathrm{cavity} = -\sum_\alpha \left(\epsilon_{0,\alpha} q_{0,\alpha} + \frac{\epsilon_{0,\alpha}^2}{K_\alpha}D_\alpha\right) \frac{\partial d_{n\alpha}}{\partial \mathbf{R}_i}

where :math:`d_{n\alpha}` is the component of molecular dipole moment in direction :math:`\alpha`.

Adaptive Timestepping
======================

Motivation
----------

Simulations with disparate timescales and sudden coupling changes require dynamic timestep adjustment:

1. **Femtosecond Rabi oscillations** vs **nanosecond structural dynamics**
2. **Instantaneous step changes** in coupling constant (sudden quench experiments)
3. **Numerical stability** during protocol switches

Adaptive timestepping dynamically adjusts :math:`\Delta t` to maintain accuracy while maximizing efficiency.

Error Estimation
----------------

**Local Truncation Error:**

Estimate error using explicit Euler step as reference. For particle :math:`i`:

**Velocity-Verlet predictor:**

.. math::

   \mathbf{R}_{i,\mathrm{VV}}(t + \Delta t) = \mathbf{R}_i(t) + \mathbf{v}_i(t)\Delta t + \frac{1}{2}\mathbf{a}_i(t)(\Delta t)^2

where :math:`\mathbf{v}_i = \mathbf{P}_i/M_i` and :math:`\mathbf{a}_i = \mathbf{F}_i/M_i`.

**Explicit Euler:**

.. math::

   \mathbf{R}_{i,\mathrm{E}}(t + \Delta t) = \mathbf{R}_i(t) + \mathbf{v}_i(t)\Delta t

**Per-Particle Error:**

.. math::

   \delta_i(t) = |\mathbf{R}_{i,\mathrm{VV}} - \mathbf{R}_{i,\mathrm{E}}| = \frac{1}{2}|\mathbf{a}_i(t)|(\Delta t)^2

**System-Wide Error:**

Root-mean-square over all :math:`N` particles:

.. math::

   \delta(t) = \sqrt{\frac{1}{N}\sum_{i=1}^N \delta_i(t)^2}

Timestep Adjustment
-------------------

**Update Rule:**

.. math::

   \Delta t' = \Delta t \sqrt{\frac{\varepsilon^*}{\delta(t)}}

where :math:`\varepsilon^*` is the target error tolerance.

**Behavior:**

- :math:`\delta > \varepsilon^*`: Decrease timestep (high forces)
- :math:`\delta < \varepsilon^*`: Increase timestep (weak forces)
- :math:`\delta \approx \varepsilon^*`: Maintain timestep (optimal)

**Practical Limits:**

Impose bounds to prevent extreme changes:

.. math::

   \Delta t_\min \leq \Delta t' \leq \Delta t_\max

Typical: :math:`\Delta t_\min = 10^{-5}` fs, :math:`\Delta t_\max = 2.0` fs.

Protocol-Aware Error Control
=============================

Dynamic Error Tolerance
-----------------------

Sudden coupling changes cause numerical instabilities even with adaptive timestepping. Implement **error tolerance ramping**:

**Reset and Ramp:**

When coupling constant changes by :math:`\Delta\lambda \geq \varepsilon_\lambda`, reset error tolerance:

.. math::

   \varepsilon(t) = \begin{cases}
   \varepsilon^*, & t < t_0 \\
   \varepsilon^*[1 - (1 - f_0)e^{-(t-t_0)/\tau^*}], & t \geq t_0
   \end{cases}

where:

- :math:`t_0`: Time of coupling change
- :math:`f_0`: Initial fraction (typically :math:`10^{-3}`)
- :math:`\tau^*`: Ramping time constant (typically 50 ps)
- :math:`\varepsilon^*`: Steady-state tolerance (typically 5.0)

**Physical Interpretation:**

- Start with very conservative timestep (:math:`\Delta t \sim 10^{-4}` fs) immediately after switch
- Gradually increase timestep as system equilibrates
- Reach normal timestep (:math:`\Delta t \sim 1.5` fs) after ramping period
- Balances stability during transient vs efficiency in stable regime

Implementation Details
----------------------

**Threshold Detection:**

Monitor coupling constant each step:

.. math::

   \text{if } |\lambda(t) - \lambda(t - \Delta t)| \geq \varepsilon_\lambda \text{ then reset } \varepsilon(t)

**Ramping Function:**

Exponential approach to steady state:

.. math::

   \varepsilon(t) = \varepsilon^* - (\varepsilon^* - f_0\varepsilon^*)e^{-(t-t_0)/\tau^*}

**Disabling Ramping:**

For well-equilibrated continuation runs, set :math:`f_0 = 0` (no ramping) and use fixed :math:`\varepsilon^*` throughout.

Practical Parameters
====================

Recommended Settings
--------------------

**Standard Equilibrium Simulations:**

- :math:`\Delta t = 1.0-1.5` fs (fixed or adaptive)
- :math:`\varepsilon^* = 5.0` (adaptive)
- :math:`\gamma_\mathrm{c} = 0.5` ps⁻¹ (cavity)
- :math:`\tau_\mathrm{b} = 1.0` ps (Bussi)

**High-Frequency Vibrations:**

- :math:`\Delta t = 0.5` fs (fixed)
- :math:`\varepsilon^* = 1.0` (adaptive)
- Necessary for :math:`\omega > 3000` cm⁻¹

**Sudden Quench Experiments:**

- Adaptive timestepping: ON
- Error ramping: ON with :math:`f_0 = 10^{-3}`, :math:`\tau^* = 50` ps
- :math:`\varepsilon^* = 5.0`

**Energy Conservation Tests:**

- Fixed timestep: :math:`\Delta t = 0.5-1.0` fs
- No thermostats (:math:`\gamma = 0`)
- Monitor :math:`\Delta E/E < 10^{-4}` over 100 ps

Timestep Selection
------------------

**Nyquist Criterion:**

Resolve fastest oscillation:

.. math::

   \Delta t < \frac{\pi}{\omega_\max}

For :math:`\omega_\mathrm{c} = 2000` cm⁻¹ :math:`\approx` 377 fs period:

.. math::

   \Delta t < \frac{377\text{ fs}}{\pi} \approx 120\text{ fs}

In practice, use :math:`\Delta t \approx T_\mathrm{min}/10` for accuracy.

**Stiffness Considerations:**

Bonded forces are stiffer than non-bonded. Use smaller timesteps when:

- High bond spring constants
- Strong coupling to cavity
- Large dipole moment gradients

Performance Considerations
==========================

Computational Cost
------------------

**Timestep vs Accuracy Tradeoff:**

- Smaller :math:`\Delta t`: Better accuracy, more steps
- Larger :math:`\Delta t`: Faster simulation, risk instability
- Adaptive: Best of both, overhead ~5%

**Adaptive Overhead:**

Error estimation adds minimal cost:

- Compute error: :math:`\mathcal{O}(N)` per step
- Typically 3-5% overhead
- Worthwhile for stability and efficiency

**Optimal Strategy:**

Use adaptive timestepping with ramping for:

- Time-varying coupling protocols
- Non-equilibrium experiments
- Systems with unknown force scales

Use fixed timesteps for:

- Well-characterized equilibrium systems
- Maximum performance on simple systems
- Benchmarking and validation

GPU Optimization
----------------

**Memory Access Patterns:**

- Coalesced reads of positions and forces
- Reduction for error calculation
- Minimal branching in adaptive logic

**Parallel Error Estimation:**

- Compute :math:`\delta_i` for each particle in parallel
- Reduce to :math:`\delta` using shared memory
- Update timestep on host

**Typical Performance:**

- Fixed timestep: ~10⁶ steps/hour (N=500, NVIDIA A100)
- Adaptive timestep: ~9.5×10⁵ steps/hour (5% overhead)

Validation and Testing
======================

Energy Conservation
-------------------

**NVE Test:**

Run without thermostats, monitor total energy:

.. math::

   \Delta E_\mathrm{rel} = \frac{|E(t) - E(0)|}{|E(0)|}

**Acceptance Criteria:**

- :math:`\Delta E_\mathrm{rel} < 10^{-4}` over 100 ps: Excellent
- :math:`\Delta E_\mathrm{rel} < 10^{-3}` over 100 ps: Acceptable
- :math:`\Delta E_\mathrm{rel} > 10^{-2}`: Investigate (timestep too large, numerical issues)

Canonical Ensemble
------------------

**Equipartition Test:**

Check kinetic energy distribution:

.. math::

   \langle E_\mathrm{kin}\rangle = \frac{1}{2}N_f k_B T

where :math:`N_f` is degrees of freedom.

**Statistical Test:**

Over long trajectory, :math:`E_\mathrm{kin}` should follow expected distribution with correct mean and variance.

Time-Reversibility
------------------

**Reversibility Test:**

1. Run forward for time :math:`t`
2. Reverse velocities: :math:`\mathbf{v} \to -\mathbf{v}`
3. Run backward for time :math:`t`
4. Check return to initial state

**Expected:**

- With thermostats: Not reversible (stochastic)
- Without thermostats: Should return to :math:`\mathbf{r}(0)` within numerical precision

Timestep Convergence
--------------------

**Convergence Study:**

Run same simulation with :math:`\Delta t`, :math:`\Delta t/2`, :math:`\Delta t/4`. Check observable :math:`O(\Delta t)` converges as :math:`\mathcal{O}(\Delta t^2)` (second-order scheme).

Best Practices
==============

Initialization
--------------

**Starting Configuration:**

- Generate initial positions on lattice or from equilibrated structure
- Assign velocities from Maxwell-Boltzmann distribution at :math:`T`
- Run short equilibration with strong thermostat

**Cavity Initialization:**

- Sample :math:`q_\alpha`, :math:`p_\alpha` from thermal distribution at :math:`T`
- For harmonic oscillator: :math:`q \sim \mathcal{N}(0, k_BT/\omega^2)`, :math:`p \sim \mathcal{N}(0, k_BT)`

Equilibration Protocol
-----------------------

**Stage 1: Thermalization (0-10 ps):**

- Strong Bussi thermostat (:math:`\tau_\mathrm{b} = 0.1` ps)
- Fixed small timestep (:math:`\Delta t = 0.5` fs)
- Monitor temperature convergence

**Stage 2: Equilibration (10-100 ps):**

- Standard thermostat (:math:`\tau_\mathrm{b} = 1.0` ps)
- Adaptive timestepping ON
- Monitor energy, pressure, temperature

**Stage 3: Production (>100 ps):**

- Begin data collection
- Continue with production settings
- Enable ramping if protocol involves coupling changes

Troubleshooting
---------------

**Simulation Explodes:**

- Decrease :math:`\Delta t` or :math:`\varepsilon^*`
- Check for particle overlaps in initial configuration
- Enable error ramping for coupling changes

**Energy Drift:**

- Reduce timestep
- Check force cutoffs and neighbor lists
- Verify thermostat parameters

**Poor Statistics:**

- Increase number of trajectories
- Longer equilibration before production
- Check convergence of time averages

Next Sections
=============

Continue to:

- :doc:`observables` for computing time averages and correlations
- :doc:`energy_conservation` for energy analysis
- :doc:`../part1_application/running_simulations` for practical setup
- :doc:`../part3_advanced/performance` for optimization strategies

**Related Theory:**

- :doc:`thermostats` for bath implementation
- :doc:`cavity_forces` for force calculations
- :doc:`time_varying_coupling` for time-dependent protocols

