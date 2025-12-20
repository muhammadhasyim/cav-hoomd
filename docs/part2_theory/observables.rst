===========
Observables
===========

This section describes the theoretical framework for computing observables in cavity-coupled molecular dynamics, including ensemble averages, time-correlation functions, and material time for aging experiments.

.. contents:: In this Section
   :local:
   :depth: 2

Overview
========

In molecular dynamics simulations, we compute ensemble averages and time-correlation functions to characterize system behavior. For non-equilibrium cavity-coupled systems, special care is needed to properly average over stochastic trajectories and handle time-dependent distributions.

**Key Concepts:**

- **Phase space**: :math:`\Gamma = (\mathbf{R}, \mathbf{P}, \mathbf{q}, \mathbf{p})` for positions, momenta, cavity coordinates, and cavity momenta
- **Distribution**: :math:`\rho^{(\lambda)}(\Gamma, t)` evolves according to Fokker-Planck equation
- **Ensemble average**: Integral over phase space weighted by :math:`\rho`
- **Time-correlation**: Relates observable at two different times

Ensemble Averages
==================

Schrodinger vs Heisenberg Picture
-----------------------------------

**Schrodinger Picture:**

Time dependence is in the phase-space distribution. For an observable :math:`A(\Gamma)`, the ensemble average at time :math:`t` is:

.. math::

   \langle A(t) \rangle_\lambda = \int d\Gamma\,A(\Gamma)\,\rho^{(\lambda)}(\Gamma, t)

**Heisenberg Picture:**

Time dependence is in the observable. Given initial condition :math:`\Gamma`, evolve the observable averaged over all realizations of noise:

.. math::

   A_t(\Gamma; \lambda) = \mathbb{E}\left[ A(\Gamma(t; \lambda)) \mid \Gamma(0) = \Gamma \right]

The ensemble average becomes:

.. math::

   \langle A(t) \rangle_\lambda = \int d\Gamma\,A_t(\Gamma; \lambda)\,\rho_\mathrm{eq}^{(\lambda=0)}(\Gamma)

Both pictures are equivalent but offer different computational advantages.

Fokker-Planck Equation
------------------------

The phase-space distribution :math:`\rho^{(\lambda)}(\Gamma, t)` evolves according to:

.. math::

   \frac{\partial \rho^{(\lambda)}}{\partial t} = \mathcal{L}_\lambda \rho^{(\lambda)}

where the Liouvillian operator is:

.. math::

   \mathcal{L}_\lambda = \mathcal{L}_0 + \mathcal{L}_\mathrm{int}(\lambda) + \mathcal{L}_\mathrm{b}

**Components:**

1. **Bare Liouvillian** :math:`\mathcal{L}_0`: Hamiltonian dynamics at zero coupling

.. math::

   \mathcal{L}_0 \rho = -\sum_i \frac{\mathbf{P}_i}{M_i} \cdot \nabla_{\mathbf{R}_i} \rho + \sum_i \nabla_{\mathbf{R}_i} V \cdot \nabla_{\mathbf{P}_i} \rho - \sum_\alpha p_\alpha \frac{\partial \rho}{\partial q_\alpha} + \sum_\alpha \omega_\mathrm{c}^2 q_\alpha \frac{\partial \rho}{\partial p_\alpha}

2. **Interaction Liouvillian** :math:`\mathcal{L}_\mathrm{int}(\lambda)`: Light-matter coupling

.. math::

   \mathcal{L}_\mathrm{int}(\lambda)\rho = \sum_i \sum_\alpha \lambda z_i e \left(\omega_\mathrm{c} q_\alpha + \lambda \mu_\alpha\right) \hat{\mathbf{e}}_\alpha \cdot \nabla_{\mathbf{P}_i} \rho - \sum_\alpha \lambda \omega_\mathrm{c} \mu_\alpha \frac{\partial \rho}{\partial p_\alpha}

3. **Bath Liouvillian** :math:`\mathcal{L}_\mathrm{b}`: Noise and dissipation from molecular and cavity baths

**Initial Condition:**

For nonthermal aging protocols starting at :math:`t = t_0`:

.. math::

   \rho^{(\lambda)}(\Gamma, t_0) = \rho_\mathrm{eq}^{(\lambda=0)}(\Gamma)

where the equilibrium distribution at zero coupling is:

.. math::

   \rho_\mathrm{eq}^{(\lambda)}(\Gamma) = \frac{1}{Z(T, \lambda)} e^{-\beta H(\Gamma; \lambda)}

Practical Computation
----------------------

**Sample Mean Over Trajectories:**

In simulations, we approximate ensemble averages by averaging over :math:`N_\mathrm{T}` independent trajectories:

.. math::

   \langle A(t) \rangle_\lambda \approx \frac{1}{N_\mathrm{T}} \sum_{n=1}^{N_\mathrm{T}} A^{(n)}(t; \lambda)

where :math:`A^{(n)}(t; \lambda) = A(\Gamma^{(n)}(t; \lambda))` is the observable evaluated on the :math:`n`-th trajectory.

**Equilibrium Averages:**

At equilibrium, we can use both ensemble and time averaging (ergodic hypothesis):

.. math::

   \langle A \rangle_\lambda = \lim_{T \to \infty} \frac{1}{T} \int_0^T dt\,A(\Gamma(t; \lambda))

In practice:

.. math::

   \langle A \rangle_\lambda \approx \frac{1}{MN_\mathrm{T}} \sum_{n=1}^{N_\mathrm{T}} \sum_{p=0}^{M-1} A(\Gamma^{(n)}(t_p; \lambda))

where we subsample at times :math:`t_p = p\Delta t`.

Static Structure Factor Example
--------------------------------

For density-mode fluctuations at wavevector :math:`\mathbf{k}`:

.. math::

   \rho_{\mathbf{k}}^{(n)}(t; \lambda) = \sum_{i=1}^N e^{i\mathbf{k} \cdot \mathbf{R}_i^{(n)}(t; \lambda)}

The static structure factor is:

.. math::

   S_{\mathbf{k}}^{(\lambda)}(t_\mathrm{w}) = \left\langle \rho_{\mathbf{k}}(t_\mathrm{w}) \rho_{\mathbf{k}}^*(t_\mathrm{w}) \right\rangle_\lambda \approx \frac{1}{N_\mathrm{T}} \sum_{n=1}^{N_\mathrm{T}} \rho_{\mathbf{k}}^{(n)}(t_\mathrm{w}; \lambda) \rho_{\mathbf{k}}^{*(n)}(t_\mathrm{w}; \lambda)

At equilibrium, this becomes time-independent:

.. math::

   S_{\mathbf{k}}^{(\lambda)} \approx \frac{1}{N_\mathrm{T}M} \sum_{n=1}^{N_\mathrm{T}} \sum_{p=0}^{M-1} \rho_{\mathbf{k}}^{(n)}(t_p; \lambda) \rho_{\mathbf{k}}^{*(n)}(t_p; \lambda)

Time-Correlation Functions
===========================

General Non-Equilibrium Case
-----------------------------

For an observable :math:`A(t)`, the time-correlation function in a non-equilibrium aging experiment is:

.. math::

   C_{AA}^{(\lambda)}(t; t_\mathrm{w}) = \langle A(t + t_\mathrm{w}) A(t_\mathrm{w}) \rangle_\lambda

where :math:`t_\mathrm{w}` is the waiting time after initiating the protocol and :math:`t` is the correlation time.

**Interpretation:**

- Depends on both :math:`t` (time difference) and :math:`t_\mathrm{w}` (absolute time)
- Reflects non-equilibrium, aging dynamics
- Converges to equilibrium correlation as :math:`t_\mathrm{w} \to \infty`

**Ensemble Average Form:**

In terms of the conditionally averaged observable and phase-space distribution:

.. math::

   C_{AA}^{(\lambda)}(t; t_\mathrm{w}) = \int d\Gamma\,A_{t+t_\mathrm{w}}(\Gamma; \lambda)\,A(\Gamma)\,\rho^{(\lambda)}(\Gamma, t_\mathrm{w})

**Trajectory Sampling:**

We approximate by generating :math:`N_\mathrm{T}` trajectories from equilibrium at zero coupling, simulating directly to :math:`t + t_\mathrm{w}`:

.. math::

   C_{AA}^{(\lambda)}(t; t_\mathrm{w}) \approx \frac{1}{N_\mathrm{T}} \sum_{n=1}^{N_\mathrm{T}} A^{(n)}(t + t_\mathrm{w}; \lambda)\,A^{(n)}(t_\mathrm{w}; \lambda)

Intermediate Scattering Function
----------------------------------

The time-correlation function of density modes is the **intermediate scattering function (ISF)**:

.. math::

   F_{\mathbf{k}}^{(\lambda)}(t; t_\mathrm{w}) = \left\langle \rho_{\mathbf{k}}(t + t_\mathrm{w}; \lambda) \rho_{\mathbf{k}}^*(t_\mathrm{w}; \lambda) \right\rangle_\lambda

**Practical Implementation:**

Average over :math:`N_{\mathbf{k}}` wavevectors of equal magnitude distributed on a sphere:

.. math::

   F_{k}^{(\lambda)}(t; t_\mathrm{w}) \approx \frac{1}{N_\mathrm{T}N_{\mathbf{k}}} \sum_{n=1}^{N_\mathrm{T}} \sum_{s=1}^{N_{\mathbf{k}}} \rho_{\mathbf{k}_s}^{(n)}(t + t_\mathrm{w}; \lambda)\,\rho_{\mathbf{k}_s}^{*(n)}(t_\mathrm{w}; \lambda)

**Normalized ISF:**

.. math::

   \phi_k(t; t_\mathrm{w}) = \frac{F_k(t; t_\mathrm{w})}{S_k(t_\mathrm{w})}

**Structural Relaxation Time:**

Define :math:`\tau_\mathrm{s}(t_\mathrm{w}, \lambda)` as the time when ISF decays to threshold :math:`a` (typically 0.1):

.. math::

   \phi_k^{(\lambda)}(\tau_\mathrm{s}; t_\mathrm{w}) = a

Located by linear interpolation between adjacent data points.

Equilibrium Case
----------------

In equilibrium (large :math:`t_\mathrm{w}`), time-correlation functions obey **time-translational invariance**:

.. math::

   C_{AA}^{(\lambda)}(t; t_\mathrm{w}) = C_{AA,\mathrm{eq}}^{(\lambda)}(t)

depending only on time difference :math:`t`.

**Proof:**

At equilibrium, :math:`\rho^{(\lambda)}(\Gamma, t_\mathrm{w}) = \rho_\mathrm{eq}^{(\lambda)}(\Gamma)` which satisfies :math:`\mathcal{L}_\lambda \rho_\mathrm{eq}^{(\lambda)} = 0`. Since the dynamics are time-translationally invariant, :math:`A_{t+t_\mathrm{w}}(\Gamma) = A_t(\Gamma)` with initial condition at :math:`t=0`.

**Ergodic Hypothesis:**

At equilibrium, ensemble average equals long-time average:

.. math::

   C_{AA,\mathrm{eq}}^{(\lambda)}(t) = \lim_{T \to \infty} \frac{1}{T} \int_0^T dt'\,A(\Gamma(t + t'; \lambda))\,A(\Gamma(t'; \lambda))

**Discrete Estimator:**

Sample at times :math:`t_p = p\Delta t`, compute autocorrelation at lag :math:`\tau_k = k\Delta t`:

.. math::

   \bar{C}_{AA}^{(\lambda)}(\tau_k) \approx \frac{1}{N_\mathrm{T}(M-k)} \sum_{n=1}^{N_\mathrm{T}} \sum_{p=0}^{M-k-1} A^{(n)}(t_p; \lambda)\,A^{(n)}(t_{p+k}; \lambda)

This is the standard "sliding time-origin" implementation.

Dipole Autocorrelation and IR Spectrum
---------------------------------------

**Centered DACF:**

For dipole moment, subtract mean to remove bias:

.. math::

   \delta\mu_\alpha(t) = \mu_\alpha(t) - \bar{\mu}_\alpha

.. math::

   \bar{C}_{\mu\mu}^{(\lambda)}(\tau_k) \approx \frac{1}{N_\mathrm{T}(M-k)} \sum_{n=1}^{N_\mathrm{T}} \sum_{p=0}^{M-k-1} \sum_{\alpha=1,2} \delta\mu_\alpha(t_p; \lambda)\,\delta\mu_\alpha(t_{p+k}; \lambda)

**IR Absorption Spectrum:**

From the fluctuation-dissipation theorem:

.. math::

   n(\omega)\alpha(\omega) = \frac{\beta\omega^2}{4\epsilon_0 Vc} \int_{-\infty}^{\infty} dt\,e^{i\omega t}\,\bar{C}_{\mu\mu}^{(\lambda)}(t)

Using time-reversibility, :math:`\bar{C}_{\mu\mu}(t) = \bar{C}_{\mu\mu}(-t)`:

.. math::

   n(\omega)\alpha(\omega) = \frac{\beta\omega^2}{4\epsilon_0 Vc} \cdot 2\int_0^{\tau_{\max}} dt\,\cos(\omega t)\,\bar{C}_{\mu\mu}^{(\lambda)}(t)

**Discrete Cosine Transform:**

Sample :math:`C_n = C_\mu(n\Delta t)` uniformly with :math:`N` points. Define frequencies:

.. math::

   \omega_k = \frac{\pi k}{(N-1)\Delta t}

Compute DCT:

.. math::

   X_k = (-1)^k C_{N-1} + 2\sum_{n=1}^{N-2} C_n \cos\left(\frac{\pi nk}{N-1}\right)

Trapezoidal rule yields:

.. math::

   n(\omega_k)\alpha(\omega_k) \approx \frac{\beta\omega_k^2}{4\epsilon_0 Vc} \cdot \frac{\Delta t}{2} X_k

Material Time
=============

Material-Time Translational Invariance
---------------------------------------

In aging systems, time-correlation functions across different protocols can be reparameterized using **material time** :math:`h_\lambda(t)`:

**MTTI Hypothesis:**

All correlation functions collapse to a universal curve when plotted against material-time differences:

.. math::

   C_{AA}^{(\lambda)}(t; t_\mathrm{w}) = \Phi\left(h_\lambda(t + t_\mathrm{w}) - h_\lambda(t_\mathrm{w})\right)

**Physical Interpretation:**

Material time measures the "internal age" of the system. Equal material-time differences correspond to equal structural evolution, regardless of absolute time or coupling strength.

Reconstruction from Relaxation Times
--------------------------------------

**Structural Relaxation Constraint:**

For each waiting time :math:`t_\mathrm{w}^{(r)}`, measure :math:`\tau_\mathrm{s}^{(r)}` where :math:`\phi_k(\tau_\mathrm{s}^{(r)}; t_\mathrm{w}^{(r)}) = 0.1`. MTTI implies:

.. math::

   h_\lambda(t_\mathrm{w}^{(r)} + \tau_\mathrm{s}^{(r)}) - h_\lambda(t_\mathrm{w}^{(r)}) = 1, \quad r = 1, \ldots, N_\mathrm{w}

**Least-Squares Formulation:**

Discretize time on a grid :math:`t_m = m\Delta t` for :math:`m = 0, \ldots, M_\mathrm{grid}-1`. Parameterize:

.. math::

   \tilde{h}_\lambda(t) = \sum_{m=0}^{M_\mathrm{grid}-1} h_\lambda^m \theta_m(t)

where :math:`\theta_m(t)` are piecewise linear basis functions:

.. math::

   \theta_m(t) = \begin{cases}
   \frac{t - t_{m-1}}{t_m - t_{m-1}}, & t \in [t_{m-1}, t_m] \\
   \frac{t_{m+1} - t}{t_{m+1} - t_m}, & t \in [t_m, t_{m+1}] \\
   0, & \text{otherwise}
   \end{cases}

**Constraint Matrix:**

Define :math:`N_\mathrm{w} \times M_\mathrm{grid}` matrix:

.. math::

   A_{rm} = \theta_m(t_\mathrm{w}^{(r)} + \tau_\mathrm{s}^{(r)}) - \theta_m(t_\mathrm{w}^{(r)})

so the constraints become linear system :math:`\mathbf{A}\mathbf{h}_\lambda = \mathbf{b}` where :math:`\mathbf{b} = [1, \ldots, 1]^\mathsf{T}`.

Regularized Fitting
--------------------

**Smoothness Penalty:**

To avoid overfitting, add curvature penalty. Define second-derivative matrix :math:`\mathbf{D}` (finite differences):

.. math::

   D_{rm} = \begin{cases}
   1, & m = r \\
   -2, & m = r+1 \\
   1, & m = r+2 \\
   0, & \text{otherwise}
   \end{cases}

**Regularized Objective:**

.. math::

   \mathbf{h}_\lambda^* = \arg\min_{\mathbf{h}_\lambda} \left\{ \|\mathbf{A}\mathbf{h}_\lambda - \mathbf{b}\|^2 + \alpha\|\mathbf{D}\mathbf{h}_\lambda\|^2 \right\}

where :math:`\alpha > 0` controls smoothness-constraint tradeoff.

**Solution:**

Taking gradient and setting to zero:

.. math::

   (\mathbf{A}^\mathsf{T}\mathbf{A} + \alpha\mathbf{D}^\mathsf{T}\mathbf{D})\mathbf{h}_\lambda^* = \mathbf{A}^\mathsf{T}\mathbf{b}

Solve using sparse linear solvers (conjugate gradient or Cholesky decomposition).

Physical Applications
---------------------

**Effective Aging Rate:**

The derivative of material time gives the local aging rate:

.. math::

   \dot{h}_\lambda(t) = \frac{d h_\lambda}{dt}

- :math:`\dot{h} > 1`: System ages faster than real time (accelerated dynamics)
- :math:`\dot{h} = 1`: System ages at normal rate (equilibrium)
- :math:`\dot{h} < 1`: System ages slower than real time (slowed dynamics)

**Cavity Effects on Aging:**

Comparing :math:`h_\lambda(t)` for different coupling strengths :math:`\lambda` reveals how cavity coupling modifies the aging process.

Best Practices
==============

Statistical Convergence
------------------------

**Trajectory Number:**

- Equilibrium properties: :math:`N_\mathrm{T} \geq 10-50` typically sufficient
- Non-equilibrium dynamics: :math:`N_\mathrm{T} \geq 100-500` recommended
- Correlations: More trajectories improve statistical quality

**Time Sampling:**

- Save frequently enough to resolve fastest timescale of interest
- For correlation functions: :math:`\Delta t \ll` shortest correlation time
- Balance storage vs resolution

Correlation Function Tips
--------------------------

**Window Size:**

- Too short: Poor statistics, noisy correlations
- Too long: Misses decorrelation, wastes computation
- **Rule of thumb**: Window should span 5-10 correlation times

**k-Averaging:**

For structure factor and ISF, average over many :math:`\mathbf{k}` vectors:

- Fibonacci sphere distribution for uniform coverage
- :math:`N_{\mathbf{k}} \geq 50` typically sufficient
- Balance isotropy vs computational cost

Material Time Reconstruction
-----------------------------

**Grid Resolution:**

- :math:`M_\mathrm{grid} = 500-1000` for smooth curves
- Finer grid near protocol changes
- Regularization :math:`\alpha \sim 10^{-3}` to :math:`10^{-1}` depending on noise

**Validation:**

- Check that reconstructed :math:`h_\lambda(t)` is monotonically increasing
- Verify constraints are satisfied within tolerance
- Compare aging rates :math:`\dot{h}_\lambda` across conditions

Next Sections
=============

Continue to:

- :doc:`integration_schemes` for numerical implementation
- :doc:`../part3_advanced/correlation_analysis` for practical correlation analysis
- :doc:`../part3_advanced/fdr_temperature` for FDR-based observables
- :doc:`../part1_application/analysis_tools` for analysis workflows

**Related Theory:**

- :doc:`introduction` for fundamental concepts
- :doc:`thermostats` for understanding bath effects on observables
- :doc:`strong_coupling` for polariton effects on dynamics

