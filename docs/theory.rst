======
Theory
======

Cavity HOOMD simulates molecules interacting with optical cavity modes. This enables study of light-matter interactions where molecular vibrations couple strongly to photons. The plugin now supports both constant and time-varying coupling for dynamic experiments.

Mathematical Formulation
========================

**General Multimode Theory**

The general formulation for cavity molecular dynamics involves multiple cavity modes and can be quite complex. In the classical limit, the equations of motion for the coupled nuclei-photonic system are:

**Nuclear Motion:**

.. math::

   M_{nj} \ddot{R}_{nj} = F_{nj}^{(0)} - \sum_{k,\lambda} \left( \tilde{\varepsilon}_{k,\lambda} \tilde{q}_{k,\lambda} + \frac{\tilde{\varepsilon}_{k,\lambda}^2}{m_{k,\lambda}\omega_{k,\lambda}^2} \sum_{l=1}^{N_{\text{sub}}} d_{lg,\lambda} \right) \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}

**Photonic Mode Dynamics:**

.. math::

   m_{k,\lambda} \ddot{\tilde{q}}_{k,\lambda} = -m_{k,\lambda}\omega_{k,\lambda}^2 \tilde{q}_{k,\lambda} - \tilde{\varepsilon}_{k,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}

Where:
- :math:`F_{nj}^{(0)}` is the cavity-free force on each nucleus
- :math:`\tilde{q}_{k,\lambda} = q_{k,\lambda}/\sqrt{N_{\text{cell}}}` is the normalized photonic coordinate
- :math:`\tilde{\varepsilon}_{k,\lambda} = \sqrt{N_{\text{cell}} m_{k,\lambda}\omega_{k,\lambda}^2}/\Omega\varepsilon_0` is the effective coupling strength
- :math:`N_{\text{cell}}` is the number of periodic simulation cells
- :math:`N_{\text{sub}}` is the number of molecules in a single simulation cell
- :math:`d_{ng,\lambda}` is the dipole moment component

Reducing to Single Mode Cavity Dynamics
=======================================

**Single Mode Approximation (κ = 0)**

When we reduce the general multimode cavity molecular dynamics to a single mode case, we make several simplifying assumptions:

1. **Single cavity mode**: We consider only one photonic mode, typically the fundamental mode (κ = 0)
2. **Specific field geometry**: The k-vector points in the z-direction
3. **Polarization considerations**: The summation over λ now represents polarizations in the x and y directions

**Simplified Nuclear Motion Equation**

For the single mode case, the nuclear motion equation becomes:

.. math::

   M_{nj}\ddot{R}_{nj} = F_{nj}^{(0)} - \tilde{\varepsilon}_{0,\lambda}\tilde{q}_{0,\lambda} + \frac{\tilde{\varepsilon}_{0,\lambda}^2}{m_{0,\lambda}\omega_{0,\lambda}^2} \sum_{l=1}^{N_{\text{sub}}} d_{lg,\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}

**Simplified Photonic Mode Dynamics**

The photonic mode equation reduces to:

.. math::

   m_{0,\lambda}\ddot{\tilde{q}}_{0,\lambda} = -m_{0,\lambda}\omega_{0,\lambda}^2 \tilde{q}_{0,\lambda} - \tilde{\varepsilon}_{0,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}

**Key Simplifications**

1. **Mode Index Reduction**
   
   - The general (k,λ) indices reduce to (0,λ)
   - Only the fundamental cavity mode (κ = 0) is retained
   - The k-vector points along z-axis: :math:`\vec{k} = k_z \hat{z}`

2. **Polarization Summation**
   
   With the field propagating in the z-direction, the polarization index λ now specifically refers to:
   
   - **λ = x**: Electric field polarized in x-direction
   - **λ = y**: Electric field polarized in y-direction
   
   The summation :math:`\sum_\lambda` effectively becomes :math:`\sum_{\lambda=x,y}`, representing the two transverse polarization modes.

3. **Coupling Strength Simplification**
   
   The effective coupling strength becomes:
   
   .. math::
   
      \tilde{\varepsilon}_{0,\lambda} = \sqrt{\frac{N_{\text{cell}} m_{0,\lambda}\omega_{0,\lambda}^2}{\Omega\varepsilon_0}}

4. **Physical Interpretation**
   
   - :math:`\tilde{q}_{0,\lambda}`: Normalized coordinate for the single cavity mode with polarization λ
   - :math:`d_{ng,\lambda}`: Dipole moment component of molecule n in direction λ (x or y)
   - The coupling now involves only the transverse components of the molecular dipole moments

**Resulting Dynamics**

This single-mode approximation captures the essential physics of cavity-molecule coupling while dramatically simplifying the computational complexity. The system now describes:

1. **Nuclear motion** influenced by a single cavity mode through the gradient of dipole-field coupling
2. **Single photonic mode** driven by the collective dipole moment of all molecules in the transverse directions
3. **Coherent coupling** between molecular vibrations and the cavity photon mode

This reduction is particularly useful for studying strong coupling phenomena, polariton formation, and cavity-enhanced molecular dynamics in systems where one cavity mode dominates the interaction.

Time-Varying Coupling Theory
=============================

**NEW: Dynamic Coupling Control**

The plugin now supports time-varying coupling strength, enabling study of dynamic cavity-molecule interactions. The coupling strength becomes a function of time:

.. math::

   \tilde{\varepsilon}_{0,\lambda}(t) = g(t) \sqrt{\frac{N_{\text{cell}} m_{0,\lambda}\omega_{0,\lambda}^2}{\Omega\varepsilon_0}}

where :math:`g(t)` is a time-dependent coupling function.

**Step Function Coupling**

The most commonly implemented time-varying coupling is a step function:

.. math::

   g(t) = \begin{cases}
   0 & \text{if } t < t_{\text{switch}} \\
   g_{\text{target}} & \text{if } t \geq t_{\text{switch}}
   \end{cases}

This simulates experiments where cavity coupling is suddenly activated.

**Modified Equations of Motion**

With time-varying coupling, the equations become:

**Nuclear Motion:**

.. math::

   M_{nj}\ddot{R}_{nj} = F_{nj}^{(0)} - g(t)\tilde{\varepsilon}_{0,\lambda}^{(0)}\tilde{q}_{0,\lambda} + \frac{g(t)^2(\tilde{\varepsilon}_{0,\lambda}^{(0)})^2}{m_{0,\lambda}\omega_{0,\lambda}^2} \sum_{l=1}^{N_{\text{sub}}} d_{lg,\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}

**Photonic Mode Dynamics:**

.. math::

   m_{0,\lambda}\ddot{\tilde{q}}_{0,\lambda} = -m_{0,\lambda}\omega_{0,\lambda}^2 \tilde{q}_{0,\lambda} - g(t)\tilde{\varepsilon}_{0,\lambda}^{(0)} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}

where :math:`\tilde{\varepsilon}_{0,\lambda}^{(0)}` is the base coupling strength.

**Energy Conservation in Time-Varying Systems**

Total energy is conserved during coupling switching:

.. math::

   E_{\text{total}}(t) = E_{\text{molecular}}(t) + \frac{1}{2}K\tilde{q}_{0,\lambda}^2(t) + g(t)\tilde{\varepsilon}_{0,\lambda}^{(0)}\tilde{q}_{0,\lambda}(t) \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}(t) + \frac{g(t)^2(\tilde{\varepsilon}_{0,\lambda}^{(0)})^2}{2K} \left(\sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}(t)\right)^2

**Finite-q Cavity Particle Displacement**

When coupling switches from 0 to :math:`g_{\text{target}}` in finite-q mode, the cavity particle position jumps to its new equilibrium:

.. math::

   \vec{q}_{\text{new}} = -\frac{g_{\text{target}}}{K} \vec{d}_{\text{total}}(\vec{r}_{\text{molecular}})

where :math:`\vec{d}_{\text{total}}` is the total molecular dipole moment and :math:`K = m_{0,\lambda}\omega_{0,\lambda}^2`.

**Dissipation in Time-Varying Systems**

Optional cavity dissipation can also be time-varying:

.. math::

   m_{0,\lambda}\ddot{\tilde{q}}_{0,\lambda} = -m_{0,\lambda}\omega_{0,\lambda}^2 \tilde{q}_{0,\lambda} - \gamma(t)\dot{\tilde{q}}_{0,\lambda} - g(t)\tilde{\varepsilon}_{0,\lambda}^{(0)} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}

where :math:`\gamma(t)` is the time-dependent damping coefficient.

**Physical Interpretation of Time-Varying Coupling**

Time-varying coupling simulates realistic experimental scenarios:

1. **Pump-probe experiments**: Cavity coupling activated by external laser pulse
2. **Cavity tuning**: Real-time adjustment of cavity-molecule resonance
3. **Switching dynamics**: Study of non-equilibrium cavity-molecule interactions
4. **Relaxation processes**: System equilibration after sudden coupling changes

**Advantages of Time-Varying Coupling**

- **Dynamic studies**: Investigate non-equilibrium cavity-molecule dynamics
- **Energy redistribution**: Study energy flow between molecular and cavity modes
- **Switching protocols**: Optimize coupling activation strategies
- **Realistic simulations**: Model experimental pump-probe scenarios

Physical Interpretation
=======================

**Why Single Mode?**

1. **Computational Efficiency**: Single-mode calculations are much faster than multimode
2. **Dominant Physics**: The lowest cavity mode often dominates the coupling
3. **Proof of Concept**: Establishes the methodology before extending to multimode

**Limitations of Single-Mode Approximation**

- Cannot capture mode-specific effects
- Misses higher-order mode coupling
- Limited to cavities where one mode dominates

**When Single-Mode is Valid**

- Strong coupling to fundamental cavity mode
- Molecular transitions resonant with lowest mode
- Initial studies of cavity effects

Strong Coupling Regime
======================

**Collective Coupling**

Even in single-mode, collective effects emerge when many molecules couple to the same cavity mode:

.. math::

   g_{\text{eff}} = g \sqrt{N}

where :math:`N` is the number of molecules.

**Energy Scales**

The system enters strong coupling when:

.. math::

   g_{\text{eff}} > \sqrt{\gamma \kappa}

where :math:`\gamma` is molecular damping and :math:`\kappa` is cavity loss rate.

**Polariton Formation**

Strong coupling creates hybrid light-matter states (polaritons) with energy splitting:

.. math::

   \Omega_R = 2g_{\text{eff}}

**Time-Varying Strong Coupling**

With time-varying coupling, the Rabi splitting becomes time-dependent:

.. math::

   \Omega_R(t) = 2g(t)\sqrt{N}

This enables study of polariton dynamics during coupling activation.

Applications and Observables
============================

**What You Can Study**

With the single-mode implementation:
- Fundamental cavity-molecule coupling effects
- Energy transfer between molecules and single cavity mode
- Collective vibrational strong coupling
- Modified molecular dynamics under cavity influence
- **NEW**: Dynamic coupling switching and non-equilibrium processes
- **NEW**: Cavity particle displacement in finite-q modes
- **NEW**: Energy conservation during coupling transitions

**Typical Parameters**

- Coupling strength: :math:`10^{-5}` to :math:`10^{-2}` (atomic units)
- Cavity frequency: 1000-3000 cm⁻¹ (molecular vibrations)
- Temperature: 50-300 K
- **NEW**: Switch times: 0.1-10 ps for dynamic experiments

**Energy Conservation**

The total energy is conserved, including during coupling switching:

.. math::

   E_{\text{total}} = E_{\text{molecular}} + E_{\text{cavity}} + E_{\text{coupling}}

**Tracked Observables**

- Individual energy components (harmonic, coupling, dipole self-energy)
- Cavity mode position and momentum
- Molecular trajectory and total dipole moment
- Energy conservation and thermodynamic quantities
- **NEW**: Time-resolved cavity particle displacement
- **NEW**: Energy redistribution during coupling switching
- **NEW**: Polariton formation dynamics

**Advanced Analysis Features**

- **Energy tracking**: Detailed monitoring of all energy components
- **Correlation functions**: F(k,t) density correlation analysis
- **Performance metrics**: Simulation efficiency and scaling
- **Adaptive timestep**: Automatic timestep optimization
- **GPU acceleration**: High-performance computing support

This enhanced single-mode framework provides a comprehensive platform for studying both equilibrium and non-equilibrium cavity-molecular coupling effects in realistic molecular systems, while maintaining computational tractability and rigorous energy conservation. 