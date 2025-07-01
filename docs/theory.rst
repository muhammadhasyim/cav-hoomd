======
Theory
======

Cavity HOOMD simulates molecules interacting with optical cavity modes. This enables study of light-matter interactions where molecular vibrations couple strongly to photons.

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
   - The k-vector points along z-axis: **k⃗ = k_z ẑ**

2. **Polarization Summation**
   
   With the field propagating in the z-direction, the polarization index λ now specifically refers to:
   
   - **λ = x**: Electric field polarized in x-direction
   - **λ = y**: Electric field polarized in y-direction
   
   The summation ∑_λ effectively becomes ∑_{λ=x,y}, representing the two transverse polarization modes.

3. **Coupling Strength Simplification**
   
   The effective coupling strength becomes:
   
   .. math::
   
      \tilde{\varepsilon}_{0,\lambda} = \sqrt{\frac{N_{\text{cell}} m_{0,\lambda}\omega_{0,\lambda}^2}{\Omega\varepsilon_0}}

4. **Physical Interpretation**
   
   - **q̃_{0,λ}**: Normalized coordinate for the single cavity mode with polarization λ
   - **d_{ng,λ}**: Dipole moment component of molecule n in direction λ (x or y)
   - The coupling now involves only the transverse components of the molecular dipole moments

**Resulting Dynamics**

This single-mode approximation captures the essential physics of cavity-molecule coupling while dramatically simplifying the computational complexity. The system now describes:

1. **Nuclear motion** influenced by a single cavity mode through the gradient of dipole-field coupling
2. **Single photonic mode** driven by the collective dipole moment of all molecules in the transverse directions
3. **Coherent coupling** between molecular vibrations and the cavity photon mode

This reduction is particularly useful for studying strong coupling phenomena, polariton formation, and cavity-enhanced molecular dynamics in systems where one cavity mode dominates the interaction.

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

Applications and Observables
============================

**What You Can Study**

With the single-mode implementation:
- Fundamental cavity-molecule coupling effects
- Energy transfer between molecules and single cavity mode
- Collective vibrational strong coupling
- Modified molecular dynamics under cavity influence

**Typical Parameters**

- Coupling strength: :math:`10^{-5}` to :math:`10^{-2}` (atomic units)
- Cavity frequency: 1000-3000 cm⁻¹ (molecular vibrations)
- Temperature: 50-300 K

**Energy Conservation**

The total energy is conserved:

.. math::

   E_{\text{total}} = E_{\text{molecular}} + E_{\text{cavity}} + E_{\text{coupling}}

**Tracked Observables**

- Individual energy components (harmonic, coupling, dipole self-energy)
- Cavity mode position and momentum
- Molecular trajectory and total dipole moment
- Energy conservation and thermodynamic quantities

This single-mode framework provides a computationally tractable approach to study cavity quantum electrodynamics effects in realistic molecular systems, while maintaining the essential physics of light-matter strong coupling. 