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

**Our Single-Mode Implementation**

Cavity HOOMD implements a **simplified single-mode approximation** focusing on the lowest cavity mode. This reduces the complexity significantly while capturing the essential physics.

**Single-Mode Equations**

For a single cavity mode (dropping indices :math:`k,\lambda`), the interaction Hamiltonian becomes:

.. math::

   H = \frac{1}{2} K q^2 + g \vec{q} \cdot \vec{d} + \frac{g^2}{2K} d^2

Where:
- :math:`q` is the cavity mode position (photon displacement)
- :math:`d` is the total molecular dipole moment  
- :math:`g` is the coupling strength
- :math:`K = m \omega_c^2` is the cavity spring constant

**Forces in Single-Mode Case**

Molecular particles experience forces:

.. math::

   F_i = F_i^{\text{mol}} - g \, q_i \, q - \frac{g^2}{K} q_i \vec{d}

The cavity mode evolves according to:

.. math::

   m \ddot{q} = -K q - g \vec{d}

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