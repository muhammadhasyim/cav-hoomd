==============
Cavity Forces
==============

This section describes the mathematical formulation of cavity-molecule coupling and the resulting forces on nuclei and photonic modes.

.. contents:: In this Section
   :local:
   :depth: 2

General Multimode Formulation
==============================

Classical Hamiltonian
---------------------

The complete classical Hamiltonian for cavity-coupled molecular dynamics is:

.. math::

   H = H_{\text{mol}} + \sum_{k,\lambda} H_{\text{cav}}^{k,\lambda}

where :math:`H_{\text{mol}}` is the molecular Hamiltonian and the sum runs over cavity modes k and polarizations λ.

**Molecular Hamiltonian:**

.. math::

   H_{\text{mol}} = \sum_{n,j} \frac{P_{nj}^2}{2M_{nj}} + V(\{R_{nj}\})

Where:

- :math:`P_{nj}`: Momentum of nucleus j in molecule n
- :math:`M_{nj}`: Mass of nucleus j in molecule n
- :math:`R_{nj}`: Position of nucleus j in molecule n
- :math:`V(\{R_{nj}\})`: Molecular potential energy (bonds, angles, LJ, Coulomb, etc.)

**Cavity Mode Hamiltonian:**

.. math::

   H_{\text{cav}}^{k,\lambda} = \frac{1}{2}m_{k,\lambda}\omega_{k,\lambda}^2 \tilde{q}_{k,\lambda}^2 + \frac{1}{2m_{k,\lambda}}\tilde{p}_{k,\lambda}^2 + \tilde{\varepsilon}_{k,\lambda} \tilde{q}_{k,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda} + \frac{\tilde{\varepsilon}_{k,\lambda}^2}{2m_{k,\lambda}\omega_{k,\lambda}^2} \left(\sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}\right)^2

Where:

- :math:`\tilde{q}_{k,\lambda}`: Normalized photonic coordinate for mode k, polarization λ
- :math:`\tilde{p}_{k,\lambda}`: Conjugate momentum
- :math:`m_{k,\lambda}`: Effective mass of photonic mode
- :math:`\omega_{k,\lambda}`: Angular frequency of cavity mode
- :math:`\tilde{\varepsilon}_{k,\lambda}`: Effective coupling strength
- :math:`d_{ng,\lambda}`: Component of molecular dipole moment n in direction λ
- :math:`N_{\text{sub}}`: Number of molecules in simulation cell

**Physical Interpretation:**

1. **First term:** Harmonic oscillator energy of cavity mode
2. **Second term:** Kinetic energy of cavity mode
3. **Third term:** Linear dipole-field coupling
4. **Fourth term:** Dipole self-energy (counter-term ensuring correct energetics)

Equations of Motion
-------------------

**Nuclear Motion:**

From Hamilton's equations:

.. math::

   M_{nj} \ddot{R}_{nj} = F_{nj}^{(0)} - \sum_{k,\lambda} \left( \tilde{\varepsilon}_{k,\lambda} \tilde{q}_{k,\lambda} + \frac{\tilde{\varepsilon}_{k,\lambda}^2}{m_{k,\lambda}\omega_{k,\lambda}^2} \sum_{l=1}^{N_{\text{sub}}} d_{lg,\lambda} \right) \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}

Where:

- :math:`F_{nj}^{(0)} = -\frac{\partial V}{\partial R_{nj}}`: Cavity-free force
- :math:`\frac{\partial d_{ng,\lambda}}{\partial R_{nj}}`: Dipole moment gradient

**Physical Interpretation:**

The nuclear force has three contributions:

1. **Molecular force** :math:`F_{nj}^{(0)}`: Standard MD forces (bonds, LJ, etc.)
2. **Linear coupling force** :math:`-\tilde{\varepsilon}_{k,\lambda} \tilde{q}_{k,\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}`: Force from cavity field acting on molecular dipole
3. **Self-energy force** :math:`-\frac{\tilde{\varepsilon}_{k,\lambda}^2}{K_{k,\lambda}} D_{\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}`: Force from collective dipole interaction (where :math:`K_{k,\lambda} = m_{k,\lambda}\omega_{k,\lambda}^2`)

**Photonic Mode Dynamics:**

.. math::

   m_{k,\lambda} \ddot{\tilde{q}}_{k,\lambda} = -m_{k,\lambda}\omega_{k,\lambda}^2 \tilde{q}_{k,\lambda} - \tilde{\varepsilon}_{k,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}

Where:

- First term: Harmonic restoring force
- Second term: Driving force from molecular dipoles

Single Mode Approximation
==========================

Simplifying to κ=0 Mode
-----------------------

For many applications, considering only the fundamental cavity mode (κ=0) suffices:

**Reduced Hamiltonian:**

.. math::

   H = H_{\text{mol}} + \sum_{\lambda=x,y} \left[ \frac{1}{2}K_\lambda \tilde{q}_{0,\lambda}^2 + \tilde{\varepsilon}_{0,\lambda} \tilde{q}_{0,\lambda} D_\lambda + \frac{\tilde{\varepsilon}_{0,\lambda}^2}{2K_\lambda} D_\lambda^2 \right]

Where:

- :math:`D_\lambda = \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda}`: Total dipole moment in direction λ
- :math:`K_\lambda = m_{0,\lambda}\omega_{0,\lambda}^2`: Cavity spring constant
- Sum over :math:`\lambda = x, y` for transverse polarizations

**Equations of Motion (Single Mode):**

Nuclear motion:

.. math::

   M_{nj}\ddot{R}_{nj} = F_{nj}^{(0)} - \sum_{\lambda} \left( \tilde{\varepsilon}_{0,\lambda}\tilde{q}_{0,\lambda} + \frac{\tilde{\varepsilon}_{0,\lambda}^2}{K_\lambda} D_\lambda \right) \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}

Photonic mode:

.. math::

   m_{0,\lambda}\ddot{\tilde{q}}_{0,\lambda} = -K_\lambda \tilde{q}_{0,\lambda} - \tilde{\varepsilon}_{0,\lambda} D_\lambda

**When is Single Mode Valid?**

1. **Frequency separation:** :math:`\omega_1 - \omega_0 \gg \Omega_R` (next mode far from resonance)
2. **Weak higher modes:** Higher mode couplings negligible
3. **Fundamental dominates:** Most photon population in κ=0 mode

Force Decomposition
===================

Total Force on Nuclei
---------------------

The total force on nucleus j of molecule n can be decomposed as:

.. math::

   F_{nj}^{\text{total}} = F_{nj}^{\text{mol}} + F_{nj}^{\text{coupling}} + F_{nj}^{\text{self}}

**1. Molecular Force:**

.. math::

   F_{nj}^{\text{mol}} = -\frac{\partial V}{\partial R_{nj}}

Standard molecular dynamics forces: bonds, angles, dihedrals, Lennard-Jones, Coulomb, etc.

**2. Coupling Force:**

.. math::

   F_{nj}^{\text{coupling}} = -\sum_{\lambda} \tilde{\varepsilon}_{0,\lambda} \tilde{q}_{0,\lambda} \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}

This force arises from the cavity electric field acting on the molecular dipole moment.

**Magnitude:** Scales as :math:`g \times q \times |\nabla d|`

**Direction:** Along dipole moment gradient (tends to align dipole with field)

**3. Self-Energy Force:**

.. math::

   F_{nj}^{\text{self}} = -\sum_{\lambda} \frac{\tilde{\varepsilon}_{0,\lambda}^2}{K_\lambda} D_\lambda \frac{\partial d_{ng,\lambda}}{\partial R_{nj}}

This force accounts for the interaction energy of the molecular dipole with the collective dipole.

**Magnitude:** Scales as :math:`g^2 \times N \times |\nabla d|`

**Physical Meaning:** Represents mean-field interaction with all other molecules mediated by the cavity

Energy Components
-----------------

**Total Energy Decomposition:**

.. math::

   E_{\text{total}} = E_{\text{kinetic}} + E_{\text{potential}} + E_{\text{cavity}} + E_{\text{coupling}} + E_{\text{self}}

**1. Kinetic Energy:**

.. math::

   E_{\text{kinetic}} = \sum_{n,j} \frac{1}{2}M_{nj}\dot{R}_{nj}^2 + \sum_{\lambda} \frac{1}{2}m_{0,\lambda}\dot{\tilde{q}}_{0,\lambda}^2

Includes both molecular and photonic kinetic energy.

**2. Potential Energy:**

.. math::

   E_{\text{potential}} = V(\{R_{nj}\})

Standard molecular potential energy.

**3. Cavity Harmonic Energy:**

.. math::

   E_{\text{cavity}} = \sum_{\lambda} \frac{1}{2}K_\lambda \tilde{q}_{0,\lambda}^2

Energy stored in cavity mode oscillations.

**4. Coupling Energy:**

.. math::

   E_{\text{coupling}} = \sum_{\lambda} \tilde{\varepsilon}_{0,\lambda} \tilde{q}_{0,\lambda} D_\lambda

Interaction energy between molecular dipoles and cavity field.

**Sign:** Can be positive or negative depending on relative orientation.

**5. Self-Energy:**

.. math::

   E_{\text{self}} = \sum_{\lambda} \frac{\tilde{\varepsilon}_{0,\lambda}^2}{2K_\lambda} D_\lambda^2

Collective dipole-dipole interaction energy.

**Always positive:** Represents energetic cost of dipole interactions.

Coupling Strength Parameter
============================

Defining the Coupling
---------------------

The effective coupling strength is:

.. math::

   \tilde{\varepsilon}_{0,\lambda} = \sqrt{\frac{N_{\text{cell}} m_{0,\lambda}\omega_{0,\lambda}^2}{\Omega\varepsilon_0}}

Where:

- :math:`N_{\text{cell}}`: Number of periodic simulation cells
- :math:`\Omega`: Volume of simulation cell
- :math:`\varepsilon_0`: Vacuum permittivity

**In practical calculations:**

.. math::

   \tilde{\varepsilon}_{0,\lambda} = g \times \text{(conversion factors)}

where g is the user-specified coupling strength parameter.

Single-Molecule vs Collective Coupling
---------------------------------------

**Single-molecule coupling:**

.. math::

   g_{\text{single}} = \frac{\tilde{\varepsilon}_{0,\lambda}}{\sqrt{N}}

**Collective coupling:**

.. math::

   g_{\text{collective}} = g_{\text{single}} \times \sqrt{N}

For N molecules, the collective enhancement is :math:`\sqrt{N}`.

**Example:**

- N = 512 molecules
- :math:`g_{\text{single}} = 10^{-3}` (atomic units)
- :math:`g_{\text{collective}} = 10^{-3} \times \sqrt{512} \approx 0.0226`

This collective enhancement enables strong coupling even with modest single-molecule coupling.

Finite-q Cavity Mode
=====================

Wave Vector Considerations
---------------------------

For finite momentum transfer :math:`\vec{k} \neq 0`:

**Modified Coupling:**

.. math::

   H_{\text{coupling}} = \sum_{\lambda} \tilde{\varepsilon}_{0,\lambda} \tilde{q}_{0,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda} e^{i\vec{k}\cdot\vec{R}_n}

where :math:`\vec{R}_n` is the position of molecule n.

**For standing wave along z-axis:**

.. math::

   \vec{k} = (0, 0, k_z), \quad k_z = \frac{2\pi}{L_z}

where :math:`L_z` is the box size in z-direction.

**Phase factor:**

.. math::

   e^{i k_z z_n} = \cos(k_z z_n) + i\sin(k_z z_n)

For real fields, only cosine term contributes:

.. math::

   H_{\text{coupling}} = \sum_{\lambda} \tilde{\varepsilon}_{0,\lambda} \tilde{q}_{0,\lambda} \sum_{n=1}^{N_{\text{sub}}} d_{ng,\lambda} \cos(k_z z_n)

Cavity Particle Representation
-------------------------------

**In Cavity HOOMD implementation:**

The cavity mode is represented by physical "cavity particles" with:

- **q=0 mode:** Particles at origin, displacement represents mode amplitude
- **Finite-q mode:** Particles positioned to represent standing wave

**Position-mode relationship:**

.. math::

   \tilde{q}_{0,\lambda}(t) = \vec{r}_{\text{cavity}}(t) \cdot \hat{e}_\lambda

where :math:`\hat{e}_\lambda` is the polarization direction unit vector.

**Advantages:**

1. Uses standard MD machinery
2. Easy visualization
3. Natural integration with HOOMD-blue
4. Straightforward implementation of time-varying coupling

Dipole Moment Calculation
==========================

For Diatomic Molecules
-----------------------

For a diatomic molecule with atoms at :math:`\vec{r}_1, \vec{r}_2`:

**Dipole moment:**

.. math::

   \vec{d} = q_1 \vec{r}_1 + q_2 \vec{r}_2 = q_{\text{eff}}(\vec{r}_2 - \vec{r}_1)

where :math:`q_{\text{eff}}` is the effective charge.

**For O₂:** Approximately non-polar, but instantaneous dipole from charge fluctuations.

**Dipole gradient:**

.. math::

   \frac{\partial d_\lambda}{\partial r_{1,j}} = -q_{\text{eff}} \delta_{\lambda j}, \quad \frac{\partial d_\lambda}{\partial r_{2,j}} = +q_{\text{eff}} \delta_{\lambda j}

For Polyatomic Molecules
-------------------------

**General formula:**

.. math::

   \vec{d} = \sum_{i=1}^{N_{\text{atoms}}} q_i \vec{r}_i

**Gradient:**

.. math::

   \frac{\partial d_\lambda}{\partial r_{i,j}} = q_i \delta_{\lambda j}

**Practical Calculation:**

In simulations, dipole moments are computed from:

1. **Partial charges:** From force field
2. **Positions:** From trajectory
3. **Periodic boundary conditions:** Must use unwrapped coordinates

Implementation Notes
====================

Numerical Considerations
------------------------

**1. Timestep Selection:**

Cavity oscillations impose constraint:

.. math::

   \Delta t \ll \frac{2\pi}{\omega_{\text{cavity}}}

**Guideline:** :math:`\Delta t < 0.1 \times T_{\text{cavity}}`

For :math:`\omega = 2000` cm⁻¹:

.. math::

   T_{\text{cavity}} \approx 16.7 \text{ fs} \implies \Delta t < 1.7 \text{ fs}

**2. Force Cutoffs:**

Dipole gradients can be large near bonds. Use:

- Smooth switching functions
- Force capping for initialization
- Gradual coupling ramp-up

**3. Energy Conservation:**

Monitor relative energy drift:

.. math::

   \Delta E_{\text{rel}} = \frac{|E(t) - E(0)|}{|E(0)|} < 10^{-4}

For well-behaved simulations.

GPU Optimization
----------------

**Memory Layout:**

- Coalesced memory access for particle positions
- Shared memory for reductions (total dipole calculation)
- Texture memory for periodic wrapping

**Kernel Structure:**

1. **Dipole calculation:** All molecules → total :math:`D_\lambda`
2. **Force calculation:** Distribute cavity force to all nuclei
3. **Integration:** Update cavity particle positions/velocities

**Performance:**

Typical speedup: 50-100× over CPU for N>1000 molecules.

Next Sections
=============

Continue to:

- :doc:`time_varying_coupling` for time-dependent coupling theory
- :doc:`energy_conservation` for energy decomposition details
- :doc:`thermostats` for temperature control
- :doc:`strong_coupling` for polariton physics

**Related practical guides:**

- :doc:`../part1_application/running_simulations` for parameter selection
- :doc:`../part1_application/analysis_tools` for energy analysis

