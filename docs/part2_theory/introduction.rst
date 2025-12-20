============
Introduction
============

This section provides the theoretical foundation for cavity-coupled molecular dynamics simulations in Cavity HOOMD.

.. contents:: In this Section
   :local:
   :depth: 2

Overview of Cavity QED in Molecular Dynamics
=============================================

**Cavity Quantum Electrodynamics (Cavity QED)** studies the interaction between light confined in an optical cavity and matter. When molecular vibrations couple strongly to cavity photon modes, hybrid light-matter states called **polaritons** form.

Classical Treatment
-------------------

Cavity HOOMD uses a **classical description** of both nuclear and photonic degrees of freedom:

- **Nuclear motion:** Classical trajectories governed by Newton's equations
- **Photonic modes:** Harmonic oscillators with classical coordinates and momenta
- **Coupling:** Bilinear dipole-field interaction

This approach is valid when:

1. Many photons populate the cavity mode (large n)
2. Thermal energy :math:`k_B T` exceeds photon energy :math:`\hbar\omega`
3. Nuclear motion is classical (heavy nuclei)

**Quantum vs Classical:**

+-------------------------+-------------------------+-------------------------+
| Aspect                  | Quantum                 | Classical (Cavity HOOMD)|
+=========================+=========================+=========================+
| Nuclear dynamics        | Wave packets            | Trajectories            |
+-------------------------+-------------------------+-------------------------+
| Photonic modes          | Fock states             | Harmonic oscillators    |
+-------------------------+-------------------------+-------------------------+
| Coupling                | :math:`\hat{H}`         | Classical Hamiltonian   |
+-------------------------+-------------------------+-------------------------+
| Validity                | T → 0, few photons      | T > 0, many photons     |
+-------------------------+-------------------------+-------------------------+

Theoretical Foundation
======================

Hamiltonian Derivation Overview
---------------------------------

The cavity-molecule Hamiltonian in Cavity HOOMD is derived from quantum electrodynamics and reduced to a classical form through several approximations. This section provides a high-level overview; see :doc:`cavity_forces` for the complete derivation.

**Starting Point: Minimal Coupling Hamiltonian**

The quantum description begins with charged particles (electrons and nuclei) coupled to the electromagnetic field through the vector potential:

.. math::

   \hat{H} = \hat{H}_\mathrm{M} + \hat{H}_\mathrm{EM}

where the matter Hamiltonian includes kinetic energy with minimal coupling:

.. math::

   \hat{H}_\mathrm{M} = \sum_i \frac{|\hat{\mathbf{P}}_i - Z_i e \hat{\mathbf{A}}(\mathbf{R}_i)|^2}{2M_i} + \sum_j \frac{|\hat{\mathbf{p}}_j + e \hat{\mathbf{A}}(\mathbf{r}_j)|^2}{2m_j} + \hat{V}

and the electromagnetic field energy is:

.. math::

   \hat{H}_\mathrm{EM} = \frac{\epsilon_0}{2} \int d^3\mathbf{r} \left( |\hat{\mathbf{E}}|^2 + c^2|\hat{\mathbf{B}}|^2 \right)

**Key Approximations:**

1. **Coulomb Gauge**: :math:`\nabla \cdot \hat{\mathbf{A}} = 0`, electromagnetic fields are purely transverse

2. **Mode Expansion**: The vector potential is expanded in cavity normal modes with frequencies :math:`\omega_\ell = \ell \pi c / L`

3. **Uniform Field Approximation**: Since molecular length scales :math:`\ll` cavity length, :math:`\hat{\mathbf{A}}(\mathbf{r}) \approx \hat{\mathbf{A}}(0)` is spatially uniform

4. **Gauge Transformations**: Two unitary transformations are applied:
   
   - Remove vector potential from kinetic energy (length gauge)
   - 90° phase rotation (coordinate displacement form)

5. **Pauli-Fierz Form**: After transformations, the Hamiltonian becomes:

.. math::

   \hat{H}_\mathrm{PF} = \hat{H}_\mathrm{M} + \frac{1}{2}\sum_{\ell,\alpha} \left[ \hat{p}_{\ell,\alpha}^2 + \omega_\ell^2 \left( \hat{q}_{\ell,\alpha} + \frac{\lambda \hat{\mu}_\alpha}{\omega_\ell} \right)^2 \right]

where :math:`\lambda = 1/\sqrt{\epsilon_0 V}` is the light-matter coupling strength and :math:`\hat{\mu}_\alpha` is the molecular dipole moment operator.

Cavity Born-Oppenheimer Approximation
---------------------------------------

**Timescale Separation:**

Electrons move much faster than nuclei, but electromagnetic modes typically match nuclear vibrational timescales. We therefore:

1. **Solve electronic structure separately** for nuclear positions :math:`\{\mathbf{R}_i\}` and cavity coordinates :math:`\{q_{\ell,\alpha}\}`:

.. math::

   \hat{H}_\mathrm{e} |\psi_g\rangle = V_\mathrm{BO}(\{\mathbf{R}_i, q_{\ell,\alpha}\}) |\psi_g\rangle

2. **Invoke adiabatic approximation**: Since cavity modes are typically far off-resonant from electronic transitions:

.. math::

   V_\mathrm{BO}(\{\mathbf{R}_i, q_{\ell,\alpha}\}) \approx V_\mathrm{BO}(\{\mathbf{R}_i\})

This avoids treating polaritonic potential energy surfaces while retaining cavity effects on nuclear dynamics.

Classical Limit and Mean-Field Approximation
----------------------------------------------

**From Quantum to Classical:**

Treating slower degrees of freedom (nuclei and photons) classically, the potential energy becomes:

.. math::

   V(\{\mathbf{R}_i, q_{\ell,\alpha}\}) = V_\mathrm{BO}(\{\mathbf{R}_i\}) + \sum_{\ell,\alpha} \left[ \frac{\omega_\ell^2 q_{\ell,\alpha}^2}{2} + \lambda \omega_\ell \langle\psi_g|\hat{\mu}_\alpha|\psi_g\rangle q_{\ell,\alpha} + \frac{\lambda^2}{2}\langle\psi_g|\hat{\mu}_\alpha^2|\psi_g\rangle \right]

**Mean-Field Treatment:**

The quantum fluctuation term :math:`\langle\psi_g|\hat{\mu}_\alpha^2|\psi_g\rangle` is approximated by its classical mean-field value:

.. math::

   \langle\psi_g|\hat{\mu}_\alpha^2|\psi_g\rangle \approx \left(\langle\psi_g|\hat{\mu}_\alpha|\psi_g\rangle\right)^2 = \mu_\alpha^2

where :math:`\mu_\alpha(\{\mathbf{R}_i\})` is the classical dipole moment.

**Final Classical Hamiltonian:**

.. math::

   H = \sum_i \frac{\mathbf{P}_i^2}{2M_i} + V_\mathrm{BO}(\{\mathbf{R}_i\}) + \sum_{\ell,\alpha} \left[ \frac{p_{\ell,\alpha}^2}{2} + \frac{\omega_\ell^2}{2}\left(q_{\ell,\alpha} + \frac{\lambda\mu_\alpha}{\omega_\ell}\right)^2 \right]

This is the starting point for cavity molecular dynamics simulations.

Single-Mode vs Multi-Mode Treatment
-------------------------------------

**Single-Mode Approximation:**

In Cavity HOOMD, we typically retain only the lowest-frequency electromagnetic mode:

.. math::

   \omega_{\ell=1} = \frac{2\pi c}{L} = \omega_\mathrm{c}

**Justification:**

- **Energy scale separation**: Higher modes have :math:`\omega_\ell > \omega_c` and are off-resonant with molecular vibrations
- **Coupling strength**: Higher modes couple more weakly due to larger frequency denominators
- **Computational efficiency**: Single mode reduces degrees of freedom from :math:`\mathcal{O}(N_\text{modes})` to :math:`\mathcal{O}(1)`

**Multi-Mode Considerations:**

Multi-mode treatment is necessary when:

- Multiple cavity modes are near-resonant with molecular transitions
- Studying mode-selective effects
- Investigating multimode interference phenomena

The formalism in :doc:`cavity_forces` naturally extends to multiple modes.

Scope and Applicability
========================

What Cavity HOOMD Can Simulate
-------------------------------

**Strong Coupling Regime:**

- Collective vibrational strong coupling
- Polariton formation and dynamics
- Cavity-modified reaction rates
- Energy transfer in cavity-coupled systems
- Non-equilibrium switching experiments

**System Types:**

- Molecular liquids and solids
- Organic semiconductors
- Biomolecular systems
- Any polarizable molecular system

**Temperature Range:**

- 10 K to 500 K (typical)
- Limited by MD force field validity
- Best suited for T > 50 K where classical treatment is valid

Limitations
-----------

**1. Classical Approximation:**

Cannot capture:

- Zero-point energy effects
- Quantum coherence
- Tunneling
- Photon number quantization

**2. Single-Mode Approximation:**

- Only one or few cavity modes
- No multimode interference effects
- No cavity dispersion

**3. No Electronic Excitations:**

- Ground state electronic structure only
- No photochemistry or excited state dynamics
- Coupling through nuclear dipole moments only

**4. Periodic Boundary Conditions:**

- Infinite periodic system
- No cavity boundaries or mirrors
- Uniform field approximation (q=0) or single standing wave (finite-q)

Units and Conventions
=====================

Atomic Units
------------

Cavity HOOMD uses **HOOMD reduced units** internally, which are related to atomic units:

**Fundamental Constants:**

- Electron mass: :math:`m_e = 1`
- Elementary charge: :math:`e = 1`
- Reduced Planck constant: :math:`\hbar = 1`
- Coulomb constant: :math:`k_e = 1`

**Derived Units:**

+-------------------------+-------------------------------+-------------------------+
| Quantity                | Atomic Unit                   | SI Equivalent           |
+=========================+===============================+=========================+
| Length                  | Bohr radius :math:`a_0`       | 0.529 Å                 |
+-------------------------+-------------------------------+-------------------------+
| Energy                  | Hartree :math:`E_h`           | 27.2 eV = 4.36×10⁻¹⁸ J  |
+-------------------------+-------------------------------+-------------------------+
| Time                    | :math:`\hbar/E_h`             | 2.42×10⁻¹⁷ s = 0.024 fs |
+-------------------------+-------------------------------+-------------------------+
| Mass                    | Electron mass :math:`m_e`     | 9.11×10⁻³¹ kg           |
+-------------------------+-------------------------------+-------------------------+
| Temperature             | :math:`E_h/k_B`               | 315,775 K               |
+-------------------------+-------------------------------+-------------------------+
| Frequency               | :math:`E_h/\hbar`             | 4.13×10¹⁶ Hz            |
+-------------------------+-------------------------------+-------------------------+

**Practical Conversions:**

For typical molecular simulations:

.. code-block:: python

   # Energy
   1 eV = 0.0367 Hartree
   1 kcal/mol = 0.00159 Hartree
   
   # Length  
   1 Å = 1.889 Bohr
   1 nm = 18.89 Bohr
   
   # Time
   1 fs = 41.34 atomic time units
   1 ps = 41,341 atomic time units
   
   # Temperature
   300 K = 0.00095 Hartree/k_B
   100 K = 0.00032 Hartree/k_B
   
   # Frequency
   2000 cm⁻¹ = 0.00916 Hartree/ℏ

User Interface Units
--------------------

**Input Parameters (User-Friendly):**

- Temperature: Kelvin (K)
- Cavity frequency: wavenumbers (cm⁻¹)
- Time: picoseconds (ps)
- Length: Angstroms (Å)
- Mass: atomic mass units (amu)

**Internal Computation:**

All converted to reduced/atomic units automatically.

**Output:**

Can be in either reduced units or converted back to SI/conventional units.

Physical Constants Table
=========================

**Useful Physical Constants:**

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Constant
     - Symbol
     - Value
   * - Boltzmann constant
     - :math:`k_B`
     - 1.381×10⁻²³ J/K = 8.617×10⁻⁵ eV/K
   * - Planck constant
     - :math:`h`
     - 6.626×10⁻³⁴ J·s
   * - Reduced Planck constant
     - :math:`\hbar`
     - 1.055×10⁻³⁴ J·s
   * - Speed of light
     - :math:`c`
     - 2.998×10⁸ m/s
   * - Elementary charge
     - :math:`e`
     - 1.602×10⁻¹⁹ C
   * - Avogadro constant
     - :math:`N_A`
     - 6.022×10²³ mol⁻¹
   * - Vacuum permittivity
     - :math:`\varepsilon_0`
     - 8.854×10⁻¹² F/m

**Molecular Constants:**

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Molecule
     - Vibrational Frequency
     - Dipole Moment
   * - O₂
     - 1580 cm⁻¹
     - 0.0 D (non-polar)
   * - CO
     - 2143 cm⁻¹
     - 0.11 D
   * - H₂O
     - 3650 cm⁻¹ (stretch)
     - 1.85 D
   * - CO₂
     - 1333 cm⁻¹ (asymmetric)
     - 0.0 D (linear)

Energy Scales and Typical Values
=================================

Molecular Energy Scales
------------------------

**At room temperature (300 K):**

.. math::

   k_B T \approx 0.026 \text{ eV} \approx 200 \text{ cm}^{-1}

**Vibrational energies:**

- O-H stretch: ~3600 cm⁻¹ ≈ 0.45 eV
- C=O stretch: ~1700 cm⁻¹ ≈ 0.21 eV
- C-H stretch: ~3000 cm⁻¹ ≈ 0.37 eV
- Bending modes: ~1000-1500 cm⁻¹

**Bonding energies:**

- Covalent bonds: 1-5 eV
- Hydrogen bonds: 0.1-0.4 eV
- van der Waals: 0.01-0.1 eV

Cavity Energy Scales
--------------------

**Typical cavity frequencies:**

- Mid-IR: 1000-4000 cm⁻¹ (0.12-0.50 eV)
- Near-IR: 4000-13,000 cm⁻¹ (0.50-1.6 eV)
- Visible: 13,000-25,000 cm⁻¹ (1.6-3.1 eV)

**Single-molecule coupling:**

.. math::

   g_{\text{single}} \sim 10^{-5} \text{ to } 10^{-3} \text{ (atomic units)}

**Collective coupling (N molecules):**

.. math::

   g_{\text{collective}} = g_{\text{single}} \sqrt{N}

For N=512 molecules:

.. math::

   g_{\text{collective}} \approx 22.6 \times g_{\text{single}}

Strong Coupling Criterion
--------------------------

**Condition for strong coupling:**

.. math::

   \Omega_R > \sqrt{\gamma_{\text{mol}} \kappa_{\text{cav}}}

where:

- :math:`\Omega_R = 2g_{\text{collective}}`: Rabi splitting
- :math:`\gamma_{\text{mol}}`: Molecular dephasing rate
- :math:`\kappa_{\text{cav}}`: Cavity loss rate

**Typical values:**

- :math:`\Omega_R \sim 50-200` cm⁻¹ (experimental)
- :math:`\gamma_{\text{mol}} \sim 10-100` cm⁻¹
- :math:`\kappa_{\text{cav}} \sim 1-50` cm⁻¹

**Strong coupling achieved when:** :math:`\Omega_R > 50` cm⁻¹

Parameter Guidelines
====================

Recommended Parameter Ranges
-----------------------------

**Coupling Strength:**

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Regime
     - g (atomic units)
     - Description
   * - Weak
     - 10⁻⁵ - 10⁻⁴
     - Perturbative, minimal modification
   * - Moderate
     - 10⁻⁴ - 10⁻³
     - Observable effects, energy exchange
   * - Strong
     - > 10⁻³
     - Polariton formation, strong coupling

**Temperature:**

- Low T (50-100 K): Reduced thermal fluctuations, clearer coupling effects
- Room T (300 K): Realistic conditions, stronger thermal noise
- High T (>400 K): Challenging for classical MD, force field limitations

**Cavity Frequency:**

- Match molecular vibrations: :math:`\omega_{\text{cav}} \approx \omega_{\text{vib}}`
- Red-detuned: :math:`\omega_{\text{cav}} < \omega_{\text{vib}}`
- Blue-detuned: :math:`\omega_{\text{cav}} > \omega_{\text{vib}}`

**Timestep:**

- Standard: 0.001 ps (1 fs)
- High-frequency vibrations: 0.0005 ps (0.5 fs)
- Adaptive: Use ``AdaptiveTimestepUpdater``

Parameter Relationships
-----------------------

**Cavity Quality Factor:**

.. math::

   Q = \frac{\omega_{\text{cav}}}{\kappa}

where :math:`\kappa` is the cavity loss rate.

**Cavity Lifetime:**

.. math::

   \tau_{\text{cav}} = \frac{1}{\kappa}

**Molecular Dephasing Time:**

.. math::

   T_2 = \frac{1}{\gamma_{\text{mol}}}

**Rabi Period:**

.. math::

   T_{\text{Rabi}} = \frac{2\pi}{\Omega_R}

Next Sections
=============

Continue to:

- :doc:`cavity_forces` for detailed force equations
- :doc:`time_varying_coupling` for time-dependent coupling theory
- :doc:`thermostats` for temperature control methods
- :doc:`energy_conservation` for energy decomposition
- :doc:`strong_coupling` for polariton physics

**Back to:**

- :doc:`../part1_application/getting_started` for practical usage
- :doc:`../index` for main documentation

