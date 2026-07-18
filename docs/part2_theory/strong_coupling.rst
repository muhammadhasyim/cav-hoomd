================
Strong Coupling
================

This section describes the strong coupling regime where cavity and molecular excitations hybridize to form polaritons.

.. contents:: In this Section
   :local:
   :depth: 2

Defining Strong Coupling
=========================

Collective Enhancement
----------------------

For N molecules coupling to a single cavity mode:

.. math::

   g_{\text{eff}} = g_{\text{single}} \sqrt{N}

Where:

- :math:`g_{\text{single}}`: Single-molecule coupling strength
- :math:`N`: Number of molecules
- :math:`g_{\text{eff}}`: Effective collective coupling

**Example:**

- N = 512 molecules
- :math:`g_{\text{single}} = 10^{-3}`
- :math:`g_{\text{eff}} = 10^{-3} \times \sqrt{512} = 0.0226`

This :math:`\sqrt{N}` enhancement is key to achieving strong coupling.

Strong Coupling Criterion
--------------------------

**Condition:**

.. math::

   \Omega_R > \sqrt{\gamma_{\text{mol}} \kappa_{\text{cav}}}

Where:

- :math:`\Omega_R = 2g_{\text{eff}}`: Rabi splitting
- :math:`\gamma_{\text{mol}}`: Molecular dephasing rate
- :math:`\kappa_{\text{cav}}`: Cavity photon loss rate

**Physical meaning:** Coherent energy exchange faster than decoherence.

Polariton Formation
===================

Hybrid Light-Matter States
---------------------------

In strong coupling, cavity photons and molecular excitations hybridize:

.. math::

   |\text{polariton}\rangle = \alpha |\text{photon}\rangle + \beta |\text{molecule}\rangle

Creating two new eigenstates:

- **Upper polariton (UP):** :math:`E_+ = E_0 + \frac{\Omega_R}{2}`
- **Lower polariton (LP):** :math:`E_- = E_0 - \frac{\Omega_R}{2}`

Where :math:`E_0` is the uncoupled resonance energy.

Rabi Splitting
--------------

**Energy splitting:**

.. math::

   \Omega_R = 2g_{\text{eff}} = 2g_{\text{single}}\sqrt{N}

**Observable in:**

- Transmission/reflection spectra
- Photoluminescence
- Time-resolved dynamics

**Typical experimental values:**

- Vibrational strong coupling: 50-200 cm⁻¹
- Electronic strong coupling: 100-1000 meV

Cavity-Modified Dynamics
=========================

Modified Energy Landscape
-------------------------

Strong coupling modifies the potential energy surface. Effective potentials include cavity contribution.

Rate Modifications
------------------

Chemical reaction rates can be modified under strong coupling through:

1. **Energy redistribution**
2. **Modified transition states**
3. **Collective effects**

This is an active area of research ("polariton chemistry").

Resonance vs Off-Resonance Effects
===================================

Polariton Formation Requirement
--------------------------------

The mechanism described above relies on **resonance** between the cavity mode and intramolecular vibrations. When the cavity frequency :math:`\omega_\mathrm{c}` matches a molecular vibrational frequency :math:`\omega_\mathrm{vib}`:

.. math::

   \omega_\mathrm{c} \approx \omega_\mathrm{vib}

polaritons form and cavity effects are maximized.

**Resonance Condition:**

For effective polariton formation:

.. math::

   |\omega_\mathrm{c} - \omega_\mathrm{vib}| < \Omega_R

where :math:`\Omega_R` is the Rabi splitting. Outside this window, polariton formation is suppressed.

Off-Resonance Behavior
-----------------------

**Frequency Dependence:**

To test whether polariton formation is essential for cavity effects, one can scan cavity frequencies around molecular vibrational frequencies. For a fixed coupling strength :math:`\lambda` and varying :math:`\omega_\mathrm{c}`:

**Expected if polaritons are essential:**

- Sharp peak in observable at :math:`\omega_\mathrm{c} = \omega_\mathrm{vib}`
- Rapid decay to baseline at off-resonant frequencies
- Effects vanish when :math:`|\omega_\mathrm{c} - \omega_\mathrm{vib}| \gg \Omega_R`

**Observed in simulations:**

Cavity effects can persist at off-resonant frequencies, exhibiting:

- Non-monotonic frequency dependence
- Effects that taper off gradually with detuning
- Different timescales for relaxation at different frequencies

**Physical Interpretation:**

Off-resonance cavity coupling can still affect dynamics through:

1. **Non-resonant energy transfer**: Even without polariton formation, cavity provides an energy reservoir
2. **Virtual excitations**: Off-resonant coupling creates virtual excitations that modify effective interactions
3. **Collective dipole effects**: Self-energy term :math:`\sim \lambda^2 D^2` persists regardless of resonance

Non-Monotonic Effects
----------------------

**Aging and Relaxation:**

In nonthermal aging experiments, the normalized structural relaxation time :math:`\tilde{\tau}_\mathrm{s}(\omega_\mathrm{c}, t_\mathrm{w})` can show:

**vs Frequency:**

- Not necessarily maximum at resonance
- Multiple local maxima possible
- Frequency-dependent aging mechanisms

**vs Waiting Time:**

- Non-monotonic evolution: initial slowdown, then recovery
- Maximum slowdown at intermediate :math:`t_\mathrm{w}`
- Frequency-dependent memory retention

**Example Observations:**

- Low frequencies: Longer memory, reaching equilibrium at longer waiting times
- High frequencies: Faster initial response, quicker equilibration
- Resonant frequency: Strongest effects but complex temporal behavior

Physical Mechanisms
-------------------

**Resonant Case** (:math:`\omega_\mathrm{c} \approx \omega_\mathrm{vib}`):

- Polariton formation
- Coherent Rabi oscillations
- Strong modification of vibrational density of states
- Energy localized in polariton modes

**Off-Resonant Case** (:math:`|\omega_\mathrm{c} - \omega_\mathrm{vib}| > \Omega_R`):

- No polariton formation
- Virtual processes dominate
- Weak but persistent modification of dynamics
- Energy distributed across modes

**Intermediate Case**:

- Partial hybridization
- Mixed character states
- Competing timescales
- Rich dynamical behavior

Implications for Experiments
-----------------------------

**Designing Cavity Experiments:**

1. **Resonance tuning**: Match :math:`\omega_\mathrm{c}` to target vibrational mode for maximum effect
2. **Detuning studies**: Scan :math:`\omega_\mathrm{c}` to isolate resonant vs non-resonant contributions
3. **Multi-mode systems**: Consider effects from multiple vibrational modes

**Interpreting Results:**

- Presence of cavity effects :math:`\neq` polariton formation
- Off-resonance effects can be significant
- Frequency dependence reveals underlying mechanisms
- Temporal behavior depends non-trivially on :math:`\omega_\mathrm{c}` and :math:`\lambda`

**Practical Guidelines:**

- For strong polariton effects: :math:`|\omega_\mathrm{c} - \omega_\mathrm{vib}| < \Omega_R/2`
- For exploring non-resonant effects: Scan :math:`\omega_\mathrm{c}` from :math:`0.5\omega_\mathrm{vib}` to :math:`2\omega_\mathrm{vib}`
- For disentangling mechanisms: Compare resonant and off-resonant responses

Computational Considerations
-----------------------------

**Frequency Scans:**

When performing :math:`\omega_\mathrm{c}` scans:

- Keep :math:`\lambda` fixed to isolate frequency effects
- Use same equilibration protocol for all frequencies
- Monitor both static and dynamic observables
- Perform statistical analysis across multiple realizations

**Convergence:**

- Off-resonance effects may be weaker but require good statistics
- Longer simulation times for subtle effects
- More trajectories for reliable error bars

Next Steps
==========

**For Users:**

See :doc:`../part1_application/running_simulations` for simulating strong coupling regimes and frequency scans.

**Related Theory:**

- :doc:`cavity_forces` for the coupling Hamiltonian
- :doc:`observables` for computing dynamic response
- :doc:`../part3_advanced/fdr_temperature` for non-equilibrium diagnostics

**Advanced Topics:**

- :doc:`../part3_advanced/correlation_analysis` for analyzing frequency-dependent dynamics
- :doc:`time_varying_coupling` for time-dependent frequency modulation

