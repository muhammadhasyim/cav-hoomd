=================
FDR Temperature
=================

Fluctuation-Dissipation Ratio (FDR) temperature measurement for mode-specific effective temperatures in non-equilibrium cavity-coupled systems.

.. contents:: In this Section
   :local:
   :depth: 2

Overview
========

The **Fluctuation-Dissipation Ratio (FDR)** provides a physics-based method for measuring effective temperatures in systems that violate thermal equilibrium, such as cavity-coupled molecular systems under strong coupling.

**Key Concept:**

In thermal equilibrium, the fluctuation-dissipation theorem relates the autocorrelation function C(t) to the linear response function χ(t):

.. math::

   C(\omega) = \frac{2k_BT}{\omega} \chi''(\omega)

For non-equilibrium systems, this relation is violated, and we can define an **effective temperature**:

.. math::

   T_{\text{eff}}(\omega_0,t) = \frac{\omega_0}{2k_B} \times \frac{S_{AA}(\omega_0,t)}{\chi''(\omega_0,t)}

Where:

- :math:`S_{AA}(\omega_0,t)`: Power spectral density (fluctuations)
- :math:`\chi''(\omega_0,t)`: Imaginary susceptibility (dissipation)
- :math:`\omega_0`: Probe frequency
- :math:`k_B`: Boltzmann constant

**Physical Interpretation:**

- :math:`T_{\text{eff}} = T`: System in thermal equilibrium
- :math:`T_{\text{eff}} \neq T`: Non-equilibrium state
- Mode-specific: Different frequencies can have different effective temperatures

Implementation
==============

FDR Temperature Estimator
--------------------------

**Streaming algorithm for real-time effective temperature:**

.. code-block:: python

   from cavitymd import FDRTemperatureEstimator
   
   # Create estimator for dipole moment at cavity frequency
   fdr_estimator = FDRTemperatureEstimator(
       observable_name='dipole',
       probe_frequency_cm=2000.0,  # Cavity frequency
       time_tracker=time_tracker,
       update_interval_ps=0.1,
       window_size_ps=10.0,         # Sliding window
       perturbation_amplitude=0.01  # External field amplitude
   )
   
   # Computes T_eff(ω_cavity, t) in real-time

**Algorithm:**

1. **Measure fluctuations**: Compute :math:`S_{\mu\mu}(\omega_0,t)` from equilibrium trajectory
2. **Measure response**: Apply perturbation, measure :math:`\chi''_{\mu}(\omega_0,t)`
3. **Compute FDR**: :math:`T_{\text{eff}} = \frac{\omega_0}{2k_B} \frac{S_{\mu\mu}}{\chi''_{\mu}}`

FDR Integration with Simulations
---------------------------------

**Seamless integration with CavityMDSimulation:**

.. code-block:: python

   from cavitymd import FDRIntegration
   
   # Add FDR analysis to simulation
   fdr_integration = FDRIntegration(
       simulation=sim,
       observable='dipole',
       probe_frequency_cm=sim.freq,  # Use cavity frequency
       enable_equilibrium_phase=True,
       equilibrium_duration_ps=100.0,
       enable_response_phase=True,
       response_duration_ps=100.0,
       perturbation_field_strength=0.01
   )
   
   # Automatically handles:
   # 1. Equilibrium trajectory collection
   # 2. Perturbation application
   # 3. Response measurement
   # 4. FDR calculation

Complete FDR Workflow
----------------------

**End-to-end analysis pipeline:**

.. code-block:: python

   from cavitymd import FDRWorkflow
   
   # Complete workflow
   workflow = FDRWorkflow(
       job_dir='fdr_analysis',
       temperature=100.0,
       coupling=0.001,
       frequency=2000.0,
       
       # Phase 1: Equilibrium
       equilibrium_duration_ps=200.0,
       
       # Phase 2: Response
       response_duration_ps=200.0,
       perturbation_amplitude=0.01,
       
       # Analysis
       output_file='fdr_results.h5'
   )
   
   # Run complete workflow
   workflow.run()
   
   # Access results
   results = workflow.get_results()
   print(f"Effective temperature: {results['T_eff']:.2f} K")
   print(f"FDR violation: {results['fdr_ratio']:.3f}")

Mathematical Details
====================

Autocorrelation Function
-------------------------

**Equilibrium dipole autocorrelation:**

.. math::

   C_{\mu}(t) = \langle \mu(t) \cdot \mu(0) \rangle

**Power spectral density:**

.. math::

   S_{\mu\mu}(\omega) = \int_{-\infty}^{\infty} C_{\mu}(t) e^{i\omega t} dt

In practice, computed via FFT of finite-time correlation.

Linear Response Function
-------------------------

**Apply external field** :math:`E(t) = E_0 \cos(\omega_0 t)`:

.. math::

   \mu_{\text{induced}}(t) = \chi(\omega_0) E(t)

**Extract imaginary part:**

.. math::

   \chi''(\omega_0) = \frac{1}{E_0} \text{Im}[\mu(\omega_0)]

Effective Temperature
---------------------

**FDR-based effective temperature:**

.. math::

   T_{\text{eff}}(\omega_0) = \frac{\omega_0 S_{\mu\mu}(\omega_0)}{2k_B \chi''_{\mu}(\omega_0)}

**Interpretation:**

- Quantifies energy distribution at specific frequency
- Different from kinetic temperature
- Mode-specific: vibrational vs translational

Practical Examples
==================

Example 1: Cavity Frequency FDR
--------------------------------

**Measure effective temperature at cavity resonance:**

.. code-block:: python

   from cavitymd import FDRTemperatureEstimator
   from cavitymd.analysis import DipoleMomentFDRTracker
   
   # Setup FDR tracker
   fdr_tracker = DipoleMomentFDRTracker(
       simulation=sim,
       time_tracker=time_tracker,
       probe_frequency_cm=2000.0,  # Cavity frequency
       reference_interval_ps=10.0,
       output_period_ps=1.0
   )
   
   # Add to simulation
   sim.operations.updaters.append(hoomd.update.CustomUpdater(
       action=fdr_tracker,
       trigger=hoomd.trigger.Periodic(1)
   ))
   
   # Run and monitor
   sim.run(100000)
   
   # Get effective temperature history
   T_eff_history = fdr_tracker.get_temperature_history()

Example 2: Multi-Frequency FDR Scan
------------------------------------

**Scan multiple frequencies:**

.. code-block:: python

   import numpy as np
   
   # Frequency range
   frequencies = np.linspace(1500, 2500, 20)  # cm^-1
   
   T_eff_vs_freq = []
   
   for freq in frequencies:
       # Create FDR estimator for this frequency
       fdr = FDRTemperatureEstimator(
           observable_name='dipole',
           probe_frequency_cm=freq,
           time_tracker=time_tracker,
           window_size_ps=50.0
       )
       
       # Run analysis
       result = fdr.analyze(equilibrium_traj, response_traj)
       T_eff_vs_freq.append(result['T_eff'])
   
   # Plot frequency-dependent effective temperature
   plt.plot(frequencies, T_eff_vs_freq)
   plt.xlabel('Frequency (cm⁻¹)')
   plt.ylabel('Effective Temperature (K)')

Example 3: Time-Resolved FDR
-----------------------------

**Monitor FDR evolution during non-equilibrium process:**

.. code-block:: python

   from cavitymd import FDRIntegration
   
   # Setup with short windows for time resolution
   fdr = FDRIntegration(
       simulation=sim,
       observable='dipole',
       probe_frequency_cm=2000.0,
       window_size_ps=5.0,      # Short window
       update_interval_ps=0.5   # Frequent updates
   )
   
   # Run simulation with coupling switch
   sim.run(initial_equilibration_steps)
   
   # Switch coupling ON
   cavity_force.couplstr = 0.001
   
   # Monitor T_eff evolution
   sim.run(observation_steps)
   
   # Get time-resolved effective temperature
   times, T_eff = fdr.get_time_series()
   
   # Observe relaxation to new equilibrium
   plt.plot(times, T_eff)
   plt.axhline(y=100.0, color='r', linestyle='--', label='Bath T')
   plt.xlabel('Time (ps)')
   plt.ylabel('T_eff (K)')

Physical Interpretation
========================

Equilibrium vs Non-Equilibrium
-------------------------------

**Equilibrium State:**

.. math::

   T_{\text{eff}}(\omega) = T_{\text{bath}} \quad \forall \omega

All modes have same effective temperature equal to bath temperature.

**Non-Equilibrium State:**

.. math::

   T_{\text{eff}}(\omega) \neq T_{\text{bath}}

Different modes have different effective temperatures:

- High-frequency modes may be "hotter" (excess energy)
- Low-frequency modes may be "colder" (energy deficit)

Mode-Selective Effects
-----------------------

**In cavity-coupled systems:**

- **On-resonance mode** (:math:`\omega \approx \omega_{\text{cavity}}`): Often shows :math:`T_{\text{eff}} > T_{\text{bath}}`
- **Off-resonance modes**: Closer to :math:`T_{\text{bath}}`

**Physical mechanism:**

Cavity coupling preferentially excites resonant vibrational modes, creating non-thermal energy distribution.

Diagnostic Tool
---------------

**FDR as diagnostic:**

1. **Equilibration check**: :math:`T_{\text{eff}} \approx T_{\text{bath}}` indicates equilibrium
2. **Strong coupling effects**: Large :math:`|T_{\text{eff}} - T_{\text{bath}}|` indicates strong non-equilibrium
3. **Mode selectivity**: Frequency-dependent :math:`T_{\text{eff}}(\omega)` shows which modes are affected

Best Practices
==============

Measurement Guidelines
----------------------

**1. Sufficient statistics:**

.. code-block:: python

   # Long equilibrium phase for good statistics
   fdr = FDRIntegration(
       ...,
       equilibrium_duration_ps=200.0,  # At least 100-200 ps
       window_size_ps=50.0             # Multiple correlation times
   )

**2. Appropriate perturbation:**

.. code-block:: python

   # Small perturbation (linear response regime)
   fdr = FDRIntegration(
       ...,
       perturbation_amplitude=0.01  # 1% of characteristic field
   )

**3. Frequency resolution:**

.. code-block:: python

   # Match to system dynamics
   fdr = FDRIntegration(
       ...,
       probe_frequency_cm=2000.0,  # Near important resonances
       frequency_resolution=10.0    # Sufficient resolution
   )

Analysis Recommendations
------------------------

**Window size selection:**

- Too small: Poor statistics, noisy :math:`T_{\text{eff}}`
- Too large: Poor time resolution, misses dynamics
- **Recommended**: 5-10 correlation times

**Update frequency:**

- Balance between time resolution and computational cost
- **Recommended**: update_interval_ps = 0.1-1.0 ps

**Perturbation strength:**

- Must be small (linear response)
- But large enough for good signal-to-noise
- **Recommended**: 0.5-2% of characteristic energy scale

Troubleshooting
===============

Common Issues
-------------

**1. Noisy T_eff measurements:**

**Causes**: Insufficient statistics, too short windows

**Solutions**:
- Increase window_size_ps
- Increase equilibrium_duration_ps
- Average over multiple realizations

**2. T_eff far from expected:**

**Causes**: System not equilibrated, perturbation too large

**Solutions**:
- Longer equilibration before measurement
- Reduce perturbation_amplitude
- Check for systematic drifts

**3. Negative T_eff:**

**Causes**: Numerical issues, population inversion

**Solutions**:
- Check calculation of :math:`\chi''` (must be positive)
- Verify FFT settings
- May indicate true population inversion (physical)

Validation
----------

**Check against known cases:**

.. code-block:: python

   # Test 1: Equilibrium system should give T_eff = T_bath
   # Run without cavity coupling
   cavity_force.couplstr = 0.0
   fdr_result = fdr.analyze(...)
   assert abs(fdr_result['T_eff'] - T_bath) < 5.0  # Within 5K
   
   # Test 2: FDT should hold at equilibrium
   fdr_ratio = fdr_result['S_AA'] / (2 * k_B * T_bath * fdr_result['chi_imag'] / omega)
   assert abs(fdr_ratio - 1.0) < 0.1  # Within 10%

Theoretical Background
======================

Fluctuation-Dissipation Theorem
--------------------------------

**Classical FDT (equilibrium):**

.. math::

   \chi''(\omega) = \frac{\omega}{2k_BT} S_{AA}(\omega)

**Generalized FDR (non-equilibrium):**

.. math::

   \frac{S_{AA}(\omega)}{\chi''(\omega)} = \frac{2k_BT_{\text{eff}}(\omega)}{\omega}

where :math:`T_{\text{eff}}(\omega)` quantifies FDT violation.

Kubo Formula
------------

**Linear response:**

.. math::

   \chi(\omega) = \frac{i}{\hbar} \int_0^{\infty} dt \, e^{i\omega t} \langle [A(t), A(0)] \rangle

Imaginary part gives dissipation.

Cavity QED Context
------------------

In cavity-molecule systems:

- **Polariton formation**: New eigenstates with modified dynamics
- **Energy redistribution**: Non-thermal energy flow between modes
- **Collective effects**: N-molecule enhancement of coupling

FDR measurements reveal these non-equilibrium dynamics.

References
==========

**Key Papers:**

1. Cugliandolo, L. F. "The effective temperature." *J. Phys. A: Math. Theor.* 44, 483001 (2011).
2. Ciliberto, S., et al. "Experimental test of the fluctuation-dissipation theorem in a nonequilibrium steady state." *Physica A* 340, 240 (2004).
3. Herrera, F. & Spano, F. C. "Cavity-controlled chemistry in molecular ensembles." *Phys. Rev. Lett.* 116, 238301 (2016).

Next Steps
==========

- :doc:`molecular_temperatures` for molecular temperature decomposition
- :doc:`controllers` for temperature control strategies
- :doc:`../part1_application/analysis_tools` for analysis framework
- :doc:`../part2_theory/strong_coupling` for strong coupling physics
