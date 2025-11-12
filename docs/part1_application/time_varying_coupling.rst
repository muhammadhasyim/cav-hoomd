=======================
Time-Varying Coupling
=======================

Comprehensive guide to dynamic coupling protocols using all available variant types.

.. contents:: In this Guide
   :local:
   :depth: 2

.. note::

   For theoretical background on time-varying coupling, see :doc:`../part2_theory/time_varying_coupling`.

Overview
========

Cavity HOOMD provides a complete **variant system** for time-dependent coupling protocols, enabling sophisticated non-equilibrium experiments and realistic modeling of pump-probe dynamics.

**Available Variant Types:**

1. **Basic Variants**: ConstantVariant, StepVariant
2. **Periodic Variants**: PeriodicVariant, SquareWaveVariant
3. **Decay Variants**: ExponentialDecayVariant, DecayingSquareWaveVariant
4. **Adaptive Variants**: AdaptiveSquareWaveVariant, ExponentialWaveVariant
5. **Scaling Variants**: LambdaScaledVariant

Why Use Time-Varying Coupling?
-------------------------------

- Model pump-probe experiments
- Study non-equilibrium relaxation dynamics
- Investigate sudden coupling activation/deactivation
- Analyze energy redistribution and thermalization
- Simulate realistic experimental protocols
- Control system temperature via coupling modulation

Basic Variants
==============

ConstantVariant
---------------

**Fixed coupling throughout simulation:**

.. code-block:: python

   from cavitymd.variants import ConstantVariant
   from cavitymd.forces import CavityForce
   
   # Constant coupling
   coupling = ConstantVariant(value=0.001)
   
   cavity_force = CavityForce(
       kvector=[0, 0, 1],
       couplstr=coupling,
       omegac=omega_c,
       phmass=1.0
   )

**Use when**: Standard equilibrium simulations with fixed coupling.

StepVariant
-----------

**Instantaneous switch at specified time:**

.. code-block:: python

   from cavitymd.variants import StepVariant
   
   # Switch coupling at t = 10 ps
   coupling = StepVariant(
       target_value=0.001,
       switch_time_ps=10.0,
       time_tracker=time_tracker
   )

**Protocol:**

.. math::

   g(t) = \begin{cases}
   0 & t < t_{\text{switch}} \\
   g_0 & t \geq t_{\text{switch}}
   \end{cases}

**Command-line usage:**

.. code-block:: bash

   python examples/05_advanced_run.py \
       --coupling 1e-3 \
       --switch-time 10.0 \
       --runtime 500

**Applications**:
- Pump-probe simulations
- Sudden quench experiments
- Coupling activation studies

Periodic Variants
=================

PeriodicVariant
---------------

**Sinusoidal modulation:**

.. code-block:: python

   from cavitymd.variants import PeriodicVariant
   
   # Sine wave modulation
   coupling = PeriodicVariant(
       amplitude=0.001,
       frequency_hz=1e12,      # 1 THz
       phase=0.0,              # Initial phase
       offset=0.0005,          # DC offset
       time_tracker=time_tracker
   )

**Protocol:**

.. math::

   g(t) = A \sin(2\pi f t + \phi) + g_{\text{offset}}

**Parameters**:
- ``amplitude``: Oscillation amplitude
- ``frequency_hz``: Modulation frequency (Hz)
- ``phase``: Initial phase (radians)
- ``offset``: DC offset (mean coupling)

**Applications**:
- AC-driven cavity coupling
- Frequency response analysis
- Floquet engineering

**Example - Temperature control via modulation:**

.. code-block:: python

   # Modulate coupling to control heating/cooling
   coupling = PeriodicVariant(
       amplitude=0.0005,
       frequency_hz=5e11,  # 0.5 THz
       offset=0.001,
       time_tracker=time_tracker
   )

SquareWaveVariant
-----------------

**Square wave modulation with controllable duty cycle:**

.. code-block:: python

   from cavitymd.variants import SquareWaveVariant
   
   # Square wave coupling
   coupling = SquareWaveVariant(
       amplitude=0.001,
       frequency_hz=1e12,      # 1 THz period
       duty_cycle=0.5,         # 50% on, 50% off
       phase_offset=0.0,
       time_tracker=time_tracker
   )

**Protocol:**

.. math::

   g(t) = \begin{cases}
   A & \text{if } (t \mod T) < T \cdot d \\
   0 & \text{otherwise}
   \end{cases}

where :math:`T = 1/f` and :math:`d` is duty cycle.

**Parameters**:
- ``amplitude``: Coupling when "on"
- ``frequency_hz``: Switching frequency
- ``duty_cycle``: Fraction of period with coupling ON (0-1)
- ``phase_offset``: Phase shift

**Applications**:
- Pulsed coupling experiments
- Intermittent cavity exposure
- Duty-cycle-dependent heating studies

**Example - Pulsed protocol:**

.. code-block:: python

   # 10% duty cycle: short pulses
   coupling = SquareWaveVariant(
       amplitude=0.01,      # Strong when on
       frequency_hz=1e12,
       duty_cycle=0.1,      # On 10% of time
       time_tracker=time_tracker
   )

Decay Variants
==============

ExponentialDecayVariant
-----------------------

**Exponential decay from initial to final value:**

.. code-block:: python

   from cavitymd.variants import ExponentialDecayVariant
   
   # Exponentially decaying coupling
   coupling = ExponentialDecayVariant(
       initial_value=0.01,         # Strong initial coupling
       final_value=0.001,          # Weak final coupling
       decay_constant_ps=20.0,     # Time constant
       start_time_ps=10.0,         # Start decay
       time_tracker=time_tracker
   )

**Protocol:**

.. math::

   g(t) = \begin{cases}
   g_0 & t < t_0 \\
   g_f + (g_0 - g_f) e^{-(t-t_0)/\tau} & t \geq t_0
   \end{cases}

**Parameters**:
- ``initial_value``: Starting coupling
- ``final_value``: Asymptotic coupling
- ``decay_constant_ps``: Time constant τ
- ``start_time_ps``: When decay begins

**Applications**:
- Gradual coupling reduction
- Adiabatic switching
- Cooling protocols

**Example - Adiabatic turn-off:**

.. code-block:: python

   # Slowly reduce coupling to final state
   coupling = ExponentialDecayVariant(
       initial_value=0.005,
       final_value=0.0,       # Turn off completely
       decay_constant_ps=50.0,  # Slow decay
       start_time_ps=100.0,
       time_tracker=time_tracker
   )

DecayingSquareWaveVariant
--------------------------

**Square wave with exponentially decaying amplitude:**

.. code-block:: python

   from cavitymd.variants import DecayingSquareWaveVariant
   
   # Decaying pulsed coupling
   coupling = DecayingSquareWaveVariant(
       initial_amplitude=0.01,
       final_amplitude=0.001,
       frequency_hz=1e12,
       duty_cycle=0.5,
       decay_constant_ps=30.0,
       start_time_ps=10.0,
       time_tracker=time_tracker
   )

**Protocol:**

.. math::

   A(t) = A_f + (A_0 - A_f) e^{-(t-t_0)/\tau}

   g(t) = \begin{cases}
   A(t) & \text{during "on" phase} \\
   0 & \text{during "off" phase}
   \end{cases}

**Applications**:
- Decaying pulse trains
- Realistic laser pulse sequences
- Progressive weakening of modulation

Adaptive Variants
=================

AdaptiveSquareWaveVariant
--------------------------

**Square wave with temperature-dependent amplitude:**

.. code-block:: python

   from cavitymd.variants import AdaptiveSquareWaveVariant
   
   # Temperature-adaptive coupling
   coupling = AdaptiveSquareWaveVariant(
       base_amplitude=0.001,
       frequency_hz=1e12,
       duty_cycle=0.5,
       temperature_tracker=temp_tracker,
       target_temperature=100.0,
       adaptation_gain=0.1,      # Feedback strength
       time_tracker=time_tracker
   )

**Protocol:**

.. math::

   A(t) = A_{\text{base}} \times \left[1 + K \cdot \frac{T_{\text{target}} - T(t)}{T_{\text{target}}}\right]

Square wave amplitude adapts based on temperature error.

**Parameters**:
- ``base_amplitude``: Nominal amplitude
- ``target_temperature``: Desired temperature
- ``adaptation_gain``: Feedback strength K
- ``temperature_tracker``: Source of T(t)

**Applications**:
- Automatic temperature control
- Self-regulating systems
- Adaptive heating/cooling

**Example - Temperature stabilization:**

.. code-block:: python

   # Coupling increases when cold, decreases when hot
   coupling = AdaptiveSquareWaveVariant(
       base_amplitude=0.001,
       target_temperature=100.0,
       adaptation_gain=0.2,      # Strong feedback
       temperature_tracker=temp_tracker,
       time_tracker=time_tracker
   )

ExponentialWaveVariant
----------------------

**Exponential modulation patterns:**

.. code-block:: python

   from cavitymd.variants import ExponentialWaveVariant
   
   # Exponentially modulated coupling
   coupling = ExponentialWaveVariant(
       amplitude=0.001,
       growth_rate=0.1,          # Growth/decay rate
       frequency_hz=1e12,
       time_tracker=time_tracker
   )

**Protocol:**

.. math::

   g(t) = A \cdot e^{\alpha t} \cdot \sin(2\pi f t)

**Applications**:
- Growing/decaying oscillations
- Chirped pulses
- Exponentially ramped protocols

Scaling Variants
================

LambdaScaledVariant
-------------------

**Automatic scaling by cavity frequency for physical units:**

.. code-block:: python

   from cavitymd.variants import LambdaScaledVariant
   
   # Lambda coupling (dimensionless)
   lambda_variant = StepVariant(target_value=0.05, switch_time_ps=10.0, time_tracker=time_tracker)
   
   # Automatically scaled: epsilon = lambda * omega_c
   coupling = LambdaScaledVariant(
       lambda_variant=lambda_variant,
       omega_c=omega_cavity  # In atomic units
   )
   
   cavity_force = CavityForce(
       kvector=[0, 0, 1],
       couplstr=coupling,  # Uses scaled value
       omegac=omega_c,
       phmass=1.0
   )

**Purpose**: Maintain physical meaning when changing cavity frequency.

**Relation**:

.. math::

   \epsilon(t) = \lambda(t) \cdot \omega_c

where λ is dimensionless coupling strength.

Composite Variants
==================

Combine Multiple Protocols
---------------------------

**Sequential composition:**

.. code-block:: python

   from cavitymd import CompositeVariant
   from cavitymd.variants import StepVariant, ExponentialDecayVariant
   
   # Phase 1: Step up
   step_up = StepVariant(target_value=0.01, switch_time_ps=10.0, time_tracker=time_tracker)
   
   # Phase 2: Decay down
   decay = ExponentialDecayVariant(
       initial_value=0.01,
       final_value=0.001,
       decay_constant_ps=20.0,
       start_time_ps=50.0,
       time_tracker=time_tracker
   )
   
   # Combine
   coupling = CompositeVariant(
       variants=[step_up, decay],
       transition_times=[10.0, 50.0]
   )

**Weighted combination:**

.. code-block:: python

   # Superpose two frequencies
   coupling1 = PeriodicVariant(amplitude=0.0005, frequency_hz=1e12, time_tracker=time_tracker)
   coupling2 = PeriodicVariant(amplitude=0.0003, frequency_hz=2e12, time_tracker=time_tracker)
   
   coupling = CompositeVariant(
       variants=[coupling1, coupling2],
       weights=[1.0, 1.0],  # Equal weight
       operation='add'
   )

Practical Examples
==================

Example 1: Pump-Probe Simulation
---------------------------------

**Protocol**: Strong pulse, then weak coupling

.. code-block:: python

   from cavitymd.variants import DecayingSquareWaveVariant
   
   # Strong initial pulse train, decaying to weak coupling
   coupling = DecayingSquareWaveVariant(
       initial_amplitude=0.01,    # Strong pump
       final_amplitude=0.0005,    # Weak probe
       frequency_hz=5e11,
       duty_cycle=0.3,
       decay_constant_ps=10.0,
       start_time_ps=5.0,
       time_tracker=time_tracker
   )
   
   # Run simulation
   sim.run(100000)  # 100 ps

Example 2: Temperature Control
-------------------------------

**Protocol**: Adaptive square wave maintains 100K

.. code-block:: python

   from cavitymd.variants import AdaptiveSquareWaveVariant
   from cavitymd.analysis import TemperatureTracker
   
   # Setup temperature tracking
   temp_tracker = TemperatureTracker(sim, time_tracker)
   
   # Adaptive coupling for temperature control
   coupling = AdaptiveSquareWaveVariant(
       base_amplitude=0.001,
       target_temperature=100.0,
       adaptation_gain=0.15,
       frequency_hz=1e12,
       duty_cycle=0.5,
       temperature_tracker=temp_tracker,
       time_tracker=time_tracker
   )
   
   # Coupling automatically adjusts to maintain T = 100K
   sim.run(500000)

Example 3: Adiabatic Switching
-------------------------------

**Protocol**: Slowly ramp coupling on, then off

.. code-block:: python

   from cavitymd.variants import ExponentialDecayVariant
   from cavitymd import CompositeVariant
   
   # Ramp up (inverse decay)
   ramp_up = ExponentialDecayVariant(
       initial_value=0.0,
       final_value=0.005,
       decay_constant_ps=30.0,
       start_time_ps=10.0,
       time_tracker=time_tracker
   )
   
   # Hold constant
   constant = ConstantVariant(value=0.005)
   
   # Ramp down
   ramp_down = ExponentialDecayVariant(
       initial_value=0.005,
       final_value=0.0,
       decay_constant_ps=30.0,
       start_time_ps=200.0,
       time_tracker=time_tracker
   )
   
   # Sequential protocol
   coupling = CompositeVariant(
       variants=[ramp_up, constant, ramp_down],
       transition_times=[10.0, 100.0, 200.0]
   )

Best Practices
==============

Choosing the Right Variant
---------------------------

**Decision tree:**

1. **Fixed coupling?** → ConstantVariant
2. **Sudden change?** → StepVariant
3. **Periodic?** 
   - Smooth: PeriodicVariant
   - Pulsed: SquareWaveVariant
4. **Gradual change?** → ExponentialDecayVariant
5. **Need adaptation?** → AdaptiveSquareWaveVariant
6. **Complex protocol?** → CompositeVariant

Time Scale Considerations
--------------------------

**Match time scales to physics:**

- **Step changes**: Instantaneous compared to molecular timescales
- **Periodic modulation**: Match molecular vibrational periods (0.01-1 ps)
- **Exponential decay**: Choose τ >> molecular relaxation time for adiabatic
- **Adaptive feedback**: Update interval ~ 0.1-1 ps

**Example:**

.. code-block:: python

   # For molecular vibration ~20 cm⁻¹ (~1.5 ps period)
   
   # Too fast (non-adiabatic)
   decay = ExponentialDecayVariant(..., decay_constant_ps=0.1)
   
   # Good (adiabatic)
   decay = ExponentialDecayVariant(..., decay_constant_ps=15.0)

Analysis and Visualization
---------------------------

**Plot coupling vs time:**

.. code-block:: python

   import matplotlib.pyplot as plt
   
   # Extract coupling history
   times = []
   couplings = []
   
   for timestep in range(sim.timestep, final_timestep):
       t_ps = timestep * sim.dt * 0.001  # Convert to ps
       g_t = coupling_variant.compute(t_ps)
       times.append(t_ps)
       couplings.append(g_t)
   
   # Plot
   plt.plot(times, couplings)
   plt.xlabel('Time (ps)')
   plt.ylabel('Coupling Strength (a.u.)')
   plt.title('Time-Varying Coupling Protocol')

**Monitor energy during protocol:**

.. code-block:: python

   # Load energy tracker data
   data = np.loadtxt('energy_tracker.txt')
   time = data[:, 0]
   total_energy = data[:, -1]
   
   # Plot energy conservation
   plt.plot(time, total_energy - total_energy[0])
   plt.xlabel('Time (ps)')
   plt.ylabel('$\\Delta E$ (a.u.)')
   plt.axvline(x=10.0, color='r', linestyle='--', label='Switch')

Troubleshooting
===============

Common Issues
-------------

**1. Energy not conserved:**

- Check if thermostat is interfering
- Verify timestep is small enough
- Ensure variant is properly configured

**2. Unphysical coupling values:**

- Check amplitude/target_value is reasonable
- Verify units (atomic units for coupling)
- Monitor for numerical overflow

**3. Adaptive variant not converging:**

- Reduce adaptation_gain
- Increase update_interval
- Check temperature_tracker is working

**4. Composite variant conflicts:**

- Verify transition times are ordered
- Check variants don't contradict
- Test each variant individually first

Next Steps
==========

- :doc:`running_simulations` for basic simulation setup
- :doc:`analysis_tools` for analyzing results
- :doc:`../part2_theory/time_varying_coupling` for theoretical background
- :doc:`../part3_advanced/controllers` for advanced control strategies
- :doc:`../api/index` for complete API reference
