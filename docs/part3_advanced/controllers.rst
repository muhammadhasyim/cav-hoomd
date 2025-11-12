===========
Controllers
===========

This section describes all advanced control systems available in Cavity HOOMD for temperature management, parameter optimization, and adaptive feedback.

.. contents:: In this Section
   :local:
   :depth: 2

.. note::
   
   Controllers are advanced features for specialized applications. For standard equilibrium simulations, use built-in thermostats (see :doc:`../part2_theory/thermostats`).

Overview
========

Cavity HOOMD provides a comprehensive suite of controllers for different control objectives:

1. **PID Control:** Classical three-term controller with auto-tuning
2. **Differential Equation Controller:** First-order dynamics-based control
3. **Simple Setpoint Controller:** Direct bath temperature control
4. **Adaptive MPC Controller:** Model predictive control with online adaptation
5. **Dual Feedback Controller:** Independent multi-mode temperature control
6. **Empirical Controller:** Data-driven temperature-energy relationship control
7. **Feedback Controller:** Base class for custom controller implementations

These controllers enable:

- Precise temperature tracking and regulation
- Smooth temperature ramps and protocols
- Parameter identification and adaptation
- Multi-objective control (multiple temperatures)
- Disturbance rejection
- Model-based optimal control

PID Control
===========

Classical proportional-integral-derivative controller for robust temperature regulation.

Theory
------

The **PID controller** implements the standard three-term control law:

.. math::

   u(t) = K_p \left( e(t) + \frac{1}{T_i} \int_0^t e(\tau) d\tau + T_d \frac{de}{dt} \right)

Where:

- :math:`e(t) = T_{\text{setpoint}} - T_{\text{measured}}`: Temperature error
- :math:`K_p`: Proportional gain (immediate response)
- :math:`T_i`: Integral time constant (eliminates steady-state error)
- :math:`T_d`: Derivative time constant (reduces overshoot)
- :math:`u(t) = T_{\text{bath}}`: Control output (bath temperature)

**Key Features:**

- Anti-windup protection for integral term
- Filtered derivative to reject measurement noise
- Output constraints and rate limiting
- Auto-tuning via Ziegler-Nichols method
- Multiple temperature signal options

Implementation
--------------

.. code-block:: python

   from cavitymd.controllers import PIDControl
   from cavitymd.analysis import TemperatureTracker, ElapsedTimeTracker
   
   # Setup tracking
   time_tracker = ElapsedTimeTracker(sim)
   temp_tracker = TemperatureTracker(sim, time_tracker)
   
   # Manual tuning
   pid = PIDControl(
       simulation=sim,
       time_tracker=time_tracker,
       temperature_tracker=temp_tracker,
       molecular_thermostat=bussi_thermostat,
       Kp=1.0,           # Proportional gain
       Ti=5.0,           # Integral time (ps)
       Td=0.5,           # Derivative time (ps)
       setpoint_K=100.0, # Target temperature (K)
       signal_choice='lj_coulombic',  # Control signal
       update_interval_ps=0.1
   )
   
   sim.operations.updaters.append(hoomd.update.CustomUpdater(
       action=pid,
       trigger=hoomd.trigger.Periodic(1)
   ))

Auto-Tuning
-----------

**Ziegler-Nichols method** from step response:

.. code-block:: python

   # Auto-tune based on system response
   pid = PIDControl(
       simulation=sim,
       time_tracker=time_tracker,
       temperature_tracker=temp_tracker,
       molecular_thermostat=bussi_thermostat,
       setpoint_K=100.0,
       signal_choice='lj_coulombic',
       auto_tune=True,              # Enable auto-tuning
       auto_tune_duration_ps=50.0,  # Tuning duration
       step_size=10.0               # Step size for identification
   )

The controller applies a step change, measures the response, fits a First-Order Plus Dead Time (FOPDT) model, and calculates optimal PID gains using Ziegler-Nichols rules.

**FOPDT Model:**

.. math::

   \frac{K_p e^{-\theta s}}{\tau s + 1}

Where :math:`K_p` is process gain, :math:`\tau` is time constant, :math:`\theta` is dead time.

**Ziegler-Nichols PID Gains:**

.. math::

   K_p &= 1.2 \frac{\tau}{K_p \theta} \\
   T_i &= 2\theta \\
   T_d &= 0.5\theta

Operating Modes
---------------

**1. Setpoint Tracking:**

Fixed target temperature with active control:

.. code-block:: python

   pid = PIDControl(
       ...,
       mode='setpoint',
       setpoint_K=100.0
   )

**2. Self-Loop:**

Track a measured signal (e.g., harmonic temperature follows LJ temperature):

.. code-block:: python

   pid = PIDControl(
       ...,
       mode='self_loop',
       signal_choice='lj_coulombic',       # Signal to track
       control_signal='harmonic_fictive'   # Signal to control
   )

Use Cases
---------

**Temperature stabilization:**

.. code-block:: python

   # Maintain 100K with tight regulation
   pid = PIDControl(..., setpoint_K=100.0, Kp=2.0, Ti=3.0, Td=0.3)

**Temperature ramping:**

.. code-block:: python

   # Ramp from 100K to 300K
   for t in np.linspace(0, 100, 1000):  # Over 100 ps
       T_target = 100.0 + (300.0 - 100.0) * (t / 100.0)
       pid.setpoint_K = T_target
       sim.run(100)  # Run 0.1 ps

**Multi-mode coupling:**

.. code-block:: python

   # LJ bath follows harmonic temperature
   pid = PIDControl(
       ...,
       mode='self_loop',
       signal_choice='harmonic_fictive',
       control_signal='lj_coulombic'
   )

Differential Equation Controller
=================================

Simple first-order dynamics-based temperature control.

Theory
------

The **DiffEqController** implements the differential equation:

.. math::

   \frac{dT_{\text{bath}}}{dt} = -\frac{1}{\tau}(T_{\text{bath}} - T_{\text{signal}})

Where:

- :math:`T_{\text{signal}}`: Measured temperature (e.g., LJ+Coulombic fictive)
- :math:`T_{\text{bath}}`: Bath temperature (control variable)
- :math:`\tau`: Time constant (controls response speed)

This creates a **first-order low-pass filter** behavior where the bath temperature smoothly follows the signal temperature.

Implementation
--------------

.. code-block:: python

   from cavitymd.controllers import DiffEqController
   
   diffeq = DiffEqController(
       simulation=sim,
       time_tracker=time_tracker,
       temperature_tracker=temp_tracker,
       molecular_thermostat=bussi_thermostat,
       signal_name='lj_coulombic_bath',  # Temperature to follow
       time_constant_ps=2.0,              # Response time
       turn_on_time_ps=10.0,              # Start time
       turn_off_time_ps=None              # Optional stop time
   )
   
   sim.operations.updaters.append(hoomd.update.CustomUpdater(
       action=diffeq,
       trigger=hoomd.trigger.Periodic(1)
   ))

Use Cases
---------

**Smooth temperature evolution:**

Bath temperature gently follows system temperature fluctuations, avoiding abrupt changes.

**Thermal equilibration:**

Gradually equilibrate cavity and molecular subsystems:

.. code-block:: python

   diffeq = DiffEqController(
       ...,
       signal_name='lj_coulombic_bath',
       time_constant_ps=5.0  # 5 ps response time
   )

**Delayed activation:**

Start control after initial equilibration:

.. code-block:: python

   diffeq = DiffEqController(
       ...,
       turn_on_time_ps=50.0  # Activate after 50 ps
   )

Simple Setpoint Controller
===========================

Direct bath temperature control without feedback.

Theory
------

The **SimpleSetpointController** directly sets the bath temperature:

.. math::

   T_{\text{bath}}(t) = T_{\text{setpoint}}

With optional linear ramping between values.

Implementation
--------------

.. code-block:: python

   from cavitymd.controllers import SimpleSetpointController
   
   # Fixed setpoint
   setpoint = SimpleSetpointController(
       simulation=sim,
       time_tracker=time_tracker,
       molecular_thermostat=bussi_thermostat,
       target_temperature_K=100.0,
       turn_on_time_ps=0.0
   )
   
   # With ramping
   setpoint = SimpleSetpointController(
       simulation=sim,
       time_tracker=time_tracker,
       molecular_thermostat=bussi_thermostat,
       initial_temperature_K=100.0,
       final_temperature_K=300.0,
       ramp_duration_ps=100.0,
       ramp_start_time_ps=10.0
   )

Use Cases
---------

**Fixed temperature:**

Simple constant temperature control when feedback is not needed.

**Linear temperature ramps:**

Controlled heating or cooling protocols:

.. code-block:: python

   # Cool from 300K to 100K over 50 ps
   setpoint = SimpleSetpointController(
       ...,
       initial_temperature_K=300.0,
       final_temperature_K=100.0,
       ramp_duration_ps=50.0
   )

Adaptive MPC Controller
========================

Advanced model predictive control with online parameter adaptation.

Theory
------

**Model Predictive Control (MPC)** solves an optimization problem at each time step:

.. math::

   \min_{u} \sum_{k=0}^{N-1} \left[ ||x(k) - x_{\text{ref}}||_Q^2 + ||u(k)||_R^2 + ||\Delta u(k)||_S^2 \right]

Subject to:

.. math::

   x(k+1) &= Ax(k) + Bu(k) \\
   u_{\min} &\leq u(k) \leq u_{\max} \\
   \Delta u_{\min} &\leq \Delta u(k) \leq \Delta u_{\max}

Where:

- :math:`x(k)`: State (temperature)
- :math:`u(k)`: Control input (bath temperature, coupling strength)
- :math:`N`: Prediction horizon
- :math:`Q, R, S`: Cost matrices

**Online Adaptation:**

System matrices :math:`A, B` are updated via **Recursive Least Squares (RLS)** based on observed system behavior.

Implementation
--------------

.. code-block:: python

   from cavitymd.controllers import AdaptiveMPCController
   
   mpc = AdaptiveMPCController(
       simulation=sim,
       time_tracker=time_tracker,
       temperature_tracker=temp_tracker,
       target_temperature=100.0,
       turn_on_time_ps=0.0,
       system_id_duration_ps=50.0,    # System identification phase
       update_interval_ps=0.1,
       prediction_horizon=10,          # Look-ahead steps
       Q_weight=100.0,                 # State cost
       R_lambda_weight=0.01,           # Lambda control cost
       R_temp_weight=0.1,              # Temperature control cost
       delta_lambda_weight=0.1,        # Lambda rate cost
       delta_temp_weight=0.01          # Temperature rate cost
   )

**Three-Phase Operation:**

1. **Pre-initialization (0 to turn_on_time_ps):** No control
2. **System Identification (turn_on_time_ps to turn_on_time_ps + system_id_duration_ps):** Random step changes to identify system model
3. **Control (after system ID):** MPC with continuous model updates via RLS

Use Cases
---------

**Multi-input control:**

Simultaneously control bath temperature AND coupling strength:

.. code-block:: python

   mpc = AdaptiveMPCController(
       ...,
       target_temperature=100.0,
       # Optimizes both T_bath and λ
   )

**Adaptive systems:**

Systems with time-varying dynamics benefit from online adaptation.

**Constrained optimization:**

Enforce physical constraints on controls:

.. code-block:: python

   mpc = AdaptiveMPCController(
       ...,
       lambda_min=0.0,
       lambda_max=0.1,
       T_bath_min=50.0,
       T_bath_max=300.0
   )

Dual Feedback Controller
=========================

Independent control of multiple temperature modes.

Theory
------

Controls multiple subsystems independently:

.. math::

   T_{\text{bath,i}}(t) = T_{\text{base}} + K_i \cdot (T_{\text{target,i}} - T_{\text{measured,i}})

Each mode (translational, rotational, vibrational, harmonic, etc.) can have independent:

- Target temperature
- Feedback gain
- Control activation time

Implementation
--------------

.. code-block:: python

   from cavitymd.controllers.dual_feedback import DualIndependentTemperatureFeedback
   
   dual = DualIndependentTemperatureFeedback(
       simulation=sim,
       time_tracker=time_tracker,
       temperature_tracker=temp_tracker,
       molecular_thermostat=bussi_thermostat,
       target_temperatures={
           'lj_coulombic': 100.0,
           'harmonic_fictive': 100.0
       },
       gains={
           'lj_coulombic': 0.1,
           'harmonic_fictive': 0.05
       },
       turn_on_time_ps=10.0
   )

Use Cases
---------

**Mode-selective control:**

Different control for different degrees of freedom:

.. code-block:: python

   dual = DualIndependentTemperatureFeedback(
       ...,
       target_temperatures={
           'translational': 100.0,
           'rotational': 100.0,
           'vibrational': 50.0  # Colder vibrations
       }
   )

**Staged equilibration:**

Activate controls sequentially:

.. code-block:: python

   # Control LJ first, then harmonic later
   dual = DualIndependentTemperatureFeedback(
       ...,
       turn_on_times={
           'lj_coulombic': 0.0,
           'harmonic_fictive': 50.0
       }
   )

Empirical Controller
====================

Data-driven control using empirical temperature-energy relationships.

Theory
------

Uses measured :math:`T_{\text{eff}}(E)` relationship to invert and find required energy for target temperature:

.. math::

   E_{\text{target}} = T_{\text{eff}}^{-1}(T_{\text{desired}})

Then controls bath temperature to achieve that energy.

**Fitting Models:**

- Harmonic: :math:`E = \frac{aT}{1 + bT}`
- LJ+Coulombic: :math:`E = E_0 + \frac{AT^{3/5}}{1 + CT^{3/5}}`

Implementation
--------------

.. code-block:: python

   from cavitymd.controllers.empirical import EmpiricalTemperatureData, EmpiricalTemperatureFeedback
   
   # Load empirical data
   emp_data = EmpiricalTemperatureData(
       filename='potential_energy_vs_T.txt',
       fit_type='harmonic'  # or 'lj_coulombic'
   )
   
   # Create controller
   emp_control = EmpiricalTemperatureFeedback(
       simulation=sim,
       time_tracker=time_tracker,
       temperature_tracker=temp_tracker,
       molecular_thermostat=bussi_thermostat,
       empirical_data=emp_data,
       target_temperature_K=100.0,
       feedback_gain=0.1
   )

Use Cases
---------

**Accurate temperature targeting:**

When you have calibration data and need precise temperature.

**Non-linear systems:**

Systems where temperature-energy relationship is complex.

**Model-free control:**

Don't need analytical model, just empirical data.

Feedback Controller Base Class
===============================

Base class for implementing custom controllers.

Implementation
--------------

.. code-block:: python

   from cavitymd.controllers.feedback import FeedbackController
   
   class MyCustomController(FeedbackController):
       def __init__(self, simulation, time_tracker, temperature_tracker, **kwargs):
           super().__init__(simulation, time_tracker, temperature_tracker, **kwargs)
           # Initialize custom parameters
       
       def compute_control(self, current_time_ps):
           # Implement custom control logic
           T_measured = self.temperature_tracker.get_temperature('lj_coulombic')
           T_bath = self.compute_custom_control(T_measured)
           return T_bath

Use Cases
---------

**Custom control algorithms:**

Implement specialized control strategies not covered by standard controllers.

**Research and development:**

Test new control concepts with existing infrastructure.

Controller Comparison
=====================

Performance Comparison
----------------------

.. list-table::
   :header-rows: 1
   :widths: 15 20 20 20 25

   * - Controller
     - Complexity
     - Tuning
     - Best Performance
     - Best Use Case
   * - PID
     - Medium
     - Auto or manual
     - Excellent
     - General temperature control
   * - DiffEq
     - Low
     - Single parameter
     - Good
     - Smooth tracking
   * - Setpoint
     - Very Low
     - None
     - N/A (open loop)
     - Fixed/ramped temperatures
   * - Adaptive MPC
     - High
     - Many parameters
     - Excellent (optimal)
     - Multi-input control
   * - Dual Feedback
     - Medium
     - Gains per mode
     - Good
     - Multi-mode systems
   * - Empirical
     - Medium
     - Requires data
     - Good
     - Non-linear systems

When to Use Each
----------------

**Use PID when:**

- Need robust, proven controller
- Want auto-tuning capability
- Single temperature control
- Good performance with moderate tuning effort

**Use DiffEq when:**

- Want simplest feedback control
- Smooth temperature evolution desired
- Single time constant sufficient
- Minimal tuning preferred

**Use Setpoint when:**

- No feedback needed
- Open-loop temperature control
- Simple ramps or fixed temperatures
- Minimal computational overhead

**Use Adaptive MPC when:**

- Multi-input control (T_bath AND λ)
- System dynamics change over time
- Constraints on controls important
- Optimal performance required
- Computational cost acceptable

**Use Dual Feedback when:**

- Multiple modes need independent control
- Different temperatures for different DOFs
- Staged activation desired
- Simple proportional control sufficient

**Use Empirical when:**

- Have calibration data available
- Accurate temperature targeting critical
- System highly non-linear
- Model-free approach preferred

Troubleshooting
===============

Common Issues
-------------

**1. Controller instability (oscillations):**

**Symptoms:** Temperature oscillates or diverges

**Solutions:**

- **PID:** Reduce Kp, increase Ti, reduce Td
- **DiffEq:** Increase time constant τ
- **MPC:** Increase R weights (control cost), reduce Q weights
- **Dual:** Reduce feedback gains
- All: Increase update interval

**2. Slow convergence:**

**Symptoms:** Takes too long to reach target

**Solutions:**

- **PID:** Increase Kp, reduce Ti
- **DiffEq:** Reduce time constant τ
- **MPC:** Increase Q weights, reduce R weights
- **Dual:** Increase feedback gains
- All: Decrease update interval

**3. Steady-state error:**

**Symptoms:** Never quite reaches target

**Solutions:**

- **PID:** Reduce Ti (increase integral action)
- **DiffEq:** Check for systematic bias
- **MPC:** Verify constraints not limiting
- **Empirical:** Check calibration data accuracy
- All: Verify target is physically achievable

**4. Noisy control signal:**

**Symptoms:** Bath temperature fluctuates rapidly

**Solutions:**

- **PID:** Increase Td filter time, reduce Kp
- All: Increase update interval
- All: Add rate limiting on control output
- All: Filter measurement signal

**5. Auto-tune failure:**

**Symptoms:** PID auto-tune doesn't converge

**Solutions:**

- Increase auto_tune_duration_ps
- Adjust step_size
- Check system is at steady state before tuning
- Verify measurement signal is responsive

Best Practices
==============

1. **Start simple:** Try DiffEq or Setpoint before PID or MPC
2. **Test in isolation:** Validate controller before production runs
3. **Monitor performance:** Track error, control effort, convergence
4. **Document settings:** Record all controller parameters
5. **Compare to baseline:** Verify improvement over open-loop
6. **Use auto-tuning:** Let PID auto-tune when possible
7. **Understand dynamics:** Know system time constants
8. **Respect constraints:** Don't request impossible performance
9. **Save data:** Log controller state and measurements
10. **Validate physically:** Ensure temperatures make physical sense

Advanced Topics
===============

Cascaded Control
----------------

Combine multiple controllers for hierarchical control:

.. code-block:: python

   # Outer loop: Setpoint provides slow temperature ramp
   setpoint = SimpleSetpointController(...)
   
   # Inner loop: PID tracks setpoint with fast regulation
   pid = PIDControl(..., setpoint_source=setpoint)

Gain Scheduling
---------------

Adjust controller gains based on operating point:

.. code-block:: python

   def update_gains(T_current):
       if T_current < 100:
           pid.Kp = 2.0
       else:
           pid.Kp = 1.0

Feedforward Control
-------------------

Add feedforward term for known disturbances:

.. code-block:: python

   T_bath = pid.compute_control(T_measured) + feedforward_term(t)

Example Workflows
=================

Complete Temperature Control
-----------------------------

.. code-block:: python

   from cavitymd import CavityMDSimulation
   from cavitymd.controllers import PIDControl
   from cavitymd.analysis import TemperatureTracker, ElapsedTimeTracker
   
   # Setup simulation
   sim = CavityMDSimulation(
       job_dir='pid_control_test',
       replica=0,
       temperature=100.0,
       freq=2000,
       couplstr=0.001,
       incavity=True,
       runtime_ps=1000
   )
   
   # Setup tracking
   time_tracker = ElapsedTimeTracker(sim.sim)
   temp_tracker = TemperatureTracker(
       sim.sim, time_tracker,
       output_period_ps=0.1
   )
   
   # Setup PID controller with auto-tuning
   pid = PIDControl(
       simulation=sim,
       time_tracker=time_tracker,
       temperature_tracker=temp_tracker,
       molecular_thermostat=sim.bussi_thermostat,
       setpoint_K=100.0,
       signal_choice='lj_coulombic',
       auto_tune=True,
       auto_tune_duration_ps=50.0
   )
   
   # Add to simulation
   sim.sim.operations.updaters.append(hoomd.update.CustomUpdater(
       action=pid,
       trigger=hoomd.trigger.Periodic(1)
   ))
   
   # Run
   sim.run()

Next Sections
=============

- :doc:`fdr_temperature` for advanced temperature measurement
- :doc:`molecular_temperatures` for temperature decomposition
- :doc:`../part2_theory/thermostats` for standard thermostat theory
- :doc:`../part1_application/running_simulations` for basic simulation setup
