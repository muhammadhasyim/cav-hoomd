===========
Controllers
===========

This section describes advanced control systems available in Cavity HOOMD for temperature management and parameter optimization.

.. contents:: In this Section
   :local:
   :depth: 2

.. note::
   
   Controllers are experimental features. For standard simulations, use built-in thermostats (see :doc:`../part2_theory/thermostats`).

Overview
========

Cavity HOOMD provides three advanced controllers:

1. **LQR Temperature Controller:** Optimal state-feedback control
2. **Extended Kalman Filter (EKF):** Parameter estimation
3. **Differential Equation Controller:** PID-like control

These enable:

- Precise temperature tracking
- Parameter identification
- Smooth temperature ramps
- Advanced control protocols

LQR Temperature Controller
===========================

Theory
------

**Linear-Quadratic Regulator (LQR)** is an optimal control method that minimizes:

.. math::

   J = \int_0^\infty (\vec{x}^T Q \vec{x} + \vec{u}^T R \vec{u}) dt

Where:

- :math:`\vec{x}`: State vector (temperatures, energies)
- :math:`\vec{u}`: Control vector (thermostat settings)
- :math:`Q`: State cost matrix
- :math:`R`: Control cost matrix

**Control law:**

.. math::

   \vec{u}(t) = -K \vec{x}(t)

where K is the optimal gain matrix.

Implementation
--------------

.. code-block:: python

   from cavitymd.controllers import LQRTemperatureController
   
   lqr = LQRTemperatureController(
       simulation=sim,
       target_temperatures={'T_s': 100.0, 'T_k': 100.0, 'T_v': 100.0},
       Q_matrix=np.diag([1.0, 1.0, 1.0]),  # State costs
       R_matrix=np.diag([0.1, 0.1]),        # Control costs
       update_period_ps=0.1
   )
   
   sim.operations.updaters.append(lqr)

**Parameters:**

- ``target_temperatures``: Desired T_translational, T_rotational, T_vibrational
- ``Q_matrix``: Penalizes state deviations
- ``R_matrix``: Penalizes control effort
- ``update_period_ps``: Control update frequency

Use Cases
---------

**1. Multi-mode temperature control:**

Independently control translational, rotational, and vibrational temperatures.

**2. Fast equilibration:**

Optimal controller achieves target faster than standard thermostats.

**3. Temperature ramps:**

Smoothly change temperature over time.

**Example:**

.. code-block:: python

   # Ramp from 100K to 300K over 100 ps
   for t in range(0, 100000, 100):  # Every 0.1 ps
       T_target = 100.0 + (300.0 - 100.0) * (t / 100000)
       lqr.set_target_temperatures({
           'T_s': T_target,
           'T_k': T_target,
           'T_v': T_target
       })
       sim.run(100)

Parameter Tuning
----------------

**Q matrix (state cost):**

- Larger values → tighter temperature control
- Smaller values → more relaxed control

**R matrix (control cost):**

- Larger values → gentler control (less aggressive)
- Smaller values → aggressive control (faster response)

**Trade-off:**

Balance between fast convergence (small R) and smooth control (large R).

Extended Kalman Filter
=======================

Theory
------

**Extended Kalman Filter (EKF)** estimates system state and parameters from noisy measurements.

**State estimation:**

.. math::

   \hat{\vec{x}}_{k+1} = \vec{f}(\hat{\vec{x}}_k, \vec{u}_k) + K_k(\vec{z}_k - \vec{h}(\hat{\vec{x}}_k))

Where:

- :math:`\vec{f}`: State transition function
- :math:`\vec{h}`: Measurement function
- :math:`K_k`: Kalman gain
- :math:`\vec{z}_k`: Measurements

Implementation
--------------

.. code-block:: python

   from cavitymd.controllers import ParameterEKF
   
   ekf = ParameterEKF(
       simulation=sim,
       parameters_to_estimate=['coupling_strength', 'cavity_frequency'],
       initial_guess={'coupling_strength': 1e-3, 'cavity_frequency': 2000.0},
       process_noise=1e-6,
       measurement_noise=1e-4,
       update_period_ps=1.0
   )
   
   sim.operations.updaters.append(ekf)
   
   # After simulation
   estimated_params = ekf.get_estimated_parameters()
   print(f"Estimated coupling: {estimated_params['coupling_strength']}")

Use Cases
---------

**1. Parameter identification:**

Estimate unknown system parameters from trajectory data.

**2. State estimation:**

Reconstruct full state from partial measurements.

**3. Adaptive control:**

Update control based on estimated parameters.

Convergence Diagnostics
-----------------------

**Monitor covariance:**

.. code-block:: python

   covariance = ekf.get_covariance_matrix()
   uncertainty = np.sqrt(np.diag(covariance))
   
   print(f"Parameter uncertainties: {uncertainty}")

**Convergence criterion:**

.. math::

   \sigma_{\text{param}} < 0.01 \times |\hat{\text{param}}|

Differential Equation Controller
=================================

Theory
------

**PID-like control** adjusts control based on error, integral, and derivative:

.. math::

   u(t) = K_p e(t) + K_i \int_0^t e(\tau) d\tau + K_d \frac{de}{dt}

Where:

- :math:`e(t) = T_{\text{target}} - T_{\text{current}}`: Error
- :math:`K_p, K_i, K_d`: Gain parameters

Implementation
--------------

.. code-block:: python

   from cavitymd.controllers import DiffEqController
   
   pid = DiffEqController(
       simulation=sim,
       target_temperature=100.0,
       Kp=1.0,   # Proportional gain
       Ki=0.1,   # Integral gain
       Kd=0.05,  # Derivative gain
       update_period_ps=0.1
   )
   
   sim.operations.updaters.append(pid)

Use Cases
---------

**1. Simple temperature control:**

Easy to implement and tune.

**2. Disturbance rejection:**

Integral term eliminates steady-state error.

**3. Smooth response:**

Derivative term reduces overshoot.

Tuning Guidelines
-----------------

**Ziegler-Nichols method:**

1. Set :math:`K_i = K_d = 0`
2. Increase :math:`K_p` until sustained oscillation
3. Use oscillation period and critical gain to set all three

**Manual tuning:**

- Start with :math:`K_p = 1.0, K_i = 0.1, K_d = 0.01`
- Increase :math:`K_p` for faster response
- Increase :math:`K_i` to eliminate steady-state error
- Increase :math:`K_d` to reduce overshoot

Controller Comparison
=====================

Performance Comparison
----------------------

.. list-table::
   :header-rows: 1
   :widths: 20 25 25 30

   * - Controller
     - Advantages
     - Disadvantages
     - Best For
   * - LQR
     - Optimal, multi-variable
     - Requires model, tuning
     - Complex control tasks
   * - EKF
     - Parameter estimation
     - Computationally expensive
     - Unknown parameters
   * - DiffEq (PID)
     - Simple, robust
     - Manual tuning needed
     - Simple setpoint tracking

When to Use Each
----------------

**Use LQR when:**

- Multiple temperatures need control
- Optimal performance required
- System model available

**Use EKF when:**

- Parameters unknown
- State estimation needed
- Adaptive control desired

**Use DiffEq (PID) when:**

- Simple temperature control
- Quick implementation needed
- Robust to model uncertainty

**Use standard thermostats when:**

- Standard equilibrium MD
- No special control requirements
- Simplicity preferred

Troubleshooting
===============

Common Issues
-------------

**1. Controller instability:**

**Symptoms:** Oscillations, divergence

**Solutions:**

- Reduce gains (K_p, K_i, K_d or Q, R matrices)
- Increase update period
- Check for numerical issues

**2. Slow convergence:**

**Symptoms:** Takes too long to reach target

**Solutions:**

- Increase proportional gain
- Reduce control cost (R matrix)
- Decrease update period

**3. Steady-state error:**

**Symptoms:** Never quite reaches target

**Solutions:**

- Add/increase integral term
- Check for systematic bias
- Verify target is achievable

Best Practices
==============

1. **Start simple:** Try standard thermostats first
2. **Test in isolation:** Validate controller before production runs
3. **Monitor performance:** Track error, control effort, convergence
4. **Document settings:** Record all controller parameters
5. **Compare to baseline:** Verify improvement over standard methods

Next Sections
=============

- :doc:`fdr_temperature` for advanced temperature measurement
- :doc:`molecular_temperatures` for temperature decomposition
- :doc:`../part2_theory/thermostats` for standard thermostat theory

