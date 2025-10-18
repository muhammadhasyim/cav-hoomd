LQR Temperature Controller
==========================

Overview
--------

The **LQR Temperature Controller** implements optimal state-feedback control with integral action for temperature regulation in cavity molecular dynamics simulations. It provides **zero steady-state error** and **minimal variance** for coupled thermal subsystems using Linear-Quadratic Regulator (LQR) theory and Kalman filtering.

Key Features
------------

* **Zero Drift Guarantee**: Integral action ensures zero mean tracking error
* **Optimal Control**: LQR minimizes variance while respecting actuator limits  
* **Robust Estimation**: Kalman filter handles noisy temperature measurements  
* **In-Situ System Identification**: Automatic parameter identification via cold quench  
* **Coupled System Control**: Handles signal, hot, and bath temperature interactions

Theory
------

State-Space Model
~~~~~~~~~~~~~~~~~

The controller models three coupled temperature deviations from target::

    x_s = T_signal - T_target    (regulated output, e.g., LJ+Coulombic)
    x_h = T_hot - T_target        (disturbance, e.g., harmonic)
    x_c = T_bath - T_target       (actuator)
    ξ = ∫ x_s dt                  (integral state for zero drift)

**Discrete-time dynamics:**

.. math::

    x[k+1] = A x[k] + B u[k] + w[k]
    
    y[k] = C x[k] + v[k]

where:

* ``A, B``: System matrices (identified from cold quench experiments)
* ``w, v``: Process and measurement noise
* ``u``: Control input (bath temperature deviation)

Control Law
~~~~~~~~~~~

.. math::

    u[k] = -K \hat{x}[k] - k_I \hat{\xi}[k]

where:

* ``K``: LQR state feedback gain
* ``k_I``: Integral gain
* :math:`\hat{x}, \hat{\xi}`: Kalman filter state estimates

Guarantees
~~~~~~~~~~

1. **Zero DC Gain**: Integral action eliminates steady-state error to deterministic disturbances
2. **Optimal Variance**: LQR minimizes :math:`J = \Sigma (x^T Q x + u^T R u)` subject to dynamics
3. **Bounded Control**: Respects T_min and T_max saturation limits

Usage
-----

Basic Usage
~~~~~~~~~~~

Example with DiffEqController followed by LQR::

    python3 18_unified_cavity_dynamics.py \\
      # Phase 1: DiffEq controller (0-510 ps)
      --enable-diffeq-controller \\
      --diffeq-temperature-method lj_coulombic \\
      --diffeq-time-constant 0.1 \\
      --diffeq-turn-on-time 10.0 \\
      --diffeq-turn-off-time 510.0 \\
      \\
      # Phase 2: LQR controller (510+ ps)
      --enable-lqr-controller \\
      --lqr-signal-method lj_coulombic \\
      --lqr-hot-method harmonic_equipartition \\
      --lqr-turn-on-time 510.0 \\
      --lqr-dynamic-target \\
      --lqr-system-id-mode step \\
      --lqr-system-id-temp 5.0 \\
      --lqr-system-id-duration 50.0

Key Parameters
~~~~~~~~~~~~~~

System Identification
^^^^^^^^^^^^^^^^^^^^^

* ``--lqr-system-id-mode``: ``step`` (cold quench), ``load`` (use saved params)
* ``--lqr-system-id-temp``: Temperature change for identification (K)
* ``--lqr-system-id-duration``: Duration of identification phase (ps)
* ``--lqr-system-id-file``: File to load/save system parameters

Control Tuning
^^^^^^^^^^^^^^

* ``--lqr-weight-signal``: Weight on signal temperature error (default: 100.0)
* ``--lqr-weight-hot``: Weight on hot temperature error (default: 1.0)
* ``--lqr-control-effort``: Weight on control action (default: 0.1)
* ``--lqr-weight-integral``: Weight on integral error (default: 10.0)

Temperature Methods
^^^^^^^^^^^^^^^^^^^

* ``--lqr-signal-method``: ``lj_coulombic``, ``kinetic``, ``harmonic_equipartition``
* ``--lqr-hot-method``: ``harmonic_equipartition``, ``kinetic``, ``lj_coulombic``

Adaptive Features
-----------------

The LQR controller supports online parameter adaptation using Extended Kalman Filtering:

* ``--lqr-use-ekf-adaptation``: Enable recursive parameter estimation
* ``--lqr-ekf-update-interval``: Update every N control steps (default: 10)
* ``--lqr-adaptive-lqr``: Redesign LQR gains when parameters change significantly
* ``--lqr-adaptive-lqr-threshold``: Threshold for parameter change (default: 0.05)

API Reference
-------------

.. autoclass:: hoomd.cavitymd.controllers.LQRTemperatureController
   :members:
   :undoc-members:
   :show-inheritance:

