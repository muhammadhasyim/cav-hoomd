Temperature Controllers
======================

The CavityMD plugin provides several temperature controllers for precise thermal regulation during simulations.

.. toctree::
   :maxdepth: 2
   
   lqr
   diffeq

Overview
--------

**DiffEqController**
   Simple PI controller with automatic gain tuning. Best for initial temperature stabilization.

**LQRTemperatureController**
   Optimal state-feedback controller with Kalman filtering. Provides zero steady-state error for coupled thermal systems.

Typical Workflow
----------------

1. **Initial Stabilization** (0-500 ps): Use DiffEqController to quickly reach target temperature
2. **Precise Control** (500+ ps): Switch to LQRTemperatureController for minimal variance tracking
3. **System Identification**: LQR automatically identifies system dynamics during transition

Example
-------

Complete temperature control workflow::

    python3 18_unified_cavity_dynamics.py \\
      # Initial stabilization with DiffEq
      --enable-diffeq-controller \\
      --diffeq-temperature-method lj_coulombic \\
      --diffeq-time-constant 0.1 \\
      --diffeq-turn-on-time 10.0 \\
      --diffeq-turn-off-time 510.0 \\
      \\
      # Precise control with LQR
      --enable-lqr-controller \\
      --lqr-signal-method lj_coulombic \\
      --lqr-hot-method harmonic_equipartition \\
      --lqr-turn-on-time 510.0 \\
      --lqr-system-id-mode step \\
      --target-temp 100.0

