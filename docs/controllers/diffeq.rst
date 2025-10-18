DiffEq Temperature Controller
=============================

Overview
--------

The **DiffEq Temperature Controller** provides first-order proportional-integral (PI) control with automatic gain tuning. It's ideal for initial temperature stabilization before switching to more sophisticated controllers like LQR.

Key Features
------------

* **Automatic Gain Tuning**: Self-tunes based on system dynamics
* **Robust Stability**: First-order response prevents oscillations
* **Smooth Transitions**: Bumpless transfer when switching controllers
* **Fast Response**: Quickly brings system to target temperature

Theory
------

Control Law
~~~~~~~~~~~

The controller implements a first-order differential equation:

.. math::

    \\tau \\frac{dT_{bath}}{dt} = k_P e + k_I \\int e dt

where:

* :math:`\\tau`: Time constant (ps)
* :math:`e = T_{target} - T_{measured}`: Temperature error
* :math:`k_P, k_I`: Proportional and integral gains (auto-tuned)

Auto-Tuning
~~~~~~~~~~~

The controller automatically computes gains based on:

1. System time constant :math:`\\tau`
2. Desired closed-loop bandwidth
3. Stability margins

Usage
-----

Basic Usage
~~~~~~~~~~~

Example for initial temperature control::

    python3 18_unified_cavity_dynamics.py \\
      --enable-diffeq-controller \\
      --diffeq-temperature-method lj_coulombic \\
      --diffeq-time-constant 0.1 \\
      --diffeq-turn-on-time 10.0 \\
      --diffeq-turn-off-time 510.0 \\
      --target-temp 100.0

Key Parameters
~~~~~~~~~~~~~~

* ``--diffeq-temperature-method``: Which temperature to control (``lj_coulombic``, ``kinetic``, ``harmonic``)
* ``--diffeq-time-constant``: Response time constant in ps (smaller = faster response)
* ``--diffeq-turn-on-time``: When to activate controller (ps)
* ``--diffeq-turn-off-time``: When to deactivate controller (ps)
* ``--diffeq-output-file``: CSV file for logging controller state

API Reference
-------------

.. autoclass:: hoomd.cavitymd.controllers.DiffEqController
   :members:
   :undoc-members:
   :show-inheritance:

