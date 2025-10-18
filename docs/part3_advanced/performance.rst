===========
Performance
===========

Optimization and performance considerations for cavity-coupled molecular dynamics.

.. contents:: In this Section
   :local:
   :depth: 2

GPU vs CPU
==========

**Performance comparison:**

.. list-table::
   :header-rows: 1

   * - System Size
     - CPU (steps/s)
     - GPU (steps/s)
     - Speedup
   * - 512 molecules
     - 500-2000
     - 50,000-150,000
     - 50-100×
   * - 2048 molecules
     - 100-500
     - 100,000-300,000
     - 200-600×

**Recommendation:** Use GPU for N>500 molecules.

Timestep Optimization
=====================

**Adaptive timestep:**

.. code-block:: python

   from cavitymd.simulation import AdaptiveTimestepUpdater
   
   adaptive = AdaptiveTimestepUpdater(
       initial_dt=0.001,
       min_dt=0.0001,
       max_dt=0.002,
       target_energy_drift=1e-4
   )
   
   sim.operations.updaters.append(adaptive)

Balances accuracy and speed automatically.

Scaling Studies
===============

**Strong scaling:** Fixed system size, increase processors

**Weak scaling:** System size scales with processors

See HOOMD documentation for parallelization strategies.

Memory Optimization
===================

**GPU memory usage:**

.. math::

   \text{Memory} \approx N_{\text{particles}} \times (48 \text{ bytes/particle})

Modern GPUs (8+ GB) handle 10,000+ molecules easily.

Back to: :doc:`../index`
