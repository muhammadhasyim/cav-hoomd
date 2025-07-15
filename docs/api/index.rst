=============
API Reference
=============

This section contains the complete API documentation for Cavity HOOMD, including all components for time-varying coupling, analysis, and simulation management.

Cavity MD Module (``hoomd.cavitymd``)
=====================================

Forces and Variants
-------------------

.. currentmodule:: hoomd.cavitymd.forces

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   CavityForce

.. currentmodule:: hoomd.cavitymd.variants

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   StepVariant
   ConstantVariant

The :class:`hoomd.cavitymd.forces.CavityForce` class now supports time-varying coupling through HOOMD variants, enabling dynamic coupling switching experiments.

Simulation Framework
--------------------

.. currentmodule:: hoomd.cavitymd.simulation

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   CavityMDSimulation
   AdaptiveTimestepUpdater

.. currentmodule:: hoomd.cavitymd.updaters

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   CavityParticleDisplacer

The :class:`hoomd.cavitymd.simulation.CavityMDSimulation` class provides a complete simulation framework with smart cavity particle handling, parameter validation, and comprehensive logging.

Analysis and Tracking
---------------------

.. currentmodule:: hoomd.cavitymd.analysis

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   Status
   ElapsedTimeTracker
   TimestepFormatter
   CavityModeTracker
   AutocorrelationTracker
   FieldAutocorrelationTracker
   EnergyTracker
   DipoleAutocorrelation
   PerformanceTracker

These classes provide comprehensive monitoring and analysis capabilities for cavity MD simulations, including energy conservation testing and correlation function analysis.

Utilities
---------

.. currentmodule:: hoomd.cavitymd.utils

.. autosummary::
   :toctree: _autosummary

   PhysicalConstants
   unwrap_positions
   get_slurm_info
   parse_replicas

Core utility functions for physical constants, coordinate transformations, and HPC integration.

Bussi Reservoir Module (``hoomd.bussi_reservoir``)
===================================================

Thermostats
-----------

.. currentmodule:: hoomd.bussi_reservoir

.. autoclass:: BussiReservoir
   :members:
   :inherited-members:
   :special-members: __init__
   :no-index:

The Bussi thermostat provides efficient canonical ensemble sampling with reservoir energy tracking.

Key Features and Changes
========================

Time-Varying Coupling
---------------------

**New in this version**: The plugin now supports time-varying coupling strength through HOOMD variants:

.. code-block:: python

   from hoomd.cavitymd.forces import CavityForce
   from hoomd.cavitymd.variants import StepVariant
   
   # Create step variant for coupling switching
   coupling_variant = StepVariant(
       target_value=0.001,
       switch_time_ps=1.0,
       time_tracker=time_tracker
   )
   
   # Create cavity force with time-varying coupling
   cavity_force = CavityForce(
       kvector=[0, 0, 1],
       couplstr=coupling_variant,
       omegac=0.01,
       phmass=1.0
   )

Energy Conservation Testing
-----------------------------

**Enhanced energy tracking**: Comprehensive energy component monitoring with conservation validation:

.. code-block:: python

   from hoomd.cavitymd.analysis import EnergyTracker
   
   # Create energy tracker
   energy_tracker = EnergyTracker(
       sim=sim,
       forces={'cavity': cavity_force},
       thermostats={'bussi': bussi_thermostat}
   )
   
   # Monitor energy components
   sim.operations.updaters.append(energy_tracker)

GPU Acceleration
----------------

**Optimized GPU kernels**: High-performance CUDA implementations for cavity force computation:

.. code-block:: python

   import hoomd
   
   # Create GPU device
   device = hoomd.device.GPU()
   sim = hoomd.Simulation(device=device)
   
   # Cavity force automatically uses GPU kernels
   cavity_force = CavityForce(
       kvector=[0, 0, 1],
       couplstr=0.001,
       omegac=0.01
   )

Adaptive Timestep Control
-------------------------

**Automatic timestep optimization**: Dynamic timestep adjustment for stability and performance:

.. code-block:: python

   from hoomd.cavitymd.simulation import AdaptiveTimestepUpdater
   
   # Create adaptive timestep updater
   adaptive_dt = AdaptiveTimestepUpdater(
       sim=sim,
       target_error=0.01,
       max_dt=0.01,
       min_dt=0.0001
   )
   
   sim.operations.updaters.append(adaptive_dt)

Simulation Framework
--------------------

**Complete simulation orchestration**: The :class:`CavityMDSimulation` class provides:

- Smart cavity particle handling
- Parameter validation and logging
- Comprehensive analysis setup
- SLURM integration
- Error handling and recovery

.. code-block:: python

   from hoomd.cavitymd.simulation import CavityMDSimulation
   
   # Create simulation with all features
   sim = CavityMDSimulation(
       job_dir="test_run",
       replica=1,
       freq=2000,
       couplstr=0.001,
       incavity=True,
       runtime_ps=1000,
       switch_time_ps=1.0,  # Enable time-varying coupling
       enable_energy_tracking=True,
       enable_fkt=True,
       device='GPU'
   )
   
   # Run complete simulation
   sim.run()

Analysis and Visualization
--------------------------

**Comprehensive analysis tools**: Built-in analysis components for:

- Energy conservation monitoring
- Correlation function analysis
- Cavity mode tracking
- Performance benchmarking
- Visualization utilities

.. code-block:: python

   from hoomd.cavitymd.analysis import CavityModeTracker, FieldAutocorrelationTracker
   
   # Track cavity mode properties
   cavity_tracker = CavityModeTracker(sim=sim)
   
   # Analyze density correlations
   fkt_tracker = FieldAutocorrelationTracker(
       sim=sim,
       kmag=1.0,
       num_wavevectors=50
   )

Migration Guide
===============

**From Previous Versions**

The simplified architecture requires minimal changes:

.. code-block:: python

   # Old import (still works)
   from hoomd.plugin.cavitymd import CavityForce
   
   # New recommended import
   from hoomd.cavitymd.forces import CavityForce
   
   # New time-varying coupling
   from hoomd.cavitymd.variants import StepVariant
   
   coupling_variant = StepVariant(
       target_value=0.001,
       switch_time_ps=1.0,
       time_tracker=time_tracker
   )

**New Features to Explore**

1. **Time-varying coupling**: Test dynamic coupling switching
2. **Energy conservation**: Monitor total energy during transitions
3. **GPU acceleration**: Use optimized CUDA kernels
4. **Adaptive timestep**: Automatic timestep optimization
5. **Comprehensive analysis**: Built-in tracking and visualization

Performance Optimization
========================

**GPU Performance Tips**

- Use GPU device for large systems (>1000 particles)
- Enable GPU-specific optimizations in build
- Monitor GPU memory usage for very large systems

**CPU Performance Tips**

- Use OpenMP threading for CPU builds
- Optimize timestep for your system
- Use efficient thermostat combinations

**Memory Management**

- Use appropriate output periods to control memory usage
- Monitor energy tracker memory for long simulations
- Use GSD compression for large trajectory files

This comprehensive API provides all the tools needed for advanced cavity molecular dynamics simulations, from basic force computation to sophisticated time-varying coupling experiments with full analysis capabilities. 