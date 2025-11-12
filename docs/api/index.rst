=============
API Reference
=============

This section contains the complete API documentation for Cavity HOOMD, including all components for time-varying coupling, analysis, control systems, and simulation management.

.. contents:: API Modules
   :local:
   :depth: 2

Cavity MD Module (``hoomd.cavitymd``)
======================================

The main cavity molecular dynamics package providing forces, variants, analysis tools, controllers, and simulation frameworks.

Forces
------

.. currentmodule:: hoomd.cavitymd.forces

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   CavityForce
   DipoleResponseForce
   PerturbationForce

The :class:`CavityForce` class implements the cavity-molecule coupling force with support for time-varying coupling through HOOMD variants. :class:`DipoleResponseForce` handles linear response measurements, and :class:`PerturbationForce` applies external perturbations for FDR analysis.

Coupling Variants
-----------------

Time-dependent coupling protocols for dynamic coupling experiments.

Basic Variants
~~~~~~~~~~~~~~

.. currentmodule:: hoomd.cavitymd.variants

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   ConstantVariant
   StepVariant

**ConstantVariant**: Fixed coupling strength throughout simulation.

**StepVariant**: Step change at specified time, enabling coupling switching experiments.

Periodic Variants
~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   PeriodicVariant
   SquareWaveVariant

**PeriodicVariant**: Sinusoidal modulation of coupling strength.

**SquareWaveVariant**: Square wave modulation with controllable duty cycle and frequency.

Decay Variants
~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   ExponentialDecayVariant
   DecayingSquareWaveVariant

**ExponentialDecayVariant**: Exponential decay from initial to target value.

**DecayingSquareWaveVariant**: Square wave with exponentially decaying amplitude.

Adaptive Variants
~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   AdaptiveSquareWaveVariant
   ExponentialWaveVariant

**AdaptiveSquareWaveVariant**: Square wave with temperature-dependent amplitude adaptation.

**ExponentialWaveVariant**: Exponential modulation patterns.

Scaling Variants
~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   LambdaScaledVariant

**LambdaScaledVariant**: Automatic scaling by cavity frequency for physical units.

Simulation Framework
--------------------

Complete simulation orchestration with comprehensive parameter management.

.. currentmodule:: hoomd.cavitymd.simulation

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   CavityMDSimulation
   AdaptiveTimestepUpdater

**CavityMDSimulation**: High-level simulation framework providing:

- Smart cavity particle handling
- Automatic thermostat configuration
- Time-varying coupling setup
- Comprehensive analysis and logging
- Parameter validation and error checking
- SLURM integration for HPC environments

**AdaptiveTimestepUpdater**: Dynamic timestep adjustment for stability and performance.

Updaters
--------

.. currentmodule:: hoomd.cavitymd.updaters

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   CavityParticleDisplacer
   HarmonicBondReset

**CavityParticleDisplacer**: Handles cavity mode particle positions and periodic boundary conditions.

**HarmonicBondReset**: One-time harmonic bond reinitialization for diatomic molecules.

Controllers
-----------

Advanced control systems for temperature management and parameter optimization.

.. currentmodule:: hoomd.cavitymd.controllers

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   DiffEqController
   SimpleSetpointController
   PIDControl
   AdaptiveMPCController

**DiffEqController**: Differential equation-based temperature control using first-order dynamics:

.. math::

   \frac{dT_{\text{bath}}}{dt} = -\frac{1}{\tau}(T_{\text{bath}} - T_{\text{signal}})

**SimpleSetpointController**: Direct bath temperature setpoint control with optional ramping.

**PIDControl**: Classical PID controller with integral windup protection and derivative filtering:

.. math::

   u(t) = K_p e(t) + K_i \int_0^t e(\tau) d\tau + K_d \frac{de}{dt}

**AdaptiveMPCController**: Model Predictive Control with online parameter adaptation via recursive least squares.

See :doc:`../part3_advanced/controllers` for detailed usage and comparison.

Analysis and Tracking
---------------------

Comprehensive monitoring and analysis capabilities for cavity MD simulations.

Timing Utilities
~~~~~~~~~~~~~~~~

.. currentmodule:: hoomd.cavitymd.analysis

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   Status
   ElapsedTimeTracker
   TimestepFormatter

**Status**: Simulation status tracking and reporting.

**ElapsedTimeTracker**: Precise elapsed time tracking for time-based output triggers.

**TimestepFormatter**: Format timestep information for display and logging.

Energy and Performance Tracking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   EnergyTracker
   PerformanceTracker
   TemperatureTracker

**EnergyTracker**: Comprehensive energy component monitoring with conservation validation:

- Molecular kinetic and potential energy
- Cavity photon energy
- Interaction energy
- Thermostat reservoir energy
- Total energy conservation testing

**PerformanceTracker**: Simulation performance metrics (TPS, ns/day, GPU utilization).

**TemperatureTracker**: Multi-mode temperature decomposition:

- Translational, rotational, vibrational temperatures
- Harmonic fictive temperature
- LJ+Coulombic fictive temperature
- Kinetic temperature

Correlation Analysis
~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   FieldAutocorrelationTracker
   AutocorrelationTracker
   CavityModeTracker

**FieldAutocorrelationTracker**: Density field autocorrelation F(k,t) for intermediate scattering function analysis.

**AutocorrelationTracker**: General-purpose autocorrelation function computation.

**CavityModeTracker**: Track cavity mode properties (amplitude, frequency, energy).

FDR Analysis
~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   DipoleMomentFDRTracker

**DipoleMomentFDRTracker**: Fluctuation-dissipation ratio analysis for dipole moment:

.. math::

   T_{\text{eff}}(\omega_0,t) = \frac{\omega_0}{2k_B} \times \frac{S_{\mu\mu}(\omega_0,t)}{\chi''_{\mu}(\omega_0,t)}

Measures violations of equilibrium fluctuation-dissipation theorem for effective temperature determination.

Data Management
---------------

HDF5-based observable output system with real-time access.

.. currentmodule:: hoomd.cavitymd.data

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   ObservableWriter
   ObservableReader

**ObservableWriter**: Unified HDF5 output with SWMR (Single Writer Multiple Reader) mode:

- Hierarchical organization by observable type
- Chunked storage with compression
- Thread-safe concurrent read access during simulation
- 10-100x smaller file sizes vs text output

**ObservableReader**: Read HDF5 observable data with convenient access patterns.

See :doc:`../part1_application/data_management` for usage examples.

FDR Temperature Analysis
------------------------

Complete workflow for fluctuation-dissipation ratio temperature measurement.

.. currentmodule:: hoomd.cavitymd

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   FDRTemperatureEstimator
   FDRIntegration
   FDRWorkflow

**FDRTemperatureEstimator**: Streaming algorithm for single-frequency effective temperature:

.. math::

   T_{\text{eff}}(\omega_0,t) = \frac{\omega_0}{2k_B} \times \frac{S_{AA}(\omega_0,t)}{\chi''(\omega_0,t)}

**FDRIntegration**: Seamless integration with CavityMDSimulation framework.

**FDRWorkflow**: Complete end-to-end FDR analysis pipeline.

See :doc:`../part3_advanced/fdr_temperature` for detailed examples.

Utilities
---------

Core utility functions for physical constants, coordinate transformations, and HPC integration.

.. currentmodule:: hoomd.cavitymd.utils

.. autosummary::
   :toctree: _autosummary

   PhysicalConstants
   unwrap_positions
   get_slurm_info
   parse_replicas

**PhysicalConstants**: Physical constants in atomic units:

- ``HARTREE_TO_CM_MINUS1``: Convert Hartree to wavenumbers
- ``BOLTZMANN_AU``: Boltzmann constant in a.u.
- ``KELVIN_TO_AU``: Convert temperature to a.u.
- ``FEMTOSEC_TO_AU``: Convert time to a.u.

**unwrap_positions**: Unwrap particle positions across periodic boundaries.

**get_slurm_info**: Extract SLURM environment information for HPC jobs.

**parse_replicas**: Parse replica indices from command-line arguments.

Validation
----------

Parameter validation and input checking.

.. currentmodule:: hoomd.cavitymd.validation

.. autosummary::
   :toctree: _autosummary

   validate_cavity_parameters
   validate_thermostat_config
   validate_coupling_variant

Input validation functions ensuring physical consistency and preventing common errors.

Experiment Management
---------------------

High-level experiment orchestration for parameter sweeps and workflows.

.. currentmodule:: hoomd.cavitymd.experiment

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   ExperimentManager
   ParameterSweep

**ExperimentManager**: Coordinate multiple simulations with shared configuration.

**ParameterSweep**: Automated parameter sweep execution with parallel job management.

Composite Variants
------------------

Combine multiple coupling protocols.

.. currentmodule:: hoomd.cavitymd

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   CompositeVariant

**CompositeVariant**: Combine multiple variants with weighted or sequential composition.

State Management
----------------

.. currentmodule:: hoomd.cavitymd.state_manager

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   StateManager

**StateManager**: Checkpoint and restart functionality for long simulations.

Bussi Reservoir Module (``hoomd.bussi_reservoir``)
===================================================

Thermostats
-----------

.. currentmodule:: hoomd.bussi_reservoir

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   BussiReservoir

**BussiReservoir**: Bussi-Parrinello thermostat with reservoir energy tracking:

- Efficient canonical ensemble sampling
- Energy exchange monitoring with reservoir
- Compatible with cavity coupling
- GPU-accelerated implementation

The Bussi thermostat provides proper canonical ensemble sampling while tracking energy exchange with the thermal reservoir, essential for energy conservation analysis in cavity-coupled systems.

Key Features Summary
====================

Time-Varying Coupling
---------------------

**Complete variant system** for dynamic coupling experiments:

.. code-block:: python

   from hoomd.cavitymd.variants import StepVariant, SquareWaveVariant
   from hoomd.cavitymd.forces import CavityForce
   
   # Step switching
   coupling_variant = StepVariant(
       target_value=0.001,
       switch_time_ps=1.0,
       time_tracker=time_tracker
   )
   
   # Or periodic modulation
   coupling_variant = SquareWaveVariant(
       amplitude=0.001,
       frequency_hz=1e12,
       duty_cycle=0.5,
       time_tracker=time_tracker
   )
   
   cavity_force = CavityForce(
       kvector=[0, 0, 1],
       couplstr=coupling_variant,
       omegac=0.01,
       phmass=1.0
   )

Advanced Control Systems
------------------------

**Multiple controller options** for different control objectives:

.. code-block:: python

   from hoomd.cavitymd.controllers import PIDControl, DiffEqController
   
   # PID for precise setpoint tracking
   pid = PIDControl(
       simulation=sim,
       target_temperature=100.0,
       Kp=1.0, Ki=0.1, Kd=0.05
   )
   
   # DiffEq for smooth temperature evolution
   diffeq = DiffEqController(
       simulation=sim,
       signal_name='lj_coulombic_bath',
       time_constant_ps=2.0
   )

Comprehensive Analysis
----------------------

**Built-in analysis tools** for monitoring and validation:

.. code-block:: python

   from hoomd.cavitymd.analysis import EnergyTracker, TemperatureTracker
   from hoomd.cavitymd.data import ObservableWriter
   
   # Energy conservation monitoring
   energy_tracker = EnergyTracker(
       sim=sim,
       forces={'cavity': cavity_force},
       thermostats={'bussi': bussi_thermostat},
       output_period_ps=0.1
   )
   
   # Temperature decomposition
   temp_tracker = TemperatureTracker(
       sim=sim,
       output_period_ps=0.1
   )
   
   # Unified HDF5 output
   writer = ObservableWriter(
       filename='observables.h5',
       observables=['energy', 'temperature', 'cavity_mode']
   )

GPU Acceleration
----------------

**Optimized CUDA kernels** for high performance:

.. code-block:: python

   import hoomd
   
   # Automatically uses GPU kernels when available
   device = hoomd.device.GPU()
   sim = hoomd.Simulation(device=device)
   
   # All cavity forces run on GPU
   cavity_force = CavityForce(...)  # GPU-accelerated

HDF5 Data Management
--------------------

**Modern data output system** with real-time access:

.. code-block:: python

   from hoomd.cavitymd.data import ObservableWriter
   
   # Single HDF5 file for all observables
   writer = ObservableWriter(
       filename='sim_data.h5',
       swmr_mode=True  # Enable live monitoring
   )
   
   # Read data while simulation runs
   from hoomd.cavitymd.data import ObservableReader
   reader = ObservableReader('sim_data.h5')
   energies = reader.get_observable('energy/total')

Migration from Previous Versions
=================================

**Module reorganization** - New import structure:

.. code-block:: python

   # Old imports (still work via compatibility layer)
   from hoomd.plugin.cavitymd import CavityForce
   
   # New recommended imports
   from hoomd.cavitymd.forces import CavityForce
   from hoomd.cavitymd.variants import StepVariant
   from hoomd.cavitymd.controllers import DiffEqController
   from hoomd.cavitymd.data import ObservableWriter

**New features to explore:**

1. Time-varying coupling with complete variant system
2. Advanced controllers (PID, MPC, DiffEq)
3. HDF5-based data management
4. FDR temperature measurement
5. Comprehensive energy tracking
6. Multi-mode temperature decomposition

Performance Tips
================

**GPU optimization:**

- Use GPU device for systems with >1000 particles
- Enable CUDA-aware MPI for multi-GPU runs
- Monitor GPU memory usage for large systems

**CPU optimization:**

- Enable OpenMP threading (set OMP_NUM_THREADS)
- Use efficient thermostat combinations
- Optimize timestep for your system

**Memory management:**

- Use appropriate output periods to control memory
- Enable HDF5 compression for large datasets
- Use GSD compression for trajectory files
- Clear autocorrelation buffers periodically for long runs

**I/O optimization:**

- Use HDF5 SWMR mode for concurrent access
- Batch output operations
- Use chunked storage for time-series data

See :doc:`../part3_advanced/performance` for detailed optimization strategies.

Cross-References
================

- **Getting Started**: :doc:`../part1_application/getting_started`
- **Running Simulations**: :doc:`../part1_application/running_simulations`
- **Time-Varying Coupling**: :doc:`../part1_application/time_varying_coupling`
- **Advanced Controllers**: :doc:`../part3_advanced/controllers`
- **FDR Temperature**: :doc:`../part3_advanced/fdr_temperature`
- **Data Management**: :doc:`../part1_application/data_management`
- **Theory Background**: :doc:`../part2_theory/introduction`
