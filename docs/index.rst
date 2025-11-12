===============================
Cavity HOOMD User Guide
===============================

**Cavity HOOMD** enables molecular dynamics simulations with optical cavity coupling using HOOMD-blue.

This package extends HOOMD-blue to simulate molecules interacting with optical cavity modes, enabling study of light-matter interactions, polariton chemistry, and cavity-modified molecular dynamics.

.. image:: _static/logo.svg
   :width: 200px
   :align: center
   :alt: Cavity HOOMD Logo

Quick Links
===========

- **New Users:** Start with :doc:`part1_application/getting_started`
- **Running Simulations:** See :doc:`part1_application/running_simulations`
- **Time-Varying Coupling:** Learn about :doc:`part1_application/time_varying_coupling`
- **Data Management:** Modern :doc:`part1_application/data_management` system
- **Controllers:** Advanced :doc:`part3_advanced/controllers`
- **FDR Analysis:** :doc:`part3_advanced/fdr_temperature` measurements
- **Theory Background:** Read :doc:`part2_theory/introduction`
- **API Reference:** Browse :doc:`api/index`
- **Installation Help:** :doc:`installation`

What Can You Do?
================

**Core Capabilities:**

- Simulate cavity-molecule coupling with time-dependent coupling strength
- Study strong coupling effects and polariton formation
- Analyze energy flow between molecular and photonic degrees of freedom
- Use advanced control systems for temperature and parameter management
- GPU-accelerated simulations for large systems
- Real-time data monitoring with HDF5-based output
- Comprehensive energy and temperature tracking
- Fluctuation-dissipation ratio analysis for non-equilibrium systems

**Example Applications:**

- Vibrational strong coupling in molecular ensembles
- Cavity-modified chemical dynamics
- Non-equilibrium switching experiments
- Polariton-mediated energy transfer
- Mode-selective chemistry under strong coupling
- Temperature control via coupling modulation
- Pump-probe dynamics simulations

Quick Example
=============

Get started in minutes:

.. code-block:: bash

   # Install
   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd
   ./build_install.sh

   # Run your first cavity simulation
   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000

**With time-varying coupling:**

.. code-block:: bash

   # Coupling switch at t=10ps
   python examples/05_advanced_run.py \
       --coupling 1e-3 \
       --switch-time 10.0 \
       --runtime 500

**With advanced control:**

.. code-block:: python

   from cavitymd import CavityMDSimulation
   from cavitymd.controllers import PIDControl
   
   # Setup simulation
   sim = CavityMDSimulation(...)
   
   # Add PID temperature controller
   pid = PIDControl(
       simulation=sim,
       setpoint_K=100.0,
       auto_tune=True
   )
   
   sim.run()

That's it! See :doc:`part1_application/getting_started` for detailed instructions.

What's New
==========

**Recent Major Updates:**

**Data Management**

- HDF5-based unified data output system
- 10-100x smaller file sizes
- Real-time monitoring with SWMR mode
- See :doc:`part1_application/data_management`

**Controllers**

- 7 controller types for different applications
- PID with auto-tuning
- Model Predictive Control (MPC)
- Adaptive and empirical controllers
- See :doc:`part3_advanced/controllers`

**Coupling Variants**

- 9 variant types for time-dependent protocols
- Step, periodic, decay, and adaptive variants
- Composite variant system
- See :doc:`part1_application/time_varying_coupling`

**FDR Analysis**

- Complete workflow for effective temperature measurement
- Streaming algorithms for real-time analysis
- Mode-specific temperature decomposition
- See :doc:`part3_advanced/fdr_temperature`

**Performance**

- Optimized GPU kernels
- Improved memory management
- Efficient autocorrelation algorithms
- See :doc:`part3_advanced/performance`

Documentation Structure
=======================

This documentation is organized into three main parts:

**Part I: The Cavity HOOMD Application Layer**

User-focused tutorials and guides for running simulations, from basic examples to advanced time-varying coupling experiments and data management.

**Key Sections:**

- :doc:`part1_application/getting_started` - First steps
- :doc:`part1_application/running_simulations` - Basic simulation workflow
- :doc:`part1_application/advanced_simulations` - Advanced features
- :doc:`part1_application/time_varying_coupling` - Dynamic coupling protocols
- :doc:`part1_application/data_management` - HDF5 data system
- :doc:`part1_application/analysis_tools` - Analysis and tracking

**Part II: Theory and Implementation**

Mathematical foundations, equations of motion, energy conservation, and physical interpretation of cavity-molecule coupling.

**Key Sections:**

- :doc:`part2_theory/introduction` - Theoretical foundations
- :doc:`part2_theory/cavity_forces` - Force computation
- :doc:`part2_theory/time_varying_coupling` - Time-dependent theory
- :doc:`part2_theory/thermostats` - Thermostat implementations
- :doc:`part2_theory/energy_conservation` - Energy conservation
- :doc:`part2_theory/strong_coupling` - Strong coupling physics

**Part III: Advanced Features and Analysis**

Advanced controllers, FDR temperature measurement, molecular temperature decomposition, correlation analysis, and performance optimization.

**Key Sections:**

- :doc:`part3_advanced/controllers` - Control systems
- :doc:`part3_advanced/fdr_temperature` - Effective temperature
- :doc:`part3_advanced/molecular_temperatures` - Temperature decomposition
- :doc:`part3_advanced/correlation_analysis` - Correlation functions
- :doc:`part3_advanced/performance` - Optimization strategies

----

.. toctree::
   :maxdepth: 2
   :caption: Part I: Application Layer
   :hidden:
   
   part1_application/getting_started
   part1_application/running_simulations
   part1_application/advanced_simulations
   part1_application/time_varying_coupling
   part1_application/data_management
   part1_application/analysis_tools

.. toctree::
   :maxdepth: 2
   :caption: Part II: Theory and Implementation
   :hidden:
   
   part2_theory/introduction
   part2_theory/cavity_forces
   part2_theory/time_varying_coupling
   part2_theory/thermostats
   part2_theory/energy_conservation
   part2_theory/strong_coupling

.. toctree::
   :maxdepth: 2
   :caption: Part III: Advanced Features
   :hidden:
   
   part3_advanced/controllers
   part3_advanced/fdr_temperature
   part3_advanced/molecular_temperatures
   part3_advanced/correlation_analysis
   part3_advanced/performance

.. toctree::
   :maxdepth: 1
   :caption: Reference
   :hidden:
   
   installation
   reference
   api/index
   license

Feature Highlights
==================

Time-Varying Coupling
----------------------

**Complete variant system** for any dynamic protocol:

.. code-block:: python

   from cavitymd.variants import (
       StepVariant,              # Sudden switching
       PeriodicVariant,          # Sinusoidal modulation
       SquareWaveVariant,        # Pulsed coupling
       ExponentialDecayVariant,  # Gradual changes
       AdaptiveSquareWaveVariant # Temperature-adaptive
   )
   
   # Example: Step switch
   coupling = StepVariant(
       target_value=0.001,
       switch_time_ps=10.0,
       time_tracker=time_tracker
   )

See :doc:`part1_application/time_varying_coupling` for all 9 variant types.

Advanced Controllers
--------------------

**7 controller types** for different objectives:

.. code-block:: python

   from cavitymd.controllers import (
       PIDControl,              # Classical PID
       DiffEqController,        # First-order dynamics
       SimpleSetpointController,# Direct control
       AdaptiveMPCController,   # Model predictive
       DualFeedbackController,  # Multi-mode
       EmpiricalController,     # Data-driven
       FeedbackController       # Base class
   )
   
   # Example: Auto-tuned PID
   pid = PIDControl(
       simulation=sim,
       setpoint_K=100.0,
       auto_tune=True
   )

See :doc:`part3_advanced/controllers` for detailed comparison.

Modern Data Management
----------------------

**HDF5-based system** with live monitoring:

.. code-block:: python

   from cavitymd.data import ObservableWriter
   
   # Single file for all observables
   writer = ObservableWriter(
       output_file='sim_data.h5',
       enable_swmr=True  # Live read access
   )
   
   # 10-100x smaller than text files
   # Hierarchical organization
   # Real-time monitoring while running

See :doc:`part1_application/data_management` for full capabilities.

FDR Temperature Analysis
------------------------

**Effective temperature** via fluctuation-dissipation:

.. code-block:: python

   from cavitymd import FDRWorkflow
   
   # Complete FDR analysis
   workflow = FDRWorkflow(
       observable='dipole',
       probe_frequency_cm=2000.0,
       equilibrium_duration_ps=200.0,
       response_duration_ps=200.0
   )
   
   results = workflow.run()
   print(f"T_eff = {results['T_eff']:.2f} K")

See :doc:`part3_advanced/fdr_temperature` for theory and examples.

GPU Acceleration
----------------

**Optimized CUDA kernels** for high performance:

.. code-block:: python

   import hoomd
   
   # Automatically uses GPU when available
   device = hoomd.device.GPU()
   sim = hoomd.Simulation(device=device)
   
   # All cavity forces run on GPU
   cavity_force = CavityForce(...)

See :doc:`part3_advanced/performance` for optimization tips.

Help and Support
================

- **GitHub Issues:** https://github.com/muhammadhasyim/cav-hoomd/issues
- **Discussions:** https://github.com/muhammadhasyim/cav-hoomd/discussions
- **Documentation:** You're reading it!
- **Examples:** Check `examples/` directory for working scripts

Common Questions
----------------

**Q: How do I start?**

A: Follow :doc:`part1_application/getting_started` for installation and first simulation.

**Q: How do I implement time-varying coupling?**

A: See :doc:`part1_application/time_varying_coupling` for all variant types.

**Q: How do I control temperature?**

A: See :doc:`part3_advanced/controllers` for 7 controller options.

**Q: How do I analyze my data?**

A: Use the :doc:`part1_application/data_management` HDF5 system and :doc:`part1_application/analysis_tools`.

**Q: What about GPU performance?**

A: See :doc:`part3_advanced/performance` for optimization strategies.

Citation
========

If you use Cavity HOOMD in your research, please cite:

.. code-block:: bibtex

   @software{cavity_hoomd_2025,
     title={Cavity HOOMD: Molecular Dynamics with Optical Cavity Coupling},
     author={Development Team},
     year={2025},
     url={https://github.com/muhammadhasyim/cav-hoomd}
   }

License
=======

Cavity HOOMD is released under the BSD 3-Clause License. See :doc:`license` for details.

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
