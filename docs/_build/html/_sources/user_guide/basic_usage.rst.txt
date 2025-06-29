===========
Basic Usage
===========

Getting Started
===============

This guide covers the fundamental usage patterns for Cavity HOOMD. We'll walk through setting up your first cavity-coupled molecular dynamics simulation.

Prerequisites
=============

Before starting, ensure you have:

* Cavity HOOMD installed (:doc:`../installation`)
* Basic familiarity with HOOMD-blue
* A molecular system prepared as a GSD file

Your First Simulation
=====================

Here's a minimal example to get you started:

.. code-block:: python

   import hoomd
   from hoomd.cavitymd import CavityForce

   # Initialize HOOMD simulation
   device = hoomd.device.CPU()  # or hoomd.device.GPU()
   sim = hoomd.Simulation(device=device)

   # Load your molecular system
   sim.create_state_from_gsd("molecular_system.gsd")

   # Set up molecular forces (example for Lennard-Jones)
   cell = hoomd.md.nlist.Cell(buffer=0.4)
   lj = hoomd.md.pair.LJ(nlist=cell)
   lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)
   lj.r_cut[('A', 'A')] = 2.5

   # Add cavity coupling
   cavity_force = CavityForce(
       kvector=[0, 0, 1],     # Cavity wave vector
       couplstr=0.001,        # Coupling strength
       omegac=0.1             # Cavity frequency (atomic units)
   )

   # Set up integrator
   integrator = hoomd.md.Integrator(dt=0.005)
   integrator.forces = [lj, cavity_force]

   # Add thermostats
   nve = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
   integrator.methods = [nve]
   sim.operations.integrator = integrator

   # Run simulation
   sim.run(10000)

Understanding the Components
============================

CavityForce
-----------

The :class:`~hoomd.cavitymd.CavityForce` is the core component that couples your molecular system to a cavity mode:

.. code-block:: python

   from hoomd.cavitymd import CavityForce

   cavity_force = CavityForce(
       kvector=[0, 0, 1],     # Direction of cavity mode
       couplstr=0.001,        # Light-matter coupling strength  
       omegac=0.1             # Cavity frequency in atomic units
   )

**Parameters:**

* ``kvector``: Direction vector for the cavity electric field
* ``couplstr``: Coupling strength parameter (typical values: 1e-4 to 1e-2)
* ``omegac``: Cavity frequency in atomic units

Physical Constants
------------------

Use the built-in physical constants for unit conversions:

.. code-block:: python

   from hoomd.cavitymd import PhysicalConstants

   # Convert cavity frequency from cm⁻¹ to atomic units
   freq_cm = 2000.0  # cm⁻¹
   freq_au = freq_cm / PhysicalConstants.HARTREE_TO_CM_MINUS1

   # Convert time from ps to atomic units  
   time_ps = 1.0  # ps
   time_au = PhysicalConstants.ps_to_atomic_units(time_ps)

Setting Up Thermostats
======================

Molecular Thermostat
--------------------

For the molecular degrees of freedom:

.. code-block:: python

   # Langevin thermostat
   molecular_thermostat = hoomd.md.methods.Langevin(
       filter=hoomd.filter.Type(['A', 'B']),  # Apply to molecular particles
       kT=1.0,          # Temperature
       default_gamma=1.0 # Friction coefficient
   )

   # Or Bussi thermostat (requires separate installation)
   from hoomd.bussi_reservoir.thermostats import BussiReservoir
   bussi = BussiReservoir(kT=1.0, tau=10.0)
   molecular_thermostat = hoomd.md.methods.ConstantVolume(
       filter=hoomd.filter.Type(['A', 'B']),
       thermostat=bussi
   )

Cavity Thermostat
-----------------

For the cavity mode (if you add a cavity particle):

.. code-block:: python

   # Langevin thermostat for cavity particle
   cavity_thermostat = hoomd.md.methods.Langevin(
       filter=hoomd.filter.Type(['L']),  # Cavity particle type
       kT=1.0,
       default_gamma=0.1  # Lower friction for cavity mode
   )

Output and Analysis
===================

GSD Output
----------

Save trajectory data:

.. code-block:: python

   # Create GSD writer
   gsd_writer = hoomd.write.GSD(
       filename="trajectory.gsd",
       trigger=hoomd.trigger.Periodic(1000),  # Every 1000 steps
       mode='wb'
   )
   sim.operations.writers.append(gsd_writer)

Logging
-------

Monitor simulation progress:

.. code-block:: python

   # Create logger
   logger = hoomd.logging.Logger()
   logger.add(sim, quantities=['timestep', 'tps'])

   # Add cavity force energy (if available)
   logger.add(cavity_force, quantities=['energy'])

   # Create table writer
   table = hoomd.write.Table(
       trigger=hoomd.trigger.Periodic(100),
       logger=logger
   )
   sim.operations.writers.append(table)

Common Patterns
===============

Parameter Scanning
------------------

.. code-block:: python

   import numpy as np

   coupling_values = np.logspace(-4, -2, 10)  # 10 coupling strengths

   for i, coupling in enumerate(coupling_values):
       # Create new simulation for each parameter
       sim = setup_simulation(coupling_strength=coupling)
       
       # Set unique output filename
       gsd_writer = hoomd.write.GSD(
           filename=f"traj_coupling_{coupling:.1e}.gsd",
           trigger=hoomd.trigger.Periodic(1000)
       )
       sim.operations.writers.append(gsd_writer)
       
       # Run simulation
       sim.run(50000)

Energy Conservation
-------------------

Monitor total energy:

.. code-block:: python

   # Add thermodynamic compute
   thermo = hoomd.md.compute.ThermodynamicQuantities(
       filter=hoomd.filter.All()
   )
   sim.operations.computes.append(thermo)

   # Log energies
   logger = hoomd.logging.Logger()
   logger.add(thermo, quantities=['kinetic_energy', 'potential_energy'])
   logger.add(cavity_force, quantities=['energy'])

Best Practices
==============

Timestep Selection
------------------

* Start with small timesteps (dt = 0.001 - 0.005 atomic units)
* Monitor energy conservation
* Use adaptive timestep control for long simulations

Coupling Strength
-----------------

* Weak coupling: 1e-4 to 1e-3
* Strong coupling: 1e-3 to 1e-2  
* Very strong coupling: > 1e-2 (may require smaller timesteps)

System Size
-----------

* Larger systems show more pronounced collective effects
* Minimum ~100-1000 particles for meaningful cavity effects
* Consider periodic boundary conditions

Troubleshooting
===============

Common Issues
-------------

**Energy not conserved**
   - Reduce timestep
   - Check force implementations
   - Verify cavity particle setup

**Simulation crashes**
   - Check particle overlaps in initial configuration
   - Reduce coupling strength initially
   - Verify all force parameters

**No cavity effects observed**
   - Increase coupling strength
   - Check cavity frequency matches molecular resonances
   - Ensure sufficient simulation time

Next Steps
==========

* Learn about :doc:`advanced_features` for sophisticated simulations
* Explore :doc:`parameter_sweeps` for systematic studies  
* Check out :doc:`analysis` for post-processing tools 