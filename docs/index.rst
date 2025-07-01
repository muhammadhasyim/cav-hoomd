===============================
Cavity HOOMD Documentation
===============================

**Cavity HOOMD** enables molecular dynamics simulations with optical cavity coupling using HOOMD-blue. 

The easiest way to get started is using the **05_advanced_run.py** script.

Quick Start
===========

**Installation:**

.. code-block:: bash

   git clone https://github.com/yourusername/cavity-hoomd.git
   cd cavity-hoomd
   cmake -B build -S .
   cmake --build build
   cmake --install build

**Run your first simulation:**

.. code-block:: bash

   # Basic cavity simulation
   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000

   # Control simulation (no cavity)
   python examples/05_advanced_run.py --no-cavity --runtime 1000

That's it! Your simulation will create output files with trajectory data and energy tracking.

Example Usage
=============

**Basic Parameters**

.. code-block:: bash

   # Run with specific parameters
   python examples/05_advanced_run.py \
       --coupling 1e-3 \
       --temperature 100 \
       --frequency 2000 \
       --runtime 1000

**Parameter Sweeps**

.. code-block:: bash

   # Sweep over coupling strengths
   python examples/05_advanced_run.py --coupling 1e-3,1e-4,1e-5 --runtime 500

   # Multiple replicas
   python examples/05_advanced_run.py --coupling 1e-3 --replicas 1-5 --runtime 1000

**Available Options**

- ``--coupling`` - Coupling strength (e.g., 1e-3)
- ``--temperature`` - Temperature in K (default: 100)
- ``--frequency`` - Cavity frequency in cm⁻¹ (default: 2000)
- ``--runtime`` - Simulation time in ps (default: 500)
- ``--no-cavity`` - Run without cavity (control simulation)
- ``--replicas`` - Run multiple replicas (e.g., "1-5")
- ``--device GPU`` - Use GPU acceleration
- ``--enable-energy-tracker`` - Detailed energy tracking

Output Files
============

Each simulation creates:

- ``prod-X.gsd`` - Trajectory file
- ``prod-X-energy.txt`` - Energy tracking data
- ``prod-X-cavity_mode.txt`` - Cavity mode properties

Jupyter Notebook
================

See ``examples/05_advanced_run.ipynb`` for an interactive example with the same functionality.

What It Does
============

Cavity HOOMD adds cavity-molecule coupling using the interaction:

.. math::

   H = \frac{1}{2} K q^2 + g \vec{q} \cdot \vec{d} + \frac{g^2}{2K} d^2

Where:
- :math:`q` is the cavity mode position
- :math:`d` is the molecular dipole moment  
- :math:`g` is the coupling strength
- :math:`K` is the cavity spring constant

Python API
===========

For programmatic access:

.. code-block:: python

   import hoomd
   from hoomd.cavitymd import CavityForce
   from hoomd.bussi_reservoir import BussiReservoir

   # Create simulation
   sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=42)
   sim.create_state_from_gsd('system.gsd')

   # Add cavity force
   cavity_force = CavityForce(
       kvector=[0, 0, 1],
       couplstr=0.001,
       omegac=0.01,
       phmass=1.0
   )

   # Setup integrator
   integrator = hoomd.md.Integrator(dt=0.001)
   integrator.forces.append(cavity_force)

   # Add thermostat
   molecular_thermostat = BussiReservoir(kT=0.01, tau=1.0)
   molecular_filter = hoomd.filter.Type(['O', 'N'])
   integrator.methods.append(
       hoomd.md.methods.ConstantVolume(
           filter=molecular_filter, thermostat=molecular_thermostat
       )
   )

   sim.operations.integrator = integrator
   sim.run(10000)

Help and Support
================

- Get all script options: ``python examples/05_advanced_run.py --help``
- **Issues**: https://github.com/yourusername/cavity-hoomd/issues
- **Questions**: https://github.com/yourusername/cavity-hoomd/discussions

Contents
========

.. toctree::
   :maxdepth: 1

   installation
   quickstart
   theory
   license

Citation
========

.. code-block:: bibtex

   @software{cavity_hoomd_2025,
     title={Cavity HOOMD: Molecular Dynamics with Optical Cavity Coupling},
     author={Development Team},
     year={2025},
     url={https://github.com/yourusername/cavity-hoomd}
   } 