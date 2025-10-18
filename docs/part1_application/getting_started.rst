===============
Getting Started
===============

This guide will help you install Cavity HOOMD and run your first cavity-coupled molecular dynamics simulation.

System Requirements
===================

**Operating System:**

- Linux (Ubuntu 20.04+, RHEL 8+, or equivalent recommended)
- macOS 10.15+ (experimental support)

**Hardware:**

- CPU: Multi-core processor (Intel/AMD x86_64 or ARM64)
- RAM: 4 GB minimum, 16 GB+ recommended for large systems
- GPU (optional): NVIDIA GPU with CUDA 11.0+ for acceleration

**Software Prerequisites:**

- Python 3.8 or later
- CMake 3.15 or later
- C++ compiler with C++17 support (GCC 7+, Clang 5+, or MSVC 2017+)
- CUDA Toolkit 11.0+ (optional, for GPU support)

Installation
============

Cavity HOOMD is built as a plugin for HOOMD-blue. You must install HOOMD-blue first, then build the Cavity HOOMD plugin.

Step 1: Install HOOMD-blue 5.2.0
---------------------------------

**Important:** Cavity HOOMD requires **exactly** HOOMD-blue version 5.2.0.

**1.1 Install Prerequisites**

Using micromamba (recommended):

.. code-block:: bash

   micromamba install cmake eigen git numpy pybind11 python scipy

Or using conda:

.. code-block:: bash

   conda install cmake eigen git numpy pybind11 python scipy

**1.2 Clone HOOMD-blue**

.. code-block:: bash

   git clone --recursive https://github.com/glotzerlab/hoomd-blue.git --branch v5.2.0
   cd hoomd-blue

**1.3 Configure and Build**

For GPU support:

.. code-block:: bash

   cmake -B build -S . -DENABLE_GPU=ON
   cd build
   make -j$(nproc)

For CPU only:

.. code-block:: bash

   cmake -B build -S . -DENABLE_GPU=OFF
   cd build
   make -j$(nproc)

**1.4 Install HOOMD-blue**

.. code-block:: bash

   make install

**1.5 Verify Installation**

.. code-block:: bash

   python3 -c "import hoomd; print('HOOMD version:', hoomd.version.version)"
   python3 -c "import hoomd; assert hoomd.version.version == '5.2.0'"

If both commands succeed, you're ready to install Cavity HOOMD.

Step 2: Install Cavity HOOMD Plugin
------------------------------------

**2.1 Clone the Repository**

.. code-block:: bash

   git clone https://github.com/muhammadhasyim/cav-hoomd.git
   cd cav-hoomd

**2.2 Build and Install**

The repository includes a convenient build script:

.. code-block:: bash

   ./build_install.sh

For CPU-only builds:

.. code-block:: bash

   ./build_install.sh --no-gpu

The script will:

- Verify HOOMD-blue 5.2.0 is installed
- Build the cavity force C++/CUDA extensions
- Build the Bussi thermostat extensions
- Install both plugins into your Python environment

**2.3 Verify Installation**

.. code-block:: bash

   python3 -c "import hoomd.cavitymd; print('Cavity HOOMD installed successfully')"
   python3 -c "import hoomd.bussi_reservoir; print('Bussi thermostat installed successfully')"

If both imports succeed, the installation is complete!

Your First Simulation
======================

Now let's run a simple cavity-coupled molecular dynamics simulation.

Example: "Hello Cavity" - Diatomic Molecules in a Cavity
---------------------------------------------------------

This example simulates a system of diatomic molecules (O₂) interacting with a cavity mode.

**Step 1: Navigate to Examples**

.. code-block:: bash

   cd /path/to/cav-hoomd/examples

**Step 2: Check Available Options**

.. code-block:: bash

   python 05_advanced_run.py --help

You'll see all available command-line options.

**Step 3: Run Your First Simulation**

.. code-block:: bash

   python 05_advanced_run.py --coupling 1e-3 --runtime 500 --temperature 100

This command will:

- Set cavity-molecule coupling strength to 1×10⁻³ (atomic units)
- Run for 500 ps of simulation time
- Maintain temperature at 100 K
- Use default cavity frequency of 2000 cm⁻¹
- Create output files in a new directory

**What Happens During the Simulation:**

1. **Initialization** (~10 seconds):
   
   - Loads equilibrated O₂ molecular system
   - Creates cavity particles for photonic mode
   - Sets up cavity force with specified coupling
   - Configures Bussi thermostat for molecules
   - Configures Langevin thermostat for cavity mode

2. **Equilibration** (~30 seconds):
   
   - System equilibrates with cavity coupling active
   - Energy is exchanged between molecules and cavity
   - Temperature stabilizes at target value

3. **Production** (~2-5 minutes):
   
   - Trajectory is recorded to GSD file
   - Energy components are tracked
   - Cavity mode position/velocity logged
   - Console shows progress updates

**Expected Console Output:**

.. code-block:: text

   Cavity HOOMD Simulation Starting
   ================================================================================
   
   System Configuration:
     - N_molecules: 512
     - Temperature: 100.0 K
     - Cavity frequency: 2000.0 cm⁻¹
     - Coupling strength: 1.000e-03
     - Cavity mode: q=0 (uniform)
   
   Thermostats:
     - Molecular: Bussi (τ=1.0 ps)
     - Cavity: Langevin (γ=0.1)
   
   Output directory: cavity_coupling_1e-03/
   
   Running equilibration (10000 steps)...
   [████████████████████] 100% Complete
   
   Running production (500000 steps)...
   [█                   ] 5% | T=100.2 K | E=-2456.3 | Step 25000/500000

Understanding the Output
------------------------

After the simulation completes, you'll find a new directory:

.. code-block:: text

   cavity_coupling_1e-03/
   ├── prod-1.gsd                    # Trajectory file
   ├── prod-1_energy_tracker.txt     # Energy time series
   ├── prod-1_cavity_mode.txt        # Cavity mode properties
   └── simulation.log                # Detailed log

**Trajectory File (prod-1.gsd)**

Contains atomic positions, velocities, and forces at regular intervals. View with:

.. code-block:: python

   import gsd.hoomd
   import numpy as np
   
   traj = gsd.hoomd.open('cavity_coupling_1e-03/prod-1.gsd')
   
   print(f"Number of frames: {len(traj)}")
   print(f"Number of particles: {traj[0].particles.N}")
   
   # Plot cavity particle trajectory
   cavity_positions = []
   for frame in traj:
       # Cavity particles are last in the system
       cavity_pos = frame.particles.position[-2:]  # 2 cavity particles for x,y polarizations
       cavity_positions.append(cavity_pos)
   
   cavity_positions = np.array(cavity_positions)

**Energy File (prod-1_energy_tracker.txt)**

Tab-separated file with columns:

- ``time_ps``: Simulation time in picoseconds
- ``kinetic_energy``: Molecular kinetic energy
- ``potential_energy``: Molecular potential energy
- ``cavity_harmonic``: Cavity mode harmonic energy
- ``coupling_energy``: Dipole-field coupling energy
- ``self_energy``: Dipole self-energy term
- ``total_energy``: Sum of all energy components

Plot the energy:

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   
   data = np.loadtxt('cavity_coupling_1e-03/prod-1_energy_tracker.txt', 
                     skiprows=1, delimiter='\t')
   
   time = data[:, 0]
   total_energy = data[:, -1]
   
   plt.figure(figsize=(10, 5))
   plt.plot(time, total_energy)
   plt.xlabel('Time (ps)')
   plt.ylabel('Total Energy (reduced units)')
   plt.title('Energy Conservation Check')
   plt.tight_layout()
   plt.show()

**Cavity Mode File (prod-1_cavity_mode.txt)**

Contains cavity particle positions and velocities for both polarizations (x and y).

Comparing with Control Simulation
----------------------------------

To understand the effect of the cavity, run a control simulation without coupling:

.. code-block:: bash

   python 05_advanced_run.py --no-cavity --runtime 500 --temperature 100

This creates a ``no_cavity/`` directory. Compare the molecular dynamics with and without the cavity to see cavity-induced effects.

Next Steps
==========

Now that you've run your first simulation, you're ready to:

1. **Explore Parameters:** Try different coupling strengths, temperatures, and cavity frequencies
2. **Learn More Features:** Read :doc:`running_simulations` for thermostat options, replica management, and advanced features
3. **Understand the Physics:** Study :doc:`../part2_theory/introduction` for the theoretical background
4. **Analyze Results:** Learn analysis tools in :doc:`analysis_tools`
5. **Try Time-Varying Coupling:** See :doc:`time_varying_coupling` for dynamic experiments

Troubleshooting
===============

Installation Issues
-------------------

**HOOMD version mismatch:**

.. code-block:: text

   ImportError: HOOMD-blue version mismatch: CavityMD requires exactly 
   version 5.2.0, but found version 5.1.0

**Solution:** You must use exactly HOOMD-blue 5.2.0. Uninstall other versions and follow Step 1 above.

**Cannot find CUDA:**

.. code-block:: text

   CMake Error: Could not find CUDA

**Solution:** Either install CUDA Toolkit 11.0+, or use CPU-only build with ``./build_install.sh --no-gpu``.

**Build fails with C++ errors:**

**Solution:** Ensure you have a C++17-compatible compiler:

.. code-block:: bash

   g++ --version  # Should be 7.0 or later
   clang++ --version  # Should be 5.0 or later

Runtime Issues
--------------

**"No module named hoomd":**

**Solution:** HOOMD-blue is not installed. Complete Step 1.

**"No module named hoomd.cavitymd":**

**Solution:** Cavity HOOMD plugin not installed. Complete Step 2.

**Simulation crashes or gives NaN energies:**

**Solution:** Reduce timestep, check initial configuration, or try different thermostat settings. See :doc:`running_simulations` for guidance.

Getting Help
============

If you encounter issues:

1. Check this documentation thoroughly
2. Search `GitHub Issues <https://github.com/muhammadhasyim/cav-hoomd/issues>`_
3. Post in `GitHub Discussions <https://github.com/muhammadhasyim/cav-hoomd/discussions>`_
4. File a new issue with:
   
   - Your exact command
   - Full error message
   - HOOMD version (``python -c "import hoomd; print(hoomd.version.version)"``)
   - Operating system and Python version

