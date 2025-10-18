==================
Running Simulations
==================

This guide provides progressively advanced examples of running cavity-coupled molecular dynamics simulations with Cavity HOOMD.

.. contents:: In this Guide
   :local:
   :depth: 2

Basic Cavity Simulations
=========================

Using the 05_advanced_run.py Script
------------------------------------

The primary interface for running simulations is the ``examples/05_advanced_run.py`` script, which provides a comprehensive command-line interface for cavity simulations.

**Basic Syntax:**

.. code-block:: bash

   python examples/05_advanced_run.py [OPTIONS]

**Essential Parameters:**

- ``--coupling VALUE``: Cavity-molecule coupling strength (atomic units)
- ``--temperature VALUE``: System temperature in Kelvin
- ``--frequency VALUE``: Cavity frequency in cm⁻¹
- ``--runtime VALUE``: Simulation time in picoseconds

Choosing Coupling Strength
---------------------------

The coupling strength determines the interaction between molecules and the cavity mode.

**Weak Coupling (g < 10⁻⁴):**

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-5 --runtime 1000

- Perturbative regime
- Minimal modification of molecular dynamics
- Useful for baseline comparisons

**Moderate Coupling (10⁻⁴ ≤ g ≤ 10⁻³):**

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 5e-4 --runtime 1000

- Observable cavity-molecule interactions
- Energy exchange between modes
- Typical for initial studies

**Strong Coupling (g > 10⁻³):**

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000

- Polariton formation regime
- Significant modification of dynamics
- Collective strong coupling effects

**Physical Interpretation:**

The effective collective coupling scales as:

.. math::

   g_{\text{eff}} = g \sqrt{N}

where N is the number of molecules. For N=512 molecules with g=1×10⁻³:

.. math::

   g_{\text{eff}} = 10^{-3} \times \sqrt{512} \approx 0.0226

This collective enhancement enables strong coupling even with modest single-molecule coupling.

Temperature and Frequency Selection
------------------------------------

**Temperature Guidelines:**

.. code-block:: bash

   # Low temperature (minimize thermal fluctuations)
   python examples/05_advanced_run.py --coupling 1e-3 --temperature 50 --runtime 1000
   
   # Room temperature
   python examples/05_advanced_run.py --coupling 1e-3 --temperature 300 --runtime 1000
   
   # High temperature (study thermal effects)
   python examples/05_advanced_run.py --coupling 1e-3 --temperature 500 --runtime 1000

**Cavity Frequency Selection:**

The cavity frequency should be chosen relative to molecular vibrational frequencies:

.. code-block:: bash

   # On-resonance with O-O stretch (~1550 cm⁻¹)
   python examples/05_advanced_run.py --coupling 1e-3 --frequency 1550 --runtime 1000
   
   # Red-detuned
   python examples/05_advanced_run.py --coupling 1e-3 --frequency 1200 --runtime 1000
   
   # Blue-detuned
   python examples/05_advanced_run.py --coupling 1e-3 --frequency 2000 --runtime 1000

**Resonance Condition:**

Maximum coupling effects occur when:

.. math::

   \omega_{\text{cavity}} \approx \omega_{\text{molecular}}

For O₂ molecules, the main vibrational mode is around 1550 cm⁻¹.

Thermostat Options
==================

Cavity HOOMD supports multiple thermostats for both molecular and cavity degrees of freedom.

Molecular Thermostats
---------------------

**Bussi Thermostat (Default, Recommended):**

.. code-block:: bash

   python examples/05_advanced_run.py --molecular-bath bussi --coupling 1e-3 --runtime 1000

**Advantages:**

- Correct canonical ensemble sampling
- Fast equilibration
- Minimal perturbation to dynamics
- Suitable for both equilibrium and dynamical studies

**Theory:** Bussi-Donadio-Parrinello stochastic velocity rescaling thermostat samples the canonical (NVT) ensemble rigorously.

**Langevin Thermostat:**

.. code-block:: bash

   python examples/05_advanced_run.py --molecular-bath langevin --coupling 1e-3 --runtime 1000

**Advantages:**

- Simple and robust
- Well-established in MD community
- Good for strongly coupled systems

**Considerations:** Introduces friction that may affect fast dynamics.

**No Thermostat (NVE):**

.. code-block:: bash

   python examples/05_advanced_run.py --molecular-bath none --coupling 1e-3 --runtime 1000

**Use Cases:**

- Microcanonical ensemble studies
- Energy conservation testing
- Short-time dynamics analysis

**Warning:** Temperature will drift if system is not at equilibrium initially.

Cavity Thermostats
------------------

**Langevin (Default, Recommended):**

.. code-block:: bash

   python examples/05_advanced_run.py --cavity-bath langevin --coupling 1e-3 --runtime 1000

**Advantages:**

- Represents cavity dissipation physically
- Prevents energy accumulation in cavity mode
- Mimics realistic cavity loss

**Friction Coefficient:** Default γ=0.1 corresponds to cavity lifetime τ ≈ 10 ps.

**Bussi Thermostat:**

.. code-block:: bash

   python examples/05_advanced_run.py --cavity-bath bussi --coupling 1e-3 --runtime 1000

**Use Cases:**

- Canonical ensemble for cavity mode
- Comparison with molecular Bussi thermostat
- Studying cavity thermal fluctuations

**No Thermostat (Conservative):**

.. code-block:: bash

   python examples/05_advanced_run.py --cavity-bath none --coupling 1e-3 --runtime 1000

**Use Cases:**

- Undamped cavity dynamics
- Energy conservation testing
- Short-time coherent dynamics

**Warning:** Cavity mode may accumulate large energy without dissipation.

Recommended Thermostat Combinations
------------------------------------

**Standard Setup (Most Common):**

.. code-block:: bash

   python examples/05_advanced_run.py \
       --molecular-bath bussi \
       --cavity-bath langevin \
       --coupling 1e-3 \
       --runtime 1000

- Molecules: Canonical ensemble (Bussi)
- Cavity: Dissipative Langevin (realistic loss)

**Pure Canonical Ensemble:**

.. code-block:: bash

   python examples/05_advanced_run.py \
       --molecular-bath bussi \
       --cavity-bath bussi \
       --coupling 1e-3 \
       --runtime 1000

- Both molecules and cavity in canonical ensemble
- Useful for equilibrium thermodynamics

**Conservative Dynamics:**

.. code-block:: bash

   python examples/05_advanced_run.py \
       --molecular-bath none \
       --cavity-bath none \
       --coupling 1e-3 \
       --runtime 100

- Microcanonical (NVE) ensemble
- Total energy conserved
- Best for short-time dynamics and energy conservation tests

Cavity Modes
============

Cavity HOOMD supports two types of cavity modes.

q=0 Mode (Uniform, Default)
----------------------------

**Description:** Uniform electric field throughout the simulation box.

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --runtime 1000
   # No additional flags needed (default mode)

**Properties:**

- Cavity mode has no spatial variation
- k-vector: :math:`\vec{k} = (0, 0, 0)`
- All molecules couple with same phase
- Maximum collective enhancement
- Two cavity particles (x and y polarizations)

**When to Use:**

- Strong coupling studies
- Maximum collective effects
- Initial explorations

**Cavity Hamiltonian:**

.. math::

   H_{\text{cavity}} = \frac{1}{2}K(q_x^2 + q_y^2) + g(q_x D_x + q_y D_y) + \frac{g^2}{2K}(D_x^2 + D_y^2)

where :math:`D_x, D_y` are total system dipole moments.

Finite-q Mode
-------------

**Description:** Standing wave with spatial modulation.

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --finite-q --runtime 1000

**Properties:**

- k-vector: :math:`\vec{k} = (0, 0, k_z)`
- Spatial phase variation: molecules at different positions couple with different phases
- Reduces collective enhancement
- Cavity particles positioned at finite coordinates

**When to Use:**

- Realistic cavity geometries
- Wave vector dependence studies
- Reduced collective coupling

**Phase Factor:**

Each molecule at position :math:`\vec{r}_n` couples with phase:

.. math::

   \phi_n = \vec{k} \cdot \vec{r}_n = k_z z_n

**Coupling Strength:**

.. math::

   H_{\text{coupling}} = g q \sum_{n=1}^{N} d_n \cos(\phi_n)

The cosine factor reduces collective enhancement compared to q=0 mode.

**Cavity Particle Positioning:**

In finite-q mode, cavity particles are positioned in the simulation box such that their displacement represents the mode amplitude and maintains energy conservation.

Multiple Replicas
=================

Running Independent Replicas
-----------------------------

For statistical analysis, run multiple independent simulations with different random seeds:

**Sequential Replica Numbering:**

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --replicas 1-5 --runtime 1000

This runs 5 replicas (replica IDs: 1, 2, 3, 4, 5) sequentially.

**Custom Replica IDs:**

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --replicas 1,5,10 --runtime 1000

This runs 3 replicas with IDs 1, 5, and 10.

**Large Replica Sets:**

.. code-block:: bash

   python examples/05_advanced_run.py --coupling 1e-3 --replicas 1-50 --runtime 1000

For comprehensive statistics (50 replicas).

Output Organization
-------------------

Each replica creates a separate output file:

.. code-block:: text

   cavity_coupling_1e-03/
   ├── prod-1.gsd
   ├── prod-1_energy_tracker.txt
   ├── prod-1_cavity_mode.txt
   ├── prod-2.gsd
   ├── prod-2_energy_tracker.txt
   ├── prod-2_cavity_mode.txt
   ├── ...
   ├── prod-5.gsd
   ├── prod-5_energy_tracker.txt
   └── prod-5_cavity_mode.txt

**Replica Differences:**

- Random number seed (determined by replica ID)
- Initial velocities (Boltzmann distribution)
- Stochastic thermostat evolution

**Replica Similarities:**

- Same initial positions
- Same parameters (coupling, temperature, frequency)
- Same forces and integration scheme

Statistical Analysis
--------------------

Analyze replica statistics using Python:

.. code-block:: python

   import numpy as np
   import glob
   
   # Load energy from all replicas
   energy_files = glob.glob('cavity_coupling_1e-03/prod-*_energy_tracker.txt')
   
   energies = []
   for file in sorted(energy_files):
       data = np.loadtxt(file, skiprows=1, delimiter='\t')
       energies.append(data[:, -1])  # Total energy column
   
   energies = np.array(energies)
   
   # Calculate statistics
   mean_energy = np.mean(energies, axis=0)
   std_energy = np.std(energies, axis=0)
   sem_energy = std_energy / np.sqrt(len(energies))
   
   print(f"Final energy: {mean_energy[-1]:.2f} ± {sem_energy[-1]:.2f}")

**Standard Error of the Mean:**

For N independent replicas:

.. math::

   \text{SEM} = \frac{\sigma}{\sqrt{N}}

where σ is the standard deviation across replicas.

Device Selection
================

CPU Execution
-------------

.. code-block:: bash

   python examples/05_advanced_run.py --device CPU --coupling 1e-3 --runtime 500

**When to Use:**

- No GPU available
- Small systems (<1000 particles)
- Debugging and testing

**Performance:** ~100-1000 steps/second for typical systems.

GPU Execution
-------------

.. code-block:: bash

   python examples/05_advanced_run.py --device GPU --coupling 1e-3 --runtime 1000

**Requirements:**

- NVIDIA GPU with CUDA 11.0+
- GPU-enabled HOOMD-blue build
- GPU-enabled Cavity HOOMD build

**When to Use:**

- Large systems (>1000 particles)
- Long simulations
- Production runs

**Performance:** ~10,000-100,000 steps/second depending on GPU.

**GPU Selection:**

For multi-GPU systems:

.. code-block:: bash

   CUDA_VISIBLE_DEVICES=0 python examples/05_advanced_run.py --device GPU --coupling 1e-3
   CUDA_VISIBLE_DEVICES=1 python examples/05_advanced_run.py --device GPU --coupling 1e-3

Output and Logging
==================

Standard Output Files
---------------------

**Trajectory File (prod-{ID}.gsd):**

Binary trajectory in GSD format containing:

- Particle positions
- Velocities
- Forces
- Box dimensions
- Particle types

**Energy Tracker (prod-{ID}_energy_tracker.txt):**

Tab-separated values with columns:

- time_ps: Simulation time
- kinetic_energy: Molecular kinetic energy
- potential_energy: Molecular potential energy
- cavity_harmonic: ½K(q_x² + q_y²)
- coupling_energy: g(q_x·D_x + q_y·D_y)
- self_energy: (g²/2K)(D_x² + D_y²)
- total_energy: Sum of all components

**Cavity Mode (prod-{ID}_cavity_mode.txt):**

Cavity particle positions and velocities:

- time_ps: Simulation time
- q_x: X-polarization cavity position
- q_y: Y-polarization cavity position
- v_x: X-polarization cavity velocity
- v_y: Y-polarization cavity velocity

Controlling Output Frequency
-----------------------------

Output frequency is controlled internally. For custom output, use Python API (see :doc:`advanced_simulations`).

Simulation Logs
---------------

Console output shows:

- Initialization parameters
- Progress bars for equilibration and production
- Current temperature, energy, and step number
- Performance metrics (steps/second)

Performance Considerations
==========================

Timestep Selection
------------------

Default timestep: 0.001 ps (1 fs)

**Smaller timestep (more stable):**

Requires modifying the Python script directly (see :doc:`advanced_simulations`).

**Adaptive timestep:**

See :doc:`advanced_simulations` for using ``AdaptiveTimestepUpdater``.

Simulation Length
-----------------

**Short runs (100-500 ps):**

- Initial testing
- Parameter exploration
- Energy conservation checks

**Medium runs (1-10 ns):**

- Equilibrium property calculations
- Structural analysis
- Standard production

**Long runs (>10 ns):**

- Rare event sampling
- Long-time correlation functions
- Slow relaxation processes

System Size
-----------

The default system contains 512 O₂ molecules (1024 atoms + 2 cavity particles).

For different system sizes, you must prepare custom initial configurations (see :doc:`advanced_simulations`).

Next Steps
==========

You now know how to:

- Run basic cavity simulations with various parameters
- Choose appropriate thermostats
- Use q=0 and finite-q cavity modes
- Run multiple replicas for statistics
- Select CPU or GPU execution

Continue to:

- :doc:`advanced_simulations` for custom Python scripts and time-varying coupling
- :doc:`analysis_tools` for analyzing simulation results
- :doc:`../part2_theory/introduction` for theoretical background

