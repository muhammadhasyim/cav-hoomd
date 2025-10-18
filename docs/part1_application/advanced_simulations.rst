=====================
Advanced Simulations
=====================

This guide covers advanced simulation techniques including time-varying coupling, custom Python scripts, and GPU acceleration.

.. contents:: In this Guide
   :local:
   :depth: 2

Time-Varying Coupling Experiments
==================================

Cavity HOOMD supports time-dependent coupling strength, enabling study of non-equilibrium cavity-molecule dynamics.

Step Function Coupling
-----------------------

**Basic Example:**

.. code-block:: bash

   python examples/05_advanced_run.py \
       --coupling 1e-3 \
       --switch-time 1.0 \
       --runtime 1000

**What Happens:**

1. **t < 1.0 ps:** System evolves without cavity coupling (g=0)
2. **t = 1.0 ps:** Coupling instantaneously switches to g=1×10⁻³
3. **t > 1.0 ps:** System evolves with full cavity coupling

**Physical Interpretation:**

This simulates a pump-probe experiment where the cavity is suddenly activated by an external laser pulse.

**Energy Redistribution:**

When coupling switches on, energy redistributes between molecular and cavity degrees of freedom:

.. math::

   E_{\text{total}}(t^-) = E_{\text{total}}(t^+)

Energy is conserved during the switch, but partitioning changes:

- Before (t<1.0 ps): All energy in molecular DOF
- After (t>1.0 ps): Energy shared between molecules and cavity

**Typical Switch Times:**

.. code-block:: bash

   # Fast switch (observes immediate response)
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 0.5 --runtime 500
   
   # Delayed switch (long equilibration before switching)
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 100.0 --runtime 500
   
   # Very delayed switch (study long-term equilibrated state before switch)
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 200.0 --runtime 500

Cavity Particle Displacement in Finite-q Mode
----------------------------------------------

**Important Behavior:** When using finite-q mode with time-varying coupling, the cavity particle position jumps to its new equilibrium at the switch time.

.. code-block:: bash

   python examples/05_advanced_run.py \
       --coupling 1e-3 \
       --switch-time 1.0 \
       --finite-q \
       --runtime 1000

**Physical Reason:**

At t<1.0 ps, with g=0, the cavity particle equilibrium position is at the origin. At t=1.0 ps, when g switches to 1×10⁻³, the new equilibrium position is:

.. math::

   \vec{q}_{\text{new}} = -\frac{g}{K} \vec{D}_{\text{total}}

where :math:`\vec{D}_{\text{total}}` is the total molecular dipole moment and :math:`K` is the cavity spring constant.

**The particle instantaneously moves to this new position** to maintain energy conservation and avoid transient oscillations.

**Observing the Jump:**

.. code-block:: python

   import gsd.hoomd
   import numpy as np
   import matplotlib.pyplot as plt
   
   traj = gsd.hoomd.open('cavity_coupling_1e-03_switch_1.0ps/prod-1.gsd')
   
   cavity_positions = []
   times = []
   for frame in traj:
       # Cavity particles are last 2 particles
       cavity_pos = frame.particles.position[-2:]
       cavity_positions.append(np.linalg.norm(cavity_pos[0]))  # X-polarization
       times.append(frame.configuration.step * 0.001)  # Convert to ps
   
   plt.figure(figsize=(10, 5))
   plt.plot(times, cavity_positions)
   plt.axvline(1.0, color='r', linestyle='--', label='Switch time')
   plt.xlabel('Time (ps)')
   plt.ylabel('Cavity Position (x-component)')
   plt.legend()
   plt.title('Cavity Particle Displacement at Switch')
   plt.show()

You'll see a clear jump at t=1.0 ps.

Energy Conservation Monitoring
-------------------------------

**Checking Energy Conservation:**

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   # Load energy data
   data = np.loadtxt('cavity_coupling_1e-03_switch_1.0ps/prod-1_energy_tracker.txt',
                     skiprows=1, delimiter='\t')
   
   time = data[:, 0]
   total_energy = data[:, -1]
   
   # Calculate energy drift
   switch_idx = np.argmin(np.abs(time - 1.0))
   
   # Before switch
   E_before = total_energy[:switch_idx]
   drift_before = (E_before[-1] - E_before[0]) / E_before[0] * 100
   
   # After switch
   E_after = total_energy[switch_idx:]
   drift_after = (E_after[-1] - E_after[0]) / E_after[0] * 100
   
   print(f"Energy drift before switch: {drift_before:.4f}%")
   print(f"Energy drift after switch: {drift_after:.4f}%")
   
   # Plot
   plt.figure(figsize=(12, 5))
   plt.plot(time, total_energy)
   plt.axvline(1.0, color='r', linestyle='--', label='Switch time')
   plt.xlabel('Time (ps)')
   plt.ylabel('Total Energy')
   plt.title('Energy Conservation During Coupling Switch')
   plt.legend()
   plt.tight_layout()
   plt.show()

**Expected Results:**

- Energy conservation better than 0.01% both before and after switch
- Small discontinuity at switch time due to cavity particle repositioning (finite-q mode)
- No long-term drift

**Troubleshooting Energy Drift:**

If energy drifts more than 0.1%:

1. Reduce timestep
2. Check initial configuration quality
3. Verify thermostat settings
4. See :doc:`../part2_theory/energy_conservation`

Custom Python Scripts
=====================

For full control, write custom Python scripts using Cavity HOOMD classes.

Minimal Working Example
------------------------

**File: minimal_cavity_simulation.py**

.. code-block:: python

   #!/usr/bin/env python3
   """Minimal cavity-coupled MD simulation."""
   
   import hoomd
   import numpy as np
   from hoomd.cavitymd.forces import CavityForce
   from hoomd.bussi_reservoir import BussiReservoir
   
   # Initialize
   device = hoomd.device.auto_select()
   sim = hoomd.Simulation(device=device, seed=42)
   
   # Load initial configuration
   sim.create_state_from_gsd('init-0.gsd')
   
   # Create cavity force
   cavity_force = CavityForce(
       kvector=[0, 0, 0],  # q=0 mode
       couplstr=0.001,      # Coupling strength
       omegac=0.01,         # Cavity frequency (reduced units)
       phmass=1.0           # Cavity particle mass
   )
   
   # Setup integrator
   integrator = hoomd.md.Integrator(dt=0.001)
   
   # Add all forces
   lj = hoomd.md.pair.LJ(nlist=hoomd.md.nlist.Cell(buffer=0.4))
   lj.params[('O', 'O')] = {'epsilon': 1.0, 'sigma': 1.0}
   lj.r_cut[('O', 'O')] = 2.5
   
   integrator.forces.append(lj)
   integrator.forces.append(cavity_force)
   
   # Add thermostats
   molecular_filter = hoomd.filter.Type(['O'])
   bussi = BussiReservoir(kT=0.01, tau=1.0)  # 100 K
   integrator.methods.append(
       hoomd.md.methods.ConstantVolume(
           filter=molecular_filter,
           thermostat=bussi
       )
   )
   
   cavity_filter = hoomd.filter.Type(['Cavity'])
   langevin = hoomd.md.methods.Langevin(
       filter=cavity_filter,
       kT=0.01,
       default_gamma=0.1
   )
   integrator.methods.append(langevin)
   
   sim.operations.integrator = integrator
   
   # Add output
   gsd_writer = hoomd.write.GSD(
       filename='trajectory.gsd',
       trigger=hoomd.trigger.Periodic(1000),
       mode='wb'
   )
   sim.operations.writers.append(gsd_writer)
   
   # Run
   print("Running simulation...")
   sim.run(100000)  # 100 ps
   print("Complete!")

**Running:**

.. code-block:: bash

   python minimal_cavity_simulation.py

Using CavityMDSimulation Class
-------------------------------

For more features, use the high-level ``CavityMDSimulation`` class:

.. code-block:: python

   #!/usr/bin/env python3
   """Advanced simulation using CavityMDSimulation class."""
   
   import hoomd
   from cavitymd.simulation import CavityMDSimulation
   from cavitymd.analysis import EnergyTracker
   
   # Create device
   device = hoomd.device.auto_select()
   
   # Initialize simulation with comprehensive setup
   sim = CavityMDSimulation(
       device=device,
       initial_state='init-0.gsd',
       coupling_strength=1e-3,
       cavity_frequency_cm=2000.0,
       temperature_K=100.0,
       molecular_thermostat='bussi',
       cavity_thermostat='langevin',
       timestep_ps=0.001,
       finite_q=False,
       seed=42
   )
   
   # Add energy tracking
   energy_tracker = EnergyTracker(
       simulation=sim.sim,
       output_period_ps=1.0,
       output_file='energy.txt'
   )
   sim.sim.operations.writers.append(energy_tracker)
   
   # Run equilibration
   print("Equilibrating...")
   sim.run_equilibration(n_steps=10000)
   
   # Run production
   print("Production...")
   sim.run_production(n_steps=100000, output_gsd='trajectory.gsd')

**Advantages of CavityMDSimulation:**

- Automatic parameter validation
- Built-in equilibration protocols
- Integrated analysis tools
- Error checking and diagnostics

Time-Varying Coupling with Python API
--------------------------------------

**Using StepVariant:**

.. code-block:: python

   from cavitymd.variants import StepVariant
   from hoomd.cavitymd.forces import CavityForce
   
   # Create time-varying coupling
   coupling_variant = StepVariant(
       initial_value=0.0,
       final_value=1e-3,
       switch_time_ps=1.0,
       dt_ps=0.001
   )
   
   # Create cavity force with variant
   cavity_force = CavityForce(
       kvector=[0, 0, 0],
       couplstr=coupling_variant,  # Pass variant instead of constant
       omegac=0.01,
       phmass=1.0
   )
   
   # ... rest of simulation setup ...

**Other Variants:**

.. code-block:: python

   from cavitymd.variants import (
       ExponentialDecayVariant,
       SquareWaveVariant,
       PeriodicVariant
   )
   
   # Exponential decay
   decay = ExponentialDecayVariant(
       initial_value=1e-3,
       final_value=1e-4,
       decay_time_ps=10.0,
       start_time_ps=1.0,
       dt_ps=0.001
   )
   
   # Square wave (on-off-on-off...)
   square_wave = SquareWaveVariant(
       amplitude=1e-3,
       period_ps=10.0,
       duty_cycle=0.5,  # 50% on, 50% off
       dt_ps=0.001
   )
   
   # Sinusoidal modulation
   periodic = PeriodicVariant(
       amplitude=1e-3,
       frequency_cm=100.0,
       phase_rad=0.0,
       offset=5e-4,
       dt_ps=0.001
   )

Custom Analysis Integration
----------------------------

**Real-Time Temperature Monitoring:**

.. code-block:: python

   from cavitymd.molecular_temperatures import DiatomicMolecularTemperatures
   from cavitymd.analysis import ElapsedTimeTracker
   
   # Create time tracker
   time_tracker = ElapsedTimeTracker(
       simulation=sim,
       runtime_ps=1000.0
   )
   
   # Create temperature tracker
   temp_tracker = DiatomicMolecularTemperatures(
       simulation=sim,
       time_tracker=time_tracker,
       output_period_ps=1.0,
       output_file='molecular_temps.csv'
   )
   
   # Add to simulation
   sim.operations.writers.append(time_tracker)
   sim.operations.writers.append(temp_tracker)
   
   # Run - temperatures will be logged automatically
   sim.run(100000)

**Custom Callback Functions:**

.. code-block:: python

   class CustomAnalyzer(hoomd.custom.Action):
       """Custom analysis performed every N steps."""
       
       def __init__(self, simulation):
           self.sim = simulation
           
       def act(self, timestep):
           """Called every trigger interval."""
           # Get current snapshot
           snapshot = self.sim.state.get_snapshot()
           
           # Perform custom analysis
           positions = snapshot.particles.position
           velocities = snapshot.particles.velocity
           
           # Calculate something interesting
           kinetic_energy = 0.5 * np.sum(velocities**2)
           
           print(f"Step {timestep}: KE = {kinetic_energy:.4f}")
   
   # Create and add custom analyzer
   analyzer = CustomAnalyzer(sim)
   custom_trigger = hoomd.trigger.Periodic(1000)
   custom_writer = hoomd.write.CustomWriter(
       action=analyzer,
       trigger=custom_trigger
   )
   sim.operations.writers.append(custom_writer)

GPU Acceleration
================

Building with GPU Support
--------------------------

**Requirements:**

1. NVIDIA GPU with CUDA Compute Capability ≥3.5
2. CUDA Toolkit 11.0 or later
3. GPU-enabled HOOMD-blue installation

**Build HOOMD-blue with GPU:**

.. code-block:: bash

   cd hoomd-blue
   cmake -B build -S . -DENABLE_GPU=ON -DCMAKE_CUDA_ARCHITECTURES=70
   cd build
   make -j$(nproc)
   make install

Replace ``70`` with your GPU's compute capability:

- RTX 30 series: 86
- RTX 20 series: 75
- GTX 10 series: 61
- Tesla V100: 70

**Build Cavity HOOMD with GPU:**

.. code-block:: bash

   cd cav-hoomd
   ./build_install.sh
   # GPU support is detected automatically if CUDA is available

Device Selection in Scripts
----------------------------

**Automatic Selection:**

.. code-block:: python

   device = hoomd.device.auto_select()

This chooses GPU if available, otherwise CPU.

**Explicit GPU:**

.. code-block:: python

   device = hoomd.device.GPU()

**Explicit CPU:**

.. code-block:: python

   device = hoomd.device.CPU()

**Multi-GPU Systems:**

.. code-block:: python

   # Use specific GPU
   device = hoomd.device.GPU(gpu_id=0)  # First GPU
   device = hoomd.device.GPU(gpu_id=1)  # Second GPU

Or via environment variable:

.. code-block:: bash

   CUDA_VISIBLE_DEVICES=0 python script.py  # Use GPU 0
   CUDA_VISIBLE_DEVICES=1 python script.py  # Use GPU 1

Performance Optimization
------------------------

**Memory Considerations:**

GPU memory limits system size. For a typical O₂ system:

- 512 molecules: ~100 MB
- 2048 molecules: ~400 MB
- 8192 molecules: ~1.6 GB

Modern GPUs (8+ GB) can handle 10,000+ molecules easily.

**Optimal Performance Tips:**

1. **Neighbor List Buffer:**
   
   .. code-block:: python
   
      nlist = hoomd.md.nlist.Cell(buffer=0.4)  # Good default
      nlist = hoomd.md.nlist.Cell(buffer=0.8)  # Less frequent rebuilds

2. **Output Frequency:**
   
   .. code-block:: python
   
      # Too frequent (slow)
      gsd_writer = hoomd.write.GSD(filename='traj.gsd', trigger=hoomd.trigger.Periodic(10))
      
      # Good balance
      gsd_writer = hoomd.write.GSD(filename='traj.gsd', trigger=hoomd.trigger.Periodic(1000))

3. **Timestep:**
   
   Larger timesteps = fewer GPU kernel launches:
   
   .. code-block:: python
   
      # Small timestep (more accurate, slower)
      integrator = hoomd.md.Integrator(dt=0.0005)
      
      # Larger timestep (less accurate, faster)
      integrator = hoomd.md.Integrator(dt=0.002)
   
   Balance accuracy vs speed with convergence tests.

Benchmarking
------------

**Simple Benchmark:**

.. code-block:: python

   import hoomd
   import time
   
   device_cpu = hoomd.device.CPU()
   device_gpu = hoomd.device.GPU()
   
   def benchmark(device, n_steps=10000):
       sim = hoomd.Simulation(device=device, seed=42)
       sim.create_state_from_gsd('init-0.gsd')
       
       # ... setup forces, integrator, etc. ...
       
       start = time.time()
       sim.run(n_steps)
       elapsed = time.time() - start
       
       tps = n_steps / elapsed
       return tps
   
   print(f"CPU: {benchmark(device_cpu):.1f} steps/second")
   print(f"GPU: {benchmark(device_gpu):.1f} steps/second")

**Expected Performance:**

- CPU (8-core): 500-2000 steps/sec
- GPU (RTX 3080): 50,000-150,000 steps/sec
- GPU (A100): 100,000-300,000 steps/sec

Speedup typically 50-100x for cavity simulations.

Next Steps
==========

You now understand:

- Time-varying coupling experiments and energy conservation
- Custom Python scripts with full API control
- Time-dependent coupling variants
- Custom analysis integration
- GPU acceleration and optimization

Continue to:

- :doc:`analysis_tools` for comprehensive analysis methods
- :doc:`time_varying_coupling` for detailed time-varying coupling guide
- :doc:`../part2_theory/time_varying_coupling` for theoretical background
- :doc:`../part3_advanced/performance` for performance optimization

