=======================
Time-Varying Coupling
=======================

This guide provides detailed practical instructions for running and analyzing time-varying coupling experiments.

.. contents:: In this Guide
   :local:
   :depth: 2

.. note::

   For theoretical background on time-varying coupling, see :doc:`../part2_theory/time_varying_coupling`.

Overview
========

Time-varying coupling allows you to study non-equilibrium cavity-molecule dynamics by changing the coupling strength during simulation.

**Why Use Time-Varying Coupling?**

- Model pump-probe experiments
- Study non-equilibrium relaxation
- Investigate sudden coupling activation
- Analyze energy redistribution dynamics
- Simulate realistic experimental protocols

**Key Concept:**

Instead of constant coupling :math:`g`, the coupling becomes time-dependent :math:`g(t)`.

Step Function Coupling
=======================

The most common protocol: coupling switches instantly from 0 to a target value.

Basic Example
-------------

**Command:**

.. code-block:: bash

   python examples/05_advanced_run.py \
       --coupling 1e-3 \
       --switch-time 10.0 \
       --runtime 500 \
       --temperature 100

**Protocol:**

.. math::

   g(t) = \begin{cases}
   0 & \text{if } t < 10 \text{ ps} \\
   10^{-3} & \text{if } t \geq 10 \text{ ps}
   \end{cases}

**Timeline:**

1. **t = 0-10 ps:** Equilibration without cavity (free molecular dynamics)
2. **t = 10 ps:** Coupling switches ON instantaneously
3. **t = 10-500 ps:** Evolution with cavity coupling active

Expected Behavior
-----------------

**Before Switch (t < 10 ps):**

- Pure molecular dynamics
- No cavity energy
- No coupling energy
- System equilibrates at target temperature

**At Switch (t = 10 ps):**

- Instantaneous parameter change
- Cavity particle may jump (finite-q mode)
- Energy redistributes between molecular and cavity DOF
- Total energy conserved

**After Switch (t > 10 ps):**

- Coupled molecular-cavity dynamics
- Energy oscillates between molecules and cavity
- System approaches new equilibrium
- Observables may show transient behavior

Analyzing Switch Dynamics
--------------------------

**Load and visualize:**

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   # Load energy data
   data = np.loadtxt('cavity_coupling_1e-03_switch_10.0ps/prod-1_energy_tracker.txt',
                     skiprows=1, delimiter='\t')
   
   time = data[:, 0]
   E_kinetic = data[:, 1]
   E_potential = data[:, 2]
   E_cavity = data[:, 3]
   E_coupling = data[:, 4]
   E_self = data[:, 5]
   E_total = data[:, 6]
   
   # Find switch point
   switch_idx = np.argmin(np.abs(time - 10.0))
   
   # Create figure
   fig, axes = plt.subplots(3, 1, figsize=(12, 10))
   
   # Plot 1: Total energy (conservation check)
   axes[0].plot(time, E_total, 'k-', linewidth=2)
   axes[0].axvline(10.0, color='r', linestyle='--', alpha=0.7, 
                   label='Switch time')
   axes[0].set_ylabel('Total Energy', fontsize=12)
   axes[0].set_title('Energy Conservation Check', fontsize=14)
   axes[0].legend()
   axes[0].grid(True, alpha=0.3)
   
   # Plot 2: Molecular vs Cavity energy
   axes[1].plot(time, E_kinetic + E_potential, label='Molecular', linewidth=2)
   axes[1].plot(time, E_cavity, label='Cavity', linewidth=2)
   axes[1].axvline(10.0, color='r', linestyle='--', alpha=0.7)
   axes[1].set_ylabel('Energy', fontsize=12)
   axes[1].set_title('Energy Partitioning', fontsize=14)
   axes[1].legend()
   axes[1].grid(True, alpha=0.3)
   
   # Plot 3: Cavity-related terms
   axes[2].plot(time, E_cavity, label='Cavity harmonic', linewidth=2)
   axes[2].plot(time, E_coupling, label='Coupling', linewidth=2)
   axes[2].plot(time, E_self, label='Self-energy', linewidth=2)
   axes[2].axvline(10.0, color='r', linestyle='--', alpha=0.7)
   axes[2].set_xlabel('Time (ps)', fontsize=12)
   axes[2].set_ylabel('Energy', fontsize=12)
   axes[2].set_title('Cavity Energy Components', fontsize=14)
   axes[2].legend()
   axes[2].grid(True, alpha=0.3)
   
   plt.tight_layout()
   plt.savefig('switch_dynamics.png', dpi=300)
   plt.show()

**Quantifying the response:**

.. code-block:: python

   # Calculate energy transfer
   E_mol_before = (E_kinetic + E_potential)[switch_idx-100:switch_idx].mean()
   E_mol_after = (E_kinetic + E_potential)[switch_idx+100:switch_idx+200].mean()
   E_cav_after = E_cavity[switch_idx+100:switch_idx+200].mean()
   
   energy_transferred = E_mol_before - E_mol_after
   cavity_fraction = E_cav_after / (E_mol_after + E_cav_after) * 100
   
   print(f"Energy transferred to cavity: {energy_transferred:.2f}")
   print(f"Cavity energy fraction: {cavity_fraction:.1f}%")
   
   # Relaxation time
   # Fit exponential decay to approach equilibrium
   from scipy.optimize import curve_fit
   
   def exponential(t, A, tau, offset):
       return A * np.exp(-t/tau) + offset
   
   # Use coupling energy as observable
   t_fit = time[switch_idx:switch_idx+1000] - time[switch_idx]
   E_fit = E_coupling[switch_idx:switch_idx+1000]
   
   try:
       popt, _ = curve_fit(exponential, t_fit, E_fit, 
                          p0=[E_fit[0], 10.0, E_fit[-100:].mean()])
       tau_relax = popt[1]
       print(f"Relaxation time: {tau_relax:.2f} ps")
   except:
       print("Could not fit relaxation time")

Cavity Particle Jump (Finite-q Mode)
-------------------------------------

**In finite-q mode, cavity particles jump at switch time.**

.. code-block:: bash

   python examples/05_advanced_run.py \
       --coupling 1e-3 \
       --switch-time 10.0 \
       --finite-q \
       --runtime 500

**Observing the jump:**

.. code-block:: python

   import gsd.hoomd
   
   # Load trajectory
   traj = gsd.hoomd.open('cavity_coupling_1e-03_switch_10.0ps/prod-1.gsd')
   
   # Extract cavity positions
   cavity_x_positions = []
   times = []
   
   for frame in traj:
       # Cavity particles are last 2 particles
       cavity_particle = frame.particles.position[-2]  # X-polarization
       cavity_x_positions.append(cavity_particle[0])
       times.append(frame.configuration.step * 0.001)  # Convert to ps
   
   cavity_x_positions = np.array(cavity_x_positions)
   times = np.array(times)
   
   # Plot
   plt.figure(figsize=(12, 6))
   plt.plot(times, cavity_x_positions, linewidth=2)
   plt.axvline(10.0, color='r', linestyle='--', linewidth=2, 
               alpha=0.7, label='Switch time')
   plt.xlabel('Time (ps)', fontsize=14)
   plt.ylabel('Cavity X Position', fontsize=14)
   plt.title('Cavity Particle Jump at Switch', fontsize=16)
   plt.legend(fontsize=12)
   plt.grid(True, alpha=0.3)
   plt.tight_layout()
   plt.savefig('cavity_jump.png', dpi=300)
   plt.show()

**Why does it jump?**

The equilibrium position changes from :math:`\vec{q}=0` (g=0) to:

.. math::

   \vec{q}_{\text{eq}} = -\frac{g}{K} \vec{D}_{\text{total}}

where :math:`\vec{D}_{\text{total}}` is the total molecular dipole moment at t=10 ps.

Multiple Switch Times
---------------------

**Comparing different switch times:**

.. code-block:: bash

   # Early switch (fast response)
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 1.0 --runtime 500
   
   # Medium switch (standard)
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 10.0 --runtime 500
   
   # Late switch (long equilibration)
   python examples/05_advanced_run.py --coupling 1e-3 --switch-time 100.0 --runtime 500

**Analysis:**

.. code-block:: python

   switch_times = [1.0, 10.0, 100.0]
   colors = ['blue', 'green', 'red']
   
   plt.figure(figsize=(12, 6))
   
   for switch_time, color in zip(switch_times, colors):
       dir_name = f'cavity_coupling_1e-03_switch_{switch_time}ps'
       data = np.loadtxt(f'{dir_name}/prod-1_energy_tracker.txt',
                        skiprows=1, delimiter='\t')
       
       time = data[:, 0]
       E_coupling = data[:, 4]
       
       plt.plot(time, E_coupling, color=color, linewidth=2,
               label=f'Switch at {switch_time} ps')
       plt.axvline(switch_time, color=color, linestyle='--', alpha=0.5)
   
   plt.xlabel('Time (ps)', fontsize=14)
   plt.ylabel('Coupling Energy', fontsize=14)
   plt.title('Effect of Switch Time on Coupling Energy', fontsize=16)
   plt.legend(fontsize=12)
   plt.grid(True, alpha=0.3)
   plt.tight_layout()
   plt.savefig('switch_time_comparison.png', dpi=300)
   plt.show()

Common Protocols
================

Exponential Decay
-----------------

**Simulates cavity mode decay or detuning.**

**Python API (required):**

.. code-block:: python

   from cavitymd.variants import ExponentialDecayVariant
   from hoomd.cavitymd.forces import CavityForce
   
   # Create decay variant
   decay_variant = ExponentialDecayVariant(
       initial_value=1e-3,      # Starting coupling
       final_value=1e-4,        # Final coupling (10x weaker)
       decay_time_ps=40.0,      # Decay timescale
       start_time_ps=10.0,      # When decay starts
       dt_ps=0.001              # Simulation timestep
   )
   
   # Create cavity force with variant
   cavity_force = CavityForce(
       kvector=[0, 0, 0],
       couplstr=decay_variant,
       omegac=0.01,
       phmass=1.0
   )
   
   # ... rest of simulation setup ...

**Decay function:**

.. math::

   g(t) = g_{\text{final}} + (g_{\text{initial}} - g_{\text{final}}) e^{-(t-t_0)/\tau}

for :math:`t \geq t_0`, where :math:`\tau` is the decay time.

Square Wave (Periodic On-Off)
------------------------------

**Alternating coupling for pulsed experiments.**

**Python API:**

.. code-block:: python

   from cavitymd.variants import SquareWaveVariant
   
   # Create square wave
   square_wave = SquareWaveVariant(
       amplitude=1e-3,        # Coupling strength when ON
       period_ps=20.0,        # Total period (ON + OFF)
       duty_cycle=0.5,        # 50% ON, 50% OFF
       dt_ps=0.001
   )
   
   cavity_force = CavityForce(
       kvector=[0, 0, 0],
       couplstr=square_wave,
       omegac=0.01,
       phmass=1.0
   )

**Pattern:**

- t=0-10 ps: g = 1e-3 (ON)
- t=10-20 ps: g = 0 (OFF)
- t=20-30 ps: g = 1e-3 (ON)
- ...

Sinusoidal Modulation
----------------------

**Smoothly varying coupling.**

**Python API:**

.. code-block:: python

   from cavitymd.variants import PeriodicVariant
   
   # Create sinusoidal modulation
   periodic = PeriodicVariant(
       amplitude=5e-4,           # Amplitude of oscillation
       frequency_cm=100.0,       # Modulation frequency
       phase_rad=0.0,            # Initial phase
       offset=5e-4,              # DC offset (average coupling)
       dt_ps=0.001
   )
   
   cavity_force = CavityForce(
       kvector=[0, 0, 0],
       couplstr=periodic,
       omegac=0.01,
       phmass=1.0
   )

**Function:**

.. math::

   g(t) = g_0 + A \sin(2\pi f t + \phi)

where:
- :math:`g_0` = offset
- :math:`A` = amplitude
- :math:`f` = frequency
- :math:`\phi` = phase

Analyzing Non-Equilibrium Data
===============================

Time-Resolved Observable
------------------------

**Computing running averages:**

.. code-block:: python

   def running_average(data, window_size):
       """Calculate running average with given window size."""
       cumsum = np.cumsum(np.insert(data, 0, 0))
       return (cumsum[window_size:] - cumsum[:-window_size]) / window_size
   
   # Load temperature data
   data = np.loadtxt('molecular_temps.csv', skiprows=1, delimiter=',')
   time = data[:, 0]
   T_vib = data[:, 3]
   
   # Apply running average (10 ps window)
   window_size = int(10.0 / (time[1] - time[0]))
   T_vib_smooth = running_average(T_vib, window_size)
   time_smooth = time[window_size//2:-window_size//2+1]
   
   # Plot
   plt.figure(figsize=(12, 6))
   plt.plot(time, T_vib, alpha=0.3, label='Raw')
   plt.plot(time_smooth, T_vib_smooth, linewidth=2, label='10 ps average')
   plt.axvline(10.0, color='r', linestyle='--', label='Switch')
   plt.xlabel('Time (ps)')
   plt.ylabel('Vibrational Temperature (K)')
   plt.legend()
   plt.grid(True, alpha=0.3)
   plt.show()

Transient Response Time
-----------------------

**Measuring equilibration time after switch:**

.. code-block:: python

   from scipy.optimize import curve_fit
   
   # Exponential approach to equilibrium
   def approach_equilibrium(t, A, tau, T_eq):
       return T_eq + A * np.exp(-t/tau)
   
   # Extract post-switch data
   switch_time = 10.0
   mask = (time > switch_time) & (time < switch_time + 100.0)
   t_fit = time[mask] - switch_time
   T_fit = T_vib[mask]
   
   # Fit
   try:
       popt, pcov = curve_fit(approach_equilibrium, t_fit, T_fit,
                             p0=[T_fit[0]-T_fit[-1], 20.0, T_fit[-1]])
       
       A, tau, T_eq = popt
       A_err, tau_err, T_eq_err = np.sqrt(np.diag(pcov))
       
       print(f"Equilibration time: {tau:.2f} ± {tau_err:.2f} ps")
       print(f"Equilibrium temperature: {T_eq:.2f} ± {T_eq_err:.2f} K")
       
       # Plot fit
       t_theory = np.linspace(0, 100, 1000)
       T_theory = approach_equilibrium(t_theory, A, tau, T_eq)
       
       plt.figure(figsize=(10, 6))
       plt.plot(t_fit, T_fit, 'o', alpha=0.5, label='Data')
       plt.plot(t_theory, T_theory, 'r-', linewidth=2, 
               label=f'Fit (τ={tau:.1f} ps)')
       plt.xlabel('Time after switch (ps)')
       plt.ylabel('Temperature (K)')
       plt.legend()
       plt.grid(True, alpha=0.3)
       plt.show()
       
   except Exception as e:
       print(f"Fit failed: {e}")

Energy Flow Analysis
--------------------

**Tracking energy flow between subsystems:**

.. code-block:: python

   # Load energy data
   data = np.loadtxt('prod-1_energy_tracker.txt', skiprows=1, delimiter='\t')
   
   time = data[:, 0]
   E_molecular = data[:, 1] + data[:, 2]  # Kinetic + potential
   E_cavity = data[:, 3]  # Cavity harmonic
   
   # Calculate energy flow rate
   dt = time[1] - time[0]
   dE_mol_dt = np.gradient(E_molecular, dt)
   dE_cav_dt = np.gradient(E_cavity, dt)
   
   # Find switch point
   switch_idx = np.argmin(np.abs(time - 10.0))
   
   # Plot around switch
   plot_range = slice(switch_idx-500, switch_idx+2000)
   
   fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
   
   # Energy
   ax1.plot(time[plot_range], E_molecular[plot_range], 
           label='Molecular', linewidth=2)
   ax1.plot(time[plot_range], E_cavity[plot_range], 
           label='Cavity', linewidth=2)
   ax1.axvline(10.0, color='r', linestyle='--', alpha=0.7)
   ax1.set_ylabel('Energy', fontsize=12)
   ax1.set_title('Energy Distribution', fontsize=14)
   ax1.legend()
   ax1.grid(True, alpha=0.3)
   
   # Energy flow rate
   ax2.plot(time[plot_range], dE_mol_dt[plot_range], 
           label='Molecular dE/dt', linewidth=2)
   ax2.plot(time[plot_range], dE_cav_dt[plot_range], 
           label='Cavity dE/dt', linewidth=2)
   ax2.axvline(10.0, color='r', linestyle='--', alpha=0.7)
   ax2.axhline(0, color='k', linestyle=':', alpha=0.5)
   ax2.set_xlabel('Time (ps)', fontsize=12)
   ax2.set_ylabel('Energy Flow Rate', fontsize=12)
   ax2.set_title('Energy Exchange Rate', fontsize=14)
   ax2.legend()
   ax2.grid(True, alpha=0.3)
   
   plt.tight_layout()
   plt.show()

Troubleshooting
===============

Energy Conservation Issues
---------------------------

**Problem:** Energy drifts after switch.

**Solutions:**

1. **Reduce timestep:**
   
   Smaller timestep improves energy conservation.

2. **Check initial configuration:**
   
   Ensure molecules are well-equilibrated before switch.

3. **Verify switch implementation:**
   
   In finite-q mode, cavity particle must jump correctly.

Unexpected Dynamics
-------------------

**Problem:** System behavior doesn't make physical sense.

**Checklist:**

1. **Verify switch time in output:**
   
   Check that coupling actually switches at intended time.

2. **Check coupling strength:**
   
   Very strong coupling (>10⁻²) may cause numerical instability.

3. **Equilibration before switch:**
   
   System should be at equilibrium before coupling switches.

4. **Thermostat settings:**
   
   Thermostats affect relaxation dynamics.

Numerical Instabilities
-----------------------

**Problem:** NaN energies or particle overlaps after switch.

**Solutions:**

1. **Smaller timestep:**
   
   Reduce from 0.001 ps to 0.0005 ps or smaller.

2. **Gradual switch:**
   
   Use exponential approach instead of step function.

3. **Check cavity frequency:**
   
   Very high frequencies require smaller timesteps.

Best Practices
==============

1. **Always run control simulations:**
   
   Compare with constant coupling at same parameters.

2. **Verify energy conservation:**
   
   Total energy should not drift more than 0.01%.

3. **Allow sufficient equilibration:**
   
   Let system equilibrate >50 ps before switch.

4. **Multiple replicas:**
   
   Run several independent simulations for statistics.

5. **Document parameters:**
   
   Record all parameters: coupling, switch time, thermostats, etc.

6. **Check phase space:**
   
   Visualize trajectories to spot anomalies.

7. **Save frequently:**
   
   Output trajectory at high frequency around switch.

Next Steps
==========

You now know how to:

- Run step function coupling experiments
- Handle cavity particle jumps in finite-q mode
- Implement various coupling protocols
- Analyze non-equilibrium dynamics
- Troubleshoot common issues

Continue to:

- :doc:`../part2_theory/time_varying_coupling` for theoretical details
- :doc:`../part2_theory/energy_conservation` for energy conservation theory
- :doc:`../part3_advanced/fdr_temperature` for non-equilibrium temperature measurement
- :doc:`analysis_tools` for general analysis techniques

