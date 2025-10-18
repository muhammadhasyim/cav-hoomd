==============
Analysis Tools
==============

This guide covers the analysis tools available in Cavity HOOMD for interpreting simulation results.

.. contents:: In this Guide
   :local:
   :depth: 2

Energy Tracking and Interpretation
===================================

Energy Components
-----------------

Cavity HOOMD tracks multiple energy components:

**Energy Decomposition:**

.. math::

   E_{\text{total}} = E_{\text{kinetic}} + E_{\text{potential}} + E_{\text{cavity}} + E_{\text{coupling}} + E_{\text{self}}

Where:

- :math:`E_{\text{kinetic}}`: Molecular kinetic energy
- :math:`E_{\text{potential}}`: Molecular potential energy (LJ, bonds, etc.)
- :math:`E_{\text{cavity}} = \frac{1}{2}K(q_x^2 + q_y^2)`: Cavity harmonic energy
- :math:`E_{\text{coupling}} = g(q_x D_x + q_y D_y)`: Dipole-field coupling
- :math:`E_{\text{self}} = \frac{g^2}{2K}(D_x^2 + D_y^2)`: Dipole self-energy

Reading Energy Files
--------------------

**Energy tracker output format:**

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   # Load energy data
   data = np.loadtxt('prod-1_energy_tracker.txt', skiprows=1, delimiter='\t')
   
   # Extract columns
   time_ps = data[:, 0]
   E_kinetic = data[:, 1]
   E_potential = data[:, 2]
   E_cavity = data[:, 3]
   E_coupling = data[:, 4]
   E_self = data[:, 5]
   E_total = data[:, 6]
   
   # Plot energy components
   fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
   
   # Total energy (conservation check)
   ax1.plot(time_ps, E_total, 'k-', linewidth=2)
   ax1.set_ylabel('Total Energy')
   ax1.set_title('Energy Conservation')
   ax1.grid(True, alpha=0.3)
   
   # Energy components
   ax2.plot(time_ps, E_kinetic, label='Kinetic')
   ax2.plot(time_ps, E_potential, label='Potential')
   ax2.plot(time_ps, E_cavity, label='Cavity')
   ax2.plot(time_ps, E_coupling, label='Coupling')
   ax2.plot(time_ps, E_self, label='Self-energy')
   ax2.set_xlabel('Time (ps)')
   ax2.set_ylabel('Energy')
   ax2.set_title('Energy Components')
   ax2.legend()
   ax2.grid(True, alpha=0.3)
   
   plt.tight_layout()
   plt.savefig('energy_analysis.png', dpi=300)
   plt.show()

Interpreting Energy Behavior
-----------------------------

**Normal Behavior:**

1. **Total energy conservation:**
   
   In NVE or during unperturbed NVT:
   
   .. math::
   
      \frac{\Delta E_{\text{total}}}{E_{\text{total}}} < 0.01\%

2. **Kinetic-potential exchange:**
   
   Kinetic and potential energies oscillate in anti-phase:
   
   .. code-block:: python
   
      correlation = np.corrcoef(E_kinetic, E_potential)[0, 1]
      print(f"KE-PE correlation: {correlation:.3f}")  # Should be negative

3. **Coupling energy:**
   
   - Magnitude: :math:`|E_{\text{coupling}}| \sim g \sqrt{N}`
   - Sign: Can be positive or negative
   - Fluctuations: Indicates molecular-cavity energy exchange

**Warning Signs:**

1. **Energy drift:**
   
   .. code-block:: python
   
      drift_percent = (E_total[-1] - E_total[0]) / E_total[0] * 100
      if abs(drift_percent) > 0.1:
           print(f"WARNING: Energy drift {drift_percent:.2f}%")

   **Causes:** Timestep too large, numerical instability, bad initial configuration.

2. **Energy spikes:**
   
   Sudden jumps indicate particle overlaps or force discontinuities.

3. **Unbounded growth:**
   
   Cavity energy continuously increasing suggests missing dissipation.

Cavity Mode Trajectory Analysis
================================

Loading Cavity Mode Data
-------------------------

**Cavity mode file format:**

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   # Load cavity mode data
   data = np.loadtxt('prod-1_cavity_mode.txt', skiprows=1, delimiter='\t')
   
   time_ps = data[:, 0]
   q_x = data[:, 1]  # X-polarization position
   q_y = data[:, 2]  # Y-polarization position
   v_x = data[:, 3]  # X-polarization velocity
   v_y = data[:, 4]  # Y-polarization velocity
   
   # Calculate total amplitude
   q_total = np.sqrt(q_x**2 + q_y**2)
   v_total = np.sqrt(v_x**2 + v_y**2)

Cavity Mode Oscillations
-------------------------

**Analyzing oscillation frequency:**

.. code-block:: python

   from scipy import signal
   
   # FFT analysis
   sampling_rate = 1.0 / (time_ps[1] - time_ps[0])  # Hz
   frequencies, power = signal.welch(q_x, fs=sampling_rate, nperseg=1024)
   
   # Find peak frequency
   peak_idx = np.argmax(power)
   peak_freq_Hz = frequencies[peak_idx]
   peak_freq_cm = peak_freq_Hz * 33.356  # Convert to cm⁻¹
   
   print(f"Cavity oscillation frequency: {peak_freq_cm:.1f} cm⁻¹")
   
   # Plot power spectrum
   plt.figure(figsize=(10, 5))
   plt.semilogy(frequencies * 33.356, power)
   plt.xlabel('Frequency (cm⁻¹)')
   plt.ylabel('Power Spectral Density')
   plt.title('Cavity Mode Spectrum')
   plt.xlim(0, 5000)
   plt.grid(True, alpha=0.3)
   plt.show()

**Phase Space Trajectory:**

.. code-block:: python

   # Plot phase space (position vs velocity)
   fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
   
   # X-polarization
   ax1.plot(q_x, v_x, alpha=0.5, linewidth=0.5)
   ax1.set_xlabel('Position q_x')
   ax1.set_ylabel('Velocity v_x')
   ax1.set_title('X-Polarization Phase Space')
   ax1.grid(True, alpha=0.3)
   
   # Y-polarization
   ax2.plot(q_y, v_y, alpha=0.5, linewidth=0.5)
   ax2.set_xlabel('Position q_y')
   ax2.set_ylabel('Velocity v_y')
   ax2.set_title('Y-Polarization Phase Space')
   ax2.grid(True, alpha=0.3)
   
   plt.tight_layout()
   plt.show()

Cavity-Molecule Energy Exchange
--------------------------------

**Correlation analysis:**

.. code-block:: python

   # Calculate cavity energy
   K = 1.0  # Cavity spring constant (depends on frequency and mass)
   E_cavity = 0.5 * K * (q_x**2 + q_y**2)
   
   # Load molecular kinetic energy from energy tracker
   E_kinetic = np.loadtxt('prod-1_energy_tracker.txt', 
                           skiprows=1, delimiter='\t', usecols=1)
   
   # Cross-correlation
   correlation = np.correlate(E_cavity - np.mean(E_cavity),
                             E_kinetic - np.mean(E_kinetic), 
                             mode='same')
   
   lags = np.arange(-len(correlation)//2, len(correlation)//2)
   lag_times = lags * (time_ps[1] - time_ps[0])
   
   plt.figure(figsize=(10, 5))
   plt.plot(lag_times, correlation / np.max(correlation))
   plt.xlabel('Lag Time (ps)')
   plt.ylabel('Normalized Cross-Correlation')
   plt.title('Cavity-Molecular Energy Exchange')
   plt.grid(True, alpha=0.3)
   plt.show()

Temperature Decomposition
=========================

For diatomic molecules, temperature can be decomposed into translational, rotational, and vibrational components.

Using DiatomicMolecularTemperatures
------------------------------------

**Setup in simulation:**

.. code-block:: python

   from cavitymd.molecular_temperatures import DiatomicMolecularTemperatures
   from cavitymd.analysis import ElapsedTimeTracker
   
   # Create trackers
   time_tracker = ElapsedTimeTracker(simulation=sim, runtime_ps=1000.0)
   temp_tracker = DiatomicMolecularTemperatures(
       simulation=sim,
       time_tracker=time_tracker,
       output_period_ps=1.0,
       output_file='molecular_temps.csv'
   )
   
   # Add to simulation
   sim.operations.writers.append(time_tracker)
   sim.operations.writers.append(temp_tracker)

**Analyzing output:**

.. code-block:: python

   import pandas as pd
   
   # Load temperature data
   df = pd.read_csv('molecular_temps.csv')
   
   # Extract columns
   time = df['time_ps']
   T_trans = df['T_translational_K']
   T_rot = df['T_rotational_K']
   T_vib = df['T_vibrational_K']
   T_total = df['T_total_K']
   
   # Plot
   fig, ax = plt.subplots(figsize=(12, 6))
   ax.plot(time, T_trans, label='Translational', linewidth=2)
   ax.plot(time, T_rot, label='Rotational', linewidth=2)
   ax.plot(time, T_vib, label='Vibrational', linewidth=2)
   ax.plot(time, T_total, label='Total (Kinetic)', 
           linestyle='--', color='black', linewidth=2)
   ax.axhline(100.0, color='red', linestyle=':', 
              label='Target Temperature')
   ax.set_xlabel('Time (ps)', fontsize=12)
   ax.set_ylabel('Temperature (K)', fontsize=12)
   ax.set_title('Molecular Temperature Decomposition', fontsize=14)
   ax.legend(fontsize=10)
   ax.grid(True, alpha=0.3)
   plt.tight_layout()
   plt.savefig('temperature_decomposition.png', dpi=300)
   plt.show()

Equipartition Check
-------------------

**Theory:** At equilibrium, each quadratic degree of freedom should have energy :math:`\frac{1}{2}k_B T`:

.. math::

   \langle E_{\text{trans}} \rangle &= \frac{3}{2} N k_B T \\
   \langle E_{\text{rot}} \rangle &= N k_B T \\
   \langle E_{\text{vib}} \rangle &= \frac{1}{2} N k_B T

**Checking equipartition:**

.. code-block:: python

   # Calculate average temperatures (after equilibration)
   equilibration_time = 100.0  # ps
   eq_mask = time > equilibration_time
   
   T_trans_avg = T_trans[eq_mask].mean()
   T_rot_avg = T_rot[eq_mask].mean()
   T_vib_avg = T_vib[eq_mask].mean()
   T_target = 100.0
   
   print(f"Target temperature: {T_target:.1f} K")
   print(f"Translational: {T_trans_avg:.1f} K ({T_trans_avg/T_target*100:.1f}%)")
   print(f"Rotational: {T_rot_avg:.1f} K ({T_rot_avg/T_target*100:.1f}%)")
   print(f"Vibrational: {T_vib_avg:.1f} K ({T_vib_avg/T_target*100:.1f}%)")
   
   # Check equipartition (should all be close to target)
   deviation_trans = abs(T_trans_avg - T_target) / T_target * 100
   deviation_rot = abs(T_rot_avg - T_target) / T_target * 100
   deviation_vib = abs(T_vib_avg - T_target) / T_target * 100
   
   if max(deviation_trans, deviation_rot, deviation_vib) < 5.0:
       print("✓ Equipartition satisfied (within 5%)")
   else:
       print("✗ Equipartition violation detected")

Correlation Functions and IR Spectra
====================================

Dipole Autocorrelation Function
--------------------------------

**Computing autocorrelation:**

.. code-block:: python

   import gsd.hoomd
   
   # Load trajectory
   traj = gsd.hoomd.open('prod-1.gsd')
   
   # Extract dipole moments
   dipole_moments = []
   for frame in traj:
       # Calculate total dipole (simplified for O2)
       positions = frame.particles.position[:-2]  # Exclude cavity particles
       types = frame.particles.typeid[:-2]
       
       # For diatomic O2, dipole ~ separation vector
       N_molecules = len(positions) // 2
       dipoles = []
       for i in range(N_molecules):
           r1 = positions[2*i]
           r2 = positions[2*i + 1]
           dipole = r2 - r1
           dipoles.append(dipole)
       
       total_dipole = np.sum(dipoles, axis=0)
       dipole_moments.append(total_dipole)
   
   dipole_moments = np.array(dipole_moments)
   
   # Compute autocorrelation
   def autocorr(x):
       result = np.correlate(x, x, mode='full')
       return result[len(result)//2:] / result[len(result)//2]
   
   # Autocorrelation for each component
   C_x = autocorr(dipole_moments[:, 0])
   C_y = autocorr(dipole_moments[:, 1])
   C_z = autocorr(dipole_moments[:, 2])
   
   # Average
   C_avg = (C_x + C_y + C_z) / 3.0
   
   # Time axis
   dt = 0.01  # ps (output frequency)
   times = np.arange(len(C_avg)) * dt
   
   # Plot
   plt.figure(figsize=(10, 5))
   plt.plot(times[:1000], C_avg[:1000])  # First 10 ps
   plt.xlabel('Time (ps)')
   plt.ylabel('Autocorrelation C(t)')
   plt.title('Dipole Moment Autocorrelation')
   plt.grid(True, alpha=0.3)
   plt.show()

IR Spectrum Calculation
------------------------

**Fourier transform of autocorrelation:**

.. code-block:: python

   from scipy import signal
   
   # Window function to reduce ringing
   window = signal.windows.blackman(len(C_avg))
   C_windowed = C_avg * window
   
   # FFT
   frequencies = np.fft.rfftfreq(len(C_windowed), d=dt)
   spectrum = np.fft.rfft(C_windowed)
   intensity = np.abs(spectrum)**2
   
   # Convert frequency to cm⁻¹
   freq_cm = frequencies * 33.356  # Conversion factor
   
   # Plot IR spectrum
   plt.figure(figsize=(12, 5))
   plt.plot(freq_cm, intensity)
   plt.xlabel('Frequency (cm⁻¹)')
   plt.ylabel('Intensity (arbitrary units)')
   plt.title('Infrared Absorption Spectrum')
   plt.xlim(0, 3000)
   plt.grid(True, alpha=0.3)
   plt.savefig('ir_spectrum.png', dpi=300)
   plt.show()

**Identifying peaks:**

.. code-block:: python

   # Find peaks in spectrum
   peaks, properties = signal.find_peaks(intensity, 
                                         height=np.max(intensity)*0.1,
                                         prominence=np.max(intensity)*0.05)
   
   print("Vibrational frequencies:")
   for i, peak_idx in enumerate(peaks):
       freq = freq_cm[peak_idx]
       height = intensity[peak_idx]
       print(f"  Peak {i+1}: {freq:.1f} cm⁻¹ (intensity: {height:.2e})")

Plotting Utilities and Visualization
=====================================

Trajectory Visualization
------------------------

**Using OVITO (external tool):**

.. code-block:: bash

   ovito prod-1.gsd

**Using MDAnalysis (Python):**

.. code-block:: python

   import MDAnalysis as mda
   from MDAnalysis.analysis import rdf
   
   # Load trajectory
   u = mda.Universe('init-0.gsd', 'prod-1.gsd')
   
   # Select atoms
   oxygens = u.select_atoms('type O')
   
   # Calculate radial distribution function
   rdf_analysis = rdf.InterRDF(oxygens, oxygens, 
                                nbins=100, range=(0.0, 10.0))
   rdf_analysis.run()
   
   # Plot
   plt.figure(figsize=(10, 5))
   plt.plot(rdf_analysis.results.bins, rdf_analysis.results.rdf)
   plt.xlabel('Distance (Å)')
   plt.ylabel('g(r)')
   plt.title('Radial Distribution Function')
   plt.grid(True, alpha=0.3)
   plt.show()

Creating Publication-Quality Figures
-------------------------------------

**Template for high-quality plots:**

.. code-block:: python

   import matplotlib.pyplot as plt
   import matplotlib as mpl
   
   # Set publication style
   mpl.rcParams['font.size'] = 12
   mpl.rcParams['font.family'] = 'sans-serif'
   mpl.rcParams['axes.linewidth'] = 1.5
   mpl.rcParams['xtick.major.width'] = 1.5
   mpl.rcParams['ytick.major.width'] = 1.5
   mpl.rcParams['xtick.major.size'] = 5
   mpl.rcParams['ytick.major.size'] = 5
   
   fig, ax = plt.subplots(figsize=(6, 4))
   
   # Plot data
   ax.plot(x_data, y_data, 'o-', linewidth=2, markersize=4, 
           label='Cavity-coupled')
   ax.plot(x_control, y_control, 's--', linewidth=2, markersize=4,
           label='Control')
   
   # Labels and formatting
   ax.set_xlabel('Time (ps)', fontsize=14)
   ax.set_ylabel('Observable', fontsize=14)
   ax.legend(frameon=True, fontsize=12, loc='best')
   ax.grid(True, alpha=0.3, linestyle='--')
   
   # Tight layout
   plt.tight_layout()
   
   # Save
   plt.savefig('figure.pdf', dpi=300, bbox_inches='tight')
   plt.savefig('figure.png', dpi=300, bbox_inches='tight')
   plt.show()

Next Steps
==========

You now know how to:

- Track and interpret energy components
- Analyze cavity mode trajectories
- Decompose molecular temperatures
- Calculate correlation functions and IR spectra
- Create publication-quality visualizations

Continue to:

- :doc:`time_varying_coupling` for time-dependent coupling analysis
- :doc:`../part3_advanced/fdr_temperature` for advanced temperature measurements
- :doc:`../part3_advanced/correlation_analysis` for detailed correlation analysis
- :doc:`../part2_theory/energy_conservation` for energy conservation theory

