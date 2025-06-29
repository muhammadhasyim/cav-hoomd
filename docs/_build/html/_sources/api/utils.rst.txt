==========
Utilities
==========

This module provides utility functions, constants, and helper classes for cavity MD simulations.

.. currentmodule:: hoomd.cavitymd

Physical Constants
==================

PhysicalConstants
-----------------

.. autoclass:: PhysicalConstants
   :members:
   :show-inheritance:

   Collection of physical constants and unit conversion utilities.

   .. rubric:: Key Constants

   .. autosummary::
      :nosignatures:

      ~PhysicalConstants.KB_HARTREE_PER_K
      ~PhysicalConstants.HARTREE_TO_CM_MINUS1
      ~PhysicalConstants.ATOMIC_TIME_TO_PS
      ~PhysicalConstants.BOHR_TO_ANGSTROM

   .. rubric:: Unit Conversions

   .. autosummary::
      :nosignatures:

      ~PhysicalConstants.atomic_units_to_ps
      ~PhysicalConstants.ps_to_atomic_units
      ~PhysicalConstants.gamma_from_tau_ps
      ~PhysicalConstants.tau_from_gamma_ps

Helper Functions
================

unwrap_positions
----------------

.. autofunction:: unwrap_positions

   Unwrap particle positions across periodic boundaries.

Usage Examples
==============

Unit Conversions
-----------------

.. code-block:: python

   from hoomd.cavitymd import PhysicalConstants

   # Convert time units
   dt_ps = 1.0  # 1 picosecond
   dt_au = PhysicalConstants.ps_to_atomic_units(dt_ps)
   print(f"{dt_ps} ps = {dt_au:.6f} a.u.")

   # Convert back
   dt_ps_back = PhysicalConstants.atomic_units_to_ps(dt_au)
   print(f"{dt_au:.6f} a.u. = {dt_ps_back:.6f} ps")

   # Frequency conversions
   freq_cm = 2000.0  # cm⁻¹
   freq_au = freq_cm / PhysicalConstants.HARTREE_TO_CM_MINUS1
   print(f"{freq_cm} cm⁻¹ = {freq_au:.6f} a.u.")

   # Temperature conversions
   temp_k = 300.0  # Kelvin
   temp_au = temp_k * PhysicalConstants.KB_HARTREE_PER_K
   print(f"{temp_k} K = {temp_au:.6f} a.u.")

Thermostat Parameters
---------------------

.. code-block:: python

   from hoomd.cavitymd import PhysicalConstants

   # Convert damping time to gamma
   tau_ps = 10.0  # ps
   gamma_au = PhysicalConstants.gamma_from_tau_ps(tau_ps)
   print(f"tau = {tau_ps} ps → gamma = {gamma_au:.6f} a.u.⁻¹")

   # Convert back
   tau_ps_back = PhysicalConstants.tau_from_gamma_ps(gamma_au)
   print(f"gamma = {gamma_au:.6f} a.u.⁻¹ → tau = {tau_ps_back:.6f} ps")

Position Unwrapping
--------------------

.. code-block:: python

   from hoomd.cavitymd import unwrap_positions
   import numpy as np

   # Example particle data
   positions = np.array([
       [1.0, 2.0, 3.0],
       [4.0, 5.0, 6.0],
       [7.0, 8.0, 9.0]
   ])

   images = np.array([
       [0, 0, 0],
       [1, 0, 0],
       [0, 1, 0]
   ])

   box_lengths = np.array([10.0, 10.0, 10.0])

   # Unwrap positions
   unwrapped = unwrap_positions(positions, images, box_lengths)
   print("Unwrapped positions:")
   print(unwrapped)

   # This accounts for particles that have moved across periodic boundaries

Working with Snapshots
-----------------------

.. code-block:: python

   from hoomd.cavitymd import unwrap_positions

   # In a simulation context
   with sim.state.cpu_local_snapshot as snap:
       # Get unwrapped positions for analysis
       unwrapped_pos = unwrap_positions(
           snap.particles.position,
           snap.particles.image,
           snap.configuration.box[:3]
       )
       
       # Calculate center of mass
       masses = snap.particles.mass
       com = np.average(unwrapped_pos, weights=masses, axis=0)
       print(f"Center of mass: {com}")
       
       # Calculate dipole moment
       charges = snap.particles.charge
       dipole = np.sum(charges[:, np.newaxis] * unwrapped_pos, axis=0)
       print(f"Dipole moment: {dipole}")

Constants Reference
===================

Physical Constants
------------------

.. list-table:: Key Physical Constants
   :header-rows: 1
   :widths: 30 20 50

   * - Constant
     - Value
     - Description
   * - KB_HARTREE_PER_K
     - 3.1668115e-6
     - Boltzmann constant (Hartree/K)
   * - HARTREE_TO_CM_MINUS1
     - 219474.63
     - Hartree to cm⁻¹ conversion
   * - ATOMIC_TIME_TO_PS
     - 2.418884e-5
     - Atomic time to picoseconds
   * - BOHR_TO_ANGSTROM
     - 0.529177
     - Bohr to Angstrom conversion

Unit Conversion Examples
------------------------

.. list-table:: Common Unit Conversions
   :header-rows: 1
   :widths: 30 30 40

   * - Quantity
     - Conversion
     - Example
   * - Energy
     - 1 Hartree = 27.211 eV
     - E_eV = E_hartree * 27.211
   * - Frequency
     - 1 Hartree = 219474.63 cm⁻¹
     - ω_cm⁻¹ = ω_hartree * 219474.63
   * - Time
     - 1 a.u. = 0.024189 fs
     - t_fs = t_au * 0.024189
   * - Length
     - 1 Bohr = 0.529177 Å
     - r_ang = r_bohr * 0.529177
   * - Temperature
     - 1 K = 3.1668e-6 Hartree
     - T_au = T_kelvin * 3.1668e-6

Practical Usage
===============

Setting Up Simulations
-----------------------

.. code-block:: python

   from hoomd.cavitymd import PhysicalConstants

   # Define simulation parameters in convenient units
   temperature_K = 300.0      # Kelvin
   cavity_freq_cm = 2000.0    # cm⁻¹
   runtime_ps = 1000.0        # picoseconds
   timestep_fs = 1.0          # femtoseconds

   # Convert to atomic units for HOOMD
   kT = temperature_K * PhysicalConstants.KB_HARTREE_PER_K
   omegac = cavity_freq_cm / PhysicalConstants.HARTREE_TO_CM_MINUS1
   dt = PhysicalConstants.ps_to_atomic_units(timestep_fs / 1000.0)

   print(f"Simulation parameters in atomic units:")
   print(f"  kT = {kT:.6f}")
   print(f"  ωc = {omegac:.6f}")
   print(f"  dt = {dt:.6f}")

Data Analysis
-------------

.. code-block:: python

   from hoomd.cavitymd import PhysicalConstants
   import numpy as np

   # Convert simulation data to physical units
   def convert_energy_data(data_au):
       """Convert energy from atomic units to eV."""
       return data_au * 27.211  # Hartree to eV

   def convert_time_data(data_au):
       """Convert time from atomic units to ps."""
       return PhysicalConstants.atomic_units_to_ps(data_au)

   def convert_frequency_data(data_au):
       """Convert frequency from atomic units to cm⁻¹."""
       return data_au * PhysicalConstants.HARTREE_TO_CM_MINUS1

   # Example usage
   energy_au = np.array([0.1, 0.2, 0.3])  # Hartree
   time_au = np.array([100, 200, 300])    # atomic time units

   energy_ev = convert_energy_data(energy_au)
   time_ps = convert_time_data(time_au)

   print(f"Energy: {energy_au} Hartree = {energy_ev} eV")
   print(f"Time: {time_au} a.u. = {time_ps} ps")

Parameter Validation
--------------------

.. code-block:: python

   from hoomd.cavitymd import PhysicalConstants

   def validate_parameters(temperature_K, freq_cm, coupling_strength):
       """Validate simulation parameters."""
       
       # Check temperature range
       if temperature_K < 0:
           raise ValueError("Temperature must be positive")
       if temperature_K > 1000:
           print(f"Warning: High temperature ({temperature_K} K)")
       
       # Check frequency range
       if freq_cm < 100 or freq_cm > 10000:
           print(f"Warning: Unusual frequency ({freq_cm} cm⁻¹)")
       
       # Check coupling strength
       if coupling_strength < 0:
           raise ValueError("Coupling strength must be positive")
       if coupling_strength > 0.1:
           print(f"Warning: Strong coupling ({coupling_strength})")
       
       # Convert to atomic units for validation
       kT = temperature_K * PhysicalConstants.KB_HARTREE_PER_K
       omegac = freq_cm / PhysicalConstants.HARTREE_TO_CM_MINUS1
       
       print(f"Parameters validated:")
       print(f"  kT = {kT:.6f} Hartree")
       print(f"  ωc = {omegac:.6f} Hartree")
       print(f"  g = {coupling_strength:.6f}")
       
       return kT, omegac, coupling_strength

   # Example usage
   kT, omegac, g = validate_parameters(300.0, 2000.0, 1e-3) 