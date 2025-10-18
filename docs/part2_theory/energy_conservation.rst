====================
Energy Conservation
====================

This section provides detailed theory of energy conservation in cavity-coupled molecular dynamics.

.. contents:: In this Section
   :local:
   :depth: 2

Total Energy Decomposition
===========================

Complete Energy Expression
--------------------------

The total energy of a cavity-coupled molecular system is:

.. math::

   E_{\text{total}} = E_{\text{kinetic}}^{\text{mol}} + E_{\text{kinetic}}^{\text{cav}} + E_{\text{potential}} + E_{\text{cavity}} + E_{\text{coupling}} + E_{\text{self}}

**Molecular kinetic energy:**

.. math::

   E_{\text{kinetic}}^{\text{mol}} = \sum_{n,j} \frac{1}{2}M_{nj}\dot{R}_{nj}^2

**Cavity kinetic energy:**

.. math::

   E_{\text{kinetic}}^{\text{cav}} = \sum_{\lambda} \frac{1}{2}m_{0,\lambda}\dot{\tilde{q}}_{0,\lambda}^2

**Potential energy:**

.. math::

   E_{\text{potential}} = V(\{R_{nj}\})

Includes bonds, angles, dihedrals, Lennard-Jones, Coulomb, etc.

**Cavity harmonic energy:**

.. math::

   E_{\text{cavity}} = \sum_{\lambda} \frac{1}{2}K_\lambda \tilde{q}_{0,\lambda}^2

**Coupling energy:**

.. math::

   E_{\text{coupling}} = \sum_{\lambda} g(t)\tilde{\varepsilon}_{0,\lambda} \tilde{q}_{0,\lambda} D_\lambda

**Self-energy:**

.. math::

   E_{\text{self}} = \sum_{\lambda} \frac{g(t)^2\tilde{\varepsilon}_{0,\lambda}^2}{2K_\lambda} D_\lambda^2

Conservation Laws
=================

Microcanonical (NVE) Ensemble
------------------------------

**Without thermostats:**

.. math::

   \frac{dE_{\text{total}}}{dt} = 0

Energy is exactly conserved (within numerical precision).

**Numerical conservation:**

.. math::

   \Delta E_{\text{rel}} = \frac{|E(t) - E(0)|}{|E(0)|} < \epsilon

Typically :math:`\epsilon < 10^{-4}` for good integrators.

Canonical (NVT) Ensemble
-------------------------

**With thermostats:**

.. math::

   \frac{dE_{\text{total}}}{dt} = P_{\text{thermostat}}

where :math:`P_{\text{thermostat}}` is power injected/removed by thermostats.

**Temperature is conserved:**

.. math::

   \langle T \rangle = T_{\text{target}}

**Energy fluctuates:**

.. math::

   \sigma_E^2 = k_B T^2 C_V

where :math:`C_V` is heat capacity.

Energy Conservation During Coupling Switch
===========================================

For time-dependent :math:`g(t)`:

**Step function switch (t = t_switch):**

.. math::

   E_{\text{total}}(t_{\text{switch}}^+) = E_{\text{total}}(t_{\text{switch}}^-)

Total energy is conserved across the switch.

**Energy redistribution:**

Before switch (g=0):

.. math::

   E(t^-) = E_{\text{kinetic}}^{\text{mol}} + E_{\text{potential}}

After switch (g=g_0):

.. math::

   E(t^+) = E_{\text{kinetic}}^{\text{mol}} + E_{\text{potential}} + E_{\text{cavity}} + E_{\text{coupling}} + E_{\text{self}}

**Conservation requires:**

.. math::

   \Delta E_{\text{cavity}} + \Delta E_{\text{coupling}} + \Delta E_{\text{self}} = 0

This is automatically satisfied by the cavity particle position jump in finite-q mode.

Numerical Energy Drift
======================

Sources of Energy Drift
-----------------------

**1. Integration error:**

Finite timestep :math:`\Delta t` introduces truncation error:

.. math::

   \Delta E \propto (\Delta t)^p

where p is integrator order (p=2 for Velocity Verlet).

**2. Force discontinuities:**

- Cutoffs in pair potentials
- Bond breaking/formation
- Domain decomposition boundaries

**3. Thermostat artifacts:**

- Improper thermostat implementation
- Too-strong coupling

**4. Numerical precision:**

- Round-off errors
- Accumulation over long simulations

Minimizing Energy Drift
------------------------

**1. Reduce timestep:**

.. math::

   \Delta t < 0.001 \text{ ps (1 fs)}

For systems with high-frequency modes.

**2. Use smooth potentials:**

- Continuous forces (no sharp cutoffs)
- Switching functions for long-range interactions

**3. Proper thermostats:**

- Use well-tested implementations (Bussi, Langevin)
- Reasonable coupling constants

**4. Monitor energy:**

Regular checks during simulation:

.. code-block:: python

   energy_drift = abs(E[-1] - E[0]) / abs(E[0])
   if energy_drift > 1e-3:
       print("WARNING: Significant energy drift!")

Energy Analysis in Practice
============================

Monitoring Energy Components
-----------------------------

**Track all components:**

.. code-block:: python

   import numpy as np
   
   data = np.loadtxt('energy_tracker.txt', skiprows=1, delimiter='\t')
   
   time = data[:, 0]
   E_kin = data[:, 1]
   E_pot = data[:, 2]
   E_cav = data[:, 3]
   E_coup = data[:, 4]
   E_self = data[:, 5]
   E_tot = data[:, 6]
   
   # Check conservation
   drift = (E_tot[-1] - E_tot[0]) / E_tot[0]
   print(f"Energy drift: {drift*100:.4f}%")

**Expected behavior:**

- NVE: Total energy constant
- NVT: Total energy fluctuates, temperature constant
- Components exchange energy dynamically

Energy Equipartition
--------------------

**At equilibrium:**

.. math::

   \langle E_{\text{kinetic}}^{\text{mol}} \rangle = \frac{3}{2}Nk_BT

**Check:**

.. code-block:: python

   T_kinetic = 2 * np.mean(E_kin) / (3 * N * k_B)
   T_target = 100.0  # K
   
   if abs(T_kinetic - T_target) / T_target < 0.05:
       print("Equipartition satisfied")

Virial and Pressure
===================

**Virial theorem:**

For bound systems at equilibrium:

.. math::

   2\langle E_{\text{kinetic}} \rangle = -\langle \vec{r} \cdot \vec{F} \rangle

**Cavity contribution to pressure:**

.. math::

   P_{\text{cavity}} = \frac{1}{3V}\left\langle \sum_\lambda K_\lambda \tilde{q}_{0,\lambda}^2 \right\rangle

Usually negligible compared to molecular pressure.

Next Sections
=============

- :doc:`../part1_application/analysis_tools` for practical energy analysis
- :doc:`time_varying_coupling` for energy during coupling switches
- :doc:`thermostats` for energy exchange with thermostats
