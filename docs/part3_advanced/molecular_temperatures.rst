=======================
Molecular Temperatures
=======================

Kinetic temperature decomposition for diatomic molecules into translational, rotational, and vibrational components.

.. contents:: In this Section
   :local:
   :depth: 2

Temperature Decomposition
==========================

For diatomic molecules:

**Translational:**

.. math::

   T_{\text{trans}} = \frac{2\langle KE_{\text{CM}}\rangle}{3Nk_B}

**Rotational:**

.. math::

   T_{\text{rot}} = \frac{\langle KE_{\text{rot}}\rangle}{Nk_B}

**Vibrational:**

.. math::

   T_{\text{vib}} = \frac{2\langle KE_{\text{vib}}\rangle}{Nk_B}

Implementation
==============

.. code-block:: python

   from cavitymd.molecular_temperatures import DiatomicMolecularTemperatures
   
   temp_tracker = DiatomicMolecularTemperatures(
       simulation=sim,
       time_tracker=time_tracker,
       output_period_ps=1.0,
       output_file='molecular_temps.csv'
   )
   
   sim.operations.writers.append(temp_tracker)

See :doc:`../molecular_temperatures` for detailed documentation.

Next: :doc:`correlation_analysis`
