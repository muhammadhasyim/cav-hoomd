=================
FDR Temperature
=================

Fluctuation-Dissipation Ratio (FDR) temperature measurement for mode-specific effective temperatures.

.. contents:: In this Section
   :local:
   :depth: 2

Overview
========

FDR provides physics-based effective temperature measurement:

.. math::

   T_{\text{eff}}(\omega_0,t) = \frac{\omega_0}{2k_B} \times \frac{S_{AA}(\omega_0,t)}{\chi''(\omega_0,t)}

Where:

- :math:`S_{AA}`: Power spectral density  
- :math:`\chi''`: Imaginary susceptibility
- Measures violations of equilibrium fluctuation-dissipation theorem

Implementation
==============

See existing documentation at :doc:`../FDR_TEMPERATURE_USAGE` for detailed implementation guide.

Next: :doc:`molecular_temperatures`
