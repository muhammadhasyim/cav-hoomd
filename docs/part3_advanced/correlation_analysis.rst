=====================
Correlation Analysis
=====================

Computing autocorrelation functions and spectral analysis for cavity-coupled systems.

.. contents:: In this Section
   :local:
   :depth: 2

Dipole Autocorrelation
=======================

**Definition:**

.. math::

   C(t) = \langle \vec{D}(t') \cdot \vec{D}(t'+t) \rangle

Used for computing IR spectra and relaxation times.

IR Spectrum Calculation
========================

**Fourier transform:**

.. math::

   I(\omega) \propto \omega^2 \int_0^\infty C(t) e^{-i\omega t} dt

See :doc:`../part1_application/analysis_tools` for implementation details.

Field Autocorrelation
=====================

For cavity mode analysis:

.. math::

   C_q(t) = \langle q(t') q(t'+t) \rangle

Reveals cavity mode dynamics and damping.

Next: :doc:`performance`
