=============
API Reference
=============

This section contains the complete API documentation for Cavity HOOMD.

Cavity MD Module (``hoomd.cavitymd``)
=====================================

Forces
------

.. currentmodule:: hoomd.cavitymd

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   CavityForce

Analysis and Tracking
---------------------

.. currentmodule:: hoomd.cavitymd

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   Status
   ElapsedTimeTracker
   TimestepFormatter
   CavityModeTracker
   AutocorrelationTracker
   FieldAutocorrelationTracker
   EnergyTracker
   DipoleAutocorrelation

Simulation Utilities
--------------------

.. currentmodule:: hoomd.cavitymd

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   AdaptiveTimestepUpdater

Utilities
---------

.. currentmodule:: hoomd.cavitymd

.. autosummary::
   :toctree: _autosummary

   PhysicalConstants
   unwrap_positions

Bussi Reservoir Module (``hoomd.bussi_reservoir``)
===================================================

Thermostats
-----------

.. currentmodule:: hoomd.bussi_reservoir

.. autosummary::
   :toctree: _autosummary
   :template: bussi_class.rst

   BussiReservoir 