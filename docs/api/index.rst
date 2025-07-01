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

   CavityForce

Analysis and Tracking
--------------------

.. autosummary::
   :toctree: _autosummary

   Status
   ElapsedTimeTracker
   TimestepFormatter
   CavityModeTracker
   AutocorrelationTracker
   FieldAutocorrelationTracker
   EnergyTracker
   DipoleAutocorrelation

Simulation Utilities
-------------------

.. autosummary::
   :toctree: _autosummary

   AdaptiveTimestepUpdater

Utilities
---------

.. autosummary::
   :toctree: _autosummary

   PhysicalConstants
   unwrap_positions

Bussi Reservoir Module (``hoomd.bussi_reservoir``)
=================================================

Thermostats
-----------

.. currentmodule:: hoomd.bussi_reservoir

.. autosummary::
   :toctree: _autosummary

   BussiReservoir

Detailed API Documentation
==========================

.. toctree::
   :maxdepth: 2

   _autosummary/hoomd.cavitymd.CavityForce
   _autosummary/hoomd.cavitymd.Status
   _autosummary/hoomd.cavitymd.ElapsedTimeTracker
   _autosummary/hoomd.cavitymd.TimestepFormatter
   _autosummary/hoomd.cavitymd.CavityModeTracker
   _autosummary/hoomd.cavitymd.AutocorrelationTracker
   _autosummary/hoomd.cavitymd.FieldAutocorrelationTracker
   _autosummary/hoomd.cavitymd.EnergyTracker
   _autosummary/hoomd.cavitymd.DipoleAutocorrelation
   _autosummary/hoomd.cavitymd.AdaptiveTimestepUpdater
   _autosummary/hoomd.cavitymd.PhysicalConstants
   _autosummary/hoomd.cavitymd.unwrap_positions
   _autosummary/hoomd.bussi_reservoir.BussiReservoir 