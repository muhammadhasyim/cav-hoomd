# Copyright (c) 2009-2025 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Extended Bussi thermostat with reservoir energy tracking."""

from hoomd.md import _md
from hoomd.md.methods import thermostats
import hoomd
from hoomd.operation import _HOOMDBaseObject
from hoomd.data.parameterdicts import ParameterDict
from hoomd.variant import Variant


class BussiReservoir(thermostats.Thermostat):
    r"""Extended Bussi thermostat that tracks reservoir energy.

    This thermostat extends the standard Bussi stochastic velocity rescaling
    thermostat to track the energy that is dumped into or taken from the thermal
    reservoir during the simulation.

    Args:
        kT (hoomd.variant.variant_like): Temperature set point for the
            thermostat :math:`[\mathrm{energy}]`.

        tau (float): Thermostat time constant :math:`[\mathrm{time}]`.
            Defaults to 0.

    The reservoir energy represents the cumulative energy exchange between the
    system and the thermal bath. Positive values indicate energy dumped to the
    reservoir (cooling), while negative values indicate energy taken from the
    reservoir (heating).

    .. rubric:: Example:

    .. code-block:: python

        import hoomd.bussi_reservoir as bussi_res
        
        bussi = bussi_res.thermostats.BussiReservoir(
            kT=1.5, tau=simulation.operations.integrator.dt * 20
        )
        simulation.operations.integrator.methods[0].thermostat = bussi

        # Run simulation
        simulation.run(1000)

        # Check reservoir energies
        print(f"Total reservoir energy: {bussi.total_reservoir_energy}")
        print(f"Translational: {bussi.reservoir_energy_translational}")
        print(f"Rotational: {bussi.reservoir_energy_rotational}")

    Attributes:
        kT (hoomd.variant.variant_like): Temperature set point for the
            thermostat :math:`[\mathrm{energy}]`.

        tau (float): Thermostat time constant :math:`[\mathrm{time}].`

        reservoir_energy_translational (float): Cumulative energy dumped to 
            reservoir from translational degrees of freedom 
            :math:`[\mathrm{energy}]` (read-only).

        reservoir_energy_rotational (float): Cumulative energy dumped to 
            reservoir from rotational degrees of freedom 
            :math:`[\mathrm{energy}]` (read-only).

        total_reservoir_energy (float): Total cumulative energy dumped to 
            reservoir :math:`[\mathrm{energy}]` (read-only).

        instantaneous_reservoir_translational (float): Energy dumped to 
            reservoir from translational DOF in the last timestep 
            :math:`[\mathrm{energy}]` (read-only).

        instantaneous_reservoir_rotational (float): Energy dumped to 
            reservoir from rotational DOF in the last timestep 
            :math:`[\mathrm{energy}]` (read-only).

        instantaneous_reservoir_total (float): Total energy dumped to 
            reservoir in the last timestep :math:`[\mathrm{energy}]` (read-only).
    """

    def __init__(self, kT, tau=0.0):
        super().__init__(kT)
        param_dict = ParameterDict(tau=float)
        param_dict["tau"] = tau
        self._param_dict.update(param_dict)

    def _attach_hook(self):
        from . import _bussi_reservoir
        group = self._simulation.state._get_group(self._filter)
        self._cpp_obj = _bussi_reservoir.BussiReservoirThermostat(
            self.kT, group, self._thermo, self._simulation.state._cpp_sys_def, self.tau
        )
        self._simulation._warn_if_seed_unset()

    @hoomd.logging.log()
    def reservoir_energy_translational(self):
        """Cumulative energy dumped to reservoir from translational DOF."""
        if not self._attached:
            return 0.0
        return self._cpp_obj.getReservoirEnergyTranslational()

    @hoomd.logging.log()
    def reservoir_energy_rotational(self):
        """Cumulative energy dumped to reservoir from rotational DOF."""
        if not self._attached:
            return 0.0
        return self._cpp_obj.getReservoirEnergyRotational()

    @hoomd.logging.log()
    def total_reservoir_energy(self):
        """Total cumulative energy dumped to reservoir."""
        if not self._attached:
            return 0.0
        return self._cpp_obj.getTotalReservoirEnergy()

    @hoomd.logging.log()
    def instantaneous_reservoir_translational(self):
        """Energy dumped to reservoir from translational DOF in last timestep."""
        if not self._attached:
            return 0.0
        return self._cpp_obj.getInstantaneousReservoirTranslational()

    @hoomd.logging.log()
    def instantaneous_reservoir_rotational(self):
        """Energy dumped to reservoir from rotational DOF in last timestep."""
        if not self._attached:
            return 0.0
        return self._cpp_obj.getInstantaneousReservoirRotational()

    @hoomd.logging.log()
    def instantaneous_reservoir_total(self):
        """Total energy dumped to reservoir in last timestep."""
        if not self._attached:
            return 0.0
        return self._cpp_obj.getInstantaneousReservoirTotal()

    def reset_reservoir_energy(self):
        """Reset all reservoir energy counters to zero.

        This is useful when you want to measure reservoir energy over
        a specific time period.

        .. rubric:: Example

        .. code-block:: python

            # Reset counters
            bussi.reset_reservoir_energy()
            
            # Run production simulation
            simulation.run(10000)
            
            # Measure reservoir energy for this period
            production_reservoir_energy = bussi.total_reservoir_energy
        """
        if self._attached:
            self._cpp_obj.resetReservoirEnergy()
        # If not attached, this is a no-op since energies are already 0.0


__all__ = [
    "BussiReservoir",
] 